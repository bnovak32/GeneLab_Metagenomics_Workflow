#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SAMTOOLS_FILTER_FASTQ as  NANO_REMOVE_CONTAMINANT} from "./samtools.nf"

// max_mem = 100e9 // 100GB

/*
 * ========================================================================================
 * PROCESS: SPADES
 * ========================================================================================
 *
 * SUMMARY:
 *   Assemble a group of (blank)reads with spades
 *
 * INPUTS:
 *   1. path: forward
 *      Cardinality: one
 *      Description: Input file: list of forward reads
 *
 *   2. path: reverse
 *      Cardinality: one
 *      Description: Input file: list of reverse reads
 *
 *   3. val: isPaired
 *      Cardinality: one
 *      Description: Parameter value: are the input reads paired?
 *
 *   4. val: type
 *      Cardinality: one
 *      Description: Parameter value: technology type i.e. one of illumina, pacbio or nanopore
 *
 * OUTPUTS:
 *   1. path: blank-assembly.fasta (emit: assembly) [OPTIONAL]
 *
 *   2. path: blank-warnings.log (emit: warnings) [OPTIONAL]
 *
 *   3. path: blank-assembly.log (emit: log)
 *
 *   4. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: SPAdes
 *   Container: [Defined in config/nextflow.config]
 *   Conda: envs/spades.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process SPADES {

    tag "Assembling contaminants..."


    input:
        path(forward)
        path(reverse)
        val(isPaired)
        val(type) // illumina or pacbio or nanopore

    output:
        path('blank-assembly.fasta'), optional:true, emit: assembly
        path('blank-warnings.log'), optional:true, emit: warnings
        path('blank-assembly.log'), emit: log
        path("versions.txt"), emit: version

    script:
        def maxmem = task.memory.toGiga()
    """
    if [ ${type} ==  "illumina" ];then

         if [ ${isPaired} == 'true' ];then

              zcat ${forward} > merged_R1.fastq
              zcat ${reverse} > merged_R2.fastq

              INPUT=' -1 merged_R1.fastq -2 merged_R2.fastq'
         else
              
              zcat ${forward} > merged.fastq
              INPUT=' -s merged.fastq'
         fi

    fi

    if [ ${type} ==  "pacbio" ];then

       zcat ${forward} > merged.fastq

       INPUT=' --pacbio merged.fastq'

    fi


    if [ ${type} ==  "nanopore" ];then
 
       zcat ${forward} > merged.fastq

       INPUT=' --nanopore merged.fastq'

    fi


    spades.py --meta --threads ${task.cpus} --memory ${maxmem} \${INPUT} -o .


    # Renaming output files
    mv scaffolds.fasta blank-assembly.fasta
    mv spades.log blank-assembly.log

    [ -f warnings.log ] && mv warnings.log blank-warnings.log


    VERSION=`spades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p'`
    echo "spades \${VERSION}" > versions.txt
    """
}
    

/*
 * ========================================================================================
 * PROCESS: FLYE
 * ========================================================================================
 *
 * SUMMARY:
 *   Assemble a group of (blank)reads with flye
 *
 * INPUTS:
 *   1. path: reads
 *      Cardinality: one
 *      Description: Input file: list of reads
 *
 * OUTPUTS:
 *   1. path: blank-assembly.fasta (emit: assembly)
 *
 *   2. path: blank-flye.log (emit: log)
 *
 *   3. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: Flye
 *   Container: [Defined in config/nanopore.config]
 *   Conda: envs/flye.yaml
 *   Labels: flye
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process FLYE {

    tag "Assembling contaminants.."
    label "flye"
 
    input:
       path(reads)

    output:
        path("blank-assembly.fasta"), emit: assembly
        path("blank-flye.log"), emit: log
        path("versions.txt"), emit: version

    script:
    """
    flye --meta \\
         --threads ${task.cpus} \\
         --out-dir . \\
         --nano-raw ${reads}

    mv assembly.fasta  blank-assembly.fasta
    mv flye.log blank-flye.log

    VERSION=`flye --version`
    echo "flye \${VERSION}" > versions.txt
    """
}

/*
 * ========================================================================================
 * PROCESS: BUILD_CONTAMINANT_DB
 * ========================================================================================
 *
 * SUMMARY:
 *   Build a contaminant Database from blank samples with bowtie2
 *
 * INPUTS:
 *   1. path: FASTA
 *      Cardinality: one
 *      Description: Input file: Blank assembly
 *
 * OUTPUTS:
 *   1. path: blank-index/ (emit: index)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/illumina.config]
 *   Conda: envs/bowtie2.yaml
 *   Labels: bowtie2
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process BUILD_CONTAMINANT_DB {

    tag "Building a contaminant Database from blank samples"
    label "bowtie2"

    input:
        path(FASTA)

    output:
        path("blank-index/"), emit: index
        path("versions.txt"), emit: version

    script:
        """
        mkdir blank-index/
        bowtie2-build ${FASTA} blank-index/blanks

        bowtie2 --version  | \\
             head -n 1 | \\
             sed -E 's/.*(bowtie2-align-s version.+)/\\1/' > versions.txt
        """
}

/*
 * ========================================================================================
 * PROCESS: BUILD_CONTAMINANT_INDEX
 * ========================================================================================
 *
 * SUMMARY:
 *   Build a contaminant index from blank samples with minimap2
 *
 * INPUTS:
 *   1. path: FASTA
 *      Cardinality: one
 *      Description: Input file: Blank assembly
 *
 * OUTPUTS:
 *   1. path: blanks.mmi (emit: index)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: Minimap2
 *   Container: [Defined in config/nanopore.config]
 *   Conda: envs/minimap2.yaml
 *   Labels: minimap2
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process BUILD_CONTAMINANT_INDEX {

    tag "Building a contaminant index from blank samples"
    label "minimap2"

    input:
        path(FASTA) // blank assembly
    output:
        path("blanks.mmi"), emit: index
        path("versions.txt"), emit: version
    script:
        """
        minimap2 -ax splice -t ${task.cpus} -d blanks.mmi ${FASTA}
        VERSION=`minimap2 --version`
        echo "minimap2 \${VERSION}" > versions.txt
        """
}

/*
 * ========================================================================================
 * PROCESS: REMOVE_CONTAMINANT
 * ========================================================================================
 *
 * SUMMARY:
 *   Remove reads mapping to blanks assembly with bowtie2
 *
 * INPUTS:
 *   1. each: path(BLANKS_DB)
 *      Cardinality: each
 *      Description: Iterates over each element. Blanks assembly index.
 *
 *   2. tuple: tuple val(sample_id), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Bolean specifying whether input reads are paired or not 
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("*.fastq.gz"), val(isPaired) (emit: reads)
 *
 *   2. tuple: tuple val(sample_id), path("${sample_id}-mapping-info.txt") (emit: info)
 *
 *   3. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/illumina.config]
 *   Conda: envs/bowtie2.yaml
 *   Labels: bowtie2
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process REMOVE_CONTAMINANT {

    tag "Removing reads mapping to the blanks assembly.."
    label "bowtie2"

    input:
        each path(BLANKS_DB)
        tuple val(sample_id), path(reads), val(isPaired)

    output:
        tuple val(sample_id), path("*.fastq.gz"), val(isPaired), emit: reads
        tuple val(sample_id), path("${sample_id}-mapping-info.txt"), emit: info
        path("versions.txt"), emit: version

    script:
        def input =  isPaired == "true" ?  
                 "-1 ${reads[0]} -2 ${reads[1]} --un-conc-gz" : 
                 "-U ${reads[0]} --un-gz"
        """
        INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\${INDEX}" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\${INDEX}" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 -p ${task.cpus} -x \${INDEX} --very-sensitive-local \\
               ${input} ${sample_id}_decontam${params.assay_suffix}.fastq.gz \\
               > ${sample_id}.sam 2> ${sample_id}-mapping-info.txt 

        # Rename Fastq Files
        if [ -f ${sample_id}_decontam${params.assay_suffix}.fastq.1.gz ]; then
            mv ${sample_id}_decontam${params.assay_suffix}.fastq.1.gz ${sample_id}_R1_decontam${params.assay_suffix}.fastq.gz
        fi

        if [ -f ${sample_id}_decontam${params.assay_suffix}.fastq.2.gz ]; then
            mv ${sample_id}_decontam${params.assay_suffix}.fastq.2.gz ${sample_id}_R2_decontam${params.assay_suffix}.fastq.gz
        fi


        # Delete sam file to free up memory
        rm -rf ${sample_id}.sam 

        bowtie2 --version  | \\
             head -n 1 | \\
             sed -E 's/.*(bowtie2-align-s version.+)/\\1/' > versions.txt
    """
}


/*
 * ========================================================================================
 * PROCESS: MAPPING_TO_CONTAMINANT
 * ========================================================================================
 *
 * SUMMARY:
 *   Map reads to blanks assembly with minimap2
 *
 * INPUTS:
 *   1. each: path(BLANKS_DB)
 *      Cardinality: each
 *      Description: Iterates over each element. Blanks minimap index.
 *
 *   2. tuple: tuple val(sample_id), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Bolean specifying whether input reads are paired or not 
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}-mapping-info.txt") (emit: sam)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/nanopore.config]
 *   Conda: envs/minimap2.yaml
 *   Labels: minimap2
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// This process builds the bowtie2 index and runs the mapping for each sample
process MAPPING_TO_CONTAMINANT {

    tag "mapping reads to the blanks assembly..."
    label "minimap2"
    

    input:
        each path(BLANKS_DB)
        tuple val(sample_id), path(reads), val(isPaired)
    output:
        tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}-mapping-info.txt"), emit: sam
        path("versions.txt"), emit: version
    script:
        """
        INDEX=`find -L ./ -name "*.mmi"`
        minimap2 -t ${task.cpus} -ax splice \${INDEX} ${reads} \\
                     > ${sample_id}.sam  2> ${sample_id}-mapping-info.txt 
        VERSION=`minimap2 --version`
        echo "minimap2 \${VERSION}" > versions.txt
        """
}



// A function to delete white spaces from an input string and covert it to lower case
def deleteWS(string){

    return string.replaceAll(/\s+/, '').toLowerCase()

}


workflow illumina_remove_contaminants {


    take:
        file_ch
        filtered_ch

    main:

        file_ch.map{row -> if(deleteWS(row.NTC) == "true") [row.sample_id] }.set{blank_samples}

        isPaired  = filtered_ch.map{sample_id, reads, paired -> paired}.first()
         
         // Retain only blank samples
         filtered_ch.join(blank_samples).set{reads_ch}

         forward   = reads_ch.map{sample_id, reads, paired -> reads[0]}.toSortedList()
         reverse   = reads_ch.map{sample_id, reads, paired -> reads[1]}.toSortedList()

         
        SPADES(forward, reverse, isPaired, Channel.of("illumina"))
        BUILD_CONTAMINANT_DB(SPADES.out.assembly)
        REMOVE_CONTAMINANT(BUILD_CONTAMINANT_DB.out.index, filtered_ch)

        // Collect software versions
       software_versions_ch = Channel.empty()
       SPADES.out.version | mix(software_versions_ch) | set{software_versions_ch}
       BUILD_CONTAMINANT_DB.out.version | mix(software_versions_ch) | set{software_versions_ch}
       REMOVE_CONTAMINANT.out.version | mix(software_versions_ch) | set{software_versions_ch}

    emit:
       clean_reads = REMOVE_CONTAMINANT.out.reads
       versions = software_versions_ch

}



workflow nano_remove_contaminants {


    take:
        file_ch
        trimmed_ch

    main:

        // Get blank samples names from input file
        file_ch.map{row -> if(deleteWS(row.NTC) == "true") [row.sample_id] }.set{blank_samples}

        // Retain only blank samples
        trimmed_ch.join(blank_samples).set{reads_ch}
        // Get only the blank reads
        blank_reads   = reads_ch.map{sample_id, reads, paired -> reads}.toSortedList()

        // Assemble contaminants / blank reads
        FLYE(blank_reads)
        // Build contaminant index for mapping
        BUILD_CONTAMINANT_INDEX(FLYE.out.assembly)
        // Map all samples reads to contaminant assembly
        MAPPING_TO_CONTAMINANT(BUILD_CONTAMINANT_INDEX.out.index, trimmed_ch)
        
       /*
       - Sort and convert sam to bam
       - Index bam
       - Collect sample mapping stats
       - Remove unmapped reads/contaminants using samtools fastq
       */
       MAPPING_TO_CONTAMINANT.out.sam.map{ sample_id, sam, mapping_info ->
                         tuple([sample_id: sample_id, suffix: "_decontam${params.assay_suffix}", isPaired: 'false'],
                                sam, mapping_info)
                  }.set{mod_sam_ch}

       NANO_REMOVE_CONTAMINANT(mod_sam_ch)

       NANO_REMOVE_CONTAMINANT.out.log.map{ sample_id, stats, logs ->
                                   tuple(sample_id, [stats,logs])
                                   }.set{log_ch}

        // Collect software versions
       software_versions_ch = Channel.empty()
       FLYE.out.version | mix(software_versions_ch) | set{software_versions_ch}
       BUILD_CONTAMINANT_INDEX.out.version | mix(software_versions_ch) | set{software_versions_ch}
       MAPPING_TO_CONTAMINANT.out.version| mix(software_versions_ch) | set{software_versions_ch}
       NANO_REMOVE_CONTAMINANT.out.version | mix(software_versions_ch) | set{software_versions_ch}

    emit:
       clean_reads = NANO_REMOVE_CONTAMINANT.out.reads
       logs        = log_ch
       versions    = software_versions_ch

}
