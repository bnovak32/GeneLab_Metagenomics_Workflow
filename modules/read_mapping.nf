#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Read mapping to contig assembly using Bowtie2 ********************
****************************************************************************************/

/*
 * ========================================================================================
 * PROCESS: MAPPING
 * ========================================================================================
 *
 * SUMMARY:
 *   Map sample reads to sample assembly with bowtie2
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(assembly), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - assembly: path to sample assembly/contigs
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Bolean specifying whether input reads are paired or not 
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}-mapping-info${params.assay_suffix}.txt") (emit: sam)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: Bowtie2
 *   Container: [Defined in config/illumina.config]
 *   Conda: envs/bowtie2.yaml
 *   Labels: mapping, bowtie2
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// This process builds the bowtie2 index and runs the mapping for each sample
process MAPPING {

    tag "Mapping ${sample_id}-s reads to its assembly ${assembly}..."
    label "mapping"
    label "bowtie2"

    input:
        tuple val(sample_id), path(assembly), path(reads), val(isPaired)
    output:
        tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}-mapping-info${params.assay_suffix}.txt"), emit: sam
        path("versions.txt"), emit: version
    script:
        """
        if [ ${isPaired}  == 'true' ]; then
            # Only running if the assembly produced anything
            if [ -s ${assembly} ]; then

                bowtie2-build ${assembly} ${sample_id}-index 
                bowtie2 --mm -q --threads ${task.cpus} \\
                        -x ${sample_id}-index  -1 ${reads[0]} -2 ${reads[1]} \\
                        --no-unal > ${sample_id}.sam  2> ${sample_id}-mapping-info${params.assay_suffix}.txt 
            rm ${sample_id}-index*
            else

                touch ${sample_id}.sam
                echo "Mapping not performed for ${sample_id} because the assembly didn't produce anything." > ${sample_id}-mapping-info${params.assay_suffix}.txt
                printf "Mapping not performed for ${sample_id} because the assembly didn't produce anything.\\n"

            fi
        # Single-end
        else

            # Only running if the assembly produced anything
            if [ -s ${assembly} ]; then

                bowtie2-build ${assembly} ${sample_id}-index 
                bowtie2 --mm -q --threads ${task.cpus} \\
                        -x ${sample_id}-index -r ${reads[0]} \\
                        --no-unal > ${sample_id}.sam  2> ${sample_id}-mapping-info${params.assay_suffix}.txt
                        
                rm ${sample_id}-index*
            else

                touch ${sample_id}.sam
                echo "Mapping not performed for ${sample_id} because the assembly didn't produce anything."  > ${sample_id}-mapping-info${params.assay_suffix}.txt
                printf "Mapping not performed for ${sample_id} because the assembly didn't produce anything.\\n"

            fi

        fi
        bowtie2 --version  | head -n 1 | sed -E 's/.*(bowtie2-align-s version.+)/\\1/' > versions.txt
        """
}


/*
 * ========================================================================================
 * PROCESS: LONG_MAPPING
 * ========================================================================================
 *
 * SUMMARY:
 *   Map sample reads to sample assembly with minimap2
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(assembly), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - assembly: path to sample assembly/contigs
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Bolean specifying whether input reads are paired or not  
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}-mapping-info${params.assay_suffix}.txt") (emit: sam)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: Minimap2
 *   Container: [Defined in config/nanopore.config]
 *   Conda: envs/minimap2.yaml
 *   Labels: minimap2, mapping
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process LONG_MAPPING {

    tag "Mapping ${sample_id}-s reads to its assembly ${assembly}..."
    label "minimap2"
    label "mapping"

    input:
        tuple val(sample_id), path(assembly), path(reads), val(isPaired)
    output:
        tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}-mapping-info${params.assay_suffix}.txt"), emit: sam
        path("versions.txt"), emit: version
    script:
        """
            # Only running if the assembly produced anything
            if [ -s ${assembly} ]; then

               minimap2 -ax map-ont -t ${task.cpus} ${assembly} ${reads[0]} \\
                     > ${sample_id}.sam  2> ${sample_id}-mapping-info${params.assay_suffix}.txt
                        
     
            else

                touch ${sample_id}.sam
                echo "Mapping not performed for ${sample_id} because the assembly didn't produce anything."  > ${sample_id}-mapping-info${params.assay_suffix}.txt
                printf "Mapping not performed for ${sample_id} because the assembly didn't produce anything.\\n"

            fi

        VERSION=`minimap2 --version`
        echo "minimap2 \${VERSION}" > versions.txt
        """
}


/*
 * ========================================================================================
 * PROCESS: SAM_TO_BAM
 * ========================================================================================
 *
 * SUMMARY:
 *   Sort and convert sample sam to bam files
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(sam), path(mapping_info)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - sam: path to sample sam file
 *                 - mapping_info: path to sample mapping info
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}${params.assay_suffix}.bam") (emit: bam)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: SAMtools
 *   Container: [Defined in config/default.config]
 *   Conda: envs/samtools.yaml
 *   Labels: samtools
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process SAM_TO_BAM {

    tag "Sorting and converting ${sample_id}-s sam to bam files..."
    label "samtools"

    input:
        tuple val(sample_id), path(sam), path(mapping_info)
    output:
        tuple val(sample_id), path("${sample_id}${params.assay_suffix}.bam"), emit: bam
        path("versions.txt"), emit: version
    script:
        """
        # Only running if the assembly produced anything
        if [ -s ${sam} ]; then

            samtools sort -@ ${task.cpus} ${sam} > ${sample_id}${params.assay_suffix}.bam 2> /dev/null

        else

            touch ${sample_id}${params.assay_suffix}.bam
            printf "Sorting and converting not performed for ${sample_id} because read mapping didn't produce anything.\\n"

        fi
        samtools --version | head -n1 > versions.txt
        """
}
