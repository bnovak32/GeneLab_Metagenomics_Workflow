#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//params.input_dir = "path/to/pod5/directory/"
// can be path sof existing config of just the model category i.e. hac, sup and fast etc
//params.model = "/path/to/dorado/mmodel/file/dna_r10.4.1_e8.2_400bps_hac@v4.3.0"
//params.kit_name = "SQK-RBK114-24" // nanopore kit used.options are:
// EXP-NBD103 EXP-NBD104 EXP-NBD114 EXP-N             
// BD114-24 EXP-NBD196 EXP-PBC001 EXP-PBC096 SQK-16S024 SQK-16S114-24 
// SQK-LWB001 SQK-MAB114-24 SQK-MLK111-96-XL SQK-MLK114-96-XL              
// SQK-NBD111-24 SQK-NBD111-96 SQK-NBD114-24 SQK-NBD114-96 SQK-PBK004 
// SQK-PCB109 SQK-PCB110 SQK-PCB111-24 SQK-PCB114-24 SQK-RAB201 
// SQK-RAB204 SQK-RBK001 SQK-RBK004 SQK-RBK110-96 SQK-RBK111-24 
// SQK-RBK111-96 SQK-RBK114-24 SQK-RBK114-96 SQK-RLB001 SQK-RPB004 
// SQK-RPB114-24 TWIST-16-UDI TWIST-96A-UDI VSK-PTC001 VSK-VMK001 VSK-VMK004 VSK-VPS001


/*
 * ========================================================================================
 * PROCESS: DORADO_BASECALLER
 * ========================================================================================
 *
 * SUMMARY:
 *   Bascall pod5 files with dorado
 *
 * INPUTS:
 *   1. path: input_dir
 *      Cardinality: one
 *      Description: Input file: input directory contaning pod5 files to basecall
 *
 *   2. val: kit_name
 *      Cardinality: one
 *      Description: Parameter value:  Name of the nanopore sequencing kit used
 *
 * OUTPUTS:
 *   1. path: basecalled.bam (emit: bam)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/nanopore.config]
 *   Conda: envs/dorado.yaml
 *   Labels: dorado
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

//  Runs the basecalling process to convert raw nanopore signal data to nucleotide sequences
// https://github.com/nanoporetech/dorado/
process DORADO_BASECALLER {


    tag "Bascalling pod5 files using dorado"
    label "dorado"

    input:
        path(input_dir) // pod5 didrectory
        val(kit_name)

    output:
        path("basecalled.bam"), emit: bam
        path("versions.txt"), emit: version


    script:
       """
       dorado basecaller hac ${input_dir} \\
             --recursive \\
             --kit-name ${kit_name} \\
             --min-qscore 8 > basecalled.bam

       VERSION=`dorado --version`
       echo "dorado \${VERSION}" > versions.txt
       """
}


/*
 * ========================================================================================
 * PROCESS: DORADO_DEMUX
 * ========================================================================================
 *
 * SUMMARY:
 *   Demultiplex basecalled bam
 *
 * INPUTS:
 *   1. path: basecalled
 *      Cardinality: one
 *      Description: Input file: basecalled bam file
 *
 *   2. val: kit_name
 *      Cardinality: one
 *      Description: Parameter value: Name of the nanopore sequencing kit used
 *
 * OUTPUTS:
 *   1. path: demultiplexed/ (emit: demux_dir)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/nanopore.config]
 *   Conda: envs/dorado.yaml
 *   Labels: dorado
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// Performs demultiplexing to separate reads based on their barcodes
process DORADO_DEMUX {


    tag "Demultiplexing basecalled bam"
    label "dorado"

    input:
        path(basecalled)
        val(kit_name)

    output:
        path("demultiplexed/"), emit: demux_dir
        path("versions.txt"), emit: version


    script:
       """
       dorado demux \\
           --output-dir demultiplexed/ \\
           --emit-fastq \\
           --emit-summary \\
           --kit-name ${kit_name} \\
           ${basecalled}

       VERSION=`dorado --version`
       echo "dorado \${VERSION}" > versions.txt
       """
}

/*
 * ========================================================================================
 * PROCESS: CAT_FASTQ_FILES
 * ========================================================================================
 *
 * SUMMARY:
 *   Concatenate sample fastq files
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(reads)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - reads: path to a list sample fastq reads to concatenate
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}.fastq.gz") (emit: reads)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: bit
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process CAT_FASTQ_FILES {

    tag "Concatenating barcode split fastqs..."
    label "bit"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.fastq.gz"), emit: reads
        path("versions.txt"), emit: version

    script:
    def readList = reads instanceof List ? reads.collect { it.toString() } : [reads.toString()]
    def command = readList[0].endsWith('.gz') ? 'zcat' : 'cat'
    """

    ${command}  ${readList.join(' ')} > ${sample_id}.fastq 

    [ -e ${sample_id}.fastq.gz ] && rm -rf ${sample_id}.fastq.gz
    gzip ${sample_id}.fastq

    VERSION=\$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    echo "cat \${VERSION}" > versions.txt
    """
}

/*
 * ========================================================================================
 * PROCESS: CAT_FASTQ_DIR
 * ========================================================================================
 *
 * SUMMARY:
 *   Concatenating barcode split fastqs...
 *
 * INPUTS:
 *   1. path: sample_to_barcode_file
 *      Cardinality: one
 *      Description: Input file: CSV file mapping sample name to barcode name
 *
 *   2. path: demux_dir
 *      Cardinality: one
 *      Description: Input file: directory containing fastq files to be concatenated
 *
 * OUTPUTS:
 *   1. path: *.fastq.gz (emit: reads)
 *
 *   2. path: runsheet.csv (emit: runsheet)
 *
 *   3. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: bit
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process CAT_FASTQ_DIR {

    tag "Concatenating barcode split fastqs..."
    label "bit"

    input:
        path(sample_to_barcode_file) // a 2-column comma separated file mapping barcode to sample id
        path(demux_dir) // demultiplexed/

    output:
        path("*.fastq.gz"), emit: reads
        path("runsheet.csv"), emit: runsheet
        path("versions.txt"), emit: version

    script:
    """
    # Keeping track of working directory
    WORK_DIR=`pwd`
    # Initiate runsheet header
    echo "sample_id,forward" > \${WORK_DIR}/runsheet.csv

    # Initiate sample name to barcode associative array
    declare -A SAMPLE_TO_BARCODE

    # Read file line-by-line into associative array
    # Handles spaces in values safely
    while IFS=',' read -r sample barcode; do
         # Skip empty lines or lines without both columns
        [[ -z "\${sample}" || -z "\${barcode}" ]] && continue
        SAMPLE_TO_BARCODE["\${sample}"]="\${barcode}"
    done < ${sample_to_barcode_file}

    # Change to directory containing split fastq files generated
    # from the split fastq process above
    cd ${demux_dir} # demultiplexed/

    # Concat separate barcode/sample fastq files into per sample fastq gzippped files
    for sample in "\${!SAMPLE_TO_BARCODE[@]}"; do


    [ -d  \${WORK_DIR}/\${sample}/ ] ||  mkdir -p \${WORK_DIR}/\${sample}/  
    cp *_\${SAMPLE_TO_BARCODE[\$sample]}*  \${WORK_DIR}/\${sample}/ 

    
    if ls -1 \${WORK_DIR}/\${sample}/| head -n1| xargs -I {} gzip -t \${WORK_DIR}/\${sample}/{}  2>/dev/null; then 
   
        # If files are gzipped
        zcat \${WORK_DIR}/\${sample}/* > \${WORK_DIR}/\${sample}.fastq 

    else

        # If files are not gzipped
        cat \${WORK_DIR}/\${sample}/* > \${WORK_DIR}/\${sample}.fastq

    fi

    gzip \${WORK_DIR}/\${sample}.fastq && \\
    rm -rf \${WORK_DIR}/\${sample}/

    # Create runsheet
    echo "\${sample},\${WORK_DIR}/\${sample}.fastq.gz" >> \${WORK_DIR}/runsheet.csv

    done
    cd \${WORK_DIR}/

    VERSION=\$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    echo "cat \${VERSION}" > versions.txt
    """
}

workflow{


     pod5_dir   = Channel.fromPath(params.input_dir, checkIfExists: true)
     config_file = Channel.fromPath(params.config_file, checkIfExists: true)
 
     DORADO_BASECALLER(pod5_dir, config_file, params.kit_name)

     DORADO_DEMUX(DORADO_BASECALLER.out.bam, params.kit_name)
}
