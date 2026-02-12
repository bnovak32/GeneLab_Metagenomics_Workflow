#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//params.isFast5 = true
//params.input_dir = "path/to/fast5/or/pod5/directory/"
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

// conver5 fast5 to pod5 file
process FAST52POD5 {

    tag "converting fast5 to post5"

    input:
        path(fast5_dir)

    output:
        path("pod5_dir/"), emit: pod5_dir
        path("versions.txt"), emit: version
        
    script:
        """   
        # Convert each fast5 to its relative converted output. The output files are written
        # into the output directory at paths relatve to the path given to the
        # --one-to-one argument. Note: This path must be a relative parent to all
        #input paths.
        #ls input/*.fast5
        #file_1.fast5 file_2.fast5 ... file_N.fast5
        pod5 convert fast5  --output pod5_dir/ --one-to-one ${fast5_dir}


        VERSION=`pod5 --version | awk '{print \$3}'`
        echo pod5 \${VERSION} > versions.txt
        """
}

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



// performs demultiplexing to separate reads based on their barcodes
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
           --kit-name ${kit_name} \\
           ${basecalled}

       VERSION=`dorado --version`
       echo "dorado \${VERSION}" > versions.txt
       """
}



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
 
     if(params.isFast5){

     input_dir   = Channel.fromPath(params.input_dir, checkIfExists: true)
     FAST52POD5(input_dir)

     pod5_dir = FAST52POD5.out.pod5_dir
     }

     DORADO_BASECALLER(pod5_dir, config_file, params.kit_name)

     DORADO_DEMUX(DORADO_BASECALLER.out.bam, params.kit_name)
}
