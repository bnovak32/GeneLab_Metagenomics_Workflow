#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * ========================================================================================
 * PROCESS: SAMTOOLS_SORT
 * ========================================================================================
 *
 * SUMMARY:
 *   Sort and convert sample sam file to bam file
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(sam), path(mapping_info)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - sam: path to sample sam file
 *                 - mapping_info: path to sample mapping 
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}.bam") (emit: bam)
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

process SAMTOOLS_SORT {

    tag "Sorting and converting ${sample_id}-s sam to bam files..."
    label "samtools"
    

    input:
        tuple val(sample_id), path(sam), path(mapping_info)
    output:
        tuple val(sample_id), path("${sample_id}.bam"), emit: bam
        path("versions.txt"), emit: version
    script:
        """
        samtools sort -@ ${task.cpus} ${sam} > ${sample_id}.bam 2> /dev/null

        samtools --version | head -n1 > versions.txt
        """

}

/*
 * ========================================================================================
 * PROCESS: SAMTOOLS_INDEX
 * ========================================================================================
 *
 * SUMMARY:
 *   Index sample bam file
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(bam)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - bam: path to sample bam file
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}.bam.bai") (emit: bai)
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

process SAMTOOLS_INDEX {

    tag "Indexing ${sample_id}-s bam file.."
    label "samtools"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}.bam.bai"), emit: bai
        path("versions.txt"), emit: version

    script:
    """
    samtools index ${bam} ${sample_id}.bam.bai  \\
            > ${sample_id}.log 2>&1

    samtools --version | head -n1 > versions.txt
    """
}


/*
 * ========================================================================================
 * PROCESS: SAMTOOLS_STATS
 * ========================================================================================
 *
 * SUMMARY:
 *   Collect sample mapping to reference statistics
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(bam), path(bai)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - bam: path to sample bam file
 *                 - bai: path to sample bam index file 
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("*.txt"), path("*.log") (emit: log)
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

process SAMTOOLS_STATS {

    tag "collecting ${sample_id}-s mapping stats..."
    label "samtools"

    input:
        tuple val(sample_id), path(bam), path(bai)

    output:
        tuple val(sample_id), path("*.txt"), path("*.log"), emit: log
        path("versions.txt"), emit: version
        
    script:
    """
    # Collect mapping stats to be used by multiqc to generate summary reports
    samtools flagstat ${bam} > \\
             ${sample_id}_flagstats.txt  \\
             2> ${sample_id}_flagstats.log

    samtools stats --remove-dups ${bam} \\
             > ${sample_id}_stats.txt   \\
             2> ${sample_id}_stats.log

    samtools idxstats ${bam}  \\
             > ${sample_id}_idxstats.txt \\
             2> ${sample_id}_idxstats.log

    samtools --version | head -n1 > versions.txt
    """
}


/*
 * ========================================================================================
 * PROCESS: SAMTOOLS_FASTQ
 * ========================================================================================
 *
 * SUMMARY:
 *   Retrieve unmapped sample reads with samtools fastq
 *
 * INPUTS:
 *   1. tuple: tuple val(meta), path(bam), path(bai)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - meta: map with the following structure [sample_id:sample_name, suffix:suffix_name, isPaired:paired]
 *                 - bam: path to sample bam file
 *                 - bai: path to sample bam index file  
 *
 * OUTPUTS:
 *   1. tuple: tuple val(meta.sample_id), path("*.fastq.gz"), val(meta.isPaired) (emit: reads)
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

// Filter Unmapped reads
process SAMTOOLS_FASTQ {
   
    tag "Retrieving ${meta.sample_id}-s unmapped reads..."
    label "samtools"

    input:
        tuple val(meta), path(bam), path(bai) // [[sample_id:sample_name, suffix:suffix_name, isPaired:paired], bam, bai]

    output:
        tuple val(meta.sample_id), path("*.fastq.gz"), val(meta.isPaired), emit: reads
        path("versions.txt"), emit: version
    script:
    """
    if [ ${meta.isPaired} == 'true' ]; then

        samtools fastq \\
              -f12 -F256 \\
              -1 ${meta.sample_id}${meta.suffix}_R1.fastq.gz \\
              -2 ${meta.sample_id}${meta.suffix}_R2.fastq.gz ${bam} 

    else

        samtools fastq -t -f 4  ${bam} \\
              > ${meta.sample_id}${meta.suffix}.fastq  && \\
        gzip ${meta.sample_id}${meta.suffix}.fastq

    fi
    
    samtools --version | head -n1 > versions.txt    
    """
}

/*
 * ========================================================================================
 * PROCESS: SAMTOOLS_FILTER_FASTQ
 * ========================================================================================
 *
 * SUMMARY:
 *   Retrieve sample unmapped reads with samtools fastq. Includes prior steps of sorting, indexing and alignment stats.
 *
 * INPUTS:
 *   1. tuple: tuple val(meta), path(sam), path(mapping_info)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - meta: map with the following structure [sample_id:sample_name, suffix:suffix_name, isPaired:paired]
 *                 - sam: path to sample sam file
 *                 - mapping_info: path to sample mapping info
 *
 * OUTPUTS:
 *   1. tuple: tuple val(meta.sample_id), path("${meta.sample_id}${meta.suffix}.fastq.gz"), val(meta.isPaired) (emit: reads)
 *
 *   2. tuple: tuple val(meta.sample_id), path("${meta.sample_id}.bam"),
 *
 *   3. path: ${meta.sample_id}.bam.bai (emit: bam)
 *
 *   4. tuple: tuple val(meta.sample_id), path("${meta.sample_id}_stats.txt"), path("*.log") (emit: log)
 *
 *   5. path: versions.txt (emit: version)
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

process SAMTOOLS_FILTER_FASTQ {

    tag "Retrieving ${meta.sample_id}-s unmapped reads..."
    label "samtools"

    input:
        tuple val(meta), path(sam), path(mapping_info) // [[sample_id:sample_name, suffix:suffix_name, isPaired:paired], sam, mapping_info]
    output:
        tuple val(meta.sample_id), path("${meta.sample_id}${meta.suffix}.fastq.gz"), val(meta.isPaired), emit: reads
        tuple val(meta.sample_id), path("${meta.sample_id}.bam"), 
              path("${meta.sample_id}.bam.bai"), emit: bam
        tuple val(meta.sample_id), path("${meta.sample_id}_stats.txt"), path("*.log"), emit: log
        path("versions.txt"), emit: version
    script:
        """
        # Sam to bam
        samtools sort -@ ${task.cpus} ${sam} > ${meta.sample_id}.bam 2> /dev/null
        # Index bam
        samtools index ${meta.sample_id}.bam ${meta.sample_id}.bam.bai  \\
                > ${meta.sample_id}.log 2>&1

        # Collect mapping stats to be used by multiqc to generate summary reports
        samtools flagstat ${meta.sample_id}.bam \\
             > ${meta.sample_id}_flagstats.txt  \\
             2> ${meta.sample_id}_flagstats.log

        samtools stats --remove-dups ${meta.sample_id}.bam \\
             > ${meta.sample_id}_stats.txt   \\
             2> ${meta.sample_id}_stats.log

        samtools idxstats ${meta.sample_id}.bam \\
             > ${meta.sample_id}_idxstats.txt \\
             2> ${meta.sample_id}_idxstats.log

        if [ ${meta.isPaired} == 'true' ]; then

            samtools fastq \\
              -f12 -F256 \\
              -1 ${meta.sample_id}${meta.suffix}_R1.fastq.gz \\
              -2 ${meta.sample_id}${meta.suffix}_R2.fastq.gz ${meta.sample_id}.bam 

        else

            samtools fastq -t -f 4  ${meta.sample_id}.bam \\
              > ${meta.sample_id}${meta.suffix}.fastq  && \\
            gzip ${meta.sample_id}${meta.suffix}.fastq

        fi
    
        samtools --version | head -n1 > versions.txt
        """

}

