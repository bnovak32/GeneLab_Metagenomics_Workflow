#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * ========================================================================================
 * PROCESS: LONG_MAP2GENOME
 * ========================================================================================
 *
 * SUMMARY:
 *   Map sample reads to a custom genome reference
 *
 * INPUTS:
 *   1. each: path(REF)
 *      Cardinality: each
 *      Description: Iterates over each element. Custom reference for genome mapping.
 *
 *   2. each: prefix
 *      Cardinality: each
 *      Description: Iterates over each element. Prefix to add to output file.
 *
 *   3. tuple: tuple val(sample_id), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Boolean specifying whether input reads are paired or not 
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id) , path("${prefix}-${sample_id}_scaffold_stats.txt") (emit: stats)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bbmap.yaml
 *   Labels: bbtools
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// Process to count the number / percentage of reads mapping to a custome genome
process LONG_MAP2GENOME {

    tag "mapping ${sample_id}-s reads to a custome reference"
    label "bbtools"

    input:
        each path(REF)
        each prefix
        tuple val(sample_id), path(reads), val(isPaired)

    output:
        tuple val(sample_id) , path("${prefix}-${sample_id}_scaffold_stats.txt"), emit: stats
        path("versions.txt"), emit: version

    script:
    def maxmem = task.memory.toGiga()
    """
    mapPacBio.sh -Xmx${maxmem}g in=${reads[0]}  \\
              ref=${REF} \\
              ambiguous=all \\
              ignorebadquality=t \\
              k=8 \\
              scafstats=${prefix}-${sample_id}_scaffold_stats.txt

    VERSION=`bbversion.sh`
    echo "bbtools \${VERSION}" > versions.txt
    """
}

/*
 * ========================================================================================
 * PROCESS: SHORT_MAP2GENOME
 * ========================================================================================
 *
 * SUMMARY:
 *   Map sample reads to a custom reference
 *
 * INPUTS:
 *   1. each: path(REF)
 *      Cardinality: each
 *      Description: Iterates over each element. Custom reference for genome mapping.
 *
 *   2. each: prefix
 *      Cardinality: each
 *      Description: Iterates over each element. refix to add to output file.
 *
 *   3. tuple: tuple val(sample_id), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Boolean specifying whether input reads are paired or not 
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id) , path("${prefix}-${sample_id}_scaffold_stats.txt") (emit: stats)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bbmap.yaml
 *   Labels: bbtools
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process SHORT_MAP2GENOME {

    tag "mapping ${sample_id}-s reads to a custome reference"
    label "bbtools"

    input:
        each path(REF)
        each prefix
        tuple val(sample_id), path(reads), val(isPaired)

    output:
        tuple val(sample_id) , path("${prefix}-${sample_id}_scaffold_stats.txt"), emit: stats
        path("versions.txt"), emit: version

    script:
    def maxmem = task.memory.toGiga()
    def input = isPaired == 'true' ? "in1=${reads[0]} in2=${reads[1]}" : "in=${reads[0]}" 
    """

    bbmap.sh -Xmx${maxmem}g ${input}  \\
              ref=${REF} \\
              ambiguous=all \\
              ignorebadquality=t \\
              k=8 \\
              scafstats=${prefix}-${sample_id}_scaffold_stats.txt

    VERSION=`bbversion.sh`
    echo "bbtools \${VERSION}" > versions.txt
    """
}
