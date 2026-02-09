#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


// process to count the number / percentage of reads mapping to a custome genome
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
