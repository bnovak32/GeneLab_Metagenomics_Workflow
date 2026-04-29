#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * ========================================================================================
 * PROCESS: KRAKEN2KRONA
 * ========================================================================================
 *
 * SUMMARY:
 *   Convert kraken2 report file to krona file
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(report)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - report: path to sample kraken2 report
 *
 * OUTPUTS:
 *   1. path: ${sample_id}.krona (emit: krona)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: Kraken2
 *   Container: [Defined in config/default.config]
 *   Conda: envs/kraken2.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process KRAKEN2KRONA {

    tag "Converting kraken file to krona file .." 

    input:
        tuple val(sample_id), path(report)

    output:
        path("${sample_id}.krona"), emit: krona
        path("versions.txt"), emit: version

    script:
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    kreport2krona.py \\
          --report-file ${report}  \\
          --output ${sample_id}.krona

    echo "krakentools ${VERSION}"  > versions.txt
    """
}

/*
 * ========================================================================================
 * PROCESS: METAPHLAN2KRONA
 * ========================================================================================
 *
 * SUMMARY:
 *   Convert metaphlan file to krona file
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(report)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - report: path to sample metaphlan taxonomy report 
 *
 * OUTPUTS:
 *   1. path: ${sample_id}.krona (emit: krona)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: Python
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: bit (Python environment)
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process METAPHLAN2KRONA {

    tag "Converting metaphlan file to krona file .."
    label "bit" //python environment

    input:
        tuple val(sample_id), path(report)

    output:
        path("${sample_id}.krona"), emit: krona
        path("versions.txt"), emit: version
    script:
    """
    metaphlan2krona.py \\
            --profile ${report} \\
            --krona ${sample_id}.krona
    python --version >  versions.txt
    """

}

/*
 * ========================================================================================
 * PROCESS: KAIJU2KRONA
 * ========================================================================================
 *
 * SUMMARY:
 *   Convert kaiju file to krona file
 *
 * INPUTS:
 *   1. each: path(DB)
 *      Cardinality: each
 *      Description: Iterates over each element. Kaiju database directory.
 *
 *   2. tuple: tuple val(sample_id), path(report)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - report: path to sample kaiju report  
 *
 * OUTPUTS:
 *   1. path: ${sample_id}.krona (emit: krona)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/kaiju.yaml
 *   Labels: kaiju
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process KAIJU2KRONA {

    tag "Converting kaiju file to krona file .."
    label "kaiju"


    input:
        each path(DB)
        tuple val(sample_id), path(report)

    output:
        path("${sample_id}.krona"), emit: krona
        path("versions.txt"), emit: version

    script:
    """
    NODES=`find -L ${DB} -name "*nodes.dmp"`
    NAMES=`find -L ${DB} -name "*names.dmp"`

    kaiju2krona -u \\
        -n \${NAMES} \\
        -t \${NODES} \\
	-i ${report} \\
	-o ${sample_id}.krona

    VERSION=`kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //'`

    echo "kaiju \${VERSION}" > versions.txt
    """
}


/*
 * ========================================================================================
 * PROCESS: KRONA_REPORT
 * ========================================================================================
 *
 * SUMMARY:
 *   Create a krona html report..
 *
 * INPUTS:
 *   1. val: prefix
 *      Cardinality: one
 *      Description: Parameter value: prefix for read classifier i.e. kaiju, kraken2 and metaphlan
 *
 *   2. path: krona_files
 *      Cardinality: one
 *      Description: Input file: krona files
 *
 * OUTPUTS:
 *   1. path: ${prefix}-report${params.assay_suffix}.html (emit: html)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/krona.yaml
 *   Labels: krona
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process KRONA_REPORT {

    tag "Creating a krona html report.."
    label "krona"


    input:
        val(prefix) // kaiju, kraken, metaphlan etc
        path(krona_files)

    output:
        path("${prefix}-report${params.assay_suffix}.html"), emit: html
        path("versions.txt"), emit: version

    script:
    """
    find -L . -type f -name "*.krona" |sort -uV > krona_files.txt

    FILES=(\$(find -L . -type f -name "*.krona"))
    basename -a -s '.krona' \${FILES[*]} | sort -uV  > sample_names.txt
    KTEXT_FILES=(\$(paste -d',' "krona_files.txt" "sample_names.txt"))
    
    ktImportText  -o ${prefix}-report${params.assay_suffix}.html \${KTEXT_FILES[*]}

    VERSION=`echo \$(ktImportText 2>&1) | sed 's/^.*KronaTools //g; s/- ktImportText.*\$//g'`
    echo "krona \${VERSION}" > versions.txt
    """
}
