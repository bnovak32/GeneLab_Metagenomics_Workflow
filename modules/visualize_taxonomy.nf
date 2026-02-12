#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

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



process KRONA_REPORT {

    tag "Creating a krona html report.."
    label "read_based_outputs"
    label "krona"


    input:
        val(prefix) // kaiju, kraken, metaphlan etc
        path(krona_files)

    output:
        path("${prefix}-report.html"), emit: html
        path("versions.txt"), emit: version

    script:
    """
    find -L . -type f -name "*.krona" |sort -uV > krona_files.txt

    FILES=(\$(find -L . -type f -name "*.krona"))
    basename -a -s '.krona' \${FILES[*]} | sort -uV  > sample_names.txt
    KTEXT_FILES=(\$(paste -d',' "krona_files.txt" "sample_names.txt"))
    
    ktImportText  -o ${prefix}-report.html \${KTEXT_FILES[*]}

    VERSION=`echo \$(ktImportText 2>&1) | sed 's/^.*KronaTools //g; s/- ktImportText.*\$//g'`
    echo "krona \${VERSION}" > versions.txt
    """
}
