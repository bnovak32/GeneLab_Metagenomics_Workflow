#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


// Classify reads using kraken2
process KRAKEN_CLASSIFY {

    tag "Classifying ${sample_id}-s reads with kraken..."
    label "kraken2"

    input:
        each path(DB)
        tuple val(sample_id), path(reads), val(isPaired)

    output:
        tuple val(sample_id), path("${sample_id}-kraken2-report.tsv"), emit: report
        path("versions.txt"), emit: version

    script:
    """
    if [ ${isPaired} == 'true' ]; then

        kraken2 --db ${DB} --gzip-compressed \\
            --threads ${task.cpus} \\
            --use-names --paired \\
            --output ${sample_id}-kraken2-output.txt \\
            --report ${sample_id}-kraken2-report.tsv  ${reads[0]} ${reads[1]}

    else

        # Single end
        kraken2 --db ${DB} --gzip-compressed \\
            --threads ${task.cpus} --use-names \\
            --output ${sample_id}-kraken2-output.txt \\
            --report ${sample_id}-kraken2-report.tsv ${reads[0]}

    fi

    VERSION=`echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//'`
    echo "kraken2 \${VERSION}"  > versions.txt
    """
}


process KRAKEN2TABLE {


    tag "Creating a species table from multiple kraken reports.."
    label "read_based_outputs"

    input:
        path(reports)

    output:
        path("${params.additional_filename_prefix}merged_kraken_table${params.assay_suffix}.tsv"), emit: table
        path ("versions.txt"), emit: version

    script:
    """
    merge_kraken_reports.R \\
               --reports-dir '.' \\
               --output-prefix '${params.additional_filename_prefix}' \\
               --assay-suffix '${params.assay_suffix}'

     Rscript -e "VERSIONS=sprintf('pavian %s\\n', packageVersion('pavian')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
    """

}



process KAIJU_CLASSIFY {

    tag "Classifying ${sample_id}-s reads with kaiju..."
    label "kaiju"


    input:
        each path(DB)
        tuple val(sample_id), path(reads), val(isPaired)
    
    output:
        tuple val(sample_id), path("${sample_id}_kaiju.out"), emit: report
        path ("versions.txt"), emit: version

    script:
        def input = isPaired ? "-i ${reads[0]}" : "-i ${reads[0]} -j ${reads[1]}"
    """
    NODES=`find -L . ${DB} -name "*nodes.dmp"`
    FMI=`find -L . ${DB} -name "*.fmi" -not -name "._*"`

    kaiju \\
        -f \${FMI} \\
        -t \${NODES} \\
        -z 10 \\
        -E 1e-05 \\
        -o ${sample_id}_kaiju.out \\
        ${input}

    VERSION=`kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //'`
    echo "kaiju \${VERSION}" > versions.txt    
    """
}


process KAIJU2TABLE {
 
    tag "Merging kaiju reports in a ${taxon_level} table.."
    label "kaiju"
    label "read_based_outputs"

    input:
        path(DB)
        val(taxon_level) // species, genus, order, class, phylum
        path(reports) // _kaiju.out

    output:
        path("${params.additional_filename_prefix}merged_kaiju_table${params.assay_suffix}.tsv"), emit: table
        path("versions.txt"), emit: version

    script:
    """
    NODES=`find -L . ${DB} -name "*nodes.dmp"`
    NAMES=`find -L . ${DB} -name "*names.dmp"`
    kaiju_out_FILES=(\$(find -L . -type f -name "*_kaiju.out")) 
    
    kaiju2table \\
            -t \${NODES} \\
            -n \${NAMES} \\
            -p  \\
            -r ${taxon_level} \\
            -o ${params.additional_filename_prefix}merged_kaiju_table${params.assay_suffix}.tsv \\
             \${kaiju_out_FILES[*]}

    # Convert the file names to sample names
    sed -i -E 's/.+\\/(.+)_kaiju\\.out/\\1/g' ${params.additional_filename_prefix}merged_kaiju_table${params.assay_suffix}.tsv && \\
    sed -i -E 's/file/sample/' ${params.additional_filename_prefix}merged_kaiju_table${params.assay_suffix}.tsv
    VERSION=`echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' )`
    echo "kaiju \${VERSION}"  > versions.txt
    """
}
