#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * ========================================================================================
 * PROCESS: KRAKEN_CLASSIFY
 * ========================================================================================
 *
 * SUMMARY:
 *   Classify sample reads with kraken2
 *
 * INPUTS:
 *   1. each: path(DB)
 *      Cardinality: each
 *      Description: Iterates over each DB element. Ensures all samples are processed not just one. 
 *
 *   2. tuple: tuple val(sample_id), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Bolean specifying whether input reads are paired or not 
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}-kraken2-report.tsv") (emit: report)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: Kraken2
 *   Container: [Defined in config/default.config]
 *   Conda: envs/kraken2.yaml
 *   Labels: kraken2
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

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
    REF_DB=`find -L ${DB} -name '*.k2d' |head -n 1|xargs -I {} dirname {}`
    if [ ${isPaired} == 'true' ]; then

        kraken2 --db \${REF_DB} --gzip-compressed \\
            --threads ${task.cpus} \\
            --use-names --paired \\
            --output ${sample_id}-kraken2-output.txt \\
            --report ${sample_id}-kraken2-report.tsv  ${reads[0]} ${reads[1]}

    else

        # Single end
        kraken2 --db \${REF_DB} --gzip-compressed \\
            --threads ${task.cpus} --use-names \\
            --output ${sample_id}-kraken2-output.txt \\
            --report ${sample_id}-kraken2-report.tsv ${reads[0]}

    fi

    VERSION=`echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//'`
    echo "kraken2 \${VERSION}"  > versions.txt
    """
}


/*
 * ========================================================================================
 * PROCESS: KRAKEN2TABLE
 * ========================================================================================
 *
 * SUMMARY:
 *   Create species table from multiple kraken2 reports with pavian
 *
 * INPUTS:
 *   1. path: reports
 *      Cardinality: one
 *      Description: Input file: Kraken2 reports
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}kraken2_species_table${params.assay_suffix}.tsv (emit: table)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: Pavian
 *   Container: [Defined in config/default.config]
 *   Conda: envs/pavian.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process KRAKEN2TABLE {


    tag "Creating a species table from multiple kraken reports.."

    input:
        path(reports)

    output:
        path("${params.additional_filename_prefix}kraken2_species_table${params.assay_suffix}.tsv"), emit: table
        path ("versions.txt"), emit: version

    script:
    """
    merge_kraken_reports.R \\
               --reports-dir '.' \\
               --output-prefix '${params.additional_filename_prefix}' \\
               --assay-suffix '${params.assay_suffix}'

     if [ -f ${params.additional_filename_prefix}merged_kraken_table${params.assay_suffix}.tsv ]; then

          mv ${params.additional_filename_prefix}merged_kraken_table${params.assay_suffix}.tsv \\
             ${params.additional_filename_prefix}kraken2_species_table${params.assay_suffix}.tsv
  
     fi
    
     Rscript -e "VERSIONS=sprintf('pavian %s\\n', packageVersion('pavian')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
    """

}


/*
 * ========================================================================================
 * PROCESS: KAIJU_CLASSIFY
 * ========================================================================================
 *
 * SUMMARY:
 *   Classify sample reads with kaiju
 *
 * INPUTS:
 *   1. each: path(DB)
 *      Cardinality: each
 *      Description: Iterates over each DB element. Ensures all samples are processed not just one. 
 *
 *   2. tuple: tuple val(sample_id), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Bolean specifying whether input reads are paired or not 
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}_kaiju.out") (emit: report)
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


/*
 * ========================================================================================
 * PROCESS: KAIJU2TABLE
 * ========================================================================================
 *
 * SUMMARY:
 *   Merge kaiju reports in a table at the specified taxon level
 *
 * INPUTS:
 *   1. path: DB
 *      Cardinality: one
 *      Description: Input file: Kaiju database directory
 *
 *   2. val: taxon_level
 *      Cardinality: one
 *      Description: Parameter value: taxon_level i.e phylum, class, order, family, genus and species
 *
 *   3. path: reports
 *      Cardinality: one
 *      Description: Input file: kaiju reports
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}merged_kaiju_table.tsv (emit: table)
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

process KAIJU2TABLE {
 
    tag "Merging kaiju reports in a ${taxon_level} table.."
    label "kaiju"

    input:
        path(DB)
        val(taxon_level) // species, genus, order, class, phylum
        path(reports) // _kaiju.out

    output:
        path("${params.additional_filename_prefix}merged_kaiju_table.tsv"), emit: table
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
            -o ${params.additional_filename_prefix}merged_kaiju_table.tsv \\
             \${kaiju_out_FILES[*]}

    # Convert the file names to sample names
    sed -i -E 's/.+\\/(.+)_kaiju\\.out/\\1/g' ${params.additional_filename_prefix}merged_kaiju_table.tsv && \\
    sed -i -E 's/file/sample/' ${params.additional_filename_prefix}merged_kaiju_table.tsv
    VERSION=`echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' )`
    echo "kaiju \${VERSION}"  > versions.txt
    """
}
