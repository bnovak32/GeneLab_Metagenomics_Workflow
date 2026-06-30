#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * ========================================================================================
 * PROCESS: METAPHLAN2COUNT
 * ========================================================================================
 *
 * SUMMARY:
 *   Convert metaphlan's default relative abundance output to a raw counts species table
 *
 * INPUTS:
 *   1. path: metaphlan_table
 *      Cardinality: one
 *      Description: Input file: default metaphlan table
 *
 *   2. path: reads_per_sample
 *      Cardinality: one
 *      Description: Input file: reads per sample file
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}metaphlan_species_table${params.assay_suffix}.tsv (emit: table)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/R_visualizations.yaml
 *   Labels: R_downstream
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process METAPHLAN2COUNT {

    tag "Processing metaphlan count table..."
    label "R_downstream"

    input:
        path(metaphlan_table) // metaphlan-taxonomy_GLmetagenomics.tsv
        path(reads_per_sample) // reads_per_sample.txt

    output:
        path("${params.additional_filename_prefix}metaphlan_species_table${params.assay_suffix}.tsv"), emit: table
        path("versions.txt"), emit: version
    script:
        """
        process_metaphlan.R \\
                  --metaphlan-table '${metaphlan_table}' \\
                  --read-count '${reads_per_sample}' \\
                  --output-prefix '${params.additional_filename_prefix}' \\
                  --assay-suffix '${params.assay_suffix}'

        Rscript -e "VERSIONS=sprintf('tibble %s\\ntidyr %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\nmagrittr %s\\nglue %s\\n', \\
                                   packageVersion('tibble'), \\
                                   packageVersion('tidyr'), \\
                                   packageVersion('dplyr'), \\
                                   packageVersion('purrr'), \\
                                   packageVersion('readr'), \\
                                   packageVersion('stringr'), \\
                                   packageVersion('magrittr'), \\
                                   packageVersion('glue')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}


/*
 * ========================================================================================
 * PROCESS: KAIJU2SPECIES_TABLE
 * ========================================================================================
 *
 * SUMMARY:
 *   Generate kaiju species count table
 *
 * INPUTS:
 *   1. path: merged_table
 *      Cardinality: one
 *      Description: Input file: merged kraken2 table generated from kraken2 reports
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}kaiju_species_table${params.assay_suffix}.tsv (emit: table)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/R_visualizations.yaml
 *   Labels: R_downstream
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process KAIJU2SPECIES_TABLE  { 

    tag "Processing kaiju species table..."
    label "R_downstream"

    input:
        path(merged_table) // merged_kaiju_table.tsv

    output:

       path("${params.additional_filename_prefix}kaiju_species_table${params.assay_suffix}.tsv"), emit: table
       path("versions.txt"), emit: version

    script:
        """
        process_kaiju_table.R \\
              --merged-table '${merged_table}' \\
              --output-prefix '${params.additional_filename_prefix}' \\
              --assay-suffix '${params.assay_suffix}'


        Rscript -e "VERSIONS=sprintf('tibble %s\\ntidyr %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\nmagrittr %s\\nglue %s\\n', \\
                                   packageVersion('tibble'), \\
                                   packageVersion('tidyr'), \\
                                   packageVersion('dplyr'), \\
                                   packageVersion('purrr'), \\
                                   packageVersion('readr'), \\
                                   packageVersion('stringr'), \\
                                   packageVersion('magrittr'), \\
                                   packageVersion('glue')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}

/*
 * ========================================================================================
 * PROCESS: FILTER_RARE
 * ========================================================================================
 *
 * SUMMARY:
 *   Filter out rare features
 *
 * INPUTS:
 *   1. val: meta
 *      Cardinality: one
 *      Description: Parameter value: a mapping with the following keys:
 *                  mode: filtering mode. 'across_samples' and 'values_sum' for read and assembly based analysis, respectively.  
 *                  filter_threshold : threshold for filtering out rare taxa
 *                  output_file: Output tsv file name
 *
 *   2. path: feature_table
 *      Cardinality: one
 *      Description: Input file: feature table to filter
 *
 * OUTPUTS:
 *   1. path: ${meta.output_file} (emit: table)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/R_visualizations.yaml
 *   Labels: R_downstream
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process FILTER_RARE {

    tag "Filtering out rare features..."
    label "R_downstream"


    input:
       val(meta) // [mode: 'max_value' //['max_value', 'across_samples', 'values_sum'],
                 //  filter_threshold : 0.1,
                 //  output_file: 'kaiju_filtered_species_table_GLlbnMetag.tsv' ] 
        path(feature_table)  // 'kaiju_species_table_GLlbnMetag.tsv'

    output:
       path("${meta.output_file}"), emit: table
       path("versions.txt"), emit: version

    script:
        """
         filter_feature_table.R \\
                  --feature-table '${feature_table}' \\
                  --mode  '${meta.mode}' \\
                  --threshold  ${meta.filter_threshold} \\
                  --output-file  '${meta.output_file}'


        Rscript -e "VERSIONS=sprintf('tibble %s\\ntidyr %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\nmagrittr %s\\nglue %s\\n', \\
                                   packageVersion('tibble'), \\
                                   packageVersion('tidyr'), \\
                                   packageVersion('dplyr'), \\
                                   packageVersion('purrr'), \\
                                   packageVersion('readr'), \\
                                   packageVersion('stringr'), \\
                                   packageVersion('magrittr'), \\
                                   packageVersion('glue')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}

/*
 * ========================================================================================
 * PROCESS: ASSEMBLY_TABLE
 * ========================================================================================
 *
 * SUMMARY:
 *   Process Contig or Gene level taxonomy or KO table
 *
 * INPUTS:
 *   1. tuple: tuple val(type), val(level)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                   - type: feature type i.e taxonomy or KO
 *                   - level: assembly level i.e. 'Gene' or 'Contig'
 *
 *   2. path: feature_table
 *      Cardinality: one
 *      Description: Input file: feature (KO or taxonomy) table
 *
 *   3. path: summary_table
 *      Cardinality: one
 *      Description: Input file: assembly summary table used to retrieve sample names
 *
 * OUTPUTS:
 *   1. path: *.tsv (emit: table)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/R_visualizations.yaml
 *   Labels: R_downstream
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process ASSEMBLY_TABLE {
  
    tag "processing your ${level} ${type} table..."
    label "R_downstream"


    input:
        tuple val(type), val(level) // ['taxonomy', 'Contig']
        path(feature_table) // 'Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv'
        path(summary_table) // 'assembly-summaries_GLmetagenomics.tsv'

    output:
       path("*.tsv"), emit: table
       path("versions.txt"), emit: version

    script:
        """

          process_assembly_table.R \\
                  --assembly-table '${feature_table}' \\
                  --assembly-summary '${summary_table}' \\
                  --level  '${level}' \\
                  --type '${type}' \\
                  --output-prefix '${params.additional_filename_prefix}' \\
                  --assay-suffix '${params.assay_suffix}'

        Rscript -e "VERSIONS=sprintf('tibble %s\\ntidyr %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\nmagrittr %s\\nglue %s\\n', \\
                                   packageVersion('tibble'), \\
                                   packageVersion('tidyr'), \\
                                   packageVersion('dplyr'), \\
                                   packageVersion('purrr'), \\
                                   packageVersion('readr'), \\
                                   packageVersion('stringr'), \\
                                   packageVersion('magrittr'), \\
                                   packageVersion('glue')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"

        """
}

/*
 * ========================================================================================
 * PROCESS: HUMANN_TABLE
 * ========================================================================================
 *
 * SUMMARY:
 *   Process Humann generated pathway, uniref or KO table
 *
 * INPUTS:
 *   1. val: type
 *      Cardinality: one
 *      Description: Parameter value: feature type. one of 'pathway', 'uniref' or  'KO'.
 *
 *   2. path: feature_table
 *      Cardinality: one
 *      Description: Input file: feature table
 *
 * OUTPUTS:
 *   1. path: *.tsv (emit: table)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/R_visualizations.yaml
 *   Labels: R_downstream
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process HUMANN_TABLE {

    tag "processing your ${type} table..."
    label "R_downstream"


    input:
        val(type) // 'pathway', 'uniref' or  'KO'
        path(feature_table) // 'Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv'

    output:
       path("*.tsv"), emit: table
       path("versions.txt"), emit: version

    script:
        """
          process_humann_table.R \\
                  --table '${feature_table}' \\
                  --type '${type}' \\
                  --output-prefix '${params.additional_filename_prefix}' \\
                  --assay-suffix '${params.assay_suffix}'

        Rscript -e "VERSIONS=sprintf('tibble %s\\ntidyr %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\nmagrittr %s\\nglue %s\\n', \\
                                   packageVersion('tibble'), \\
                                   packageVersion('tidyr'), \\
                                   packageVersion('dplyr'), \\
                                   packageVersion('purrr'), \\
                                   packageVersion('readr'), \\
                                   packageVersion('stringr'), \\
                                   packageVersion('magrittr'), \\
                                   packageVersion('glue')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """
}


/*
 * ========================================================================================
 * PROCESS: DECONTAM
 * ========================================================================================
 *
 * SUMMARY:
 *   Statistical feature table decontamination with decontam
 *
 * INPUTS:
 *   1. val: meta
 *      Cardinality: one
 *      Description: Parameter value: A map with the following keys:
 *                  feature:  feature column name e.g. 'Species'
 *                  samples:  samples column name e.g. 'Sample_ID'
 *                  prevalence: prevalence column name e.g. 'NTC'
 *                  frequency: frequency column name e.g. 'concentration'
 *                  decontam_threshold: decontam's threshold for identifying contaminants. Default: 0.5
 *                  method: classification method e.g. 'kaiju'
 *                  ntc_name: name of ntc in prevalence column e.g. 'true'
 *
 *   2. path: metadata
 *      Cardinality: one
 *      Description: Input file: samples metadata file
 *
 *   3. path: feature_table
 *      Cardinality: one
 *      Description: Input file: feature (species, KO etc) table
 *
 * OUTPUTS:
 *   1. path: *_decontam_results*.tsv (emit: result)
 *
 *   2. path: *_decontam_*_table*.tsv (emit: table) [OPTIONAL]
 *
 *   3. path: *_decontam_failure.txt (emit: failure) [OPTIONAL] # A failure text file generated if the values in both prevalence and frequency columns are not different between samples within each column. i.e no difference between negative control(s) and other samples
 *
 *   4. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/R_visualizations.yaml
 *   Labels: R_downstream
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process DECONTAM  { 

    tag "Decontaminating ${feature_table} with decontam..."
    label "R_downstream"

    input:
        val(meta) // [feature: 'Species', 
              // samples: 'Sample_ID',
              // prevalence: 'Sample_or_Control',
              // frequency: 'concentration',
              // decontam_threshold: 0.1,
              // method: 'kaiju',
              // ntc_name: 'Control_Sample']
        path(metadata) // mapping/metadata.csv
        path(feature_table) // kaiju_species_table_GLlbnMetag.csv

    output:
        path("*_decontam_results*.tsv"), emit: result // decontam's primary results
        path("*_decontam_*_table*.tsv"), optional: true, emit: table // decontaminated feature table
        /* A failure text file generated if the values in both prevalence and frequency columns 
           are not different between samples within each column. 
           i.e no difference between negative control(s) and other samples */
        path("*_decontam_failure.txt"), optional: true, emit: failure 
        path("versions.txt"), emit: version

    script:
        """
        run_decontam.R \\
                  --feature-table '${feature_table}' \\
                  --feature-column '${meta.feature}' \\
                  --metadata-table '${metadata}' \\
                  --samples-column '${meta.samples}' \\
                  --prevalence-column '${meta.prevalence}' \\
                  --frequency-column '${meta.frequency}' \\
                  --ntc-name  '${meta.ntc_name}' \\
                  --threshold ${meta.decontam_threshold} \\
                  --classification-method '${meta.method}' \\
                  --output-prefix '${params.additional_filename_prefix}' \\
                  --assay-suffix '${params.assay_suffix}'

        Rscript -e "VERSIONS=sprintf('tibble %s\\ntidyr %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\nmagrittr %s\\nglue %s\\nphyloseq %s\\ndecontam %s\\n', \\
                                   packageVersion('tibble'), \\
                                   packageVersion('tidyr'), \\
                                   packageVersion('dplyr'), \\
                                   packageVersion('purrr'), \\
                                   packageVersion('readr'), \\
                                   packageVersion('stringr'), \\
                                   packageVersion('magrittr'), \\
                                   packageVersion('glue'), \\
                                   packageVersion('phyloseq'), \\
                                   packageVersion('decontam')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}

/*
 * ========================================================================================
 * PROCESS: BARPLOT
 * ========================================================================================
 *
 * SUMMARY:
 *   Make relative abundance static and interactive bar plots
 *
 * INPUTS:
 *   1. val: meta
 *      Cardinality: one
 *      Description: Parameter value:  A map with the following keys:
 *                  feature:  feature column name e.g. 'Species'
 *                  samples:  samples column name e.g. 'Sample_ID'
 *                  prefix: output file prefix e.g. 'kaiju_filtered_species'
 *
 *   2. path: feature_table
 *      Cardinality: one
 *      Description: Input file: read-based feature table to plot
 *
 *   3. path: metadata
 *      Cardinality: one
 *      Description: Input file: samples metadata file
 *
 * OUTPUTS:
 *   1. path: ${meta.prefix}_barplot${params.assay_suffix}.png (emit: plot)
 *
 *   2. path: ${meta.prefix}_barplot${params.assay_suffix}.html (emit: html)
 *
 *   3. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/R_visualizations.yaml
 *   Labels: R_downstream
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process BARPLOT {

    tag "Making your bar plot..."
    label "R_downstream"

    input:
      val(meta)
      path(feature_table) // 'kaiju_species_table_GLlbnMetag'
      path(metadata)  // 'mapping/metadata.txt'
        
    output:
         path("${meta.prefix}_barplot${params.assay_suffix}.png"), emit: plot
         path("${meta.prefix}_barplot${params.assay_suffix}.html"), emit: html
         path("versions.txt"), emit: version

    script:
        """
        # To fix Fontconfig error: No writable cache directories
        #mkdir -p cache/fontconfig/ && export FONTCONFIG_CACHE=cache/fontconfig/
        make_barplot.R \\
                  --metadata-table '${metadata}' \\
                  --feature-table '${feature_table}' \\
                  --group-column '${meta.group}' \\
                  --samples-column '${meta.samples}'  \\
                  --output-prefix  '${meta.prefix}' \\
                  --assay-suffix '${params.assay_suffix}'


        Rscript -e "VERSIONS=sprintf('tibble %s\\ntidyr %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\nmagrittr %s\\nglue %s\\nplotly %s\\nhtmlwidgets %s\\nggplot2 %s\\n', \\
                                   packageVersion('tibble'), \\
                                   packageVersion('tidyr'), \\
                                   packageVersion('dplyr'), \\
                                   packageVersion('purrr'), \\
                                   packageVersion('readr'), \\
                                   packageVersion('stringr'), \\
                                   packageVersion('magrittr'), \\
                                   packageVersion('glue'), \\
                                   packageVersion('plotly'), \\
                                   packageVersion('htmlwidgets'),
                                   packageVersion('ggplot2')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}


/*
 * ========================================================================================
 * PROCESS: HEATMAP
 * ========================================================================================
 *
 * SUMMARY:
 *   Make feature heatmap
 *
 * INPUTS:
 *   1. val: meta
 *      Cardinality: one
 *      Description: Parameter value: A map with the following keys:
 *                  group:  sample grouping column name e.g. 'group'
 *                  samples:  samples column name e.g. 'sample_id'
 *                  prefix: output file prefix e.g. 'Gene-families-uniref_unfiltered'
 *
 *   2. path: feature_table
 *      Cardinality: one
 *      Description: Input file: feature table
 *
 *   3. path: metadata
 *      Cardinality: one
 *      Description: Input file: samples metadata file
 *
 * OUTPUTS:
 *   1. path: *_heatmap${params.assay_suffix}.png (emit: plot)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/R_visualizations.yaml
 *   Labels: R_downstream
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process HEATMAP {

    tag "Making your heatmap..."
    label "R_downstream"


    input:
      val(meta) // [group: "group", prefix: "filtered_gene_functions"]
      path(feature_table) // 'kaiju_species_table_GLlbnMetag'
      path(metadata)  // 'mapping/metadata.txt'

        
    output:
       path("*_heatmap${params.assay_suffix}.png"), emit: plot
       path("versions.txt"), emit: version


    script:
        """
        make_heatmap.R \\
                  --metadata-table '${metadata}' \\
                  --feature-table '${feature_table}' \\
                  --samples-column '${meta.samples}' \\
                  --group-column '${meta.group}' \\
                  --output-prefix '${meta.prefix}' \\
                  --assay-suffix '${params.assay_suffix}'

        Rscript -e "VERSIONS=sprintf('tibble %s\\ntidyr %s\\ndplyr %s\\npurrr %s\\nreadr %s\\nstringr %s\\nmagrittr %s\\nglue %s\\npheatmap %s\\n', \\
                                   packageVersion('tibble'), \\
                                   packageVersion('tidyr'), \\
                                   packageVersion('dplyr'), \\
                                   packageVersion('purrr'), \\
                                   packageVersion('readr'), \\
                                   packageVersion('stringr'), \\
                                   packageVersion('magrittr'), \\
                                   packageVersion('glue'), \\
                                   packageVersion('pheatmap')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}

