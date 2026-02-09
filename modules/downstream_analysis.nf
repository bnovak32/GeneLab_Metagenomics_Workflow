#!/usr/bin/env nextflow
nextflow.enable.dsl = 2



process METAPHLAN2COUNT {

    tag "Processing metaphlan count table..."
    label "R_downstream"
    label "read_based_outputs"

    input:
        path(metaphlan_table) // Metaphlan-taxonomy_GLmetagenomics.tsv
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

        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\n', \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}



process KAIJU2SPECIES_TABLE  { 

    tag "Processing kaiju species table..."
    label "R_downstream"
    label "read_based_outputs"


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

        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\n', \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}


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


        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\n', \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}


process ASSEMBLY_TABLE {
  
    tag "processing your ${level} ${type} table..."
    label "R_downstream"
    label "combine_outputs"


    input:
        tuple val(type), val(level) // ['taxonomy', 'Contig']
        path(feature_table) // 'Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv'
        path(summary_table) // 'assembly-summaries_GLmetagenomics.tsv'

    output:
       path("${params.additional_filename_prefix}${level}-level*.tsv"), emit: table
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

          Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\n',  \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """
}


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
        path("*_decontam_*_results*.tsv"), emit: result // decontam's primary results
        path("*_decontam_*_table*.tsv"), optional: true, emit: table // decontaminated feature table
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



        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\nphyloseq %s\\ndecontam %s\\n', \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue'), \\
                                    packageVersion('phyloseq'), \\
                                    packageVersion('decontam')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}



process BARPLOT {

    tag "Making your bar plot..."
    label "R_downstream"
    label "read_based_outputs"

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
        make_barplot.R \\
                  --metadata-table '${metadata}' \\
                  --feature-table '${feature_table}' \\
                  --group-column '${meta.group}' \\
                  --feature-column '${meta.feature}' \\
                  --samples-column '${meta.samples}'  \\
                  --output-prefix  '${meta.prefix}' \\
                  --assay-suffix '${params.assay_suffix}'

        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\nplotly %s\\nhtmlwidgets %s\\n', \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue'), \\
                                    packageVersion('plotly'), \\
                                    packageVersion('htmlwidgets')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}


process HEATMAP {

    tag "Making your heatmap..."
    label "R_downstream"
    label "combine_outputs"


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


        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\npheatmap %s\\n', \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue'), \\
                                    packageVersion('pheatmap')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}

