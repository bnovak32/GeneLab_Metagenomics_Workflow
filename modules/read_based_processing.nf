#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
//params.additional_filename_prefix = ""
//params.assay_suffix = "_GLmetagenomics"

/**************************************************************************************** 
********************* Read-based processing using Humann3 *******************************
****************************************************************************************/
// Make Kraken, Kaiju and Humann3 databases
include { SETUP_KAIJU; SETUP_KRAKEN; make_humann_db } from "./database_creation.nf"

if(params.technology == "illumina"){
    // Metaphlan
    include { METAPHLAN2KRONA; KRONA_REPORT as METAPHLAN_REPORT } from "./visualize_taxonomy.nf"
    include { METAPHLAN2COUNT; BARPLOT as METAPHLAN_UNFILTERED_BARPLOT } from "./downstream_analysis.nf"
    include { FILTER_RARE as METAPHLAN_FILTER_RARE; BARPLOT as METAPHLAN_FILTERED_BARPLOT } from "./downstream_analysis.nf"
}

// Kraken2
include { KRAKEN_CLASSIFY; KRAKEN2TABLE } from "./assign_taxonomy.nf"
include { MULTIQC as KRAKEN_MULTIQC; ZIP_MULTIQC as KRAKEN_ZIP_MULTIQC } from "./quality_assessment.nf"
include { BARPLOT as KRAKEN_UNFILTERED_BARPLOT } from "./downstream_analysis.nf"
include { FILTER_RARE as KRAKEN_FILTER_RARE; BARPLOT as KRAKEN_FILTERED_BARPLOT } from "./downstream_analysis.nf"
// Kaiju
include { KAIJU_CLASSIFY; KAIJU2TABLE } from "./assign_taxonomy.nf"
include { KAIJU2SPECIES_TABLE; BARPLOT as KAIJU_UNFILTERED_BARPLOT } from "./downstream_analysis.nf"
include { FILTER_RARE as KAIJU_FILTER_RARE; BARPLOT as KAIJU_FILTERED_BARPLOT } from "./downstream_analysis.nf"
// Krona plots
include { KRONA_REPORT as KRAKEN_REPORT } from "./visualize_taxonomy.nf"
include { KRAKEN2KRONA; KAIJU2KRONA } from "./visualize_taxonomy.nf"
include { KRONA_REPORT as KAIJU_REPORT } from "./visualize_taxonomy.nf"

// Functional analysis

// Gene families UNIREF90
include { HUMANN_TABLE as GFU_HUMANN_TABLE; HEATMAP as GFU_UNFILTERED_HEATMAP } from "./downstream_analysis.nf"
include { FILTER_RARE as GFU_FILTER_RARE; HEATMAP as GFU_FILTERED_HEATMAP } from "./downstream_analysis.nf"

// Gene families KO
include { HUMANN_TABLE as GKO_HUMANN_TABLE; HEATMAP as GKO_UNFILTERED_HEATMAP } from "./downstream_analysis.nf"
include { FILTER_RARE as GKO_FILTER_RARE; HEATMAP as GKO_FILTERED_HEATMAP } from "./downstream_analysis.nf"

// Pathway Abundance
include { HUMANN_TABLE as PATH_HUMANN_TABLE; HEATMAP as PATH_UNFILTERED_HEATMAP } from "./downstream_analysis.nf"
include { FILTER_RARE as PATH_FILTER_RARE; HEATMAP as PATH_FILTERED_HEATMAP } from "./downstream_analysis.nf"


if(params.sample_type == "low_biomass"){

    include { DECONTAM as METAPHLAN_DECONTAM; BARPLOT as METAPHLAN_DECONTAM_BARPLOT } from "./downstream_analysis.nf"
    include { DECONTAM as KRAKEN_DECONTAM; BARPLOT as KRAKEN_DECONTAM_BARPLOT } from "./downstream_analysis.nf"
    include { DECONTAM as KAIJU_DECONTAM; BARPLOT as KAIJU_DECONTAM_BARPLOT } from "./downstream_analysis.nf"
    include { DECONTAM as GFU_DECONTAM; HEATMAP as GFU_DECONTAM_HEATMAP } from "./downstream_analysis.nf"
    include { DECONTAM as GKO_DECONTAM; HEATMAP as GKO_DECONTAM_HEATMAP } from "./downstream_analysis.nf"
    include { DECONTAM as PATH_DECONTAM; HEATMAP as PATH_DECONTAM_HEATMAP } from "./downstream_analysis.nf"
}

/*
    This process runs humann3 and metaphlan4 on each individual sample generating the
    read-based functional annotations and taxonomic classifications.
*/

process HUMANN {

    tag "Running humann on ${sample_id}-s reads..."
    label "read_based"
    

    input:
        tuple val(sample_id), path(reads), val(isPaired)
        path(chocophlan_dir)
        path(uniref_dir)
        path(metaphlan_dir)
        
    output:
        path("${sample_id}-humann3-out-dir/${sample_id}_genefamilies.tsv"), emit: genefamilies
        path("${sample_id}-humann3-out-dir/${sample_id}_pathabundance.tsv"), emit: pathabundance
        path("${sample_id}-humann3-out-dir/${sample_id}_pathcoverage.tsv"), emit: pathcoverage
        tuple val(sample_id), path("${sample_id}-humann3-out-dir/${sample_id}_metaphlan_bugs_list.tsv"), emit: metaphlan_bugs_list
        path("versions.txt"), emit: version
    script:
        """
        zcat ${reads} > ${sample_id}-reads.tmp.fq

        humann --input ${sample_id}-reads.tmp.fq \\
                    --output ${sample_id}-humann3-out-dir/ \\
                    --threads ${task.cpus} \\
                    --output-basename ${sample_id} \\
                    --metaphlan-options "--bowtie2db ${metaphlan_dir} --unclassified_estimation --add_viruses --sample_id ${sample_id}" \\
                    --nucleotide-database ${chocophlan_dir} \\
                    --protein-database ${uniref_dir} \\
                    --bowtie-options "--sensitive --mm" && \\
        mv ${sample_id}-humann3-out-dir/${sample_id}_humann_temp/${sample_id}_metaphlan_bugs_list.tsv \\
               ${sample_id}-humann3-out-dir/${sample_id}_metaphlan_bugs_list.tsv

        humann3 --version  > versions.txt
        """
}


/*
    This process combines the read-based humann3 output functional 
    tables from indiviual samples into single tables across the GLDS dataset.
*/

process COMBINE_READ_BASED_PROCESSING_TABLES {

    tag "Combining the read based processing tables..."
    label "read_based"

    input:
        path(gene_families)
        path(path_abundances)
        path(path_coverages)
        path(utilities_path)
    output:
        path("${params.additional_filename_prefix}gene-families.tsv"), emit: gene_families 
        path("${params.additional_filename_prefix}pathway-abundances.tsv"), emit: path_abundances
        path("${params.additional_filename_prefix}pathway-coverages.tsv"), emit: path_coverages
        path("versions.txt"), emit: version
    script:
        """
        if [ ${params.use_conda} == true ]; then
            # Setting humann3 utilities location (can be off if we pointed to
            # a previously installed database, and doesn't hurt to reset if it was already good-to-go)
            humann_config --update database_folders utility_mapping ${utilities_path} > /dev/null 2>&1
        fi

        # they each need to be in the same directories to be merged
        mkdir -p gene-family-results/ path-abundance-results/ path-coverage-results/
        cp ${gene_families} gene-family-results/ 
        cp ${path_abundances} path-abundance-results/
        cp ${path_coverages} path-coverage-results/

        humann_join_tables -i gene-family-results/ -o ${params.additional_filename_prefix}gene-families.tsv > /dev/null 2>&1
        humann_join_tables -i path-abundance-results/ -o ${params.additional_filename_prefix}pathway-abundances.tsv > /dev/null 2>&1
        humann_join_tables -i path-coverage-results/ -o ${params.additional_filename_prefix}pathway-coverages.tsv > /dev/null 2>&1

        humann3 --version  > versions.txt
        """
}


/*
    The read-based functional annotation tables have taxonomic info and non-taxonomic info mixed
    together initially. Humann comes with utility scripts to split these. This process does that,
    generating non-taxonomically grouped functional info files and taxonomically grouped ones.
*/

process SPLIT_READ_BASED_PROCESSING_TABLES {

    tag "Splitting humann stratified tables..."
    label "read_based"

    input:
        path(gene_families)
        path(path_abundances)
        path(path_coverages)
    output:
        path("${params.additional_filename_prefix}Gene-families${params.assay_suffix}.tsv"), emit: gene_families
        path("${params.additional_filename_prefix}Gene-families-grouped-by-taxa${params.assay_suffix}.tsv"), emit: gene_families_grouped
        path("${params.additional_filename_prefix}Pathway-abundances${params.assay_suffix}.tsv"), emit: path_abundances
        path("${params.additional_filename_prefix}Pathway-abundances-grouped-by-taxa${params.assay_suffix}.tsv"), emit: path_abundances_grouped
        path("${params.additional_filename_prefix}Pathway-coverages${params.assay_suffix}.tsv"), emit: path_coverages 
        path("${params.additional_filename_prefix}Pathway-coverages-grouped-by-taxa${params.assay_suffix}.tsv"), emit: path_coverages_grouped
        path("versions.txt"), emit: version
    script:
        """
        [ -d temp_processing/ ] && rm -rf temp_processing/
        mkdir temp_processing/

        # Gene Families
        humann_split_stratified_table -i ${gene_families} -o temp_processing/ > /dev/null 2>&1
        mv temp_processing/${params.additional_filename_prefix}gene-families_stratified.tsv \\
           ${params.additional_filename_prefix}Gene-families-grouped-by-taxa${params.assay_suffix}.tsv
        
        mv temp_processing/${params.additional_filename_prefix}gene-families_unstratified.tsv \\
           ${params.additional_filename_prefix}Gene-families${params.assay_suffix}.tsv

        # Pathway Abundance
        humann_split_stratified_table -i ${path_abundances} -o temp_processing/ > /dev/null 2>&1
        mv temp_processing/${params.additional_filename_prefix}pathway-abundances_stratified.tsv \\
           ${params.additional_filename_prefix}Pathway-abundances-grouped-by-taxa${params.assay_suffix}.tsv

        mv temp_processing/${params.additional_filename_prefix}pathway-abundances_unstratified.tsv \\
           ${params.additional_filename_prefix}Pathway-abundances${params.assay_suffix}.tsv

        # Pathway Coverage
        humann_split_stratified_table -i ${path_coverages} -o temp_processing/ > /dev/null 2>&1
        mv temp_processing/${params.additional_filename_prefix}pathway-coverages_stratified.tsv \\
           ${params.additional_filename_prefix}Pathway-coverages-grouped-by-taxa${params.assay_suffix}.tsv

        mv temp_processing/${params.additional_filename_prefix}pathway-coverages_unstratified.tsv \\
           ${params.additional_filename_prefix}Pathway-coverages${params.assay_suffix}.tsv

        humann3 --version  > versions.txt
        """
}


/*
    This process generates some normalized tables of the read-based functional outputs from
    humann that are more readily suitable for across sample comparisons.
*/

process GEN_NORMALIZED_READ_BASED_PROCESSING_TABLES {

    tag "Generating normalized humann tables..."
    label "read_based"

    input:
       path(gene_families)
       path(path_abundances)
    output:
        path("${params.additional_filename_prefix}Gene-families-cpm${params.assay_suffix}.tsv"), emit: gene_families
        path("${params.additional_filename_prefix}Pathway-abundances-cpm${params.assay_suffix}.tsv"), emit: path_abundances
        path("versions.txt"), emit: version
    script:
        """
        humann_renorm_table \\
               -i ${gene_families} \\
               -o ${params.additional_filename_prefix}Gene-families-cpm${params.assay_suffix}.tsv \\
               --update-snames > /dev/null 2>&1

        humann_renorm_table \\
               -i ${path_abundances} \\
               -o ${params.additional_filename_prefix}Pathway-abundances-cpm${params.assay_suffix}.tsv \\
               --update-snames > /dev/null 2>&1

        humann3 --version  > versions.txt
        """
}


/*
    This process summarizes the read-based humann annotations based on Kegg Orthlogy terms.
*/

process GEN_READ_BASED_PROCESSING_KO_TABLE {

    tag "Retrieving Kegg Orthologs..."
    label "read_based"
    
    input:
        path(gene_families)
    output:
        path("${params.additional_filename_prefix}Gene-families-KO-cpm${params.assay_suffix}.tsv"), emit: gene_families
        path("versions.txt"), emit: version
    script:
        """
        humann_regroup_table \\
              -i ${gene_families} \\
              -g uniref90_ko 2> /dev/null | \\
        humann_rename_table \\
               -n kegg-orthology 2> /dev/null | \\
        humann_renorm_table \\
               -o ${params.additional_filename_prefix}Gene-families-KO-cpm${params.assay_suffix}.tsv \\
               --update-snames > /dev/null 2>&1

        humann3 --version  > versions.txt
        """
}



//This process merges the taxonomy tables generated by metaphlan
process COMBINE_READ_BASED_PROCESSING_TAXONOMY {

    tag "Merging metaphlan taxonomy tables..."
    label "read_based"

    input:
        path(metaphlan_bugs_list_files)
    output:
        path("${params.additional_filename_prefix}metaphlan-taxonomy${params.assay_suffix}.tsv"), emit: taxonomy
        path("versions.txt"), emit: version
    script:
        """
        merge_metaphlan_tables.py ${metaphlan_bugs_list_files} \\
                         > ${params.additional_filename_prefix}metaphlan-taxonomy${params.assay_suffix}.tsv 2> /dev/null

        # Removing redundant text from headers 
        sed -i 's/_metaphlan_bugs_list//g' ${params.additional_filename_prefix}metaphlan-taxonomy${params.assay_suffix}.tsv

        metaphlan --version > versions.txt
        """
}



workflow read_based {

    take:
        reads_per_sample
        metadata
        filtered_reads
        krakendb_dir
        kaijudb_dir
        chocophlan_dir
        uniref_dir
        metaphlan_dir
        utilities_dir

    main:
        
        software_versions_ch = Channel.empty()

        
        // ------------------------ Kraken
        if(krakendb_dir){
            KRAKEN_CLASSIFY(krakendb_dir, filtered_reads)
        }else{
            SETUP_KRAKEN(params.krakendb_url)
            KRAKEN_CLASSIFY(SETUP_KRAKEN.out.krakendb_dir, filtered_reads)
            SETUP_KRAKEN.out.version | mix(software_versions_ch) | set{software_versions_ch}
        }
        kraken_reports = KRAKEN_CLASSIFY.out.report.map{sample_id, report -> report}.collect()
        KRAKEN_MULTIQC(Channel.of('kraken2'), params.multiqc_config, kraken_reports)
        KRAKEN_ZIP_MULTIQC(Channel.of('kraken2'), KRAKEN_MULTIQC.out.data)
        KRAKEN2TABLE(kraken_reports)
        KRAKEN2KRONA(KRAKEN_CLASSIFY.out.report)
        // Unfiltered
        unfilt_kraken_barplot_meta = Channel.of([group: "group",
                               feature: 'Species',
                               samples: 'sample_id',
                               prefix:  'kraken2_unfiltered_species'])
        KRAKEN_UNFILTERED_BARPLOT(unfilt_kraken_barplot_meta, KRAKEN2TABLE.out.table, metadata)
        // Filtered - drop species with relative abundance less than 0.5% across samples
        filt_kraken_meta = Channel.of([mode: 'across_samples', filter_threshold : 0.5,
                            output_file: "kraken2_filtered_species_table${params.assay_suffix}.tsv"])
        KRAKEN_FILTER_RARE(filt_kraken_meta, KRAKEN2TABLE.out.table)
        filt_kraken_barplot_meta = Channel.of([group: "group",
                               feature: 'Species',
                               samples: 'sample_id',
                               prefix:  'kraken2_filtered_species'])
        KRAKEN_FILTERED_BARPLOT(filt_kraken_barplot_meta, KRAKEN_FILTER_RARE.out.table, metadata)
        KRAKEN_REPORT("kraken2", KRAKEN2KRONA.out.krona.collect())

        // -------------------- Kaiju
        if(kaijudb_dir){
           kaijuDB = kaijudb_dir
           KAIJU_CLASSIFY(kaijuDB, filtered_reads)
        }else{
            SETUP_KAIJU(params.kaijudb_name)
            KAIJU_CLASSIFY(SETUP_KAIJU.out.kaijudb_dir, filtered_reads)
            kaijuDB = SETUP_KAIJU.out.kaijudb_dir
            SETUP_KAIJU.out.version | mix(software_versions_ch) | set{software_versions_ch}
        }
        kaiju_reports = KAIJU_CLASSIFY.out.report.map{sample_id, report -> report}.collect()
        KAIJU2TABLE(kaijuDB, "species", kaiju_reports)
        KAIJU2KRONA(kaijuDB, KAIJU_CLASSIFY.out.report)
        KAIJU_REPORT("kaiju", KAIJU2KRONA.out.krona.collect())

        // Unfiltered
        KAIJU2SPECIES_TABLE(KAIJU2TABLE.out.table)
        unfilt_kaiju_barplot_meta = Channel.of([group: "group",
                               feature: 'Species',
                               samples: 'sample_id',
                               prefix:  'kaiju_unfiltered_species'])
        KAIJU_UNFILTERED_BARPLOT(unfilt_kaiju_barplot_meta, KAIJU2SPECIES_TABLE.out.table, metadata)
        // Filtered - drop species with relative abundance less than 0.5% across samples
        filt_kaiju_meta = Channel.of([mode: 'across_samples', filter_threshold : 0.5,
                            output_file: "kaiju_filtered_species_table${params.assay_suffix}.tsv"])
        KAIJU_FILTER_RARE(filt_kaiju_meta, KAIJU2SPECIES_TABLE.out.table)
        filt_kaiju_barplot_meta = Channel.of([group: "group",
                               feature: 'Species',
                               samples: 'sample_id',
                               prefix:  'kaiju_filtered_species'])
        KAIJU_FILTERED_BARPLOT(filt_kaiju_barplot_meta, KAIJU_FILTER_RARE.out.table, metadata)


        // Humann
        if(chocophlan_dir && uniref_dir && metaphlan_dir && utilities_dir){
            HUMANN(filtered_reads, chocophlan_dir, uniref_dir, metaphlan_dir)
        }else{
            make_humann_db()
            HUMANN(filtered_reads, make_humann_db.out.chocophlan_dir, 
                   make_humann_db.out.uniref_dir,
                   make_humann_db.out.metaphlan_db_dir) 
            make_humann_db.out.versions | mix(software_versions_ch) | set{software_versions_ch} 
        } 

        gene_families_ch = HUMANN.out.genefamilies.collect()
        pathabundance_ch = HUMANN.out.pathabundance.collect()
        pathcoverage_ch = HUMANN.out.pathcoverage.collect()
        metaphlan_bugs_list_ch = HUMANN.out.metaphlan_bugs_list.map{sample_id, bug_list -> bug_list}.collect()
        
        if(chocophlan_dir && uniref_dir && metaphlan_dir && utilities_dir){
            COMBINE_READ_BASED_PROCESSING_TABLES(gene_families_ch, pathabundance_ch, pathcoverage_ch, utilities_dir)
        }else{
            COMBINE_READ_BASED_PROCESSING_TABLES(gene_families_ch, pathabundance_ch,
                                              pathcoverage_ch, make_humann_db.out.utilities_dir)
        }
        SPLIT_READ_BASED_PROCESSING_TABLES(COMBINE_READ_BASED_PROCESSING_TABLES.out.gene_families,
                                           COMBINE_READ_BASED_PROCESSING_TABLES.out.path_abundances,
                                           COMBINE_READ_BASED_PROCESSING_TABLES.out.path_coverages)

        GEN_NORMALIZED_READ_BASED_PROCESSING_TABLES(SPLIT_READ_BASED_PROCESSING_TABLES.out.gene_families,
                                                    SPLIT_READ_BASED_PROCESSING_TABLES.out.path_abundances)

        GEN_READ_BASED_PROCESSING_KO_TABLE(SPLIT_READ_BASED_PROCESSING_TABLES.out.gene_families)
        ko_table_ch      = GEN_READ_BASED_PROCESSING_KO_TABLE.out.gene_families
        uniref_table_ch  =  GEN_NORMALIZED_READ_BASED_PROCESSING_TABLES.out.gene_families
        pathway_table_ch =  GEN_NORMALIZED_READ_BASED_PROCESSING_TABLES.out.path_abundances
       

        // ------------------------- Gene families UNIREF 90
        // Unfiltered
        GFU_HUMANN_TABLE('uniref', uniref_table_ch)
        unfilt_uniref_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                               prefix:  'Gene-families-uniref_unfiltered']) 
        GFU_UNFILTERED_HEATMAP(unfilt_uniref_heatmap_meta, GFU_HUMANN_TABLE.out.table, metadata)
        // Filtered - filter out uniref less than 1000 CPM across samples
        filt_uniref_meta = Channel.of([mode: 'values_sum', filter_threshold : 1000,
                            output_file: "Gene-families-uniref_filtered${params.assay_suffix}.tsv"])
        GFU_FILTER_RARE(filt_uniref_meta, GFU_HUMANN_TABLE.out.table)
        filt_uniref_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                             prefix:  'Gene-families-uniref_filtered'])
        GFU_FILTERED_HEATMAP(filt_uniref_heatmap_meta, GFU_FILTER_RARE.out.table, metadata)

        // ------------------------- Gene families KO
        // Unfiltered
        GKO_HUMANN_TABLE('KO', ko_table_ch)
        unfilt_KO_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                               prefix:  'Gene-families-KO_unfiltered'])
        GKO_UNFILTERED_HEATMAP(unfilt_KO_heatmap_meta, GKO_HUMANN_TABLE.out.table, metadata)
        // Filtered - filter out KO less than 1000 CPM across samples
        filt_KO_meta = Channel.of([mode: 'values_sum', filter_threshold : 1000,
                            output_file: "Gene-families-KO_filtered${params.assay_suffix}.tsv"])
        GKO_FILTER_RARE(filt_KO_meta, GKO_HUMANN_TABLE.out.table)
        filt_KO_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                             prefix:  'Gene-families-KO_filtered'])
        GKO_FILTERED_HEATMAP(filt_KO_heatmap_meta, GKO_FILTER_RARE.out.table, metadata)

        // ------------------------- Pathway abundances
        // Unfiltered
        PATH_HUMANN_TABLE('pathway', pathway_table_ch)
        unfilt_pathway_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                               prefix:  'Pathway-abundances_unfiltered'])
        PATH_UNFILTERED_HEATMAP(unfilt_pathway_heatmap_meta, PATH_HUMANN_TABLE.out.table, metadata)
        // Filtered - filter out pathways less than 1000 CPM across samples
        filt_pathway_meta = Channel.of([mode: 'values_sum', filter_threshold : 1000,
                            output_file: "Pathway-abundances_filtered${params.assay_suffix}.tsv"])
        PATH_FILTER_RARE(filt_pathway_meta, PATH_HUMANN_TABLE.out.table)
        filt_pathway_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                             prefix:  'Pathway-abundances_filtered'])
        PATH_FILTERED_HEATMAP(filt_pathway_heatmap_meta, PATH_FILTER_RARE.out.table, metadata)

 
        COMBINE_READ_BASED_PROCESSING_TAXONOMY(metaphlan_bugs_list_ch)
        taxonomy_ch = COMBINE_READ_BASED_PROCESSING_TAXONOMY.out.taxonomy

        if(params.technology == "illumina"){

       // ------------------- Metaphlan
        METAPHLAN2KRONA(HUMANN.out.metaphlan_bugs_list)
        METAPHLAN_REPORT("metaphlan", METAPHLAN2KRONA.out.krona.collect())

        // Unfiltered
        // Create raw table
        METAPHLAN2COUNT(taxonomy_ch, reads_per_sample)
        unfilt_metaphlan_barplot_meta = Channel.of([group: "group",
                               feature: 'Species',
                               samples: 'sample_id',
                               prefix:  'metaplan_unfiltered_species'])
        METAPHLAN_UNFILTERED_BARPLOT(unfilt_metaphlan_barplot_meta, METAPHLAN2COUNT.out.table, metadata)
        // Filtered - drop species with relative abundance less than 0.5% across samples
        filt_metaphlan_meta = Channel.of([mode: 'across_samples', filter_threshold : 0.5,
                            output_file: "metaphlan_filtered_species_table${params.assay_suffix}.tsv"])
        METAPHLAN_FILTER_RARE(filt_metaphlan_meta, METAPHLAN2COUNT.out.table)
        filt_metaphlan_barplot_meta = Channel.of([group: "group",
                               feature: 'Species',
                               samples: 'sample_id',
                               prefix:  'metaplan_filtered_species'])
        METAPHLAN_FILTERED_BARPLOT(filt_metaphlan_barplot_meta, METAPHLAN_FILTER_RARE.out.table, metadata)
        
        if(params.sample_type == "low_biomass"){
            // Decontaminate with decontam
            decontam_metaphlan_meta = Channel.of([feature: 'Species', samples: 'sample_id',
                                   prevalence: 'NTC', frequency: 'concentration',
                                   decontam_threshold: 0.5, method: 'metaphlan',
                                   ntc_name: 'true'])
            METAPHLAN_DECONTAM(decontam_metaphlan_meta, metadata, METAPHLAN_FILTER_RARE.out.table)
            decontam_metaphlan_barplot_meta = Channel.of([group: "group",
                               feature: 'Species',
                               samples: 'sample_id',
                               prefix:  'metaplan_decontam_species'])
            METAPHLAN_DECONTAM_BARPLOT(decontam_metaphlan_barplot_meta, METAPHLAN_DECONTAM.out.table, metadata)
        }
         
        METAPHLAN_REPORT.out.version | mix(software_versions_ch) | set{software_versions_ch}
        }


        // Decontaminate with decontam
        if(params.sample_type == "low_biomass"){

        // Kraken2
        decontam_kraken_meta = Channel.of([feature: 'Species', samples: 'sample_id',
                                   prevalence: 'NTC', frequency: 'concentration',
                                   decontam_threshold: 0.5, method: 'kraken2',
                                   ntc_name: 'true'])
        KRAKEN_DECONTAM(decontam_kraken_meta, metadata, KRAKEN_FILTER_RARE.out.table)
        decontam_kraken_barplot_meta = Channel.of([group: "group",
                               feature: 'Species',
                               samples: 'sample_id',
                               prefix:  'kraken2_decontam_species'])
        KRAKEN_DECONTAM_BARPLOT(decontam_kraken_barplot_meta, KRAKEN_DECONTAM.out.table, metadata)


        // Kaiju
        decontam_kaiju_meta = Channel.of([feature: 'Species', samples: 'sample_id',
                                   prevalence: 'NTC', frequency: 'concentration',
                                   decontam_threshold: 0.5, method: 'kaiju',
                                   ntc_name: 'true'])
        KAIJU_DECONTAM(decontam_kaiju_meta, metadata, KAIJU_FILTER_RARE.out.table)
        decontam_kaiju_barplot_meta = Channel.of([group: "group",
                               feature: 'Species',
                               samples: 'sample_id',
                               prefix:  'kaiju_decontam_species'])
        KAIJU_DECONTAM_BARPLOT(decontam_kaiju_barplot_meta, KAIJU_DECONTAM.out.table, metadata)


        // Gene families UNIREF90
        decontam_uniref_meta = Channel.of([feature: 'Uniref90', samples: 'sample_id',
                                   prevalence: 'NTC', frequency: 'concentration',
                                   decontam_threshold: 0.5, method: 'Gene-families-uniref',
                                   ntc_name: 'true'])

        GFU_DECONTAM(decontam_uniref_meta, metadata, GFU_FILTER_RARE.out.table)
        decontam_uniref_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                                  prefix:  'Gene-families-uniref_decontam'])        
        GFU_DECONTAM_HEATMAP(decontam_uniref_heatmap_meta, GFU_DECONTAM.out.table, metadata)

        // Gene families KO
        decontam_KO_meta = Channel.of([feature: 'KO', samples: 'sample_id',
                                   prevalence: 'NTC', frequency: 'concentration',
                                   decontam_threshold: 0.5, method: 'Gene-families-KO',
                                   ntc_name: 'true'])

        GKO_DECONTAM(decontam_KO_meta, metadata, GKO_FILTER_RARE.out.table)
        decontam_KO_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                                  prefix:  'Gene-families-KO_decontam'])
        GKO_DECONTAM_HEATMAP(decontam_KO_heatmap_meta, GKO_DECONTAM.out.table, metadata)


        // Pathway
        decontam_pathway_meta = Channel.of([feature: 'Pathway', samples: 'sample_id',
                                   prevalence: 'NTC', frequency: 'concentration',
                                   decontam_threshold: 0.5, method: 'Pathway-abundances',
                                   ntc_name: 'true'])

        PATH_DECONTAM(decontam_pathway_meta, metadata, PATH_FILTER_RARE.out.table)
        decontam_pathway_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                                  prefix:  'Pathway-abundances_decontam'])
        PATH_DECONTAM_HEATMAP(decontam_pathway_heatmap_meta, PATH_DECONTAM.out.table, metadata)



        KAIJU_DECONTAM.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KAIJU_DECONTAM_BARPLOT.out.version | mix(software_versions_ch) | set{software_versions_ch} 

        }


        KRAKEN_CLASSIFY.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KRAKEN_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KRAKEN_ZIP_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KRAKEN2TABLE.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KRAKEN2KRONA.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KRAKEN_REPORT.out.version  | mix(software_versions_ch) | set{software_versions_ch}
        KAIJU_CLASSIFY.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KAIJU2TABLE.out.version  | mix(software_versions_ch) | set{software_versions_ch}
        KAIJU2KRONA.out.version  | mix(software_versions_ch) | set{software_versions_ch}
        KAIJU_REPORT.out.version | mix(software_versions_ch) | set{software_versions_ch}
        HUMANN.out.version | mix(software_versions_ch) | set{software_versions_ch}
        
        COMBINE_READ_BASED_PROCESSING_TABLES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SPLIT_READ_BASED_PROCESSING_TABLES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        GEN_NORMALIZED_READ_BASED_PROCESSING_TABLES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        GEN_READ_BASED_PROCESSING_KO_TABLE.out.version | mix(software_versions_ch) | set{software_versions_ch}
        COMBINE_READ_BASED_PROCESSING_TAXONOMY.out.version | mix(software_versions_ch) | set{software_versions_ch}

    emit:
        gene_families   = uniref_table_ch
        path_abundances = pathway_table_ch
        ko_table        = ko_table_ch
        taxonomy        = taxonomy_ch 
        versions        = software_versions_ch
}


workflow {

     read_based(reads_per_sample,
                metadata,
                filtered_ch, 
                params.krakendb_dir,
                params.kaijudb_dir,
                params.chocophlan_dir,
                params.uniref_dir,
                params.metaphlan_db_dir,
                params.utilities_dir)
}
