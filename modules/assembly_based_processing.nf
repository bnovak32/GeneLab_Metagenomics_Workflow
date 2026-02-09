#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
***************************  Assembly-based processing workflow *************************
****************************************************************************************/

// Assembly-based workflow
// Processes to create the required database(s) if not provided
include { SETUP_CAT_DB; SETUP_KOFAMSCAN_DB; SETUP_GTDBTK_DB;
          SETUP_CHOCOPHLAN} from "./database_creation.nf"

if(params.technology == "illumina"){

include { ASSEMBLE } from "./assembly.nf"
include { MAPPING as SHORT_MAPPING; SAM_TO_BAM } from "./read_mapping.nf"

}else{

    include { FLYE ; POLISH_ASSEMBLY } from "./assembly.nf"
    include { LONG_MAPPING } from "./long_read_mapping.nf"
    include { SAMTOOLS_SORT as SAM_TO_BAM } from "./samtools.nf"
}

include { RENAME_HEADERS; SUMMARIZE_ASSEMBLIES } from "./assembly.nf"
include { CALL_GENES; REMOVE_LINEWRAPS } from "./assembly_annotation.nf"
include { KO_ANNOTATION; FILTER_KFAMSCAN } from "./assembly_annotation.nf"
include { TAX_CLASSIFICATION } from "./assembly_annotation.nf"
include { GET_COV_AND_DET } from "./coverage.nf"
include { COMBINE_GENE_ANNOTS_TAX_AND_COVERAGE; MAKE_COMBINED_GENE_LEVEL_TABLES } from "./combine_contig_annotation.nf"

// Gene taxonomy
include { ASSEMBLY_TABLE as GT_ASSEMBLY_TABLE; HEATMAP as GT_UNFILTERED_HEATMAP } from "./downstream_analysis.nf"
include { FILTER_RARE as GT_FILTER_RARE; HEATMAP as GT_FILTERED_HEATMAP } from "./downstream_analysis.nf"
include { DECONTAM as GT_DECONTAM; HEATMAP as GT_DECONTAM_HEATMAP } from "./downstream_analysis.nf"

// Gene functions/KO
include { ASSEMBLY_TABLE as GF_ASSEMBLY_TABLE; HEATMAP as GF_UNFILTERED_HEATMAP } from "./downstream_analysis.nf"
include { FILTER_RARE as GF_FILTER_RARE; HEATMAP as GF_FILTERED_HEATMAP } from "./downstream_analysis.nf"
include { DECONTAM as GF_DECONTAM; HEATMAP as GF_DECONTAM_HEATMAP } from "./downstream_analysis.nf"

// Contig taxonomy
include { COMBINE_CONTIG_TAX_AND_COVERAGE; MAKE_COMBINED_CONTIG_TAX_TABLES } from "./combine_contig_annotation.nf"
include { ASSEMBLY_TABLE as CONTIG_ASSEMBLY_TABLE; HEATMAP as CONTIG_UNFILTERED_HEATMAP } from "./downstream_analysis.nf"
include { FILTER_RARE as CONTIG_FILTER_RARE; HEATMAP as CONTIG_FILTERED_HEATMAP } from "./downstream_analysis.nf"
include { DECONTAM as CONTIG_DECONTAM; HEATMAP as CONTIG_DECONTAM_HEATMAP } from "./downstream_analysis.nf"

include { METABAT_BINNING } from "./binning.nf"
include { summarize_bins } from "./summarize_bins.nf"
include { summarize_mags } from "./summarize_MAG.nf"
include { GENERATE_ASSEMBLY_PROCESSING_OVERVIEW_TABLE } from "./summarize_assembly-based_processing.nf"


workflow assembly_based {

    take:
        metadata
        file_ch
        filtered_ch
        ko_db_dir
        cat_db
        gtdbtk_db_dir
        use_gtdbtk_scratch_location

    main:
        /*****************************************************
        ************* Assembly-based analysis ****************
        *****************************************************/

        software_versions_ch = Channel.empty()

        

        if(params.technology == "illumina"){
            // Assemble reads to contigs
            ASSEMBLE(filtered_ch)
            // Rename assembly header
            ASSEMBLE.out.contigs | RENAME_HEADERS
            ASSEMBLE.out.version | mix(software_versions_ch) | set{software_versions_ch}
        }else{
            // Assemble reads to contigs
            FLYE(filtered_ch)          
            // Polish assembly then rename assembly header
            raw_assembly_ch = FLYE.out.contigs
            POLISH_ASSEMBLY(raw_assembly_ch.join(filtered_ch))
            POLISH_ASSEMBLY.out.contigs | RENAME_HEADERS
            FLYE.out.version | mix(software_versions_ch) | set{software_versions_ch}
            POLISH_ASSEMBLY.out.version | mix(software_versions_ch) | set{software_versions_ch}
        }


        assembly_ch = RENAME_HEADERS.out.contigs
        assemblies_ch = assembly_ch.map{
                                      sample_id, assembly -> file("${assembly}")
                                      }.collect()
        SUMMARIZE_ASSEMBLIES(assemblies_ch)

        // Write failed assemblies to a Failed assemblies file
        failed_assemblies = RENAME_HEADERS.out.failed_assembly
        failed_assemblies
              .map{ it.text }
              .collectFile(name: "${params.assemblies_dir}/Failed-assemblies.tsv", cache: false)
        
        // Map reads to assembly
        if(params.technology == "illumina"){
            SHORT_MAPPING(assembly_ch.join(filtered_ch))
            SHORT_MAPPING.out.sam | SAM_TO_BAM
            SHORT_MAPPING.out.version | mix(software_versions_ch) | set{software_versions_ch}
        }else{
            LONG_MAPPING(assembly_ch.join(filtered_ch))
            LONG_MAPPING.out.sam | SAM_TO_BAM
            LONG_MAPPING.out.version | mix(software_versions_ch) | set{software_versions_ch}
        }
        read_mapping_ch = SAM_TO_BAM.out.bam

        // Annotate assembly
        CALL_GENES(assembly_ch)
        CALL_GENES.out.genes | REMOVE_LINEWRAPS
        genes_ch = REMOVE_LINEWRAPS.out.genes

        if (ko_db_dir){

            KO_ANNOTATION(assembly_ch.join(genes_ch), ko_db_dir)
            KO_ANNOTATION.out.temp_table | FILTER_KFAMSCAN
            annotations_ch = FILTER_KFAMSCAN.out.ko_annotation

        }else{

            SETUP_KOFAMSCAN_DB()
            SETUP_KOFAMSCAN_DB.out.version | mix(software_versions_ch) | set{software_versions_ch}
            KO_ANNOTATION(assembly_ch.join(genes_ch), SETUP_KOFAMSCAN_DB.out.ko_db_dir)
            KO_ANNOTATION.out.temp_table | FILTER_KFAMSCAN
            annotations_ch = FILTER_KFAMSCAN.out.ko_annotation

        }

        if (cat_db){         

            TAX_CLASSIFICATION(assembly_ch.join(genes_ch), cat_db)
            taxonomy_ch = TAX_CLASSIFICATION.out.taxonomy

        }else{

            SETUP_CAT_DB(params.CAT_DB_LINK)
            SETUP_CAT_DB.out.version | mix(software_versions_ch) | set{software_versions_ch}
            TAX_CLASSIFICATION(assembly_ch.join(genes_ch), SETUP_CAT_DB.out.cat_db)
            taxonomy_ch = TAX_CLASSIFICATION.out.taxonomy
        }

        // Calculate gene coverage and depth
        GET_COV_AND_DET(read_mapping_ch
                           .join(assembly_ch)
                           .join(genes_ch))
        coverage_ch = GET_COV_AND_DET.out.coverages

        // Combine contig annotations
        tax_and_cov_ch = COMBINE_GENE_ANNOTS_TAX_AND_COVERAGE(coverage_ch
                                                                .join(annotations_ch)
                                                                .join(taxonomy_ch)
                                                                .join(genes_ch)
                                                                .join(assembly_ch))
                                                                
        gene_coverage_annotation_and_tax_files_ch = tax_and_cov_ch.map{
                                                         sample_id, coverage -> file("${coverage}")
                                                         }.collect()

        // Gene level
        MAKE_COMBINED_GENE_LEVEL_TABLES(gene_coverage_annotation_and_tax_files_ch)

        // ------------------ Gene taxonomy
        // Unfiltered
        gene_taxonomy_ch = Channel.of(['taxonomy', 'Gene'])
        GT_ASSEMBLY_TABLE(gene_taxonomy_ch, MAKE_COMBINED_GENE_LEVEL_TABLES.out.norm_taxonomy_coverages, 
                          SUMMARIZE_ASSEMBLIES.out.summary)
        unfilt_gene_taxonomy_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                               prefix:  'NoFilter-Combined-gene-level-taxonomy']) 
        GT_UNFILTERED_HEATMAP(unfilt_gene_taxonomy_heatmap_meta, GT_ASSEMBLY_TABLE.out.table, metadata)
        // Filtered - filter out species less than 1000 CPM across samples
        filt_gene_taxonomy_meta = Channel.of([mode: 'values_sum', filter_threshold : 1000,
                            output_file: "Filter-Combined-gene-level-taxonomy${params.assay_suffix}.tsv"])
        GT_FILTER_RARE(filt_gene_taxonomy_meta, GT_ASSEMBLY_TABLE.out.table)
        filt_gene_taxonomy_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                             prefix:  'Filter-Combined-gene-level-taxonomy'])
        GT_FILTERED_HEATMAP(filt_gene_taxonomy_heatmap_meta, GT_FILTER_RARE.out.table, metadata)
        // Decontaminate with decontam
        decontam_gene_taxonomy_meta = Channel.of([feature: 'species', samples: 'sample_id',
                                   prevalence: 'NTC', frequency: 'concentration',
                                   decontam_threshold: 0.5, method: 'gene-taxonomy',
                                   ntc_name: 'true'])

        GT_DECONTAM(decontam_gene_taxonomy_meta, metadata, GT_FILTER_RARE.out.table)
        decontam_gene_taxonomy_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                                  prefix:  'Decontam-Combined-gene-level-taxonomy'])        
        GT_DECONTAM_HEATMAP(decontam_gene_taxonomy_heatmap_meta, GT_FILTER_RARE.out.table, metadata)

       // ---------------------- Gene functions
       gene_function_ch = Channel.of(['KO', 'Gene'])
       // Unfiltered
       GF_ASSEMBLY_TABLE(gene_function_ch, MAKE_COMBINED_GENE_LEVEL_TABLES.out.norm_function_coverages,
                         SUMMARIZE_ASSEMBLIES.out.summary)
       unfilt_gene_function_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                        prefix:  'NoFilter-Combined-gene-level-KO'])
       GF_UNFILTERED_HEATMAP(unfilt_gene_function_heatmap_meta, GF_ASSEMBLY_TABLE.out.table, metadata)

       // Filtered - filter out KO_ID less than 1000 CPM across samples
       filt_gene_function_meta = Channel.of([mode: 'values_sum', filter_threshold : 1000,
                            output_file: "Filter-Combined-gene-level-KO${params.assay_suffix}.tsv"])
       GF_FILTER_RARE(filt_gene_function_meta, GF_ASSEMBLY_TABLE.out.table)
       filt_gene_function_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                             prefix:  'Filter-Combined-gene-level-KO'])
       GF_FILTERED_HEATMAP(filt_gene_function_heatmap_meta, GF_FILTER_RARE.out.table, metadata)
       // Decontaminate with decontam
       decontam_gene_function_meta = Channel.of([feature: 'KO_ID', samples: 'sample_id',
                                   prevalence: 'NTC', frequency: 'concentration',
                                   decontam_threshold: 0.5, method: 'gene-function',
                                   ntc_name: 'true']) 
       GF_DECONTAM(decontam_gene_function_meta, metadata, GF_FILTER_RARE.out.table)
       decontam_gene_function_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                                  prefix:  'Decontam-Combined-gene-level-KO'])
       GF_DECONTAM_HEATMAP(decontam_gene_function_heatmap_meta, GF_FILTER_RARE.out.table, metadata)

   
        combined_cov_ch = COMBINE_CONTIG_TAX_AND_COVERAGE(coverage_ch
                                                            .join(taxonomy_ch)
                                                            .join(genes_ch)
                                                            .join(assembly_ch))

        // ------------------------ Contig taxonomy
        MAKE_COMBINED_CONTIG_TAX_TABLES(combined_cov_ch.map{
                                               sample_id, coverage -> file("${coverage}")
                                               }.collect())

       // Unfiltered
       contig_ch = Channel.of(['taxonomy', 'Contig'])
       CONTIG_ASSEMBLY_TABLE(contig_ch, MAKE_COMBINED_CONTIG_TAX_TABLES.out.norm_taxonomy, SUMMARIZE_ASSEMBLIES.out.summary)
       unfilt_contig_heatmap_meta = Channel.of([group: "group", samples: 'sample_id', 
                                        prefix:  'NoFilter-Combined-contig-level-taxonomy'])
       CONTIG_UNFILTERED_HEATMAP(unfilt_contig_heatmap_meta, CONTIG_ASSEMBLY_TABLE.out.table, metadata)
       // Filtered - filter at 1000 CPM
       filt_contig_meta = Channel.of([mode: 'values_sum', filter_threshold : 1000, 
                            output_file: "Filter-Combined-contig-level-taxonomy${params.assay_suffix}.tsv"])
       CONTIG_FILTER_RARE(filt_contig_meta, CONTIG_ASSEMBLY_TABLE.out.table)
       filt_contig_heatmap_meta = Channel.of([group: "group", samples: 'sample_id',
                                             prefix:  'Filter-Combined-contig-level-taxonomy'])
       CONTIG_FILTERED_HEATMAP(filt_contig_heatmap_meta, CONTIG_FILTER_RARE.out.table, metadata)
       // Decontaminated - using decontam threshold of 0.5 sice it was more likely to detect contaminats than 0.1.
       decontam_contig_meta = Channel.of([feature: 'species', samples: 'sample_id',
                                   prevalence: 'NTC', frequency: 'concentration',
                                   decontam_threshold: 0.5, method: 'contig-taxonomy', 
                                   ntc_name: 'true'])
       CONTIG_DECONTAM(decontam_contig_meta, metadata, CONTIG_FILTER_RARE.out.table)
       decontam_contig_heatmap_meta = Channel.of([group: "group", samples: 'sample_id', 
                                                  prefix:  'Decontam-Combined-contig-level-taxonomy'])
       CONTIG_DECONTAM_HEATMAP(decontam_contig_heatmap_meta, CONTIG_DECONTAM.out.table, metadata)


        // Assembly binning
        METABAT_BINNING(assembly_ch.join(read_mapping_ch))
        binning_ch = METABAT_BINNING.out.bins
        binning_ch | summarize_bins
        depth_ch = METABAT_BINNING.out.depth
        metabat_assembly_depth_files_ch = depth_ch.map{
                                           sample_id, depth -> file("${depth}")
                                           }.collect()
        bins_ch = binning_ch.map{
                         sample_id, bins -> bins instanceof List ? bins.each{it}: bins 
                         }.flatten().collect()
   
         
        // Check Bins and Summarize MAGs
        if(gtdbtk_db_dir){
            summarize_mags(summarize_bins.out.bins_checkm_results,
                      bins_ch,
                      gtdbtk_db_dir, use_gtdbtk_scratch_location,
                      gene_coverage_annotation_and_tax_files_ch)
        }else{
            SETUP_GTDBTK_DB(params.GTDBTK_LINK)
            SETUP_GTDBTK_DB.out.version | mix(software_versions_ch) | set{software_versions_ch}
            summarize_mags(summarize_bins.out.bins_checkm_results,
                      bins_ch, 
                      SETUP_GTDBTK_DB.out.gtdbtk_db_dir, use_gtdbtk_scratch_location,
                      gene_coverage_annotation_and_tax_files_ch)          
        }

        // Get the predicted amino acids for all the samples
        genes_aa_ch = genes_ch.map{sample_id, aa, nt -> file("${aa}")}.collect() 
    
        // Generating a file with sample ids on a new line
        file_ch.map{row -> "${row.sample_id}"}
              .collectFile(name: "${launchDir}/unique-sample-IDs.txt", newLine: true)
              .set{sample_ids_ch}


        bam_files = read_mapping_ch.map{sample_id, bam -> file("${bam}")}.collect()
        // Summarize Assembly-based analysis
        GENERATE_ASSEMBLY_PROCESSING_OVERVIEW_TABLE(sample_ids_ch, summarize_mags.out.MAGs_overview,
                                                    summarize_mags.out.MAGs_dir, assemblies_ch,
                                                    genes_aa_ch,
                                                    metabat_assembly_depth_files_ch,
                                                    bins_ch,
                                                    bam_files)
   
        // Capture software versions
        CONTIG_DECONTAM.out.version | mix(software_versions_ch) | set{software_versions_ch}
        CONTIG_DECONTAM_HEATMAP.out.version | mix(software_versions_ch) | set{software_versions_ch} 
        RENAME_HEADERS.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SUMMARIZE_ASSEMBLIES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SAM_TO_BAM.out.version | mix(software_versions_ch) | set{software_versions_ch}
        CALL_GENES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        REMOVE_LINEWRAPS.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KO_ANNOTATION.out.version | mix(software_versions_ch) | set{software_versions_ch}
        FILTER_KFAMSCAN.out.version | mix(software_versions_ch) | set{software_versions_ch}
        TAX_CLASSIFICATION.out.version | mix(software_versions_ch) | set{software_versions_ch}
        GET_COV_AND_DET.out.version | mix(software_versions_ch) | set{software_versions_ch}
        MAKE_COMBINED_GENE_LEVEL_TABLES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        MAKE_COMBINED_CONTIG_TAX_TABLES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        METABAT_BINNING.out.version | mix(software_versions_ch) | set{software_versions_ch}
        summarize_bins.out.versions | mix(software_versions_ch) | set{software_versions_ch}
        summarize_mags.out.versions | mix(software_versions_ch) | set{software_versions_ch}


    emit:
        versions = software_versions_ch

}
