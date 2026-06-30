#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Summarize Meta assembled genomes (MAGs) **************************
****************************************************************************************/
include { ZIP_FASTA as ZIP_MAGS } from "./zip_fasta.nf"

//params.min_est_comp = 90
//params.max_est_redund = 10
//params.max_est_strain_het = 50

/*
Scratch directory for gtdb-tk, if wanting to use disk space instead of RAM, can be memory intensive;
see https://ecogenomics.github.io/GTDBTk/faq.html#gtdb-tk-reaches-the-memory-limit-pplacer-crashes
leave empty if wanting to use memory, the default, put in quotes the path to a directory that 
already exists if wanting to use disk space
*/

//params.gtdb_tk_scratch_location = ""

/*
 * ========================================================================================
 * PROCESS: FILTER_CHECKM_RESULTS_AND_COPY_MAGS
 * ========================================================================================
 *
 * SUMMARY:
 *   Filter checkm results to retrieve MAGs
 *
 * INPUTS:
 *   1. path: bins_checkm_results
 *      Cardinality: one
 *      Description: Input file: bins checkm results
 *
 *   2. path: bins
 *      Cardinality: one
 *      Description: Input file: bins
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}MAGs-checkm-out.tsv (emit: MAGs_checkm_out)
 *
 *   2. path: MAGs_dir/ (emit: MAGs_dir)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: CheckM
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: bit
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

/*  Retrieve MAGS.
    Filters checkm results based on estimate completion, redundancy, and 
    strain heterogeneity. Defaults are conservatively 90, 10, and 50  
*/

process FILTER_CHECKM_RESULTS_AND_COPY_MAGS {

    tag "Filtering checkm-s results..."
    label "bit"

    input:
        path(bins_checkm_results)
        path(bins) 
    output:
        path("${params.additional_filename_prefix}MAGs-checkm-out.tsv"), emit: MAGs_checkm_out
        path("MAGs_dir/"), emit: MAGs_dir
    script:
        """        
        mkdir MAGs_dir/
        # Only running if there were bins recovered
        if [ `find -L . -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

            cat <( printf "Bin Id\\tMarker lineage\\t# genomes\\t# markers\\t# marker sets\\t0\\t1\\t2\\t3\\t4\\t5+\\tCompleteness\\tContamination\\tStrain heterogeneity\\n" ) \\
                <( awk -F '\\t' ' \$12 >= ${params.min_est_comp} && \$13 <= ${params.max_est_redund} && \$14 <= ${params.max_est_strain_het} ' ${bins_checkm_results} ) \\
                > MAGs-checkm-out.tmp

            sed 's/-bin\\./-MAG-/' MAGs-checkm-out.tmp > ${params.additional_filename_prefix}MAGs-checkm-out.tsv
            
            for MAG in `cut -f 1 MAGs-checkm-out.tmp | tail -n +2`
            do
                new_ID=`echo \$MAG | sed 's/-bin\\./-MAG-/'`
                cp \$MAG.fasta MAGs_dir/\${new_ID}.fasta
            done

        else

            printf "There were no MAGs recovered.\\n" > ${params.additional_filename_prefix}MAGs-checkm-out.tsv

        fi
        """
}


/*
 * ========================================================================================
 * PROCESS: GET_MAGS
 * ========================================================================================
 *
 * SUMMARY:
 *   Collect list of MAGs
 *
 * INPUTS:
 *   1. path: MAGs_dir
 *      Cardinality: one
 *      Description: Input file: Directory containing MAGs
 *
 * OUTPUTS:
 *   1. path: ${MAGs_dir}/*fasta [OPTIONAL]
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: bit
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process GET_MAGS {

   tag "Collecting MAGs list"
   label "bit"

   input:
       path(MAGs_dir)

    output:
       path("${MAGs_dir}/*fasta"), optional: true

    script:
       " echo collecting mags"

}


/*
 * ========================================================================================
 * PROCESS: GTDBTK_ON_MAG
 * ========================================================================================
 *
 * SUMMARY:
 *   Assign taxonomy to MAG with gtdb-tk
 *
 * INPUTS:
 *   1. path: MAGs_checkm_out
 *      Cardinality: one
 *      Description: Input file: MAGs checkm output
 *
 *   2. path: MAG
 *      Cardinality: one
 *      Description: Input file: MAG
 *
 *   3. path: gtdbtk_db_dir
 *      Cardinality: one
 *      Description: Input file:  path to GTDBTK database. Dummy here ensures dependency on process that creates the database
 *
 *   4. val: use_gtdbtk_scratch_location
 *      Cardinality: one
 *      Description: Parameter value: should a scratch location be use for gtdbtk 
 *
 *   5. env: GTDBTK_DATA_PATH
 *      Cardinality: one
 *      Description: Parameter value: environmental variable holding the path to your GTDBTK database
 *
 * OUTPUTS:
 *   1. path: *.summary.tsv (emit: summary)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: GTDB-Tk
 *   Container: [Defined in config/default.config]
 *   Conda: envs/gtdb-tk.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// Assign taxonomy to MAGs with gtdb-tk
process  GTDBTK_ON_MAG {
   
    tag "Assigning taxonomy to your MAGs with gtdb-tk..." 

    input:
        path(MAGs_checkm_out)
        path(MAG)
        path(gtdbtk_db_dir)
        val(use_gtdbtk_scratch_location)
        env('GTDBTK_DATA_PATH')
           
    output:
        path("*.summary.tsv"), emit: summary
        path("versions.txt"), emit: version
    script:
        """
        MAG_NAME=`basename -s ".fasta" ${MAG}`
        mkdir MAG_dir/ && cp ${MAG} MAG_dir/
        # Only running if any MAGs were recovered
        if [ `find -L MAG_dir -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

            if [ ${use_gtdbtk_scratch_location} == 'true' ]; then

              [ -d gtdbtk_scratch_location/ ] || mkdir gtdbtk_scratch_location/

                gtdbtk classify_wf \\
                        --scratch_dir gtdbtk_scratch_location/ \\
                        --genome_dir MAG_dir/ \\
                        -x fasta \\
                        --out_dir gtdbtk-out/ \\
                        --cpus ${task.cpus} \\
                        --skip_ani_screen 

            else

                gtdbtk classify_wf \\
                       --genome_dir MAG_dir/ \\
                       -x fasta \\
                       --out_dir gtdbtk-out/ \\
                       --cpus ${task.cpus} \\
                       --skip_ani_screen

            fi

        cat gtdbtk-out/*summary.tsv > \${MAG_NAME}.summary.tsv
        else

            mkdir -p gtdbtk-out/
            printf "There were no MAGs recovered.\\n" \\
                   > gtdbtk-out/No-MAGs-recovered.txt
                   
            printf "\\n\\nThere were no MAGs recovered, so GTDB-tk was not run.\\n\\n"


            touch failed.summary.tsv 

        fi
        gtdbtk -h |grep "GTDB-Tk" | sed -E 's/.+\\s+(GTDB-Tk v.+)\\s+.+/\\1/' > versions.txt

        """
}


/*
 * ========================================================================================
 * PROCESS: COMBINE_GTDBTK
 * ========================================================================================
 *
 * SUMMARY:
 *   Combine GTDBTK summaries from all MAGs
 *
 * INPUTS:
 *   1. path: summaries
 *      Cardinality: one
 *      Description: Input file: list of GTDBTK summaries
 *
 * OUTPUTS:
 *   1. path: gtdbtk_summary.tsv (emit: summary)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: GTDB-Tk
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: bit
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process COMBINE_GTDBTK {

    tag "Combining GTDBTK summaries..."
    label "bit"


    input:
        path(summaries)

    output:
        path("gtdbtk_summary.tsv"), emit: summary
  

    script:
       """
       cat ${summaries} > temp_summary.tsv

       (grep "^user_genome" temp_summary.tsv | sort -u; \\
         grep -v "^user_genome" temp_summary.tsv | sort -uV) \\
           >  gtdbtk_summary.tsv
       """
}


/*
 * ========================================================================================
 * PROCESS: SUMMARIZE_MAG_ASSEMBLIES
 * ========================================================================================
 *
 * SUMMARY:
 *   Summarize MAG assemblies
 *
 * INPUTS:
 *   1. path: MAGs_dir
 *      Cardinality: one
 *      Description: Input file: Directory of MAGs
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}MAG-assembly-summaries.tsv (emit: summary)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: bit
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// Summarize MAG assemblies
process  SUMMARIZE_MAG_ASSEMBLIES {

    tag "Summarizing MAG assemblies..."
    label "bit"

    input:
        path(MAGs_dir)
    output:
        path("${params.additional_filename_prefix}MAG-assembly-summaries.tsv"), emit: summary
        path("versions.txt"), emit: version
    script:
        """
        # Only running if any MAGs were recovered
        if [ `find -L ${MAGs_dir} -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then
            
            # Remove fasta index if already exists
            rm -rf ${MAGs_dir}/*.fxi
            bit-summarize-assembly ${MAGs_dir}/*.fasta -o MAG-summaries.tmp -t

            # Slimming down the output
            cut -f 1,2,3,5,6,8,11,18,19,20 MAG-summaries.tmp \\
               > ${params.additional_filename_prefix}MAG-assembly-summaries.tsv

        else

            printf "There were no MAGs recovered.\\n" \\
               > ${params.additional_filename_prefix}MAG-assembly-summaries.tsv

        fi
        bit-version |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt
        """
}

/*
 * ========================================================================================
 * PROCESS: GENERATE_MAGS_OVERVIEW_TABLE
 * ========================================================================================
 *
 * SUMMARY:
 *   Generate an overview table of all MAGs
 *
 * INPUTS:
 *   1. path: MAG_assembly_summaries
 *      Cardinality: one
 *      Description: Input file: MAG assembly summaries
 *
 *   2. path: MAGs_checkm_out
 *      Cardinality: one
 *      Description: Input file: MAGs checkm output
 *
 *   3. path: gtdbtk_summary
 *      Cardinality: one
 *      Description: Input file: gtdbtk summary
 *
 *   4. path: MAGs_dir
 *      Cardinality: one
 *      Description: Input file: directory of MAGs
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}MAGs-overview${params.assay_suffix}.tsv
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: mags, bit
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process  GENERATE_MAGS_OVERVIEW_TABLE {

    tag "Generating an overview table of all MAGs..."
    label "mags"
    label "bit"

    input:
        path(MAG_assembly_summaries)
        path(MAGs_checkm_out)
        path(gtdbtk_summary)
        path(MAGs_dir)
    output:
        path("${params.additional_filename_prefix}MAGs-overview${params.assay_suffix}.tsv")

    script:
        """
        # Only running if any MAGs were recovered
        if [ `find -L ${MAGs_dir}  -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

        #--------------------------- get_MAGs_estimates_and_taxonomy.sh ------------------------------------#
        get_MAGs_estimates_and_taxonomy.sh ${MAGs_dir} ${MAG_assembly_summaries} ${MAGs_checkm_out} ${gtdbtk_summary}
        #----------------------------------------------------------------------------------------------------#

            # Adding headers
            cat <(printf "est. completeness\\test. redundancy\\test. strain heterogeneity\\n") \\
                checkm-estimates.tmp > checkm-estimates-with-headers.tmp

            cat <(printf "domain\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\n") \\
                gtdb-taxonomies.tmp > gtdb-taxonomies-with-headers.tmp

            paste ${MAG_assembly_summaries} \\
                  checkm-estimates-with-headers.tmp \\
                  gtdb-taxonomies-with-headers.tmp \\
                  > MAGs-overview.tmp

            # Ordering by taxonomy
            head -n 1 MAGs-overview.tmp > MAGs-overview-header.tmp

            tail -n +2 MAGs-overview.tmp | \\
            sort -t \$'\\t' -k 14,20 > MAGs-overview-sorted.tmp

            cat MAGs-overview-header.tmp MAGs-overview-sorted.tmp \\
               > ${params.additional_filename_prefix}MAGs-overview${params.assay_suffix}.tsv

        else

            printf "There were no MAGs recovered.\\n" \\
              > ${params.additional_filename_prefix}MAGs-overview${params.assay_suffix}.tsv

        fi
        """
}

/*
 * ========================================================================================
 * PROCESS: SUMMARIZE_MAG_LEVEL_KO_ANNOTATIONS
 * ========================================================================================
 *
 * SUMMARY:
 *   Parse MAG KO annotations
 *
 * INPUTS:
 *   1. path: MAGs_overview
 *      Cardinality: one
 *      Description: Input file: MAGs overview
 *
 *   2. path: gene_coverage_annotation_and_tax_files
 *      Cardinality: one
 *      Description: Input file: gene coverage annotation and taxonomy files
 *
 *   3. path: MAGs_dir
 *      Cardinality: one
 *      Description: Input file: directory of MAGs
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}MAG-level-KO-annotations${params.assay_suffix}.tsv (emit: summary)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: mags, bit
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process SUMMARIZE_MAG_LEVEL_KO_ANNOTATIONS {

    tag "Parsing MAG KO annotations..."
    label "mags"
    label "bit"

    input:
        path(MAGs_overview)
        path(gene_coverage_annotation_and_tax_files)
        path(MAGs_dir)
    output:
        path("${params.additional_filename_prefix}MAG-level-KO-annotations${params.assay_suffix}.tsv"), emit: summary
        path("versions.txt"), emit: version 
    script:
        """
        # Only running if any MAGs were recovered
        if [ `find -L ${MAGs_dir} -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

            for MAG in `cut -f 1 ${MAGs_overview} | tail -n +2`
            do

                sample_ID=`echo \$MAG | sed 's/-MAG-[0-9]*\$//'`
                grep "^>" ${MAGs_dir}/\$MAG.fasta | tr -d ">" > curr-contig-ids.tmp

                parse-MAG-annots.py \\
                        -i \${sample_ID}-gene-coverage-annotation-and-tax.tsv \\
                        -w curr-contig-ids.tmp \\
                        -M \$MAG \\
                        -o ${params.additional_filename_prefix}MAG-level-KO-annotations${params.assay_suffix}.tsv

            done

        else

            printf "There were no MAGs recovered.\\n" \\
               > ${params.additional_filename_prefix}MAG-level-KO-annotations${params.assay_suffix}.tsv

        fi
        python --version > versions.txt
        """
}

/*
 * ========================================================================================
 * PROCESS: SUMMARIZE_MAG_KO_ANNOTS_WITH_KEGG_DECODER
 * ========================================================================================
 *
 * SUMMARY:
 *   Summarize MAG KO annotations with kegg decoder
 *
 * INPUTS:
 *   1. path: MAG_level_KO_annotations
 *      Cardinality: one
 *      Description: Input file: MAG level KO annotations
 *
 *   2. path: MAGs_dir
 *      Cardinality: one
 *      Description: Input file: directory of MAGs
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}MAG-KEGG-Decoder-out${params.assay_suffix}.* (emit: summary)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/keggdecoder.yaml
 *   Labels: mags
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process SUMMARIZE_MAG_KO_ANNOTS_WITH_KEGG_DECODER {

    tag "Summarizing MAG KO annotations using kegg decoder..."
    label "mags"


    input:
        path(MAG_level_KO_annotations)
        path(MAGs_dir)
    output:
        path("${params.additional_filename_prefix}MAG-KEGG-Decoder-out${params.assay_suffix}.*"), emit: summary
        path("versions.txt"), emit: version
    script:
        """
        # Getting number of MAGs recovered
        num_mags_recovered=`find -L ${MAGs_dir}/ -name '*.fasta' | wc -l | sed 's/^ *//'`
        # Only running if any MAGs were recovered
        if [ \$num_mags_recovered -gt 0 ]; then

            # KEGGDecoder splits on the first underscore to identify unique genome/MAG IDs
            # this can be problematic with how things are named, so we are swapping them all to not have
            # any "_" first, then afterwards we are changing the output table back to the original names so 
            # they match elsewhere (they will still be slightly different in the html output, but that is
            # only manually explored anyway)

            # Making version of input for KEGGDecoder with no underscores
            tr "_" "-" < ${MAG_level_KO_annotations} > mod-MAG-level-KO-annotations.tmp

            # Making mapping file
            paste <( cut -f 1 ${MAG_level_KO_annotations} ) \\
                  <( cut -f 1 mod-MAG-level-KO-annotations.tmp ) \\
                  > MAG-ID-map.tmp

            # Running KEGGDecoder
            # can only create html output if there are more than 1
            if [ \$num_mags_recovered -gt 1 ]; then
                KEGG-decoder -v interactive -i mod-MAG-level-KO-annotations.tmp -o MAG-KEGG-Decoder-out.tmp
                
                ## adding additional prefix to html output if there is one
               if [ -f MAG-KEGG-Decoder-out.html ]; then


                [ -f ${params.additional_filename_prefix}MAG-KEGG-Decoder-out${params.assay_suffix}.html ] || \\
                mv MAG-KEGG-Decoder-out.html  ${params.additional_filename_prefix}MAG-KEGG-Decoder-out${params.assay_suffix}.html
                
               fi

            else

                KEGG-decoder -i mod-MAG-level-KO-annotations.tmp -o MAG-KEGG-Decoder-out.tmp

            fi

            # Swapping MAG IDs back in output tsv from KEGGDecoder
            swap-MAG-IDs.py -i MAG-KEGG-Decoder-out.tmp -m MAG-ID-map.tmp -o MAG-KEGG-Decoder-out.tsv && \\
            
            [ -f ${params.additional_filename_prefix}MAG-KEGG-Decoder-out${params.assay_suffix}.tsv ] || \\
            mv MAG-KEGG-Decoder-out.tsv \\
              ${params.additional_filename_prefix}MAG-KEGG-Decoder-out${params.assay_suffix}.tsv


        else

            printf "There were no MAGs recovered.\\n" \\
               >  ${params.additional_filename_prefix}MAG-KEGG-Decoder-out${params.assay_suffix}.tsv

        fi
        python --version > versions.txt
        """
}



workflow summarize_mags {
    take:
        bins_checkm_results_ch
        bins_ch
        gtdbtk_db_dir
        use_gtdbtk_scratch_location
        gene_coverage_annotation_and_tax_files_ch 


    main:
        FILTER_CHECKM_RESULTS_AND_COPY_MAGS(bins_checkm_results_ch, bins_ch) 
        MAGs_checkm_out_ch = FILTER_CHECKM_RESULTS_AND_COPY_MAGS.out.MAGs_checkm_out
        MAGs_dir_ch = FILTER_CHECKM_RESULTS_AND_COPY_MAGS.out.MAGs_dir
        ZIP_MAGS(Channel.of("MAG"), MAGs_dir_ch)

        
        MAGs_ch = GET_MAGS(MAGs_dir_ch).flatten()
        GTDBTK_ON_MAG(MAGs_checkm_out_ch, MAGs_ch, gtdbtk_db_dir, use_gtdbtk_scratch_location, gtdbtk_db_dir)
        gtdbtk_summaries_ch = GTDBTK_ON_MAG.out.summary.collect()

        COMBINE_GTDBTK(gtdbtk_summaries_ch)

        SUMMARIZE_MAG_ASSEMBLIES(MAGs_dir_ch)
        MAG_assembly_summaries_ch = SUMMARIZE_MAG_ASSEMBLIES.out.summary

        MAGs_overview_ch = GENERATE_MAGS_OVERVIEW_TABLE(MAG_assembly_summaries_ch,
                                                        MAGs_checkm_out_ch,
                                                        COMBINE_GTDBTK.out.summary,
                                                        MAGs_dir_ch)

        SUMMARIZE_MAG_LEVEL_KO_ANNOTATIONS(MAGs_overview_ch, 
                                           gene_coverage_annotation_and_tax_files_ch, 
                                           MAGs_dir_ch)
        MAG_level_KO_annotations_ch = SUMMARIZE_MAG_LEVEL_KO_ANNOTATIONS.out.summary

        SUMMARIZE_MAG_KO_ANNOTS_WITH_KEGG_DECODER(MAG_level_KO_annotations_ch, MAGs_dir_ch)

        // Capture software versions
        software_versions_ch = Channel.empty()
        ZIP_MAGS.out.version | mix(software_versions_ch) | set{software_versions_ch}
        GTDBTK_ON_MAG.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SUMMARIZE_MAG_ASSEMBLIES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SUMMARIZE_MAG_LEVEL_KO_ANNOTATIONS.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SUMMARIZE_MAG_KO_ANNOTS_WITH_KEGG_DECODER.out.version | mix(software_versions_ch) | set{software_versions_ch}

    emit:
        MAGs_overview = MAGs_overview_ch
        MAGs_dir = MAGs_dir_ch
        versions = software_versions_ch
    
}
