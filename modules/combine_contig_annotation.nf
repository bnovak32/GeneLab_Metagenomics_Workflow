#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


/**************************************************************************************** 
**************************  Combine Contig Annotation ***********************************
****************************************************************************************/

/*
 * ========================================================================================
 * PROCESS: COMBINE_GENE_ANNOTS_TAX_AND_COVERAGE
 * ========================================================================================
 *
 * SUMMARY:
 *   Combine gene and taxonomy annotations for sample
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(gene_coverages), path(contig_coverages), path(annotations),
 *                   path(gene_tax), path(contig_tax), path(aa), path(nt), path(assembly)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - gene_coverages: path to sample gene coverages
 *                 - contig_coverages: path to sample contig coverages
 *                 - gene_tax: path to sample gene taxonomy
 *                 - contig_tax: path to sample contig taxonomy
 *                 - aa: path to sample amino acids file
 *                 - nt: path to sample nucleotide file
 *                 - assembly: path to sample assembly/contigs
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}-gene-coverage-annotation-and-tax${params.assay_suffix}.tsv")
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: bit, contig_annotation
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

/*
This process combines the gene-level functional annotations, taxonomic classifications,
 and coverage information for each individual sample.
*/
process COMBINE_GENE_ANNOTS_TAX_AND_COVERAGE {

    tag "Combining gene and taxonomy annotations for ${sample_id}"
    label "bit"
    label "contig_annotation"

    input:
        tuple val(sample_id), path(gene_coverages), path(contig_coverages),
               path(annotations), path(gene_tax), path(contig_tax),
               path(aa), path(nt), path(assembly)
    output:
        tuple val(sample_id), path("${sample_id}-gene-coverage-annotation-and-tax${params.assay_suffix}.tsv")
    script:
        """
        # Only running if the assembly produced anything 
        # and genes were identified (they are required for this)
        if [ -s ${assembly} ] && [ -s ${aa} ]; then

            paste <( tail -n +2 ${gene_coverages} | \\
            sort -V -k 1 ) <( tail -n +2 ${annotations} | \\
            sort -V -k 1 | cut -f 2- ) <( tail -n +2 ${gene_tax} | \\
            sort -V -k 1 | \\
            cut -f 2- ) > ${sample_id}-gene.tmp

            paste <( head -n 1 ${gene_coverages} ) \\
                  <( head -n 1  ${annotations} | \\
                  cut -f 2- ) <( head -n 1 ${gene_tax} | \\
                  cut -f 2- ) > ${sample_id}-gene-header.tmp

            cat ${sample_id}-gene-header.tmp ${sample_id}-gene.tmp \\
                > ${sample_id}-gene-coverage-annotation-and-tax${params.assay_suffix}.tsv

            rm -rf ${sample_id}-gene.tmp ${sample_id}-gene-header.tmp
           

        else

            printf "gene_ID\\tcoverage\\tKO_ID\\tKO_function\\ttaxid\\tdomain\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\n" \\
            > ${sample_id}-gene-coverage-annotation-and-tax${params.assay_suffix}.tsv
            
        fi
        """
}


/*
 * ========================================================================================
 * PROCESS: MAKE_COMBINED_GENE_LEVEL_TABLES
 * ========================================================================================
 *
 * SUMMARY:
 *   Combine all gene level annotations
 *
 * INPUTS:
 *   1. path: gene_coverage_annotation_and_tax_files
 *      Cardinality: one
 *      Description: Input file: gene coverage annotation and taxonomy files
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}Combined-gene-level-KO-function-coverages${params.assay_suffix}.tsv (emit: raw_function_coverages)
 *
 *   2. path: ${params.additional_filename_prefix}Combined-gene-level-KO-function-coverages-CPM${params.assay_suffix}.tsv (emit: norm_function_coverages)
 *
 *   3. path: ${params.additional_filename_prefix}Combined-gene-level-taxonomy-coverages${params.assay_suffix}.tsv (emit: raw_taxonomy_coverages)
 *
 *   4. path: ${params.additional_filename_prefix}Combined-gene-level-taxonomy-coverages-CPM${params.assay_suffix}.tsv (emit: norm_taxonomy_coverages)
 *
 *   5. path: versions.txt (emit: version)
 *
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


process MAKE_COMBINED_GENE_LEVEL_TABLES {

    tag "Combining all gene level annotations...."
    label "bit"

    input:
        path(gene_coverage_annotation_and_tax_files)
    output:
        path("${params.additional_filename_prefix}Combined-gene-level-KO-function-coverages${params.assay_suffix}.tsv"), emit: raw_function_coverages
        path("${params.additional_filename_prefix}Combined-gene-level-KO-function-coverages-CPM${params.assay_suffix}.tsv"), emit: norm_function_coverages
        path("${params.additional_filename_prefix}Combined-gene-level-taxonomy-coverages${params.assay_suffix}.tsv"), emit: raw_taxonomy_coverages
        path("${params.additional_filename_prefix}Combined-gene-level-taxonomy-coverages-CPM${params.assay_suffix}.tsv"), emit: norm_taxonomy_coverages
        path("versions.txt"), emit: version
    script:
        """
        bit-GL-combine-KO-and-tax-tables ${gene_coverage_annotation_and_tax_files} -o ${params.additional_filename_prefix}Combined

        # Renaming to have GL assay-specific suffix
        [ -f "${params.additional_filename_prefix}Combined-gene-level-KO-function-coverages${params.assay_suffix}.tsv" ] || \\
        mv "${params.additional_filename_prefix}Combined-gene-level-KO-function-coverages.tsv" \\
           "${params.additional_filename_prefix}Combined-gene-level-KO-function-coverages${params.assay_suffix}.tsv"

        [ -f "${params.additional_filename_prefix}Combined-gene-level-KO-function-coverages-CPM${params.assay_suffix}.tsv" ] || \\
        mv "${params.additional_filename_prefix}Combined-gene-level-KO-function-coverages-CPM.tsv" \\
           "${params.additional_filename_prefix}Combined-gene-level-KO-function-coverages-CPM${params.assay_suffix}.tsv"

        [ -f "${params.additional_filename_prefix}Combined-gene-level-taxonomy-coverages${params.assay_suffix}.tsv" ] || \\
        mv "${params.additional_filename_prefix}Combined-gene-level-taxonomy-coverages.tsv" \\
           "${params.additional_filename_prefix}Combined-gene-level-taxonomy-coverages${params.assay_suffix}.tsv"

        [ -f "${params.additional_filename_prefix}Combined-gene-level-taxonomy-coverages-CPM${params.assay_suffix}.tsv" ] || \\
        mv "${params.additional_filename_prefix}Combined-gene-level-taxonomy-coverages-CPM.tsv" \\
           "${params.additional_filename_prefix}Combined-gene-level-taxonomy-coverages-CPM${params.assay_suffix}.tsv"
        bit-version |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt
        """
}


/*
 * ========================================================================================
 * PROCESS: COMBINE_CONTIG_TAX_AND_COVERAGE
 * ========================================================================================
 *
 * SUMMARY:
 *   Combine taxonomy and coverage for sample
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(gene_coverages), path(contig_coverages),
 *                   path(gene_tax), path(contig_tax), path(aa), path(nt), path(assembly)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - gene_coverages: path to sample gene coverages
 *                 - contig_coverages: path to sample contig coverages
 *                 - gene_tax: path to sample gene taxonomy
 *                 - contig_tax: path to sample contig taxonomy
 *                 - aa: path to sample amino acids file
 *                 - nt: path to sample nucleotide file
 *                 - assembly: path to sample assembly/contigs
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}-contig-coverage-and-tax${params.assay_suffix}.tsv")
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: bit, contig_annotation
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// This process combines the contig-level taxonomic and 
// coverage information for each individual sample. 
process COMBINE_CONTIG_TAX_AND_COVERAGE {

    tag "Combining taxonomy and coverage for ${sample_id}...."
    label "bit"
    label "contig_annotation"

    input:
        tuple val(sample_id), path(gene_coverages), path(contig_coverages),
               path(gene_tax), path(contig_tax),
               path(aa), path(nt), path(assembly)        
    output:
        tuple val(sample_id), path("${sample_id}-contig-coverage-and-tax${params.assay_suffix}.tsv")
    script:
        """
        # Only running if the assembly produced anything
        if [ -s ${assembly} ]; then

            # If there were no genes called, there is no contig-level taxonomy, so dealing with that here
            if [ -s ${aa} ]; then

                paste <( tail -n +2 ${contig_coverages} | \\
                sort -V -k 1 ) <( tail -n +2 ${contig_tax} | \\
                sort -V -k 1 | cut -f 2- ) > ${sample_id}-contig.tmp
                paste <( head -n 1 ${contig_coverages} ) \\
                      <( head -n 1 ${contig_tax} | cut -f 2- ) \\
                      > ${sample_id}-contig-header.tmp

                cat ${sample_id}-contig-header.tmp ${sample_id}-contig.tmp \\
                > ${sample_id}-contig-coverage-and-tax${params.assay_suffix}.tsv

                rm -rf ${sample_id}-contig.tmp ${sample_id}-contig-header.tmp

            else

                paste <( tail -n +2 ${contig_coverages} | sort -V -k 1 ) > ${sample_id}-contig-p1.tmp
                
                sed 's/.*/NA/g' ${sample_id}-contig-p1.tmp > ${sample_id}-tax-col.tmp

                paste ${sample_id}-contig-p1.tmp ${sample_id}-tax-col.tmp ${sample_id}-tax-col.tmp \\
                      ${sample_id}-tax-col.tmp ${sample_id}-tax-col.tmp ${sample_id}-tax-col.tmp \\
                      ${sample_id}-tax-col.tmp ${sample_id}-tax-col.tmp ${sample_id}-tax-col.tmp \\
                      > ${sample_id}-contig.tmp

                cat <( printf "contig_ID\\tcoverage\\ttaxid\\tdomain\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\n" ) \\
                    ${sample_id}-contig.tmp > ${sample_id}-contig-coverage-and-tax${params.assay_suffix}.tsv
                rm -rf ${sample_id}-contig-p1.tmp ${sample_id}-tax-col.tmp ${sample_id}-contig.tmp

            fi

        else

            printf "contig_ID\\tcoverage\\ttaxid\\tdomain\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\n" \\
              > ${sample_id}-contig-coverage-and-tax${params.assay_suffix}.tsv

        fi
        """
}

/*
 * ========================================================================================
 * PROCESS: MAKE_COMBINED_CONTIG_TAX_TABLES
 * ========================================================================================
 *
 * SUMMARY:
 *   Make summary contig taxonomy table
 *
 * INPUTS:
 *   1. path: contig_coverage_and_tax_files
 *      Cardinality: one
 *      Description: Input file: contig coverage and tax files
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}Combined-contig-level-taxonomy-coverages${params.assay_suffix}.tsv (emit: raw_taxonomy)
 *
 *   2. path: ${params.additional_filename_prefix}Combined-contig-level-taxonomy-coverages-CPM${params.assay_suffix}.tsv (emit: norm_taxonomy)
 *
 *   3. path: versions.txt (emit: version)
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

process MAKE_COMBINED_CONTIG_TAX_TABLES {

    tag "Making a summary contig taxonomy table...."
    label "bit"

    input:
        path(contig_coverage_and_tax_files)
    output:
        path("${params.additional_filename_prefix}Combined-contig-level-taxonomy-coverages${params.assay_suffix}.tsv"), emit: raw_taxonomy
        path("${params.additional_filename_prefix}Combined-contig-level-taxonomy-coverages-CPM${params.assay_suffix}.tsv"), emit: norm_taxonomy
        path("versions.txt"), emit: version
    script:
        """
        bit-GL-combine-contig-tax-tables ${contig_coverage_and_tax_files} -o ${params.additional_filename_prefix}Combined

        # Renaming to have GL assay-specific suffix
        [ -f "${params.additional_filename_prefix}Combined-contig-level-taxonomy-coverages${params.assay_suffix}.tsv" ] || \\
        mv "${params.additional_filename_prefix}Combined-contig-level-taxonomy-coverages.tsv" \\
           "${params.additional_filename_prefix}Combined-contig-level-taxonomy-coverages${params.assay_suffix}.tsv"

        [ -f "${params.additional_filename_prefix}Combined-contig-level-taxonomy-coverages-CPM${params.assay_suffix}.tsv" ] || \\
        mv "${params.additional_filename_prefix}Combined-contig-level-taxonomy-coverages-CPM.tsv" \\
           "${params.additional_filename_prefix}Combined-contig-level-taxonomy-coverages-CPM${params.assay_suffix}.tsv"
        bit-version |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt
        """
}




workflow combine_contig_coverage {
    take:
        coverages_ch 
        annotations_ch 
        taxonomy_ch
        gene_call_ch 
        assembly_ch

    main:
        tax_and_cov_ch = COMBINE_GENE_ANNOTS_TAX_AND_COVERAGE(coverages_ch
                                                                 .join(annotations_ch) 
                                                                 .join(taxonomy_ch)
                                                                 .join(gene_call_ch)
                                                                .join(assembly_ch))

        MAKE_COMBINED_GENE_LEVEL_TABLES(tax_and_cov_ch.map{sample_id, coverage -> file("${coverage}")}.collect())

        combined_cov_ch = COMBINE_CONTIG_TAX_AND_COVERAGE(coverages_ch
                                                             .join(taxonomy_ch) 
                                                             .join(gene_call_ch)
                                                             .join(assembly_ch))

        MAKE_COMBINED_CONTIG_TAX_TABLES(combined_cov_ch.map{sample_id, coverage -> file("${coverage}")}.collect())
}

