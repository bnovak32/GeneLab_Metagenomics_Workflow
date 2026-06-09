#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
**************************  Sequence Assembly Annotation *******************************
****************************************************************************************/

/*
 * ========================================================================================
 * PROCESS: CALL_GENES
 * ========================================================================================
 *
 * SUMMARY:
 *   Predict genes in sample assembly with prodigal
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(assembly)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - assembly: path to sample assembly/contigs
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}-genes.faa"), path("${sample_id}-genes.fasta"), path("${sample_id}-genes${params.assay_suffix}.gff") (emit: genes)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: Prodigal
 *   Container: [Defined in config/default.config]
 *   Conda: envs/prodigal.yaml
 *   Labels: call_genes
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// This process calls genes on each assembly file.
process CALL_GENES {

    tag "Predicting genes for ${sample_id}-s assembly"
    label "call_genes"
    
    input:
        tuple val(sample_id), path(assembly) 
    output:
        // Amino acids, nucleotides and gff
        tuple val(sample_id), path("${sample_id}-genes.faa"), path("${sample_id}-genes.fasta"), path("${sample_id}-genes${params.assay_suffix}.gff"), emit: genes
        path("versions.txt"), emit: version
    script:
        """
        # Only running if assembly produced any contigs
        if [ -s ${assembly} ]; then

            prodigal -q -c -p meta -a ${sample_id}-genes.faa \\
                     -d ${sample_id}-genes.fasta \\
                     -f gff -o ${sample_id}-genes${params.assay_suffix}.gff \\
                     -i ${assembly} 
        else

            touch ${sample_id}-genes.faa ${sample_id}-genes.fasta ${sample_id}-genes${params.assay_suffix}.gff
            printf "Gene-calling not performed because the assembly didn't produce anything.\\n"

        fi
        prodigal -v 2>&1 | grep Prodigal > versions.txt
        """
}        

/*
 * ========================================================================================
 * PROCESS: REMOVE_LINEWRAPS
 * ========================================================================================
 *
 * SUMMARY:
 *   Remove line wraps in sample nucleotide and amino acid files
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(aa), path(nt), path(gff)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - aa: path to sample amino acids file
 *                 - nt: path to sample nucleotide file
 *                 - gff: path to sample gene feature file
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}-genes${params.assay_suffix}.faa"), path("${sample_id}-genes${params.assay_suffix}.fasta") (emit: genes)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: call_genes, bit
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// Removing line-wraps using bit
process REMOVE_LINEWRAPS {

    tag "Remove line wraps in ${sample_id}-s nucleotide and amino acid file..."
    label "call_genes"
    label "bit"

    input:
        tuple val(sample_id), path(aa), path(nt), path(gff)
   
    output:
        tuple val(sample_id), path("${sample_id}-genes${params.assay_suffix}.faa"), path("${sample_id}-genes${params.assay_suffix}.fasta"), emit: genes
        path("versions.txt"), emit: version
    script:
        """
         if [ -s ${aa} ] && [ -s ${nt} ]; then
            # Removing line-wraps
            bit-remove-wraps ${aa} > ${sample_id}-genes.faa.tmp 2> /dev/null && \\
            mv ${sample_id}-genes.faa.tmp ${sample_id}-genes${params.assay_suffix}.faa
            
            bit-remove-wraps ${nt} > ${sample_id}-genes.fasta.tmp 2> /dev/null && \\
            mv ${sample_id}-genes.fasta.tmp ${sample_id}-genes${params.assay_suffix}.fasta
        else

            touch ${sample_id}-genes${params.assay_suffix}.faa ${sample_id}-genes${params.assay_suffix}.fasta
            printf "Line wrapping not performed because gene-calling wasn't performed on ${sample_id}.\\n"
        fi 
        bit-version |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt
        """
}


/*
 * ========================================================================================
 * PROCESS: KO_ANNOTATION
 * ========================================================================================
 *
 * SUMMARY:
 *   KO annotation of sample predicted amino acids with kofamscan.
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(assembly), path(aa), path(nt)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - assembly: path to sample assembly/contigs
 *                 - aa: path to sample amino acids file
 *                 - nt: path to sample nucleotide file
 *
 *   2. path: ko_db_dir
 *      Cardinality: one
 *      Description: Input file: KO Database directory
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}-KO-tab.tmp") (emit: temp_table)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/kofamscan.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// This process runs the gene-level (KO) functional annotation for each sample.
// KO annotation of the predicted amino acids
process KO_ANNOTATION {

    tag "Running KO annotation of ${sample_id}-s predicted amino acids.."
    
    input:
       tuple val(sample_id), path(assembly), path(aa), path(nt)
       path(ko_db_dir)
    output:
        tuple val(sample_id), path("${sample_id}-KO-tab.tmp"), emit: temp_table
        path("versions.txt"), emit: version
    script:
        """
        # only running if assembly produced any contigs and genes were identified (they are required for this)
        if [ -s ${assembly} ] && [ -s ${aa} ]; then

            exec_annotation -p ${ko_db_dir}/profiles/ \\
                            -k ${ko_db_dir}/ko_list \\
                            --cpu ${task.cpus} -f detail-tsv \\
                            -o ${sample_id}-KO-tab.tmp --tmp-dir ${sample_id}-tmp-KO-dir \\
                            --report-unannotated ${aa}

        else

            touch ${sample_id}-KO-tab.tmp
            printf "Functional annotations not performed because the assembly didn't produce anything and/or no genes were identified.\\n"

        fi
        exec_annotation -v > versions.txt
        """
}


/*
 * ========================================================================================
 * PROCESS: FILTER_KFAMSCAN
 * ========================================================================================
 *
 * SUMMARY:
 *   Filter sample KO annotation results
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(KO_tab_tmp)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - KO_tab_tmp: path to temporary KO annotation table produced by the KO_ANNOTATION process.
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}-annotations.tsv") (emit: ko_annotation)
 *
 *   2. path: versions.txt (emit: version)
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

process FILTER_KFAMSCAN {

    tag "Filtering ${sample_id}-s KO annotation results..."
    label "bit"
    label "contig_annotation"

    input:
       tuple val(sample_id), path(KO_tab_tmp)   
    output:
        tuple val(sample_id), path("${sample_id}-annotations.tsv"), emit: ko_annotation
        path("versions.txt"), emit: version
    script:
        """
        if [ -s ${KO_tab_tmp} ]; then

            bit-filter-KOFamScan-results -i ${KO_tab_tmp} -o ${sample_id}-annotations.tsv

        else 

            touch ${sample_id}-annotations.tsv
            printf "Nothing to filter since functional annotation was not performed.\\n"

        fi
        bit-version |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt
        """

}

/*
 * ========================================================================================
 * PROCESS: TAX_CLASSIFICATION
 * ========================================================================================
 *
 * SUMMARY:
 *   Taxonomy classification of sample contigs with CAT (Contig Annotation Tool) 
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(assembly), path(aa), path(nt)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - assembly: path to sample assembly/contigs
 *                 - aa: path to sample amino acids file
 *                 - nt: path to sample nucleotide file
 *
 *   2. path: cat_db
 *      Cardinality: one
 *      Description: Input file: CAT database directory
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}-gene-tax.tsv"), path("${sample_id}-contig-tax.tsv") (emit: taxonomy)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: CAT
 *   Container: [Defined in config/default.config]
 *   Conda: envs/cat.yaml
 *   Labels: contig_annotation
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// This process runs the gene- and contig-level taxonomic classifications for each assembly.
process TAX_CLASSIFICATION {

    tag "Taxonomy classification of ${sample_id}-s "
    label "contig_annotation"

    input:
       tuple val(sample_id), path(assembly), path(aa), path(nt)
       path(cat_db)
    output:
        // Gene and contig taxonomy
        tuple val(sample_id), path("${sample_id}-gene-tax.tsv"), path("${sample_id}-contig-tax.tsv"), emit: taxonomy
        path("versions.txt"), emit: version
    script:
        """
        # Only running if assembly produced any contigs and 
        # genes were identified (they are required for this)
        if [ -s ${assembly} ] && [ -s ${aa} ]; then

            CAT contigs -d ${cat_db}/${params.cat_db_sub_dir} -t ${cat_db}/${params.cat_taxonomy_dir} \\
                        -n ${task.cpus} -r 3 --top 4 \\
                        --I_know_what_Im_doing -c ${assembly} \\
                        -p ${aa} -o ${sample_id}-tax-out.tmp \\
                        --no_stars --block_size ${params.block_size} \\
                        --index_chunks 2 --force 

            # Adding names to gene classifications
            CAT add_names -i ${sample_id}-tax-out.tmp.ORF2LCA.txt \\
                          -o ${sample_id}-gene-tax.tmp -t ${cat_db}/${params.cat_taxonomy_dir} \\
                          --only_official --exclude_scores

            # Formatting gene classifications
            bash format-gene-tax-classifications.sh \\
                       ${sample_id}-gene-tax.tmp ${sample_id}-gene-tax.tsv

            # Adding names to contig classifications
            CAT add_names -i ${sample_id}-tax-out.tmp.contig2classification.txt \\
                          -o ${sample_id}-contig-tax.tmp -t ${cat_db}/${params.cat_taxonomy_dir} \\
                          --only_official --exclude_scores

            # Formatting contig classifications
            bash format-contig-tax-classifications.sh \\
                  ${sample_id}-contig-tax.tmp ${sample_id}-contig-tax.tsv

        else

            touch ${sample_id}-gene-tax.tsv ${sample_id}-contig-tax.tsv
            printf "Assembly-based taxonomic classification not performed because the assembly didn't produce anything and/or no genes were identified.\\n" 

        fi
        CAT --version | sed -E 's/(CAT v.+)\\s\\(.+/\\1/'  > versions.txt
        """
}

workflow annotate_assembly {
    take:
        assembly_ch
        ko_db_dir
        cat_db
 
    main:
        CALL_GENES(assembly_ch)
        genes_ch = CALL_GENES.out.genes | REMOVE_LINEWRAPS.out.genes
        KO_ANNOTATION(assembly_ch.join(genes_ch), ko_db_dir)
        KO_ANNOTATION.out.temp_table | FILTER_KFAMSCAN.out.ko_annotation
        TAX_CLASSIFICATION(assembly_ch.join(genes_ch), cat_db)

}
