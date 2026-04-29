#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * ========================================================================================
 * PROCESS: CLEAN_FASTQC_PATHS
 * ========================================================================================
 *
 * SUMMARY:
 *   Purge genelab paths from MultiQC zip files
 *
 * INPUTS:
 *   1. path: FastQC_Outputs_dir
 *      Cardinality: one
 *      Description: Input file: FastQC Outputs dir
 *
 * OUTPUTS:
 *   1. path: ${OUT_DIR} (emit: clean_dir)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/post_processing.config]
 *   Conda: envs/genelab.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process CLEAN_FASTQC_PATHS {
    tag "Purging genelab paths from MultiQC zip files in ${params.FastQC_Outputs}"
    input:
        path(FastQC_Outputs_dir)
    output:
        path("${OUT_DIR}"), emit: clean_dir
    script:
        OUT_DIR = "${FastQC_Outputs_dir.baseName}"
        """
        WORKDIR=`pwd`
        mv ${FastQC_Outputs_dir} FastQC_Outputs_dir

        [ -d ${OUT_DIR}/ ] || mkdir  ${OUT_DIR}/ && \\
        cp -r FastQC_Outputs_dir/*  ${OUT_DIR}/
        
        [ -f ${OUT_DIR}/versions.txt ] && rm -rf ${OUT_DIR}/versions.txt

        cat `which clean-paths.sh` > \${WORKDIR}/clean-paths.sh
        chmod +x \${WORKDIR}/clean-paths.sh

        echo "Purging paths from multiqc outputs"
        cd \${WORKDIR}/${OUT_DIR}/
        echo "Cleaning raw multiqc files with path info"
        unzip raw_multiqc${params.assay_suffix}_report.zip && rm raw_multiqc${params.assay_suffix}_report.zip
        cd raw_multiqc_report/raw_multiqc_data/

        # No reason not to just run it on all
        echo "Purging paths in all raw QC files..."
        find . -type f -exec bash \${WORKDIR}/clean-paths.sh '{}' ${params.baseDir} \\;
        cd \${WORKDIR}/${OUT_DIR}/

        echo "Re-zipping up raw multiqc"
        zip -r raw_multiqc${params.assay_suffix}_report.zip raw_multiqc_report/ && rm -rf raw_multiqc_report/

        echo "Cleaning filtered multiqc files with path info..."
        unzip filtered_multiqc${params.assay_suffix}_report.zip && rm filtered_multiqc${params.assay_suffix}_report.zip
        cd filtered_multiqc_report/filtered_multiqc_data/


        # No reason not to just run it on all
        echo "Purging paths in all filtered QC files..."
        find . -type f -exec bash \${WORKDIR}/clean-paths.sh '{}' ${params.baseDir} \\;
        cd \${WORKDIR}/${OUT_DIR}/


        echo "Re-zipping up filtered multiqc..."
        zip -r filtered_multiqc${params.assay_suffix}_report.zip filtered_multiqc_report/ && rm -rf filtered_multiqc_report/
        cd \${WORKDIR}

        echo "Purging paths from multiqc outputs completed successfully..."

        echo "Done! Paths purged successfully."
        """

}

/*
 * ========================================================================================
 * PROCESS: PACKAGE_PROCESSING_INFO
 * ========================================================================================
 *
 * SUMMARY:
 *   Purge file paths and zip processing info
 *
 * INPUTS:
 *   1. val: files_and_dirs
 *      Cardinality: one
 *      Description: Parameter value: files and directories to zip
 *
 * OUTPUTS:
 *   1. path: processing_info${params.assay_suffix}.zip (emit: zip)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/post_processing.config]
 *   Conda: envs/genelab.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process PACKAGE_PROCESSING_INFO {

    tag "Purging file paths and zipping processing info"

    input:
        val(files_and_dirs) 
    output:
        path("processing_info${params.assay_suffix}.zip"), emit: zip

    script:
        """
        cat `which clean-paths.sh` > clean-paths.sh
        chmod +x ./clean-paths.sh

        [ -d processing_info/ ] || mkdir processing_info/ && \\
        cp -r ${files_and_dirs.join(" ")} processing_info/

        echo "Purging file paths"
        find processing_info/ -type f -exec bash ./clean-paths.sh '{}' ${params.baseDir} \\;
        
        # Purge file paths and then zip
        zip -r processing_info${params.assay_suffix}.zip processing_info/
        """
} 


/*
 * ========================================================================================
 * PROCESS: GENERATE_README
 * ========================================================================================
 *
 * SUMMARY:
 *   Generate README for an OSD accession
 *
 * INPUTS:
 *   1. tuple: tuple val(name), val(email), val(output_prefix), val(OSD_accession),
                     val(protocol_id), val(FastQC_Outputs), val(Filtered_Sequence_Data),
                     val(Read_Based_Processing), val(Assembly_Based_Processing),
                     val(Assemblies), val(Genes), val(Annotations_And_Tax), val(Mapping), val(Combined_Output)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *
 *   2. path: processing_info
 *      Cardinality: one
 *      Description: Input file: processing info zip file
 *
 *   8. path: Bins
 *      Cardinality: one
 *      Description: Input file: Bins directory
 *
 *   9. path: MAGS
 *      Cardinality: one
 *      Description: Input file: MAGS directory
 *
 * OUTPUTS:
 *   1. path: README${params.assay_suffix}.txt (emit: readme)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/post_processing.config]
 *   Conda: envs/genelab.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process GENERATE_README {

    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating README for ${OSD_accession}"
    input:
        tuple val(name), val(email), val(output_prefix),
              val(OSD_accession), val(protocol_id),
              val(FastQC_Outputs), val(Filtered_Sequence_Data),
              val(Read_Based_Processing), val(Assembly_Based_Processing),
              val(Assemblies), val(Genes), val(Annotations_And_Tax),
              val(Mapping), val(Combined_Output)
        path(processing_info)
        path(Bins)
        path(MAGS)
    output:
        path("README${params.assay_suffix}.txt"), emit: readme

    script:
        """    
        GL-gen-processed-metagenomics-readme \\
             --output 'README${params.assay_suffix}.txt' \\
             --GLDS-ID '${OSD_accession}' \\
             --output-prefix '${output_prefix}' \\
             --name '${name}' \\
             --email '${email}' \\
             --protocol_ID '${protocol_id}' \\
             --assay_suffix '${params.assay_suffix}' \\
             --processing_zip_file '${processing_info}' \\
             --fastqc_dir '${FastQC_Outputs}' \\
             --filtered_reads_dir '${Filtered_Sequence_Data}' \\
             --read_based_dir '${Read_Based_Processing}' \\
             --assembly_based_dir '${Assembly_Based_Processing}' \\
             --assemblies_dir '${Assemblies}' \\
             --genes_dir   '${Genes}' \\
             --annotations_and_tax_dir '${Annotations_And_Tax}' \\
             --mapping_dir '${Mapping}' \\
             --bins_dir '${Bins}' \\
             --MAGs_dir '${MAGS}' \\
             --combined_output_dir  '${Combined_Output}' ${params.readme_extra}
        """

}

/*
 * ========================================================================================
 * PROCESS: VALIDATE_PROCESSING
 * ========================================================================================
 *
 * SUMMARY:
 *   Automated validation and verification
 *
 * INPUTS:
 *   1. tuple: tuple val(GLDS_accession), val(V_V_guidelines_link), val(output_prefix),
 *                   val(target_files), val(assay_suffix), val(log_dir_basename),
 *                   val(raw_suffix), val(raw_R1_suffix), val(raw_R2_suffix),
 *                   val(filtered_suffix), val(filtered_R1_suffix), val(filtered_R2_suffix)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *
 *   2. tuple: tuple path(Filtered_Sequence_Data),  path(Read_Based),
 *                   path(Assembly_Based), path(Assemblies), path(Mapping),
 *                   path(Genes), path(Annotation_And_Tax), path(Bins),
 *                   path(MAGS), path(Combined_Output), path(FastQC_Outputs)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *
 *   3. path: sample_ids_file
 *      Cardinality: one
 *      Description: Input file: one column sample ids file
 *
 *   4. path: README
 *      Cardinality: one
 *      Description: Input file: README file
 *
 *   5. path: processing_info
 *      Cardinality: one
 *      Description: Input file: processing info zip file
 *
 * OUTPUTS:
 *   1. path: ${GLDS_accession}_${output_prefix}metagenomics-validation.log (emit: log)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/post_processing.config]
 *   Conda: envs/genelab.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process VALIDATE_PROCESSING {

    tag "Running automated validation and verification...."

    input:
        // Labels
        tuple val(GLDS_accession), val(V_V_guidelines_link), val(output_prefix),
               val(target_files), val(assay_suffix), val(log_dir_basename),
               val(raw_suffix), val(raw_R1_suffix), val(raw_R2_suffix),
               val(filtered_suffix), val(filtered_R1_suffix), val(filtered_R2_suffix)
        // Directory paths
        tuple path(Filtered_Sequence_Data),  path(Read_Based),
              path(Assembly_Based), path(Assemblies), path(Mapping),
              path(Genes), path(Annotation_And_Tax), path(Bins), 
              path(MAGS), path(Combined_Output), path(FastQC_Outputs)
        // File paths
        path(sample_ids_file)
        path(README)
        path(processing_info) 

    output:
        path("${GLDS_accession}_${output_prefix}metagenomics-validation.log"), emit: log

    script:
        """
        GL-validate-processed-metagenomics-data \\
             --output '${GLDS_accession}_${output_prefix}metagenomics-validation.log' \\
             --GLDS-ID '${GLDS_accession}' \\
             --readme '${README}' \\
             --sample-IDs-file '${sample_ids_file}' \\
             --V_V_guidelines_link '${V_V_guidelines_link}' \\
             --processing_zip_file '${processing_info}' \\
             --output-prefix '${output_prefix}' \\
             --zip_targets '${target_files}' \\
             --assay_suffix '${assay_suffix}' \\
             --raw_suffix '${raw_suffix}' \\
             --raw_R1_suffix '${raw_R1_suffix}' \\
             --raw_R2_suffix '${raw_R2_suffix}' \\
             --filtered_suffix '${filtered_suffix}' \\
             --filtered_R1_suffix '${filtered_R1_suffix}' \\
             --filtered_R2_suffix '${filtered_R2_suffix}' \\
             --logs_dir_basename '${log_dir_basename}' \\
             --fastqc_dir ${FastQC_Outputs} \\
             --filtered_reads_dir ${Filtered_Sequence_Data} \\
             --read_based_dir ${Read_Based} \\
             --assembly_based_dir ${Assembly_Based} \\
             --assemblies_dir ${Assemblies} \\
             --genes_dir ${Genes} \\
             --annotations_and_tax_dir ${Annotation_And_Tax} \\
             --mapping_dir ${Mapping} \\
             --bins_dir ${Bins} \\
             --MAGs_dir ${MAGS} \\
             --combined_output_dir ${Combined_Output} ${params.validation_extra}
        """
}

/*
 * ========================================================================================
 * PROCESS: GENERATE_CURATION_TABLE
 * ========================================================================================
 *
 * SUMMARY:
 *   Generate a file association table for curation
 *
 * INPUTS:
 *   1. tuple: tuple val(GLDS_accession), val(output_prefix), val(assay_suffix),
 *                   val(raw_suffix), val(raw_R1_suffix), val(raw_R2_suffix),
 *                   val(filtered_suffix), val(filtered_R1_suffix), val(filtered_R2_suffix)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *
 *   2. tuple: tuple val(processing_zip_file), val(readme)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *
 *   3. tuple: tuple path(raw_reads_dir), path(filtered_reads_dir), path(read_based_dir),
 *                   path(assembly_based_dir), path(annotation_and_tax_dir), path(combined_output_dir)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *
 *   4. tuple: tuple path(Assemblies), path(Genes), path(Mapping), path(Bins), path(MAGS), path(FastQC_Outputs)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *
 *   5. path: input_table
 *      Cardinality: one
 *      Description: Input file: input_table
 *
 *   6. path: runsheet
 *      Cardinality: one
 *      Description: Input file: runsheet
 *
 * OUTPUTS:
 *   1. path: ${GLDS_accession}_${output_prefix}-associated-file-names.tsv (emit: curation_table)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/post_processing.config]
 *   Conda: envs/genelab.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process GENERATE_CURATION_TABLE {

    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating a file association table for curation..."

    input:
        // GeneLab accession and Suffixes
        tuple val(GLDS_accession), val(output_prefix), val(assay_suffix),
               val(raw_suffix), val(raw_R1_suffix), val(raw_R2_suffix),
               val(filtered_suffix), val(filtered_R1_suffix), val(filtered_R2_suffix)
        // File labels
        tuple val(processing_zip_file), val(readme)
        // Directory labels as paths - these paths are utilized as mere labels by the script
        tuple path(raw_reads_dir), path(filtered_reads_dir), path(read_based_dir),
              path(assembly_based_dir), path(annotation_and_tax_dir), path(combined_output_dir)
        // Directory paths
        tuple path(Assemblies), path(Genes), path(Mapping),
              path(Bins), path(MAGS), path(FastQC_Outputs) 
        path(input_table)
        path(runsheet)
        
    output:
        path("${GLDS_accession}_${output_prefix}-associated-file-names.tsv"), emit: curation_table

    script:
        def INPUT_TABLE = params.assay_table ? "--assay-table ${input_table}" : "--isa-zip  ${input_table}"
        """
        GL-gen-metagenomics-file-associations-table ${INPUT_TABLE} \\
                    --runsheet '${runsheet}' \\
                    --output '${GLDS_accession}_${output_prefix}-associated-file-names.tsv' \\
                    --GLDS-ID  '${GLDS_accession}' \\
                    --output-prefix '${output_prefix}' \\
                    --assay_suffix '${assay_suffix}' \\
                    --raw_suffix '${raw_suffix}' \\
                    --raw_R1_suffix '${raw_R1_suffix}' \\
                    --raw_R2_suffix '${raw_R2_suffix}' \\
                    --filtered_suffix '${filtered_suffix}' \\
                    --filtered_R1_suffix '${filtered_R1_suffix}' \\
                    --filtered_R2_suffix '${filtered_R2_suffix}' \\
                    --processing_zip_file '${processing_zip_file}' \\
                    --readme '${readme}' \\
                    --fastqc_dir '${FastQC_Outputs}' \\
                    --assemblies_dir '${Assemblies}' \\
                    --genes_dir '${Genes}' \\
                    --mapping_dir '${Mapping}' \\
                    --bins_dir '${Bins}' \\
                    --MAGs_dir '${MAGS}' \\
                    --raw_reads_dir '${raw_reads_dir}' \\
                    --filtered_reads_dir '${filtered_reads_dir}' \\
                    --read_based_dir  '${read_based_dir}' \\
                    --assembly_based_dir '${assembly_based_dir}' \\
                    --annotations_and_tax_dir '${annotation_and_tax_dir}' \\
                    --combined_output_dir '${combined_output_dir}' ${params.file_association_extra}
        """
}

/*
 * ========================================================================================
 * PROCESS: GENERATE_MD5SUMS
 * ========================================================================================
 *
 * SUMMARY:
 *   Generate md5sums for the files to be released on OSDR
 *
 * INPUTS:
 *   1. path: processing_info
 *      Cardinality: one
 *      Description: Input file: processing info zip file
 *
 *   2. path: README
 *      Cardinality: one
 *      Description: Input file: README files
 *
 *   3. val: dirs
 *      Cardinality: one
 *      Description: Parameter value: directory with files to generate md5sums for.
 *
 * OUTPUTS:
 *   1. path: processed_md5sum${params.assay_suffix}.tsv (emit: md5sum)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/post_processing.config]
 *   Conda: envs/genelab.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process GENERATE_MD5SUMS {
    
    tag "Generating md5sums for the files to be released on OSDR..."
 
    input:
        path(processing_info)
        path(README)
        val(dirs)

    output:
        path("processed_md5sum${params.assay_suffix}.tsv"), emit: md5sum
    script:
        """
        mkdir processing/ && \\
        cp -r ${dirs.join(" ")} ${processing_info} ${README} \\
              processing/

        # Generate md5sums
        find -L processing/ -type f -exec md5sum '{}' \\; |
        awk -v OFS='\\t' 'BEGIN{OFS="\\t"; printf "File Path\\tFile Name\\tmd5\\n"} \\
                {N=split(\$2,a,"/"); sub(/processing\\//, "", \$2); print \$2,a[N],\$1}' \\
                | grep -v "versions.txt" > processed_md5sum${params.assay_suffix}.tsv
        """
}

/*
 * ========================================================================================
 * PROCESS: GENERATE_PROTOCOL
 * ========================================================================================
 *
 * SUMMARY:
 *   Generate analysis protocol
 *
 * INPUTS:
 *   1. path: software_versions
 *      Cardinality: one
 *      Description: Input file: software versions file
 *
 *   2. val: protocol_id
 *      Cardinality: one
 *      Description: Parameter value: GeneLab pipeline id
 *
 * OUTPUTS:
 *   1. path: protocol.txt
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/post_processing.config]
 *   Conda: envs/genelab.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process GENERATE_PROTOCOL {

    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating your analysis protocol..."

    input:
        path(software_versions)
        val(protocol_id)
    output:
        path("protocol.txt")
    script:
        """
        generate_protocol.sh ${software_versions} ${protocol_id} > protocol.txt
        """
}
