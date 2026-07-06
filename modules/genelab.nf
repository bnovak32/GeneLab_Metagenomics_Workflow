#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * ========================================================================================
 * PROCESS: PACKAGE_PROCESSING_INFO
 * ========================================================================================
 *
 * SUMMARY:
 *   Purge file paths and zip processing info
 *
 * INPUTS:
 *   1. val: files
 *      Cardinality: one
 *      Description: Parameter value: files to zip
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
        val(files) 
    output:
        path("processing_info${params.assay_suffix}.zip"), emit: zip

    script:
        """
        cat `which clean-paths.sh` > clean-paths.sh
        chmod +x ./clean-paths.sh

        [ -d processing_info/ ] || mkdir processing_info/ && \\
        cp -r ${files.join(" ")} processing_info/

        echo "Purging file paths"
        find processing_info/ -type f -exec bash ./clean-paths.sh '{}' ${params.baseDir} \\;
        
        # Purge file paths and then zip
        zip -r processing_info${params.assay_suffix}.zip processing_info/
        """
} 


/*
 * ========================================================================================
 * PROCESS: CLEAN_MULTIQC_PATHS
 * ========================================================================================
 *
 * SUMMARY:
 *   Purge genelab paths from MultiQC zip file
 *
 * INPUTS:
 *   1. path:  zip_file
 *      Cardinality: one
 *      Description: Input file: MultiQC zip file
 *
 * OUTPUTS:
 *   1. path: "*.zip", emit: clean_zip
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

process CLEAN_MULTIQC_PATHS {
    tag "Purging genelab paths from MultiQC zip file..."
    input:
        path(zip_file) // e.g. raw_multiqc_data.zip
    output:
        path("*.zip"), emit: clean_zip
    script:
        """ 
        dirname=`basename -s .zip ${zip_file}`
        unzip ${zip_file} && rm ${zip_file}
        clean_multiqc_paths.py ${dirname} .
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
 *   1. val: meta
 *      Cardinality: one
 *      Description: Metadata object containing multiple channel elements
 *
 *
 *   2. path: runsheet
 *      Cardinality: one
 *      Description: Input file: run sheet csv file with sample names in the first column.
 *
 *
 *   3. path: processing_info
 *      Cardinality: one
 *      Description: Input file: processing info zip file
 *
 * OUTPUTS:
 *   1. path: ${meta.GLDS_accession}_${meta.output_prefix}metagenomics-validation.log (emit: log)
 *   2. path: ${meta.GLDS_accession}_${meta.output_prefix}metagenomics-validation.manifest.json (emit: json)
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
        // Labeling and suffixes
        val(meta)
        // File paths
        path(runsheet)
        path(processing_info) 

    output:
        path("${meta.GLDS_accession}_${meta.output_prefix}metagenomics-validation.log"), emit: log
        path("${meta.GLDS_accession}_${meta.output_prefix}metagenomics-validation.manifest.json"), emit: json

    script:
      def single_end_flag = params.single_end ? "--single-ended" : ""
      def host_removed_flag = ""
      def decontam_flag = ""

      if (params.host_removed) {


        if ( (params.technology == "illumina"  && params.single_end) || (params.technology == "nanopore") ) {
          // illumina single-end or nanopore
         host_removed_flag = " --host-removed  --host-suffix ${meta.host_suffix}"

        }else if (params.technology == "illumina") {
         // illumina paired-end
         host_removed_flag = " --host-removed  --host-suffix ${meta.host_suffix} --host_R1-suffix  ${meta.host_R1_suffix} --host_R2-suffix  ${meta.host_R2_suffix}"

        } 


      }
       

      if (params.sample_type == "low_biomass") {

        if ( (params.technology == "illumina"  && params.single_end) ||  params.technology == "nanopore" ) {
            // illumina single-end or nanopore
            decontam_flag = "--decontam-suffix ${meta.decontam_suffix}"

        }else if (params.technology == "illumina") {

            // illumina paired-end
            decontam_flag = "--decontam-suffix ${meta.decontam_suffix} --decontam_R1-suffix ${meta.decontam_R1_suffix} --decontam_R2-suffix ${meta.decontam_R2_suffix}"

       }

       }

        """
        if [ ${meta.technology} == "illumina" ]; then

        # Illumina
        GL-validate-processed-metagenomics-data \\
             --technology '${meta.technology}' \\
             --sample-type '${meta.sample_type}' \\  
             --output '${meta.glds_accession}_${meta.output_prefix}metagenomics-validation.log' \\
             --manifest '${meta.glds_accession}_${meta.output_prefix}metagenomics-validation.manifest.json' \\
             --glds-id '${meta.glds_accession}' \\
             --runsheet '${runsheet}' \\
             --V-V-guidelines_link '${meta.v_v_guidelines_link}' \\
             --processing-zip-file '${processing_info}' \\
             --output-prefix '${meta.output_prefix}' \\
             --assay-suffix '${meta.assay_suffix}' \\
             --raw-suffix '${meta.raw_suffix}' \\
             --raw-R1-suffix '${meta.raw_R1_suffix}' \\
             --raw-R2-suffix '${meta.raw_R2_suffix}' \\
             --filtered-suffix '${meta.filtered_suffix}' \\
             --filtered-R1-suffix '${meta.filtered_R1_suffix}' \\
             --filtered-R2-suffix '${meta.filtered_R2_suffix}'  \\
             ${decontam_flag} ${host_removed_flag} ${single_end_flag}

        else

        # Nanopore
        GL-validate-processed-metagenomics-data \\
             --technology '${meta.technology}' \\
             --sample-type '${meta.sample_type}' \\  
             --output '${meta.glds_accession}_${meta.output_prefix}metagenomics-validation.log' \\
             --manifest '${meta.glds_accession}_${meta.output_prefix}metagenomics-validation.manifest.json"' \\
             --glds-id '${meta.glds_accession}' \\
             --runsheet '${runsheet}' \\
             --V-V-guidelines_link '${meta.V_V_guidelines_link}' \\
             --processing-zip-file '${processing_info}' \\
             --output-prefix '${meta.output_prefix}' \\
             --assay-suffix '${meta.assay_suffix}' \\
             --raw-suffix '${meta.raw_suffix}' \\
             --filtered-suffix '${meta.filtered_suffix}' \\
             --trimmed-suffix '${meta.trimmed_suffix}' \\
             --human-suffix '${meta.human_suffix}' \\
             --host-suffix '${meta.host_suffix}' \\
             ${decontam_flag} ${host_removed_flag}

        fi
        """
}


/*
 * ========================================================================================
 * PROCESS: GENERATE_READ_STATS
 * ========================================================================================
 *
 * SUMMARY:
 *   Generate read statistics for the processed data
 *
 * INPUTS:
 *   1. path: raw_zip 
 *      Cardinality: one
 *      Description: Input file: raw multiqc data zip file.
 *
 *   2. path: human_summary
 *      Cardinality: one
 *      Description: Input file: human removed summary file with read counts for raw and human removed data.
 *
 *   3. path: filtered_zip
 *      Cardinality: one
 *      Description: Input file: filtered multiqc data zip file.
 *
 *   4. path: trimmed_zip
 *      Cardinality: one
 *      Description: Input file: trimmed multiqc data zip file.
 *
 *   5. path: human_zip
 *      Cardinality: one
 *      Description: Input file: human removed multiqc data zip file.
 *
 *   6. path: decontam_zip
 *      Cardinality: one
 *      Description: Input file: decontam multiqc data zip file.
 *
 *   7. path: host_zip
 *      Cardinality: one
 *      Description: Input file: host removed multiqc data zip file.
 *
 * OUTPUTS:
 *   1. path: reads_statistics.tsv (emit: stats)
 *
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


process GENERATE_READ_STATS {
    tag "Generating read statistics..."
    input:
        path(raw_zip) // raw_multiqc_data.zip
        path(human_summary) // human_removed_summary.tsv
        path(filtered_zip) // filtered_multiqc_data.zip
        path(trimmed_zip) //  trimmed_multiqc_data.zip
        path(human_zip) // HRrm_multiqc_data.zip
        path(decontam_zip) // decontam_multiqc_data.zip
        path(host_zip)  // HostRm_multiqc_data.zip
    output:
        path("reads_statistics.tsv"), emit: stats
    script:
        
        blank_flag = params.decontam_zip ? "--blank-removed ${decontam_zip}" : ""
        host_flag  = params.host_removed ? "--host-removed ${host_zip}" : ""

        if (params.technology == "nanopore") {
            // Nanopore
            raw_flag = "--raw ${raw_zip} --trimmed ${trimmed_zip} --human-removed ${human_zip} "
        } else if (params.human_summary ){
            // Illumina with human summary file provided - it contains the read counts for raw data and human removed
            //  data so can be used to generate statistics without needing the full multiqc report for the raw data
            raw_flag = "--human-removed-summary ${human_summary}"
        } else {
            // Illumina with no human summary file provided - use the raw multiqc report to get the read counts for
            //  raw data and human removed data and generate statistics that way
            raw_flag = "--raw ${raw_zip} --human-removed ${human_zip}"
        }

        """
        generate_reads_statistics.py --output reads_statistics.tsv \\
                                     --sample-type ${params.sample_type} \\
                                     --technology ${params.technology} \\
                                     --filtered ${filtered_zip} \\
                                     ${raw_flag} ${blank_flag} ${host_flag}
                                     
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
 *   1. Map:  meta
 *      Cardinality: one
 *      Description: Map input combining multiple channel elements
 *
 *   2. path: processing_dir
 *      Cardinality: one
 *      Description: Input directory: processing directory with files to generate md5sums for
 *
 *   3. path: input_table
 *      Cardinality: one
 *      Description: Input file: input_table
 *
 *   4. path: runsheet
 *      Cardinality: one
 *      Description: Input file: processing runsheet
 *
 *   5. path: read_statistics
 *      Cardinality: one
 *      Description: Input file: read_statistics
 *
 *   6. path: validation_manifest
 *      Cardinality: one
 *      Description: Input file: validation manifest JSON file   
 * 
 *
 * OUTPUTS:
 *   1. path: ${meta.glds_accession}_${meta.output_prefix}-associated-file-names.tsv (emit: curation_table)
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
        val(meta)
        path(processing_dir)
        path(input_table)
        path(runsheet)
        path(human_summary) // human_removed_summary.tsv
        path(read_statistics)
        path(validation_manifest)
        
    output:
        path("${meta.glds_accession}_${meta.output_prefix}-associated-file-names.tsv"), emit: curation_table

    script:
        def INPUT_TABLE = params.assay_table ? "--assay-table ${input_table}" : "--isa-zip  ${input_table}"
        def hum_summary = params.human_summary ? "--hrrm_stats ${human_summary}" : ""
        """
        update_assay_table.py ${INPUT_TABLE} \\
                    --technology ${params.technology} \\
                    --sample-type ${params.sample_type} \\
                    --runsheet '${runsheet}' \\
                    --output '${meta.glds_accession}_${meta.output_prefix}-associated-file-names.tsv' \\ # it doesn't have the output agument so I wonder what it will be
                    --processed_file_path ${processing_dir} \\
                    --glds_accession  '${meta.glds_accession}' \\
                    --output-prefix '${meta.output_prefix}' \\
                    --assay_suffix '${meta.assay_suffix}' \\
                    --read_stats_file '${read_statistics}' \\
                    --validation_output_json '${meta.validation_manifest}' ${hum_summary}
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
 *   1. Map:  meta
 *      Cardinality: one
 *      Description: Map input combining multiple channel elements
 *
 *   2. path: processing_info
 *      Cardinality: one
 *      Description: Input file: processing info zip file
 *
 *   3. path: runsheet
 *      Cardinality: one
 *      Description: Input file: processing runsheet
 *
 *   4. path: validation_manifest
 *      Cardinality: one
 *      Description: Input file: validation manifest JSON file
 *
 *
 * OUTPUTS:
 *   1. path: ${meta.output_prefix}README${meta.assay_suffix}.txt (emit: readme)
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
    tag "Generating README for ${meta.osd_accession}"
    input:
        val(meta)
        path(processing_info)
        path(runsheet)
        path(validation_manifest)
    output:
        path("${meta.output_prefix}README${meta.assay_suffix}.txt"), emit: readme

    script:
        def runsheet_name = runsheet.getName()
        def host_removed_flag = params.host_removed ? "--host-removed" : ""  
        """    
        GL-gen-processed-metagenomics-readme \\
             --output '${meta.output_prefix}README${meta.assay_suffix}.txt' \\
             --glds-id '${meta.osd_accession}' \\
             --name '${meta.name}' \\
             --email '${meta.email}' \\
             --protocol-id '${meta.protocol_id}' \\
             --assay_suffix '${meta.assay_suffix}' \\
             --runsheet '${runsheet_name}' \\
             --validation_output_json '${meta.validation_manifest}' ${host_removed_flag}
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
 *   1. path: processing_dir
 *      Cardinality: one
 *      Description: Input directory: processing directory with files to generate md5sums for
 *
 *   2. path: processing_info
 *      Cardinality: one
 *      Description: Input file: processing info zip file
 *
 *   3. path: README
 *      Cardinality: one
 *      Description: Input file: README files
 *
 *
 * OUTPUTS:
 *   1. path: ${params.output_prefix}processed_md5sum${params.assay_suffix}.tsv (emit: md5sum)
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
        path(processing_dir)
        path(processing_info)
        path(README)

    output:
        path("${params.output_prefix}processed_md5sum${params.assay_suffix}.tsv"), emit: md5sum
    script:
        """
        # Generate md5sums
        generate_md5sums.py --outdir ${processing_dir} \\
                            --output_prefix '${params.output_prefix}' \\
                            --assay-suffix '${params.assay_suffix}'  
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
 *   1. val: meta
 *      Cardinality: one
 *      Description: Map input combining multiple channel elements [protocol_id, sample_type, technology]
 *   2. path: software_versions
 *      Cardinality: one
 *      Description: Input file: software versions file
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
        val(meta)
        path(software_versions)
    output:
        path("protocol.txt")
    script:
        """
        generate_protocol.py --versions-file ${software_versions} \\
                             --protocol-id ${meta.protocol_id} \\
                             --sample-type ${meta.sample_type} \\
                             --technology ${meta.technology}  > protocol.txt
        """
}


