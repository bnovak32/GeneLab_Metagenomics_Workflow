#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Terminal text color definitions
c_back_bright_red = "\u001b[41;1m";
c_bright_green    = "\u001b[32;1m";
c_blue            = "\033[0;34m";
c_reset           = "\033[0m";

params.help = false


/**************************************************
* HELP MENU  **************************************
**************************************************/
if(params.help){

  println()
  println("GeneLab Post Processing Pipeline: $workflow.manifest.version")
  println("USAGE:")
  println("Example: Submit and run jobs with slurm in singularity containers.")
  println("   > nextflow -C post_processing.config run post_processing.nf -resume -profile slurm,singularity")
  println()
  println("Required Parameters:")
  println("  --technology [STRING]  Sequencing technology. Options are 'illumina' or 'nanopore'. Default: null.")
  println("  --sample_type [STRING] Sample type. Options are 'standard' or 'low_biomass'. Default: null.")
  println("""-profile [STRING] Specifies the profile to be used to run the workflow. Options are [slurm, singularity, docker, and  conda].
	                    singularity, docker and conda will run the workflow locally using singularity, docker, and conda, respectively.
                      To combine profiles, separate two or more profiles with a comma. 
                      For example, to combine slurm and singularity profiles, pass 'slurm,singularity' as argument. """)	
  println("  --publishDir_mode [STRING]  Specifies how nextflow handles output file publishing. Options can be found here https://www.nextflow.io/docs/latest/process.html#publishdir Default: link.")
  println("  --glds_accession [STRING]  A Genelab GLDS accession number. Example GLDS-574. Default: empty string")
  println("  --osd_accession [STRING]  A Genelab OSD accession number. Example OSD-574. Default: empty string")
  println("  --name [STRING] The analyst's full name. E.g. 'FirstName A. LastName'.  Default: FirstName A. LastName")
  println("  --email [STRING] The analyst's email address. E.g. 'mail@nasa.gov'.  Default: mail@nasa.gov")
  println("  --assay_suffix [STRING]  Genelab's assay suffix. Default: _GLmetagenomics.")
  println("  --output_prefix [STRING] Unique name to tag onto output files. Default: empty string.")
  println("  --v_v_guidelines_link [URL] Genelab metagenomics data validation and verification guidelines link. Default: https://genelab-tools.arc.nasa.gov/confluence/pages/viewpage.action?pageId=8225175.")
  println("File Suffixes:")
  println("      --raw_suffix [STRING]  Raw reads suffix for datasets during processing. Default: _HRrm.fastq.gz.")  
  println("      --raw_R1_suffix [STRING]  Raw forward reads suffix for illumina datasets. Default: _R1_HRrm.fastq.gz.")
  println("      --raw_R2_suffix [STRING]  Raw reverse reads suffix for illumina datasets. Default: _R2_HRrm.fastq.gz.")
  println("      --filtered_suffix [STRING]  Filtered reads suffix for datasets during processing. Default: _filtered.fastq.gz.")
  println("      --filtered_R1_suffix [STRING]  Filtered forward reads suffix for illumina datasets. Default: _R1_filtered.fastq.gz.")
  println("      --filtered_R2_suffix [STRING]  Filtered reverse reads suffix for illumina datasets. Default: _R2_filtered.fastq.gz.")
  println("      --trimmed_suffix [STRING]  Trimmed reads suffix for nanopore datasets. Default: _trimmed.fastq.gz.")
  println("      --human_suffix [STRING]  Human removed reads suffix for datasets during processing. Default: _HRrm.fastq.gz.")
  println("      --decontam_suffix [STRING]  Decontaminated reads suffix for datasets during processing. Default: _decontam.fastq.gz.")
  println("      --decontam_R1_suffix [STRING]  Decontaminated forward reads suffix for illumina datasets. Default: _R1_decontam.fastq.gz.")
  println("      --decontam_R2_suffix [STRING]  Decontaminated reverse reads suffix for illumina datasets. Default: _R2_decontam.fastq.gz.")
  println("      --host_removed [BOOLEAN]  Indicates if host sequences were removed from the dataset. Default: false")
  println("      --host_suffix [STRING]  Host removed reads suffix for datasets during processing. Only used if --host-removed is true. Default: _HostRm.fastq.gz.")
  println("      --host_R1_suffix [STRING]  Host removed forward reads suffix for illumina datasets. Only used if --host-removed is true. Default: _R1_HostRm.fastq.gz.")
  println("      --host_R2_suffix [STRING]  Host removed reverse reads suffix for illumina datasets. Only used if --host-removed is true. Default: _R2_HostRm.fastq.gz.")
  println("Files:")
  println("    --human_removed_summary  [PATH]  Path to the human reads removed summary file. Default: empty string")
  println("    --run_command  [PATH] File containing the nextflow run command used in processing. Default: ./processing_scripts/command.txt")
  println("    --processing_commands  [PATH] File containing all the process names and scripts used during processing. Default: ./processing_scripts/nextflow_processing_info_GLAmpliseq.txt")
  println("    --assay_table  [PATH] GLDS assay table generated after running the processing pipeline with accession number as input.")
  println("                   Example, ../Genelab/a_OSD-574_metagenomic-sequencing_whole-genome-shotgun-sequencing_illumina.txt. Default: empty string")
  println("    --isa_zip  [PATH] Genelab ISA zip files containing an assay atable for the OSD accession. This is only required if --files.assay_table is not set.")
  println("                   Example, ../Genelab/OSD-574_metadata_OSD-574-ISA.zip. Default: empty string")
  println("    --runsheet  [PATH] Input csv file used to run the processing pipeline with the first column containing sample names. This is the value set to the paremater --input_file when run the processing pipeline with a csv file as input otherwise it is the GLfile.csv in the GeneLab directory if --GLDS_accession was used as input. Example '../GeneLab/GLfile.csv'.  Default: null")
  println("    --software_versions  [PATH] A file generated after running the processing pipeline listing the software versions used. Default: ../Metadata/software_versions.txt")
  println()
  println("Directories:")
  println("    --output_dir  [PATH] Specifies the directory where outputs of this post-processing workflow will be published. Default: ../Post_Processing/")
  println("    --processing_dir  [PATH] Specifies the directory where outputs of the main processing workflow were published. Default: ../../processing/")
  println()
  println("Optional arguments:")  
  println("    --help [BOOLEAN] Print this help message and exit")
  println("    --debug [BOOLEAN] Set to true if you'd like to see the values of your set parameters printed to the terminal. Default: false.")
  println()
  println("Paths to existing conda environments to use otherwise a new one will be created using the yaml file in envs/.")
  println("    --conda_genelab [PATH] Path to a conda environment containing genelab-utils. Default: null.")
  exit 0
}


/************************************************
*********** Show pipeline parameters ************
*************************************************/
if(params.debug){

log.info """${c_blue}
         GeneLab Post Processing Pipeline: $workflow.manifest.version
         
         You have set the following parameters:
         Technology: ${params.technology}
         Sample Type: ${params.sample_type}
         Profile: ${workflow.profile} 
         Analyst's Name: ${params.name}
         Analyst's Email: ${params.email}
         GLDS Accession: ${params.glds_accession}
         OSD Accession: ${params.osd_accession}
         Assay Suffix: ${params.assay_suffix}
         Output Prefix: ${params.output_prefix}
         V & V Link: ${params.v_v_guidelines_link}
         Human Genome Reference: ${params.genome}
         Nextflow directory publishing mode: ${params.publishDir_mode}
        

         Suffixes:
         Raw Suffix: ${params.raw_suffix}
         Raw R1 suffix: ${params.raw_R1_suffix}
         Raw R2 suffix: ${params.raw_R2_suffix}          
         Filtered Suffix: ${params.filtered_suffix}
         Filtered R1 suffix: ${params.filtered_R1_suffix}
         Filtered R2 suffix: ${params.filtered_R2_suffix}
         Trimmed Suffix: ${params.trimmed_suffix}
         Human Removed Suffix: ${params.human_suffix}
         Decontaminated Suffix: ${params.decontam_suffix}
         Decontaminated R1 Suffix: ${params.decontam_R1_suffix} 
         Decontaminated R2 Suffix: ${params.decontam_R2_suffix}
         Host Removed Suffix: ${params.host_suffix}
         Host Removed R1 Suffix: ${params.host_R1_suffix}
         Host Removed R2 Suffix: ${params.host_R2_suffix} 


         Files:
         Human Removed Summary: ${params.human_summary}
         Nextflow Command: ${params.run_command}
         Processing Commands: ${params.processing_commands}
         Assay Table: ${params.assay_table}
         ISA Zip: ${params.isa_zip}
         Input Runsheet: ${params.runsheet}
         Software Versions: ${params.software_versions}
        
         Boolean Flags:
         Host Removed: ${params.host_removed}
         Single Ended: ${params.single_end}
         Concated Reads: ${params.concated}

         Directories:
         Pipeline Outputs: ${params.output_dir}
         Processing : ${params.processing_dir}
         Jinja Templates: ${params.templates}
         """
}


include { CLEAN_MULTIQC_PATHS as CLEAN_RAW_PATHS; 
          CLEAN_MULTIQC_PATHS as CLEAN_FILTERED_PATHS;
          CLEAN_MULTIQC_PATHS as CLEAN_TRIMMED_PATHS;
          CLEAN_MULTIQC_PATHS as CLEAN_HUMAN_PATHS;
          CLEAN_MULTIQC_PATHS as CLEAN_DECONTAM_PATHS;
          CLEAN_MULTIQC_PATHS as CLEAN_HOST_PATHS} from './modules/genelab.nf'

include { PACKAGE_PROCESSING_INFO;  GENERATE_READ_STATS; VALIDATE_PROCESSING;
           GENERATE_CURATION_TABLE; GENERATE_README; 
           GENERATE_MD5SUMS; GENERATE_PROTOCOL} from './modules/genelab.nf'

workflow {

        // Make sure accession numbers are set
        if(!params.glds_accession || !params.osd_accession){
           error("""${c_back_bright_red}ACCESSION ERROR!. 
                    Please supply both --glds_accession and --osd_accession.
                    They can be any string you choose but must be set.
                 ${c_reset}""")
        }

        // Make sure technology is set
        if(!params.technology ){
           error("""${c_back_bright_red}PARAMETER ERROR!. 
                    Please supply the --technology.
                    One of illumina or nanopore must be set.
                 ${c_reset}""")
        }

        // Make sure  sample type is set
        if(!params.sample_type){
           error("""${c_back_bright_red}PARAMETER ERROR!. 
                    Please supply the --sample_type.
                    One of standard or low_biomass must be set.
                 ${c_reset}""")
        }

       // ---------------------- Input channels -------------------------------- //
       // Input Value channels
       meta_ch   =  channel.of([name: params.name, email: params.email, output_prefix: params.output_prefix, genome: params.genome,
                                protocol_id: params.protocol_id, technology: params.technology, sample_type: params.sample_type,
                                osd_accession: params.osd_accession, glds_accession: params.glds_accession, concated: params.concated,
                                v_v_guidelines_link: params.v_v_guidelines_link, assay_suffix: params.assay_suffix,  
                                raw_suffix: params.raw_suffix, raw_R1_suffix: params.raw_R1_suffix, raw_R2_suffix: params.raw_R2_suffix,
                                filtered_suffix: params.filtered_suffix, filtered_R1_suffix: params.filtered_R1_suffix,
                                filtered_R2_suffix: params.filtered_R2_suffix,
                                trimmed_suffix: params.trimmed_suffix, human_suffix: params.human_suffix,
                                decontam_suffix: params.decontam_suffix, decontam_R1_suffix: params.decontam_R1_suffix,
                                decontam_R2_suffix: params.decontam_R2_suffix,
                                host_removed: params.host_removed, host_suffix: params.host_suffix, host_R1_suffix: params.host_R1_suffix,
                                host_R2_suffix: params.host_R2_suffix])

       // Input files
       software_versions   =  channel.fromPath(params.software_versions, checkIfExists: true)

       // Processing directory containing files to be validated and packaged for OSDR release
       processing_dir = channel.fromPath(params.processing_dir, type: 'dir', checkIfExists: true)
       // Jinja templates for protocol generation
       templates = channel.fromPath(params.templates, type: 'dir', checkIfExists: true)

       // If the assay table is provided use it as the input table otherwise use the isa_zip
       assay_table_ch = channel.fromPath(params.assay_table ?  params.assay_table : params.isa_zip,
                                          checkIfExists: true)

      // Runsheet used to execute the processing workflow
      runsheet_ch = channel.fromPath(params.runsheet, checkIfExists: true)

      // Human removed reads summary generated after running the human reads removal workflow
      human_summary_ch = params.human_summary ? channel.fromPath(params.human_summary, checkIfExists: true) : file('hsempty.txt')

      // Files to be packaged in processing_info.zip
      files_ch = channel.of(params.run_command, params.processing_commands,
                            params.software_versions, params.runsheet)
                                      .collect()
                                      .map{ run_command, processing_commands, software_versions, runsheet ->
                                            tuple( file(run_command, checkIfExists: true),
                                                   file(processing_commands, checkIfExists: true),
                                                   file(software_versions, checkIfExists: true),
                                                   file(runsheet, checkIfExists: true)
                                                 ) }

        // ---------------------- Post-processing begins ---------------------------------//
        PACKAGE_PROCESSING_INFO(files_ch)
       
        // Clean paths in mutiqc reports  
        raw_multiqc      =  channel.fromPath(params.raw_zip,  checkIfExists: true)
        filtered_multiqc =  channel.fromPath(params.filtered_zip,  checkIfExists: true)
        trimmed_multiqc  =  params.trimmed_zip ? channel.fromPath(params.trimmed_zip,  checkIfExists: true) : channel.empty()
        human_multiqc    =  params.human_zip ? channel.fromPath(params.human_zip,  checkIfExists: true) : channel.empty()
        decontam_multiqc =  params.decontam_zip ? channel.fromPath(params.decontam_zip,  checkIfExists: true) : channel.empty()
        host_multiqc     =  params.host_zip ? channel.fromPath(params.host_zip,  checkIfExists: true) : channel.empty()

        CLEAN_RAW_PATHS(raw_multiqc)
        CLEAN_FILTERED_PATHS(filtered_multiqc)
        if (params.technology == "nanopore") {
          CLEAN_TRIMMED_PATHS(trimmed_multiqc)
          CLEAN_HUMAN_PATHS(human_multiqc)
        }

        if (params.sample_type == "low_biomass"){ CLEAN_DECONTAM_PATHS(decontam_multiqc) } 
        
        if (params.host_removed) { CLEAN_HOST_PATHS(host_multiqc) }

         
        // Automatic verification and validation
        VALIDATE_PROCESSING(meta_ch, processing_dir, runsheet_ch, PACKAGE_PROCESSING_INFO.out.zip) 

        GENERATE_READ_STATS(raw_multiqc, human_summary_ch, filtered_multiqc,
                            trimmed_multiqc.ifEmpty(file('trim_empty.txt')), 
                            human_multiqc.ifEmpty(file('human_empty.txt')), 
                            decontam_multiqc.ifEmpty(file('decontam_empty.txt')), 
                            host_multiqc.ifEmpty(file('host_empty.txt')))
        
        /*
        // Generate curation file association table
        GENERATE_CURATION_TABLE(meta_ch, processing_dir, assay_table_ch, runsheet_ch,
                                human_summary_ch, GENERATE_READ_STATS.out.stats, 
                                VALIDATE_PROCESSING.out.json)
        */

        // Generate README file
        GENERATE_README(meta_ch, PACKAGE_PROCESSING_INFO.out.zip, 
                      runsheet_ch, VALIDATE_PROCESSING.out.json)
        
        // Generate md5sums
        GENERATE_MD5SUMS(processing_dir, PACKAGE_PROCESSING_INFO.out.zip,
                            GENERATE_README.out.readme)
        
        // Write methods
        GENERATE_PROTOCOL(meta_ch, software_versions, templates)
}


workflow.onComplete {

    println("${c_bright_green}Pipeline completed at: $workflow.complete")
    println("""Execution status: ${ workflow.success ? 'OK' : "${c_back_bright_red}failed" }""")
    log.info ( workflow.success ? "\nDone! Workflow completed without any error\n" : "Oops .. something went wrong${c_reset}" )

    if ( workflow.success ) {

    println("Post-processing Outputs: ${params.output_dir} ${c_reset}")
    println()

    }
}
