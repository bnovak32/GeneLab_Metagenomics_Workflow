#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Terminal text color defintions
c_back_bright_red = "\u001b[41;1m";
c_bright_green    = "\u001b[32;1m";
c_blue            = "\033[0;34m";
c_reset           = "\033[0m";

params.help = false

/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println()
  println("Nextflow Metagenomics Consensus Workflow: $workflow.manifest.version")
  println("USAGE:")
  println("Example 1: Illumina Standard Sample.")
  println("   > nextflow run main.nf -resume -profile slurm,singularity --technology illumina --sample_type standard --input_file PE_samples.csv  --errorStrategy 'ignore'")
  println()
  println("Example 2: Nanopore Standard Sample. Start with pod5 files as input, and run the pipeline on slurm in conda environments.")
  println("   > nextflow run main.nf -resume -profile slurm,mamba --technology nanopore --sample_type standard --input_file input_dir_barcodes.csv --input_type 'directory' --kit_name 'SQK-RPB114-24' --input_dir /path/to/pod5/directory/ --errorStrategy 'ignore' ")
  println()
  println("Example 3: Illumina Low-Biomass Sample.")
  println("   > nextflow run main.nf -resume -profile slurm,singularity --technology illumina --sample_type low_biomass --input_file PE_samples.csv  --errorStrategy 'ignore' ")
  println()
  println("Example 4: : Nanopore Low-Biomass Sample. Start with multiple FASTQ files per sample as input, and run the pipeline on slurm in singularity environments.")
  println("   > nextflow run main.nf -resume -profile slurm,singularity --technology nanopore --sample_type low_biomass --input_file multiple.csv --input_type 'multiple' --errorStrategy 'ignore' ")
  println()
  println("Example 5: Run jobs locally in conda environments, supply a GLDS accession, and specify the path to an existing conda environment.")
  println("   > nextflow run main.nf -resume -profile mamba --technology illumina --sample_type standard --accession OSD-574 --conda_megahit <path/to/existing/conda/environment>")
  println()
  println("Required arguments:")
  println("""-profile [STRING] Specifies the profile to be used to run the workflow. Options are [slurm, singularity, docker, mamba, and  conda].
	                    singularity, docker and conda will run the pipeline locally using singularity, docker, and conda, respectively.
                      To combine profiles, separate two or more profiles with comma. For example, to combine slurm and singularity profiles, pass 'slurm,singularity' as argument. """)			 
  println("--input_file  [PATH] A 3-column (single-end) or 4-column (paired-end) csv input file (sample_id, forward, [reverse,] paired). Required only if a GLDS accession is not provided. Default : null")
  println("   Please see the files: SE_file.csv and PE_file.csv for single-end and paired-end examples, respectively.")
  println("   The sample_id column should contain unique sample ids.")
  println("   The forward and reverse columns should contain the absolute or relative path to the sample's forward and reverse reads.")
  println("   The paired column should be true for paired-end or anything else for single-end reads.")
  println()
  println("Optional arguments:")  
  println("  --help  Print this help message and exit")
  println("  --workflow [STRING] Specifies that workflow to be run. Options are one of [read-based, assembly-based, both]. Default: both.")
  println("  --publishDir_mode [STRING]  Specifies how nextflow handles output file publishing. Options can be found here https://www.nextflow.io/docs/latest/process.html#publishdir Default: link.")
  println("  --errorStrategy [STRING] Specifies how nextflow handles errors. Options can be found here https://www.nextflow.io/docs/latest/process.html#errorstrategy. Default: ignore")
  println("  --multiqc_config [PATH] Path to a custom multiqc config file. Default: config/multiqc.config.")
  println("  --use_gtdbtk_scratch_location [BOOLEAN] Should a scratch location be used to store GTDBTK temp files? true or false.")
  println("    Scratch directory for gtdb-tk, if wanting to use disk space instead of RAM, can be memory intensive;")
  println("    see https://ecogenomics.github.io/GTDBTk/faq.html#gtdb-tk-reaches-the-memory-limit-pplacer-crashes")
  println("    leave empty if wanting to use memory, the default, put in quotes the path to a directory that")
  println("    already exists if wanting to use disk space. Default: false.")
  println()
  println("MAG parameters: MAG filtering cutoffs based on checkm quality assessments (in percent); see https://github.com/Ecogenomics/CheckM/wiki/Reported-Statistics.")
  println("	 --min_est_comp [INT] Minimum estimated completion. Default: 90.") 
  println("	 --max_est_redund [INT] Minimum estimated redundancy. Default: 10.") 
  println("	 --max_est_strain_het [INT] Minimum estimated strain heterogeneity. Default: 50.")
  println("	 --reduced_tree [STRING] reduced_tree option for checkm, limits the RAM usage to 16GB; https://github.com/Ecogenomics/CheckM/wiki/Genome-Quality-Commands#tree.")
  println("    'True' for yes, anything else will be considered 'False' and the default full tree will be used. Default: 'True'. ")
  println("	 --max_mem [INT] Maximum memory allowed, passed to megahit assembler. Can be set either by proportion of available on system, e.g. 0.5")
  println("    or by absolute value in bytes, e.g. 100e9 would be 100 GB. Default: 100e9.")
  println()
  println("	 --pileup_mem [STRING] pileup.sh paramater for calculating contig coverage and depth. Memory used by bbmap's pileup.sh (within the GET_COV_AND_DET process). ")
  println("	   passed as the -Xmx parameter, 20g means 20 gigs of RAM, 20m means 20 megabytes.")
  println("	   5g should be sufficient for most assemblies, but if that rule is failing, this may need to be increased.Default: '5g' .")
  println("	 --block_size [int] Block size variable for CAT/diamond, lower value means less RAM usage; see https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options. Default: 4.")
  println()
  println("Output directories:")
  println("      --merged_dir [PATH] Specifies where the raw/merged reads will be published. Default: ../Merged_Sequence_Data/.")
  println("      --filtered_reads_dir [PATH] Specifies where filtered reads will be published.  Default: ../Filtered_Sequence_Data/.")
  println("      --trimmed_dir [PATH] Specifies where trimmed reads will be published. Default: ../Trimmed_Sequence_Data/.")
  println("      --human_removed_dir [PATH] Specifies where human-removed reads will be published. Default: ../HR-removed_Sequence_Data/.")
  println("      --decontaminated_dir [PATH] Specifies where decontaminated reads will be published.  Default: ../Decontaminated_Sequence_Data/.")
  println("      --host_removed_dir [PATH] Specifies where host-removed reads will be published.  Default: ../HostRM-removed_Sequence_Data/.")
  println("      --genelab_dir [PATH] Specifies where Genelab outputs will be published.  Default: ../GeneLab/.")
  println("      --logs_dir [PATH] Specifies where tool log outputs will be published.  Default: ../Logs/.")
  println("      --metadata_dir [PATH] Specifies where metadata outputs (e.g software versions) will be published.  Default: ../Metadata/.")
  println("      --read_based_dir [PATH] Read-based analysis outputs directory.  Default: ../Read-based_Processing/.")
  println("      --assembly_based_dir [PATH] Specifies where the results of assembly-based analysis will be published. Default: ../Assembly-based_Processing/.")
  println()
  println("Genelab specific arguements:")
  println("      --accession [STRING]  A Genelab accession number if the --input_file parameter is not set. If this parameter is set, it will ignore the --input_file parameter. Default: null.")
  println("      --RawFilePattern [STRING]  If we do not want to download all files (which we often won't), we can specify a pattern here to subset the total files.")
  println("                                 For example, if we know we want to download just the fastq.gz files, we can say 'fastq.gz'. We can also provide multiple patterns")
  println("                                 as a comma-separated list. For example, If we want to download the fastq.gz files that also have 'NxtaFlex', 'metagenomics', and 'raw' in") 
  println("                                 their filenames, we can provide '-p fastq.gz,NxtaFlex,metagenomics,raw'. Default: null.")
  println("      --assay_suffix [STRING]  Genelab's assay suffix. Default: _GLmetagenomics.")
  println("      --additional_filename_prefix [STRING] additional prefix to add to output files that describe more than one sample (to make them unique compared to other datasets).")
  println("      include separator at end if adding one, e.g. Swift1S_ if wanted. Default: empty string .")
  println()
  println("Paths to existing databases and database links.")
  println("        --DB_ROOT [PATH]   FULL PATH to root directory where the databases will be downloaded if they don't exist.") 
  println("                  Relative paths such as '~/' and '../' will fail, please don't use them. Default: ../Reference_DBs/ ")
  println("CAT database directory strings:")
  println("    The strings below will be added to the end of the --database.cat_db path arguement provided below.")
  println("         --cat_taxonomy_dir [PATH] CAT taxonomy database directory. Default: 2021-01-07_taxonomy/.")
  println("         --cat_db_sub_dir [PATH] CAT database sub directory. Default: 2021-01-07_CAT_database/.")
  println("         --CAT_DB_LINK [URL] CAT database online download link. Default: https://tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20210107.tar.gz.")
  println("CAT database ")
  println("         --cat_db [PATH] Path to CAT database. Example, /path/to/Reference_DBs/CAT_prepare_20210107/. Default: null.")
  println("Humann database:")
  println("      --metaphlan_db_dir [PATH] Path to metaphlan database. Example, /path/to/Reference_DBs/metaphlan4-db/. Default: null.")
  println("      --chocophlan_dir [PATH] Path to Humann's chocophlan nucleotide database. Example, /path/to/Reference_DBs/humann3-db/chocophlan/. Default: null.")
  println("      --uniref_dir [PATH] Path to Humann's Uniref protein database. Example, /path/to/Reference_DBs/humann3-db/uniref/. Default: null.")
  println("      --utilities_dir [PATH] Path to Humann's untilities database. Example, /path/to/Reference_DBs/humann3-db/utility_mapping/.  Default: null.")
  println("GTDBTK database:")
  println("      --GTDBTK_LINK [URL] GTDBTK database online download link. Default: https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz.")
  println("      --gtdbtk_db_dir  [PATH] Path to GTDBTK database. Example, /path/Reference_DBs/GTDB-tk-ref-db/. Default: null.")
  println("kofam scan database database:")
  println("      --ko_db_dir  [PATH] Path to kofam scan database. Example, /path/to/Reference_DBs/kofamscan_db/. Default: null.")
  println()
  println("Paths to existing conda environments to use, otherwise, new ones will be created using the yaml files in envs/.")
  println("      --conda_qc [PATH] Path to a conda environment containing fastqc, multiqc, zip and python. Default: null.")
  println("      --conda_humann3 [PATH] Path to a conda environment with humann3 installed. Default: null.")
  println("      --conda_cat  [PATH] Path to a conda environment containing CAT (Contig annotation tool). Default: null.")
  println("      --conda_prodigal [PATH] Path to a conda environment with prodigal installed. Default: null.")
  println("      --conda_metabat [PATH] Path to a conda environment containing metabat. Default: null.")
  println("      --conda_gtdbtk [PATH] Path to a conda environment containing gtdbtk. Default: null.")
  println("      --conda_kegg_decoder [PATH] Path to a conda environment with kegg_decoder installed. Default: null.")
  println("      --conda_megahit  [PATH] Path to a conda environment containing megahit. Default: null.")
  println("      --conda_bit [PATH] Path to a conda environment with bit installed. Default: null.")
  println("      --conda_kofamscan [PATH] Path to a conda environment containing KOFAM SCAN. Default: null.")
  println("      --conda_mapping [PATH] Path to a conda environment with bowtie and samtools installed. Default: null.")
  println("      --conda_checkm [PATH] Path to a conda environment with checkm installed. Default: null.")
  println()
  print("Advanced users can edit the nextflow.config file for more control over default settings such container choice, number of cpus, memory per task etc.")
  exit 0
  }

/************************************************
*********** Show pipeline parameters ************
*************************************************/

if (params.debug) {
log.info """${c_blue}
         Nextflow Metagenomics Illumina Consensus Pipeline: $workflow.manifest.version
         
         You have set the following parameters:
         Profile: ${workflow.profile} 
         Input csv file : ${params.input_file}
         GLDS or OSD Accession : ${params.accession}
         GLDS Raw File Pattern: ${params.RawFilePattern}         
         Workflow : ${params.workflow}
         Nextflow Directory publishing mode: ${params.publishDir_mode}
         Nextflow Error strategy: ${params.errorStrategy}
         Use GTDBTK Scratch Location: ${params.use_gtdbtk_scratch_location}
         MultiQC configuration file: ${params.multiqc_config}
         Megahit Maximum Memory: ${params.max_mem}
         Pile-up Memory: ${params.pileup_mem}
         CAT block size: ${params.block_size}

         MAG Parameters:
         Minimum completion: ${params.min_est_comp}
         Maximum redundancy: ${params.max_est_redund}
         Maximum strain heterogeneity: ${params.max_est_strain_het}
         Use Reduced Tree: ${params.reduced_tree}
 
         Output Directories:
         Raw reads: ${params.merged_dir}
         Filtered Reads: ${params.filtered_reads_dir}
         Assembly-based Analysis: ${params.assembly_based_dir}
         Read-based Analysis: ${params.read_based_dir}

         Genelab Assay Suffix: ${params.assay_suffix}
         Additional Filename Prefix: ${params.additional_filename_prefix}

         Conda Environments:
         humann3: ${params.conda_humann3}
         CAT: ${params.conda_cat}
         prodigal: ${params.conda_prodigal}
         metabat: ${params.conda_metabat}
         gtdbtk: ${params.conda_gtdbtk}
         kegg decoder: ${params.conda_kegg_decoder}
         megahit: ${params.conda_megahit}
         bit: ${params.conda_bit}
         kofamscan: ${params.conda_kofamscan}
         checkm: ${params.conda_checkm}         

         Databases:
         CAT Taxonomy: ${params.cat_taxonomy_dir}
         CAT DB sub directory: ${params.cat_db_sub_dir}
         CAT URL: ${params.CAT_DB_LINK}
         CAT DB: ${params.cat_db}
         KOFAM Scan: ${params.ko_db_dir}
         Metaphlan: ${params.metaphlan_db_dir}
         Chocophlan: ${params.chocophlan_dir}
         Uniref: ${params.uniref_dir}
         Utilities: ${params.utilities_dir}
         GTDBTK URL: ${params.GTDBTK_LINK}
         GTDBTK DB: ${params.gtdbtk_db_dir}
         ${c_reset}"""
}

// Illumina workflow
include { illumina } from "./workflows/illumina.nf"

// Nanopore workflow
include { nanopore } from "./workflows/nanopore.nf"


// Read-based workflow
include { read_based } from "./modules/read_based_processing.nf"

// Assembly-based workflow
include { assembly_based } from "./modules/assembly_based_processing.nf"


// Workflow to perform read-based analysis
workflow run_read_based_analysis {

    take:
        reads_per_sample
        metadata
        filtered_ch

    main:

        software_versions_ch = Channel.empty()    
        read_based(reads_per_sample, metadata, filtered_ch, 
                    params.krakendb_dir,
                    params.kaijudb_dir,
                    params.chocophlan_dir,
                    params.uniref_dir,
                    params.metaphlan_db_dir,
                    params.utilities_dir)

         read_based.out.versions  | mix(software_versions_ch) | set{software_versions_ch}

    emit:
        versions =  software_versions_ch

}


// Workflow to perform assembly-based analysis
workflow run_assembly_based_analysis {

    take:
        metadata
        file_ch
        filtered_ch


    main:
        software_versions_ch = Channel.empty()

        kofam_db = params.ko_db_dir
        cat_db = params.cat_db
        gtdbtk_db_dir = params.gtdbtk_db_dir

        // Run assembly based workflow 
        assembly_based(metadata, file_ch, filtered_ch, kofam_db, 
                        cat_db, gtdbtk_db_dir, params.use_gtdbtk_scratch_location)


        assembly_based.out.versions | mix(software_versions_ch) | set{software_versions_ch}


    emit:
        versions =  software_versions_ch

}

// A function to delete white spaces from an input string and convert it to lower case 
def deleteWS(string){

    return string.replaceAll(/\s+/, '').toLowerCase()

}



workflow {


    // Sanity check : Test input requirement
    if (!params.accession &&  !params.input_file){
     
       error("""${c_back_bright_red}INPUT ERROR!
              Please supply either an accession (OSD or Genelab number) or an input CSV file
              by passing either to the --accession or --input_file parameter, respectively.
              ${c_reset}""")
    } 
        
     // Software Version Capturing - runsheet
     software_versions_ch = Channel.empty()

     // Parse file input
       if(params.accession){

       GET_RUNSHEET(params.accession)
       GET_RUNSHEET.out.input_file
           .splitCsv(header:true)
           .set{file_ch}

       GET_RUNSHEET.out.version | mix(software_versions_ch) | set{software_versions_ch}
      }else{
 
       Channel.fromPath(params.input_file, checkIfExists: true)
           .splitCsv(header:true)
           .set{file_ch}
      }



    if( params.technology == "illumina" ) {

         illumina(params.sample_type, file_ch)

         clean_reads      = illumina.out.clean_reads
         reads_per_sample = illumina.out.reads_per_sample
         metadata         = illumina.out.metadata
         illumina.out.software_versions | mix(software_versions_ch) | set{software_versions_ch}


    }else if(  params.technology == "nanopore" ) {


         nanopore(params.input_type, params.sample_type, file_ch)
 
         clean_reads      = nanopore.out.clean_reads
         reads_per_sample = nanopore.out.reads_per_sample
         metadata         = nanopore.out.metadata
         nanopore.out.software_versions | mix(software_versions_ch) | set{software_versions_ch}

    }else{

        error("""${c_back_bright_red}INPUT ERROR!
              You must specify a recognized technology by passing one of
              'illumina' or  'nanopore' to  the --technology parameter.
              'illumina' - Illumina short reads.
              'nanopore' - Oxford nanopore long reads.
              ${c_reset}""")

    }


    // Run the analysis based on selection i.e, read-based, assembly-based or both
    // it will run both by default
    if(params.workflow == 'read-based'){
           
          run_read_based_analysis(reads_per_sample, metadata, clean_reads)
          
          run_read_based_analysis.out.versions | mix(software_versions_ch) | set{software_versions_ch}
          
    }else if(params.workflow == 'assembly-based') {

          run_assembly_based_analysis(metadata, file_ch, clean_reads)
          run_assembly_based_analysis.out.versions | mix(software_versions_ch) | set{software_versions_ch}

    }else{

          run_read_based_analysis(reads_per_sample, metadata, clean_reads)
          run_assembly_based_analysis(metadata, file_ch, clean_reads)

          run_read_based_analysis.out.versions | mix(software_versions_ch) | set{software_versions_ch}
          run_assembly_based_analysis.out.versions | mix(software_versions_ch) | set{software_versions_ch}
    }


     // Software Version Capturing - combining all captured sofware versions
     nf_version = "Nextflow Version ".concat("${nextflow.version}")
     nextflow_version_ch = Channel.value(nf_version)
     workflow_version = "Metagenomics ".concat("${workflow.manifest.version}")
     workflow_version_ch =  Channel.value(workflow_version)

     //  Write software versions to file
     software_versions_ch | map { it.text.strip() }
                          | unique
                          | mix(nextflow_version_ch)
                          | mix(workflow_version_ch)
                          | collectFile(name: "${params.metadata_dir}/software_versions.txt", newLine: true, cache: false)
                          | set{final_software_versions_ch}

}

workflow.onComplete {

    println("${c_bright_green}Pipeline completed at: $workflow.complete")
    println("""Execution status: ${ workflow.success ? 'OK' : "${c_back_bright_red}failed" }""")
    log.info ( workflow.success ? "\nDone! Workflow completed without any error\n" : "Oops .. something went wrong${c_reset}" )

    if ( workflow.success ) {

    println("Merged/Raw outputs location: ${params.merged_dir}")
    println("Filtered outputs location: ${params.filtered_dir}")
    println("Read-based Analysis: ${params.read_based_dir}")
    println("Assembly-based Analysis: ${params.assembly_based_dir}")
    println("Software versions location: ${params.metadata_dir}")
    println("Pipeline tracing/visualization files location:  ../Resource_Usage${c_reset}")
    println()
    }

}

