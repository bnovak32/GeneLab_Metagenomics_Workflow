#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Create GLDS runsheet
include { GET_RUNSHEET } from "../modules/create_runsheet.nf"

// Demultiplexing
include {DORADO_BASECALLER; DORADO_DEMUX} from "../modules/demultiplexing.nf"
// Split fastq concatenation
include {CAT_FASTQ_FILES; CAT_FASTQ_DIR} from "../modules/demultiplexing.nf"

// Read quality check and filtering
include { nano_quality_check as raw_qc; FILTLONG ; PORECHOP} from "../modules/quality_assessment.nf"
include { nano_quality_check as filtered_qc } from "../modules/quality_assessment.nf"
include { nano_quality_check as trimmed_qc } from "../modules/quality_assessment.nf"
include { nano_quality_check as nohost_qc } from "../modules/quality_assessment.nf"

// Remove human reads
include { remove_host as remove_human } from "../modules/remove_host.nf"
include { nano_quality_check as nohum_qc } from "../modules/quality_assessment.nf"

if(params.sample_type == "low_biomass"){

   // Remove contaminant
   include { nano_remove_contaminants as remove_contaminants } from "../modules/remove_contaminant.nf"
   include { nano_quality_check as noblank_qc } from "../modules/quality_assessment.nf"
}

// Remove host
include { remove_host } from "../modules/remove_host.nf"

// Custom genome mapping
include { LONG_MAP2GENOME  as FILTERED_MAP2GENOME} from "../modules/genome_mapping.nf"
include { LONG_MAP2GENOME  as DECONTAMED_MAP2GENOME} from "../modules/genome_mapping.nf"

// A function to delete white spaces from an input string and convert it to lower case 
def deleteWS(string){

    return string.replaceAll(/\s+/, '').toLowerCase()

}

// Main workflow
workflow nanopore {


     take:
      input_type
      sample_type
      file_ch



     main:

     // Software Version Capturing - runsheet
     software_versions_ch = Channel.empty()

    //---------------------------- Parse input file -------------------------------------------//
     if(input_type == 'single'){

        // One fastq file per sample
        file_ch.map{
                row -> tuple("${row.sample_id}", [file("${row.forward}", checkIfExists: true)], deleteWS(row.paired))
                }.set{reads_ch}
     
     }else if(input_type == 'multiple'){

         // Multiple fastq files per sample i.e Fastq files have been split but not concatenated per sample

        // Read multiple fastq files per sample to keep unique sample ids and a list of fastq files per sample
        file_ch.map{ row ->
            def meta = [id: row.sample_id]
            [meta, [file(row.forward, checkIfExists: true)]]
        }.groupTuple() // group by map ids i.e sample_id
         .map { meta, reads -> [ meta.id, reads.flatten() ] } // [sample_id, [reads]]
         .set{read_ch}
         // Concatenate multiple fastq to one fastq per sample
         CAT_FASTQ_FILES(read_ch)


        if(sample_type == "low_biomass"){
         // Get distinct sample metadata
         file_ch.map{
                row -> tuple( "${row.sample_id}", deleteWS(row.group),
                                deleteWS(row.NTC), deleteWS(row.concentration),
                                deleteWS(row.paired) )
                }
                .distinct()
                .set{meta_ch}
          }else{
         file_ch.map{
                row -> tuple( "${row.sample_id}", deleteWS(row.group), deleteWS(row.paired))
                }
                .distinct()
                .set{meta_ch}

          }
   


      CAT_FASTQ_FILES.out.reads
            .join(meta_ch)
            .set{concated_files_ch}


        if(sample_type == "low_biomass"){

        header = Channel.of(["sample_id", "forward", "group", "NTC", "concentration", "paired"]) 

        header.concat(concated_files_ch)
              .map{ sample_id, forward, group, NTC, concentration, paired ->
                   "${sample_id},${forward},${group},${deleteWS(NTC)},${concentration},${paired}"
              }
              .collectFile(name: "${launchDir}/samples_file.csv", newLine: true, sort:false)
              .splitCsv(header:true)
              .set{InFile_ch}

        }else{

        header = Channel.of(["sample_id", "forward", "group", "paired"])

        header.concat(concated_files_ch)
              .map{ sample_id, forward, group, paired ->
                   "${sample_id},${forward},${group},${paired}"
              }
              .collectFile(name: "${launchDir}/samples_file.csv", newLine: true, sort:false)
              .splitCsv(header:true)
              .set{InFile_ch}


        }

        InFile_ch.map{ row -> 
                    tuple( row.sample_id, [file(row.forward, checkIfExists: true)], deleteWS(row.paired))
                }.set{reads_ch}

        CAT_FASTQ_FILES.out.version | mix(software_versions_ch) | set{software_versions_ch}

     }else if(input_type == 'directory'){

         // Directory containing POD5 files
         pod5_dir   = Channel.fromPath(params.input_dir, checkIfExists: true)
 
         // Basecall, demultiplex and concatenate demultiplexed fastq files
         DORADO_BASECALLER(pod5_dir, params.kit_name)
         DORADO_DEMUX(DORADO_BASECALLER.out.bam, params.kit_name)

         // Create 2-column file to rename barcode names to sample names
         // sample_id to barcode_id column in --input_file
        
        file_ch.map{  row ->
                   "${row.sample_id},${row.barcode_id}"
              }
              .collectFile(name: "sample2barcode_file.csv", newLine: true, sort:false)
              .set{sample2barcode}

         CAT_FASTQ_DIR(sample2barcode, DORADO_DEMUX.out.demux_dir)
         
        // Read-in runsheet generated fromm concatenating fastq files above
        CAT_FASTQ_DIR.out.runsheet.splitCsv(header:true)
           .map{
                row -> tuple( "${row.sample_id}", [file("${row.forward}", checkIfExists: true)])
                }.set{runsheet_ch}
         
         if(sample_type == "low_biomass"){

         // Read input file into tuples
         file_ch.map{
                row -> tuple( "${row.sample_id}", deleteWS(row.group),
                                deleteWS(row.NTC), deleteWS(row.concentration),
                                deleteWS(row.paired) )
                }.set{InFile_ch}

        // Merge the generated fastq files with their corresponding metadata
        runsheet_ch.join(InFile_ch)
                    .map{sample_id, forward, group, NTC, concentration, paired -> 
                    tuple(sample_id, forward, paired)
                }.set{reads_ch}

        }else{

         // Read input file into tuples
         file_ch.map{
                row -> tuple( "${row.sample_id}", deleteWS(row.group), deleteWS(row.paired) )
                }.set{InFile_ch}

        // Merge the generated fastq files with their corresponding metadata
        runsheet_ch.join(InFile_ch)
                    .map{sample_id, forward, group, paired ->
                    tuple(sample_id, forward, paired)
                }.set{reads_ch}


        }

        DORADO_BASECALLER.out.version | mix(software_versions_ch) | set{software_versions_ch}
        DORADO_DEMUX.out.version | mix(software_versions_ch) | set{software_versions_ch}
        CAT_FASTQ_DIR.out.version | mix(software_versions_ch) | set{software_versions_ch}

     }else{

        error("""${c_back_bright_red}INPUT ERROR!
              You must specify a recognized input type by passing one of
              'single', 'multiple' or 'directory' to --input_type parameter.
              'single' -  On fastq file per sample.
              'multiple' - Multiple fastq files per sample.
              'directory' - Directory containing POD5 files.
              ${c_reset}""")

     }
    
    // Quality check  input reads
    raw_qc(Channel.of("raw"), params.multiqc_config,reads_ch,Channel.empty())

    // Filter input reads based on length and quality then quality check the filtered reads
    FILTLONG(reads_ch)
    filtered_qc(Channel.of("filtered"), params.multiqc_config, FILTLONG.out.reads, FILTLONG.out.log)

    // Trim off primers and adapters using porechop
    PORECHOP(FILTLONG.out.reads)
    trimmed_ch = PORECHOP.out.reads
    trimmed_qc(Channel.of("trimmed"), params.multiqc_config, trimmed_ch, PORECHOP.out.log)
    // Map trimmed reads to a custom genome with bbmap
    //FILTERED_MAP2GENOME(params.custom_genome, Channel.of("trimmed"), trimmed_ch)
    
    // Quality check software capturing
    raw_qc.out.versions | mix(software_versions_ch) | set{software_versions_ch}
    //FILTERED_MAP2GENOME.out.version | mix(software_versions_ch) | set{software_versions_ch}
    FILTLONG.out.version | mix(software_versions_ch) | set{software_versions_ch}
    filtered_qc.out.versions | mix(software_versions_ch) | set{software_versions_ch}
    trimmed_qc.out.versions | mix(software_versions_ch) | set{software_versions_ch}

    // Get the number of reads per sample after read filtering
    // relative abundance to count calculation for metaphlan results
    reads_per_sample = trimmed_qc.out.reads_per_sample

    // Remove human reads and quality check
    remove_human("HRrm", "human", params.human_db_url, null, params.human_db_dir, trimmed_ch)
    nohum_qc(Channel.of("HRrm"), params.multiqc_config, remove_human.out.clean_reads,remove_human.out.logs)
    remove_human.out.versions | mix(software_versions_ch) | set{software_versions_ch}

    // By default trimmed reads are reads after quality filtering, trimming with filtlong 
    // and porechop and human reads removed with kraken2
    trimmed_reads = remove_human.out.clean_reads

    if(sample_type == "low_biomass"){

    // Remove contaminants and quality check
    remove_contaminants(file_ch, remove_human.out.clean_reads)
    noblank_qc(Channel.of("decontam"), params.multiqc_config,
               remove_contaminants.out.clean_reads, remove_contaminants.out.logs)
    remove_contaminants.out.versions | mix(software_versions_ch) | set{software_versions_ch}
    noblank_qc.out.versions | mix(software_versions_ch) | set{software_versions_ch}

   
    // Prepare metadata
    meta_header = Channel.of(["sample_id", "group", "NTC", "concentration"])
    file_ch.map{
                row -> tuple( "${row.sample_id}", row.group, deleteWS(row.NTC), row.concentration )
                }
                .distinct()
                .set{body}

    meta_header.concat(body)
              .map{ sample_id, group, NTC, concentration ->
                   "${sample_id},${group},${NTC},${concentration}"
              }
              .collectFile(name: "metadata_file.txt", newLine: true, sort:false)
              .set{metadata}

   // If the dataset is a low biomass dataset then the trimmed reads
   // are reads after contaminants (negative control) have been removed
   trimmed_reads = remove_contaminants.out.clean_reads

   // Get the number of reads per sample after removing contaminants for
   // relative abundance to count calculation for metaphlan results
   reads_per_sample = noblank_qc.out.reads_per_sample

   }else{
 
    // Prepare metadata
    meta_header = Channel.of(["sample_id", "group"])
    file_ch.map{
                row -> tuple( "${row.sample_id}", row.group )
                }
                .distinct()
                .set{body}

    meta_header.concat(body)
              .map{ sample_id, group->
                   "${sample_id},${group}"
              }
              .collectFile(name: "metadata_file.txt", newLine: true, sort:false)
              .set{metadata}

  }



   // If host database or host link or host fats or host name is provided then remove 
   // host sequences
   if( params.host_db_dir || params.host_url || params.host_fasta || params.host_name ){

    // Remove host
    remove_host("HostRm", params.host_name, params.host_url, params.host_fasta,
                params.host_db_dir, trimmed_reads)

    //DECONTAMED_MAP2GENOME(params.custom_genome, Channel.of("decontamed"), clean_reads)
    nohost_qc(Channel.of("HostRm"), params.multiqc_config, remove_host.out.clean_reads, remove_host.out.logs)

    //DECONTAMED_MAP2GENOME.out.version | mix(software_versions_ch) | set{software_versions_ch}
    remove_host.out.versions | mix(software_versions_ch) | set{software_versions_ch}

    clean_reads = remove_host.out.clean_reads
    // Get the number of reads per sample after removing contaminants and host reads for
    // relative abundance to count calculation for metaphlan results
    reads_per_sample = nohost_qc.out.reads_per_sample

   }else{

    clean_reads = trimmed_reads

   }

     emit:
         clean_reads       = clean_reads
         reads_per_sample  = reads_per_sample
         metadata          = metadata
         software_versions = software_versions_ch

}



workflow {


    Channel.fromPath(params.input_file, checkIfExists: true)
           .splitCsv(header:true)
           .set{file_ch}

    nanopore(params.input_type, params.sample_type, file_ch)
}
