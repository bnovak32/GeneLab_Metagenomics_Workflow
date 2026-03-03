#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Create GLDS runsheet
include { GET_RUNSHEET } from "../modules/create_runsheet.nf"

// Read quality check and filtering
include { quality_check as raw_qc; FASTP } from "../modules/quality_assessment.nf"
include { FASTP as POLYG_FASTP } from "../modules/quality_assessment.nf"
include { MULTIQC as FASTP_MULTIQC } from "../modules/quality_assessment.nf" 
include { quality_check as filtered_qc } from "../modules/quality_assessment.nf"
include { quality_check as noblank_qc } from "../modules/quality_assessment.nf"
include { quality_check as nohost_qc } from "../modules/quality_assessment.nf"

if(params.sample_type == "low_biomass"){

   // Remove contaminant
   include { illumina_remove_contaminants as remove_contaminants } from "../modules/remove_contaminant.nf"
}

// Remove host
include { remove_host } from "../modules/remove_host.nf"

// Custome genome mapping
include { SHORT_MAP2GENOME  as NODECONTAM_MAP2GENOME } from "../modules/genome_mapping.nf"
include { SHORT_MAP2GENOME  as FILTERED_MAP2GENOME } from "../modules/genome_mapping.nf"
include { SHORT_MAP2GENOME  as DECONTAMED_MAP2GENOME } from "../modules/genome_mapping.nf"


// A function to delete white spaces from an input string and covert it to lower case 
def deleteWS(string){

    return string.replaceAll(/\s+/, '').toLowerCase()

}



// Main workflow
workflow illumina {
   take:
     sample_type
     file_ch


    main:

     // Software Version Capturing - runsheet
     software_versions_ch = Channel.empty()

    file_ch.map{
                row -> deleteWS(row.paired) == 'true'  ? tuple( "${row.sample_id}", [file("${row.forward}", checkIfExists: true), file("${row.reverse}", checkIfExists: true)], deleteWS(row.paired)) : 
                                                         tuple( "${row.sample_id}", [file("${row.forward}", checkIfExists: true)], deleteWS(row.paired))
                }.set{reads_ch}



    // Quality check and trim the input reads
    raw_qc(Channel.of("raw"), params.multiqc_config,reads_ch)

    //NODECONTAM_MAP2GENOME(params.custome_genome, Channel.of("no_decontam"), reads_ch)
    FASTP(Channel.of('false'), reads_ch) // no ployG trimming
    POLYG_FASTP(Channel.of('true'), FASTP.out.reads) // ployg trimming
    filtered_ch = POLYG_FASTP.out.reads
    json_ch = POLYG_FASTP.out.json.map{sample_id, json -> json}.collect()
    //FILTERED_MAP2GENOME(params.custome_genome, Channel.of("filtered"), filtered_ch)
    filtered_qc(Channel.of("filtered"), params.multiqc_config, filtered_ch)
    FASTP_MULTIQC(Channel.of('fastp'), params.multiqc_config, json_ch)

    // Quality check software capturing
    raw_qc.out.versions | mix(software_versions_ch) | set{software_versions_ch}
    //NODECONTAM_MAP2GENOME.out.version | mix(software_versions_ch) | set{software_versions_ch}
    //FILTERED_MAP2GENOME.out.version | mix(software_versions_ch) | set{software_versions_ch}
    filtered_qc.out.versions | mix(software_versions_ch) | set{software_versions_ch}

    // By default filtered reads are reads after quality filtering with fastp
    filtered_reads = filtered_ch
    // Get the number of reads per sample after read filtering
    // relative abundance to count calcultaion for metaphalan results
   reads_per_sample = filtered_qc.out.reads_per_sample


    if(sample_type == "low_biomass"){

    // Remove contaminants and quality ckeck
    remove_contaminants(file_ch, filtered_ch)
    noblank_qc(Channel.of("decontam"), params.multiqc_config,remove_contaminants.out.clean_reads)
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
 
   // If the dataset is a low biomass dataset then the filtered reads
   // are reads after contaminants (negative control) have been removed
   filtered_reads = remove_contaminants.out.clean_reads
   // Get the number of reads per sample after removing contaminats for
   // relative abundance to count calcultaion for metaphalan results
   reads_per_sample = noblank_qc.out.reads_per_sample


   } else{


    // Prepare metadata
    meta_header = Channel.of(["sample_id", "group"])
    file_ch.map{
                row -> tuple( "${row.sample_id}", row.group)
                }
                .distinct()
                .set{body}

    meta_header.concat(body)
              .map{ sample_id, group ->
                   "${sample_id},${group}"
              }
              .collectFile(name: "metadata_file.txt", newLine: true, sort:false)
              .set{metadata}
 
   }    


   // if host database or host link or host fats or host name is provided then remove 
   // host sequences
   if( params.host_db_dir || params.host_url || params.host_fasta || params.host_name ){

    // Remove host
    remove_host("HostRM", params.host_name, params.host_url, params.host_fasta,
                params.host_db_dir, filtered_reads)

    //DECONTAMED_MAP2GENOME(params.custome_genome, Channel.of("decontamed"), clean_reads)
    nohost_qc(Channel.of("HostRM"), params.multiqc_config, remove_host.out.clean_reads)

    //DECONTAMED_MAP2GENOME.out.version | mix(software_versions_ch) | set{software_versions_ch}
    remove_host.out.versions | mix(software_versions_ch) | set{software_versions_ch}

    clean_reads = remove_host.out.clean_reads

     // Get the number of reads per sample after removing contaminats and host reads for
     // relative abundance to count calcultaion for metaphalan results
     reads_per_sample = nohost_qc.out.reads_per_sample


   }else{

    clean_reads = filtered_reads


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

       illumina(params.sample_type, file_ch)


}
