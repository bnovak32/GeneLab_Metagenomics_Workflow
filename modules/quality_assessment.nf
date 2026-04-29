#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Sequence quality assessment and control processes ****************
****************************************************************************************/
//params.multiqc_config = "config/multiqc.config"

/*
 * ========================================================================================
 * PROCESS: FASTQC
 * ========================================================================================
 *
 * SUMMARY:
 *   Sample quality assessment with fastqc
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Bolean specifying whether input reads are paired or not  
 *
 * OUTPUTS:
 *   1. tuple: tuple path("*.html"), path("*.zip") (emit: html)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: FastQC
 *   Container: [Defined in config/illumina.config]
 *   Conda: envs/fastqc.yaml
 *   Labels: quality_check
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process FASTQC {

  // FastQC performed on reads
  tag "Running fastqc on ${sample_id}"
  label "quality_check"
  beforeScript "chmod +x ${projectDir}/bin/*"

  input:
    tuple val(sample_id), path(reads), val(isPaired)
  output:
    tuple path("*.html"), path("*.zip"), emit: html
    path("versions.txt"), emit: version
  script:
    """
    fastqc -o . \\
     -t ${task.cpus} \\
      ${reads}

    fastqc --version > versions.txt
    """
}

/*
 * ========================================================================================
 * PROCESS: NANOPLOT
 * ========================================================================================
 *
 * SUMMARY:
 *   Sample quality assessment with nanoplot
 *
 * INPUTS:
 *   1. each: prefix
 *      Cardinality: each
 *      Description: Iterates over each element. Prefix to add to output file.
 *
 *   2. tuple: tuple val(sample_id), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *
 * OUTPUTS:
 *   1. tuple: tuple path("*.log"), path("*_NanoStats.txt"), path("*.html") (emit: html)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: NanoPlot
 *   Container: [Defined in config/nanopore.config]
 *   Conda: envs/nanoplot.yaml
 *   Labels: quality_check
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// Quality control on Nanopore data
process NANOPLOT {

    tag "Running nanoplot on ${sample_id}"
    label "quality_check"
    beforeScript "chmod +x ${projectDir}/bin/*"

    input:
        each prefix
        tuple val(sample_id), path(reads), val(isPaired)

    output:
        tuple path("*.log"), path("*_NanoStats.txt"), path("*.html"), emit: html
        path("versions.txt"), emit: version

    script:
    """
    NanoPlot \\
        --prefix ${sample_id}_${prefix}_ \\
        -t ${task.cpus} \\
        --fastq ${reads[0]} \\
        -o .
    # Rename html report
    mv ${sample_id}_${prefix}_NanoPlot-report.html \\
       ${sample_id}_${prefix}_NanoPlot-report${params.assay_suffix}.html

    VERSION=`NanoPlot --version 2>&1 | sed 's/^.*NanoPlot //; s/ .*\$//'`
    echo "NanPlot \${VERSION}"  > versions.txt
    """
}

/*
 * ========================================================================================
 * PROCESS: MULTIQC
 * ========================================================================================
 *
 * SUMMARY:
 *   Summarize QC outputs with multiqc
 *
 * INPUTS:
 *   1. val: prefix
 *      Cardinality: one
 *      Description: Parameter value: prefix to add to output file name
 *
 *   2. path: multiqc_config
 *      Cardinality: one
 *      Description: Input file: multiqc config file to customize output
 *
 *   3. path: files
 *      Cardinality: one
 *      Description: Input file: QC files to summarize
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}${prefix}_multiqc${params.assay_suffix}.html (emit: html)
 *
 *   2. path: ${params.additional_filename_prefix}${prefix}_multiqc_report/${prefix}_multiqc_data/ (emit: data)
 *
 *   3. path: ${params.additional_filename_prefix}${prefix}_reads_per_sample.tsv (emit: reads_per_sample)
 *
 *   4. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: MultiQC
 *   Container: [Defined in config/default.config]
 *   Conda: envs/multiqc.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process MULTIQC {

  tag "Running multiqc on the ${prefix} files.."

  input:
    val(prefix)   
    path(multiqc_config)
    path(files)
  output:
    path("${params.additional_filename_prefix}${prefix}_multiqc${params.assay_suffix}.html"), emit: html
    path("${params.additional_filename_prefix}${prefix}_multiqc_report/${prefix}_multiqc_data/"), emit: data
    path("${params.additional_filename_prefix}${prefix}_reads_per_sample.tsv"), emit: reads_per_sample
    path("versions.txt"), emit: version
  script:
    """
      multiqc -q --filename ${params.additional_filename_prefix}${prefix}_multiqc \\
              --force --cl-config 'max_table_rows: 99999999' \\
              --interactive --config ${multiqc_config} \\
              --outdir ${params.additional_filename_prefix}${prefix}_multiqc_report  ${files} > /dev/null 2>&1


      # Copy and rename html file 
      cp ${params.additional_filename_prefix}${prefix}_multiqc_report/${prefix}_multiqc.html  \\
          ${params.additional_filename_prefix}${prefix}_multiqc${params.assay_suffix}.html
  

       
      if [ `find -type f  -name 'multiqc_nanostat.txt' | wc -l` -gt 0 ]; then

      # Nanopore dataset - Nanoplot

      FILENAME=`find  -type f -name multiqc_nanostat.txt`
      # Write out the number of reads per sample to file
      awk 'BEGIN{print "Sample_ID\\tReads"} NR>1{printf "%s\\t%s\\n", \$1,\$6 }' \${FILENAME} \\
           > ${params.additional_filename_prefix}${prefix}_reads_per_sample.tsv
    
      else

      # Illumina dataset - fastqc
      FILENAME=`find  -type f -name multiqc_general_stats.txt`
      # Write out the number of reads per sample to file
      awk 'BEGIN{print "Sample_ID\\tReads"} NR>1{printf "%s\\t%s\\n", \$1,\$NF }' \${FILENAME} \\
           > ${params.additional_filename_prefix}${prefix}_reads_per_sample.tsv
          
      fi  

      multiqc --version > versions.txt
    """
  }

/*
 * ========================================================================================
 * PROCESS: ZIP_MULTIQC
 * ========================================================================================
 *
 * SUMMARY:
 *   Zip directory
 *
 * INPUTS:
 *   1. val: prefix
 *      Cardinality: one
 *      Description: Parameter value: prefix to add to output zip file.
 *
 *   2. path: multiqc_dir
 *      Cardinality: one
 *      Description: Input file: directory to zip
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}${prefix}_multiqc${params.assay_suffix}_data.zip (emit: data)
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: MultiQC
 *   Container: [Defined in config/default.config]
 *   Conda: envs/zip.yaml
 *   Labels: zip
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process ZIP_MULTIQC {

    tag "Zipping ${prefix} multiqc.."
    label "zip"
 
    input:
        val(prefix)
        path(multiqc_dir)

    output:
        path("${params.additional_filename_prefix}${prefix}_multiqc${params.assay_suffix}_data.zip"), emit: data
        path("versions.txt"), emit: version

    script:
        """
        # zipping and removing unzipped dir
        zip -q -r \\
           ${params.additional_filename_prefix}${prefix}_multiqc${params.assay_suffix}_data.zip \\
           ${multiqc_dir}

        zip -h | grep "Zip" | sed -E 's/(Zip.+\\)).+/\\1/' > versions.txt
        """
}



/*
 * ========================================================================================
 * PROCESS: FASTP
 * ========================================================================================
 *
 * SUMMARY:
 *   Quality filter sample reads with fastp
 *
 * INPUTS:
 *   1. each: trimPolyG
 *      Cardinality: each
 *      Description: Iterates over each element
 *
 *   2. tuple: tuple val(sample_id), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Bolean specifying whether input reads are paired or not 
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("*fastq.gz"), val(isPaired) (emit: reads)
 *
 *   2. tuple: tuple val(sample_id), path("${sample_id}.fastp.json") (emit: json)
 *
 *   3. tuple: tuple val(sample_id), path("${sample_id}.fastp.html") (emit: html)
 *
 *   4. tuple: tuple val(sample_id), path("${sample_id}-fastp.log") (emit: log)
 *
 *   5. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/illumina.config]
 *   Conda: envs/fastp.yaml
 *   Labels: fastp
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

//Adapter trimming and filtering for short reads 
// (quality_score>20;min_length=20;low complexity filter;threading;
// automatically detect adapters;json report)
process FASTP {

    tag "Quality filtering ${sample_id}-s reads.."
    beforeScript "chmod +x ${projectDir}/bin/*"
    label "fastp"
    
    input:
    each trimPolyG
    tuple val(sample_id), path(reads), val(isPaired)
    
    output:
    tuple val(sample_id), path("*fastq.gz"), val(isPaired), emit: reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"), emit: html
    tuple val(sample_id), path("${sample_id}-fastp.log"), emit: log
    path("versions.txt"), emit: version   
    
    script:
        def polyG = trimPolyG == 'true' ? "--trim_poly_g": "-G"
        def out_prefix = trimPolyG == 'true' ? "" : "temp_"
    """
    if [ ${isPaired} == true ]; then
    
        fastp --in1 ${reads[0]} --out1 ${out_prefix}${sample_id}_R1_filtered${params.assay_suffix}.fastq.gz \\
          --in2 ${reads[1]} --out2 ${out_prefix}${sample_id}_R2_filtered${params.assay_suffix}.fastq.gz \\
          --qualified_quality_phred  20 \\
          --length_required 50 \\
          --thread ${task.cpus} \\
          --detect_adapter_for_pe \\
          --json ${sample_id}.fastp.json \\
          --html ${sample_id}.fastp.html \\
          --trim_front1 0 \\
          --trim_tail1 0 \\
          --trim_front2 0 \\
          --trim_tail2 0  ${polyG}  2> ${sample_id}-fastp.log
   
    else

        fastp --in1 ${reads[0]} --out1 ${out_prefix}${sample_id}_filtered${params.assay_suffix}.fastq.gz \\
          --qualified_quality_phred  20 \\
          --length_required 50 \\
          --thread ${task.cpus} \\
          --json ${sample_id}.fastp.json \\
          --html ${sample_id}.fastp.html \\
          --trim_front1 0 \\
          --trim_tail1 0 ${polyG}  2> ${sample_id}-fastp.log

    fi

    VERSION=`fastp --version 2>&1 | sed -e "s/fastp //g"`
    echo "fastp \${VERSION}" > versions.txt
    """

}

/*
 * ========================================================================================
 * PROCESS: FILTLONG
 * ========================================================================================
 *
 * SUMMARY:
 *   Quality filter sample reads by length with filtlong
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Bolean specifying whether input reads are paired or not  
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), val(isPaired) (emit: reads)
 *
 *   2. tuple: tuple val(sample_id), path("${sample_id}-filtlong.log") (emit: log)
 *
 *   3. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/nanopore.config]
 *   Conda: envs/filtlong.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process FILTLONG {

    tag "Quality filtering ${sample_id}-s reads.."
    beforeScript "chmod +x ${projectDir}/bin/*"

    input:
        tuple val(sample_id), path(reads), val(isPaired)

    output:
        tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), val(isPaired), emit: reads
        tuple val(sample_id), path("${sample_id}-filtlong.log"), emit: log
        path("versions.txt"), emit: version

    script:
    """
    filtlong \\
        --min_length 200 \\
        --min_mean_q 8 \\
        ${reads[0]} 2> >(tee ${sample_id}-filtlong.log >&2) \\
        | gzip -n > ${sample_id}_filtered.fastq.gz

    VERSION=\$(filtlong --version | sed -e "s/Filtlong v//g")
    echo "filtlong \${VERSION}"  > versions.txt
    """
}


/*
 * ========================================================================================
 * PROCESS: PORECHOP
 * ========================================================================================
 *
 * SUMMARY:
 *   Trim adapters from sample reads with porechop
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(reads), val(isPaired)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - reads: path to sample fastq reads
 *                 - isPaired: Bolean specifying whether input reads are paired or not  
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), val(isPaired) (emit: reads)
 *
 *   2. tuple: tuple val(sample_id), path("${sample_id}-porechop.log") (emit: log)
 *
 *   3. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/nanopore.config]
 *   Conda: envs/porechop.yaml 
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// Adapter trimming for Nanopore
// Porechop is a tool for removing adapter sequences from nanopore reads
process PORECHOP {
    tag "Trimming ${sample_id}-s reads...."
    beforeScript "chmod +x ${projectDir}/bin/*"

    input:
        tuple val(sample_id), path(reads), val(isPaired)

    output:
        tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), val(isPaired), emit: reads
        tuple val(sample_id), path("${sample_id}-porechop.log"), emit: log
        path("versions.txt"), emit: version


    script:
    """
    porechop \\
        -i ${reads[0]} \\
        -t ${task.cpus} \\
        --discard_middle \\
        -o ${sample_id}_trimmed.fastq.gz \\
        > ${sample_id}-porechop.log

    VERSION=`porechop --version` 
    echo "porechop \${VERSION}" > versions.txt
    """
}


workflow nano_quality_check {

    take:
        prefix_ch
        multiqc_config
        reads_ch
        logs_ch


    main:
        NANOPLOT(prefix_ch,reads_ch)
        reads = reads_ch.map{ sample_id, reads, paired -> reads}
        logs = logs_ch.map{ sample_id, log -> log}.flatten()
        nanoplot_ch = NANOPLOT.out.html.flatten()
                              .mix(reads)
                              .mix(logs)
                              .flatten()
                              .collect()
        MULTIQC(prefix_ch, multiqc_config, nanoplot_ch)
        ZIP_MULTIQC(prefix_ch, MULTIQC.out.data)

        software_versions_ch = Channel.empty()
        NANOPLOT.out.version | mix(software_versions_ch) | set{software_versions_ch}
        MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
        ZIP_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}


    emit:
        reads_per_sample = MULTIQC.out.reads_per_sample
        versions = software_versions_ch
}



workflow quality_check {

    take:
        prefix_ch
        multiqc_config
        reads_ch
    

    main:
        FASTQC(reads_ch)
        fastqc_ch = FASTQC.out.html.flatten().collect()
        MULTIQC(prefix_ch, multiqc_config, fastqc_ch)
        ZIP_MULTIQC(prefix_ch, MULTIQC.out.data)

        software_versions_ch = Channel.empty()
        FASTQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
        MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
        ZIP_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}

    emit:
        reads_per_sample = MULTIQC.out.reads_per_sample 
        versions = software_versions_ch
}

workflow {

        Channel.fromPath(params.input_file)
               .splitCsv()
               .map{ row -> row.paired == 'true' ? tuple( "${row.sample_id}", [file("${row.forward}", checkIfExists: true), file("${row.reverse}", checkIfExists: true)], row.paired) : 
                                                   tuple( "${row.sample_id}", [file("${row.forward}", checkIfExists: true)], row.paired)}
               .set{reads_ch}   

    res_ch = quality_check(Channel.of(params.prefix), params.multiqc_config, reads_ch)
}
