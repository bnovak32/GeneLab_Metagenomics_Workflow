#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ========================================================================================
 * PROCESS: ZIP_FASTA
 * ========================================================================================
 *
 * SUMMARY:
 *   Zip bins or MAGs
 *
 * INPUTS:
 *   1. val: TYPE
 *      Cardinality: one
 *      Description: Parameter value: bin or MAG
 *
 *   2. path: DIR
 *      Cardinality: one
 *      Description: Input file: Directory to be zipped
 *
 * OUTPUTS:
 *   1. path: *.zip (emit: zip_files) [OPTIONAL]
 *
 *   2. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
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

process ZIP_FASTA {

    tag "Zipping up your ${TYPE}s..."
    label "zip"


    input:
        val(TYPE)
        path(DIR)

    output:
        path("*.zip"), emit: zip_files, optional: true
        path("versions.txt"), emit: version

    script:
        """
        function zip_sample() {

            local SAMPLE=\$1
            local TYPE=\$2

            mkdir -p \${SAMPLE}-\${TYPE}s && \\
            cp -f \${SAMPLE}-\${TYPE}*.fasta \${SAMPLE}-\${TYPE}s && \\
            zip -r \${SAMPLE}-\${TYPE}s${params.assay_suffix}.zip \${SAMPLE}-\${TYPE}s

           }


        export -f zip_sample

        if [ ${TYPE} == 'bin' ]; then 
       
              WORKDIR=`pwd`
        else

              WORKDIR=${DIR}
        fi

        if [ `find -L \${WORKDIR} -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then


             if [ ${TYPE} == 'MAG' ]; then

                 find -L \${WORKDIR} -name '*.fasta' | xargs -I {} cp {} .

             fi

             SAMPLES=(`ls -1 *.fasta | sed -E 's/(.+)-${TYPE}.*.fasta/\\1/g'`)

             for SAMPLE in \${SAMPLES[*]};do
                  
                  zip_sample \${SAMPLE} ${TYPE}

             done

        fi

        zip -h | grep "Zip" | sed -E 's/(Zip.+\\)).+/\\1/' > versions.txt
        """

}
