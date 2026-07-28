#!/bin/bash

#SBATCH --job-name="nf_master" ## Replace job_name with the name of the job you are running ##
#SBATCH --output=nf_master.o.%j ## Replace job_name with the name of the job you are running ##
#SBATCH --error=nf_master.e.%j ## Replace job_name with the name of the job you are running ##
#SBATCH --partition=priority ## Specifies the job queue to use, for urgent jobs change normal to priority ##
#SBATCH --mem=20G ## Memory required to run the job in MB, this example is showing 10,000 MB or 10GB, change this number based on how much RAM you need ##
#SBATCH --cpus-per-task=1 ## Number of CPUs to run the job, this example is showing 5 CPUs, change this number based on how many CPUs you need ##
#SBATCH --mail-user=user@domain.com ## Specifies the e-mail address to e-mail when the job is complete, replace this e-mail address with your NASA e-mail address ##
#SBATCH --mail-type=END ## Tells slurm to e-mail the address above when the job has completed ##

. ~/.profile


echo "nf_master" ## Replace job_name with the name of the job you are running ##
echo ""


## Add a time-stamp at the start of the job ##
start=$(date +%s)
echo "start time: $start"

## Print the name of the compute node executing the job ##
echo $HOSTNAME

WORKFLOW_DIR='/path/to/workflow_code/'

# ---------------------------------------------------------------------------- Nanopore ---------------------------------------------------------------------------------------- #

# ---- Low biomass
# Processing
# Input dir
bash ./launch.sh processing ${WORKFLOW_DIR}/main.nf ${WORKFLOW_DIR}/nextflow.config '--sample_type "low_biomass" --input_type "directory"  --input_file input_dir_barcodes.csv --errorStrategy "ignore" --technology "nanopore" --assay_suffix "_GLlblMetag" --input_dir "/path/to/pod5/directory/" --kit_name "SQK-RPB114-24"'

# Single
bash ./launch.sh processing ${WORKFLOW_DIR}/main.nf ${WORKFLOW_DIR}/nextflow.config '--sample_type "low_biomass" --input_type "single"  --input_file single.csv  --errorStrategy "ignore" --technology "nanopore" --assay_suffix "_GLlblMetag"'

# Multiple
bash ./launch.sh processing ${WORKFLOW_DIR}/main.nf ${WORKFLOW_DIR}/nextflow.config '--sample_type "low_biomass" --input_type "multiple" --input_file multiple.csv  --errorStrategy "ignore" --technology "nanopore" --assay_suffix "_GLlblMetag"'

# Post-proscessing
bash ./launch.sh post_processing  ${WORKFLOW_DIR}/post_processing.nf  ${WORKFLOW_DIR}/config/post_processing.config \
           '--name "First M. Last" --email "username@nasa.com" --glds_accession "OSD-944" --osd_accession "OSD-944" 
           --assay_table assay_table.txt --runsheet "multiple.csv" 
           --technology "nanopore" --sample_type "low_biomass" --protocol_id "GL-DPPD-7116" --raw_suffix ".fastq.gz" 
           --raw_zip "../Merged_Sequence_Data/MultiQC_Reports/raw_multiqc_data.zip" 
           --filtered_zip "../Filtered_Sequence_Data/MultiQC_Reports/filtered_multiqc_data.zip" 
           --trimmed_zip "../Trimmed_Sequence_Data/MultiQC_Reports/trimmed_multiqc_data.zip" 
           --human_zip "../HR-removed_Sequence_Data/MultiQC_Reports/HRrm_multiqc_data.zip" 
           --decontam_zip "../Decontaminated_Sequence_Data/MultiQC_Reports/decontam_multiqc_data.zip" ' 

# ----- Standard
# Processing
# Input dir
bash ./launch.sh processing ${WORKFLOW_DIR}/main.nf ${WORKFLOW_DIR}/nextflow.config '--sample_type "standard" --input_type "directory"  --input_file input_dir_barcodes.csv --errorStrategy "ignore" --technology "nanopore" --assay_suffix "_GLlMetag" --input_dir "/path/to/pod5/directory/" --kit_name "SQK-RPB114-24"'

# Single
bash ./launch.sh processing ${WORKFLOW_DIR}/main.nf ${WORKFLOW_DIR}/nextflow.config '--sample_type "standard" --input_type "single"  --input_file single.csv  --errorStrategy "ignore" --technology "nanopore" --assay_suffix "_GLlMetag"'

# Multiple
bash ./launch.sh processing ${WORKFLOW_DIR}/main.nf ${WORKFLOW_DIR}/nextflow.config '--sample_type "standard" --input_type "multiple" --input_file multiple.csv  --errorStrategy "ignore" --technology "nanopore" --assay_suffix "_GLlMetag"'

# Post processing 
bash ./launch.sh post_processing  ${UPDATED_WORKFLOW_DIR}/post_processing.nf  ${UPDATED_WORKFLOW_DIR}/config/post_processing.config \
           '--name "First M. Last" --email "username@nasa.com" --glds_accession "OSD-944" --osd_accession "OSD-944" 
           --assay_table assay_table.txt --runsheet "multiple.csv" --technology "nanopore" --sample_type "standard" 
           --protocol_id "GL-DPPD-7116" --raw_suffix ".fastq.gz" 
           --raw_zip "../Merged_Sequence_Data/MultiQC_Reports/raw_multiqc_data.zip" 
           --filtered_zip "../Filtered_Sequence_Data/MultiQC_Reports/filtered_multiqc_data.zip" 
           --trimmed_zip "../Trimmed_Sequence_Data/MultiQC_Reports/trimmed_multiqc_data.zip" 
           --human_zip "../HR-removed_Sequence_Data/MultiQC_Reports/HRrm_multiqc_data.zip" ' 

# ------------------------------------------------------------------------------- Illumina ----------------------------------------------------------------------------------------- #

# ---- Low biomass
# Processing
bash ./launch.sh processing ${WORKFLOW_DIR}/main.nf ${WORKFLOW_DIR}/nextflow.config '--sample_type "low_biomass" --input_file sub_PE_file.csv --errorStrategy "ignore"  --technology "illumina" --assay_suffix "_GLlbsMetag"'

# Post processing 
bash ./launch.sh post_processing  ${WORKFLOW_DIR}/post_processing.nf  ${WORKFLOW_DIR}/config/post_processing.config \
           '--name "First M. Last" --email "username@nasa.com" --glds_accession "OSD-944" --osd_accession "OSD-944" 
           --assay_table assay_table.txt  --runsheet "PE_file.csv" --technology "illumina" --sample_type "low_biomass" 
           --protocol_id "GL-DPPD-7117" --raw_suffix ".fastq.gz" --raw_R1_suffix "_R1.fastq.gz" --raw_R2_suffix "_R2.fastq.gz"  
           --filtered_zip "../Filtered_Sequence_Data/MultiQC_Reports/filtered_multiqc_data.zip" 
           --raw_zip "../Merged_Sequence_Data/MultiQC_Reports/HRrm_multiqc_data.zip" 
           --human_summary "../../01-HRrm/results/human-read-removal-summary.tsv" 
           --decontam_zip "../Decontaminated_Sequence_Data/MultiQC_Reports/decontam_multiqc_data.zip" '


# ---- Standard
# Processing
bash ./launch.sh processing ${WORKFLOW_DIR}/main.nf ${WORKFLOW_DIR}/nextflow.config '--sample_type "standard" --input_file sub_PE_file.csv --errorStrategy "ignore"  --technology "illumina" --assay_suffix "_GLsMetag"'

# Post processing 
bash ./launch.sh post_processing  ${WORKFLOW_DIR}/post_processing.nf  ${WORKFLOW_DIR}/config/post_processing.config \
           '--name "First M. Last" --email "username@nasa.com" --glds_accession "OSD-944" --osd_accession "OSD-944" 
           --assay_table assay_table.txt --runsheet "PE_file.csv" --technology "illumina" --sample_type "standard" 
           --protocol_id "GL-DPPD-7117" --raw_suffix ".fastq.gz" --raw_R1_suffix "_R1.fastq.gz" --raw_R2_suffix "_R2.fastq.gz"  
           --filtered_zip "../Filtered_Sequence_Data/MultiQC_Reports/filtered_multiqc_data.zip" 
           --raw_zip "../Merged_Sequence_Data/MultiQC_Reports/HRrm_multiqc_data.zip" 
           --human_summary "../../01-HRrm/results/human-read-removal-summary.tsv" '

## Add a time-stamp at the end of the job then calculate how long the job took to run in seconds, minutes, and hours ##
echo ""
end=$(date +%s)
echo "end time: $end"
runtime_s=$(echo $(( end - start )))
echo "total run time(s): $runtime_s"
sec_per_min=60
sec_per_hr=3600
runtime_m=$(echo "scale=2; $runtime_s / $sec_per_min;" | bc)
echo "total run time(m): $runtime_m"
runtime_h=$(echo "scale=2; $runtime_s / $sec_per_hr;" | bc)
echo "total run time(h): $runtime_h"
echo ""


## Print the slurm job ID so you have it recorded and can view slurm job statistics if needed ##
echo "slurm job ID: ${SLURM_JOB_ID}"
