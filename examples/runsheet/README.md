# Runsheet File Specification

## Description

* The runsheet is a comma-separated file that contains the metadata required for processing metagenomics sequence datasets through the GeneLab Metagenomics data processing pipelines.

* The runsheet format is dependent on the workflow and approach used for processing. 

<br> 

## Jump To:

- [Low Biomass Long Read Runsheet Info](#low-biomass-long-read-examples)
- [Low Biomass Short Read Runsheet Info](#low-biomass-short-read-examples)

<br> 

## Low Biomass Long Read Examples

1. Example runsheet for Approach 1: Start with pod5 or fast5 files as input | [input_dir_barcodes.csv](input_dir_barcodes.csv)
2. Example runsheet for Approach 2: Start with multiple FASTQ files per sample as input | [multiple.csv](multiple.csv)
3. Example runsheet for Approach 3: Start with one FASTQ file per sample as input | [single.csv](single.csv)

<br> 

## Low Biomass Long Read Required Columns

### Approach 1: Start with pod5 or fast5 files as input
> *Note: There should be one row per sample. See the [input_dir_barcodes.csv](input_dir_barcodes.csv) as an example.*

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| sample_id | string | Unique sample name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | sample-1 |
| barcode_id | string | Unique barcode ID associated with the respective sample. Should not include spaces or weird characters. | barcode01 |
| group | string | Name of the treatment group that the sample belongs to. Should not include spaces or weird characters. | cleanroom_floor |
| NTC | bool | Set to TRUE if the sample is a non-template control (NTC), and FALSE if it is not. | FALSE |
| concentration (ng) | float | The amount of input DNA used for library preparation. | 0.005 |
| paired | bool | Set to TRUE if the samples were sequenced as paired-end. If set to FALSE, samples are assumed to be single-end. | FALSE |

<br> 

### Approach 2: Start with multiple FASTQ files per sample as input
> *Note: There should be one row per fastq file. So there will be multiple rows for each sample depending on the number of fastq files associated with that sample. See the [multiple.csv](multiple.csv) as an example.*

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| sample_id | string | Unique sample name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | sample-1 |
| barcode_id | string | Unique barcode ID associated with the respective sample. Should not include spaces or weird characters. | barcode01 |
| forward | string (local path or URL) | Location of one of the raw FASTQ files associated with the respective sample. *Note: If more than one FASTQ file is associated with the sample, add an additional row for each FASTQ file.* | /fastq_pass/barcode01/PBC13394_pass_barcode01_6db14c65_7b9d01fc_51.fastq.gz |
| group | string | Name of the treatment group that the sample belongs to. Should not include spaces or weird characters. | cleanroom_floor |
| NTC | bool | Set to TRUE if the sample is a non-template control (NTC), and FALSE if it is not. | FALSE |
| concentration (ng) | float | The amount of input DNA used for library preparation. | 0.005 |
| paired | bool | Set to TRUE if the samples were sequenced as paired-end. If set to FALSE, samples are assumed to be single-end. | FALSE |

<br> 

### Approach 3: Start with one FASTQ file per sample as input
> *Note: This approach assumes that all FASTQ files associated with each sample have been concatenated so that there is only one FASTQ files per sample. Therefore, there should be one row per sample. See the [single.csv](single.csv) as an example.*

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| sample_id | string | Unique sample name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | sample-1 |
| barcode_id | string | Unique barcode ID associated with the respective sample. Should not include spaces or weird characters. | barcode01 |
| forward | string (local path or URL) | Location of the raw FASTQ file associated with the respective sample. | /raw_fastqs/sample1.fastq.gz |
| group | string | Name of the treatment group that the sample belongs to. Should not include spaces or weird characters. | cleanroom_floor |
| NTC | bool | Set to TRUE if the sample is a non-template control (NTC), and FALSE if it is not. | FALSE |
| concentration (ng) | float | The amount of input DNA used for library preparation. | 0.005 |
| paired | bool | Set to TRUE if the samples were sequenced as paired-end. If set to FALSE, samples are assumed to be single-end. | FALSE |

<br>

## Low Biomass Short Read Examples

1. Example runsheet for Approach 1: Start with paired-end FASTQ files as input | [PE_file.csv](PE_file.csv)
2. Example runsheet for Approach 2: Start with single-end FASTQ files as input | [SE_file.csv](SE_file.csv)

<br> 

## Low Biomass Short Read Required Columns

### Approach 1: Start with paired-end FASTQ files as input
> *Note: There should be one row per sample. See the [PE_file.csv](PE_file.csv) as an example.*

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| sample_id | string | Unique sample name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | sample-1 |
| forward | string (local path or URL) | Location of the raw FASTQ file associated with the respective sample. For paired-end data, this specifies the forward reads fastq.gz file. | /raw_fastqs/sample1_R1_raw.fastq.gz |
| reverse | string (local path or URL) | Location of the raw FASTQ file associated with the respective sample. For paired-end data, this specifies the reverse reads fastq.gz file. | /raw_fastqs/sample1_R2_raw.fastq.gz |
| group | string | Name of the treatment group that the sample belongs to. Should not include spaces or weird characters. | cleanroom_floor |
| NTC | bool | Set to TRUE if the sample is a non-template control (NTC), and FALSE if it is not. | FALSE |
| concentration (ng) | float | The amount of input DNA used for library preparation. | 0.005 |
| paired | bool | Set to TRUE if the samples were sequenced as paired-end. If set to FALSE, samples are assumed to be single-end. | TRUE |

<br> 

### Approach 2: Start with single-end FASTQ files as input
> *Note: There should be one row per sample. See the [SE_file.csv](SE_file.csv) as an example.*

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| sample_id | string | Unique sample name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | sample-1 |
| forward | string (local path or URL) | Location of the raw FASTQ file associated with the respective sample. | /raw_fastqs/sample1_raw.fastq.gz |
| group | string | Name of the treatment group that the sample belongs to. Should not include spaces or weird characters. | cleanroom_floor |
| NTC | bool | Set to TRUE if the sample is a non-template control (NTC), and FALSE if it is not. | FALSE |
| concentration (ng) | float | The amount of input DNA used for library preparation. | 0.005 |
| paired | bool | Set to TRUE if the samples were sequenced as paired-end. If set to FALSE, samples are assumed to be single-end. | TRUE |

