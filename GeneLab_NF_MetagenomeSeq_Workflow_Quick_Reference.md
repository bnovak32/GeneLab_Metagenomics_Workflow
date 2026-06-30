# GeneLab Metagenomics Workflow - Quick Reference Guide <!-- omit in toc -->

This is a quick start guide and basic reference for key parameter selection and troubleshooting. 
For more detailed information, see: [**GeneLab_Metagenomics_Workflow_Documentation.md**](GeneLab_Metagenomics_Workflow_Documentation.md)

- [Quick Start Commands](#quick-start-commands)
- [Input Sample File Templates](#input-sample-file-templates)
- [Essential Parameters Cheat Sheet](#essential-parameters-cheat-sheet)
- [Common Nextflow Profile Combinations](#common-nextflow-profile-combinations)
- [Database Parameters Quick Reference](#database-parameters-quick-reference)
- [Host Removal Options (Choose ONE)](#host-removal-options-choose-one)
- [Performance Tuning](#performance-tuning)
- [Key Output File Locations](#key-output-file-locations)
- [Troubleshooting Quick Fixes](#troubleshooting-quick-fixes)
- [MAG Quality Thresholds](#mag-quality-thresholds)
- [Workflow Modes](#workflow-modes)
- [Common Use Cases](#common-use-cases)
- [Monitoring \& Logs](#monitoring--logs)
- [Help Command](#help-command)
- [Version Information](#version-information)
- [Example Directory Structure Before Running](#example-directory-structure-before-running)
- [Memory Requirements Summary](#memory-requirements-summary)
- [Time Estimates (Approximate)](#time-estimates-approximate)
- [Key Differences: Standard vs Low Biomass](#key-differences-standard-vs-low-biomass)
- [Getting Help](#getting-help)
- [Best Practices](#best-practices)
- [Quick Validation Checklist](#quick-validation-checklist)


## Quick Start Commands

### Illumina Standard Sample
```bash
nextflow run main.nf \
  -resume \
  -profile slurm,singularity \
  --technology illumina \
  --sample_type standard \
  --input_file PE_samples.csv
```

### Low Biomass Nanopore
```bash
nextflow run main.nf \
  -resume \
  -profile slurm,singularity \
  --technology nanopore \
  --sample_type low_biomass \
  --input_file nanopore_samples.csv
```

### Low Biomass Illumina
```bash
nextflow run main.nf \
  -resume \
  -profile slurm,singularity \
  --technology illumina \
  --sample_type low_biomass \
  --input_file low_biomass_samples.csv
```

### Using GeneLab Accession
```bash
nextflow run main.nf \
  -resume \
  -profile slurm,singularity \
  --technology illumina \
  --sample_type standard \
  --accession OSD-574
```

---

## Input Sample File Templates

### Illumina Standard Paired-End
```csv
sample_id,forward,reverse,paired,group
sample1,reads/s1_R1.fq.gz,reads/s1_R2.fq.gz,TRUE,control
sample2,reads/s2_R1.fq.gz,reads/s2_R2.fq.gz,TRUE,treatment
```

### Illumina Standard Single-End
```csv
sample_id,forward,paired,group
sample1,reads/s1.fq.gz,FALSE,control
sample2,reads/s2.fq.gz,FALSE,treatment
```

### Illumina Low Biomass Paired-End
```csv
sample_id,forward,reverse,paired,group,NTC,concentration
sample1,reads/s1_R1.fq.gz,reads/s1_R2.fq.gz,TRUE,treatment,FALSE,15.2
blank1,reads/blank1_R1.fq.gz,reads/blank1_R2.fq.gz,TRUE,control,TRUE,0.0
```

### Nanopore Low Biomass Single FASTQ
```csv
sample_id,reads,group,NTC,concentration,paired
sample1,reads/s1.fastq.gz,control,TRUE,0.0,FALSE
sample2,reads/s2.fastq.gz,treatment,FALSE,15.2,FALSE
```

### Nanopore Low Biomass Multiple FASTQs
```csv
sample_id,barcode_id,forward,group,NTC,concentration,paired
sample1,barcode01,reads/barcode01/barcode01_0.fastq.gz,treatment,FALSE,0.005,FALSE
sample1,barcode01,reads/barcode01/barcode01_1.fastq.gz,treatment,FALSE,0.005,FALSE
sample2,barcode02,reads/barcode02/barcode02_0.fastq.gz,treatment,FALSE,0.005,FALSE
sample2,barcode02,reads/barcode02/barcode02_1.fastq.gz,treatment,FALSE,0.005,FALSE
ntc1,barcode09,reads/barcode09/barcode09_0.fastq.gz,control,TRUE,0,FALSE
ntc1,barcode09,reads/barcode09/barcode09_1.fastq.gz,control,TRUE,0,FALSE
```

### Nanopore Low Biomass POD5
```csv
sample_id,barcode_id,group,NTC,concentration,paired
sample1,barcode01,treatment,FALSE,0.1,FALSE
sample2,barcode02,control,TRUE,0,FALSE
```

---

## Essential Parameters Cheat Sheet

| Parameter       | Values                                 | Description                       |
| --------------- | -------------------------------------- | --------------------------------- |
| `--technology`  | `illumina`, `nanopore`                 | **REQUIRED**                      |
| `--sample_type` | `standard`, `low_biomass`              | **REQUIRED**                      |
| `--input_file`  | CSV path                               | **REQUIRED** (or use --accession) |
| `--accession`   | OSD-### / GLDS-###                     | Alternative to --input_file       |
| `--workflow`    | `both`, `read-based`, `assembly-based` | Default: `both`                   |
| `-profile`      | `slurm,singularity`                    | **REQUIRED**                      |
| `-resume`       | -                                      | Resume from last checkpoint       |

---

## Common Nextflow Profile Combinations

| Execution           | Container/Env | Command                      | Best for                                                                   |
| ------------------- | ------------- | ---------------------------- | -------------------------------------------------------------------------- |
| Slurm + Singularity | Singularity   | `-profile slurm,singularity` | HPC environments with a job scheduler                                      |
| Slurm + Conda       | Conda         | `-profile slurm,conda`       | HPC environments with a job scheduler                                      |
| Local + Docker      | Docker        | `-profile docker`            | Single machine with sufficient resources and docker container support      |
| Local + Conda       | Conda         | `-profile conda`             | Single machine with sufficient resources (without container support)       |
| Local + Singularity | Singularity   | `-profile singularity`       | Single machine with sufficient resources and singularity container support |

---

## Database Parameters Quick Reference

### Auto-Download (Leave as null)
```bash
# Workflow will download automatically:
--krakendb_dir null
--kaijudb_dir null
--metaphlan_db_dir null
--cat_db null
--gtdbtk_db_dir null
--ko_db_dir null
```

### Manual Database Paths
```bash
--krakendb_dir /databases/kraken2_pluspfp/
--kaijudb_dir /databases/kaiju_nr_euk/
--metaphlan_db_dir /databases/metaphlan4-db/
--chocophlan_dir /databases/humann3-db/chocophlan/
--uniref_dir /databases/humann3-db/uniref/
--utilities_dir /databases/humann3-db/utility_mapping/
--cat_db /databases/CAT_prepare_20210107/
--gtdbtk_db_dir /databases/GTDB-tk-ref-db/
--ko_db_dir /databases/kofamscan_db/
```

---

## Host Removal Options (Choose ONE)

### Option 1: Use Host Name
Supports all reference libraries available to the kraken2-build command (see full 
list in the [Kraken2 documentation](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases))
```bash
--host_name human
```

### Option 2: Download Pre-built Database
```bash
--host_url https://zenodo.org/records/8339700/files/k2_Human_20230629.tar.gz
```

### Option 3: Build from FASTA
```bash
--host_fasta /path/to/host_genome.fna
```

### Option 4: Existing Database
```bash
--host_db_dir /databases/kraken2-host-db/
```

---

## Performance Tuning

### High-Memory Datasets
```bash
--max_mem 200e9              # MEGAHIT: 200GB
--pileup_mem 10g             # BBMap: 10GB
--block_size 2               # CAT: Lower = less RAM
--use_gtdbtk_scratch_location true
```

### Low-Memory Systems
```bash
--max_mem 50e9
--pileup_mem 2g
--block_size 6
--reduced_tree True
```

---

## Key Output File Locations
For a complete list of output files, see [**GeneLab_Metagenomics_Workflow_Documentation.md**](GeneLab_Metagenomics_Workflow_Documentation.md)

```
../
├── Filtered_Sequence_Data/          # Clean reads
├── Read-based_Processing/
│   ├── Kraken2-Outputs/
│   ├── Kaiju-Outputs/
│   ├── Metaphlan-Outputs/           # Illumina only
│   └── HUMAnN-Outputs/
├── Assembly-based_Processing/
│   ├── assemblies/
│   ├── annotations-and-taxonomy/
│   ├── bins/
│   ├── combined-outputs/
│   ├── predicted-genes/
│   ├── read-mapping/
│   └── MAGs/                        # High-quality only
├── Metadata/
│   └── software_versions.txt
└── Resource_Usage/
    └── execution_report_*.html
```

---

## Troubleshooting Quick Fixes

### Problem: Out of Memory
```bash
# Reduce memory for key processes
--max_mem 50e9 --pileup_mem 2g --block_size 6 --reduced_tree True
```

### Problem: GTDB-Tk RAM Issues
```bash
--use_gtdbtk_scratch_location true
```

### Problem: Database Download Fails
```bash
# Download manually and provide path
--cat_db /path/to/CAT_prepare_20210107/
```

### Problem: Process Keeps Failing
```bash
# Check logs
less work/<process_hash>/.command.log

# Check resource usage
firefox ../Resource_Usage/execution_report_*.html
```

### Problem: Pipeline Stops Unexpectedly
```bash
# Resume from checkpoint
nextflow run main.nf -resume [other parameters]
```

---

## MAG Quality Thresholds

### Default Settings
```bash
--min_est_comp 90         # ≥90% complete
--max_est_redund 10       # ≤10% redundant
--max_est_strain_het 50   # ≤50% strain heterogeneity
```

### More Stringent (High-Quality Only)
```bash
--min_est_comp 95
--max_est_redund 5
--max_est_strain_het 25
```

### More Permissive (Medium-Quality MAGs)
```bash
--min_est_comp 50
--max_est_redund 20
--max_est_strain_het 100
```

---

## Workflow Modes

### Read-Based Only (Faster)
```bash
--workflow read-based
```
**Use when**: Only need taxonomic/functional profiles, no binning needed

### Assembly-Based Only
```bash
--workflow assembly-based
```
**Use when**: Only need MAGs, assemblies, or contig-level analysis

### Both (Default)
```bash
--workflow both
```
**Use when**: Comprehensive analysis needed

---

## Common Use Cases

### Case 1: NASA GeneLab Dataset
```bash
nextflow run main.nf -resume -profile slurm,singularity \
  --technology illumina --sample_type standard \
  --accession OSD-574
```

### Case 2: Mouse Microbiome (Remove Host)
```bash
nextflow run main.nf -resume -profile slurm,singularity \
  --technology illumina --sample_type standard \
  --input_file samples.csv \
  --host_name mouse
```

### Case 3: Environmental Low Biomass
```bash
nextflow run main.nf -resume -profile slurm,singularity \
  --technology illumina --sample_type low_biomass \
  --input_file low_biomass.csv
```

### Case 4: Environmental Low Biomass Nanopore from Multiple FASTQs
```bash
nextflow run main.nf -resume -profile slurm,singularity \
  --technology nanopore --sample_type low_biomass \
  --input_file multiple.csv
```
---

## Monitoring & Logs

### Check Progress
```bash
# Terminal output
tail -f .nextflow.log

# Resource usage (Open after completion. Using firefox, for example.)
firefox ../Resource_Usage/execution_timeline_*.html
firefox ../Resource_Usage/execution_report_*.html
```

### Find Failed Process
```bash
# Check trace file
less ../Resource_Usage/execution_trace_*.txt

# Navigate to work directory
cd work/<hash>/
less .command.log
less .command.err
```

---

## Help Command

```bash
nextflow run main.nf --help
```

Displays all parameters with descriptions.

---

## Version Information

```bash
# Workflow version
grep version nextflow.config

# Check Nextflow version
nextflow -version

# Software versions used in run
cat ../Metadata/software_versions.txt
```

---

## Example Directory Structure Before Running

```
project/
├── GeneLab_Metagenomics_Workflow/  # Cloned repo
│   ├── main.nf
│   ├── nextflow.config
│   ├── modules/
│   ├── workflows/
│   └── config/
├── reads/                          # Your FASTQ files
│   ├── sample1_R1.fastq.gz
│   └── sample1_R2.fastq.gz
├── samples.csv                     # Your input CSV
└── databases/                      # Optional: pre-downloaded DBs
    ├── kraken2_pluspfp/
    ├── metaphlan4-db/
    └── CAT_prepare_20210107/
```

### Run from Workflow Directory
```bash
cd GeneLab_Metagenomics_Workflow/
nextflow run main.nf [parameters]
```

---

## Memory Requirements Summary

| Process      | Default Memory | Configurable Via                |
| ------------ | -------------- | ------------------------------- |
| MEGAHIT      | 100 GB         | `--max_mem`                     |
| GTDB-Tk      | 100-200 GB     | `--use_gtdbtk_scratch_location` |
| HUMAnN3      | 40-80 GB       | Fixed (per process config)      |
| CAT          | 20-100 GB      | `--block_size`                  |
| BBMap Pileup | 5 GB           | `--pileup_mem`                  |
| MetaBAT2     | 32 GB          | Fixed                           |
| CheckM       | 16-32 GB       | `--reduced_tree`                |

---

## Time Estimates (Approximate)

| Dataset Size          | Workflow Mode       | Time (Standard) | Time (Low Biomass) |
| --------------------- | ------------------- | --------------- | ------------------ |
| 1 sample, 10M reads   | Both                | 4-8 hours       | 6-10 hours         |
| 5 samples, 10M reads  | Both                | 12-24 hours     | 18-30 hours        |
| 10 samples, 20M reads | Both                | 24-48 hours     | 36-60 hours        |
| 1 sample, 10M reads   | Read-based only     | 2-4 hours       | 3-6 hours          |
| 1 sample, 10M reads   | Assembly-based only | 3-6 hours       | 4-8 hours          |

*Times vary significantly based on system resources and database sizes*

Database download/build can also take a substantial amount of time. Kaiju database download, for example, can take up to two weeks. Timing will depend on database size and network speed.

---

## Key Differences: Standard vs Low Biomass

| Feature             | Standard                                   | Low Biomass                         |
| ------------------- | ------------------------------------------ | ----------------------------------- |
| Input CSV           | sample_id, forward, reverse, paired, group | + NTC, concentration columns        |
| Contaminant Removal | No                                         | Yes (statistical decontamination)   |
| Blank/NTC Samples   | Not used                                   | Required for decontamination        |
| Output Directories  | Standard outputs                           | + Decontaminated outputs            |
| Processing Time     | Faster                                     | Slower (additional steps)           |
| Recommended For     | Normal samples                             | Cleanroom, space, skin, air samples |

---

## Getting Help

1. **Check documentation**: This guide + main documentation
2. **View nextflow logs**: `nextflow log`
3. **Check resources**: `../Resource_Usage/execution_report_*.html`
4. **Enable debug**: `--debug true`
5. **GitHub issues**: https://github.com/olabiyi/GeneLab_Metagenomics_Workflow/issues

---

## Best Practices

✅ **DO**:
- Use `-resume` to restart from failures
- Provide absolute paths for databases
- Check input CSV format carefully
- Monitor resource usage reports
- Keep databases in persistent storage

❌ **DON'T**:
- Use relative paths (`~/`, `../`) for `--DB_ROOT`
- Mix PE and SE samples in same CSV
- Forget to specify `--technology` and `--sample_type`
- Run without sufficient disk space for databases
- Delete `work/` directory until analysis is complete

---

## Quick Validation Checklist

Before running, verify:
- [ ] Input CSV format matches technology/sample type
- [ ] All FASTQ files in CSV exist and are readable
- [ ] Sufficient disk space (databases + outputs)
- [ ] Correct `--technology` parameter (`illumina` or `nanopore`)
- [ ] Correct `--sample_type` parameter (`standard` or `low_biomass`)
- [ ] Profile specified (`-profile` parameter)
- [ ] Input method chosen (`--input_file` OR `--accession`)
  
---