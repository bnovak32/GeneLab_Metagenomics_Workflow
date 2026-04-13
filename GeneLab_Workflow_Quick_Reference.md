# GeneLab Metagenomics Workflow - Quick Reference Guide

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

### Nanopore Standard Sample
```bash
nextflow run main.nf \
  -resume \
  -profile slurm,singularity \
  --technology nanopore \
  --sample_type standard \
  --input_file nanopore_samples.csv
```

### Low-Biomass Illumina
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

## Input File Templates

### Illumina Paired-End
```csv
sample_id,forward,reverse,paired,group
sample1,reads/s1_R1.fq.gz,reads/s1_R2.fq.gz,true,control
sample2,reads/s2_R1.fq.gz,reads/s2_R2.fq.gz,true,treatment
```

### Illumina Single-End
```csv
sample_id,forward,paired,group
sample1,reads/s1.fq.gz,false,control
sample2,reads/s2.fq.gz,false,treatment
```

### Low-Biomass Template
```csv
sample_id,forward,reverse,paired,group,NTC,concentration
sample1,reads/s1_R1.fq.gz,reads/s1_R2.fq.gz,true,treatment,false,15.2
blank1,reads/blank1_R1.fq.gz,reads/blank1_R2.fq.gz,true,control,true,0.0
```

### Nanopore Single FASTQ
```csv
sample_id,reads,group
sample1,reads/s1.fastq.gz,control
sample2,reads/s2.fastq.gz,treatment
```

### Nanopore Multiple FASTQs
```csv
sample_id,reads,group
sample1,reads/s1_pass1.fastq.gz,reads/s1_pass2.fastq.gz,control
sample2,reads/s2_pass1.fastq.gz,reads/s2_pass2.fastq.gz,treatment
```

---

## Essential Parameters Cheat Sheet

| Parameter | Values | Description |
|-----------|--------|-------------|
| `--technology` | `illumina`, `nanopore` | **REQUIRED** |
| `--sample_type` | `standard`, `low_biomass` | **REQUIRED** |
| `--input_file` | CSV path | **REQUIRED** (or use --accession) |
| `--accession` | OSD-### / GLDS-### | Alternative to --input_file |
| `--workflow` | `both`, `read-based`, `assembly-based` | Default: `both` |
| `-profile` | `slurm,singularity` | **REQUIRED** |
| `-resume` | - | Resume from last checkpoint |

---

## Profile Combinations

| Execution | Container/Env | Command |
|-----------|---------------|---------|
| Slurm + Singularity | Singularity | `-profile slurm,singularity` |
| Slurm + Conda | Conda | `-profile slurm,conda` |
| Local + Docker | Docker | `-profile docker` |
| Local + Conda | Conda | `-profile conda` |
| Local + Singularity | Singularity | `-profile singularity` |

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

## Key Output Locations

```
../
├── Filtered_Sequence_Data/          # Clean reads
├── Read-based_Processing/
│   ├── Kraken2-Outputs/
│   ├── Kaiju-Outputs/
│   ├── Metaphlan-Outputs/          # Illumina only
│   └── HUMAnN-Outputs/
├── Assembly-based_Processing/
│   ├── assemblies/
│   ├── annotations-and-taxonomy/
│   ├── bins/
│   └── MAGs/                        # High-quality only
├── Metadata/
│   └── software_versions.txt
└── Resource_Usage/
    └── execution_report_*.html
```

---

## Important Output Files

### Read-Based
```
Read-based_Processing/
├── Kraken2-Outputs/Filtered-species-abundance-table.tsv
├── Kaiju-Outputs/Filtered-species-abundance-table.tsv
├── Metaphlan-Outputs/Filtered-species-abundance-table.tsv
└── HUMAnN-Outputs/
    ├── Gene-Families/Filtered-gene-families-grouped-uniref90.tsv
    ├── Gene-Families-KO/Filtered-gene-families-grouped-KO.tsv
    └── Pathway-Abundance/Filtered-pathway-abundance.tsv
```

### Assembly-Based
```
Assembly-based_Processing/
├── assemblies/*_assembly.fasta
├── annotations-and-taxonomy/
│   ├── Gene-taxonomy-grouped/Filtered-grouped-gene-taxonomy.tsv
│   ├── Gene-KO-annotation/Filtered-grouped-gene-KO-annotation.tsv
│   └── Contig-taxonomy/Filtered-grouped-contig-taxonomy.tsv
├── bins/Bins-overview.tsv
└── MAGs/
    ├── MAGs-overview.tsv
    └── GTDB-tk-classification.tsv
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

### Case 2: Human Microbiome (Remove Host)
```bash
nextflow run main.nf -resume -profile slurm,singularity \
  --technology illumina --sample_type standard \
  --input_file samples.csv \
  --host_name human
```

### Case 3: Environmental Low-Biomass
```bash
nextflow run main.nf -resume -profile slurm,singularity \
  --technology illumina --sample_type low_biomass \
  --input_file low_biomass.csv
```

### Case 4: Nanopore Basecalling from Pod5
```bash
nextflow run main.nf -resume -profile slurm,singularity \
  --technology nanopore --sample_type standard \
  --input_type directory \
  --input_dir /path/to/pod5_files/ \
  --input_file barcodes.csv \
  --kit_name SQK-RPB114-24
```

### Case 5: Swift 1S Library Prep
```bash
nextflow run main.nf -resume -profile slurm,singularity \
  --technology illumina --sample_type standard \
  --input_file samples.csv \
  --swift_1S true
```

---

## Execution Environment

### SLURM Cluster
```bash
-profile slurm,singularity
```
**Best for**: HPC environments with job scheduler

### Local Workstation
```bash
-profile singularity
```
**Best for**: Single machine with sufficient resources

### Cloud/Container
```bash
-profile docker
```
**Best for**: AWS, Google Cloud, Docker-enabled systems

### Conda (Any Environment)
```bash
-profile conda
```
**Best for**: Systems without container support

---

## Monitoring & Logs

### Check Progress
```bash
# Terminal output
tail -f .nextflow.log

# Resource usage (open after completion)
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

## Directory Structure Before Running

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

| Process | Default Memory | Configurable Via |
|---------|----------------|------------------|
| MEGAHIT | 100 GB | `--max_mem` |
| GTDB-Tk | 100-200 GB | `--use_gtdbtk_scratch_location` |
| HUMAnN3 | 40-80 GB | Fixed (per process config) |
| CAT | 20-100 GB | `--block_size` |
| BBMap Pileup | 5 GB | `--pileup_mem` |
| MetaBAT2 | 32 GB | Fixed |
| CheckM | 16-32 GB | `--reduced_tree` |

---

## Time Estimates (Approximate)

| Dataset Size | Workflow Mode | Time (Standard) | Time (Low-Biomass) |
|--------------|---------------|-----------------|---------------------|
| 1 sample, 10M reads | Both | 4-8 hours | 6-10 hours |
| 5 samples, 10M reads | Both | 12-24 hours | 18-30 hours |
| 10 samples, 20M reads | Both | 24-48 hours | 36-60 hours |
| 1 sample, 10M reads | Read-based only | 2-4 hours | 3-6 hours |
| 1 sample, 10M reads | Assembly-based only | 3-6 hours | 4-8 hours |

*Times vary significantly based on system resources and database sizes*

---

## Key Differences: Standard vs Low-Biomass

| Feature | Standard | Low-Biomass |
|---------|----------|-------------|
| Input CSV | sample_id, forward, reverse, paired, group | + NTC, concentration columns |
| Contaminant Removal | No | Yes (statistical decontamination) |
| Blank/NTC Samples | Not used | Required for decontamination |
| Output Directories | Standard outputs | + Decontaminated outputs |
| Processing Time | Faster | Slower (additional steps) |
| Recommended For | Normal samples | Space, skin, air samples |

---

## Getting Help

1. **Check documentation**: This guide + main documentation
2. **View logs**: `work/<hash>/.command.log`
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

## End of Quick Reference

For detailed information, see: **GeneLab_Metagenomics_Workflow_Documentation.md**
