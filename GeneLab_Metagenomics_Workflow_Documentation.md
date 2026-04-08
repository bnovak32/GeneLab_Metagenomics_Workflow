# GeneLab Metagenomics Workflow Documentation

## Overview

The **GeneLab Metagenomics Sequencing Data Processing Workflow** is a comprehensive Nextflow DSL2 pipeline developed by NASA's GeneLab (part of the Open Science Data Repository - OSDR) for processing both Illumina short-read and Oxford Nanopore long-read metagenomics sequencing data. This workflow implements NASA's standardized metagenomics data processing pipelines for both standard and low-biomass samples.

**Repository**: https://github.com/olabiyi/GeneLab_Metagenomics_Workflow  
**Branch**: DEV  
**Version**: 1.0.0  
**Author**: Olabiyi Aderemi Obayomi  
**Nextflow Version**: >=24.04.4

---

## Workflow Capabilities

### Supported Technologies
- **Illumina**: Short-read sequencing (paired-end and single-end)
- **Nanopore**: Long-read sequencing (Oxford Nanopore Technologies)

### Sample Types
- **Standard**: Traditional metagenomics samples
- **Low Biomass**: Samples with limited microbial material requiring specialized processing (e.g., contaminant removal from negative controls)

### Analysis Modes
1. **Read-based analysis**: Taxonomic and functional profiling directly from reads
2. **Assembly-based analysis**: Genome assembly, binning, and annotation
3. **Both**: Combined read-based and assembly-based workflows (default)

---

## Workflow Architecture

### Main Entry Point
- **`main.nf`**: Primary workflow script that orchestrates the entire pipeline

### Technology-Specific Workflows
1. **`workflows/illumina.nf`**: Illumina short-read processing
2. **`workflows/nanopore.nf`**: Oxford Nanopore long-read processing
3. **`workflows/post_processing.nf`**: Post-processing and reporting

### Core Processing Modules (modules/)

#### Quality Control & Preprocessing
- **`quality_assessment.nf`**: FastQC, FASTP trimming/filtering, MultiQC reporting
- **`remove_contaminant.nf`**: Low-biomass contaminant removal from negative controls
- **`remove_host.nf`**: Host sequence removal using Kraken2
- **`demultiplexing.nf`**: Barcode demultiplexing for Nanopore data

#### Read-Based Analysis (`read_based_processing.nf`)
- **Taxonomic classification**:
  - Metaphlan4 (Illumina only)
  - Kraken2 (standard and PlusPFP databases)
  - Kaiju (multiple database options)
- **Functional profiling**:
  - HUMAnN3 pathway analysis
  - Gene family annotation (UniRef90)
  - KEGG Ortholog (KO) annotation
- **Visualization**: Krona plots, heatmaps, barplots

#### Assembly-Based Analysis (`assembly_based_processing.nf`)
- **Assembly**:
  - MEGAHIT (Illumina)
  - Flye + Medaka polishing (Nanopore)
- **Gene prediction**: Prodigal
- **Annotation**:
  - KOFAM Scan (KEGG function annotation)
  - CAT (Contig Annotation Tool) for taxonomy
- **Binning**: MetaBAT2
- **Quality assessment**: CheckM
- **MAG characterization**: GTDB-Tk taxonomy
- **Coverage analysis**: BBMap pileup.sh

#### Supporting Modules
- **`database_creation.nf`**: Automatic database download and setup
- **`create_runsheet.nf`**: GeneLab accession metadata retrieval
- **`downstream_analysis.nf`**: Statistical filtering, visualization, decontamination
- **`visualize_taxonomy.nf`**: Krona report generation
- **`genome_mapping.nf`**: Custom genome alignment
- **`genelab.nf`**: GeneLab-specific formatting and outputs

---

## Input Requirements

### Input File Formats

The workflow accepts CSV input files with different formats depending on sequencing type:

#### Illumina Paired-End (PE_file.csv)
```csv
sample_id,forward,reverse,paired,group
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,true,control
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,true,treatment
```

#### Illumina Single-End (SE_file.csv)
```csv
sample_id,forward,paired,group
sample1,/path/to/sample1.fastq.gz,false,control
sample2,/path/to/sample2.fastq.gz,false,treatment
```

#### Low-Biomass Samples (additional columns)
```csv
sample_id,forward,reverse,paired,group,NTC,concentration
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,true,treatment,false,10.5
blank1,/path/to/blank1_R1.fastq.gz,/path/to/blank1_R2.fastq.gz,true,control,true,0.0
```

#### Nanopore Inputs
- **`single.csv`**: One FASTQ file per sample
- **`multiple.csv`**: Multiple FASTQ files per sample (will be concatenated)
- **`input_dir_barcodes.csv`**: Directory containing pod5/fast5 files with barcode information

### GeneLab Accession Support
Instead of providing an input file, you can directly specify a GeneLab accession:
```bash
--accession OSD-574
```

---

## Key Parameters

### Required Parameters

| Parameter | Description | Options |
|-----------|-------------|---------|
| `--technology` | Sequencing platform | `illumina`, `nanopore` |
| `--sample_type` | Sample biomass level | `standard`, `low_biomass` |
| `-profile` | Execution environment | `slurm`, `singularity`, `docker`, `conda` (can combine: `slurm,singularity`) |

### Input Selection (Choose ONE)

| Parameter | Description |
|-----------|-------------|
| `--input_file` | Path to CSV input file |
| `--accession` | GeneLab/OSD accession number |

### Workflow Control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--workflow` | `both` | Analysis type: `read-based`, `assembly-based`, or `both` |
| `--errorStrategy` | `ignore` | Nextflow error handling strategy |
| `--publishDir_mode` | `link` | How outputs are published (`link`, `copy`, `symlink`) |

### Assembly & Binning Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_mem` | `100e9` | Maximum memory for MEGAHIT (100GB) |
| `--pileup_mem` | `5g` | Memory for BBMap pileup coverage calculation |
| `--block_size` | `4` | CAT/DIAMOND block size (lower = less RAM) |
| `--reduced_tree` | `True` | Use CheckM reduced tree (limits RAM to ~16GB) |

### MAG Quality Thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_est_comp` | `90` | Minimum estimated MAG completion (%) |
| `--max_est_redund` | `10` | Maximum estimated MAG redundancy (%) |
| `--max_est_strain_het` | `50` | Maximum estimated strain heterogeneity (%) |

### Host/Contaminant Removal

| Parameter | Description |
|-----------|-------------|
| `--host_name` | Host species name for database building (e.g., `human`) |
| `--host_url` | URL to pre-built Kraken2 host database |
| `--host_fasta` | Path to host genome FASTA file |
| `--host_db_dir` | Path to existing host Kraken2 database |
| `--swift_1S` | Enable Swift 1S library prep trimming (`true`/`false`) |

---

## Database Management

The workflow can **automatically download and set up** required databases if not provided. Database paths are specified via parameters:

### Read-Based Analysis Databases

| Database | Parameter | Auto-download |
|----------|-----------|---------------|
| Kraken2 PlusPFP | `--krakendb_dir` | ✅ Yes |
| Kaiju | `--kaijudb_dir` | ✅ Yes (specify `--kaijudb_name`) |
| Metaphlan4 | `--metaphlan_db_dir` | ✅ Yes |
| HUMAnN3 Chocophlan | `--chocophlan_dir` | ✅ Yes |
| HUMAnN3 UniRef | `--uniref_dir` | ✅ Yes |
| HUMAnN3 Utilities | `--utilities_dir` | ✅ Yes |

### Assembly-Based Analysis Databases

| Database | Parameter | Auto-download |
|----------|-----------|---------------|
| CAT | `--cat_db` | ✅ Yes (via `--CAT_DB_LINK`) |
| KOFAM Scan | `--ko_db_dir` | ✅ Yes |
| GTDB-Tk | `--gtdbtk_db_dir` | ✅ Yes (via `--GTDBTK_LINK`) |

### Database Root Directory
```bash
--DB_ROOT /full/path/to/Reference_DBs/
```
⚠️ **Important**: Use absolute paths, not relative paths (`~/` or `../` will fail)

---

## Usage Examples

### Example 1: Illumina Paired-End Standard Sample (Slurm + Singularity)
```bash
nextflow run main.nf \
  -resume \
  -profile slurm,singularity \
  --technology illumina \
  --sample_type standard \
  --input_file PE_samples.csv \
  --workflow both
```

### Example 2: Low-Biomass Illumina with Conda (Local Execution)
```bash
nextflow run main.nf \
  -resume \
  -profile conda \
  --technology illumina \
  --sample_type low_biomass \
  --input_file low_biomass_PE.csv \
  --workflow both
```

### Example 3: GeneLab Accession with Pre-existing Databases
```bash
nextflow run main.nf \
  -resume \
  -profile slurm,singularity \
  --technology illumina \
  --sample_type standard \
  --accession OSD-574 \
  --krakendb_dir /databases/kraken2_pluspfp_20251015/ \
  --metaphlan_db_dir /databases/metaphlan4-db/ \
  --cat_db /databases/CAT_prepare_20210107/
```

### Example 4: Nanopore Long-Reads with Host Removal
```bash
nextflow run main.nf \
  -resume \
  -profile docker \
  --technology nanopore \
  --sample_type standard \
  --input_file nanopore_single.csv \
  --host_name human \
  --workflow assembly-based
```

### Example 5: Assembly-Based Only with Custom Resources
```bash
nextflow run main.nf \
  -resume \
  -profile slurm,singularity \
  --technology illumina \
  --sample_type standard \
  --input_file samples.csv \
  --workflow assembly-based \
  --max_mem 200e9 \
  --pileup_mem 10g
```

---

## Output Structure

The workflow generates outputs in the following directory structure:

```
../
├── Raw_Sequence_Data/                 # Raw read FastQC reports
├── Filtered_Sequence_Data/            # Quality-filtered reads
│   ├── *_filtered.fastq.gz
│   └── FastQC/
├── Decontaminated_Sequence_Data/      # Low-biomass decontaminated reads
├── HostRM-removed_Sequence_Data/      # Host-removed reads
├── Read-based_Processing/
│   ├── Kraken2-Outputs/
│   │   ├── *_kraken_report.txt
│   │   ├── Filtered-species-abundance-table.tsv
│   │   └── krona.html
│   ├── Kaiju-Outputs/
│   │   ├── *_kaiju.out
│   │   ├── Filtered-species-abundance-table.tsv
│   │   └── krona.html
│   ├── Metaphlan-Outputs/             # Illumina only
│   │   ├── *_metaphlan_bugs_list.tsv
│   │   ├── Filtered-species-abundance-table.tsv
│   │   └── metaphlan-combined-barplot.pdf
│   └── HUMAnN-Outputs/
│       ├── Gene-Families/
│       │   ├── *_genefamilies.tsv
│       │   ├── Filtered-gene-families-grouped-uniref90.tsv
│       │   └── heatmap.pdf
│       ├── Gene-Families-KO/
│       │   ├── *_genefamilies_ko.tsv
│       │   └── Filtered-gene-families-grouped-KO.tsv
│       └── Pathway-Abundance/
│           ├── *_pathabundance.tsv
│           └── Filtered-pathway-abundance.tsv
├── Assembly-based_Processing/
│   ├── assemblies/
│   │   ├── *_assembly.fasta
│   │   └── Assembly-summaries.tsv
│   ├── predicted-genes/
│   │   ├── *_genes.faa
│   │   └── *_genes.fna
│   ├── annotations-and-taxonomy/
│   │   ├── *_gene_annotations.tsv
│   │   ├── *_contig_taxonomy.tsv
│   │   ├── Gene-taxonomy-grouped/
│   │   │   ├── Filtered-grouped-gene-taxonomy.tsv
│   │   │   └── heatmap.pdf
│   │   ├── Gene-KO-annotation/
│   │   │   ├── Filtered-grouped-gene-KO-annotation.tsv
│   │   │   └── heatmap.pdf
│   │   └── Contig-taxonomy/
│   │       ├── Filtered-grouped-contig-taxonomy.tsv
│   │       └── heatmap.pdf
│   ├── read-mapping/
│   │   ├── *_mapped_sorted.bam
│   │   └── *_coverage.tsv
│   ├── bins/
│   │   ├── sample1.*.fa
│   │   └── Bins-overview.tsv
│   ├── MAGs/                          # High-quality bins only
│   │   ├── sample1.*.fa
│   │   ├── MAGs-overview.tsv
│   │   └── GTDB-tk-classification.tsv
│   └── combined-outputs/
│       └── Assembly-based-processing-overview.tsv
├── Metadata/
│   ├── software_versions.txt
│   └── metadata_file.txt
└── Resource_Usage/
    ├── execution_timeline_*.html
    ├── execution_report_*.html
    └── execution_trace_*.txt
```

---

## Key Output Files

### Read-Based Analysis

#### Taxonomic Classification
- **Kraken2**: `Filtered-species-abundance-table.tsv`, `krona.html`
- **Kaiju**: `Filtered-species-abundance-table.tsv`, `krona.html`
- **Metaphlan**: `Filtered-species-abundance-table.tsv`, barplots

#### Functional Analysis
- **Gene Families (UniRef90)**: `Filtered-gene-families-grouped-uniref90.tsv`
- **Gene Families (KO)**: `Filtered-gene-families-grouped-KO.tsv`
- **Pathways**: `Filtered-pathway-abundance.tsv`

### Assembly-Based Analysis

#### Assemblies & Genes
- **Assemblies**: `*_assembly.fasta`
- **Predicted genes**: `*_genes.faa` (protein), `*_genes.fna` (nucleotide)

#### Annotation & Taxonomy
- **Gene annotations**: `*_gene_annotations.tsv` (KO functions + taxonomy + coverage)
- **Contig taxonomy**: `*_contig_taxonomy.tsv` (CAT results + coverage)
- **Combined tables**: Filtered, grouped abundance tables for genes and contigs

#### Binning & MAGs
- **All bins**: Individual FASTA files in `bins/`
- **MAGs**: High-quality bins (>90% complete, <10% redundant) in `MAGs/`
- **GTDB-Tk**: `GTDB-tk-classification.tsv` (phylogenetic placement)
- **CheckM**: Quality metrics in overview files

---

## Workflow Logic

### Main Workflow Execution Flow

```
main.nf
  │
  ├─→ Input validation (accession OR input_file)
  │
  ├─→ Technology-specific workflow
  │   ├─→ illumina.nf
  │   │   ├─→ Quality control (FastQC, FASTP)
  │   │   ├─→ [Low biomass] Contaminant removal
  │   │   ├─→ [Optional] Host removal
  │   │   └─→ Returns: clean_reads, metadata
  │   │
  │   └─→ nanopore.nf
  │       ├─→ [Optional] Demultiplexing (pod5/fast5)
  │       ├─→ Quality filtering
  │       ├─→ [Low biomass] Contaminant removal
  │       ├─→ [Optional] Host removal
  │       └─→ Returns: clean_reads, metadata
  │
  ├─→ Analysis mode selection
  │   ├─→ read-based → run_read_based_analysis()
  │   ├─→ assembly-based → run_assembly_based_analysis()
  │   └─→ both → run both workflows
  │
  └─→ Software version collection
```

### Illumina Workflow Details

```
illumina.nf
  │
  ├─→ Parse CSV input (PE/SE detection)
  │
  ├─→ QC: Raw FastQC
  │
  ├─→ Filtering: FASTP
  │   ├─→ Adapter trimming
  │   ├─→ Quality filtering
  │   └─→ [Optional] PolyG trimming (NextSeq/NovaSeq)
  │
  ├─→ QC: Filtered FastQC + MultiQC
  │
  ├─→ [Low biomass only] Contaminant removal
  │   ├─→ Statistical decontamination from blanks
  │   └─→ QC: Decontaminated FastQC
  │
  ├─→ [Optional] Host removal (Kraken2)
  │   ├─→ Database setup if needed
  │   ├─→ Classification & extraction of non-host reads
  │   └─→ QC: Host-removed FastQC
  │
  └─→ Output: clean_reads channel
```

### Nanopore Workflow Details

```
nanopore.nf
  │
  ├─→ Input type handling
  │   ├─→ [pod5/fast5] Dorado basecalling + demultiplexing
  │   ├─→ [multiple fastq] Concatenation
  │   └─→ [single fastq] Direct processing
  │
  ├─→ QC: Raw NanoPlot
  │
  ├─→ Filtering: Chopper
  │   ├─→ Length filtering
  │   └─→ Quality score filtering
  │
  ├─→ QC: Filtered NanoPlot
  │
  ├─→ [Low biomass] Contaminant removal
  │
  ├─→ [Optional] Host removal
  │
  └─→ Output: clean_reads channel
```

### Read-Based Processing Details

```
read_based_processing.nf
  │
  ├─→ Database setup (if needed)
  │   ├─→ Kraken2
  │   ├─→ Kaiju
  │   ├─→ Metaphlan
  │   └─→ HUMAnN3 (Chocophlan, UniRef, Utilities)
  │
  ├─→ [Illumina] HUMAnN3 + Metaphlan4
  │   ├─→ Per-sample profiling
  │   ├─→ Merge tables
  │   ├─→ Normalize (CPM, relab)
  │   ├─→ Group/stratify
  │   ├─→ Filtering (rare taxa/functions)
  │   └─→ Visualization
  │
  ├─→ Kraken2
  │   ├─→ Classify reads
  │   ├─→ Generate count tables
  │   ├─→ Filter rare species
  │   └─→ Krona plots + barplots
  │
  ├─→ Kaiju
  │   ├─→ Classify reads
  │   ├─→ Species tables
  │   ├─→ Filter rare species
  │   └─→ Krona plots + barplots
  │
  └─→ [Low biomass] Statistical decontamination
      └─→ Batch correction based on blanks
```

### Assembly-Based Processing Details

```
assembly_based_processing.nf
  │
  ├─→ Database setup
  │   ├─→ CAT
  │   ├─→ KOFAM Scan
  │   └─→ GTDB-Tk
  │
  ├─→ Assembly
  │   ├─→ [Illumina] MEGAHIT
  │   └─→ [Nanopore] Flye + Medaka polishing
  │
  ├─→ Header renaming (standardization)
  │
  ├─→ Gene prediction (Prodigal)
  │
  ├─→ Functional annotation (KOFAM Scan)
  │
  ├─→ Contig taxonomy (CAT)
  │
  ├─→ Read mapping
  │   ├─→ [Illumina] Bowtie2
  │   └─→ [Nanopore] Minimap2
  │
  ├─→ Coverage calculation (BBMap pileup)
  │
  ├─→ Combine annotations
  │   ├─→ Gene-level: KO + taxonomy + coverage
  │   └─→ Contig-level: taxonomy + coverage
  │
  ├─→ Binning (MetaBAT2)
  │
  ├─→ Bin quality (CheckM)
  │
  ├─→ MAG filtering
  │   └─→ Completion ≥90%, Redundancy ≤10%
  │
  ├─→ MAG taxonomy (GTDB-Tk)
  │
  ├─→ Summary tables
  │   ├─→ Gene taxonomy/function abundance
  │   ├─→ Contig taxonomy abundance
  │   ├─→ Bin overview
  │   └─→ MAG overview
  │
  └─→ [Low biomass] Decontamination
      └─→ Statistical correction
```

---

## Configuration Files

### `config/params.config`
Global parameters including:
- Technology and sample type
- Input/output directories
- Database paths and URLs
- Analysis thresholds
- Conda environment paths

### `config/default.config`
Default process settings:
- CPU/memory allocations
- Container images
- Error strategies
- Publishing directories

### `config/illumina.config`
Illumina-specific process configurations:
- MEGAHIT assembly settings
- Bowtie2 mapping parameters
- Illumina-specific tool containers

### `config/nanopore.config`
Nanopore-specific process configurations:
- Flye assembly parameters
- Medaka polishing settings
- Minimap2 mapping options
- Dorado basecalling models

### `config/profiles.config`
Execution profiles:
- `standard`: Local execution
- `slurm`: SLURM job scheduler
- `conda`: Conda environment management
- `singularity`: Singularity containers
- `docker`: Docker containers

---

## Container & Environment Management

### Singularity Containers (Recommended)
The workflow automatically pulls Singularity images from:
- Docker Hub
- Quay.io
- BioContainers

### Docker Support
All processes support Docker execution via `-profile docker`

### Conda Environments
Pre-defined environments in `envs/`:
- `humann3.yml`
- `cat.yml`
- `metabat.yml`
- `gtdbtk.yml`
- `megahit.yml`
- And more...

You can provide paths to existing conda environments:
```bash
--conda_megahit /path/to/existing/megahit/env
```

---

## Resource Requirements

### Minimum Recommendations
- **CPU**: 10 cores per process
- **Memory**: 32 GB minimum, 300 GB maximum per process
- **Storage**: Variable (depends on dataset size)
  - Raw data: ~original dataset size
  - Databases: ~200-500 GB (all databases combined)
  - Outputs: ~2-5x raw data size

### High-Memory Processes
- **MEGAHIT**: Up to 100 GB (configurable via `--max_mem`)
- **GTDB-Tk**: 100-200 GB (use `--use_gtdbtk_scratch_location` to offload to disk)
- **HUMAnN3**: 40-80 GB per sample
- **CAT**: Depends on `--block_size` (lower = less RAM)

---

## Advanced Features

### GeneLab Integration
- Direct data retrieval from GeneLab/OSDR via accession numbers
- Automatic runsheet generation
- GeneLab-specific file naming conventions (`--assay_suffix`, `--additional_filename_prefix`)

### Low-Biomass Sample Processing
Specialized statistical decontamination:
1. Identifies negative controls (NTC = true)
2. Performs batch correction
3. Removes contaminant signals from samples
4. Applies to both read-based and assembly-based results

### Custom Genome Mapping
Map reads to custom reference genomes:
```bash
--custome_genome /path/to/reference.fna
```

### Monitoring with Seqera Platform
Enable tower integration in `config/profiles.config`:
```groovy
tower {
    accessToken = 'your-token'
    enabled = true
}
```

---

## Troubleshooting

### Common Issues

#### 1. Database Download Failures
**Solution**: Manually download databases and provide paths:
```bash
--cat_db /full/path/to/CAT_prepare_20210107/
--gtdbtk_db_dir /full/path/to/GTDB-tk-ref-db/
```

#### 2. Out of Memory Errors
**Solutions**:
- Reduce `--max_mem` for MEGAHIT
- Use `--reduced_tree True` for CheckM
- Enable GTDB-Tk scratch: `--use_gtdbtk_scratch_location true`
- Lower `--block_size` for CAT
- Reduce `--pileup_mem`

#### 3. GTDB-Tk Memory Issues
**Solution**: Use scratch directory to offload RAM:
```bash
--use_gtdbtk_scratch_location true
```

#### 4. Process Failures
**Solutions**:
- Check logs in `../Logs/` and `work/` directories
- Review resource usage in `../Resource_Usage/`
- Adjust `--errorStrategy` (default: `ignore`)
- Enable debug mode: `--debug true`

---

## Software Dependencies

The workflow integrates the following major tools:

### Quality Control
- FastQC
- MultiQC
- NanoPlot (Nanopore)
- FASTP (Illumina)
- Chopper (Nanopore)

### Taxonomic Classification
- Metaphlan4 (Illumina)
- Kraken2
- Kaiju
- CAT (assembly)
- GTDB-Tk (MAGs)

### Functional Annotation
- HUMAnN3
- KOFAM Scan
- KEGG Decoder

### Assembly & Binning
- MEGAHIT (Illumina)
- Flye (Nanopore)
- Medaka (Nanopore polishing)
- Prodigal (gene calling)
- MetaBAT2 (binning)
- CheckM (bin quality)

### Mapping & Coverage
- Bowtie2 (Illumina)
- Minimap2 (Nanopore)
- Samtools
- BBMap

### Visualization
- Krona
- R (ggplot2, heatmaps, barplots)

---

## Citation & References

If you use this workflow, please cite:
- **GeneLab Metagenomics Workflow**: https://github.com/olabiyi/GeneLab_Metagenomics_Workflow
- **NASA OSDR**: https://osdr.nasa.gov/
- Individual tool citations (see `software_versions.txt` in outputs)

### Pipeline Documents
- GL-DPPD-7107-A: Standard Illumina metagenomics
- GL-DPPD-7116: Nanopore metagenomics
- GL-DPPD-7117: Low-biomass metagenomics

---

## Support & Contact

- **Issues**: https://github.com/olabiyi/GeneLab_Metagenomics_Workflow/issues
- **GeneLab**: https://genelab.nasa.gov/
- **OSDR**: https://osdr.nasa.gov/

---

## Changelog

### Version 1.0.0 (Current - DEV branch)
- Initial Nextflow DSL2 implementation
- Support for Illumina and Nanopore
- Standard and low-biomass sample types
- Read-based and assembly-based workflows
- Automatic database management
- GeneLab accession integration
- Comprehensive quality control and reporting

---

## License

This workflow is developed by NASA GeneLab and is part of the Open Science Data Repository (OSDR) initiative. Please refer to the repository for license information.
