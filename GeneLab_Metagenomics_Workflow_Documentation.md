# GeneLab Metagenomics Workflow Documentation

## Overview

The **GeneLab Metagenomics Sequencing Data Processing Workflow** is a comprehensive Nextflow DSL2 pipeline developed by NASA's GeneLab (part of the Open Science Data Repository - OSDR) for processing both Illumina short-read and Oxford Nanopore long-read metagenomics sequencing data. This workflow implements NASA's standardized metagenomics data processing pipelines for both standard and low-biomass samples.

**Repository**: https://github.com/nasa/GeneLab_Metagenomics_Workflow  
**Branch**: main  
**Version**: 1.0.0  
**Author**: Olabiyi Aderemi Obayomi  
**Nextflow Version**: 24.10.5-0

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
- **`quality_assessment.nf`**: FastQC, FASTP trimming/filtering, Filtlong filtering, PoreChop trimming,  MultiQC reporting etc
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
- **Visualization**: Heatmaps

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
|-- Assembly-based_Processing
|   |-- Assembly-based-processing-overview.tsv
|   |-- MAGs # High-quality bins only
|   |   |-- *-MAGs.zip
|   |   |-- MAG-KEGG-Decoder-out.html
|   |   |-- MAG-KEGG-Decoder-out.tmp
|   |   |-- MAG-KEGG-Decoder-out.tsv
|   |   |-- MAG-level-KO-annotations.tsv
|   |   |-- MAGs-overview.tsv
|   |-- annotations-and-taxonomy
|   |   |-- *-annotations.tsv
|   |   |-- *-contig-coverage-and-tax.tsv
|   |   |-- *-contig-tax.tsv
|   |   |-- *-gene-coverage-annotation-and-tax.tsv
|   |   |-- *-gene-tax.tsv
|   |-- assemblies
|   |   |-- *-assembly.fasta
|   |   |-- Failed-assemblies.tsv # Only if some sample assemblies fail
|   |   |-- assembly-summaries.tsv
|   |-- bins
|   |   |-- *-bins.zip
|   |   |-- bin-assembly-summaries.tsv
|   |   |-- bins-checkm-out.tsv
|   |   |-- bins-overview.tsv
|   |-- combined-outputs
|   |   |-- Contig-level
|   |   |   |-- Combined-contig-level-taxonomy-coverages-CPM.tsv
|   |   |   |-- Combined-contig-level-taxonomy-coverages.tsv
|   |   |   |-- Combined-contig-level-taxonomy.tsv
|   |   |   |-- Combined-contig-level-taxonomy_decontam_failure.txt # Only if decontam fails
|   |   |   |-- Combined-contig-level-taxonomy_decontam_results.tsv
|   |   |   |-- Combined-contig-level-taxonomy_filtered.tsv
|   |   |   |-- Combined-contig-level-taxonomy_filtered_heatmap.png
|   |   |   |-- Combined-contig-level-taxonomy_filtered_top_50_heatmap.png
|   |   |   |-- Combined-contig-level-taxonomy_unfiltered_heatmap.png
|   |   |   |-- Combined-contig-level-taxonomy_unfiltered_top_50_heatmap.png
|   |   |-- Gene-level
|   |       |-- KO
|   |       |   |-- Combined-gene-level-KO-function-coverages-CPM.tsv
|   |       |   |-- Combined-gene-level-KO-function-coverages.tsv
|   |       |   |-- Combined-gene-level-KO-function.tsv
|   |       |   |-- Combined-gene-level-KO-function_decontam_failure.txt # Only if decontam fails
|   |       |   |-- Combined-gene-level-KO-function_decontam_results.tsv
|   |       |   |-- Combined-gene-level-KO-function_filtered.tsv
|   |       |   |-- Combined-gene-level-KO-function_filtered_heatmap.png
|   |       |   |-- Combined-gene-level-KO-function_filtered_top_50_heatmap.png
|   |       |   |-- Combined-gene-level-KO-function_unfiltered_heatmap.png
|   |       |   |-- Combined-gene-level-KO-function_unfiltered_top_50_heatmap.png
|   |       |-- Taxonomy
|   |           |-- Combined-gene-level-taxonomy-coverages-CPM.tsv
|   |           |-- Combined-gene-level-taxonomy-coverages.tsv
|   |           |-- Combined-gene-level-taxonomy.tsv
|   |           |-- Combined-gene-level-taxonomy_decontam_failure.txt
|   |           |-- Combined-gene-level-taxonomy_decontam_results.tsv
|   |           |-- Combined-gene-level-taxonomy_filtered.tsv
|   |           |-- Combined-gene-level-taxonomy_filtered_heatmap.png
|   |           |-- Combined-gene-level-taxonomy_filtered_top_50_heatmap.png
|   |           |-- Combined-gene-level-taxonomy_unfiltered_heatmap.png
|   |           |-- Combined-gene-level-taxonomy_unfiltered_top_50_heatmap.png
|   |-- predicted-genes
|   |   |-- *-genes.faa
|   |   |-- *-genes.fasta
|   |   |-- *-genes.gff
|   |-- read-mapping
|       |-- *-contig-coverages.tsv
|       |-- *-gene-coverages.tsv
|       |-- *-mapping-info.txt
|       |-- *-metabat-assembly-depth.tsv
|       |-- *.bam
|-- Decontaminated_Sequence_Data # Low-biomass decontaminated reads
|   |-- FastQC | Nanoplots
|   |   |-- *_decontam_*.html
|   |   |-- *_decontam_fastqc.zip # Illumina only
|   |-- Fastq
|   |   |-- *_decontam.fastq.gz
|   |-- Logs
|   |   |-- blank-assembly.log
|   |-- MultiQC_Reports
|       |-- decontam_multiqc.html
|       |-- decontam_multiqc_data.zip
|-- Filtered_Sequence_Data # Quality-filtered reads
|   |-- FastQC
|   |   |-- *_filtered_*.html
|   |   |-- *_filtered_fastqc.zip # Illumina only
|   |-- Fastq
|   |   |-- `Tool Name` [Fastp or Filtlong]
|   |   |   |-- *_filtered.fastq.gz
|   |   |-- Fastp_temp  # First round of Fastp with PollG trimmed - illumina only
|   |       |-- temp_*_filtered.fastq.gz
|   |-- MultiQC_Reports
|       |-- filtered_multiqc.html
|       |-- filtered_multiqc_data.zip
|-- HostRM-removed_Sequence_Data # Host-removed reads - Optional
|   |-- FastQC | Nanoplots 
|   |   |-- *_HostRM_*.html
|   |-- Fastq
|   |   |-- *_HostRM_*.fastq.gz
|   |-- MultiQC_Reports
|       |-- HostRM_multiqc.html
|       |-- HostRM_multiqc_data.zip
|-- HR-removed_Sequence_Data # Human-removed reads - Nanopore data only
|   |-- Nanoplots
|   |   |-- *_HRrm_*.html
|   |-- Fastq
|   |   |-- *_HRrm_*.fastq.gz
|   |-- MultiQC_Reports
|       |-- HRrm_multiqc.html
|       |-- HRrm_multiqc_data.zip
|-- Trimmed_Sequence_Data  # PoreChop Quality-trimmed reads - Nanopore data only
|   |-- Nanoplots
|   |   |-- *_trimmed_*.html
|   |-- Fastq
|   |   |-- *_trimmed.fastq.gz
|   |-- MultiQC_Reports
|       |-- trimmed_multiqc.html
|       |-- trimmed_multiqc_data.zip
|-- Logs
|   |-- <TOOL NAME>
|   |   |-- *.log
|-- Merged_Sequence_Data # Raw read or Merged reads
|   |-- FastQC | Nanoplots
|   |   |-- *.html
|   |   |-- *_fastqc.zip # Illumina only
|   |-- MultiQC_Reports
|       |-- HRrm_multiqc.html
|       |-- HRrm_multiqc_data.zip
|-- Metadata
|   |-- software_versions.txt
|-- Read-based_Processing
|   |-- Humann_Outputs
|   |   |-- Gene_Families
|   |   |   |-- Gene-families-KO-cpm.tsv
|   |   |   |-- Gene-families-KO.tsv
|   |   |   |-- Gene-families-KO_decontam_results.tsv
|   |   |   |-- Gene-families-KO_filtered.tsv
|   |   |   |-- Gene-families-KO_filtered_heatmap.png
|   |   |   |-- Gene-families-KO_filtered_top_50_heatmap.png
|   |   |   |-- Gene-families-KO_unfiltered_heatmap.png
|   |   |   |-- Gene-families-KO_unfiltered_top_50_heatmap.png
|   |   |   |-- Gene-families-cpm.tsv
|   |   |   |-- Gene-families-grouped-by-taxa.tsv
|   |   |   |-- Gene-families-uniref.tsv
|   |   |   |-- Gene-families-uniref_decontam_results.tsv
|   |   |   |-- Gene-families-uniref_filtered.tsv
|   |   |   |-- Gene-families-uniref_filtered_heatmap.png
|   |   |   |-- Gene-families-uniref_filtered_top_50_heatmap.png
|   |   |   |-- Gene-families-uniref_unfiltered_heatmap.png
|   |   |   |-- Gene-families-uniref_unfiltered_top_50_heatmap.png
|   |   |   |-- Gene-families.tsv
|   |   |-- Pathway_Abundances
|   |   |   |-- Pathway-abundances-cpm.tsv
|   |   |   |-- Pathway-abundances-grouped-by-taxa.tsv
|   |   |   |-- Pathway-abundances.tsv
|   |   |   |-- Pathway-abundances_decontam_results.tsv
|   |   |   |-- Pathway-abundances_filtered.tsv
|   |   |   |-- Pathway-abundances_filtered_heatmap.png
|   |   |   |-- Pathway-abundances_unfiltered_heatmap.png
|   |   |   |-- Pathway-abundances_unfiltered_top_50_heatmap.png
|   |   |-- Pathway_Coverage
|   |       |-- Pathway-coverages-grouped-by-taxa.tsv
|   |       |-- Pathway-coverages.tsv
|   |-- Kaiju_Outputs
|   |   |-- Barplots
|   |   |   |-- kaiju_filtered_species_barplot.html
|   |   |   |-- kaiju_filtered_species_barplot.png
|   |   |   |-- kaiju_unfiltered_species_barplot.html
|   |   |   |-- kaiju_unfiltered_species_barplot.png
|   |   |-- Count_tables
|   |   |   |-- kaiju_decontam_results.tsv
|   |   |   |-- kaiju_filtered_species_table.tsv
|   |   |   |-- kaiju_species_table.tsv
|   |   |-- Krona_Reports
|   |       |-- kaiju-report.html
|   |-- Kraken2_Outputs
|   |   |-- Barplots
|   |   |   |-- kraken2_filtered_species_barplot.html
|   |   |   |-- kraken2_filtered_species_barplot.png
|   |   |   |-- kraken2_unfiltered_species_barplot.html
|   |   |   |-- kraken2_unfiltered_species_barplot.png
|   |   |-- Count_tables
|   |   |   |-- kraken2_decontam_results.tsv
|   |   |   |-- kraken2_filtered_species_table.tsv
|   |   |   |-- kraken2_species_table.tsv
|   |   |-- Krona_Reports
|   |   |   |-- kraken2-report.html
|   |   |-- MultiQC_Reports
|   |       |-- kraken2_multiqc.html
|   |       |-- kraken2_multiqc_data.zip
|   |-- Metaphlan_Outputs # Illumina only
|       |-- Barplots
|       |   |-- metaplan_filtered_species_barplot.html
|       |   |-- metaplan_filtered_species_barplot.png
|       |   |-- metaplan_unfiltered_species_barplot.html
|       |   |-- metaplan_unfiltered_species_barplot.png
|       |-- Krona_Reports
|       |   |-- metaphlan-report.html
|       |-- Taxonomy
|           |-- metaphlan-taxonomy.tsv
|           |-- metaphlan_decontam_results.tsv
|           |-- metaphlan_filtered_species_table.tsv
|           |-- metaphlan_species_table.tsv
|-- Resource_Usage
|   |-- execution_report_*.html
|   |-- execution_timeline_*.html
|   |-- execution_trace_*.txt```

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
  │       ├─→ Adapter trimming  
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
  │   └─→ PolyG trimming (NextSeq/NovaSeq)
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
  ├─→ QC: Raw NanoPlot  + MultiQC
  │
  ├─→ Filtering: 
  │   ├─→ Length filtering
  │   └─→ Quality score filtering
  |
  ├─→ QC: Filtered NanoPlot  + MultiQC
  │   
  ├─→ Adapter Trimming: 
  |
  ├─→ QC: Trimmed NanoPlot + MultiQC
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
- `humann3.yaml`
- `cat.yaml`
- `metabat.yaml`
- `gtdbtk.yaml`
- `megahit.yaml`
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
- **Kaiju DB Setup**: 500 GB

---

## Advanced Features

### GeneLab Integration
- Direct data retrieval from GeneLab/OSDR via accession numbers
- Automatic runsheet generation
- GeneLab-specific file naming conventions (`--assay_suffix`, `--additional_filename_prefix`)

### Low-Biomass Sample Processing
Specialized statistical decontamination:
1. Identifies negative controls (NTC = true)
2. Assembles negative controls and maps reads to the assembly
3. Removes contaminant signals/reads from samples by retaining only unmapped reads
4. Statistically identifies and removes contaminant features (taxa or functions) using decontam
5. Applies to both read-based and assembly-based results

### Custom Genome Mapping
Map reads to custom reference genomes:
>Note: You'll need to uncomment relavent lines in main.nf to use this functionality
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
- Filtlong (Nanopore)
- PoreChop (Nanopore)

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
- **GeneLab Metagenomics Workflow**: https://github.com/nasa/GeneLab_Metagenomics_Workflow
- **NASA OSDR**: https://osdr.nasa.gov/
- Individual tool citations (see `software_versions.txt` in outputs)

### Pipeline Documents
- [GL-DPPD-7107-B](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Illumina/Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-B.md) : Standard Illumina metagenomics
- [GL-DPPD-7116](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Pipeline_GL-DPPD-7116_Versions/GL-DPPD-7116.md): Low-biomass Nanopore metagenomics
- [GL-DPPD-7117](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Pipeline_GL-DPPD-7117_Versions/GL-DPPD-7117.md): Low-biomass Illumina metagenomics

---

## Support & Contact

- **Issues**: https://github.com/nasa/GeneLab_Metagenomics_Workflow/issues
- **GeneLab**: https://genelab.nasa.gov/
- **OSDR**: https://osdr.nasa.gov/

---

## Changelog

### Version 1.0.0 (Current - main branch)
- Initial Nextflow DSL2 implementation
- Support for Illumina and Nanopore
- Standard and low-biomass sample types
- Read-based and assembly-based workflows
- Automatic database management
- GeneLab accession integration for traditional metagenomics
- Comprehensive quality control and reporting

---

## License

This workflow is developed by NASA GeneLab and is part of the Open Science Data Repository (OSDR) initiative. Please refer to the repository for license information.
