#!/usr/bin/env bash

# Generate protocol according to a pipeline document

# USAGE:
# generate_protocol.sh <software_versions> <protocol_id> <sample_type> <technology>
# EXAMPLE
# generate_protocol.sh ../Metadata/software_versions.txt GL-DPPD-7107-A low_biomass illumina

FASTQC=`grep -i 'fastqc' $1 | awk '{print $2}' |sed -E 's/v//'`
MULTIQC=`grep -i 'multiqc' $1 | awk '{print $3}'`
BBMAP=`grep -i 'bbtools' $1 | awk '{print $2}'`
HUMANN=`grep -i 'humann' $1 | awk '{print $2}'|sed -E 's/v//'`
MEGAHIT=`grep -i 'megahit' $1 | awk '{print $2}'|sed -E 's/v//'`
PRODIGAL=`grep -i 'prodigal' $1 | awk '{print $2}'|sed -E 's/[vV:]//g'`
CAT=`grep 'CAT' $1 | awk '{print $2}'|sed -E 's/v//'`
KOFAMSCAN=`grep 'exec_annotation' $1 | awk '{print $2}'`
BOWTIE2=`grep -i 'bowtie' $1 | awk '{print $3}'`
SAMTOOLS=`grep -i 'samtools' $1 | awk '{print $2}'`
METABAT2=`grep -i 'metabat' $1 | awk '{print $2}'`
BIT=`grep -i 'bioinformatics tools' $1 | awk '{print $3}' | sed 's/v//' | sed -E 's/.+([0-9]+.[0-9]+.[0-9]+).+/\1/'`
CHECKM=`grep -i 'checkm' $1 | awk '{print $2}' |sed -E 's/v//'`
GTDBTK=`grep -i '^GTDB' $1 | awk '{print $2}' |sed -E 's/v//' | head -n2` # If 2 versions are used, choose the second
FASTP=`grep -i 'fastp' $1 | awk '{print $2}' |sed -E 's/v//'`
METAPHLAN=`grep -i 'metaphlan' $1 | awk '{print $2}' |sed -E 's/v//'`
KRAKEN2=`grep -i 'kraken2' $1 | awk '{print $2}' |sed -E 's/v//'`
KAIJU=`grep -i 'kaiju' $1 | awk '{print $2}' |sed -E 's/v//'`
KRONA=`grep -i 'krona' $1 | awk '{print $2}' |sed -E 's/v//'`
GGPPLOT2=`grep -i 'ggplot2' $1 | awk '{print $2}' |sed -E 's/v//'`
PHEATMAP=`grep -i 'pheatmap' $1 | awk '{print $2}' |sed -E 's/v//'`
DECONTAM=`grep -i 'decontam' $1 | awk '{print $2}' |sed -E 's/v//'`
NANOPLOT=`grep -i 'NanoPlot' $1 | awk '{print $2}' |sed -E 's/v//'`
FILTLONG=`grep -i 'Filtlong' $1 | awk '{print $2}' |sed -E 's/v//'`
PORECHOP=`grep -i 'Porechop' $1 | awk '{print $2}' |sed -E 's/v//'`
FLYE=`grep -i 'Flye' $1 | awk '{print $2}' |sed -E 's/v//'`
MINIMAP=`grep -i 'minimap2' $1 | awk '{print $2}' |sed -E 's/v//'`


PROTOCOL_ID=$2
SAMPLE_TYPE=$3
TECHNOLOGY=$4

if [[ ${SAMPLE_TYPE} == "low_biomass" && ${TECHNOLOGY} == "illumina" ]]; then

# Illumina low biomass
PROTOCOL="Data were processed as described in ${PROTOCOL_ID} (https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Pipeline_GL-DPPD-7117_Versions/${PROTOCOL_ID}.md) \
using workflow GeneLab_Metagenomics_Workflow v1.0.0_beta (https://github.com/nasa/GeneLab_Metagenomics_Workflow/blob/main/README.md). \
Quality assessment: Quality assessment of raw (human-removed), filtered, and blanks-removed reads was performed with FastQC v${FASTQC} and reports summarized with MultiQC v${MULTIQC}. \
Quality control: Two rounds of FASTP v${FASTP} were run to quality filter and remove adapters from the human removed raw reads. First, raw reads were filtered by quality (Phred quality >= 20) and length (>=50bp), \
and adapters auto detected then removed using FASTP v${FASTP}. Next, PolyG adapters were detected and removed by running FASTP again on the quality filtered reads after setting FASTP's --trim_poly_g flag. \
Blanks reads removal: negative control samples were assembled using Spades v${SPADES}. Quality controlled reads were decontaminated by mapping them to the assembled blank contigs then \
unmapped/uncontaminated reads filtered out using bowtie2 v${BOWTIE2} by setting the parameters --very-sensitive-local and --un-conc-gz to retain only unmapped reads. \
Read-based processing: taxonomy assignment of decontaminated reads was performed with Metaphlan v${METAPHLAN}, Kraken2 v${KRAKEN2} and Kaiju v${KAIJU}. \
The resulting assignments were visualized as krona plots using krona v${KRONA} and as barplots in R using tidyverse v${TIDYVERSE}. \
Assembly-based analysis: decontaminated reads were assembled with megahit v${MEGAHIT}. Genes were called with prodigal v${PRODIGAL}. Taxonomic classification of genes and contigs was performed with CAT v${CAT}. \
Functional annotation was done with KOFamScan v${KOFAMSCAN}. Reads were mapped to assemblies with bowtie2 v${BOWTIE2}, and coverage information was extracted for reads and contigs with samtools v${SAMTOOLS} and bbmap v${BBMAP}. \
Binning of contigs was performed with metabat2 v${METABAT2}. Bins were summarized with bit v${BIT} and estimates of bin quality were generated with checkm v${CHECKM}. \
High-quality bins (greater than 90% est. completeness and less than 10% est. redundancy) were taxonomically classified with gtdb-tk v${GTDBTK} as MAGs (Meta assembled genomes). \
Feature table decontamination and visualizations: All graphical displays were conducted using R. Bar plots and heatmaps were generated using R packages ggplot2 v${GGPPLOT2} and pheatmap v${PHEATMAP}, respectively. \
Taxonomy and function annotation tables from both read and assembly-based processing steps were decontaminated with decontam v${DECONTAM} an R package designed to statistically identify contaminant features in a feature table."

elseif [[ ${SAMPLE_TYPE} == "standard" && ${TECHNOLOGY} == "illumina" ]]; then

# Illumina standard
PROTOCOL="Data were processed as described in ${PROTOCOL_ID} (https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Pipeline_GL-DPPD-7117_Versions/${PROTOCOL_ID}.md) \
using workflow GeneLab_Metagenomics_Workflow v1.0.0_beta (https://github.com/nasa/GeneLab_Metagenomics_Workflow/blob/main/README.md). \
Quality assessment: Quality assessment of raw (human-removed), filtered, and blanks-removed reads was performed with FastQC v${FASTQC} and reports summarized with MultiQC v${MULTIQC}. \
Quality control: Two rounds of FASTP v${FASTP} were run to quality filter and remove adapters from the human removed raw reads. First, raw reads were filtered by quality (Phred quality >= 20) and length (>=50bp), \
and adapters auto detected then removed using FASTP v${FASTP}. Next, PolyG adapters were detected and removed by running FASTP again on the quality filtered reads after setting FASTP's --trim_poly_g flag. \
Read-based processing: taxonomy assignment of quality controlled reads was performed with Metaphlan v${METAPHLAN}, Kraken2 v${KRAKEN2} and Kaiju v${KAIJU}. \
The resulting assignments were visualized as krona plots using krona v${KRONA} and as barplots in R using tidyverse v${TIDYVERSE}. \
Assembly-based analysis: Quality controlled reads were assembled with megahit v${MEGAHIT}. Genes were called with prodigal v${PRODIGAL}. Taxonomic classification of genes and contigs was performed with CAT v${CAT}. \
Functional annotation was done with KOFamScan v${KOFAMSCAN}. Reads were mapped to assemblies with bowtie2 v${BOWTIE2}, and coverage information was extracted for reads and contigs with samtools v${SAMTOOLS} and bbmap v${BBMAP}. \
Binning of contigs was performed with metabat2 v${METABAT2}. Bins were summarized with bit v${BIT} and estimates of bin quality were generated with checkm v${CHECKM}. \
High-quality bins (greater than 90% est. completeness and less than 10% est. redundancy) were taxonomically classified with gtdb-tk v${GTDBTK} as MAGs (Meta assembled genomes). \
Feature table visualizations: All graphical displays were conducted using R. Bar plots and heatmaps were generated using R packages ggplot2 v${GGPPLOT2} and pheatmap v${PHEATMAP}, respectively."

elseif [[ ${SAMPLE_TYPE} == "low_biomass" && ${TECHNOLOGY} == "nanopore" ]]; then

# Nanopore low biomass
PROTOCOL="Data were processed as described in ${PROTOCOL_ID} (https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Pipeline_GL-DPPD-7117_Versions/${PROTOCOL_ID}.md) 
using workflow GeneLab_Metagenomics_Workflow v1.0.0_beta (https://github.com/nasa/GeneLab_Metagenomics_Workflow/blob/main/README.md). \
Quality assessment: Quality assessment of raw, filtered, trimmed, blanks-removed and host-removed reads was performed with Nanoplot v${NANOPLOT} and reports summarized with MultiQC v${MULTIQC}. \
Reads concatenation: demultiplexed reads were concatenated by barcode to form raw per sample fastq files.\
Quality control: raw concatenated reads were filtered by quality (Phred quality >= 8) and length (>=200bp) using filtlong v${FILTLONG}. Adapters were then removed from the filtered reads using porechop v${PORECHOP}. \
Human reads removal: Human reads were removed from the quality controlled reads using kraken2 v${KRAKEN2}. \
In short, human reads were identified and removed from the raw reads using kraken2 v${KRAKEN2} against a human genome reference database that was constructed from NCBI's RefSeq (GCF_000001405.39) GRCh38.p13. \
The database was constructed by running kraken2 build command with the following parameter set --no-masking, kmer-length 35 and minimizer-length 31. \
Blanks reads removal: Blank / negative control samples were assembled using Flye v${FLYE}. Human removed reads were decontaminated by mapping them to the assembled blank contigs using minimap2 v${MINIMAP} \
then unmapped/uncontaminated reads filtered out using samtools v${SAMTOOLS} fastq by setting the -f parameter to 4 to retain only unmapped reads. 
Read-based processing: taxonomy assignment of decontaminated reads was performed with Kraken2 v${KRAKEN2} and Kaiju v${KAIJU}. \
The resulting assignments were visualized as krona plots using krona v${KRONA} and as barplots in R using tidyverse v${TIDYVERSE}.
Assembly-based analysis: decontaminated reads were assembled with Flye v${FLYE}. Genes were called with prodigal v${PRODIGAL}. \
Taxonomic classification of genes and contigs was performed with CAT v${CAT}. Functional annotation was done with KOFamScan v${KOFAMSCAN}. \
Reads were mapped to assemblies with minimap2 v${MINIMAP}, and coverage information was extracted for reads and contigs with samtools v${SAMTOOLS} and bbmap v${BBMAP}. \
Binning of contigs was performed with metabat2 v${METABAT}. Bins were summarized with bit v${BIT} and estimates of bin quality were generated with checkm v${CHECKM}. \
High-quality bins (greater than 90% est. completeness and less than 10% est. redundancy) were taxonomically classified with gtdb-tk v${GTDB_TK} as MAGs (Meta assembled genomes). \
Feature table decontamination and visualizations: All graphical displays were conducted using R. Bar plots and heatmaps were generated using R packages ggplot2 v${GGPLOT2} and pheatmap v${PHEATMAP}, respectively. \
Taxonomy and function annotation tables from both read and assembly-based processing steps were decontaminated with decontam v${DECONTAM} an R package designed to statistically identify contaminant features in a feature table."

else

# Nanopore standard
PROTOCOL="Data were processed as described in ${PROTOCOL_ID} (https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/Metagenomics/Low_Biomass/Pipeline_GL-DPPD-7116_Versions/${PROTOCOL_ID}.md) 
using workflow GeneLab_Metagenomics_Workflow v1.0.0_beta (https://github.com/nasa/GeneLab_Metagenomics_Workflow/blob/main/README.md). \
Quality assessment: Quality assessment of raw, filtered, trimmed, blanks-removed and host-removed reads was performed with Nanoplot v${NANOPLOT} and reports summarized with MultiQC v${MULTIQC}. \
Reads concatenation: Pre-demultiplexed reads were concatenated by barcode to form raw per sample fastq files. \
Quality control: raw concatenated reads were filtered by quality (Phred quality >= 8) and length (>=200bp) using filtlong v${FILTLONG}. Adapters were then removed from the filtered reads using porechop v${PORECHOP}. \
Human reads removal: Human reads were removed from the quality controlled reads using kraken2 v${KRAKEN2}. \
In short, human reads were identified and removed from the raw reads using kraken2 v${KRAKEN2} against a human genome reference database that was constructed from NCBI's RefSeq (GCF_000001405.39) GRCh38.p13. \
The database was constructed by running kraken2 build command with the following parameter set --no-masking, kmer-length 35 and minimizer-length 31. \
Read-based processing: taxonomy assignment of human-removed reads was performed with Kraken2 v${KRAKEN2} and Kaiju v${KAIJU}. \
The resulting assignments were visualized as krona plots using krona v${KRONA} and as barplots in R using tidyverse v${TIDYVERSE}. \
Assembly-based analysis: human-removed reads were assembled with Flye v${FLYE}. Genes were called with prodigal v${PRODIGAL}. Taxonomic classification of genes and contigs was performed with CAT v${CAT}. \
Functional annotation was done with KOFamScan v${KOFAMSCAN}. Reads were mapped to assemblies with minimap2 v${MINIMAP}, and coverage information was extracted for reads and contigs with samtools v${SAMTOOLS} and bbmap v${BBMAP}. \
Binning of contigs was performed with metabat2 v${METABAT}. Bins were summarized with bit v${BIT} and estimates of bin quality were generated with checkm v${CHECKM}. \
High-quality bins (greater than 90% est. completeness and less than 10% est. redundancy) were taxonomically classified with gtdb-tk v${GTDB_TK} as MAGs (Meta assembled genomes). \
Feature table visualizations: All graphical displays were conducted using R. Bar plots and heatmaps were generated using R packages ggplot2 v${GGPLOT2} and pheatmap v${PHEATMAP}, respectively."

fi

echo ${PROTOCOL}