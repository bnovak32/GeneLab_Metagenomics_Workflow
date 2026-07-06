#!/usr/bin/env python3

import argparse
import re
import sys

def parse_versions(versions_file):
    """Parse software versions from the specified versions file."""
    # Initialize expected tool keys to avoid potential KeyError issues
    tools = [
        'fastqc', 'multiqc', 'bbmap', 'humann', 'megahit', 'prodigal', 'cat',
        'kofamscan', 'bowtie2', 'samtools', 'metabat2', 'bit', 'checkm',
        'gtdbtk', 'fastp', 'metaphlan', 'kraken2', 'kaiju', 'krona',
        'ggplot2', 'pheatmap', 'decontam', 'nanoplot', 'filtlong',
        'porechop', 'flye', 'minimap', 'spades', 'tidyverse'
    ]
    versions = {t: "" for t in tools}
    gtdb_list = []

    try:
        with open(versions_file) as f:
            for line in f:
                line_strip = line.strip()
                if not line_strip:
                    continue
                line_lower = line_strip.lower()
                parts = line_strip.split()
                if len(parts) < 2:
                    continue

                # Mimic the exact awk/sed filtering logic from the shell script
                if 'fastqc' in line_lower:
                    versions['fastqc'] = re.sub(r'^[vV]', '', parts[1])
                elif 'multiqc' in line_lower and len(parts) > 2:
                    versions['multiqc'] = parts[2]
                elif 'bbtools' in line_lower:
                    versions['bbmap'] = parts[1]
                elif 'humann' in line_lower:
                    versions['humann'] = re.sub(r'^[vV]', '', parts[1])
                elif 'megahit' in line_lower:
                    versions['megahit'] = re.sub(r'^[vV]', '', parts[1])
                elif 'prodigal' in line_lower:
                    versions['prodigal'] = re.sub(r'[vV:]', '', parts[1])
                elif 'CAT' in line:  # Case-sensitive check to mirror grep 'CAT'
                    versions['cat'] = re.sub(r'^[vV]', '', parts[1])
                elif 'exec_annotation' in line_lower:
                    versions['kofamscan'] = parts[1]
                elif 'bowtie' in line_lower and len(parts) > 2:
                    versions['bowtie2'] = parts[2]
                elif 'samtools' in line_lower:
                    versions['samtools'] = re.sub(r'^[vV]', '', parts[1])
                elif 'metabat' in line_lower:
                    versions['metabat2'] = parts[1]
                elif 'bioinformatics tools' in line_lower and len(parts) > 2:
                    match = re.search(r'(\d+\.\d+\.\d+)', parts[2])
                    if match:
                        versions['bit'] = match.group(1)
                elif 'checkm' in line_lower:
                    versions['checkm'] = re.sub(r'^[vV]', '', parts[1])
                elif line_strip.startswith('GTDB'):
                    gtdb_list.append(re.sub(r'^[vV]', '', parts[1]))
                elif 'fastp' in line_lower:
                    versions['fastp'] = re.sub(r'^[vV]', '', parts[1])
                elif 'metaphlan' in line_lower:
                    versions['metaphlan'] = re.sub(r'^[vV]', '', parts[1])
                elif 'kraken2' in line_lower:
                    versions['kraken2'] = re.sub(r'^[vV]', '', parts[1])
                elif 'kaiju' in line_lower:
                    versions['kaiju'] = re.sub(r'^[vV]', '', parts[1])
                elif 'krona' in line_lower:
                    versions['krona'] = re.sub(r'^[vV]', '', parts[1])
                elif 'ggplot2' in line_lower:
                    versions['ggplot2'] = re.sub(r'^[vV]', '', parts[1])
                elif 'pheatmap' in line_lower:
                    versions['pheatmap'] = re.sub(r'^[vV]', '', parts[1])
                elif 'decontam' in line_lower:
                    versions['decontam'] = re.sub(r'^[vV]', '', parts[1])
                elif 'nanplot' in line_lower:
                    versions['nanoplot'] = re.sub(r'^[vV]', '', parts[1])
                elif 'filtlong' in line_lower:
                    versions['filtlong'] = re.sub(r'^[vV]', '', parts[1])
                elif 'porechop' in line_lower:
                    versions['porechop'] = re.sub(r'^[vV]', '', parts[1])
                elif 'flye' in line_lower:
                    versions['flye'] = re.sub(r'^[vV]', '', parts[1])
                elif 'minimap2' in line_lower:
                    versions['minimap'] = re.sub(r'^[vV]', '', parts[1])
                elif 'spades' in line_lower:
                    versions['spades'] = re.sub(r'^[vV]', '', parts[1])
                elif 'tidyverse' in line_lower:
                    versions['tidyverse'] = re.sub(r'^[vV]', '', parts[1])

    except FileNotFoundError:
        sys.exit(f"Error: versions file '{versions_file}' not found")

    if gtdb_list:
        # If 2 versions are used, choose the second (replicating bash head -n2 intent)
        versions['gtdbtk'] = gtdb_list[1] if len(gtdb_list) > 1 else gtdb_list[0]

    return versions


def generate_protocol(protocol_id, sample_type, technology, v):
    """Build the protocol paragraphs cleanly by reducing text duplication."""
    
    # 1. Overview
    pipeline_num = "7116" if (technology == "nanopore" and sample_type == "standard") else "7117"
    overview = (
        f"Data were processed as described in {protocol_id} "
        f"(https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_Metagenomics_low_biomass/"
        f"Metagenomics/Low_Biomass/Pipeline_GL-DPPD-{pipeline_num}_Versions/{protocol_id}.md) "
        f"using workflow GeneLab_Metagenomics_Workflow v1.0.0_beta "
        f"(https://github.com/nasa/GeneLab_Metagenomics_Workflow/blob/main/README.md)."
    )

    # 2. Quality Assessment
    if technology == "illumina":
        qa = f"Quality assessment: Quality assessment of raw (human-removed), filtered, and blanks-removed reads was performed with FastQC v{v['fastqc']} and reports summarized with MultiQC v{v['multiqc']}."
    else:  # nanopore
        qa = f"Quality assessment: Quality assessment of raw, filtered, trimmed, blanks-removed and host-removed reads was performed with Nanoplot v{v['nanoplot']} and reports summarized with MultiQC v{v['multiqc']}."

    # 3. Reads Concatenation (Nanopore only)
    reads_concat = ""
    if technology == "nanopore":
        prefix = "Pre-" if sample_type == "standard" else ""
        reads_concat = f"Reads concatenation: {prefix}demultiplexed reads were concatenated by barcode to form raw per sample fastq files."

    # 4. Quality Control
    if technology == "illumina":
        qc = (
            f"Quality control: Two rounds of FASTP v{v['fastp']} were run to quality filter and remove adapters from the human removed raw reads. "
            f"First, raw reads were filtered by quality (Phred quality >= 20) and length (>=50bp), and adapters auto detected then removed using FASTP v{v['fastp']}. "
            f"Next, PolyG adapters were detected and removed by running FASTP again on the quality filtered reads after setting FASTP's --trim_poly_g flag."
        )
    else:  # nanopore
        qc = (
            f"Quality control: raw concatenated reads were filtered by quality (Phred quality >= 8) and length (>=200bp) using filtlong v{v['filtlong']}. "
            f"Adapters were then removed from the filtered reads using porechop v{v['porechop']}."
        )

    # 5. Human Reads Removal (Nanopore only)
    human_removal = ""
    if technology == "nanopore":
        human_removal = (
            f"Human reads removal: Human reads were removed from the quality controlled reads using kraken2 v{v['kraken2']}. "
            f"In short, human reads were identified and removed from the raw reads using kraken2 v{v['kraken2']} against a human genome reference database "
            f"that was constructed from NCBI's RefSeq (GCF_000001405.39) GRCh38.p13. The database was constructed by running kraken2 build command "
            f"with the following parameter set --no-masking, kmer-length 35 and minimizer-length 31."
        )

    # 6. Blanks Reads Removal (Low Biomass only)
    blanks_removal = ""
    if sample_type == "low_biomass":
        if technology == "illumina":
            blanks_removal = (
                f"Blanks reads removal: negative control samples were assembled using Spades v{v['spades']}. "
                f"Quality controlled reads were decontaminated by mapping them to the assembled blank contigs then "
                f"unmapped/uncontaminated reads filtered out using bowtie2 v{v['bowtie2']} by setting the parameters "
                f"--very-sensitive-local and --un-conc-gz to retain only unmapped reads."
            )
        else:  # nanopore
            blanks_removal = (
                f"Blanks reads removal: Blank / negative control samples were assembled using Flye v{v['flye']}. "
                f"Human removed reads were decontaminated by mapping them to the assembled blank contigs using minimap2 v{v['minimap']} "
                f"then unmapped/uncontaminated reads filtered out using samtools v{v['samtools']} fastq by setting the -f parameter to 4 to retain only unmapped reads."
            )

    # Determine read types used for subsequent processing steps
    if sample_type == "low_biomass":
        read_type = "decontaminated reads"
    else:  # standard
        read_type = "quality controlled reads" if technology == "illumina" else "human-removed reads"

    # 7. Read-Based Processing
    if technology == "illumina":
        rb_tools = f"Metaphlan v{v['metaphlan']}, Kraken2 v{v['kraken2']} and Kaiju v{v['kaiju']}"
    else:  # nanopore
        rb_tools = f"Kraken2 v{v['kraken2']} and Kaiju v{v['kaiju']}"
        
    read_processing = (
        f"Read-based processing: taxonomy assignment of {read_type} was performed with {rb_tools}. "
        f"The resulting assignments were visualized as krona plots using krona v{v['krona']} and as barplots in R using tidyverse v{v['tidyverse']}."
    )

    # 8. Assembly-Based Analysis
    assembly_read_type = read_type
    if sample_type == "standard" and technology == "illumina":
        assembly_read_type = "Quality controlled reads"  # Handles exact capitalization discrepancy

    assembly_tool = f"megahit v{v['megahit']}" if technology == "illumina" else f"Flye v{v['flye']}"
    mapping_tool = f"bowtie2 v{v['bowtie2']}" if technology == "illumina" else f"minimap2 v{v['minimap']}"

    assembly_analysis = (
        f"Assembly-based analysis: {assembly_read_type} were assembled with {assembly_tool}. "
        f"Genes were called with prodigal v{v['prodigal']}. Taxonomic classification of genes and contigs was performed with CAT v{v['cat']}. "
        f"Functional annotation was done with KOFamScan v{v['kofamscan']}. Reads were mapped to assemblies with {mapping_tool}, "
        f"and coverage information was extracted for reads and contigs with samtools v{v['samtools']} and bbmap v{v['bbmap']}. "
        f"Binning of contigs was performed with metabat2 v{v['metabat2']}. Bins were summarized with bit v{v['bit']} and "
        f"estimates of bin quality were generated with checkm v{v['checkm']}. High-quality bins (greater than 90% est. completeness "
        f"and less than 10% est. redundancy) were taxonomically classified with gtdb-tk v{v['gtdbtk']} as MAGs (Meta assembled genomes)."
    )

    # 9. Feature Table Visualizations / Decontamination
    if sample_type == "low_biomass":
        visualizations = (
            f"Feature table decontamination and visualizations: All graphical displays were conducted using R. "
            f"Bar plots and heatmaps were generated using R packages ggplot2 v{v['ggplot2']} and pheatmap v{v['pheatmap']}, respectively. "
            f"Taxonomy and function annotation tables from both read and assembly-based processing steps were decontaminated with "
            f"decontam v{v['decontam']} an R package designed to statistically identify contaminant features in a feature table."
        )
    else:  # standard
        visualizations = (
            f"Feature table visualizations: All graphical displays were conducted using R. "
            f"Bar plots and heatmaps were generated using R packages ggplot2 v{v['ggplot2']} and pheatmap v{v['pheatmap']}, respectively."
        )

    # Combine parts and strip extra internal whitespace cleanly
    parts = [overview, qa, reads_concat, qc, human_removal, blanks_removal, read_processing, assembly_analysis, visualizations]
    return " ".join([p for p in parts if p])


def main():
    parser = argparse.ArgumentParser(description="Generate metagenomics processing protocol text.")
    parser.add_argument('--versions-file', help='Path to software_versions.txt file')
    parser.add_argument('--protocol-id', help='Protocol ID (e.g. GL-DPPD-7107-A)')
    parser.add_argument('--sample-type', choices=['low_biomass', 'standard'], help='Sample structural biomass type')
    parser.add_argument('--technology', choices=['illumina', 'nanopore'], help='Sequencing technology framework')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    versions = parse_versions(args.versions_file)
    protocol = generate_protocol(args.protocol_id, args.sample_type, args.technology, versions)
    
    # Print to standard output to match the original echo behavior
    print(protocol)


if __name__ == "__main__":
    main()
