#!/usr/bin/env python3

import argparse
import re
import sys
import jinja2

def parse_versions(versions_file):
    """Parse software versions from the specified versions file."""
    # Initialize expected tool keys to avoid potential KeyError issues
    tools = [
        'fastqc', 'multiqc', 'bbmap', 'humann', 'megahit', 'prodigal', 'cat',
        'kofamscan', 'bowtie2', 'samtools', 'metabat2', 'bit', 'checkm',
        'gtdbtk', 'fastp', 'metaphlan', 'kraken2', 'kaiju', 'krona',
        'ggplot2', 'pheatmap', 'decontam', 'nanoplot', 'filtlong',
        'porechop', 'flye', 'minimap', 'spades', 'tidyverse', 'r', 'workflow'
    ]
    versions = {t: "" for t in tools}
    gtdb_list = []

    versions['r'] = "4.5.3"

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
                elif 'nanplot' in line_lower or 'nanoplot' in line_lower:
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
                elif 'metagenomics' in line_lower:
                    versions['workflow'] = re.sub(r'^[vV]', '', parts[1])

    except FileNotFoundError:
        sys.exit(f"Error: versions file '{versions_file}' not found")

    if gtdb_list:
        # If 2 versions are used, choose the second (replicating bash head -n2 intent)
        versions['gtdbtk'] = gtdb_list[1] if len(gtdb_list) > 1 else gtdb_list[0]

    return versions


def generate_protocol_jinja(protocol_id: str, sample_type: str, technology: str, v: dict, 
                            kraken2_genome_reference: str = "(GCF_000001405.39) GRCh38.p13",
                            concat:bool = True):
    """Build the protocol paragraphs cleanly by reducing text duplication."""
    
    # Determine read types used for subsequent processing steps
    if sample_type == "low_biomass":
        read_type = "decontaminated"
    else:  # standard
        read_type = "quality controlled" if technology == "illumina" else "human-removed"


    jinja_env = jinja2.Environment(loader=jinja2.FileSystemLoader("./templates/"), 
                                   trim_blocks=True, lstrip_blocks=True, keep_trailing_newline=False)

    viz_template = jinja_env.get_template("Visualizations_template.txt")

    protocol = ""

    if technology == "illumina":
        if sample_type == "standard":
            pipeline_number = "7107"
            pipeline_version = "-B"
            workflow_type = "illumina"
        else:
            pipeline_number = "7117"
            pipeline_version = ""
            workflow_type = "Low_Biomass"

        template = jinja_env.get_template("Protocol_template_illumina.txt")
        protocol = template.render(
            sample_type=sample_type,
            read_type=read_type,
            workflow_type=workflow_type,
            pipeline_number=pipeline_number,
            pipeline_version=pipeline_version,
            workflow_version=v['workflow'],
            fastqc_version=v['fastqc'],
            multiqc_version=v['multiqc'],
            fastp_version=v['fastp'],
            spades_version=v['spades'],
            bowtie2_version=v['bowtie2'],
            metaphlan_version=v['metaphlan'],
            kraken2_version=v['kraken2'],
            kaiju_version=v['kaiju'],
            krona_version=v['krona'],
            r_version=v['r'],
            tidyverse_version=v['tidyverse'],
            megahit_version=v['megahit'],
            prodigal_version=v['prodigal'],
            cat_version=v['cat'],
            kofamscan_version=v['kofamscan'],
            samtools_version=v['samtools'],
            bbtools_version=v['bbmap'],
            metabat_version=v['metabat2'],
            bit_version=v['bit'],
            checkm_version=v['checkm'],
            gtdbtk_version=v['gtdbtk'],
        )

    else:
        pipeline_version=""
        workflow_type=""
        pipeline_number=""
        if sample_type == "low_biomass":
            pipeline_number = "7116"
            workflow_type = "Low_Biomass"

        template = jinja_env.get_template("Protocol_template_nanopore.txt")
        
        protocol = template.render(
            sample_type=sample_type,
            read_type=read_type,
            workflow_type=workflow_type,
            pipeline_number=pipeline_number,
            pipeline_version=pipeline_version,
            workflow_version=v['workflow'],
            nanoplot_version=v['nanoplot'],
            multiqc_version=v['multiqc'],
            kraken2_version=v['kraken2'],
            minimap2_version=v['minimap'],
            flye_version=v['flye'],
            kaiju_version=v['kaiju'],
            krona_version=v['krona'],
            r_version=v['r'],
            tidyverse_version=v['tidyverse'],
            megahit_version=v['megahit'],
            prodigal_version=v['prodigal'],
            cat_version=v['cat'],
            kofamscan_version=v['kofamscan'],
            samtools_version=v['samtools'],
            bbtools_version=v['bbmap'],
            metabat_version=v['metabat2'],
            bit_version=v['bit'],
            checkm_version=v['checkm'],
            gtdbtk_version=v['gtdbtk'],
            ggplot2_version=v['ggplot2'],
            pheatmap_version=v['pheatmap'],
            decontam_version=v['decontam'],
            human_genome_version=kraken2_genome_reference,
            concat=concat
        )

    visualizations = viz_template.render(
        sample_type=sample_type,
        r_version=v['r'],
        ggplot2_version=v['ggplot2'],
        pheatmap_version=v['pheatmap'],
        decontam_version=v['decontam']
    )

    return " ".join([protocol, visualizations]).replace('\n', '')


def main():
    parser = argparse.ArgumentParser(description="Generate metagenomics processing protocol text.")
    parser.add_argument('--versions-file', help='Path to software_versions.txt file')
    parser.add_argument('--protocol-id', help='Protocol ID (e.g. GL-DPPD-7107-A)')
    parser.add_argument('--sample-type', choices=['low_biomass', 'standard'], help='Sample structural biomass type')
    parser.add_argument('--technology', choices=['illumina', 'nanopore'], help='Sequencing technology framework')
    parser.add_argument('--concat-reads', type=bool, default=False, help='Flag indicating that reads were concatenated before processing (specific to Nanopore workflow).')
    parser.add_argument('--host-removed', type=bool, default=False, help='Flag indicating if host reads were removed.')
    parser.add_argument('--kraken2-genome-reference', type=str, default="(GCF_000001405.39) GRCh38.p13", help='Human genome reference information')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    versions = parse_versions(args.versions_file)
    protocol = generate_protocol_jinja(args.protocol_id, args.sample_type, args.technology, versions, 
                                       kraken2_genome_reference=args.kraken2_genome_reference,
                                       concat=args.concat_reads, )
    
    # Print to standard output to match the original echo behavior
    print(protocol)


if __name__ == "__main__":
    main()
