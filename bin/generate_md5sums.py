#!/usr/bin/env python3

from email import parser
import os
import sys
import hashlib
import argparse

def calculate_md5(filepath):
    """Calculate MD5 hash for a file."""
    md5_hash = hashlib.md5()
    
    # Follow symlinks to get the actual file
    actual_path = os.path.realpath(filepath) if os.path.islink(filepath) else filepath
    
    try:
        with open(actual_path, "rb") as f:
            # Read in chunks in case of large files
            for chunk in iter(lambda: f.read(4096), b""):
                md5_hash.update(chunk)
        return md5_hash.hexdigest()
    except Exception as e:
        sys.stderr.write(f"Error calculating MD5 for {filepath}: {str(e)}\n")
        return "ERROR"

def is_raw_file(filepath):
    """Check if file is a raw FASTQ or raw multiqc file (zip or html)."""
    if "/Merged_Sequence_Data/" in filepath:
        # Match raw fastq files
        if filepath.endswith(".fastq.gz"):
            return True
        # Match raw multiqc reports (zip or html)
        elif "raw_multiqc" in filepath and (filepath.endswith(".zip") or filepath.endswith(".html")):
            return True
    return False

def should_include(filepath, allowed_dirs):
    """Check if file should be included in MD5 calculation."""
    
    # Skip any files not in the list of workflow output parent directories
    if not any(allowed_dir in filepath for allowed_dir in allowed_dirs):
        return False
    
    # Skip versions.txt, *decontam_failure.txt files from all directories
    if filepath.endswith("versions.txt") or filepath.endswith("decontam_failure.txt"):
        return False

    # Skip any individual log files
    if ("/Logs/" in filepath and filepath.endswith(".log")):
        return False
    
    # Skip any files with 'fastqc' in the path or filename (case-sensitive to make sure Post_Processing/FastQC_Outputs/ is still included) 
    if 'fastqc' in filepath:
        return False

    # Skip all .png files in alpha Barplot outputs
    if ("/Barplots/" in filepath and filepath.endswith(".png")):
        return False

    # Skip ISA.zip and GLfile.csv files
    if filepath.endswith("ISA.zip") or filepath.endswith("GLfile.csv"):
        return False

    # Skip files in Post_Processing except for README and processing_info.zip
    if "/Post_Processing/" in filepath and not (filepath.endswith("README" + args.assay_suffix + ".txt") or filepath.endswith("processing_info" + args.assay_suffix + ".zip")):
        return False
    
    return True

def trim_to_allowed_dir(filepath, allowed_dirs):
    """Trim filepath to start from the first matching allowed_dir segment."""
    for allowed_dir in allowed_dirs:
        # Strip leading slash for the match so the output starts like: FastQC_Outputs/...
        marker = allowed_dir.strip("/")
        if marker in filepath:
            idx = filepath.index(marker)
            return filepath[idx:]
    return os.path.basename(filepath)  # fallback: just the filename

def main():
    parser = argparse.ArgumentParser(description='Generate MD5 sum files for GeneLab data.')
    parser.add_argument('--outdir', required=True, help='Output directory containing files to process')
    parser.add_argument('--assay_suffix', default='', help="GeneLab assay suffix. one of '_GLmetagenomics', '_GLlbsMetag', '_GLlblMetag', and '_GLlongMetag' (default: empty)")
    parser.add_argument('--output_prefix', default='', help='Prefix for output MD5 sum files (default: empty)')
    parser.add_argument('--generate_raw_md5sums', action='store_true', help='Generate raw MD5sums')
    
    global args
    args = parser.parse_args()
    
    # Make sure outdir is absolute path
    outdir = os.path.abspath(args.outdir)
    
    # Create output files and initialize them without headers
    if args.generate_raw_md5sums:
        raw_md5_file = f"{args.output_prefix}raw_md5sum{args.assay_suffix}.tsv"
        with open(raw_md5_file, 'w') as f:
            f.write("File Path\tFile Name\tmd5\n")
    
    processed_md5_file = f"{args.output_prefix}processed_md5sum{args.assay_suffix}.tsv"
    with open(processed_md5_file, 'w') as f:
        f.write("File Path\tFile Name\tmd5\n")
    
    
    # Track processed files for reporting
    raw_count = 0
    processed_count = 0
    allowed_dirs = [
            "/Metadata/", 
            "/GeneLab/", 
            "/Merged_Sequence_Data/", 
            "/Filtered_Sequence_Data/",
            "Trimmed_Sequence_Data/",
            "/Decontaminated_Sequence_Data/",
            "/HR-removed_Sequence_Data/",
            "HostRm-removed_Sequence_Data/",
            "/Assembly-based_Processing/",
            "/Read-based_Processing/",  
            "/Post_Processing/"]
 
    
    # Walk through all files recursively
    print(f"Scanning directory: {outdir}")
    for root, _, files in os.walk(outdir):
        for filename in files:
            filepath = os.path.join(root, filename)
            
            # Skip files that shouldn't be included
            if not should_include(filepath, allowed_dirs):
                continue
            
            # Get just the filename (basename)
            basename = os.path.basename(filepath)
            trimmed_path = trim_to_allowed_dir(filepath, allowed_dirs)
            
            # Generate and save md5sum of raw files only if enabled
            if is_raw_file(filepath):
                if not args.generate_raw_md5sums:
                    continue

                md5sum = calculate_md5(filepath)
                with open(raw_md5_file, 'a') as f:
                    f.write(f"{trimmed_path}\t{basename}\t{md5sum}\n")
                raw_count += 1
            else:
                #Generate and save md5sum of processed files always
                md5sum = calculate_md5(filepath)
                with open(processed_md5_file, 'a') as f:
                    f.write(f"{trimmed_path}\t{basename}\t{md5sum}\n")
                processed_count += 1
    
    if args.generate_raw_md5sums:
        print(f"Added {raw_count} files to {raw_md5_file}")
    print(f"Added {processed_count} files to {processed_md5_file}")

    def dedup_file(filename):
        seen = set()
        lines = []
        with open(filename, 'r') as f:
            for line in f:
                key = line.split('\t', 1)[0]  # dedup by basename
                if key not in seen:
                    seen.add(key)
                    lines.append(line)
        with open(filename, 'w') as f:
            f.writelines(lines)

    if args.generate_raw_md5sums:
        dedup_file(raw_md5_file)
    dedup_file(processed_md5_file)

if __name__ == "__main__":
    main()
