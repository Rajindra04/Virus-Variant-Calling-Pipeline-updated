#!/usr/bin/env python
import sys
import os
import argparse

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(description='Generate sample sheet from fastq files')
    parser.add_argument('directory', type=str, help='Directory containing fastq files')
    parser.add_argument('sample_sheet_file', type=str, help='Output sample sheet file')
    parser.add_argument('--sample_names', type=str, default=None,
                        help='Comma-separated list of sample names (must match number of FASTQ pairs)')
    args = parser.parse_args(argv)

    directory = args.directory
    sample_sheet_file = args.sample_sheet_file

    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        sys.exit(1)

    if not os.access(os.path.dirname(sample_sheet_file) or '.', os.W_OK):
        print(f"Error: No write permission for {os.path.dirname(sample_sheet_file)}")
        sys.exit(1)

    sample_names = None
    if args.sample_names:
        sample_names = [s.strip() for s in args.sample_names.split(',')]

    generate_sample_sheet(directory, sample_sheet_file, sample_names)

def get_fastq_files(directory):
    fastq_files = []
    for file in os.listdir(directory):
        if file.endswith((".fastq.gz", ".fq.gz")):
            fastq_files.append(os.path.join(directory, file))
    return fastq_files

def find_pair(fastq, fastq_files):
    base = os.path.basename(fastq)
    for suffix in ["_R1_", "_R1.", "_1_", "_1.", "_forward_"]:
        if suffix in base:
            read1 = fastq
            read2_base = base.replace(suffix, suffix.replace("1", "2").replace("forward", "reverse"))
            read2 = os.path.join(os.path.dirname(fastq), read2_base)
            if read2 in fastq_files:
                return read1, read2
    return None, None

def generate_sample_sheet(directory, sample_sheet_file, sample_names=None):
    abs_directory = os.path.abspath(directory)
    fastq_files = get_fastq_files(abs_directory)
    if not fastq_files:
        print(f"Error: No FASTQ files found in {directory}")
        sys.exit(1)

    # Collect all pairs first
    pairs = []
    processed = set()
    for fastq in sorted(fastq_files):
        if fastq in processed:
            continue
        read1, read2 = find_pair(fastq, fastq_files)
        if read1 and read2:
            pairs.append((read1, read2))
            processed.add(read1)
            processed.add(read2)

    if not pairs:
        print(f"Error: No valid FASTQ pairs found in {directory}")
        sys.exit(1)

    # Validate user-provided names if given
    if sample_names:
        if len(sample_names) != len(pairs):
            print(f"Error: {len(sample_names)} sample names provided but {len(pairs)} FASTQ pairs found")
            sys.exit(1)

    with open(sample_sheet_file, "w") as f:
        f.write("sample_name\tread1\tread2\n")
        for i, (read1, read2) in enumerate(pairs):
            if sample_names:
                name = sample_names[i]
            else:
                name = os.path.basename(read1).split("_")[0]
            f.write(f"{name}\t{read1}\t{read2}\n")

if __name__ == "__main__":
    main()
