import sys
import argparse
import glob
import os
import subprocess
import logging

from virus_pipeline.config import load_config

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command):
    logging.info(f"Running command: {command}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command execution failed with return code {process.returncode}, stderr: {stderr.decode('utf-8')}")
    return stdout.decode("utf-8"), stderr.decode("utf-8")

def validate_bam(bam_file):
    logging.info(f"Validating BAM file: {bam_file}")
    run_command(f"samtools quickcheck {bam_file}")
    stdout, stderr = run_command(f"samtools view -c {bam_file}")
    alignments = int(stdout.strip())
    logging.info(f"BAM file {bam_file} validated with {alignments} alignments")
    return alignments > 0

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True, help='Path to input directory containing SAM files.')
    parser.add_argument('--reference_fasta', type=str, required=True, help='Path to reference FASTA file.')
    parser.add_argument('--output_dir', type=str, required=True, help='Path to output directory for BAM files.')
    parser.add_argument('--config', type=str, required=True, help='Path to virus config YAML file')
    args = parser.parse_args(argv)

    config = load_config(args.config)
    af = config['alignment_filtering']
    dedup = config['deduplication']

    input_dir = os.path.abspath(args.input_dir)
    reference_fasta = os.path.abspath(args.reference_fasta)
    output_dir = os.path.abspath(args.output_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sam_files = glob.glob(os.path.join(input_dir, "*_aln.sam"))
    logging.info(f"Found {len(sam_files)} SAM files in {input_dir}: {sam_files}")
    if not sam_files:
        logging.error(f"No SAM files found in {input_dir}")
        raise FileNotFoundError(f"No SAM files found in {input_dir}")

    for sam_file in sam_files:
        logging.info(f"Processing: {sam_file}")
        sample_name = os.path.splitext(os.path.basename(sam_file))[0].replace('_aln', '')

        # Convert SAM to BAM with MAPQ and flag filtering
        bam_file = os.path.join(output_dir, f"{sample_name}.bam")
        sam_to_bam_command = (
            f"samtools view -S -b -T {reference_fasta} "
            f"-q {af['min_mapping_quality']} -F {af['exclude_flags']} "
            f"{sam_file} > {bam_file}"
        )
        run_command(sam_to_bam_command)
        validate_bam(bam_file)

        if dedup['enabled']:
            # Add mate information for duplicate marking
            fixmate_bam = os.path.join(output_dir, f"{sample_name}.fixmate.bam")
            fixmate_command = f"samtools fixmate -m {bam_file} {fixmate_bam}"
            run_command(fixmate_command)

            # Sort BAM file
            sorted_bam_file = os.path.join(output_dir, f"{sample_name}.sorted.bam")
            sort_command = f"samtools sort -T {os.path.join(output_dir, f'temp_{sample_name}')} -o {sorted_bam_file} {fixmate_bam}"
            run_command(sort_command)

            # Count pre-dedup reads
            pre_dedup_stdout, _ = run_command(f"samtools view -c {sorted_bam_file}")
            pre_dedup_count = int(pre_dedup_stdout.strip())

            # Mark and remove duplicates
            dedup_bam = os.path.join(output_dir, f"{sample_name}.dedup.bam")
            markdup_flags = "-r -s" if dedup['remove_duplicates'] else "-s"
            markdup_command = f"samtools markdup {markdup_flags} {sorted_bam_file} {dedup_bam}"
            run_command(markdup_command)
            os.rename(dedup_bam, sorted_bam_file)
            validate_bam(sorted_bam_file)

            # Count post-dedup reads
            post_dedup_stdout, _ = run_command(f"samtools view -c {sorted_bam_file}")
            post_dedup_count = int(post_dedup_stdout.strip())
            duplicates_removed = pre_dedup_count - post_dedup_count
            pct_duplicates = round(duplicates_removed / pre_dedup_count * 100, 2) if pre_dedup_count > 0 else 0
            logging.info(f"Deduplication: {pre_dedup_count} -> {post_dedup_count} reads "
                         f"({duplicates_removed} duplicates removed, {pct_duplicates}%)")

            # Write dedup stats file
            dedup_stats_file = os.path.join(output_dir, f"{sample_name}_dedup_stats.txt")
            with open(dedup_stats_file, 'w') as ds:
                ds.write(f"pre_dedup_reads\t{pre_dedup_count}\n")
                ds.write(f"post_dedup_reads\t{post_dedup_count}\n")
                ds.write(f"duplicates_removed\t{duplicates_removed}\n")
                ds.write(f"pct_duplicates\t{pct_duplicates}\n")

            # Clean up intermediate files
            for tmp in [bam_file, fixmate_bam]:
                if os.path.exists(tmp):
                    os.remove(tmp)
        else:
            # No dedup: just sort
            sorted_bam_file = os.path.join(output_dir, f"{sample_name}.sorted.bam")
            sort_command = f"samtools sort -T {os.path.join(output_dir, f'temp_{sample_name}')} -o {sorted_bam_file} {bam_file}"
            run_command(sort_command)
            validate_bam(sorted_bam_file)

            if os.path.exists(bam_file):
                os.remove(bam_file)

        # Index BAM file
        index_command = f"samtools index {sorted_bam_file}"
        run_command(index_command)

        logging.info(f"Sorted BAM file generated and indexed: {sorted_bam_file}")

if __name__ == '__main__':
    main()
