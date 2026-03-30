import sys
import json
import argparse
import pandas as pd
import subprocess
import os
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

def samplesheet_verify(samplesheet):
    try:
        samples = pd.read_csv(samplesheet, sep='\t')
        if not all(col in samples.columns for col in ['sample_name', 'read1', 'read2']):
            raise ValueError("Sample sheet must have columns 'sample_name', 'read1', 'read2'.")
        if samples.empty:
            raise ValueError("Sample sheet contains no samples.")
    except Exception as e:
        raise ValueError(f"Error reading sample sheet: {e}")
    return samples


def check_qc_gate(fastp_json, sample_name, base_dir):
    """Parse fastp JSON and enforce QC thresholds.

    Returns:
        str: "PASS", "WARN", or "FAIL"
    """
    qc_fail_log = os.path.join(base_dir, "qc_fail_log.txt")
    qc_warn_log = os.path.join(base_dir, "qc_warn_log.txt")

    with open(fastp_json, "r") as f:
        data = json.load(f)

    after = data["summary"]["after_filtering"]
    q30_rate = after["q30_rate"]
    total_reads = after["total_reads"]

    # Duplication rate is at top level in fastp JSON
    duplication_rate = data.get("duplication", {}).get("rate", 0.0)

    status = "PASS"
    reasons = []

    # FAIL conditions
    if q30_rate < 0.70:
        reasons.append(f"Q30 rate {q30_rate:.3f} < 0.70")
        status = "FAIL"
    if total_reads < 10000:
        reasons.append(f"Total reads {total_reads} < 10000")
        status = "FAIL"

    # WARN conditions (only if not already FAIL)
    if status != "FAIL":
        if q30_rate < 0.80:
            reasons.append(f"Q30 rate {q30_rate:.3f} < 0.80")
            status = "WARN"
        if duplication_rate > 0.80:
            reasons.append(f"Duplication rate {duplication_rate:.3f} > 0.80")
            status = "WARN"

    # Log results
    if status == "FAIL":
        with open(qc_fail_log, "a") as f:
            f.write(f"{sample_name}\tFAIL\t{'; '.join(reasons)}\t"
                    f"q30={q30_rate:.3f}\treads={total_reads}\tdup={duplication_rate:.3f}\n")
        logging.error(f"QC FAIL for {sample_name}: {'; '.join(reasons)}")
    elif status == "WARN":
        with open(qc_warn_log, "a") as f:
            f.write(f"{sample_name}\tWARN\t{'; '.join(reasons)}\t"
                    f"q30={q30_rate:.3f}\treads={total_reads}\tdup={duplication_rate:.3f}\n")
        logging.warning(f"QC WARN for {sample_name}: {'; '.join(reasons)}")
    else:
        logging.info(f"QC PASS for {sample_name}: q30={q30_rate:.3f}, reads={total_reads}, dup={duplication_rate:.3f}")

    return status


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('--samplesheet', type=str, required=True, help='Path to sample sheet with columns "sample_name", "read1", "read2".')
    parser.add_argument('--reference', type=str, required=True, help='Path to reference FASTA file')
    parser.add_argument('--config', type=str, required=True, help='Path to virus config YAML file')
    args = parser.parse_args(argv)

    config = load_config(args.config)
    fp = config['fastp']

    samples = samplesheet_verify(args.samplesheet)
    base_dir = os.path.abspath(os.path.dirname(args.samplesheet))
    sam_files_dir = os.path.join(base_dir, "sam_files")
    os.makedirs(sam_files_dir, exist_ok=True)

    # Reference indexing with bwa-mem2
    index_command = f"bwa-mem2 index {args.reference}"
    logging.info(f"Indexing reference {args.reference}...")
    stdout, stderr = run_command(index_command)
    logging.info(stdout)
    logging.info(stderr)

    for i, row in samples.iterrows():
        sample_name = row['sample_name']
        read1 = row['read1']
        read2 = row['read2']

        sample_folder = os.path.join(base_dir, f"{sample_name}_output")
        os.makedirs(sample_folder, exist_ok=True)

        # fastp with parameters from config
        trimmed_r1 = os.path.join(sample_folder, f'{sample_name}_trimmed_R1.fastq.gz')
        trimmed_r2 = os.path.join(sample_folder, f'{sample_name}_trimmed_R2.fastq.gz')
        fastp_json = os.path.join(sample_folder, f'{sample_name}_fastp.json')
        fastp_html = os.path.join(sample_folder, f'{sample_name}_fastp.html')

        fastp_parts = [
            f"fastp -i {read1} -I {read2}",
            f"-o {trimmed_r1} -O {trimmed_r2}",
            f"--qualified_quality_phred {fp['qualified_quality_phred']}",
            f"--length_required {fp['length_required']}",
            f"--cut_window_size {fp['cut_window_size']}",
            f"--cut_mean_quality {fp['cut_mean_quality']}",
            f"--overlap_len_require {fp['overlap_len_require']}",
            f"--json {fastp_json}",
            f"--html {fastp_html}",
            f"--thread {fp['threads']}",
        ]
        if fp['cut_front']:
            fastp_parts.append("--cut_front")
        if fp['cut_tail']:
            fastp_parts.append("--cut_tail")
        if fp['detect_adapter_for_pe']:
            fastp_parts.append("--detect_adapter_for_pe")
        if fp['correction']:
            fastp_parts.append("--correction")

        fastp_command = " ".join(fastp_parts)
        logging.info(f"Running fastp for sample {sample_name}...")
        stdout, stderr = run_command(fastp_command)
        logging.info(stdout)
        logging.info(stderr)

        # QC Gate: check fastp results before proceeding
        qc_status = check_qc_gate(fastp_json, sample_name, base_dir)
        if qc_status == "FAIL":
            logging.error(f"Skipping sample {sample_name} due to QC failure")
            continue

        # FastQC
        fastqc_command = (
            f"fastqc {trimmed_r1} {trimmed_r2} --outdir={sample_folder}"
        )
        logging.info(f"Running FastQC for sample {sample_name}...")
        stdout, stderr = run_command(fastqc_command)
        logging.info(stdout)
        logging.info(stderr)

        # Mapping with bwa-mem2
        sam_file = os.path.join(sample_folder, f"{sample_name}_aln.sam")
        bwa_command = f"bwa-mem2 mem {args.reference} {trimmed_r1} {trimmed_r2} > {sam_file}"
        logging.info(f"Running bwa-mem2 for sample {sample_name}...")
        stdout, stderr = run_command(bwa_command)
        logging.info(stdout)
        logging.info(stderr)

        # Move SAM files to sam_files_dir
        run_command(f"mv {sample_folder}/{sample_name}*.sam {sam_files_dir}")
        logging.info(f"SAM files moved to {sam_files_dir}")

if __name__ == '__main__':
    main()
