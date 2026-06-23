#!/usr/bin/env python3
import sys
import json
import argparse
import pandas as pd
import subprocess
import os
import logging
import shutil

from virus_pipeline.config import load_config

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def run_command(command):
    logging.info(f"Running command: {command}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    stdout_text = stdout.decode("utf-8") if stdout else ""
    stderr_text = stderr.decode("utf-8") if stderr else ""
    if process.returncode != 0:
        raise Exception(f"Command execution failed with return code {process.returncode}, stderr: {stderr_text}")
    return stdout_text, stderr_text


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


def run_host_removal(trimmed_r1, trimmed_r2, host_reference, threads, sample_folder, sample_name):
    """
    Map reads to combined host reference (human+mosquito) and extract reads that DO NOT map to host.
    Returns paths to host-depleted R1 and R2 (gzipped FASTQ).
    Requires bwa-mem2 and samtools in PATH.
    """
    host_sam = os.path.join(sample_folder, f"{sample_name}_host.sam")
    host_bam = os.path.join(sample_folder, f"{sample_name}_host.bam")
    host_depleted_bam = os.path.join(sample_folder, f"{sample_name}_host_depleted.bam")
    depleted_r1 = os.path.join(sample_folder, f"{sample_name}_host_depleted_1.fastq.gz")
    depleted_r2 = os.path.join(sample_folder, f"{sample_name}_host_depleted_2.fastq.gz")

    logging.info(f"Removing host reads (ref={host_reference}) for {sample_name}...")

    # 1) Map to host
    cmd_map = f"bwa-mem2 mem -t {threads} {host_reference} {trimmed_r1} {trimmed_r2} > {host_sam}"
    run_command(cmd_map)

    # 2) SAM -> BAM
    run_command(f"samtools view -bS {host_sam} -o {host_bam}")

    # 3) Extract pairs where BOTH mates are unmapped (-f 12)
    # Use -f 12 to require both mates unmapped; this is conservative and helps keep paired reads
    run_command(f"samtools view -b -f 12 -o {host_depleted_bam} {host_bam}")

    # 4) BAM -> paired FASTQ (compressed)
    # -n keeps read names, -1/-2 specify paired output, -0 /dev/null drop singletons, -s /dev/null drop singles
    run_command(f"samtools fastq -1 {depleted_r1} -2 {depleted_r2} -0 /dev/null -s /dev/null -n {host_depleted_bam}")

    # Cleanup intermediate files to save space (optional)
    for p in (host_sam, host_bam, host_depleted_bam):
        try:
            if os.path.exists(p):
                os.remove(p)
        except Exception:
            logging.debug(f"Could not remove intermediate file {p}")

    logging.info(f"Host-depleted FASTQs: {depleted_r1}, {depleted_r2}")
    return depleted_r1, depleted_r2


def run_denovo_and_classify(depleted_r1, depleted_r2, sample_folder, sample_name, threads, meta_cfg):
    """
    Run de-novo assembler (megahit by default) on host-depleted reads and classify contigs (kraken2 by default).
    Writes results into sample_folder/assembly_{sample_name}/ and classification files.
    """
    assembler = meta_cfg.get('assembler', 'megahit')
    kraken_db = meta_cfg.get('kraken2_db')
    outdir = os.path.join(sample_folder, f"assembly_{sample_name}")
    os.makedirs(outdir, exist_ok=True)

    contigs = None
    if assembler == 'megahit':
        logging.info(f"Running MEGAHIT for {sample_name}...")
        megahit_cmd = f"megahit -1 {depleted_r1} -2 {depleted_r2} -o {outdir} --min-contig-len {meta_cfg.get('min_contig_len',200)} -t {threads}"
        run_command(megahit_cmd)
        contigs = os.path.join(outdir, "final.contigs.fa")
    elif assembler in ('metaspades', 'metaSPAdes'):
        logging.info(f"Running metaSPAdes for {sample_name}...")
        spades_out = os.path.join(outdir, "spades_output")
        spades_cmd = f"spades.py --meta -1 {depleted_r1} -2 {depleted_r2} -o {spades_out} -t {threads}"
        run_command(spades_cmd)
        contigs = os.path.join(spades_out, "contigs.fasta")
    else:
        logging.warning(f"Assembler {assembler} not recognized; skipping assembly.")

    if contigs and os.path.isfile(contigs) and kraken_db:
        logging.info(f"Classifying contigs with kraken2 for {sample_name}...")
        kraken_report = os.path.join(outdir, "kraken2_contigs.report")
        kraken_out = os.path.join(outdir, "kraken2_contigs.out")
        kraken_cmd = f"kraken2 --db {kraken_db} --threads {threads} --report {kraken_report} --output {kraken_out} {contigs}"
        run_command(kraken_cmd)
        logging.info(f"Kraken2 report saved to {kraken_report}")
    else:
        if not kraken_db:
            logging.info("No kraken2 DB provided; skipping classification.")
        else:
            logging.info("No contigs produced; skipping classification.")

    return contigs


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('--samplesheet', type=str, required=True, help='Path to sample sheet with columns "sample_name", "read1", "read2".')
    parser.add_argument('--reference', type=str, required=True, help='Path to reference FASTA file')
    parser.add_argument('--config', type=str, required=True, help='Path to virus config YAML file')
    parser.add_argument('--threads', type=int, default=2, help='Number of parallel sample processing threads (default: 2)')
    parser.add_argument('--meta', action='store_true', help='Enable metagenomic mode (host removal + assembly + classification).')
    args = parser.parse_args(argv)

    config = load_config(args.config)
    fp = config.get('fastp', {})
    meta_cfg = config.get('meta', {})

    samples = samplesheet_verify(args.samplesheet)
    base_dir = os.path.abspath(os.path.dirname(args.samplesheet))
    sam_files_dir = os.path.join(base_dir, "sam_files")
    os.makedirs(sam_files_dir, exist_ok=True)

    # Reference indexing with bwa-mem2 (target virus reference)
    index_command = f"bwa-mem2 index {args.reference}"
    logging.info(f"Indexing reference {args.reference}...")
    stdout, stderr = run_command(index_command)
    logging.info(stdout)
    logging.info(stderr)

    # If meta mode: verify host reference is present and index it
    if args.meta:
        host_ref = meta_cfg.get('host_reference')
        if not host_ref:
            logging.error("Meta mode requested but meta.host_reference not set in config.")
            sys.exit(1)
        logging.info(f"Indexing host reference {host_ref} for host-removal step...")
        run_command(f"bwa-mem2 index {host_ref}")

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
            f"--qualified_quality_phred {fp.get('qualified_quality_phred',20)}",
            f"--length_required {fp.get('length_required',50)}",
            f"--cut_window_size {fp.get('cut_window_size',4)}",
            f"--cut_mean_quality {fp.get('cut_mean_quality',20)}",
            f"--overlap_len_require {fp.get('overlap_len_require',30)}",
            f"--json {fastp_json}",
            f"--html {fastp_html}",
            f"--thread {fp.get('threads',1)}",
        ]
        if fp.get('cut_front'):
            fastp_parts.append("--cut_front")
        if fp.get('cut_tail'):
            fastp_parts.append("--cut_tail")
        if fp.get('detect_adapter_for_pe'):
            fastp_parts.append("--detect_adapter_for_pe")
        if fp.get('correction'):
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
        fastqc_command = f"fastqc {trimmed_r1} {trimmed_r2} --outdir={sample_folder}"
        logging.info(f"Running FastQC for sample {sample_name}...")
        stdout, stderr = run_command(fastqc_command)
        logging.info(stdout)
        logging.info(stderr)

        # If meta mode enabled: remove host, assemble & classify, then map host-depleted reads.
        try:
            if args.meta:
                depleted_r1, depleted_r2 = run_host_removal(trimmed_r1, trimmed_r2, meta_cfg['host_reference'], args.threads, sample_folder, sample_name)

                # optional: run denovo assembly and classify contigs
                contigs = run_denovo_and_classify(depleted_r1, depleted_r2, sample_folder, sample_name, args.threads, meta_cfg)

                # If host depletion produced outputs, use them for downstream mapping
                if os.path.exists(depleted_r1) and os.path.exists(depleted_r2):
                    trimmed_r1, trimmed_r2 = depleted_r1, depleted_r2
                else:
                    logging.warning("Host depletion did not produce paired FASTQs; proceeding with original trimmed reads.")

            # Mapping with bwa-mem2 to virus reference
            sam_file = os.path.join(sample_folder, f"{sample_name}_aln.sam")
            bwa_command = f"bwa-mem2 mem {args.reference} {trimmed_r1} {trimmed_r2} > {sam_file}"
            logging.info(f"Running bwa-mem2 for sample {sample_name}...")
            stdout, stderr = run_command(bwa_command)
            logging.info(stdout)
            logging.info(stderr)

            # Move SAM files to sam_files_dir
            run_command(f"mv {sample_folder}/{sample_name}*.sam {sam_files_dir}")
            logging.info(f"SAM files moved to {sam_files_dir}")

        except Exception as e:
            logging.error(f"Error encountered while processing sample {sample_name}: {e}")


if __name__ == '__main__':
    main()
