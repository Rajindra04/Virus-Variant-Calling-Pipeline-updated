#!/usr/bin/env python

import argparse
import os
import sys
import logging
import shutil
import subprocess

from virus_pipeline import (
    create_samplesheet,
    map_reads,
    samtobamdenv,
    create_snpeff_database,
    variant_calling_consensus,
    summarize_result,
    summarize_snpEff,
)
from virus_pipeline.config import load_config
from virus_pipeline.extract_proteins import run_extraction as extract_proteins
from virus_pipeline.provenance import ProvenanceTracker

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def check_file_exists(file_path, description):
    if not file_path:
        return
    if not os.path.exists(file_path):
        logging.error(f"{description} not found: {file_path}")
        sys.exit(1)

def check_write_permission(directory):
    try:
        test_file = os.path.join(directory, '.test_write')
        with open(test_file, 'w') as f:
            f.write('test')
        os.remove(test_file)
    except (PermissionError, OSError) as e:
        logging.error(f"No write permission for {directory}: {e}")
        sys.exit(1)

def detect_java_versions():
    """
    Detect Java 8 for GATK and Java 21/17/11 for SnpEff.
    Returns (java8_path, snpeff_java_path)
    """
    candidates = [
        "/usr/lib/jvm/java-21-openjdk-amd64/bin/java",
        "/usr/lib/jvm/java-17-openjdk-amd64/bin/java",
        "/usr/lib/jvm/java-11-openjdk-amd64/bin/java",
        "/usr/lib/jvm/java-8-openjdk-amd64/bin/java",
        "java",
    ]

    java8 = None
    java_modern = None

    for path in candidates:
        exe = shutil.which(path)
        if not exe:
            continue

        try:
            out = subprocess.check_output([exe, "-version"], stderr=subprocess.STDOUT).decode()
        except Exception:
            continue

        # Detect Java 8 (specifically for GATK)
        if ("1.8" in out or " 8." in out) and not java8:
            java8 = exe
        
        # Detect Modern Java (Prioritizing 21 for SnpEff)
        if "21." in out:
            java_modern = exe  # Highest priority
        elif ("17." in out or "11." in out) and not java_modern:
            java_modern = exe

    return java8, java_modern

def check_tools(annotation_mode='snpeff', gatk_java='auto', snpeff_java='auto'):
    """Check if required binaries are available."""
    tools = ['bwa-mem2', 'samtools', 'fastp', 'fastqc', 'ivar', 'bcftools']

    for j_path, label in [(gatk_java, "GATK Java"), (snpeff_java, "SnpEff Java")]:
        if j_path not in ('auto', 'java') and not os.path.exists(j_path):
            logging.error(f"Custom {label} path not found: {j_path}")
            sys.exit(1)

    if annotation_mode == 'snpeff':
        snpeff_found = False
        for tool in ['snpeff', 'snpEff']:
            if shutil.which(tool):
                snpeff_found = True
                break
        if not snpeff_found and snpeff_java in ('auto', 'java'):
            logging.error("Required tool not found: snpeff/snpEff and no custom --snpeff_java provided.")
            sys.exit(1)

    for tool in tools:
        if not shutil.which(tool):
            logging.error(f"Required tool not found: {tool}")
            sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Automate variant calling pipeline.')
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--reference_fasta', type=str, required=True)
    parser.add_argument('--genbank_file', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--config', type=str, required=True, help='Path to config YAML.')
    parser.add_argument('--sample_names', type=str, default=None)
    parser.add_argument('--primer_bed', type=str, default=None)
    
    # New parameter to configure parallel pipeline execution threads
    parser.add_argument('--threads', type=int, default=2, help='Number of parallel sample processing threads (default: 2)')
    
    # These flags are now properly recognized by the parser
    parser.add_argument('--gatk_java', type=str, default='auto', help='Java path for GATK (Java 8)')
    parser.add_argument('--snpeff_java', type=str, default='auto', help='Java path for SnpEff (Java 21/17/11)')
    
    parser.add_argument('--annotation_mode', type=str, default='snpeff', choices=['snpeff', 'config'])
    parser.add_argument('--gatk_memory', type=str, default=None)

    args = parser.parse_args()

    # Initial Validations
    check_file_exists(args.config, "Config file")
    config = load_config(args.config)
    database_name = config['database_name']

    check_file_exists(args.input_dir, "Input directory")
    check_file_exists(args.reference_fasta, "Reference FASTA")

    if args.annotation_mode == 'snpeff':
        check_file_exists(args.genbank_file, "GenBank file")

    # Handle Java Version Auto-Detection
    java8_path, java_modern_path = detect_java_versions()

    if args.gatk_java == 'auto':
        if java8_path:
            args.gatk_java = java8_path
            logging.info(f"Auto-detected Java 8 for GATK: {java8_path}")
        else:
            args.gatk_java = 'java' # Fallback
            logging.warning("Java 8 not found. Falling back to default 'java'.")

    if args.snpeff_java == 'auto':
        if java_modern_path:
            args.snpeff_java = java_modern_path
            logging.info(f"Auto-detected Modern Java for SnpEff: {java_modern_path}")
        else:
            logging.error("No suitable Java (11, 17, or 21) found for SnpEff.")
            sys.exit(1)

    check_tools(args.annotation_mode, args.gatk_java, args.snpeff_java)

    # Setup Directories
    os.makedirs(args.output_dir, exist_ok=True)
    sam_files_dir = os.path.join(args.output_dir, "sam_files")
    os.makedirs(sam_files_dir, exist_ok=True)
    sample_sheet = os.path.join(args.output_dir, "samplesheet.tsv")

    # Provenance Tracking
    tracker = ProvenanceTracker(args.output_dir)
    tracker.set_pipeline_args(vars(args))
    tracker.set_config(config)
    tracker.detect_all_tool_versions()

    # --- Execution Steps ---
    
    try:
        logging.info("Step 1: Creating Samplesheet")
        ss_args = [args.input_dir, sample_sheet]
        if args.sample_names: ss_args.extend(['--sample_names', args.sample_names])
        create_samplesheet(ss_args)

        logging.info("Step 2: Mapping Reads")
        map_reads([
            '--samplesheet', sample_sheet, 
            '--reference', args.reference_fasta, 
            '--config', args.config,
            '--threads', str(args.threads)  # Passed threads value downstream
        ])

        logging.info("Step 3: SAM to BAM conversion")
        samtobamdenv([
            '--input_dir', sam_files_dir, 
            '--reference_fasta', args.reference_fasta, 
            '--output_dir', args.output_dir, 
            '--config', args.config,
            '--threads', str(args.threads)  # Passed threads value downstream
        ])

        if args.annotation_mode == 'snpeff':
            logging.info(f"Step 4: Building SnpEff database with {args.snpeff_java}")
            create_snpeff_database([
                '--genbank_file', args.genbank_file,
                '--reference_fasta', args.reference_fasta,
                '--output_dir', args.output_dir,
                '--database_name', database_name,
                '--snpeff_java', args.snpeff_java
            ])

        logging.info("Step 5: Variant Calling & Consensus")
        vcc_args = [
            '--input_dir', args.output_dir,
            '--reference_fasta', args.reference_fasta,
            '--output_dir', args.output_dir,
            '--database_name', database_name,
            '--config', args.config,
            '--annotation_mode', args.annotation_mode,
            '--gatk_java', args.gatk_java,
            '--snpeff_java', args.snpeff_java,
            '--threads', str(args.threads),  # Passed threads value downstream
        ]
        if args.primer_bed: vcc_args.extend(['--primer_bed', args.primer_bed])
        if args.gatk_memory: vcc_args.extend(['--gatk_memory', args.gatk_memory])
        variant_calling_consensus(vcc_args)

        logging.info("Finalizing: Extracting Proteins and Summarizing")
        extract_proteins(consensus_dir=args.output_dir, config_source=args.config, reference=args.reference_fasta, output_dir=os.path.join(args.output_dir, 'proteins'))
        summarize_result(['--input_dir', args.output_dir, '--output_dir', args.output_dir, '--database_name', database_name])
        
        if args.annotation_mode == 'snpeff':
            summarize_snpEff(['--input_dir', args.output_dir, '--output_dir', args.output_dir, '--config', args.config])

    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        sys.exit(1)

    tracker.write_json()
    tracker.write_report()
    logging.info("Pipeline completed successfully!")

if __name__ == "__main__":
    main()
