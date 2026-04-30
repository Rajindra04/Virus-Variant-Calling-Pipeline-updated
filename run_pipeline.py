#!/usr/bin/env python

import argparse
import os
import sys
import logging
import shutil

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

def check_tools(annotation_mode='snpeff', gatk_java='java', snpeff_java='java'):
    """
    Check if required binaries are available. 
    If a custom path is provided for Java, we verify that path exists.
    """
    tools = ['bwa-mem2', 'samtools', 'fastp', 'fastqc', 'ivar', 'bcftools']
    
    # Check custom Java paths if they aren't the default 'java'
    for j_path, label in [(gatk_java, "GATK Java"), (snpeff_java, "SnpEff Java")]:
        if j_path != 'java' and not os.path.exists(j_path):
            logging.error(f"Custom {label} path not found: {j_path}")
            sys.exit(1)

    if annotation_mode == 'snpeff':
        snpeff_found = False
        for tool in ['snpeff', 'snpEff']:
            if shutil.which(tool):
                snpeff_found = True
                break
        if not snpeff_found and snpeff_java == 'java':
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
    parser.add_argument('--config', type=str, required=True,
                        help='Path to virus config YAML file (e.g., configs/denv1.yaml)')
    parser.add_argument('--sample_names', type=str, default=None,
                        help='Comma-delimited list of sample names.')
    
    # Feature 1: Optional Primer Trimming
    parser.add_argument('--primer_bed', type=str, default=None,
                        help='BED file for ivar trim. If missing/not provided, trimming is skipped.')
    
    # Feature 2: Multi-Java Version Support
    parser.add_argument('--gatk_java', type=str, default='java',
                        help='Path to Java binary for GATK (e.g., /usr/bin/java8)')
    parser.add_argument('--snpeff_java', type=str, default='java',
                        help='Path to Java binary for SnpEff (e.g., /usr/bin/java11)')
    
    parser.add_argument('--annotation_mode', type=str, default='snpeff',
                        choices=['snpeff', 'config'])
    parser.add_argument('--gatk_memory', type=str, default=None)
    
    args = parser.parse_args()

    # Initial Validations
    check_file_exists(args.config, "Config file")
    config = load_config(args.config)
    database_name = config['database_name']

    check_file_exists(args.input_dir, "Input directory")
    check_file_exists(args.reference_fasta, "Reference FASTA")
    
    # Validate primer BED: If provided but doesn't exist, set to None to trigger skip logic
    if args.primer_bed:
        if not os.path.exists(args.primer_bed):
            logging.warning(f"Primer BED file '{args.primer_bed}' not found. Skipping trimming.")
            args.primer_bed = None

    if args.annotation_mode == 'snpeff':
        check_file_exists(args.genbank_file, "GenBank file")
    
    check_tools(args.annotation_mode, args.gatk_java, args.snpeff_java)

    # Setup Directories
    os.makedirs(args.output_dir, exist_ok=True)
    sam_files_dir = os.path.join(args.output_dir, "sam_files")
    os.makedirs(sam_files_dir, exist_ok=True)
    sample_sheet = os.path.join(args.output_dir, "samplesheet.tsv")

    # Provenance Tracking
    tracker = ProvenanceTracker(args.output_dir)
    pipeline_args = vars(args).copy()
    tracker.set_pipeline_args(pipeline_args)
    tracker.set_config(config)
    tracker.detect_all_tool_versions()

    # --- Step 1: Samplesheet ---
    try:
        logging.info("Starting create_samplesheet")
        samplesheet_args = [args.input_dir, sample_sheet]
        if args.sample_names:
            samplesheet_args.extend(['--sample_names', args.sample_names])
        create_samplesheet(samplesheet_args)
    except Exception as e:
        logging.error(f"Failed in create_samplesheet: {e}")
        sys.exit(1)

    # --- Step 2: Mapping ---
    try:
        logging.info("Starting map_reads")
        map_reads(['--samplesheet', sample_sheet, '--reference', args.reference_fasta, '--config', args.config])
    except Exception as e:
        logging.error(f"Failed in map_reads: {e}")
        sys.exit(1)

    # --- Step 3: SAM to BAM ---
    try:
        logging.info("Starting samtobamdenv")
        samtobamdenv(['--input_dir', sam_files_dir, '--reference_fasta', args.reference_fasta,
                      '--output_dir', args.output_dir, '--config', args.config])
    except Exception as e:
        logging.error(f"Failed in samtobamdenv: {e}")
        sys.exit(1)

    # --- Step 4: SnpEff Database ---
    if args.annotation_mode == 'snpeff':
        try:
            logging.info(f"Building SnpEff database using: {args.snpeff_java}")
            create_snpeff_database(['--genbank_file', args.genbank_file,
                                    '--reference_fasta', args.reference_fasta,
                                    '--output_dir', args.output_dir,
                                    '--database_name', database_name,
                                    '--java_path', args.snpeff_java])
        except Exception as e:
            logging.error(f"Failed in create_snpeff_database: {e}")
            sys.exit(1)

    # --- Step 5: Variant Calling & Consensus (Trimming happens here) ---
    vcc_args = [
        '--input_dir', args.output_dir,
        '--reference_fasta', args.reference_fasta,
        '--output_dir', args.output_dir,
        '--database_name', database_name,
        '--config', args.config,
        '--annotation_mode', args.annotation_mode,
        '--gatk_java', args.gatk_java,    
        '--snpeff_java', args.snpeff_java 
    ]
    
    if args.primer_bed:
        vcc_args.extend(['--primer_bed', args.primer_bed])
        logging.info(f"Primer trimming active: {args.primer_bed}")
    else:
        logging.info("Primer trimming SKIPPED (no valid BED provided).")

    if args.gatk_memory:
        vcc_args.extend(['--gatk_memory', args.gatk_memory])

    try:
        logging.info("Starting variant_calling_consensus")
        variant_calling_consensus(vcc_args)
        
        # Log status to provenance
        pt_status = "completed" if args.primer_bed else "skipped"
        tracker.record_step("Primer trimming", "ivar trim", {'bed': args.primer_bed}, status=pt_status)
    except Exception as e:
        logging.error(f"Failed in variant_calling_consensus: {e}")
        sys.exit(1)

    # --- Final Steps (Proteins & Summaries) ---
    try:
        extract_proteins(consensus_dir=args.output_dir, config_source=args.config, 
                         reference=args.reference_fasta, output_dir=os.path.join(args.output_dir, 'proteins'))
        summarize_result(['--input_dir', args.output_dir, '--output_dir', args.output_dir, '--database_name', database_name])
        
        if args.annotation_mode == 'snpeff':
            summarize_snpEff(['--input_dir', args.output_dir, '--output_dir', args.output_dir, '--config', args.config])
    except Exception as e:
        logging.warning(f"Final summary steps encountered issues: {e}")

    tracker.write_json()
    tracker.write_report()
    logging.info("Pipeline completed successfully!")

if __name__ == "__main__":
    main()
