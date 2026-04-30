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

def check_tools(annotation_mode='snpeff'):
    tools = ['bwa-mem2', 'samtools', 'fastp', 'fastqc', 'gatk', 'ivar', 'bcftools']
    if annotation_mode == 'snpeff':
        snpeff_found = False
        for tool in ['snpeff', 'snpEff']:
            if shutil.which(tool):
                snpeff_found = True
                break
        if not snpeff_found:
            logging.error("Required tool not found: snpeff/snpEff")
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
                        help='Comma-delimited list of sample names. Must match number of FASTQ pairs.')
    
    # UPDATE 1: Conditional Primer Trimming logic
    parser.add_argument('--primer_bed', type=str, default=None,
                        help='BED file with primer coordinates for ivar trim. '
                             'If not provided or file is missing, primer trimming is skipped.')
    
    # UPDATE 2: Java Version arguments
    parser.add_argument('--gatk_java', type=str, default='java',
                        help='Path to Java executable for GATK (e.g., /usr/bin/java8)')
    parser.add_argument('--snpeff_java', type=str, default='java',
                        help='Path to Java executable for SnpEff (e.g., /usr/bin/java11)')
    
    parser.add_argument('--annotation_mode', type=str, default='snpeff',
                        choices=['snpeff', 'config'],
                        help='Annotation mode: snpeff (default) or config (lightweight, no snpEff)')
    parser.add_argument('--gatk_memory', type=str, default=None,
                        help='Max Java heap size for GATK (e.g., 2g, 4g).')
    
    args = parser.parse_args()

    # Load virus-specific config
    check_file_exists(args.config, "Config file")
    config = load_config(args.config)
    database_name = config['database_name']

    # Check inputs and tools
    check_file_exists(args.input_dir, "Input directory")
    check_file_exists(args.reference_fasta, "Reference FASTA")
    
    # Verify primer BED only if path was provided
    if args.primer_bed:
        if not os.path.exists(args.primer_bed):
            logging.warning(f"Primer BED file '{args.primer_bed}' not found. Trimming will be skipped.")
            args.primer_bed = None

    if args.annotation_mode == 'snpeff':
        check_file_exists(args.genbank_file, "GenBank file")
    check_tools(args.annotation_mode)

    # Create and check output directories
    os.makedirs(args.output_dir, exist_ok=True)
    check_write_permission(args.output_dir)
    sam_files_dir = os.path.join(args.output_dir, "sam_files")
    os.makedirs(sam_files_dir, exist_ok=True)
    check_write_permission(sam_files_dir)

    sample_sheet = os.path.join(args.output_dir, "samplesheet.tsv")

    # Initialize provenance tracker
    tracker = ProvenanceTracker(args.output_dir)
    pipeline_args = vars(args).copy()
    pipeline_args['virus_name'] = config['virus_name']
    pipeline_args['ploidy'] = config['ploidy']
    pipeline_args['database_name'] = database_name
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
        check_file_exists(sample_sheet, "Sample sheet")
        tracker.record_step("Create sample sheet", "create_samplesheet",
                            {'input_dir': args.input_dir, 'output': sample_sheet})
        logging.info("Completed create_samplesheet")
    except Exception as e:
        tracker.record_step("Create sample sheet", "create_samplesheet",
                            {'input_dir': args.input_dir}, status="failed", notes=str(e))
        logging.error(f"Failed in create_samplesheet: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    # --- Step 2: Mapping ---
    try:
        logging.info("Starting map_reads")
        map_reads(['--samplesheet', sample_sheet, '--reference', args.reference_fasta,
                   '--config', args.config])
        tracker.record_step("Mapping (bwa-mem2)", "bwa-mem2", {'reference': args.reference_fasta})
        logging.info("Completed map_reads")
    except Exception as e:
        tracker.record_step("Read Mapping", "bwa-mem2", {}, status="failed", notes=str(e))
        logging.error(f"Failed in map_reads: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    # --- Step 3: SAM to BAM ---
    try:
        logging.info("Starting samtobamdenv")
        samtobamdenv(['--input_dir', sam_files_dir, '--reference_fasta', args.reference_fasta,
                      '--output_dir', args.output_dir, '--config', args.config])
        tracker.record_step("SAM to BAM", "samtools", {'deduplication': config['deduplication']['enabled']})
        logging.info("Completed samtobamdenv")
    except Exception as e:
        tracker.record_step("SAM to BAM", "samtools", {}, status="failed", notes=str(e))
        logging.error(f"Failed in samtobamdenv: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    # --- Step 4: SnpEff Database ---
    if args.annotation_mode == 'snpeff':
        try:
            logging.info(f"Starting create_snpeff_database using Java: {args.snpeff_java}")
            create_snpeff_database(['--genbank_file', args.genbank_file,
                                    '--reference_fasta', args.reference_fasta,
                                    '--output_dir', args.output_dir,
                                    '--database_name', database_name,
                                    '--java_path', args.snpeff_java]) # Updated
            tracker.record_step("Build SnpEff database", "snpEff", {'java': args.snpeff_java})
            logging.info("Completed create_snpeff_database")
        except Exception as e:
            tracker.record_step("Build SnpEff database", "snpEff", {}, status="failed", notes=str(e))
            logging.error(f"Failed in create_snpeff_database: {e}")
            tracker.write_json()
            tracker.write_report()
            sys.exit(1)

    # --- Step 5: Variant Calling & Consensus ---
    vcc_args = [
        '--input_dir', args.output_dir,
        '--reference_fasta', args.reference_fasta,
        '--output_dir', args.output_dir,
        '--database_name', database_name,
        '--config', args.config,
        '--annotation_mode', args.annotation_mode,
        '--gatk_java', args.gatk_java,    # Updated
        '--snpeff_java', args.snpeff_java # Updated
    ]
    
    # Conditional Primer Trimming
    if args.primer_bed:
        vcc_args.extend(['--primer_bed', args.primer_bed])
        logging.info(f"Primer trimming enabled with {args.primer_bed}")
    else:
        logging.info("No primer BED provided; primer trimming will be skipped.")

    if args.gatk_memory:
        vcc_args.extend(['--gatk_memory', args.gatk_memory])

    try:
        logging.info("Starting variant_calling_consensus")
        variant_calling_consensus(vcc_args)
        
        # Provenance Tracking logic
        pt_status = "completed" if args.primer_bed else "skipped"
        tracker.record_step("Primer trimming", "ivar trim", {'bed': args.primer_bed}, status=pt_status)
        tracker.record_step("Variant calling", "gatk HaplotypeCaller", {'java': args.gatk_java})
        tracker.record_step("Annotation", "snpEff", {'java': args.snpeff_java})

        logging.info("Completed variant_calling_consensus")
    except Exception as e:
        tracker.record_step("Variant calling + consensus", "multiple", {}, status="failed", notes=str(e))
        logging.error(f"Failed in variant_calling_consensus: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    # --- Step 6: Protein Extraction ---
    try:
        logging.info("Starting extract_proteins")
        proteins_dir = os.path.join(args.output_dir, 'proteins')
        extract_proteins(consensus_dir=args.output_dir, config_source=args.config, 
                         reference=args.reference_fasta, output_dir=proteins_dir)
        logging.info("Completed extract_proteins")
    except Exception as e:
        logging.warning(f"Protein extraction failed (non-fatal): {e}")

    # --- Step 7: Summarization ---
    try:
        logging.info("Starting summarize_result")
        summarize_result(['--input_dir', args.output_dir, '--output_dir', args.output_dir,
                          '--database_name', database_name])
        
        if args.annotation_mode == 'snpeff':
            summarize_snpEff(['--input_dir', args.output_dir, '--output_dir', args.output_dir,
                              '--config', args.config])
        else:
            from virus_pipeline import summarize_annotations
            summarize_annotations(['--input_dir', args.output_dir, '--output_dir', args.output_dir])
            
        logging.info("Summarization complete.")
    except Exception as e:
        logging.error(f"Summarization failed: {e}")

    tracker.write_json()
    tracker.write_report()
    logging.info("Pipeline completed successfully!")

if __name__ == "__main__":
    main()
