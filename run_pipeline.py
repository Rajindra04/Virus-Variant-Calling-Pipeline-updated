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
        # Check for both snpeff and snpEff
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
    parser.add_argument('--primer_bed', type=str, default=None,
                        help='BED file with primer coordinates for ivar trim. '
                             'If not provided, primer trimming is skipped.')
    parser.add_argument('--annotation_mode', type=str, default='snpeff',
                        choices=['snpeff', 'config'],
                        help='Annotation mode: snpeff (default) or config (lightweight, no snpEff)')
    parser.add_argument('--gatk_memory', type=str, default=None,
                        help='Max Java heap size for GATK (e.g., 2g, 4g). '
                             'Overrides the gatk_memory value in the config YAML. '
                             'Default from config is 4g. Reduce on low-memory systems.')
    args = parser.parse_args()

    # Load virus-specific config
    check_file_exists(args.config, "Config file")
    config = load_config(args.config)

    # Config provides database_name
    database_name = config['database_name']

    # Check inputs and tools
    check_file_exists(args.input_dir, "Input directory")
    check_file_exists(args.reference_fasta, "Reference FASTA")
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

    # Run pipeline steps
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
                           {'input_dir': args.input_dir}, status="failed",
                           notes=str(e))
        logging.error(f"Failed in create_samplesheet: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    try:
        logging.info("Starting map_reads")
        map_reads(['--samplesheet', sample_sheet, '--reference', args.reference_fasta,
                   '--config', args.config])
        tracker.record_step("Read trimming (fastp) + QC (FastQC) + Mapping (bwa-mem2)",
                           "fastp, fastqc, bwa-mem2", {
            'fastp_qualified_quality_phred': config['fastp']['qualified_quality_phred'],
            'fastp_length_required': config['fastp']['length_required'],
            'fastp_cut_front': config['fastp']['cut_front'],
            'fastp_cut_tail': config['fastp']['cut_tail'],
            'fastp_cut_window_size': config['fastp']['cut_window_size'],
            'fastp_cut_mean_quality': config['fastp']['cut_mean_quality'],
            'fastp_detect_adapter_for_pe': config['fastp']['detect_adapter_for_pe'],
            'fastp_correction': config['fastp']['correction'],
            'fastp_overlap_len_require': config['fastp']['overlap_len_require'],
            'bwa_mem2_reference': args.reference_fasta,
        })
        logging.info("Completed map_reads")
    except Exception as e:
        tracker.record_step("Read trimming + QC + Mapping", "fastp, fastqc, bwa-mem2",
                           {}, status="failed", notes=str(e))
        logging.error(f"Failed in map_reads: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    try:
        logging.info("Starting samtobamdenv")
        samtobamdenv(['--input_dir', sam_files_dir, '--reference_fasta', args.reference_fasta,
                      '--output_dir', args.output_dir, '--config', args.config])
        tracker.record_step("SAM to BAM conversion + filtering + dedup", "samtools", {
            'samtools_view_min_mapq': config['alignment_filtering']['min_mapping_quality'],
            'samtools_view_exclude_flags': config['alignment_filtering']['exclude_flags'],
            'samtools_fixmate': True,
            'samtools_markdup_remove': config['deduplication']['remove_duplicates'],
            'deduplication_enabled': config['deduplication']['enabled'],
        })
        logging.info("Completed samtobamdenv")
    except Exception as e:
        tracker.record_step("SAM to BAM conversion", "samtools",
                           {}, status="failed", notes=str(e))
        logging.error(f"Failed in samtobamdenv: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    if args.annotation_mode == 'snpeff':
        try:
            logging.info("Starting create_snpeff_database")
            create_snpeff_database(['--genbank_file', args.genbank_file,
                                    '--reference_fasta', args.reference_fasta,
                                    '--output_dir', args.output_dir,
                                    '--database_name', database_name])
            tracker.record_step("Build SnpEff database", "snpEff", {
                'genbank_file': args.genbank_file,
                'database_name': database_name,
            })
            logging.info("Completed create_snpeff_database")
        except Exception as e:
            tracker.record_step("Build SnpEff database", "snpEff",
                               {}, status="failed", notes=str(e))
            logging.error(f"Failed in create_snpeff_database: {e}")
            tracker.write_json()
            tracker.write_report()
            sys.exit(1)
    else:
        logging.info("Skipping snpEff database creation (annotation_mode=config)")
        tracker.record_step("Build SnpEff database", "snpEff", {},
                           status="skipped", notes="Using config-based annotation")

    # Build variant_calling_consensus args
    vcc_args = [
        '--input_dir', args.output_dir,
        '--reference_fasta', args.reference_fasta,
        '--output_dir', args.output_dir,
        '--database_name', database_name,
        '--config', args.config,
        '--annotation_mode', args.annotation_mode,
    ]
    if args.primer_bed:
        vcc_args.extend(['--primer_bed', args.primer_bed])
    if args.gatk_memory:
        vcc_args.extend(['--gatk_memory', args.gatk_memory])

    try:
        logging.info("Starting variant_calling_consensus")
        variant_calling_consensus(vcc_args)

        # Record primer trimming status (after step has actually run)
        pt = config['primer_trimming']
        if args.primer_bed:
            tracker.record_step("Primer trimming", "ivar trim", {
                'primer_bed': args.primer_bed,
                'min_quality': pt['min_quality'],
                'min_length': pt['min_length'],
                'sliding_window': pt['sliding_window'],
                'include_reads_no_primer': pt['include_reads_no_primer'],
            })
        else:
            tracker.record_step("Primer trimming", "ivar trim", {
                'primer_bed': 'NOT PROVIDED',
            }, status="skipped",
               notes="No --primer_bed argument provided. Primer sequences were NOT "
                     "removed from aligned reads. Variants at primer binding sites "
                     "may be masked by primer sequence. Provide a BED file with "
                     "--primer_bed to enable this step.")

        cov = config['coverage']
        tracker.record_step("Coverage calculation", "samtools depth", {
            'min_mapping_quality': cov['min_mapping_quality'],
            'min_base_quality': cov['min_base_quality'],
        })

        cons = config['consensus']
        tracker.record_step("Consensus generation", "samtools mpileup + ivar consensus", {
            'mpileup_max_depth': cons['mpileup_max_depth'],
            'mpileup_min_base_quality': cons['mpileup_min_base_quality'],
            'mpileup_min_mapping_quality': cons['mpileup_min_mapping_quality'],
            'ivar_min_quality': cons['ivar_min_quality'],
            'ivar_min_frequency_threshold': cons['ivar_min_frequency'],
            'ivar_min_depth': cons['ivar_min_depth'],
            'ivar_ambiguous_char': cons['ivar_ambiguous_char'],
        })

        vc = config['variant_calling']
        tracker.record_step("Variant calling", "gatk HaplotypeCaller", {
            'ploidy': vc['ploidy'],
            'standard_min_confidence': vc['standard_min_confidence'],
            'min_base_quality_score': vc['min_base_quality_score'],
        })

        vf = config['vcf_filtering']
        tracker.record_step("VCF filtering", "gatk VariantFiltration + SelectVariants", {
            'filters': vf['filters'],
            'select_pass_only': vf['select_pass_only'],
        })

        tracker.record_step("Variant annotation", "snpEff + SnpSift", {
            'database': database_name,
            'snpeff_memory': config['annotation']['snpeff_memory'],
            'input': 'filtered PASS variants',
        })
        logging.info("Completed variant_calling_consensus")
    except Exception as e:
        tracker.record_step("Variant calling + consensus", "multiple",
                           {}, status="failed", notes=str(e))
        logging.error(f"Failed in variant_calling_consensus: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    try:
        logging.info("Starting extract_proteins")
        proteins_dir = os.path.join(args.output_dir, 'proteins')
        extract_proteins(
            consensus_dir=args.output_dir,
            config_source=args.config,
            reference=args.reference_fasta,
            output_dir=proteins_dir
        )
        tracker.record_step("Protein extraction", "extract_proteins", {
            'output_dir': proteins_dir,
        })
        logging.info("Completed extract_proteins")
    except Exception as e:
        tracker.record_step("Protein extraction", "extract_proteins",
                           {}, status="failed", notes=str(e))
        logging.warning(f"Protein extraction failed (non-fatal): {e}")

    try:
        logging.info("Starting summarize_result")
        summarize_result(['--input_dir', args.output_dir, '--output_dir', args.output_dir,
                          '--database_name', database_name])
        tracker.record_step("Coverage summary", "summarize_result", {
            'metrics': 'mean, median, min, max, std, %genome at 1x/10x/30x/100x',
        })
        logging.info("Completed summarize_result")
    except Exception as e:
        tracker.record_step("Coverage summary", "summarize_result",
                           {}, status="failed", notes=str(e))
        logging.error(f"Failed in summarize_result: {e}")
        tracker.write_json()
        tracker.write_report()
        sys.exit(1)

    if args.annotation_mode == 'snpeff':
        try:
            logging.info("Starting summarize_snpEff")
            summarize_snpEff(['--input_dir', args.output_dir, '--output_dir', args.output_dir,
                              '--config', args.config])
            tracker.record_step("SnpEff summary", "summarize_snpEff", {})
            logging.info("Completed summarize_snpEff")
        except Exception as e:
            tracker.record_step("SnpEff summary", "summarize_snpEff",
                               {}, status="failed", notes=str(e))
            logging.error(f"Failed in summarize_snpEff: {e}")
            tracker.write_json()
            tracker.write_report()
            sys.exit(1)
    else:
        try:
            logging.info("Starting summarize_annotations (config mode)")
            from virus_pipeline import summarize_annotations
            summarize_annotations(['--input_dir', args.output_dir,
                                   '--output_dir', args.output_dir])
            tracker.record_step("Annotation summary", "summarize_annotations", {})
            logging.info("Completed summarize_annotations")
        except Exception as e:
            tracker.record_step("Annotation summary", "summarize_annotations",
                               {}, status="failed", notes=str(e))
            logging.error(f"Failed in summarize_annotations: {e}")
            tracker.write_json()
            tracker.write_report()
            sys.exit(1)

    # Write provenance reports
    tracker.write_json()
    tracker.write_report()

    logging.info("Pipeline completed successfully!")
    logging.info(f"Provenance report: {os.path.join(args.output_dir, 'provenance_report.txt')}")
    logging.info(f"Provenance JSON:   {os.path.join(args.output_dir, 'provenance.json')}")

if __name__ == "__main__":
    main()
