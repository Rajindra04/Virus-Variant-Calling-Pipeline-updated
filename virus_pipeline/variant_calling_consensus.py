import sys
import argparse
import glob
import os
import subprocess
import matplotlib.pyplot as plt
import re
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

def validate_fasta(fasta_file):
    valid_bases = set("ACGTNacgtnRYKMSWBDHVrykmswbvh")
    with open(fasta_file, "r") as f:
        sequence = ""
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    invalid_bases = set(sequence) - valid_bases
    if invalid_bases:
        raise ValueError(f"Non-IUPAC bases found in {fasta_file}: {invalid_bases}")
    logging.info(f"FASTA file {fasta_file} validated: contains only IUPAC bases.")

def validate_bam(bam_file):
    """Validate BAM file using samtools quickcheck and check alignment count."""
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file {bam_file} does not exist")
    validate_command = f"samtools quickcheck {bam_file}"
    stdout, stderr = run_command(validate_command)
    if stderr:
        raise Exception(f"BAM validation failed for {bam_file}: {stderr}")
    # Check alignment count
    count_command = f"samtools view -c {bam_file}"
    stdout, stderr = run_command(count_command)
    alignment_count = int(stdout.strip())
    if alignment_count == 0:
        raise ValueError(f"BAM file {bam_file} contains no alignments")
    logging.info(f"BAM file {bam_file} validated successfully with {alignment_count} alignments")

def plot_coverage(coverage_file, output_dir, sample_name=None):
    positions = []
    depths = []
    with open(coverage_file, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            position = int(columns[1])
            depth = int(columns[2])
            positions.append(position)
            depths.append(depth)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.fill_between(positions, depths, alpha=0.4, color='steelblue')
    ax.plot(positions, depths, linewidth=0.5, color='steelblue')

    # 20X threshold line (ivar consensus min depth)
    ax.axhline(y=20, color='red', linestyle='--', linewidth=1, label='20X (consensus threshold)')

    ax.set_yscale('log')
    ax.set_ylim(bottom=0.5)  # So 0-depth positions show near the bottom
    ax.set_xlabel('Genome Position (bp)')
    ax.set_ylabel('Read Depth (log scale)')
    title = f'Coverage: {sample_name}' if sample_name else 'Read Coverage'
    ax.set_title(title)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    coverage_plot_file = os.path.join(output_dir,
        os.path.splitext(os.path.basename(coverage_file))[0] + '.png')
    fig.savefig(coverage_plot_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    logging.info(f"Coverage plot generated: {coverage_plot_file}")

def prepare_reference(reference_fasta, output_dir):
    validate_fasta(reference_fasta)
    if not os.path.exists(f"{reference_fasta}.fai"):
        logging.info("Indexing reference FASTA...")
        run_command(f"samtools faidx {reference_fasta}")

    dict_file = os.path.splitext(reference_fasta)[0] + ".dict"
    if not os.path.exists(dict_file):
        logging.info("Creating sequence dictionary for reference...")
        run_command(f"gatk CreateSequenceDictionary -R {reference_fasta} -O {dict_file}")

def add_read_groups(bam_file, sample_name, output_dir):
    rg_bam = os.path.join(output_dir, f"{sample_name}_rg.bam")
    rg_command = (
        f"samtools addreplacerg -r 'ID:{sample_name}\tSM:{sample_name}\tLB:lib1\tPL:ILLUMINA' "
        f"-o {rg_bam} {bam_file}"
    )
    stdout, stderr = run_command(rg_command)
    print("Add read groups output:", stdout)
    print("Add read groups error:", stderr)

    # Validate and sort the BAM file
    run_command(f"samtools sort -T {os.path.join(output_dir, f'temp_{sample_name}')} -o {rg_bam}.sorted {rg_bam}")
    os.rename(f"{rg_bam}.sorted", rg_bam)
    run_command(f"samtools index {rg_bam}")
    return rg_bam

def run_variant_calling(bam_file, reference_fasta, sample_name, output_dir, config):
    validate_bam(bam_file)
    raw_vcf = os.path.join(output_dir, f"{sample_name}.vcf")
    vc = config['variant_calling']
    gatk_command = (
        f"gatk HaplotypeCaller "
        f"-R {reference_fasta} "
        f"-I {bam_file} "
        f"-O {raw_vcf} "
        f"-ploidy {vc['ploidy']} "
        f"--standard-min-confidence-threshold-for-calling {vc['standard_min_confidence']} "
        f"--min-base-quality-score {vc['min_base_quality_score']}"
    )
    stdout, stderr = run_command(gatk_command)
    logging.info(f"GATK HaplotypeCaller output: {stdout}")
    logging.info(f"GATK HaplotypeCaller error: {stderr}")
    return raw_vcf

def filter_vcf(raw_vcf, reference_fasta, sample_name, output_dir, config):
    """Filter VCF for quality, depth, and strand bias."""
    vf = config['vcf_filtering']
    filtered_vcf = os.path.join(output_dir, f"{sample_name}_filtered.vcf")

    # Build filter expressions from config
    filter_parts = []
    for filter_name, expression in vf['filters'].items():
        filter_parts.append(f"--filter-expression '{expression}' --filter-name '{filter_name}'")

    filter_command = (
        f"gatk VariantFiltration "
        f"-R {reference_fasta} "
        f"-V {raw_vcf} "
        f"{' '.join(filter_parts)} "
        f"-O {filtered_vcf}"
    )
    run_command(filter_command)

    # Select only PASS variants
    pass_vcf = os.path.join(output_dir, f"{sample_name}_pass.vcf")
    if vf['select_pass_only']:
        select_command = (
            f"gatk SelectVariants "
            f"-R {reference_fasta} "
            f"-V {filtered_vcf} "
            f"--exclude-filtered "
            f"-O {pass_vcf}"
        )
        run_command(select_command)
    else:
        pass_vcf = filtered_vcf

    logging.info(f"VCF filtering complete: {pass_vcf}")
    return pass_vcf

def trim_primers(bam_file, primer_bed, sample_name, output_dir, config):
    """Trim primer sequences from aligned reads using ivar trim."""
    pt = config['primer_trimming']
    trimmed_prefix = os.path.join(output_dir, f"{sample_name}_trimmed")
    trimmed_bam = f"{trimmed_prefix}.bam"
    trimmed_sorted = os.path.join(output_dir, f"{sample_name}_trimmed.sorted.bam")

    trim_parts = [
        f"ivar trim -b {primer_bed}",
        f"-p {trimmed_prefix}",
        f"-i {bam_file}",
        f"-q {pt['min_quality']} -m {pt['min_length']} -s {pt['sliding_window']}",
    ]
    if pt['include_reads_no_primer']:
        trim_parts.append("-e")

    trim_command = " ".join(trim_parts)
    run_command(trim_command)

    # Sort and index the trimmed BAM
    run_command(f"samtools sort -o {trimmed_sorted} {trimmed_bam}")
    run_command(f"samtools index {trimmed_sorted}")

    # Clean up unsorted trimmed BAM
    if os.path.exists(trimmed_bam):
        os.remove(trimmed_bam)

    logging.info(f"Primer trimming complete: {trimmed_sorted}")
    return trimmed_sorted

def run_snpeff_annotation(raw_vcf, sample_name, output_dir, config, reference_name="denv1"):
    ann = config['annotation']
    annotated_vcf = os.path.join(output_dir, f"{sample_name}_annotated.vcf")
    summary_html = os.path.join(output_dir, f"{sample_name}_snpEff_summary.html")
    summary_csv = os.path.join(output_dir, f"{sample_name}_snpEff_summary.csv")
    snpeff_command = (
        f"snpEff -Xmx{ann['snpeff_memory']} "
        f"-c {os.path.join(output_dir, 'snpEff.config')} "
        f"-v {reference_name} "
        f"-s {summary_html} "
        f"-csvStats {summary_csv} "
        f"{raw_vcf} > {annotated_vcf}"
    )
    stdout, stderr = run_command(snpeff_command)
    logging.info(f"SnpEff annotation output: {stdout}")
    logging.info(f"SnpEff annotation error: {stderr}")
    return annotated_vcf, summary_html, summary_csv, os.path.join(output_dir, f"{sample_name}_snpEff_summary.genes.txt")


def write_low_coverage_positions(coverage_file, output_dir, sample_name, min_depth=20):
    """Write a file listing positions with coverage below the consensus threshold."""
    low_cov_file = os.path.join(output_dir, f"{sample_name}_low_coverage.tsv")
    low_cov_count = 0
    total_positions = 0

    with open(coverage_file, 'r') as fin, open(low_cov_file, 'w') as fout:
        fout.write("CHROM\tPOSITION\tDEPTH\n")
        for line in fin:
            columns = line.strip().split('\t')
            chrom = columns[0]
            pos = int(columns[1])
            depth = int(columns[2])
            total_positions += 1
            if depth < min_depth:
                fout.write(f"{chrom}\t{pos}\t{depth}\n")
                low_cov_count += 1

    logging.info(f"Low coverage file: {low_cov_file} "
                 f"({low_cov_count}/{total_positions} positions below {min_depth}X)")
    return low_cov_file

def create_annotation_tsv(annotated_vcf, sample_name, output_dir, config):
    """Parse annotated VCF and create user-friendly TSV with one row per variant."""
    transcript_map = config.get('transcript_annotations', {})
    output_file = os.path.join(output_dir, f"{sample_name}_annotations.tsv")

    header = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
        "Total_Depth", "Allele_Frequency", "Strand_Bias", "Allelic_Depths",
        "EFFECT", "PUTATIVE_IMPACT", "GENE_NAME", "GENE_ID",
        "FEATURE_TYPE", "FEATURE_ID", "TRANSCRIPT_TYPE",
        "HGVSc", "HGVSp",
        "cDNA_POSITION_AND_LENGTH", "CDS_POSITION_AND_LENGTH",
        "PROTEIN_POSITION_AND_LENGTH", "ERROR"
    ]

    rows = []
    with open(annotated_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom = fields[0]
            pos = fields[1]
            vid = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filt = fields[6]
            info = fields[7]

            # Parse INFO field
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v

            total_depth = info_dict.get('DP', '')
            allele_freq = info_dict.get('AF', '')
            strand_bias = info_dict.get('FS', '')

            # Get AD from FORMAT field (ref,alt depths)
            allelic_depths = ''
            if len(fields) > 9:
                fmt_keys = fields[8].split(':')
                fmt_vals = fields[9].split(':')
                fmt_dict = dict(zip(fmt_keys, fmt_vals))
                allelic_depths = fmt_dict.get('AD', '')

            # Parse ANN field - take the FIRST annotation (highest impact)
            ann_str = info_dict.get('ANN', '')
            effect = impact = gene_name = gene_id = ''
            feature_type = feature_id = transcript_type = ''
            hgvsc = hgvsp = cdna_pos = cds_pos = protein_pos = error = ''

            if ann_str:
                annotations = ann_str.split(',')

                # Score impacts for comparison
                impact_rank = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}

                # Find the best non-Polyprotein annotation (highest impact)
                best_specific = None
                best_specific_rank = 0
                polyprotein_ann = None

                for ann in annotations:
                    parts = ann.split('|')
                    if len(parts) >= 16:
                        feat_id = parts[6]
                        mapped_name = transcript_map.get(feat_id, '')
                        rank = impact_rank.get(parts[2], 0)

                        # Save the polyprotein annotation as fallback
                        if mapped_name == 'Polyprotein' and polyprotein_ann is None:
                            polyprotein_ann = parts

                        # Track the best non-Polyprotein annotation by impact
                        if mapped_name and mapped_name != 'Polyprotein' and rank > best_specific_rank:
                            best_specific = parts
                            best_specific_rank = rank

                # Decision logic:
                # 1. If specific protein has impact >= polyprotein's, use it
                # 2. If polyprotein has better impact, use polyprotein annotation
                poly_rank = impact_rank.get(polyprotein_ann[2], 0) if polyprotein_ann else 0

                if best_specific and best_specific_rank >= poly_rank:
                    best_ann = best_specific
                elif polyprotein_ann:
                    best_ann = polyprotein_ann
                else:
                    best_ann = annotations[0].split('|') if annotations else None

                if best_ann and len(best_ann) >= 16:
                    effect = best_ann[1]
                    impact = best_ann[2]
                    gene_name_raw = best_ann[3]
                    gene_id = best_ann[4]
                    feature_type = best_ann[5]
                    feature_id = best_ann[6]
                    transcript_type = best_ann[7]
                    hgvsc = best_ann[9]
                    hgvsp = best_ann[10]
                    cdna_pos = best_ann[11]
                    cds_pos = best_ann[12]
                    protein_pos = best_ann[13]
                    error = best_ann[15] if len(best_ann) > 15 else ''
                    gene_name = transcript_map.get(feature_id, gene_name_raw)

            rows.append([
                chrom, pos, vid, ref, alt, qual, filt,
                total_depth, allele_freq, strand_bias, allelic_depths,
                effect, impact, gene_name, gene_id,
                feature_type, feature_id, transcript_type,
                hgvsc, hgvsp, cdna_pos, cds_pos, protein_pos, error
            ])

    with open(output_file, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')

    logging.info(f"Annotation TSV created: {output_file} ({len(rows)} variants)")
    return output_file

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing BAM files.')
    parser.add_argument('--reference_fasta', type=str, help='Path to reference FASTA file.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for consensus FASTA and VCF files.')
    parser.add_argument('--database_name', type=str, default="denv1", help='Name of the SnpEff database (default: denv1).')
    parser.add_argument('--config', type=str, required=True, help='Path to virus config YAML file')
    parser.add_argument('--primer_bed', type=str, default=None,
                        help='BED file with primer coordinates for ivar trim. '
                             'If not provided, primer trimming is skipped.')
    args = parser.parse_args(argv)

    config = load_config(args.config)
    cov = config['coverage']
    cons = config['consensus']

    input_dir = args.input_dir
    reference_fasta = args.reference_fasta
    output_dir = args.output_dir
    database_name = args.database_name
    primer_bed = args.primer_bed

    if primer_bed and not os.path.exists(primer_bed):
        logging.error(f"Primer BED file not found: {primer_bed}")
        sys.exit(1)
    if primer_bed:
        logging.info(f"Primer trimming enabled with BED file: {primer_bed}")
    else:
        logging.warning("No primer BED file provided -- skipping primer trimming")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bam_files = glob.glob(os.path.join(input_dir, "*.sorted.bam"))
    if not bam_files:
        logging.error(f"No sorted BAM files found in {input_dir}")
        sys.exit(1)

    prepare_reference(reference_fasta, output_dir)
    bam_files = glob.glob(os.path.join(input_dir, "*sorted.bam"))
    for bam_file in bam_files:
        print("Processing:", bam_file)
        sample_name = os.path.splitext(os.path.basename(bam_file))[0].replace('.sorted', '')
        try:
            # Add read groups before GATK
            rg_bam = add_read_groups(bam_file, sample_name, output_dir)

            # Primer trimming (if BED file provided)
            analysis_bam = rg_bam
            if primer_bed:
                analysis_bam = trim_primers(rg_bam, primer_bed, sample_name, output_dir, config)

            # Coverage calculation with quality filters from config
            coverage_file = os.path.join(output_dir, f"{sample_name}_coverage.txt")
            coverage_command = (
                f"samtools depth -a -q {cov['min_base_quality']} -Q {cov['min_mapping_quality']} "
                f"{analysis_bam} > {coverage_file}"
            )
            stdout, stderr = run_command(coverage_command)
            logging.info(f"Coverage command output: {stdout}")
            logging.info(f"Coverage command error: {stderr}")
            plot_coverage(coverage_file, output_dir, sample_name=sample_name)
            write_low_coverage_positions(coverage_file, output_dir, sample_name)

            # Flagstat for on-target read metrics
            flagstat_file = os.path.join(output_dir, f"{sample_name}_flagstat.txt")
            flagstat_command = f"samtools flagstat {analysis_bam} > {flagstat_file}"
            stdout, stderr = run_command(flagstat_command)
            logging.info(f"Flagstat output saved to {flagstat_file}")

            # ivar consensus with parameters from config
            pileup_command = (
                f"samtools mpileup -aa -A -d {cons['mpileup_max_depth']} "
                f"-Q {cons['mpileup_min_base_quality']} "
                f"-q {cons['mpileup_min_mapping_quality']} "
                f"{analysis_bam} | "
                f"ivar consensus -p {sample_name} "
                f"-q {cons['ivar_min_quality']} "
                f"-t {cons['ivar_min_frequency']} "
                f"-m {cons['ivar_min_depth']} "
                # -n omitted: default is N
            )
            stdout, stderr = run_command(pileup_command)
            logging.info(f"Consensus command output: {stdout}")
            logging.info(f"Consensus command error: {stderr}")
            consensus_fasta = f"{sample_name}.fa"
            os.rename(consensus_fasta, os.path.join(output_dir, consensus_fasta))

            # Variant calling (ploidy from config)
            raw_vcf = run_variant_calling(analysis_bam, reference_fasta, sample_name, output_dir, config)

            # VCF filtering (thresholds from config)
            filtered_vcf = filter_vcf(raw_vcf, reference_fasta, sample_name, output_dir, config)

            # Annotate filtered variants
            annotated_vcf, summary_html, summary_csv, summary_txt = run_snpeff_annotation(
                filtered_vcf, sample_name, output_dir, config, database_name)
            annotation_tsv = create_annotation_tsv(annotated_vcf, sample_name, output_dir, config)
            logging.info(
                f"Processing complete for {sample_name}: "
                f"consensus={os.path.join(output_dir, consensus_fasta)}, "
                f"coverage={coverage_file}, filtered_vcf={filtered_vcf}, "
                f"annotated_vcf={annotated_vcf}"
            )
        except Exception as e:
            print(f"Error occurred during processing {sample_name}: {str(e)}")
            continue

if __name__ == '__main__':
    main()
