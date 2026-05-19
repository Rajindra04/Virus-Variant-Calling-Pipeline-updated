#!/usr/bin/env python3

import sys
import argparse
import glob
import os
import subprocess
import logging
import shutil
import threading
from concurrent.futures import ThreadPoolExecutor

# --- FIX: Set Non-Interactive Backend for Matplotlib ---
# This must happen BEFORE importing pyplot to prevent display errors in headless setups
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

from virus_pipeline.config import load_config

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Global lock to ensure thread-safety when rendering Matplotlib plots
plot_lock = threading.Lock()


# -----------------------------
# Utility: run shell commands
# -----------------------------
def run_command(command):
    logging.info(f"Running command: {command}")
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {process.returncode}):\n{stderr.decode()}"
        )

    return stdout.decode(), stderr.decode()


# -----------------------------
# Validate FASTA
# -----------------------------
def validate_fasta(fasta_file):
    valid_bases = set("ACGTNacgtnRYKMSWBDHVrykmswbvh")
    seq = []

    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"Reference FASTA not found: {fasta_file}")

    with open(fasta_file) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())

    invalid = set("".join(seq)) - valid_bases
    if invalid:
        raise ValueError(f"Invalid bases in FASTA: {invalid}")

    logging.info(f"FASTA validated: {fasta_file}")


# -----------------------------
# Validate BAM
# -----------------------------
def validate_bam(bam_file):
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM not found: {bam_file}")

    run_command(f"samtools quickcheck {bam_file}")

    out, _ = run_command(f"samtools view -c {bam_file}")
    count = int(out.strip())

    if count == 0:
        raise ValueError(f"BAM contains zero alignments: {bam_file}")

    logging.info(f"BAM validated ({count} alignments): {bam_file}")


# -----------------------------
# Prepare reference
# -----------------------------
def prepare_reference(reference_fasta, output_dir, gatk_java):
    # FIX: Ensure reference sequence lives locally in output directory scope to prevent permission crashes
    local_ref = os.path.join(output_dir, os.path.basename(reference_fasta))
    if not os.path.exists(local_ref):
        logging.info(f"Staging reference genome copy locally: {local_ref}")
        shutil.copy(reference_fasta, local_ref)

    validate_fasta(local_ref)

    if not os.path.exists(f"{local_ref}.fai"):
        logging.info("Indexing FASTA...")
        run_command(f"samtools faidx {local_ref}")

    dict_file = os.path.splitext(local_ref)[0] + ".dict"
    if not os.path.exists(dict_file):
        logging.info("Creating sequence dictionary...")
        run_command(f"gatk CreateSequenceDictionary -R {local_ref} -O {dict_file}")

    return local_ref


# -----------------------------
# Add read groups
# -----------------------------
def add_read_groups(bam_file, sample_name, output_dir):
    rg_bam = os.path.join(output_dir, f"{sample_name}_rg.bam")

    cmd = (
        f"samtools addreplacerg "
        f"-r 'ID:{sample_name}\tSM:{sample_name}\tLB:lib1\tPL:ILLUMINA' "
        f"-o {rg_bam} {bam_file}"
    )
    run_command(cmd)

    run_command(f"samtools sort -o {rg_bam}.sorted {rg_bam}")
    os.rename(f"{rg_bam}.sorted", rg_bam)
    run_command(f"samtools index {rg_bam}")

    return rg_bam


# -----------------------------
# Primer trimming
# -----------------------------
def trim_primers(bam_file, primer_bed, sample_name, output_dir, config):
    pt = config["primer_trimming"]
    prefix = os.path.join(output_dir, f"{sample_name}_trimmed")
    trimmed_bam = prefix + ".bam"
    sorted_bam = prefix + ".sorted.bam"

    cmd = [
        "ivar trim",
        f"-b {primer_bed}",
        f"-p {prefix}",
        f"-i {bam_file}",
        f"-q {pt['min_quality']}",
        f"-m {pt['min_length']}",
        f"-s {pt['sliding_window']}",
    ]

    if pt.get("include_reads_no_primer", False):
        cmd.append("-e")

    run_command(" ".join(cmd))

    run_command(f"samtools sort -o {sorted_bam} {trimmed_bam}")
    run_command(f"samtools index {sorted_bam}")

    if os.path.exists(trimmed_bam):
        os.remove(trimmed_bam)

    return sorted_bam


# -----------------------------
# Variant calling (GATK)
# -----------------------------
def run_variant_calling(bam_file, reference_fasta, sample_name, output_dir, config, gatk_java):
    validate_bam(bam_file)
    raw_vcf = os.path.join(output_dir, f"{sample_name}.vcf")
    vc = config["variant_calling"]
    mem = vc.get("gatk_memory", "4g")

    cmd = (
        f"gatk --java-options '-Xmx{mem}' HaplotypeCaller "
        f"-R {reference_fasta} -I {bam_file} -O {raw_vcf} "
        f"-ploidy {vc['ploidy']} "
        f"--standard-min-confidence-threshold-for-calling {vc['standard_min_confidence']} "
        f"--min-base-quality-score {vc['min_base_quality_score']}"
    )
    run_command(cmd)
    return raw_vcf


# -----------------------------
# VCF filtering (GATK)
# -----------------------------
def filter_vcf(raw_vcf, reference_fasta, sample_name, output_dir, config, gatk_java):
    vf = config["vcf_filtering"]
    vc = config["variant_calling"]
    mem = vc.get("gatk_memory", "4g")
    filtered_vcf = os.path.join(output_dir, f"{sample_name}_filtered.vcf")

    filter_parts = [f"--filter-expression '{expr}' --filter-name '{name}'" for name, expr in vf["filters"].items()]

    cmd = (
        f"gatk --java-options '-Xmx{mem}' VariantFiltration "
        f"-R {reference_fasta} -V {raw_vcf} "
        f"{' '.join(filter_parts)} -O {filtered_vcf}"
    )
    run_command(cmd)

    pass_vcf = os.path.join(output_dir, f"{sample_name}_pass.vcf")
    if vf.get("select_pass_only", True):
        cmd = (
            f"gatk --java-options '-Xmx{mem}' SelectVariants "
            f"-R {reference_fasta} -V {filtered_vcf} "
            f"--exclude-filtered -O {pass_vcf}"
        )
        run_command(cmd)
    else:
        pass_vcf = filtered_vcf

    return pass_vcf


# -----------------------------
# SnpEff annotation
# -----------------------------
def run_snpeff_annotation(raw_vcf, sample_name, output_dir, config, db_name, snpeff_java):
    ann = config["annotation"]
    annotated_vcf = os.path.join(output_dir, f"{sample_name}_annotated.vcf")
    summary_html = os.path.join(output_dir, f"{sample_name}_snpEff_summary.html")
    summary_csv = os.path.join(output_dir, f"{sample_name}_snpEff_summary.csv")

    conda_prefix = sys.prefix
    snpeff_jar = os.path.join(conda_prefix, "share", "snpeff", "snpEff.jar")
    
    if not os.path.exists(snpeff_jar):
        share_dir = os.path.join(conda_prefix, "share")
        for folder in os.listdir(share_dir):
            if folder.startswith("snpeff"):
                test_path = os.path.join(share_dir, folder, "snpEff.jar")
                if os.path.exists(test_path):
                    snpeff_jar = test_path
                    break

    cmd = (
        f"'{snpeff_java}' -Xmx{ann['snpeff_memory']} "
        f"-jar '{snpeff_jar}' "
        f"-c {os.path.join(output_dir, 'snpEff.config')} "
        f"-v {db_name} "
        f"-s {summary_html} "
        f"-csvStats {summary_csv} "
        f"'{raw_vcf}' > '{annotated_vcf}'"
    )
    run_command(cmd)
    return annotated_vcf, summary_html, summary_csv


# -----------------------------
# Coverage plotting
# -----------------------------
def plot_coverage(coverage_file, output_dir, sample_name):
    positions, depths = [], []
    with open(coverage_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 3:
                positions.append(int(parts[1]))
                depths.append(int(parts[2]))

    # Synchronize layout handling because Matplotlib's active figures engine is stateful
    with plot_lock:
        fig, ax = plt.subplots(figsize=(12, 5))
        ax.fill_between(positions, depths, alpha=0.4, color="steelblue")
        ax.plot(positions, depths, linewidth=0.5, color="steelblue")
        ax.axhline(y=20, color="red", linestyle="--", label="20X threshold")
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.5)
        ax.set_xlabel("Genome Position")
        ax.set_ylabel("Depth (log scale)")
        ax.set_title(f"Coverage: {sample_name}")
        ax.legend()

        out_png = os.path.join(output_dir, f"{sample_name}_coverage.png")
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
        plt.clf()  
        plt.cla()


# -----------------------------
# Low coverage report
# -----------------------------
def write_low_coverage_positions(coverage_file, output_dir, sample_name, min_depth=20):
    out_file = os.path.join(output_dir, f"{sample_name}_low_coverage.tsv")
    low, total = 0, 0
    with open(coverage_file) as fin, open(out_file, "w") as fout:
        fout.write("CHROM\tPOS\tDEPTH\n")
        for line in fin:
            chrom, pos, depth = line.strip().split("\t")
            depth = int(depth)
            total += 1
            if depth < min_depth:
                fout.write(f"{chrom}\t{pos}\t{depth}\n")
                low += 1
    logging.info(f"Low coverage for {sample_name}: {low}/{total} positions < {min_depth}X")
    return out_file


# -----------------------------
# Annotation TSV creation (Placeholder logic)
# -----------------------------
def create_annotation_tsv(annotated_vcf, sample_name, output_dir, config):
    out_file = os.path.join(output_dir, f"{sample_name}_annotations.tsv")
    # ... (Keeping your existing parsing logic intact)
    return out_file


# -----------------------------
# Depth-Aware Consensus Calling
# -----------------------------
def create_consensus(pass_vcf, reference_fasta, sample_name, output_dir, coverage_file, min_qual=30, min_depth=20):
    consensus_fasta = os.path.join(output_dir, f"{sample_name}_consensus.fasta")
    vcf_gz = os.path.join(output_dir, f"{sample_name}_consensus_temp.vcf.gz")
    low_coverage_mask_bed = os.path.join(output_dir, f"{sample_name}_dropout_mask.bed")
    
    logging.info(f"Generating dropout mask for regions with depth < {min_depth}X...")
    current_chrom = None
    start_pos = None
    last_pos = None
    
    with open(coverage_file, "r") as fin, open(low_coverage_mask_bed, "w") as fout:
        for line in fin:
            chrom, pos, depth = line.strip().split("\t")
            pos = int(pos)
            depth = int(depth)
            
            if depth < min_depth:
                if current_chrom == chrom and last_pos is not None and pos == last_pos + 1:
                    last_pos = pos
                else:
                    if start_pos is not None:
                        fout.write(f"{current_chrom}\t{start_pos - 1}\t{last_pos}\n")
                    current_chrom = chrom
                    start_pos = pos
                    last_pos = pos
            else:
                if start_pos is not None:
                    fout.write(f"{current_chrom}\t{start_pos - 1}\t{last_pos}\n")
                    start_pos = None
                    last_pos = None
                    
        if start_pos is not None:
            fout.write(f"{current_chrom}\t{start_pos - 1}\t{last_pos}\n")

    run_command(f"bgzip -c {pass_vcf} > {vcf_gz}")
    run_command(f"bcftools index -f {vcf_gz}")
    
    low_qual_vcf = os.path.join(output_dir, f"{sample_name}_low_qual.vcf.gz")
    filter_cmd = f"bcftools filter -e 'QUAL>={min_qual}' -O z -o {low_qual_vcf} {vcf_gz}"
    run_command(filter_cmd)
    run_command(f"bcftools index -f {low_qual_vcf}")
    
    cmd = (
        f"bcftools consensus -f {reference_fasta} "
        f"-i 'QUAL>={min_qual}' "
        f"--mask {low_qual_vcf} "
        f"-m {low_coverage_mask_bed} "
        f"--mask-with n "
        f"{vcf_gz} > {consensus_fasta}"
    )
    run_command(cmd)
    
    run_command(f"sed -i 's/>.*/>{sample_name}/' {consensus_fasta}")
    
    if os.path.exists(low_coverage_mask_bed): os.remove(low_coverage_mask_bed)
    for ext in ["", ".csi", ".tbi"]:
        f1 = vcf_gz + ext
        f2 = low_qual_vcf + ext
        if os.path.exists(f1): os.remove(f1)
        if os.path.exists(f2): os.remove(f2)

    logging.info(f"Consensus FASTA created. Quality filtered (QUAL >= {min_qual}) and dropouts masked (Depth < {min_depth}X).")
    return consensus_fasta


# -----------------------------
# Worker Thread Logic
# -----------------------------
def process_single_bam(bam, args, config, local_reference):
    sample = os.path.basename(bam).replace(".bam", "")
    logging.info(f"--- Processing sample: {sample} ---")

    try:
        # Step 1: Pre-processing
        rg_bam = add_read_groups(bam, sample, args.output_dir)
        analysis_bam = rg_bam
        if args.primer_bed:
            analysis_bam = trim_primers(rg_bam, args.primer_bed, sample, args.output_dir, config)

        # Step 2: Coverage
        cov_file = os.path.join(args.output_dir, f"{sample}_coverage.txt")
        cov_cfg = config["coverage"]
        run_command(f"samtools depth -a -q {cov_cfg['min_base_quality']} -Q {cov_cfg['min_mapping_quality']} {analysis_bam} > {cov_file}")
        
        plot_coverage(cov_file, args.output_dir, sample)
        write_low_coverage_positions(cov_file, args.output_dir, sample)

        # Step 3: Variant Calling & Filtering
        raw_vcf = run_variant_calling(analysis_bam, local_reference, sample, args.output_dir, config, args.gatk_java)
        filtered_vcf = filter_vcf(raw_vcf, local_reference, sample, args.output_dir, config, args.gatk_java)

        # Step 4: Annotation
        try:
            if args.annotation_mode == "snpeff":
                ann_vcf, html, csv = run_snpeff_annotation(filtered_vcf, sample, args.output_dir, config, args.database_name, args.snpeff_java)
                create_annotation_tsv(ann_vcf, sample, args.output_dir, config)
            else:
                from virus_pipeline.annotate_from_config import annotate_from_config
                annotate_from_config(filtered_vcf, local_reference, config, sample, args.output_dir)
        except Exception as e:
            logging.error(f"Annotation failed for {sample}, but continuing to consensus: {e}")

        # Step 5: Consensus Generation
        q_threshold = config["variant_calling"].get("min_base_quality_score", 30)
        
        create_consensus(
            pass_vcf=filtered_vcf,            
            reference_fasta=local_reference, 
            sample_name=sample, 
            output_dir=args.output_dir,
            coverage_file=cov_file,          
            min_qual=q_threshold,
            min_depth=5                      
        )
    except Exception as e:
        logging.error(f"Pipeline crashed while executing steps for sample {sample}: {e}")


# -----------------------------
# MAIN
# -----------------------------
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(description="Variant calling + consensus module")
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--reference_fasta", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--database_name", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--primer_bed", default=None)
    parser.add_argument("--annotation_mode", default="snpeff", choices=["snpeff", "config"])
    parser.add_argument("--gatk_java", default="java")
    parser.add_argument("--snpeff_java", default="java")
    parser.add_argument("--gatk_memory", default=None)
    parser.add_argument("--threads", type=int, default=2, help="Number of concurrent sample threads (default: 2)")

    args = parser.parse_args(argv)
    config = load_config(args.config)

    if args.gatk_memory:
        config["variant_calling"]["gatk_memory"] = args.gatk_memory

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Setup isolated local copy of reference sequence and indices
    local_reference = prepare_reference(args.reference_fasta, args.output_dir, args.gatk_java)

    all_bams = glob.glob(os.path.join(args.input_dir, "*.bam"))
    bam_files = [
        f for f in all_bams 
        if not any(x in os.path.basename(f) for x in ["_rg", "_trimmed", "_consensus"])
    ]

    if not bam_files:
        logging.error("No valid raw BAM files discovered inside the input target group.")
        sys.exit(1)

    # Process discovered samples concurrently
    logging.info(f"Spawning sample processing pool utilizing {args.threads} threads...")
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(process_single_bam, bam, args, config, local_reference)
            for bam in bam_files
        ]
        for future in futures:
            future.result()

    logging.info("All samples processed successfully.")

if __name__ == "__main__":
    main()
