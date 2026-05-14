#!/usr/bin/env python3

import sys
import argparse
import glob
import os
import subprocess
import logging
import shutil
import matplotlib.pyplot as plt

from virus_pipeline.config import load_config

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


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

    # samtools quickcheck
    run_command(f"samtools quickcheck {bam_file}")

    # alignment count
    out, _ = run_command(f"samtools view -c {bam_file}")
    count = int(out.strip())

    if count == 0:
        raise ValueError(f"BAM contains zero alignments: {bam_file}")

    logging.info(f"BAM validated ({count} alignments): {bam_file}")


# -----------------------------
# Prepare reference (FASTA index + dict)
# -----------------------------
def prepare_reference(reference_fasta, output_dir, gatk_java):
    validate_fasta(reference_fasta)

    if not os.path.exists(f"{reference_fasta}.fai"):
        logging.info("Indexing FASTA...")
        run_command(f"samtools faidx {reference_fasta}")

    dict_file = os.path.splitext(reference_fasta)[0] + ".dict"
    if not os.path.exists(dict_file):
        logging.info("Creating sequence dictionary...")
        # Use the 'gatk' wrapper directly
        run_command(f"gatk CreateSequenceDictionary -R {reference_fasta} -O {dict_file}")


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

    # sort + index
    sorted_bam = rg_bam + ".sorted"
    run_command(f"samtools sort -o {sorted_bam} {rg_bam}")
    os.rename(sorted_bam, rg_bam)
    run_command(f"samtools index {rg_bam}")

    return rg_bam


# -----------------------------
# Primer trimming (optional)
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

    # sort + index
    run_command(f"samtools sort -o {sorted_bam} {trimmed_bam}")
    run_command(f"samtools index {sorted_bam}")

    if os.path.exists(trimmed_bam):
        os.remove(trimmed_bam)

    logging.info(f"Primer trimming complete: {sorted_bam}")
    return sorted_bam


# -----------------------------
# Variant calling (GATK)
# -----------------------------
def run_variant_calling(bam_file, reference_fasta, sample_name, output_dir, config, gatk_java):
    validate_bam(bam_file)

    raw_vcf = os.path.join(output_dir, f"{sample_name}.vcf")
    vc = config["variant_calling"]
    mem = vc.get("gatk_memory", "4g")

    # Use the 'gatk' wrapper with --java-options
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

    filter_parts = []
    for name, expr in vf["filters"].items():
        filter_parts.append(f"--filter-expression '{expr}' --filter-name '{name}'")

    # Use 'gatk' wrapper
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

    # FIX: Dynamically locate the JAR in the Conda environment
    conda_prefix = sys.prefix
    snpeff_jar = os.path.join(conda_prefix, "share", "snpeff", "snpEff.jar")
    
    # Fallback search for versioned folders (e.g., snpeff-5.2)
    if not os.path.exists(snpeff_jar):
        share_dir = os.path.join(conda_prefix, "share")
        if os.path.exists(share_dir):
            for folder in os.listdir(share_dir):
                if folder.startswith("snpeff"):
                    test_path = os.path.join(share_dir, folder, "snpEff.jar")
                    if os.path.exists(test_path):
                        snpeff_jar = test_path
                        break

    if not os.path.exists(snpeff_jar):
        raise FileNotFoundError(
            f"SnpEff JAR not found in Conda environment: {snpeff_jar}"
        )

    # Construct command using absolute path to the modern Java and the JAR
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
    positions = []
    depths = []

    with open(coverage_file) as f:
        for line in f:
            chrom, pos, depth = line.strip().split("\t")
            positions.append(int(pos))
            depths.append(int(depth))

    import matplotlib.pyplot as plt

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

    logging.info(f"Coverage plot saved: {out_png}")


# -----------------------------
# Low coverage report
# -----------------------------
def write_low_coverage_positions(coverage_file, output_dir, sample_name, min_depth=20):
    out_file = os.path.join(output_dir, f"{sample_name}_low_coverage.tsv")
    low = 0
    total = 0

    with open(coverage_file) as fin, open(out_file, "w") as fout:
        fout.write("CHROM\tPOS\tDEPTH\n")
        for line in fin:
            chrom, pos, depth = line.strip().split("\t")
            depth = int(depth)
            total += 1
            if depth < min_depth:
                fout.write(f"{chrom}\t{pos}\t{depth}\n")
                low += 1

    logging.info(f"Low coverage: {low}/{total} positions < {min_depth}X")
    return out_file


# -----------------------------
# Annotation TSV creation
# -----------------------------
def create_annotation_tsv(annotated_vcf, sample_name, output_dir, config):
    transcript_map = config.get("transcript_annotations", {})
    out_file = os.path.join(output_dir, f"{sample_name}_annotations.tsv")

    header = [
        "CHROM","POS","ID","REF","ALT","QUAL","FILTER",
        "DP","AF","FS","AD",
        "EFFECT","IMPACT","GENE","GENE_ID",
        "FEATURE_TYPE","FEATURE_ID","TRANSCRIPT_TYPE",
        "HGVSc","HGVSp","cDNA","CDS","PROTEIN","ERROR"
    ]

    rows = []

    with open(annotated_vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]

            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v

            dp = info_dict.get("DP", "")
            af = info_dict.get("AF", "")
            fs = info_dict.get("FS", "")

            # FORMAT field
            ad = ""
            if len(cols) > 9:
                fmt_keys = cols[8].split(":")
                fmt_vals = cols[9].split(":")
                fmt = dict(zip(fmt_keys, fmt_vals))
                ad = fmt.get("AD", "")

            # ANN parsing
            ann = info_dict.get("ANN", "")
            effect = impact = gene = gene_id = ""
            feature_type = feature_id = transcript_type = ""
            hgvsc = hgvsp = cdna = cds = protein = err = ""

            if ann:
                first = ann.split(",")[0].split("|")
                if len(first) >= 16:
                    effect = first[1]
                    impact = first[2]
                    gene = first[3]
                    gene_id = first[4]
                    feature_type = first[5]
                    feature_id = first[6]
                    transcript_type = first[7]
                    hgvsc = first[9]
                    hgvsp = first[10]
                    cdna = first[11]
                    cds = first[12]
                    protein = first[13]
                    err = first[15] if len(first) > 15 else ""

                    # map gene name if provided
                    gene = transcript_map.get(feature_id, gene)

            rows.append([
                chrom, pos, vid, ref, alt, qual, filt,
                dp, af, fs, ad,
                effect, impact, gene, gene_id,
                feature_type, feature_id, transcript_type,
                hgvsc, hgvsp, cdna, cds, protein, err
            ])

    with open(out_file, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")

    logging.info(f"Annotation TSV written: {out_file}")
    return out_file


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

    args = parser.parse_args(argv)

    config = load_config(args.config)

    # Override memory if provided
    if args.gatk_memory:
        config["variant_calling"]["gatk_memory"] = args.gatk_memory

    # Prepare reference
    prepare_reference(args.reference_fasta, args.output_dir, args.gatk_java)

    # Process BAM files
    bam_files = glob.glob(os.path.join(args.input_dir, "*.bam"))
    if not bam_files:
        logging.error("No BAM files found")
        sys.exit(1)

    for bam in bam_files:
        sample = os.path.basename(bam).replace(".bam", "")
        logging.info(f"Processing sample: {sample}")

        # Add read groups
        rg_bam = add_read_groups(bam, sample, args.output_dir)

        # Primer trimming
        analysis_bam = rg_bam
        if args.primer_bed:
            analysis_bam = trim_primers(rg_bam, args.primer_bed, sample, args.output_dir, config)

        # Coverage
        cov_file = os.path.join(args.output_dir, f"{sample}_coverage.txt")
        cov_cfg = config["coverage"]

        cmd = (
            f"samtools depth -a -q {cov_cfg['min_base_quality']} "
            f"-Q {cov_cfg['min_mapping_quality']} "
            f"{analysis_bam} > {cov_file}"
        )
        run_command(cmd)

        plot_coverage(cov_file, args.output_dir, sample)
        write_low_coverage_positions(cov_file, args.output_dir, sample)

        # Variant calling
        raw_vcf = run_variant_calling(
            analysis_bam, args.reference_fasta, sample, args.output_dir, config, args.gatk_java
        )

        # Filtering
        filtered_vcf = filter_vcf(
            raw_vcf, args.reference_fasta, sample, args.output_dir, config, args.gatk_java
        )

        # Annotation
        if args.annotation_mode == "snpeff":
            annotated_vcf, html, csv = run_snpeff_annotation(
                filtered_vcf, sample, args.output_dir, config, args.database_name, args.snpeff_java
            )
            create_annotation_tsv(annotated_vcf, sample, args.output_dir, config)
        else:
            from virus_pipeline.annotate_from_config import annotate_from_config
            annotate_from_config(filtered_vcf, args.reference_fasta, config, sample, args.output_dir)

        logging.info(f"Sample complete: {sample}")


if __name__ == "__main__":
    main()
