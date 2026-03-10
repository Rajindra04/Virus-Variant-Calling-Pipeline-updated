#!/usr/bin/env python3
"""Read an approved config TSV and generate a pipeline YAML config."""
import argparse
import sys


def parse_config_tsv(filepath):
    """Parse the [GENERAL] and [TRANSCRIPT_ANNOTATIONS] sections from a TSV."""
    general = {}
    annotations = {}  # feature_id -> protein_name
    gene_coords = []  # list of dicts with full annotation info
    current_section = None
    ann_header = None

    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            if line.startswith("[GENERAL]"):
                current_section = "general"
                continue
            elif line.startswith("[TRANSCRIPT_ANNOTATIONS]"):
                current_section = "annotations"
                continue
            elif line.startswith("["):
                current_section = None
                continue

            parts = line.split("\t")

            if current_section == "general":
                if parts[0] == "Field":
                    continue
                if len(parts) >= 2:
                    general[parts[0]] = parts[1]

            elif current_section == "annotations":
                if parts[0] == "feature_id":
                    ann_header = parts
                    continue
                if len(parts) >= 2:
                    annotations[parts[0]] = parts[1]
                    # Also capture full row for gene_coordinates
                    if ann_header:
                        row = dict(zip(ann_header, parts))
                        location = row.get("location", "")
                        if ".." in location:
                            start_str, end_str = location.split("..")
                            gene_coords.append({
                                "feature_id": row.get("feature_id", ""),
                                "protein_name": row.get("protein_name", ""),
                                "feature_type": row.get("feature_type", ""),
                                "start": int(start_str),
                                "end": int(end_str),
                            })

    return general, annotations, gene_coords


def build_yaml(general, annotations, gene_coords=None):
    """Generate YAML config string matching the pipeline's expected structure."""
    virus_name = general.get("virus_name", "unknown")
    genome_type = general.get("genome_type", "ssRNA(+)")
    ploidy = int(general.get("ploidy", 1))
    genome_size = int(general.get("expected_genome_size", 0))
    database_name = general.get("database_name", "custom_db")

    lines = []
    lines.append(f'virus_name: "{virus_name}"')
    lines.append(f'genome_type: "{genome_type}"')
    lines.append(f'ploidy: {ploidy}')
    lines.append(f'expected_genome_size: {genome_size}')
    lines.append(f'database_name: "{database_name}"')
    lines.append("")

    lines.append("fastp:")
    lines.append("  qualified_quality_phred: 20")
    lines.append("  length_required: 50")
    lines.append("  cut_front: true")
    lines.append("  cut_tail: true")
    lines.append("  cut_window_size: 4")
    lines.append("  cut_mean_quality: 20")
    lines.append("  detect_adapter_for_pe: true")
    lines.append("  correction: true")
    lines.append("  overlap_len_require: 30")
    lines.append("  threads: 4")
    lines.append("")

    lines.append("alignment_filtering:")
    lines.append("  min_mapping_quality: 20")
    lines.append('  exclude_flags: "0x904"')
    lines.append("")

    lines.append("deduplication:")
    lines.append("  enabled: true")
    lines.append("  remove_duplicates: true")
    lines.append("")

    lines.append("primer_trimming:")
    lines.append("  min_quality: 20")
    lines.append("  min_length: 50")
    lines.append("  sliding_window: 4")
    lines.append("  include_reads_no_primer: true")
    lines.append("")

    lines.append("coverage:")
    lines.append("  min_base_quality: 20")
    lines.append("  min_mapping_quality: 20")
    lines.append("")

    lines.append("consensus:")
    lines.append("  mpileup_max_depth: 10000")
    lines.append("  mpileup_min_base_quality: 0")
    lines.append("  mpileup_min_mapping_quality: 0")
    lines.append("  ivar_min_quality: 20")
    lines.append("  ivar_min_frequency: 0.5")
    lines.append("  ivar_min_depth: 20")
    lines.append('  ivar_ambiguous_char: "N"')
    lines.append("")

    lines.append("variant_calling:")
    lines.append(f"  ploidy: {ploidy}")
    lines.append("  standard_min_confidence: 30")
    lines.append("  min_base_quality_score: 20")
    lines.append("")

    lines.append("vcf_filtering:")
    lines.append("  filters:")
    lines.append('    LowQD: "QD < 2.0"')
    lines.append('    StrandBias: "FS > 60.0"')
    lines.append('    LowMQ: "MQ < 40.0"')
    lines.append('    LowDepth: "DP < 20"')
    lines.append("  select_pass_only: true")
    lines.append("")

    lines.append("annotation:")
    lines.append('  snpeff_memory: "4g"')
    lines.append("")

    lines.append("transcript_annotations:")
    for feat_id, prot_name in annotations.items():
        lines.append(f'  "{feat_id}": "{prot_name}"')

    if gene_coords:
        lines.append("")
        lines.append("gene_coordinates:")
        for gc in gene_coords:
            pname = gc["protein_name"]
            lines.append(f"  {pname}:")
            lines.append(f"    start: {gc['start']}")
            lines.append(f"    end: {gc['end']}")
            lines.append(f'    feature_id: "{gc["feature_id"]}"')
            lines.append(f'    feature_type: "{gc["feature_type"]}"')

    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description="Build pipeline YAML from approved config TSV")
    parser.add_argument("--tsv", required=True, help="Input config TSV file")
    parser.add_argument("--output", required=True, help="Output YAML config file")
    args = parser.parse_args()

    general, annotations, gene_coords = parse_config_tsv(args.tsv)

    if not general.get("database_name"):
        print("ERROR: No database_name found in TSV", file=sys.stderr)
        sys.exit(1)

    yaml_content = build_yaml(general, annotations, gene_coords)

    with open(args.output, "w") as f:
        f.write(yaml_content)

    print(f"YAML config written to {args.output}")
    print(f"  Virus: {general.get('virus_name', 'unknown')}")
    print(f"  Database: {general.get('database_name', 'unknown')}")
    print(f"  Transcript annotations: {len(annotations)}")


if __name__ == "__main__":
    main()
