#!/usr/bin/env python3
"""Lightweight annotation summary (replaces summarize_snpEff when using config mode)."""
import argparse
import csv
import json
import logging
import os
import sys
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

EFFECT_COLUMNS = [
    'missense_variant', 'synonymous_variant', 'frameshift_variant',
    'stop_gained', 'inframe_insertion', 'inframe_deletion',
    'upstream_gene_variant', 'downstream_gene_variant',
]


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(description='Summarize config-based annotation TSVs')
    parser.add_argument('--input_dir', required=True, help='Directory with *_annotations.tsv files')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args(argv)

    os.makedirs(args.output_dir, exist_ok=True)

    # Find all annotation TSV files
    ann_files = sorted([
        f for f in os.listdir(args.input_dir)
        if f.endswith('_annotations.tsv')
    ])

    if not ann_files:
        logging.warning("No *_annotations.tsv files found")
        return

    # Aggregate: per-protein counts of each effect type across all samples
    # Also per-sample breakdown
    protein_effect_counts = defaultdict(lambda: defaultdict(int))
    sample_data = []

    for ann_file in ann_files:
        sample_name = ann_file.replace('_annotations.tsv', '')
        filepath = os.path.join(args.input_dir, ann_file)

        sample_effects = defaultdict(lambda: defaultdict(int))
        total_variants = 0

        with open(filepath) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                effect = row.get('EFFECT', '')
                protein = row.get('GENE_NAME', '')
                if not effect:
                    continue
                total_variants += 1
                protein_effect_counts[protein][effect] += 1
                sample_effects[protein][effect] += 1

        sample_data.append({
            'sample': sample_name,
            'total_variants': total_variants,
            'effects': dict(sample_effects),
        })

    # Build summary table: one row per protein
    summary_rows = []
    for protein in sorted(protein_effect_counts.keys()):
        counts = protein_effect_counts[protein]
        row = {'Protein': protein, 'Samples': len(ann_files)}
        total = 0
        for eff in EFFECT_COLUMNS:
            row[eff] = counts.get(eff, 0)
            total += counts.get(eff, 0)
        # Add any effects not in the standard list
        for eff, count in counts.items():
            if eff not in EFFECT_COLUMNS:
                row[eff] = count
                total += count
        row['total'] = total
        summary_rows.append(row)

    # Write summary CSV
    output_file = os.path.join(args.output_dir, 'summary_table.csv')
    all_columns = ['Protein', 'Samples'] + EFFECT_COLUMNS + ['total']
    # Add any extra effect columns
    extra_effects = set()
    for row in summary_rows:
        for k in row:
            if k not in all_columns:
                extra_effects.add(k)
    all_columns = all_columns[:-1] + sorted(extra_effects) + ['total']

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=all_columns, extrasaction='ignore')
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    logging.info(f"Summary table saved to {output_file}")

    # Frameshift QC warning
    total_frameshifts = sum(
        protein_effect_counts[p].get('frameshift_variant', 0)
        for p in protein_effect_counts
    )
    if total_frameshifts > 0:
        logging.warning(
            "QC ALERT: %d frameshift variant(s) detected. "
            "For single-polyprotein viruses, frameshifts may indicate "
            "alignment artifacts or assembly errors. Review the annotations TSV.",
            total_frameshifts
        )

    # Write chart data JSON (Chart.js-ready format)
    chart_data = {
        'labels': [row['Protein'] for row in summary_rows],
        'datasets': []
    }
    for eff in EFFECT_COLUMNS:
        chart_data['datasets'].append({
            'label': eff,
            'data': [row.get(eff, 0) for row in summary_rows],
        })

    chart_file = os.path.join(args.output_dir, 'chart_data.json')
    with open(chart_file, 'w') as f:
        json.dump(chart_data, f, indent=2)

    logging.info(f"Chart data saved to {chart_file}")
    logging.info(f"Processed {len(ann_files)} samples, "
                 f"{len(protein_effect_counts)} proteins")


if __name__ == '__main__':
    main()
