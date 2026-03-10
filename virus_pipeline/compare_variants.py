#!/usr/bin/env python3
"""Compare Pass 1 (RefSeq) and Pass 2 (self-reference) variants to classify them."""
import argparse
import csv
import os
import sys
from collections import defaultdict


def parse_annotation_tsv(filepath):
    """Parse Pass 1 annotation TSV. Returns list of variant dicts."""
    variants = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            variants.append({
                'CHROM': row.get('CHROM', ''),
                'POS': int(row.get('POS', 0)),
                'REF': row.get('REF', ''),
                'ALT': row.get('ALT', ''),
                'AF': float(row.get('Allele_Frequency', 0) or 0),
                'DP': int(row.get('Total_Depth', 0) or 0),
                'EFFECT': row.get('EFFECT', ''),
                'PROTEIN': row.get('GENE_NAME', ''),
                'HGVSp': row.get('HGVSp', ''),
                'IMPACT': row.get('PUTATIVE_IMPACT', ''),
            })
    return variants


def parse_ivar_tsv(filepath):
    """Parse ivar variants TSV. Returns list of variant dicts."""
    variants = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('PASS', '') != 'TRUE':
                continue
            variants.append({
                'POS': int(row.get('POS', 0)),
                'REF': row.get('REF', ''),
                'ALT': row.get('ALT', ''),
                'AF': float(row.get('ALT_FREQ', 0) or 0),
                'DP': int(row.get('TOTAL_DP', 0) or 0),
            })
    return variants


def main():
    parser = argparse.ArgumentParser(description='Compare Pass 1 and Pass 2 variants')
    parser.add_argument('--pass1_base', required=True, help='Base directory containing pass1_denv1/, pass1_denv2/, pass1_denv3/ output dirs')
    parser.add_argument('--pass2_dir', required=True, help='Pass 2 output directory')
    parser.add_argument('--sample_list', required=True, help='Sample list TSV')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Read sample list
    samples = []
    with open(args.sample_list) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('SAMPLE'):
                continue
            parts = line.split('\t')
            samples.append({
                'name': parts[0],
                'serotype': parts[1],
            })

    # Group samples by serotype for shared-variant analysis
    serotype_samples = defaultdict(list)
    for s in samples:
        serotype_samples[s['serotype']].append(s['name'])

    # Load all Pass 1 variants to build the shared-variant matrix
    pass1_data = {}  # sample -> list of variants
    for s in samples:
        ann_file = os.path.join(args.pass1_base, f"pass1_{s['serotype']}", 'output', f"{s['name']}_annotations.tsv")
        if os.path.exists(ann_file):
            pass1_data[s['name']] = parse_annotation_tsv(ann_file)
        else:
            print(f"WARNING: No Pass 1 annotations for {s['name']}: {ann_file}", file=sys.stderr)
            pass1_data[s['name']] = []

    # Build shared-variant counts per serotype
    # Key: (serotype, pos, ref, alt) -> set of samples
    shared_counts = defaultdict(set)
    for s in samples:
        for v in pass1_data[s['name']]:
            key = (s['serotype'], v['POS'], v['REF'], v['ALT'])
            shared_counts[key].add(s['name'])

    # Process each sample
    summary_rows = []
    for s in samples:
        sample = s['name']
        serotype = s['serotype']
        n_serotype_samples = len(serotype_samples[serotype])

        # Load Pass 2 variants
        ivar_file = os.path.join(args.pass2_dir, sample, f"{sample}_ivar_variants.tsv")
        pass2_variants = []
        if os.path.exists(ivar_file):
            pass2_variants = parse_ivar_tsv(ivar_file)

        # Index Pass 2 by position for lookup
        pass2_by_pos = {}
        for v in pass2_variants:
            pass2_by_pos[(v['POS'], v['REF'], v['ALT'])] = v

        # Classify each Pass 1 variant
        classified = []
        counts = {'FIXED_SHARED': 0, 'FIXED_UNIQUE': 0, 'INTRAHOST': 0, 'MIXED': 0}

        for v in pass1_data[sample]:
            key = (serotype, v['POS'], v['REF'], v['ALT'])
            n_shared = len(shared_counts.get(key, set()))
            pct_shared = n_shared / n_serotype_samples if n_serotype_samples > 0 else 0

            p2_key = (v['POS'], v['REF'], v['ALT'])
            p2 = pass2_by_pos.pop(p2_key, None)

            if v['AF'] >= 0.95 and pct_shared > 0.5:
                classification = 'FIXED_SHARED'
            elif v['AF'] >= 0.95:
                classification = 'FIXED_UNIQUE'
            elif v['AF'] < 0.95 and p2:
                classification = 'MIXED'
            else:
                classification = 'INTRAHOST'

            counts[classification] += 1
            classified.append({
                'CHROM': v['CHROM'],
                'POS': v['POS'],
                'REF': v['REF'],
                'ALT': v['ALT'],
                'PASS1_AF': round(v['AF'], 4),
                'PASS1_DEPTH': v['DP'],
                'PASS2_AF': round(p2['AF'], 4) if p2 else '',
                'PASS2_DEPTH': p2['DP'] if p2 else '',
                'EFFECT': v['EFFECT'],
                'PROTEIN': v['PROTEIN'],
                'HGVSp': v['HGVSp'],
                'CLASS': classification,
            })

        # Add Pass 2-only variants (not seen in Pass 1)
        for p2_key, p2 in pass2_by_pos.items():
            counts['INTRAHOST'] += 1
            classified.append({
                'CHROM': '',
                'POS': p2['POS'],
                'REF': p2['REF'],
                'ALT': p2['ALT'],
                'PASS1_AF': '',
                'PASS1_DEPTH': '',
                'PASS2_AF': round(p2['AF'], 4),
                'PASS2_DEPTH': p2['DP'],
                'EFFECT': '',
                'PROTEIN': '',
                'HGVSp': '',
                'CLASS': 'INTRAHOST',
            })

        # Write per-sample classification
        out_file = os.path.join(args.output_dir, f"{sample}_variant_classification.tsv")
        header = ['CHROM', 'POS', 'REF', 'ALT', 'PASS1_AF', 'PASS1_DEPTH',
                  'PASS2_AF', 'PASS2_DEPTH', 'EFFECT', 'PROTEIN', 'HGVSp', 'CLASS']
        with open(out_file, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for row in sorted(classified, key=lambda x: x['POS']):
                f.write('\t'.join(str(row[h]) for h in header) + '\n')

        total = sum(counts.values())
        summary_rows.append({
            'SAMPLE': sample,
            'SEROTYPE': serotype,
            'TOTAL_PASS1': len(pass1_data[sample]),
            'FIXED_SHARED': counts['FIXED_SHARED'],
            'FIXED_UNIQUE': counts['FIXED_UNIQUE'],
            'INTRAHOST': counts['INTRAHOST'],
            'MIXED': counts['MIXED'],
        })

        print(f"{sample}: {total} variants "
              f"(shared={counts['FIXED_SHARED']}, unique={counts['FIXED_UNIQUE']}, "
              f"intrahost={counts['INTRAHOST']}, mixed={counts['MIXED']})")

    # Write summary
    summary_file = os.path.join(args.output_dir, 'variant_summary.tsv')
    summary_header = ['SAMPLE', 'SEROTYPE', 'TOTAL_PASS1', 'FIXED_SHARED',
                      'FIXED_UNIQUE', 'INTRAHOST', 'MIXED']
    with open(summary_file, 'w') as f:
        f.write('\t'.join(summary_header) + '\n')
        for row in summary_rows:
            f.write('\t'.join(str(row[h]) for h in summary_header) + '\n')

    print(f"\nSummary written to {summary_file}")
    print(f"Per-sample reports in {args.output_dir}/")


if __name__ == '__main__':
    main()
