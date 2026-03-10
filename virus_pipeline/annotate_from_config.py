#!/usr/bin/env python3
"""Lightweight VCF annotator using gene coordinates from the pipeline YAML config."""
import argparse
import logging
import os
import sys

import yaml

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

AA_3LETTER = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'Ter',
}


def read_fasta(filepath):
    """Read a single-sequence FASTA file. Returns (header, sequence)."""
    header = None
    seq_parts = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
            else:
                seq_parts.append(line)
    return header, ''.join(seq_parts).upper()


def translate_codon(codon):
    """Translate a 3-letter codon to amino acid."""
    return CODON_TABLE.get(codon.upper(), 'X')


def find_gene(pos, gene_coords):
    """Find the best gene for a genomic position.

    Priority: mat_peptide over CDS, shorter (more specific) over longer.
    If position falls outside all genes, find nearest gene for upstream/downstream.
    """
    mat_peptide_hits = []
    cds_hits = []

    for name, info in gene_coords.items():
        if info['start'] <= pos <= info['end']:
            gene_len = info['end'] - info['start'] + 1
            if info['feature_type'] == 'mat_peptide':
                mat_peptide_hits.append((name, info, gene_len))
            else:
                cds_hits.append((name, info, gene_len))

    # Prefer mat_peptide, then shortest (most specific) gene
    if mat_peptide_hits:
        mat_peptide_hits.sort(key=lambda x: x[2])
        return mat_peptide_hits[0][0], mat_peptide_hits[0][1]
    if cds_hits:
        cds_hits.sort(key=lambda x: x[2])
        return cds_hits[0][0], cds_hits[0][1]

    # Position not in any gene — find the nearest one for upstream/downstream
    nearest_name = None
    nearest_info = None
    nearest_dist = float('inf')
    for name, info in gene_coords.items():
        dist_start = abs(pos - info['start'])
        dist_end = abs(pos - info['end'])
        dist = min(dist_start, dist_end)
        if dist < nearest_dist:
            nearest_dist = dist
            nearest_name = name
            nearest_info = info

    if nearest_info:
        return nearest_name, nearest_info

    return None, None


def classify_variant(pos, ref, alt, gene_name, gene_info, ref_seq):
    """Classify a variant and compute effect, impact, HGVSc, HGVSp."""
    result = {
        'effect': '', 'impact': '', 'hgvsc': '', 'hgvsp': '',
        'cdna_pos': '', 'cds_pos': '', 'protein_pos': '',
        'feature_type': '', 'gene_name': gene_name or '',
    }

    if gene_info is None:
        result['effect'] = 'intergenic_region'
        result['impact'] = 'MODIFIER'
        result['gene_name'] = 'intergenic'
        return result

    gene_start = gene_info['start']
    gene_end = gene_info['end']
    gene_len = gene_end - gene_start + 1
    result['feature_type'] = gene_info['feature_type']

    # Check if position is upstream or downstream
    if pos < gene_start:
        result['effect'] = 'upstream_gene_variant'
        result['impact'] = 'MODIFIER'
        return result
    if pos > gene_end:
        result['effect'] = 'downstream_gene_variant'
        result['impact'] = 'MODIFIER'
        return result

    # Position is within the gene
    # CDS-relative position (1-based)
    cds_pos = pos - gene_start + 1
    protein_len = gene_len // 3

    # Handle indels
    if len(ref) != len(alt):
        indel_len = abs(len(alt) - len(ref))
        if indel_len % 3 == 0:
            result['effect'] = 'inframe_insertion' if len(alt) > len(ref) else 'inframe_deletion'
            result['impact'] = 'MODERATE'
        else:
            result['effect'] = 'frameshift_variant'
            result['impact'] = 'HIGH'
        result['hgvsc'] = f"c.{cds_pos}{'ins' if len(alt) > len(ref) else 'del'}"
        result['cdna_pos'] = f"{cds_pos}/{gene_len}"
        result['cds_pos'] = f"{cds_pos}/{gene_len}"
        codon_num = (cds_pos - 1) // 3 + 1
        result['protein_pos'] = f"{codon_num}/{protein_len}"
        return result

    # SNV — compute codon change
    codon_num = (cds_pos - 1) // 3 + 1
    codon_offset = (cds_pos - 1) % 3

    # Extract reference codon (0-based indexing into ref_seq)
    codon_start_genome = gene_start + (codon_num - 1) * 3 - 1  # 0-based
    if codon_start_genome + 3 > len(ref_seq):
        result['effect'] = 'incomplete_terminal_codon'
        result['impact'] = 'LOW'
        return result

    ref_codon = list(ref_seq[codon_start_genome:codon_start_genome + 3])
    alt_codon = ref_codon.copy()
    alt_codon[codon_offset] = alt[0] if len(alt) == 1 else alt

    ref_codon_str = ''.join(ref_codon)
    alt_codon_str = ''.join(alt_codon)

    ref_aa = translate_codon(ref_codon_str)
    alt_aa = translate_codon(alt_codon_str)

    result['hgvsc'] = f"c.{cds_pos}{ref[0]}>{alt[0]}"
    result['cdna_pos'] = f"{cds_pos}/{gene_len}"
    result['cds_pos'] = f"{cds_pos}/{gene_len}"
    result['protein_pos'] = f"{codon_num}/{protein_len}"

    ref_aa3 = AA_3LETTER.get(ref_aa, ref_aa)
    alt_aa3 = AA_3LETTER.get(alt_aa, alt_aa)

    if ref_aa == alt_aa:
        result['effect'] = 'synonymous_variant'
        result['impact'] = 'LOW'
        result['hgvsp'] = f"p.{ref_aa3}{codon_num}{alt_aa3}"
    elif alt_aa == '*':
        result['effect'] = 'stop_gained'
        result['impact'] = 'HIGH'
        result['hgvsp'] = f"p.{ref_aa3}{codon_num}{alt_aa3}"
    elif ref_aa == '*':
        result['effect'] = 'stop_lost'
        result['impact'] = 'HIGH'
        result['hgvsp'] = f"p.{ref_aa3}{codon_num}{alt_aa3}"
    else:
        result['effect'] = 'missense_variant'
        result['impact'] = 'MODERATE'
        result['hgvsp'] = f"p.{ref_aa3}{codon_num}{alt_aa3}"

    return result


def annotate_from_config(vcf_file, reference, config, sample_name, output_dir):
    """Annotate a filtered VCF using gene coordinates from config YAML."""
    gene_coords = config.get('gene_coordinates', {})
    if not gene_coords:
        logging.error("No gene_coordinates found in config YAML")
        return None

    _, ref_seq = read_fasta(reference)

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
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
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

            # Get AD from FORMAT field
            allelic_depths = ''
            if len(fields) > 9:
                fmt_keys = fields[8].split(':')
                fmt_vals = fields[9].split(':')
                fmt_dict = dict(zip(fmt_keys, fmt_vals))
                allelic_depths = fmt_dict.get('AD', '')

            # Find gene and classify variant
            gene_name, gene_info = find_gene(pos, gene_coords)
            result = classify_variant(pos, ref, alt, gene_name, gene_info, ref_seq)

            rows.append([
                chrom, str(pos), vid, ref, alt, qual, filt,
                total_depth, allele_freq, strand_bias, allelic_depths,
                result['effect'], result['impact'],
                result['gene_name'], result['gene_name'],
                result['feature_type'], 'config_annotation', 'protein_coding',
                result['hgvsc'], result['hgvsp'],
                result['cdna_pos'], result['cds_pos'], result['protein_pos'],
                ''
            ])

    with open(output_file, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')

    logging.info(f"Config-based annotation TSV: {output_file} ({len(rows)} variants)")
    return output_file


def main():
    parser = argparse.ArgumentParser(description='Annotate VCF from config gene coordinates')
    parser.add_argument('--vcf', required=True, help='Filtered VCF file')
    parser.add_argument('--reference', required=True, help='Reference FASTA')
    parser.add_argument('--config', required=True, help='Pipeline YAML config')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    os.makedirs(args.output_dir, exist_ok=True)
    annotate_from_config(args.vcf, args.reference, config, args.sample_name, args.output_dir)


if __name__ == '__main__':
    main()
