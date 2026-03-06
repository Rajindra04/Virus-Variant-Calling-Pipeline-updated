#!/usr/bin/env python3
"""Replace Ns in consensus FASTA with reference bases."""
import argparse
import sys


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
    return header, ''.join(seq_parts)


def write_fasta(filepath, header, sequence, line_width=80):
    """Write a FASTA file with wrapped sequence lines."""
    with open(filepath, 'w') as f:
        f.write(f'>{header}\n')
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i + line_width] + '\n')


def main():
    parser = argparse.ArgumentParser(description='Fill Ns in consensus with reference bases')
    parser.add_argument('--consensus', required=True, help='Consensus FASTA file')
    parser.add_argument('--reference', required=True, help='Reference FASTA file')
    parser.add_argument('--output', required=True, help='Output filled FASTA file')
    parser.add_argument('--sample_name', required=True, help='Sample name for FASTA header')
    args = parser.parse_args()

    _, cons_seq = read_fasta(args.consensus)
    _, ref_seq = read_fasta(args.reference)

    if len(cons_seq) != len(ref_seq):
        diff = len(ref_seq) - len(cons_seq)
        print(f"WARNING: Consensus length ({len(cons_seq)}) != reference length ({len(ref_seq)}), "
              f"difference={diff}bp", file=sys.stderr)
        if abs(diff) > 10:
            print("ERROR: Length difference >10bp, likely a major assembly issue. Aborting.",
                  file=sys.stderr)
            sys.exit(1)
        # Pad consensus with Ns or truncate to match reference length
        if len(cons_seq) < len(ref_seq):
            cons_seq = cons_seq + 'N' * diff
            print(f"  Padded consensus with {diff} Ns at 3' end", file=sys.stderr)
        else:
            cons_seq = cons_seq[:len(ref_seq)]
            print(f"  Truncated consensus by {-diff}bp at 3' end", file=sys.stderr)

    filled = []
    n_count = 0
    for i in range(len(cons_seq)):
        if cons_seq[i] in ('N', 'n'):
            filled.append(ref_seq[i])
            n_count += 1
        else:
            filled.append(cons_seq[i])

    filled_seq = ''.join(filled)
    write_fasta(args.output, args.sample_name, filled_seq)

    pct = round(n_count / len(cons_seq) * 100, 2)
    print(f"Filled {n_count} of {len(cons_seq)} Ns ({pct}% of genome)")


if __name__ == '__main__':
    main()
