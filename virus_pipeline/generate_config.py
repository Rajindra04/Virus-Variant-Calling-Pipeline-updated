#!/usr/bin/env python3
"""Parse a GenBank file and generate a TSV config for user review."""
import argparse
import re
import sys


def parse_genbank(filepath):
    """Parse a GenBank file without BioPython."""
    organism = ""
    accession_version = ""
    seq_length = 0
    features = []

    with open(filepath) as f:
        lines = f.readlines()

    # --- Header parsing ---
    features_line_idx = None
    for i, line in enumerate(lines):
        if line.startswith("LOCUS"):
            parts = line.split()
            seq_length = int(parts[2])

        elif line.startswith("VERSION"):
            accession_version = line.split()[1].strip()

        elif line.strip().startswith("ORGANISM"):
            organism = line.strip().replace("ORGANISM", "").strip()

        elif line.startswith("FEATURES"):
            features_line_idx = i
            break

    if features_line_idx is None:
        print("ERROR: No FEATURES section found", file=sys.stderr)
        sys.exit(1)

    # --- Feature parsing ---
    current_feature = None
    last_qual_key = None

    for line in lines[features_line_idx + 1:]:
        if line.startswith("ORIGIN") or line.startswith("//"):
            if current_feature and current_feature["type"] in ("CDS", "mat_peptide"):
                features.append(current_feature)
            break

        # Feature key line: 5 spaces then non-space at column 5
        if len(line) > 5 and line[:5] == "     " and line[5] != " ":
            # Save previous feature
            if current_feature and current_feature["type"] in ("CDS", "mat_peptide"):
                features.append(current_feature)

            stripped = line.strip()
            parts = stripped.split(None, 1)
            if len(parts) == 2:
                feat_type, location = parts
                current_feature = {
                    "type": feat_type,
                    "location": location.strip(),
                    "qualifiers": {}
                }
                last_qual_key = None
            else:
                current_feature = None
                last_qual_key = None

        elif line.startswith(" " * 21) and current_feature:
            stripped = line.strip()
            if stripped.startswith("/"):
                if "=" in stripped:
                    key, val = stripped[1:].split("=", 1)
                    val = val.strip('"')
                    current_feature["qualifiers"][key] = val
                    last_qual_key = key
                else:
                    key = stripped[1:]
                    current_feature["qualifiers"][key] = True
                    last_qual_key = None
            elif last_qual_key:
                # Continuation of multi-line qualifier
                prev = current_feature["qualifiers"].get(last_qual_key, "")
                if isinstance(prev, str):
                    current_feature["qualifiers"][last_qual_key] = prev + stripped.strip('"')

    return organism, accession_version, seq_length, features


def suggest_short_name(product):
    """Suggest a short protein name from GenBank product field."""
    if not product:
        return "unknown"
    product_lower = product.lower()

    if product_lower == "polyprotein":
        return "Polyprotein"

    # Check for NS proteins
    ns_match = re.search(r'NS[1-5][AB]?', product, re.IGNORECASE)
    if ns_match:
        return ns_match.group(0).upper()

    # Check for 2K
    if "2k" in product_lower or "2K" in product:
        return "2K"

    # Specific patterns
    patterns = [
        (r'anchored capsid.*ancC', 'ancC'),
        (r'capsid protein C', 'Capsid'),
        (r'envelope protein E', 'E'),
        (r'membrane glycoprotein precursor.*prM', 'prM'),
        (r'protein pr$', 'pr'),
        (r'membrane glycoprotein M', 'M'),
        (r'RNA-dependent RNA polymerase', 'NS5'),
    ]
    for pattern, name in patterns:
        if re.search(pattern, product, re.IGNORECASE):
            return name

    return product


def main():
    parser = argparse.ArgumentParser(description="Parse GenBank file and generate config TSV")
    parser.add_argument("--genbank", required=True, help="Path to GenBank file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    organism, accession_version, seq_length, features = parse_genbank(args.genbank)

    with open(args.output, "w") as f:
        # General section
        f.write("[GENERAL]\n")
        f.write("Field\tValue\n")
        f.write(f"virus_name\t{organism}\n")
        f.write(f"database_name\t{accession_version}\n")
        f.write(f"expected_genome_size\t{seq_length}\n")
        f.write(f"genome_type\tssRNA(+)\n")
        f.write(f"ploidy\t1\n")
        f.write("\n")

        # Transcript annotations section
        f.write("[TRANSCRIPT_ANNOTATIONS]\n")
        f.write("feature_id\tprotein_name\tfeature_type\tproduct\tlocation\n")

        for feat in features:
            q = feat["qualifiers"]
            product = q.get("product", "")
            protein_id = q.get("protein_id", "")
            locus_tag = q.get("locus_tag", "")

            if feat["type"] == "CDS":
                feature_id = locus_tag if locus_tag else protein_id
            else:  # mat_peptide
                feature_id = protein_id if protein_id else f"{locus_tag}_mp"

            if not feature_id:
                print(f"WARNING: No feature_id for {feat['type']} at {feat['location']}", file=sys.stderr)
                continue

            short_name = suggest_short_name(product)
            f.write(f"{feature_id}\t{short_name}\t{feat['type']}\t{product}\t{feat['location']}\n")

    print(f"Config TSV written to {args.output}")
    print(f"  Organism: {organism}")
    print(f"  Accession: {accession_version}")
    print(f"  Genome size: {seq_length} bp")
    print(f"  Features: {len(features)} (CDS + mat_peptide)")


if __name__ == "__main__":
    main()
