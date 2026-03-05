import sys
import pandas as pd
import glob
import os
import json
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_transcript_map(config_path):
    """Load transcript-to-protein mapping from config YAML if available."""
    if not config_path:
        return {}
    try:
        import yaml
        with open(config_path) as f:
            config = yaml.safe_load(f)
        return config.get('transcript_annotations', {})
    except Exception as e:
        logging.warning(f"Could not load transcript map from {config_path}: {e}")
        return {}

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True, help='Path to input directory containing SnpEff summary files.')
    parser.add_argument('--output_dir', type=str, required=True, help='Path to output directory for summary files.')
    parser.add_argument('--config', type=str, default=None, help='Path to virus config YAML (for transcript annotations)')
    args = parser.parse_args(argv)

    input_dir = args.input_dir
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load transcript-to-protein mapping
    transcript_map = load_transcript_map(args.config)

    # File pattern
    file_pattern = "*_snpEff_summary.genes.txt"

    # Expected columns to keep
    columns_to_keep = [
        "TranscriptId",
        "variants_impact_HIGH",
        "variants_impact_LOW",
        "variants_impact_MODERATE",
        "variants_impact_MODIFIER",
        "variants_effect_downstream_gene_variant",
        "variants_effect_frameshift_variant",
        "variants_effect_missense_variant",
        "variants_effect_splice_region_variant",
        "variants_effect_synonymous_variant"
    ]

    # Initialize data list
    all_data = []

    # Check if directory exists
    if not os.path.exists(input_dir):
        logging.error(f"Directory '{input_dir}' does not exist.")
        sys.exit(1)

    # Get list of files
    files = glob.glob(os.path.join(input_dir, file_pattern))
    if not files:
        logging.error(f"No files found matching pattern '{file_pattern}' in '{input_dir}'.")
        logging.error(f"Files in directory: {os.listdir(input_dir)}")
        sys.exit(1)

    # Iterate over files
    for file_path in files:
        sample_name = os.path.basename(file_path).replace("_snpEff_summary.genes.txt", "")
        logging.info(f"Processing file: {file_path} (Sample: {sample_name})")

        try:
            # Read with tab delimiter; line 0 is a comment, line 1 is the header
            df = pd.read_csv(file_path, sep="\t", header=1)

            # Clean column names (strip whitespace, handle leading #)
            df.columns = [c.strip().lstrip('#') for c in df.columns]

            logging.info(f"Columns found in {file_path}: {list(df.columns)} ({len(df.columns)} columns)")

            # Select only the columns we want that actually exist
            available_cols = [col for col in columns_to_keep if col in df.columns]

            if not available_cols:
                logging.error(f"No expected columns found in {file_path}. Skipping.")
                continue

            df = df[available_cols].copy()

            # Fill missing expected columns with 0
            for col in columns_to_keep:
                if col not in df.columns:
                    df[col] = 0

            df["Sample"] = sample_name
            all_data.append(df)
        except Exception as e:
            logging.error(f"Error reading {file_path}: {str(e)}")
            continue

    if not all_data:
        logging.error("No valid dataframes to concatenate. Check file contents or column names.")
        sys.exit(1)

    combined_df = pd.concat(all_data, ignore_index=True)
    summary_df = combined_df.groupby("TranscriptId").agg({
        "variants_impact_HIGH": "sum",
        "variants_impact_LOW": "sum",
        "variants_impact_MODERATE": "sum",
        "variants_impact_MODIFIER": "sum",
        "variants_effect_downstream_gene_variant": "sum",
        "variants_effect_frameshift_variant": "sum",
        "variants_effect_missense_variant": "sum",
        "variants_effect_splice_region_variant": "sum",
        "variants_effect_synonymous_variant": "sum",
        "Sample": lambda x: ";".join(x)
    }).reset_index()

    summary_df.columns = [
        "TranscriptId",
        "Total_variants_impact_HIGH",
        "Total_variants_impact_LOW",
        "Total_variants_impact_MODERATE",
        "Total_variants_impact_MODIFIER",
        "Total_variants_effect_downstream_gene_variant",
        "Total_variants_effect_frameshift_variant",
        "Total_variants_effect_missense_variant",
        "Total_variants_effect_splice_region_variant",
        "Total_variants_effect_synonymous_variant",
        "Samples"
    ]

    # Add Protein column from transcript mapping
    if transcript_map:
        summary_df.insert(1, "Protein", summary_df["TranscriptId"].map(transcript_map).fillna(""))
    else:
        summary_df.insert(1, "Protein", "")

    output_file = os.path.join(output_dir, "summary_table.csv")
    summary_df.to_csv(output_file, index=False)
    logging.info(f"Summary table saved to {output_file}")

    chart_data = {
        "labels": summary_df["TranscriptId"].tolist(),
        "datasets": [
            {
                "label": "High Impact",
                "data": summary_df["Total_variants_impact_HIGH"].tolist(),
                "backgroundColor": "rgba(255, 99, 132, 0.7)"
            },
            {
                "label": "Low Impact",
                "data": summary_df["Total_variants_impact_LOW"].tolist(),
                "backgroundColor": "rgba(54, 162, 235, 0.7)"
            },
            {
                "label": "Moderate Impact",
                "data": summary_df["Total_variants_impact_MODERATE"].tolist(),
                "backgroundColor": "rgba(255, 206, 86, 0.7)"
            },
            {
                "label": "Modifier Impact",
                "data": summary_df["Total_variants_impact_MODIFIER"].tolist(),
                "backgroundColor": "rgba(75, 192, 192, 0.7)"
            }
        ]
    }

    chart_json_file = os.path.join(output_dir, "chart_data.json")
    with open(chart_json_file, "w") as f:
        json.dump(chart_data, f, indent=2)
    logging.info(f"Chart data for external visualization saved to {chart_json_file}")

if __name__ == '__main__':
    main()
