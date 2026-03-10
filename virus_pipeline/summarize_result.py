import sys
import os
import json
import argparse
import pandas as pd
import logging
from openpyxl import load_workbook
from openpyxl.workbook.views import BookView
pd.set_option('future.no_silent_downcasting', True)

def fix_excel_window(filepath):
    """Set reasonable Excel window dimensions."""
    try:
        wb = load_workbook(filepath)
        wb.views = [BookView(windowWidth=19200, windowHeight=12000)]
        wb.save(filepath)
    except Exception as e:
        logging.warning(f'Could not fix Excel window for {filepath}: {e}')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def summarize_fastp(input_dir, output_file):
    """Parse fastp JSON files for QC metrics."""
    sample_names = []
    reads_before = []
    reads_after = []
    q30_before = []
    q30_after = []
    duplication_rates = []
    adapter_pcts = []

    # Look for fastp JSON in sample output folders
    for item in sorted(os.listdir(input_dir)):
        item_path = os.path.join(input_dir, item)
        if os.path.isdir(item_path) and item.endswith('_output'):
            for fname in os.listdir(item_path):
                if fname.endswith('_fastp.json'):
                    json_path = os.path.join(item_path, fname)
                    sample_name = fname.replace('_fastp.json', '')
                    try:
                        with open(json_path, 'r') as f:
                            data = json.load(f)

                        before = data['summary']['before_filtering']
                        after = data['summary']['after_filtering']
                        dup = data.get('duplication', {}).get('rate', 0.0)
                        adapter = data.get('adapter_cutting', {})
                        adapter_pct = 0.0
                        if adapter and before['total_reads'] > 0:
                            adapter_trimmed = adapter.get('adapter_trimmed_reads', 0)
                            adapter_pct = round(adapter_trimmed / before['total_reads'] * 100, 2)

                        sample_names.append(sample_name)
                        reads_before.append(before['total_reads'])
                        reads_after.append(after['total_reads'])
                        q30_before.append(round(before.get('q30_rate', 0) * 100, 2))
                        q30_after.append(round(after.get('q30_rate', 0) * 100, 2))
                        duplication_rates.append(round(dup * 100, 2))
                        adapter_pcts.append(adapter_pct)
                    except Exception as e:
                        logging.error(f"Error parsing fastp JSON {json_path}: {e}")

    if sample_names:
        df = pd.DataFrame({
            'Sample': sample_names,
            'Reads_Before_Filtering': reads_before,
            'Reads_After_Filtering': reads_after,
            'Pct_Q30_Before': q30_before,
            'Pct_Q30_After': q30_after,
            'Pct_Duplication': duplication_rates,
            'Pct_Adapter': adapter_pcts,
        })
        df.to_excel(output_file, index=False)
        fix_excel_window(output_file)
        logging.info(f"Fastp QC summary saved to {output_file}")
    else:
        logging.warning("No fastp JSON files found")


def summarize_flagstat(input_dir, output_file):
    """Parse samtools flagstat output files."""
    sample_names = []
    total_reads_list = []
    mapped_reads_list = []
    pct_mapped_list = []
    properly_paired_list = []
    pct_properly_paired_list = []

    for fname in sorted(os.listdir(input_dir)):
        if fname.endswith('_flagstat.txt'):
            flagstat_path = os.path.join(input_dir, fname)
            sample_name = fname.replace('_flagstat.txt', '')
            try:
                total_reads = 0
                mapped_reads = 0
                properly_paired = 0
                with open(flagstat_path, 'r') as f:
                    for line in f:
                        if 'in total' in line:
                            total_reads = int(line.split()[0])
                        elif 'mapped (' in line and 'primary mapped' not in line:
                            mapped_reads = int(line.split()[0])
                        elif 'properly paired' in line:
                            properly_paired = int(line.split()[0])

                pct_mapped = round(mapped_reads / total_reads * 100, 2) if total_reads > 0 else 0
                pct_pp = round(properly_paired / total_reads * 100, 2) if total_reads > 0 else 0

                sample_names.append(sample_name)
                total_reads_list.append(total_reads)
                mapped_reads_list.append(mapped_reads)
                pct_mapped_list.append(pct_mapped)
                properly_paired_list.append(properly_paired)
                pct_properly_paired_list.append(pct_pp)
            except Exception as e:
                logging.error(f"Error parsing flagstat {flagstat_path}: {e}")

    if sample_names:
        df = pd.DataFrame({
            'Sample': sample_names,
            'Total_Reads': total_reads_list,
            'Mapped_Reads': mapped_reads_list,
            'Pct_Mapped': pct_mapped_list,
            'Properly_Paired': properly_paired_list,
            'Pct_Properly_Paired': pct_properly_paired_list,
        })
        df.to_excel(output_file, index=False)
        fix_excel_window(output_file)
        logging.info(f"Flagstat summary saved to {output_file}")
    else:
        logging.warning("No flagstat files found")



def summarize_dedup(input_dir, output_file):
    """Parse dedup stats files for deduplication metrics."""
    sample_names = []
    pre_dedup_list = []
    post_dedup_list = []
    pct_dup_list = []

    for fname in sorted(os.listdir(input_dir)):
        if fname.endswith('_dedup_stats.txt'):
            stats_path = os.path.join(input_dir, fname)
            sample_name = fname.replace('_dedup_stats.txt', '')
            try:
                stats = {}
                with open(stats_path, 'r') as f:
                    for line in f:
                        key, val = line.strip().split('\t')
                        stats[key] = val

                sample_names.append(sample_name)
                pre_dedup_list.append(int(stats.get('pre_dedup_reads', 0)))
                post_dedup_list.append(int(stats.get('post_dedup_reads', 0)))
                pct_dup_list.append(float(stats.get('pct_duplicates', 0)))
            except Exception as e:
                logging.error(f"Error parsing dedup stats {stats_path}: {e}")

    if sample_names:
        df = pd.DataFrame({
            'Sample': sample_names,
            'Pre_Dedup_Reads': pre_dedup_list,
            'Post_Dedup_Reads': post_dedup_list,
            'Pct_Duplicates': pct_dup_list,
        })
        df.to_excel(output_file, index=False)
        fix_excel_window(output_file)
        logging.info(f"Dedup summary saved to {output_file}")
    else:
        logging.warning("No dedup stats files found")

def summarize_coverage(input_dir, output_file):
    """Fix 8: Comprehensive coverage metrics instead of average-only."""
    import numpy as np

    sample_names = []
    avg_cov = []
    median_cov = []
    min_cov = []
    max_cov = []
    std_cov = []
    pct_1x = []
    pct_10x = []
    pct_20x = []
    pct_100x = []
    evenness = []

    for file_name in sorted(os.listdir(input_dir)):
        if file_name.endswith('_coverage.txt'):
            coverage_file = os.path.join(input_dir, file_name)
            sample_name = file_name.replace('_coverage.txt', '')
            sample_names.append(sample_name)

            coverage_values = []
            try:
                with open(coverage_file, 'r') as file:
                    for line in file:
                        columns = line.strip().split('\t')
                        depth = int(columns[2])
                        coverage_values.append(depth)
                if coverage_values:
                    arr = np.array(coverage_values)
                    total = len(arr)
                    avg_cov.append(round(float(np.mean(arr)), 2))
                    median_cov.append(round(float(np.median(arr)), 2))
                    min_cov.append(int(np.min(arr)))
                    max_cov.append(int(np.max(arr)))
                    std_cov.append(round(float(np.std(arr)), 2))
                    pct_1x.append(round(float(np.sum(arr >= 1)) / total * 100, 2))
                    pct_10x.append(round(float(np.sum(arr >= 10)) / total * 100, 2))
                    pct_20x.append(round(float(np.sum(arr >= 20)) / total * 100, 2))
                    pct_100x.append(round(float(np.sum(arr >= 100)) / total * 100, 2))
                    # Evenness: 1 - Gini coefficient (1=perfectly even, 0=all reads at one position)
                    sorted_arr = np.sort(arr)
                    n = len(sorted_arr)
                    cumsum = np.cumsum(sorted_arr)
                    gini = (2.0 * np.sum((np.arange(1, n+1) * sorted_arr)) / (n * np.sum(sorted_arr))) - (n + 1) / n if np.sum(sorted_arr) > 0 else 0
                    evenness.append(round(1.0 - gini, 4))
                else:
                    for lst in [avg_cov, median_cov, min_cov, max_cov, std_cov,
                                pct_1x, pct_10x, pct_20x, pct_100x]:
                        lst.append(0)
            except Exception as e:
                logging.error(f"Error processing coverage file {coverage_file}: {e}")
                for lst in [avg_cov, median_cov, min_cov, max_cov, std_cov,
                            pct_1x, pct_10x, pct_20x, pct_100x]:
                    lst.append(0)

    coverage_data = {
        'Sample': sample_names,
        'Average_Coverage': avg_cov,
        'Median_Coverage': median_cov,
        'Min_Coverage': min_cov,
        'Max_Coverage': max_cov,
        'Std_Coverage': std_cov,
        'Pct_Genome_1x': pct_1x,
        'Pct_Genome_10x': pct_10x,
        'Pct_Genome_20x': pct_20x,
        'Pct_Genome_100x': pct_100x,
        'Coverage_Evenness': evenness,
    }
    coverage_df = pd.DataFrame(coverage_data)
    coverage_df.to_excel(output_file, index=False)
    fix_excel_window(output_file)
    logging.info(f"Coverage summary saved to {output_file}")

def summarize_fasta(input_dir, output_file, database_name):
    sample_names = []
    n_counts = []
    total_bases = []

    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fa') and not file_name.endswith('.fasta'):
            fasta_file = os.path.join(input_dir, file_name)
            sample_name = file_name.replace('.fa', '')
            sample_names.append(sample_name)

            n_count = 0
            total_base_count = 0
            try:
                with open(fasta_file, 'r') as file:
                    for line in file:
                        if not line.startswith('>'):
                            sequence = line.strip()
                            n_count += sequence.count('N')
                            total_base_count += len(sequence)
                n_counts.append(n_count)
                total_bases.append(total_base_count)
            except Exception as e:
                logging.error(f"Error processing FASTA file {fasta_file}: {e}")
                n_counts.append(0)
                total_bases.append(0)

    fasta_data = {
        'Sample': sample_names,
        'Number of Ns': n_counts,
        'Total Bases': total_bases
    }
    fasta_df = pd.DataFrame(fasta_data)
    fasta_df.to_excel(output_file, index=False)
    fix_excel_window(output_file)
    logging.info(f"FASTA summary saved to {output_file}")

def merge_excel_files(coverage_file, fasta_file, output_file):
    try:
        coverage_df = pd.read_excel(coverage_file)
        fasta_df = pd.read_excel(fasta_file)
        merged_df = pd.merge(coverage_df, fasta_df, on='Sample', how='outer')
        merged_df.fillna(0, inplace=True)
        merged_df.to_excel(output_file, index=False)
        logging.info(f"Merged summary saved to {output_file}")
    except Exception as e:
        logging.error(f"Error merging Excel files: {e}")
        sys.exit(1)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, help='Path to input directory containing coverage and FASTA files.')
    parser.add_argument('--output_dir', type=str, help='Path to output directory for summarized information.')
    parser.add_argument('--database_name', type=str, default='denv1', help='Name of the SnpEff database (default: denv1).')
    args = parser.parse_args(argv)

    input_dir = args.input_dir
    output_dir = args.output_dir
    database_name = args.database_name

    if not os.path.exists(input_dir):
        logging.error(f"Input directory {input_dir} does not exist")
        sys.exit(1)

    # Fastp QC summary
    fastp_output_file = os.path.join(output_dir, 'qc_summary.xlsx')
    summarize_fastp(input_dir, fastp_output_file)

    # Flagstat summary
    flagstat_output_file = os.path.join(output_dir, 'on_target_summary.xlsx')
    summarize_flagstat(input_dir, flagstat_output_file)

    # Dedup summary
    dedup_output_file = os.path.join(output_dir, 'dedup_summary.xlsx')
    summarize_dedup(input_dir, dedup_output_file)

    # Coverage summary
    coverage_output_file = os.path.join(output_dir, 'coverage_summary.xlsx')
    summarize_coverage(input_dir, coverage_output_file)

    # FASTA summary
    fasta_output_file = os.path.join(output_dir, 'fasta_summary.xlsx')
    summarize_fasta(input_dir, fasta_output_file, database_name)

    # Merge all summaries
    merged_output_file = os.path.join(output_dir, 'merged_summary.xlsx')
    dfs = []
    for f in [fastp_output_file, flagstat_output_file, dedup_output_file, coverage_output_file, fasta_output_file]:
        if os.path.exists(f):
            try:
                dfs.append(pd.read_excel(f))
            except Exception as e:
                logging.warning(f"Could not read {f}: {e}")
    if dfs:
        merged = dfs[0]
        for df in dfs[1:]:
            merged = pd.merge(merged, df, on='Sample', how='outer')
        merged.fillna(0, inplace=True)
        merged.to_excel(merged_output_file, index=False)
        fix_excel_window(merged_output_file)
        logging.info(f"Merged summary saved to {merged_output_file}")
    else:
        logging.error("No summary files to merge")

if __name__ == '__main__':
    main()
