# Pipeline Output Documentation

## Summary Files

### coverage_summary.xlsx
| Column | Definition |
|--------|-----------|
| Sample | User-defined sample name |
| Average_Coverage | Mean read depth across all genome positions |
| Median_Coverage | Median read depth across all genome positions |
| Min_Coverage | Minimum read depth at any position |
| Max_Coverage | Maximum read depth at any position |
| Std_Coverage | Standard deviation of read depth |
| Pct_Genome_1x | % of genome positions with >= 1 read |
| Pct_Genome_10x | % of genome positions with >= 10 reads |
| Pct_Genome_20x | % of genome positions with >= 20 reads (matches ivar consensus threshold) |
| Pct_Genome_100x | % of genome positions with >= 100 reads |
| Coverage_Evenness | 1 - Gini coefficient. Ranges 0-1. 1 = perfectly uniform coverage across all positions. 0 = all reads at one position. Lower values indicate uneven coverage, common with amplicon sequencing. |

### qc_summary.xlsx
| Column | Definition |
|--------|-----------|
| Sample | User-defined sample name |
| Reads_Before_Filtering | Total raw read count from sequencer |
| Reads_After_Filtering | Reads remaining after fastp quality/adapter trimming |
| Pct_Q30_Before | % of bases with Phred score >= 30 before trimming |
| Pct_Q30_After | % of bases with Phred score >= 30 after trimming |
| Pct_Duplication | Estimated duplication rate from fastp (pre-alignment, k-mer frequency heuristic). Always lower than Pct_Duplicates. |
| Pct_Adapter | % of reads with adapter contamination |

### dedup_summary.xlsx
| Column | Definition |
|--------|-----------|
| Sample | User-defined sample name |
| Pre_Dedup_Reads | Aligned reads before duplicate removal (after MAPQ/flag filtering) |
| Post_Dedup_Reads | Aligned reads after samtools markdup removes PCR duplicates |
| Pct_Duplicates | Actual PCR duplicate rate based on mapping position (samtools markdup). Higher than Pct_Duplication because position-based detection is more sensitive than k-mer estimation. |

### on_target_summary.xlsx
| Column | Definition |
|--------|-----------|
| Sample | User-defined sample name |
| Total_Reads | Total reads in the final analysis BAM (post-dedup, post-primer-trim) |
| Mapped_Reads | Reads that aligned to the reference genome |
| Pct_Mapped | % of reads that mapped |
| Properly_Paired | Read pairs where both mates mapped in expected orientation/distance |
| Pct_Properly_Paired | % of reads that are properly paired |

### fasta_summary.xlsx
| Column | Definition |
|--------|-----------|
| Sample | User-defined sample name |
| Number of Ns | Count of ambiguous bases (N) in consensus. Positions with < 20X coverage are called as N by ivar consensus. |
| Total Bases | Total length of the consensus sequence |

### merged_summary.xlsx
Combined view of all the above summaries, joined on Sample name.

## Per-Sample Files

### {sample}_annotations.tsv
User-friendly variant annotation table. One row per variant, showing the most specific
protein annotation (mature peptide, not polyprotein).
| Column | Definition |
|--------|-----------|
| CHROM | Reference chromosome/accession |
| POS | Genomic position |
| ID | Variant ID (usually ".") |
| REF | Reference allele |
| ALT | Alternate allele |
| QUAL | GATK quality score |
| FILTER | PASS or filter reason (LowQD, StrandBias, LowMQ, LowDepth) |
| Total_Depth | Total read depth at this position (INFO DP) |
| Allele_Frequency | Frequency of ALT allele (INFO AF) |
| Strand_Bias | Fisher strand bias score (INFO FS). Higher = more bias. |
| Allelic_Depths | Ref and alt read counts (FORMAT AD, e.g., "5,245") |
| EFFECT | Variant effect (missense_variant, synonymous_variant, etc.) |
| PUTATIVE_IMPACT | Impact level: HIGH, MODERATE, LOW, MODIFIER |
| GENE_NAME | Protein name from transcript annotations config |
| GENE_ID | Gene identifier |
| FEATURE_TYPE | Feature type (transcript) |
| FEATURE_ID | Transcript/protein accession |
| TRANSCRIPT_TYPE | Biotype (protein_coding) |
| HGVSc | HGVS coding DNA notation (e.g., c.5A>T) |
| HGVSp | HGVS protein notation (e.g., p.Asn2Ile) |
| cDNA_POSITION_AND_LENGTH | Position within cDNA / total cDNA length |
| CDS_POSITION_AND_LENGTH | Position within CDS / total CDS length |
| PROTEIN_POSITION_AND_LENGTH | Amino acid position / total protein length |
| ERROR | snpEff warnings (e.g., WARNING_TRANSCRIPT_NO_STOP_CODON) |

### {sample}_low_coverage.tsv
Positions with read depth below the consensus threshold (20X).
| Column | Definition |
|--------|-----------|
| CHROM | Reference chromosome/accession |
| POSITION | Genomic position |
| DEPTH | Read depth at this position |

### {sample}_coverage.png
Coverage plot showing read depth across the genome on a log scale, with a red
dashed line at 20X (the ivar consensus minimum depth threshold).

### {sample}_coverage.txt
Raw per-position coverage (samtools depth output). Tab-separated: CHROM, POS, DEPTH.
Generated from the deduplicated, primer-trimmed BAM.

### {sample}_flagstat.txt
samtools flagstat output from the deduplicated, primer-trimmed BAM.

### {sample}_dedup_stats.txt
Tab-separated deduplication statistics: pre_dedup_reads, post_dedup_reads,
duplicates_removed, pct_duplicates.

### VCF Files (in order of generation)
1. **{sample}.vcf** -- Raw GATK HaplotypeCaller output (all called variants)
2. **{sample}_filtered.vcf** -- After GATK VariantFiltration (filter tags added but all variants present)
3. **{sample}_pass.vcf** -- After GATK SelectVariants --exclude-filtered (PASS variants only)
4. **{sample}_annotated.vcf** -- After snpEff annotation (ANN field added to PASS variants)

### {sample}.fa
Consensus FASTA generated by ivar consensus. Positions with <20X coverage are called as N.

### {sample}_snpEff_summary.html / .csv / .genes.txt
snpEff summary reports. The .genes.txt contains per-gene variant impact counts.

### {sample}_snpSift.txt
Raw SnpSift field extraction (all transcript annotations per variant, wide format).
Use {sample}_annotations.tsv for a cleaner view.

### snpEff.config
Auto-generated snpEff configuration for this run. Database named by accession number.

## Pipeline Processing Order
1. fastp (read trimming + QC) -> qc_summary.xlsx
2. FastQC (quality reports)
3. QC gate (fail/warn based on fastp metrics)
4. bwa-mem2 (read mapping)
5. samtools view (MAPQ >= 20, exclude unmapped/secondary/supplementary)
6. samtools fixmate + sort + markdup -r (deduplication, removes PCR duplicates) -> dedup_stats.txt
7. samtools addreplacerg (add read groups)
8. ivar trim (primer trimming)
9. samtools depth -a (coverage) -> coverage.txt, coverage.png, low_coverage.tsv
10. samtools flagstat -> flagstat.txt, on_target_summary.xlsx
11. ivar consensus -> {sample}.fa, fasta_summary.xlsx
12. GATK HaplotypeCaller -> {sample}.vcf
13. GATK VariantFiltration -> {sample}_filtered.vcf
14. GATK SelectVariants -> {sample}_pass.vcf
15. snpEff annotation -> {sample}_annotated.vcf
16. SnpSift extract -> {sample}_snpSift.txt
17. create_annotation_tsv -> {sample}_annotations.tsv
18. summarize_result -> coverage_summary.xlsx, merged_summary.xlsx
19. summarize_snpEff -> summary_table.csv
