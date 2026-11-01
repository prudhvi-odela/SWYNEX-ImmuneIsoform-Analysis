# SWYNEX ImmuneIsoform Analysis

## Project objective
Develop a reproducible bioinformatics workflow to investigate how immune activation affects transcript isoform abundance in human monocytes.

## Dataset
- GEO accession: GSE297994
- Organism: Homo sapiens
- Biological material: Human monocytes
- Conditions: Ex vivo and LPS-treated (3 hours)
- Technologies: Short-read RNA-seq and native long-read RNA-seq
- Samples: 22 total; 16 short-read and 6 long-read

## Completed
- Created project directory structure.
- Downloaded and inspected GEO SOFT metadata.
- Identified four processed expression tables.
- Downloaded the four processed tables and passed gzip integrity checks.
- Validated table column counts, row widths, and identifier uniqueness where applicable.
- Extracted sample metadata into metadata/sample_metadata.tsv.
- Verified the expected sample counts for each technology and condition.
- Initialized the local Git repository on branch main.

## Data policy
- Prefer processed data for the initial workflow.
- Do not download raw sequencing archives by default.
- Preserve source metadata and document processing decisions.
- Investigate repeated gene identifiers in the long-read raw-count table before analysis; do not remove duplicates blindly.

## Next steps
- Inspect transcript and gene identifier relationships.
- Define suitable analysis methods for each sequencing technology.
- Build and validate the analysis workflow.
- Document results and limitations before public release.

## Task 2 — Data cleaning and basic analysis (2026-10-24)
- Ran `scripts/clean_and_summarize.py` against all four checked-in processed GEO tables using Python 3.12.14.
- Wrote whitespace-normalized, uncompressed copies to `data/cleaned/` and descriptive quality metrics to `results/task2_data_summary.json`; the original compressed files in `data/processed/` were unchanged.
- All tables had consistent row widths, no blank rows or empty cells, and no nonnumeric or negative values in quantitative columns.
- The long-read raw-count table has repeated `gene_id` annotations across transcript rows. Its separate `ids` column was used as the row-level identifier; no rows were dropped or aggregated.
- No normalization, differential testing, or biological interpretation was performed.

## Dataset source
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297994

## Task 3 — Visualization and interpretation (2026-11-01)
- Inspected and used the Task 2 cleaned short-read and native long-read raw-count tables, the existing `metadata/sample_metadata.tsv`, and `results/task2_data_summary.json` as validation context. Salmon abundance tables were kept separate from raw counts and were not plotted.
- Added `scripts/visualize_task3.py`, a Python 3.9+ standard-library script that maps sample columns to existing metadata and computes raw-count totals and the number of nonzero rows per sample.
- Added separate SVG sample-summary figures for short-read (16 samples) and native long-read (6 samples), plus the machine-readable per-sample metrics TSV.
- The figures show descriptive raw-count totals and nonzero-row counts only. No normalization, statistical testing, differential expression, isoform switching analysis, or biological conclusions were performed.
- Task 3 methods, observed ranges, interpretation boundaries, and limitations are documented in `docs/task3_visualization_interpretation.md`.
- Ran `scripts/visualize_task3.py` with Python 3.12.14 and regenerated the SVG figures and per-sample TSV from the checked-in tables and metadata.
- Revised Task 3 visuals after review: grouped condition bars, consistent condition colors, readable donor labels, and tidy axis intervals; long-read samples are condition-grouped.
- Added a paired short-read top-100 gene-count share figure. In the source matrix, the per-sample top-100 share is 24.9–32.0% ex vivo and 42.6–50.8% after LPS; all eight donor pairs are higher in LPS. This is a descriptive composition pattern, not a differential-expression result, and it flags a normalization consideration for later analysis.
- Final figure polish: added an explicit condition legend to the paired top-100 chart, moved value labels clear of donor-pair connector lines, and placed larger donor labels with colored condition captions beneath both sample-summary panels.
- Extended Task 3 with an expressed-feature CPM versus median-of-ratios boxplot, PCA on the 500 most variable CPM-filtered features, and Pearson correlation heatmaps on all CPM-filtered features. The current Task 3 report presents four figures; the earlier sample-summary-only figures were retired.
- Verified the reported filter sizes and PCA variance from the cleaned tables: 12,428 short-read genes (PC1 95.45%) and 23,386 long-read transcript rows (PC1 94.86%). The short-read filtered-expressed CPM medians range 22.2–24.4 ex vivo versus 7.6–10.0 LPS; median normalized log2 values range 6.70–7.00. Long-read HC3 ex vivo correlations with HC1 and HC2 are 0.725 and 0.785; HC1–HC2 is 0.901. These are descriptive only.
