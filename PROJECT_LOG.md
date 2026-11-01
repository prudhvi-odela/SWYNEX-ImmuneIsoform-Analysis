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
- Ran `scripts/visualize_task3.py` with Python 3.12.14 and regenerated both SVG figures and the per-sample TSV from the checked-in tables and metadata.
