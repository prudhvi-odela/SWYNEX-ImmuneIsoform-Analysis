# Task 2 — Data cleaning and basic analysis

## Scope

This workflow uses the four processed GSE297994 tables in [data/processed/](../data/processed/). It performs structural checks, conservative text cleanup, and descriptive summaries. Short-read and long-read tables, and count and abundance measures, remain separate; no cross-assay merging or statistical comparison is performed.

## Run

From the repository root:

```bash
python scripts/clean_and_summarize.py
```

The script uses Python 3.9+ and its standard library only. It writes normalized TSV copies under `data/cleaned/` and a machine-readable report to `results/task2_data_summary.json`. Original compressed source files remain unchanged.

## Cleaning decisions

- If the first header cell is blank (a common representation of row names), label that identifier column `feature_id`; reject other blank or duplicate column names.
- Decode UTF-8 (including an optional byte-order mark) and parse tab-separated fields.
- Trim surrounding whitespace from headers and cells; omit fully blank rows.
- Preserve missing cells as empty fields; do not impute values.
- Preserve all nonblank rows, including repeated gene annotations; do not drop or aggregate them. For the native long-read raw-count table, use its row-level `ids` column as the feature identifier because `gene_id` repeats across transcript rows.
- Stop if a nonblank row has a different number of fields from the header, or if a header is blank or duplicated.
- For each non-ID column, report nonempty-cell count, numeric-value count, mean, minimum, maximum, zero and negative counts, and nonnumeric-value count. Missing and non-finite values are excluded from numeric summaries.

This is formatting and integrity cleanup only. It does not normalize counts, infer sample metadata, remove outliers, perform differential testing, or establish biological findings.

## Results

Run completed on 2026-10-24 with Python 3.12.14. All four source tables had consistent row widths, zero blank rows, zero empty cells, zero missing feature IDs, zero duplicate row-level feature IDs, and zero nonnumeric or negative quantitative values. The cleaned tables and full per-column summary are in `data/cleaned/` and `results/task2_data_summary.json`.

| Processed table | Rows | Columns | Numeric sample columns | Numeric cells | Zero values | Maximum |
|---|---:|---:|---:|---:|---:|---:|
| Native long-read raw counts | 24,894 | 8 | 6 | 149,364 | 23,785 | 37,572 |
| Native long-read Salmon abundance | 24,939 | 7 | 6 | 149,634 | 12,613 | 31,261.6322 |
| Short-read raw counts | 63,140 | 17 | 16 | 1,010,240 | 712,389 | 553,043 |
| Short-read Salmon abundance | 24,939 | 17 | 16 | 399,024 | 184,585 | 217,117.4459 |

The long-read raw-count table contains 16,439 repeated values in its first (`gene_id`) column, while the separate `ids` row identifier is unique. These are repeated gene annotations across transcript rows; the workflow preserves them and does not mistake them for duplicate transcript records. Counts and Salmon abundance remain separate measures, and no biological comparison is claimed. See the [GSE297994 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297994) for dataset provenance.
