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
- Preserve duplicate feature IDs and report their number; do not drop or aggregate them.
- Stop if a nonblank row has a different number of fields from the header, or if a header is blank or duplicated.
- For each non-ID column, report nonempty-cell count, numeric-value count, mean, minimum, maximum, zero and negative counts, and nonnumeric-value count. Missing and non-finite values are excluded from numeric summaries.

This is formatting and integrity cleanup only. It does not normalize counts, infer sample metadata, remove outliers, perform differential testing, or establish biological findings.

## Results

Run the command above to generate observed counts and summaries from the checked-in tables. Values are not copied from the publication. See the [GSE297994 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297994) for dataset provenance.
