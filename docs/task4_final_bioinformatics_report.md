# Task 4 — Final Bioinformatics Report

**Project:** Immune Isoform Analysis  
**Internship:** SWYNEX Technologies  
**Dataset:** GSE297994  
**Report scope:** Tasks 1–3 completed in this repository

## Executive summary

This project uses the public GSE297994 dataset to explore how short-read and native long-read RNA-sequencing data can describe human monocytes before and after 3-hour LPS stimulation. The repository documents dataset selection, conservative cleaning of four processed expression tables, and descriptive visualizations of their count patterns.

The Task 3 figures show two notable patterns in the processed counts: the 100 most abundant short-read gene rows account for a larger share of counts in every LPS sample than in its ex vivo counterpart, and samples separate by condition along the first PCA axis in both sequencing technologies. Median-of-ratios scaling aligns the central short-read distributions after CPM medians differ markedly. These results are exploratory summaries. They are not differential-expression or differential-isoform tests and do not establish which genes or transcripts respond to LPS.

The appropriate conclusion is that these processed data show condition-associated structure and a count-composition shift that should be considered before formal inference. The project's broader isoform-level question remains open in this repository.

## Research question

How do the available processed short-read gene counts and native long-read transcript counts vary between ex vivo and 3-hour LPS-treated human monocytes, and what patterns or data limitations should be considered before testing gene-expression or isoform changes?

This scoped question reflects the analyses actually performed here. The repository does not independently reproduce the original study's differential-expression, isoform-switching, splicing, or protein analyses.

## Dataset and study design

GSE297994 contains processed expression tables for primary human CD14-positive monocytes sampled ex vivo and after 3 hours of in-vitro LPS activation. The documented design includes 16 short-read samples (8 per condition) and 6 native long-read samples (3 per condition). Short-read sequencing used Illumina NextSeq 500; native long-read sequencing used Oxford Nanopore MinION. The sample inventory is in [`metadata/sample_metadata.tsv`](../metadata/sample_metadata.tsv).

The repository contains four processed matrices: short-read raw gene counts, short-read Salmon isoform abundance, long-read raw transcript-row counts, and long-read Salmon isoform abundance. Task 2 retained these measurement types separately. In the long-read raw-count matrix, repeated `gene_id` annotations occur across transcript rows; the separate `ids` field is the row-level identifier and those rows were preserved.

## Methods

### Data preparation

Task 2 trimmed surrounding whitespace, checked table widths and numeric values, preserved all nonblank rows and original compressed inputs, and wrote cleaned copies and a JSON summary. It did not impute values, remove features, normalize counts, or conduct biological testing. The method and table-level results are described in [`docs/task2_data_cleaning.md`](task2_data_cleaning.md).

### Descriptive visualization

Task 3 used the cleaned raw-count matrices and sample metadata. Short-read and long-read data were analyzed separately; Salmon abundance tables were not mixed with raw counts.

- **Count composition:** for each short-read sample, the share of all counts represented by that sample's own 100 highest-count gene rows was calculated. The top-100 set may differ between samples.
- **Normalization display:** short-read CPM values were compared with counts divided by median-of-ratios size factors, then displayed as `log2(value + 1)`. The implementation follows the median-of-ratios idea described for DESeq2, but this plot is not a DESeq2 differential-expression analysis.
- **PCA:** within each technology, features with CPM ≥ 1 in at least the smallest condition group were retained (8 samples for short-read; 3 for long-read). PCA was performed on the 500 features with the highest variance in normalized `log2(count + 1)` values.
- **Sample correlation:** Pearson correlations were calculated between samples using all features passing the CPM filter and the normalized log-count values.

The four figures and their full interpretive notes are in [`docs/task3_visualization_interpretation.md`](task3_visualization_interpretation.md). They can be regenerated with `python scripts/visualize_task3.py` from the repository root (Python 3.9+, standard library only).

## Results and interpretation

### 1. Data integrity and scope

The Task 2 workflow reported consistent row widths and no blank rows, empty cells, missing row identifiers, nonnumeric values, or negative quantitative values in the four processed tables. It retained repeated gene annotations in the long-read raw-count table because its transcript-level row IDs are unique. These checks support using the supplied processed matrices for the descriptive work; they do not validate the original sequencing, quantification, or biological sample identity.

### 2. Short-read count composition

The 100 highest-count rows account for **24.9–32.0%** of total counts in ex vivo samples and **42.6–50.8%** in LPS samples. The share is higher for LPS in all eight short-read sample pairs, with reported paired differences of 15.1–18.9 percentage points.

This indicates that observed counts are more concentrated among the top-count rows in the LPS samples. Because each sample's top-100 set is selected independently, the comparison is about concentration, not a fixed list of genes. The chart does not identify which rows drive the change, show differential expression, prove that LPS caused it, or establish statistical significance. It does show why total-count/CPM scaling alone may give a different picture from composition-aware scaling.

### 3. Short-read normalization display

Among expressed features passing the stated filter, median CPM values were **22.2–24.4** for ex vivo and **7.6–10.0** for LPS samples. After division by median-of-ratios size factors, the per-sample median `log2(normalized count + 1)` values fell between **6.70 and 7.00**.

The plot demonstrates alignment of the central distributions under this normalization for these data. It does not establish that the normalization assumption is correct: median alignment is not external validation, and it cannot rule out broad biological shifts. The top-100 composition pattern makes this assumption especially important to revisit before inferential analysis.

### 4. PCA and sample correlations

In the selected 500-feature PCA, PC1 explained **95.45%** of the variance in the short-read data and **94.86%** in the long-read data. Samples separated by condition along PC1 in both technology-specific panels. These are exploratory visual patterns; they do not identify the features responsible, quantify statistical evidence for a treatment effect, or establish generalizability beyond the sampled donors.

The correlation heatmaps also show condition-associated grouping. Within the long-read ex vivo samples, HC3 correlated **0.725** with HC1 and **0.785** with HC2, while HC1 and HC2 correlated **0.901**. HC3 had the smallest long-read raw-count total (**329,056**). This is a sample-level observation only: the data do not establish that sequencing depth caused the lower correlations or that HC3 should be excluded.

### 5. Relationship to the published study

The publication associated with GSE297994 reports a much broader analysis, including differential gene and isoform expression and isoform switching. Those are the authors' published findings, not results generated by this repository. This project currently reproduces data preparation and descriptive summaries of supplied processed matrices; it does not claim to reproduce or independently validate the paper's inferential analyses.

## Limitations

- The work uses supplied processed matrices; raw sequencing reads were not reprocessed in this project.
- Task 2 performs integrity checks and formatting cleanup only. It does not verify upstream alignment, quantification, annotation quality, or sample identity.
- Task 3 is descriptive. It calculates no differential-expression or differential-usage model, p-value, confidence interval, or multiple-testing-adjusted result.
- The median-of-ratios method depends on assumptions about the distribution of feature changes. The observed count-composition skew means those assumptions require scrutiny; median alignment alone is not proof of correct normalization.
- The CPM threshold and variable-feature selection are exploratory choices. PCA and correlation patterns depend on those choices, normalization, and log transformation.
- There are only three long-read samples per condition, limiting how strongly patterns can support population-level conclusions. Short-read has eight samples per condition.
- Sequencing technologies were analyzed independently. Their quantification approaches and feature definitions differ, so the displayed PCA/correlation patterns are not direct cross-platform comparisons.
- Sample titles suggest donor numbering, but cross-technology donor matching must be verified from source metadata before paired cross-platform analysis.
- No external validation or functional experiment was performed by this repository.

## Conclusion

Across the available processed matrices, Task 2 confirms that the tables are structurally usable for downstream analysis, and Task 3 reveals clear descriptive condition-associated structure plus a marked short-read count-composition difference. The results justify careful normalization and sample-aware modeling in a future inferential analysis. They do not, by themselves, establish differential genes, isoform switching, causality, or biological mechanism.

## Reproducibility and project files

From the repository root:

```bash
python scripts/clean_and_summarize.py
python scripts/visualize_task3.py
```

The scripts use Python 3.9+ and the standard library only. Task 2 writes cleaned data and `results/task2_data_summary.json`; Task 3 writes four SVG figures and `results/task3_sample_metrics.tsv`. See [`README.md`](../README.md), [`PROJECT_LOG.md`](../PROJECT_LOG.md), and the Task 2 and Task 3 reports for details.

## References

1. Bodelón A, van Haaren MJH, Sobrevals Alcaraz P, et al. Native long-read RNA sequencing of human monocytes reveals activation-induced alternative splicing toward functional isoforms. *Nature Communications*. 2026;17:6982. [doi:10.1038/s41467-026-73661-5](https://doi.org/10.1038/s41467-026-73661-5).
2. Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. *Nucleic Acids Research*. 2002;30(1):207–210. [GEO record GSE297994](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297994).
3. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014;15:550. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).
