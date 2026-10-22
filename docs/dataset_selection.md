# Task 1 — Dataset selection and rationale

## Selected study

- **Accession:** [GSE297994](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297994)
- **Title:** Native long-read RNA sequencing of human monocytes reveals activation-induced alternative splicing toward functional isoforms
- **Organism:** *Homo sapiens*
- **GEO experiment type:** Expression profiling by high throughput sequencing
- **GEO status:** Public on 2026-05-04; record last updated 2026-09-24 (checked 2026-10-22)
- **BioProject:** [PRJNA1267092](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1267092)

## Why this dataset fits the project

The study profiles human monocytes in ex vivo and 3-hour LPS-stimulated conditions, directly connecting innate immune activation with transcript abundance. It combines short-read bulk RNA-seq and native Oxford Nanopore long-read RNA-seq, making it a suitable candidate for a project concerned with immune activation and isoform-level signal. The GEO record lists processed gene/isoform counts and Salmon abundance tables, allowing a later workflow to begin from compact processed data if raw-read reprocessing is outside the available compute budget.

This is a selection rationale based on the public record and publication description. It is not a reanalysis or independent validation of the paper’s conclusions.

## Sample design recorded by GEO

| Data type | Ex vivo | 3-hour LPS in vitro | Instrument |
|---|---:|---:|---|
| Short-read bulk RNA-seq | 8 healthy volunteers | 8 healthy volunteers | Illumina NextSeq 500; single-end, 75 bp |
| Native long-read RNA-seq | 3 healthy volunteers | 3 healthy volunteers | Oxford Nanopore MinION |
| **Total** | **11 samples** | **11 samples** | **22 samples** |

The series-level design describes short-read data from 8 volunteers and long-read data from 3 volunteers across the two conditions. The sample titles number volunteers within groups. Confirm donor identity, pairing, and any missingness from the full sample metadata before specifying a statistical model. Keep assay types and their available sample sizes explicit in downstream comparisons.

## Feasibility and limits

- The four compressed processed tables listed on GEO total **4,662,174 bytes** (~4.45 MB decimal).
- The ENA run manifest lists 22 compressed FASTQ files totaling **16,631,011,185 bytes** (~16.63 GB decimal / 15.49 GiB).
- Raw-read processing will need materially more disk space than the compressed input footprint for temporary files, alignments, indexes, and outputs. No raw data were downloaded for this selection task.
- The GEO record provides short-read raw gene counts, native long-read raw gene/isoform counts, and Salmon isoform abundance tables. Their meanings and quantification pipelines differ; do not combine the tables without defining and justifying the measurement scale.
- This project has not performed quality control, normalization, differential testing, isoform switching analysis, or biological interpretation.

## Selection decision

**Retain GSE297994 as the Task 1 candidate** because it contains the relevant human immune-cell context, explicit LPS perturbation, both short- and long-read assays, and compact processed data for a later reproducible workflow. Revisit the decision if the intended analysis requires more than three long-read biological replicates per condition or fully matched multi-platform donor measurements.

## References

1. NCBI GEO, [GSE297994 record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297994).
2. NCBI SRA, [BioProject PRJNA1267092](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1267092).
3. Bodelón A, et al. *Native long-read RNA sequencing of human monocytes reveals activation-induced alternative splicing toward functional isoforms*. *Nature Communications*. 2026;17:6982. [doi:10.1038/s41467-026-73661-5](https://doi.org/10.1038/s41467-026-73661-5).
