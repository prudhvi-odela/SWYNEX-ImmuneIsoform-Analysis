# SWYNEX Immune Isoform Analysis

A reproducible bioinformatics project built around NCBI GEO **GSE297994**, a human monocyte RNA-sequencing study of ex vivo and lipopolysaccharide (LPS)-activated conditions. This repository currently documents **Task 1: dataset selection and feasibility**. It does not claim that differential expression, isoform, or any other analyses have been run.

## Selected dataset

**Native long-read RNA sequencing of human monocytes reveals activation-induced alternative splicing toward functional isoforms** ([GSE297994](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297994)) studies primary human monocytes before and after 3-hour in-vitro LPS stimulation. GEO lists 22 samples across two sequencing platforms:

| Assay | Ex vivo | 3 h LPS | Platform |
|---|---:|---:|---|
| Short-read bulk RNA-seq | 8 healthy volunteers | 8 healthy volunteers | Illumina NextSeq 500, single-end 75 bp |
| Native long-read RNA-seq | 3 healthy volunteers | 3 healthy volunteers | Oxford Nanopore MinION |

The combination is relevant to immune activation and transcript isoform abundance: short reads provide broader sample coverage while long reads support transcript-level characterization. The study’s paper reports activation-associated isoform switching and describes the deposited data; those published findings are background, not results produced by this repository.

## Data footprint (snapshot: 2026-10-22)

- **Processed GEO supplementary tables:** 4 files, **4,662,174 bytes total** (about 4.45 MB decimal; 4.44 MiB), compressed.
- **Raw FASTQ files listed by ENA:** 22 runs, **16,631,011,185 bytes total** (about 16.63 GB decimal; 15.49 GiB), compressed. No raw reads were downloaded for Task 1.

See [docs/data_access.md](docs/data_access.md) for the files, exact sizes, sources, and access notes. The GEO processed tables are small enough for a later workflow; raw reads require substantially more storage and transfer time.

## Task 1 deliverables

- [Dataset selection rationale and design](docs/dataset_selection.md)
- [Sample group inventory](docs/sample_manifest.md)
- [Data files, sizes, and access plan](docs/data_access.md)

## Scope and limitations

This repository records metadata-based dataset selection only. No reads or count tables have been downloaded, quality checked, normalized, statistically tested, or analyzed here. No biological conclusions are asserted. Before a downstream comparison, inspect the sample-level metadata and define the analysis population, donor structure, and assay-specific comparisons. Do not assume that the short- and long-read samples are all matched across donors.

## Sources

- NCBI GEO Series record: [GSE297994](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297994)
- NCBI SRA Run Selector: [PRJNA1267092](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1267092)
- ENA run manifest and FASTQ byte counts: [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1267092&result=read_run&fields=run_accession%2Csample_accession%2Cexperiment_accession%2Cinstrument_platform%2Clibrary_strategy%2Cfastq_bytes&format=tsv)
- Bodelón et al. *Nature Communications* (2026), [article and DOI](https://doi.org/10.1038/s41467-026-73661-5)

Metadata and file sizes can change. Recheck the linked NCBI/ENA records before downloading or beginning analysis.
