# Data access and file-size snapshot

**Checked:** 2026-10-22. Sizes below are the byte counts displayed by the live NCBI GEO record for processed files and the ENA run report for raw reads. They describe compressed files. They are metadata snapshots, not downloads or checksums.

## Processed supplementary files at GEO

| File | Contents as described by GEO | Bytes | GEO display |
|---|---|---:|---:|
| GSE297994_Native_long-read_RNASeq_raw_count.txt.gz | Native long-read raw counts (gene/isoform; STAR/FLAIR annotation in GEO description) | 732,958 | 715.8 KB |
| GSE297994_Native_long-read_RNASeq_salmon_abundance.txt.gz | Native long-read Salmon abundance | 903,345 | 882.2 KB |
| GSE297994_short-read_RNASeq_raw_count.txt.gz | Short-read raw counts (gene-level) | 694,687 | 678.4 KB |
| GSE297994_short-read_RNASeq_salmon_abundance.txt.gz | Short-read Salmon abundance | 2,331,184 | 2.2 MB |
| **Total** | **4 compressed tables** | **4,662,174** | **about 4.45 MB decimal / 4.44 MiB** |

File links are served from NCBI GEO’s [supplementary directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE297nnn/GSE297994/suppl/). The linked directory and exact filenames are also listed on the [GSE297994 record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297994).

## Raw sequence data at ENA/SRA

The NCBI SRA project linked from GEO is [PRJNA1267092](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1267092). Its ENA run report returned 22 FASTQ records. The sum of their reported fastq_bytes is **16,631,011,185 bytes** (16.63 GB decimal; 15.49 GiB):

| Sequencing platform | Runs/samples | Reported compressed FASTQ bytes |
|---|---:|---:|
| Oxford Nanopore | 6 | 9,398,292,595 (~9.40 GB) |
| Illumina | 16 | 7,232,718,590 (~7.23 GB) |
| **Total** | **22** | **16,631,011,185 (~16.63 GB)** |

The manifest can be rechecked through the [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1267092&result=read_run&fields=run_accession%2Csample_accession%2Cexperiment_accession%2Cinstrument_platform%2Clibrary_strategy%2Cfastq_bytes&format=tsv). ENA’s byte counts are the source for the table; the NCBI GEO page directs users to SRA for raw data. No raw FASTQ or SRA files were downloaded for Task 1.

## Access plan

1. For a first downstream prototype, retrieve only the processed supplementary tables needed for the chosen question and record their source URLs and checksums.
2. Inspect table headers and sample annotations before deciding which assay-specific tables are comparable.
3. Download raw reads only if the planned workflow requires reprocessing or raw-read QC and available storage/compute can accommodate the compressed inputs plus expanded intermediate files.
4. Save a dated run manifest and checksums at download time; archive version, command, and software details in the analysis record.

GEO/ENA records and file sizes may change. Recheck before a transfer. This document does not assert that files have been downloaded, inspected, or analyzed.
