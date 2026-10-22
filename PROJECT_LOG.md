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

## Dataset source
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297994
