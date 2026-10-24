🧬 SWYNEX Bioinformatics Internship

Biological Dataset Selection & Immune Isoform Analysis

Organization: SWYNEX Technologies
Domain: Bioinformatics
Current Task: Task 1 — Biological Dataset Selection
Project Status: Dataset selection and documentation

---

📌 1. Project Overview

This repository documents my work during the SWYNEX Technologies Bioinformatics Internship.

The first task focuses on selecting a publicly available biological dataset, documenting its source and size, understanding its experimental design, and defining a biologically relevant research question.

The selected dataset is GSE297994, a human monocyte RNA-sequencing study investigating activation-associated changes in transcript isoforms.

This repository currently documents dataset selection and feasibility. Future internship tasks will be added as separate sections or folders according to their assigned requirements.

🔬 2. Selected Biological Dataset

Attribute| Details
Dataset accession| GSE297994
Repository| NCBI Gene Expression Omnibus (GEO)
Organism| Homo sapiens
Biological material| Human monocytes
Research area| RNA sequencing and immune activation
Short-read platform| Illumina NextSeq 500
Long-read platform| Oxford Nanopore MinION
Associated BioProject| PRJNA1267092

Official Data Sources

- "NCBI GEO — GSE297994" (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297994)
- "NCBI BioProject — PRJNA1267092" (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1267092)

📊 3. Experimental Design

The documented study contains 22 samples across two RNA-sequencing approaches.

Assay| Ex vivo| 3-hour LPS activation| Platform
Short-read bulk RNA-seq| 8 samples| 8 samples| Illumina NextSeq 500, single-end 75 bp
Native long-read RNA-seq| 3 samples| 3 samples| Oxford Nanopore MinION
Total| 11| 11| Two sequencing approaches

The experimental conditions allow investigation of changes associated with lipopolysaccharide (LPS)-induced immune activation.

Important: The sample counts describe the study design. Donor matching across sequencing approaches must be verified from the original metadata before conducting paired comparisons.

💾 4. Dataset Size and Data Availability

Processed supplementary files

The four documented processed GEO supplementary files have a combined size of 4,662,174 bytes (approximately 4.66 MB decimal).

Data category| File type| Number of files
Raw count tables| Compressed text (".txt.gz")| 2
Salmon abundance tables| Compressed text (".txt.gz")| 2
Total| | 4

Raw sequencing archive

The associated sequencing archive lists 22 runs with a combined reported size of approximately 16.63 GB decimal.

The raw sequencing reads have not been downloaded for Task 1. The reported archive size is based on available run metadata, not on a local download or completed analysis.

For individual file sizes and access details, see ""docs/data_access.md"" (docs/data_access.md).

❓ 5. Research Question

How does LPS-induced immune activation affect RNA transcript isoform abundance in human monocytes, and how can short-read and native long-read RNA sequencing provide complementary information about these changes?

Biological Significance

Immune activation can influence gene expression and RNA processing. Investigating transcript abundance and isoform usage can help characterize these changes.

Short-read RNA sequencing provides gene- and transcript-level measurements, while long-read sequencing can help characterize longer transcript structures. Together, these approaches offer complementary perspectives on transcriptome changes.

This research question motivates potential downstream analysis. It is not a claim that new biological findings have been generated in this repository.

🎯 6. Why This Dataset Was Selected

The dataset was selected because:

- It is publicly available through NCBI GEO.
- It addresses a relevant biological question involving immune activation.
- It includes short-read and native long-read RNA sequencing.
- Its experimental design and sample counts can be documented.
- Its available processed data and archive metadata provide a basis for planning future computational work.

🛠️ 7. Potential Downstream Applications

Depending on the requirements of subsequent internship tasks, future work may include:

1. Reviewing sample metadata and experimental design.
2. Assessing sequencing data quality.
3. Preparing and processing RNA-seq data.
4. Examining gene expression and transcript isoform abundance.
5. Creating visualizations and interpreting results.
6. Documenting methods, findings, and limitations.

These are proposed future activities. They should not be interpreted as analyses already completed.

⚠️ 8. Limitations and Considerations

- The native long-read dataset contains three samples per condition.
- Donor matching between sequencing approaches must be verified before paired comparisons.
- Short-read and long-read quantification methods may produce different measurements.
- Raw sequencing data have not been downloaded or independently processed for Task 1.
- No new quality control, normalization, differential expression analysis, statistical testing, or independent biological validation has been performed as part of this task.
- Published study findings are background information and are not results generated by this repository.

📁 9. Repository Structure

SWYNEX-ImmuneIsoform-Analysis/
│
├── README.md
├── PROJECT_LOG.md
├── .gitignore
│
├── docs/
│   ├── dataset_selection.md
│   ├── data_access.md
│   └── sample_manifest.md
│
├── metadata/
├── data/
│   └── processed/
└── scripts/

The existing documentation provides additional details about dataset selection, data access, and sample metadata. New task-specific files will be added as the internship progresses and requirements become available.

📚 10. Task Progress

Deliverable| Status
Select a public biological dataset| Completed
Document the dataset source and accession| Documented
Document sample counts and data availability| Documented
Define a research question| Completed
Publish LinkedIn explanation post| Pending confirmation
Submit Task 1 for review| Pending confirmation
Receive reviewer approval| Awaiting confirmation
Subsequent internship tasks| Not started

This table reflects the current project documentation. Official task completion remains subject to SWYNEX reviewer approval.

👨‍💻 11. About Me

Prudhvi Teja Odela
Bioinformatics Undergraduate | Vignan University

My interests include genomics, transcriptomics, RNA sequencing, and computational biology. Through this internship, I aim to strengthen my practical skills and maintain a structured, reproducible record of my bioinformatics work.

- "GitHub Profile" (https://github.com/prudhvi-odela)
- "LinkedIn Profile" (https://www.linkedin.com/in/prudhvi-teja-odela/)

---

Disclaimer: This repository documents dataset selection and planning for Task 1. Potential future analyses are not presented as completed work or original research findings.
