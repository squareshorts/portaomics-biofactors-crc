# Orthogonal Data Availability Audit (A2)

## 1. Methylation Data
**Candidate Datasets:** TCGA-COAD and TCGA-READ GDC DNA Methylation (Illumina HumanMethylation450).
**Local Paths Found:** NONE. No TCGA methylation resources, matrices, or downloading scripts currently exist in the repository.
**Sample Types Expected:** Primary tumor versus solid-tissue normal.
**Paired Tumor/Normal Data:** Yes, a subset of patients in TCGA have paired primary tumor and adjacent solid-tissue normal samples.
**Annotation Availability:** To be downloaded (Illumina 450K manifest or GDC annotations required to map probes to TSS200/TSS1500/5'UTR/1stExon).
**Expected Download Size:** ~1-5 GB for pre-compiled beta-value matrices; potentially >50 GB if raw IDATs are required.

## 2. Proteomics Data
**Candidate Datasets:**
1. CPTAC prospective/confirmatory colon global proteome (Primary priority, referenced as PDC000116 in `datasets.csv`).
2. CPTAC retrospective TCGA_Colon_Cancer_Proteome (Fallback).
**Local Paths Found:** NONE. `config/datasets.csv` contains a metadata entry for `PDC000116` but no local data, processed matrices, or scripts exist.
**Sample Types Expected:** Primary tumor versus adjacent normal.
**Paired Tumor/Normal Data:** Yes, CPTAC typically provides paired tumor/normal tissues.
**Annotation Availability:** To be downloaded (Protein to gene symbol mapping required).
**Expected Download Size:** ~50-200 MB for tabular protein expression matrices.

**Conclusion:** All external orthogonal validation datasets (TCGA and CPTAC) are currently absent from the local repository and must be acquired before evaluation can proceed.
