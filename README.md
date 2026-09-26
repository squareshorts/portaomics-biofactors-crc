# PORTA-OMICS BioFactors CRC portability study

Stage-0 reproducibility repository for the planned study:

**Portability-aware feature selection for cross-cohort colorectal cancer biomarkers with multi-omics validation**

This repository is intentionally frozen at the **data-integrity / preprocessing-audit stage**. It downloads public data, constructs sample manifests, verifies prespecified inclusion counts, maps probes to genes, checks cross-platform gene coverage, and creates single-sample rank representations. It does **not** fit a classifier, compute disease AUCs, rank disease-associated genes, tune PAFS, or inspect external predictive performance.

## Study separation

The previously submitted microbiome manuscript asks whether biomarker signals are portable across cohorts. This study asks a distinct question: whether portability can be built into **host-tissue biomarker discovery** itself. It uses new CRC transcriptomic cohorts and independent methylation/proteomic validation rather than the submitted PD/CRC microbiome datasets.

## Frozen transcriptomic cohorts

Development: GSE41258, GSE39582, GSE9348, GSE23878, GSE44861, GSE103512.

Locked external validation: GSE156451, GSE106582, GSE44076.

Orthogonal validation, not used for feature selection: GSE101764, GSE131013, CPTAC/PDC PDC000116.

See `config/datasets.csv` and `PROTOCOL_v1.0.md`.


## Windows note

For Windows, use **standard CPython 3.13**, not the MSYS2/UCRT64 Python interpreter under `C:\msys64`. The latter can force source builds of NumPy and related scientific packages because standard Windows wheels do not match its ABI. See `WINDOWS_SETUP.md`.

## Stage-0 run

Create a Python environment and install:

```bash
pip install -r requirements-stage0.txt
```

Then run:

```bash
python scripts/run_stage0.py
```

The script will:

1. download GEO Series Matrix files for the eight array cohorts;
2. download the processed RPKM supplement for GSE156451;
3. download GEO platform annotation tables for the array platforms;
4. build a sample-level manifest using accession-specific rules;
5. verify expected tumor/normal counts and paired subsets;
6. convert probe-level expression to one value per unambiguous gene using within-sample probe medians;
7. compute the cross-platform measurable-gene intersection;
8. create the prespecified single-sample percentile-rank / normal-score representation;
9. write a Stage-0 audit report without performing disease-association or classifier analyses.

Large downloaded data are excluded by `.gitignore`.

## Expected outputs

```text
results/manifests/sample_manifest.tsv
results/qc/cohort_sample_counts.csv
results/qc/pairing_summary.csv
results/qc/expression_dimensions.csv
results/qc/platform_mapping_summary.csv
results/qc/value_summary.csv
results/qc/gene_intersection.txt
results/qc/file_checksums.tsv
reports/stage0_audit.md
data/processed/gene_expression/*.parquet
data/processed/rank_normal_scores/*.parquet
```

If any frozen sample count does not match, Stage 0 fails rather than silently proceeding.

## Data sources

GEO accessions are downloaded from NCBI GEO programmatic endpoints. GSE156451 uses its per-sample processed RPKM TXT supplement rather than sequencing reads. Proteomic and methylation datasets are registered now but intentionally not pulled into feature selection.

## No-performance-peeking rule

Stage 0 must be completed and archived before `classifier_analysis_authorized` is changed. No tumor-normal AUC, differential-expression ranking, PAFS score, or external predictive result belongs in this stage.
