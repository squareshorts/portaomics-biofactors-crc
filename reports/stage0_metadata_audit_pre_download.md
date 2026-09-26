# Stage-0 metadata audit before data retrieval

Date: 2026-09-22

This audit was completed from public accession metadata before any disease-classification analysis was run.

## Cohort inventory

| Accession | Role | Frozen primary samples | Platform / processing issue to preserve |
|---|---|---:|---|
| GSE41258 | Development | 186 primary CRC + 54 normal colon | Affymetrix U133A; series also contains adenoma, metastasis, non-colon normals and cell lines, all excluded |
| GSE39582 | Development | 566 CRC + 19 non-tumoral mucosa | Affymetrix U133 Plus 2; submitted processed values include RMA and ComBat; no additional cross-cohort batch correction will be applied |
| GSE9348 | Development | 70 CRC + 12 healthy controls | Affymetrix U133 Plus 2; submitted values use MAS5/global scaling, materially different from several other cohorts |
| GSE23878 | Development | 35 CRC + 24 normal | Affymetrix U133 Plus 2 |
| GSE44861 | Development | 56 CRC + 55 adjacent noncancerous | Affymetrix HT U133A |
| GSE103512 | Development | 57 CRC + 12 matched normal | Affymetrix HT U133+ PM; FFPE; colorectal subset must be isolated from a four-cancer series |
| GSE156451 | Locked external | 72 paired tumor/native tissue patients | RNA-seq, NovaSeq; processed supplement is per-sample RPKM; strongest platform-shift test |
| GSE106582 | Locked external | 68 complete tumor-mucosa pairs | Illumina HumanHT-12 v4; entire series contains 77 tumors and 117 mucosae |
| GSE44076 | Locked external | 98 complete tumor-adjacent-normal pairs | Affymetrix U219; additionally contains 50 healthy donor mucosae reserved for a stress test |

## Orthogonal validation inventory

- GSE101764: Illumina 450K methylation; 112 tumors and 149 mucosa samples, including 105 tumor-mucosa pairs. It is linked to GSE106582; transcript-methylation overlap will be established by patient/sample identifiers only after the transcriptomic signature is frozen.
- GSE131013: Illumina 450K methylation; paired tumor/adjacent mucosa from 96 individuals plus 48 healthy donors. It is paired conceptually with the COLONOMICS expression resource GSE44076.
- PDC000116: CPTAC prospective colon proteome, TMT10 global proteomics; registered for post-signature protein validation only.

## Important methodological consequence

The public cohorts do **not** share a common upstream normalization: examples include RMA, RMA+ComBat, MAS5/global scaling, quantile-normalized Illumina intensities, and RNA-seq RPKM. Therefore pooled expression values are not a defensible primary feature space.

The frozen primary representation remains the per-sample gene percentile-rank / normal-score transform after probe-to-gene mapping. It uses no test-cohort means, variances, class labels or pooled batch estimate. This is specifically intended to make the representation insensitive to cohort-specific monotone scale differences.

A later raw-data sensitivity analysis may be performed for array cohorts with usable raw files, but it is a sensitivity analysis and may not be used to redefine the primary method after external performance is known.

## Stage-0 gates

1. Frozen tumor/normal counts must match exactly.
2. GSE156451 must recover 72 complete tumor/native pairs.
3. GSE106582 must recover at least the prespecified 68 complete pairs and primary inclusion must be exactly 68+68.
4. GSE44076 must recover exactly 98 primary pairs and separately identify the 50 healthy controls.
5. GSE103512 must contain exactly 57 CRC and 12 colorectal-normal samples after exclusion of breast, prostate and lung samples.
6. Probe-to-gene mapping must be accession/platform-specific and ambiguous multi-gene probes are discarded.
7. The gene intersection is determined without disease labels.
8. No disease association, AUC, classifier, PAFS score or external predictive evaluation is allowed before this audit passes.
