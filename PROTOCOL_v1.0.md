# Frozen protocol v1.0

Locked: 2026-09-22

## Scientific question

Does transcriptomic feature selection that explicitly rewards cross-cohort consistency and penalizes cohort dependence improve generalization to unseen colorectal-cancer cohorts, while yielding a final signature with independent methylation and protein support?

## Primary phenotype

Primary colorectal carcinoma tissue versus histologically non-neoplastic colorectal mucosa. Metastases, adenomas, recurrent lesions, unrelated cancer types, and cell lines are excluded.

The study concerns **tissue-state molecular biomarkers**. It does not claim population screening performance.

## Cohort roles

Six cohorts are development data: GSE41258, GSE39582, GSE9348, GSE23878, GSE44861, GSE103512.

Three cohorts are locked external validation data: GSE156451, GSE106582, GSE44076.

For GSE106582, the primary external analysis is restricted to the 68 complete tumor-mucosa pairs. For GSE44076, the primary external analysis is restricted to the 98 complete tumor-adjacent-normal pairs; the 50 healthy-donor mucosae are reserved for the prespecified control-source stress test.

GSE101764, GSE131013, and PDC000116 are orthogonal validation resources only. Their data may not influence transcriptomic feature selection or hyperparameter tuning.

## Preprocessing lock

No cross-cohort ComBat or other pooled batch correction is permitted in the primary analysis.

Array data are represented at gene level by mapping probes through the corresponding GEO platform annotation. A probe is retained only when it maps unambiguously to one gene symbol. When multiple retained probes map to one gene, the within-sample median is used.

RNA-seq GSE156451 uses the submitted per-sample RPKM files. Gene identifiers are harmonized to gene symbols before the cross-platform intersection is determined.

The primary feature universe is the intersection of genes technically measurable in all nine transcriptomic cohorts. This intersection is determined without disease labels.

Each sample is independently transformed to percentile ranks across the retained gene universe and then to normal scores:

r_ij = (rank(x_ij) - 0.5) / p

z_ij = Phi^{-1}(r_ij)

The transformation does not use test-cohort means, variances, disease labels, or batch estimates.

## Weighting lock for later modelling

Within any training partition, cohorts have equal total weight and cases/controls have equal weight within cohort:

w_i = 1 / (2 K n_{k,y_i})

where K is the number of training cohorts.

## PAFS lock for later modelling

For gene j and training cohort k:

e_jk = 2(AUC_jk - 0.5)

Disease magnitude:

M_j = |median_k(e_jk)|

Directional concordance:

S_j = K^{-1} sum_k I[sign(e_jk) = sign(median_k(e_jk))]

Heterogeneity:

H_j = MAD_k(e_jk)

Cohort dependence is the incremental partial R^2 of cohort after disease status in a training-only model.

All four quantities are converted to percentile ranks oriented so that one is desirable. The locked PAFS score is their unweighted mean. No learned component weights are permitted in v1.0.

Candidate signature sizes: 5, 10, 20, 50, 100 genes.

## Validation lock for later modelling

Outer validation is leave-one-development-cohort-out. Inner leave-one-cohort-out tuning selects signature size and ridge penalty. Feature selection is recomputed from scratch inside every inner fold.

The primary comparator ranks genes by disease magnitude alone and uses the same downstream ridge-logistic classifier. Additional comparators are pooled cohort-adjusted association ranking and full-universe ElasticNet.

External cohorts remain untouched until all method choices are frozen from development data.

## Stage-0 prohibition

The current repository stage is limited to data retrieval, metadata parsing, inclusion verification, feature annotation, cross-platform gene coverage, single-sample rank preprocessing, and non-outcome-based QC.

Do not compute differential expression, disease AUC, classifier predictions, PAFS, model calibration, or external predictive performance during Stage 0.
