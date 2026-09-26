# External-validation protocol amendment A1

Date: 2026-09-23

## Scope

This amendment was made after completion and audit of the six-cohort development freeze and before any external expression matrix or external outcome was loaded or evaluated.

No fitted model, selected feature set, hyperparameter, calibration parameter, or development result is changed.

## Reason for amendment

Adaptive PAFS selected `lambda_portability = 0`, with `m = 20` and `C = 0.03162277660168379`. The frozen M_ONLY model uses the same `m` and `C`, and the two frozen fits have identical fitted and calibration parameters.

Therefore, under the frozen specification, PAFS is mathematically identical to M_ONLY. The originally specified primary external contrast, `AUROC(PAFS) - AUROC(M_ONLY)`, is consequently deterministic and cannot serve as an inferential external-validation endpoint.

## Amended primary external estimand

The primary external endpoint is absolute transport performance of the frozen PAFS model:

- cohort-specific AUROC in each prespecified primary external analysis;
- 2,000-replicate participant-cluster bootstrap uncertainty within each cohort;
- unweighted macro-mean AUROC across the three primary external cohort analyses.

The PAFS-minus-M_ONLY contrast is retained only as an implementation identity check and is expected to equal exactly zero in every cohort if scoring is implemented correctly.

## Secondary predictive benchmarks

COHORT_ADJ and full-universe ElasticNet remain prespecified secondary predictive benchmarks. Their cohort-specific AUROC and unweighted macro-average will be reported descriptively. No new post-freeze superiority hypothesis is introduced.

The previously frozen secondary endpoints remain unchanged:

- AUPRC;
- balanced Brier score;
- calibration intercept;
- calibration slope;
- sensitivity at 0.5;
- specificity at 0.5.

The GSE44076 healthy-control analysis remains a separate secondary stress test and must be performed only after the prespecified primary tumor/adjacent-normal analysis.

Orthogonal methylation/protein validation remains restricted to the final adaptive-PAFS gene set and cannot feed back into transcriptomic model selection.

## Authorization

External expression loading: **AUTHORIZED**

External performance evaluation: **AUTHORIZED**

The authorization applies only to the one-shot evaluation of the already frozen models and endpoints described above. Retuning, model refitting using external outcomes, feature-selection changes, or feedback from external performance into the transcriptomic models remains prohibited.