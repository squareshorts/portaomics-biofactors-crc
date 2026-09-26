# External Validation Execution Incident Report A1

**Date:** 2026-09-24

**Original Script SHA-256:** `DB547F0C1D4A8CDBD093C2335F40B8205D75D7A1970A1AC3CBEFDFD14B81AE84`

### Failure Description

- **Exact Exception Type:** `KeyError`
- **Exact Failure Location:** `expr_sub = expr.loc[sample_ids]` in `scripts/14_one_shot_external_validation.py` (line ~233).
- **Missing Samples:** 50 missing samples (`GSM1077598` through `GSM1077647`).
- **Identity Class:** These samples correspond to the `normal_healthy` stress-test group from GSE44076.
- **Analysis Stage Failed:** The failure occurred *only* when entering the prespecified GSE44076 healthy-control secondary stress test.
- **Computational State:** The primary analyses for GSE156451, GSE106582, and the primary tumor/adjacent-normal analysis for GSE44076 had already been successfully processed computationally in memory.
- **Result Files Written:** None. No partial output files were written to disk, as file writing is safely deferred to the end of the script after all computational stages complete. 

### Analytical Validation

- **Confirmation:** No analytical choice, hyperparameter, threshold, phenotype definition, or model specification was changed in response to this incident.

### Mechanical Correction

- **Exact Mechanical Correction:** 
  1. Implemented execution atomicity by separating a structural validation phase (which verifies the availability of all required expression matrices, samples, and genes across *all* analysis paths) that must strictly pass before *any* metrics or predictions are computed.
  2. The 50 healthy samples were absent from `data/processed/rank_normal_scores/GSE44076.pkl.gz` because Stage-0 processing explicitly restricted output to samples where `primary_include` is true. The recovery dynamically loads the missing IDs directly from the GSE44076 raw matrix (`data/raw/geo/GSE44076_series_matrix.txt.gz`).
  3. Added a fallback to treat independent healthy samples lacking an explicit `participant_id` (empty or `NaN`) as independent participant units using their `sample_id` for valid bootstrap clustering logic.
- **Justification:** The correction concerns solely data-addressing and ensures preprocessing path consistency. The exact frozen Stage-0 processing logic (`collapse_probes_to_genes` and `rank_normal_scores`) is invoked on the raw matrix. Because `rank_normal_scores` operates independently within each sample, processing these 50 samples exactly reproduces the frozen normalization logic without affecting the normalization or results of any other sample.

**Repaired Script SHA-256:** `0D18035BB92610BDE60835B95902853D21AEB9D688C9588AAB1CBB49D265BB2E`
