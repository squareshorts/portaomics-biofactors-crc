# PORTA-OMICS CRC cross-cohort transcriptomic study

Reproducibility repository for a staged colorectal-cancer (CRC) transcriptomic study with locked external validation and orthogonal methylation/proteomic validation.

The repository preserves analysis locks, executable scripts, model artifacts, derived results, execution manifests, audit reports, and cryptographic hashes. Raw public molecular datasets are not redistributed.

## Study design

Six GEO cohorts were used for development:

- GSE41258
- GSE39582
- GSE9348
- GSE23878
- GSE44861
- GSE103512

Three cohorts were withheld from model development and used for locked external validation:

- GSE156451: 72 tumor / 72 normal in the primary paired analysis
- GSE106582: 68 tumor / 68 normal in the primary paired analysis
- GSE44076: 98 tumor / 98 adjacent-normal in the primary paired analysis

GSE44076 also supplied 50 independent healthy controls for a prespecified secondary stress test.

Orthogonal validation used paired TCGA-COAD/READ methylation data and CPTAC/PDC study PDC000116 proteomics.

## Analysis chronology

The repository is intentionally staged. Earlier protocol files are retained as provenance and should be read together with the later locks and amendments.

1. **Stage 0 — data integrity and harmonization.** Public transcriptomic cohorts were retrieved, sample inclusion was audited, probe-level arrays were collapsed to unambiguous gene symbols, and the common measurable feature universe was frozen at 10,998 genes. Each sample was transformed independently to percentile-rank normal scores. No disease-association or predictive performance was used at this stage.
2. **Stage 1A — initial development-only nested LOCO analysis.** The original equal-component portability-aware feature-selection formulation was evaluated using only the six development cohorts.
3. **Stage 1B — adaptive portability analysis.** A development-only amendment introduced a portability weight `lambda` selected strictly inside nested leave-one-cohort-out validation. The frozen adaptive score was
   `rank_M + lambda * mean(rank_S, rank_low_H, rank_low_C)`.
4. **Stage 1C — ElasticNet comparator.** A full-universe ElasticNet benchmark was tuned using development data only.
5. **Final development freeze.** The selected PAFS model used 20 genes, ridge `C = 0.0316227766`, and `lambda = 0`. Under this frozen specification, PAFS became mathematically identical to the disease-magnitude-only model. The selected development LOCO macro-AUROC was 0.9805.
6. **External protocol amendment A1.** Before any external expression matrix or external outcome was loaded, the primary external estimand was changed from the now-deterministic PAFS-minus-M_ONLY contrast to absolute AUROC of the frozen 20-gene PAFS model in each external cohort, with a 2,000-replicate participant-cluster bootstrap and unweighted macro-average.
7. **Locked external validation.** Primary AUROCs were 0.9782 (GSE156451), 0.9922 (GSE106582), and 0.9999 (GSE44076), giving an unweighted macro-AUROC of 0.9901. The separate GSE44076 healthy-control stress test yielded AUROC 0.9998.
8. **Orthogonal multi-omics validation A2.** All 20 signature genes were testable for promoter methylation; 11/20 followed the expected inverse direction. Nine genes were testable in CPTAC proteomics; all 9/9 followed the transcriptomically expected direction and remained significant after Benjamini-Hochberg correction.

## External-validation execution incident

The first external-validation execution reached the prespecified GSE44076 healthy-control stress-test path and stopped because those 50 samples had intentionally not been written to the Stage-0 primary-only processed matrix. No result files were written by the failed run.

The correction was mechanical: a preflight/atomicity check was added, the 50 healthy samples were processed using the already frozen sample-wise transformation, and missing participant identifiers for independent healthy samples were replaced by their sample IDs for bootstrap clustering. No model, feature, hyperparameter, threshold, phenotype definition, or external endpoint was changed. The full incident record is retained in `reports/external_validation_execution_incident_A1.md`.

Accordingly, this repository does not describe the external evaluation as an uninterrupted "one-shot" execution; it documents the failed execution and the constrained mechanical repair explicitly.

## Key frozen artifacts

- `config/final_development_lock.json`
- `config/final_external_lock_pre_A1.json`
- `config/final_external_lock.json`
- `reports/final_external_protocol_amendment_A1.md`
- `reports/final_development_freeze.md`
- `reports/external_validation_execution_incident_A1.md`
- `results/final_development/`
- `results/external_validation/`
- `config/orthogonal_source_lock_A2.json`
- `config/orthogonal_analysis_lock_A2.json`
- `results/orthogonal_validation/`

## Main analysis scripts

- `scripts/00_download_geo.py`
- `scripts/01_build_manifest.py`
- `scripts/02_preprocess_expression.py`
- `scripts/03_stage0_report.py`
- `scripts/10_stage1_nested_loco.py`
- `scripts/11_stage1b_adaptive_loco.py`
- `scripts/12_stage1c_elasticnet.py`
- `scripts/13_final_development_freeze.py`
- `scripts/14_one_shot_external_validation.py`
- `scripts/15_orthogonal_multiomics_validation.py`

The filename `14_one_shot_external_validation.py` is retained for provenance; the documented execution incident above should be consulted when interpreting the execution history.

## Reproducing the data layer

Raw GEO, TCGA/GDC, and PDC inputs are intentionally excluded from Git. Stage-0 acquisition and harmonization instructions are retained in `PROTOCOL_v1.0.md`, `WINDOWS_SETUP.md`, `STAGE0_COLAB.ipynb`, and `run_stage0.ps1`.

The root `PROTOCOL_v1.0.md` is the original Stage-0 protocol and is not a complete description of the later adaptive development, external-amendment, or orthogonal-validation stages. Later locks and reports listed above supersede it where explicitly documented.

## Software environment

Exact software versions used in the completed analyses are recorded in the execution/provenance manifests. In particular, the development and external-validation manifests record Python 3.14.6, NumPy 2.4.6, pandas 3.0.3, and SciPy 1.18.0 for the completed downstream analyses.

The historical `requirements-stage0.txt` and `WINDOWS_SETUP.md` describe the earlier Stage-0 setup path. A release-level environment specification should therefore be interpreted together with the execution manifests rather than assuming that the Stage-0 dependency file describes the complete downstream environment.

## Data sources

Transcriptomic data are publicly available from NCBI GEO under GSE41258, GSE39582, GSE9348, GSE23878, GSE44861, GSE103512, GSE156451, GSE106582, and GSE44076.

Methylation validation uses TCGA-COAD and TCGA-READ data obtained through the Genomic Data Commons. Proteomic validation uses Proteomic Data Commons study PDC000116.

## Scope

The frozen model is a research classifier for tumor-versus-non-neoplastic colorectal tissue. The case-control external cohorts do not represent clinical prevalence, and the frozen probability values should not be interpreted as absolute clinical CRC risk.

## Citation and archival release

The manuscript-associated reproducibility archive is frozen as GitHub release `v1.0.0` and archived on Zenodo.

- Zenodo DOI: [10.5281/zenodo.22979172](https://doi.org/10.5281/zenodo.22979172)
- GitHub release: [v1.0.0](https://github.com/squareshorts/portaomics-biofactors-crc/releases/tag/v1.0.0)
- Frozen release commit: `51810d06ffe6aa9e94548b6b91981d2d3d0ea791`

Suggested software citation:

> Gama, Jessica Silva, and Antonio Pereira. (2026). *PORTA-OMICS CRC: Cross-cohort colorectal cancer transcriptomic validation with multi-omics support* (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.22979172


## Funding

This work was supported by the Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq), Grant 309589/2023-1 to A.P., and by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES).
