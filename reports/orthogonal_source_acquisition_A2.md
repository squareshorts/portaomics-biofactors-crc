# Orthogonal Source Acquisition A2

## TCGA Methylation Source
- **API endpoint:** `https://api.gdc.cancer.gov/files`
- **Query/Filter specification:** `data_category` = "DNA Methylation", `data_format` = "txt", `platform` = "Illumina Human Methylation 450", `cases.project.project_id` in ["TCGA-COAD", "TCGA-READ"], `cases.samples.sample_type` in ["Primary Tumor", "Solid Tissue Normal"].
- **GDC Release/API Status:** Active, queried via /files endpoint.
- **COAD HM450 participant counts:** 38 paired participants (out of larger total).
- **READ HM450 participant counts:** 7 paired participants.
- **Combined paired participant count:** 45.
- **Selected beta-file count:** 90 (1 primary tumor and 1 solid tissue normal per participant).
- **Downloaded bytes:** 0 (download deferred pending annotation audit).
- **Hashes:** N/A (downloads pending).
- **Annotation source and hash:** `HM450.hg38.manifest.gencode.v36.tsv.gz`. Expected MD5: `e163fc110043abb5a7ef623816383bb9`. Currently unable to download from public GDC endpoints.

## CPTAC Proteomics Source
- **PDC Study:** PDC000116 (Prospective Colon PNNL Proteome Qeplus).
- **Exact study version UUID:** `bbc1441e-57b8-11e8-b07a-00a098d917f8`.
- **Current study version metadata:** Experiment Type: TMT10, Analytical Fraction: Proteome, Cases Count: 102.
- **Quantitation type locked:** Protein Assembly (Log2 Ratio - Only Unshared Peptides).
- **Tumor count:** 97.
- **Normal count:** 100.
- **Paired participant count:** 96.
- **Downloaded processed quantitation file(s):** Pending direct download (candidate file: `CPTAC2_Colon_Prospective_Collection_PNNL_Proteome.tmt10.tsv`).
- **Hashes:** Expected MD5: `6cc8208af6da2554ade54de2afe090d8`.
- **Protein-to-gene annotation fields available:** To be determined upon file extraction.
