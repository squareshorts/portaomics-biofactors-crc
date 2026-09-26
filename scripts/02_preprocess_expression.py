from __future__ import annotations

import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    PROCESSED,
    RAW_GEO,
    RAW_PLATFORM,
    RESULTS,
    collapse_probes_to_genes,
    load_registry,
    parse_rnaseq_txt_bytes,
    rank_normal_scores,
    read_platform_annotation,
    read_series_expression,
)


def load_manifest() -> pd.DataFrame:
    p = RESULTS / "manifests" / "sample_manifest.tsv"
    if not p.exists():
        raise SystemExit("Missing sample manifest. Run 01_build_manifest.py first.")
    return pd.read_csv(p, sep="\t", dtype=str).fillna("")


def read_gse156451_tar(path: Path, manifest: pd.DataFrame) -> pd.DataFrame:
    sample_to_title = dict(zip(manifest["sample_id"], manifest.get("Sample_title", "")))
    title_to_sample = {v: k for k, v in sample_to_title.items() if v}
    series = {}
    with tarfile.open(path, "r") as tar:
        for member in tar.getmembers():
            if not member.isfile() or not member.name.endswith(".txt.gz"):
                continue
            base = Path(member.name).name
            # Example: GSM4731764_N26-RNA.txt.gz
            stem = base[:-7] if base.endswith(".txt.gz") else Path(base).stem
            parts = stem.split("_", 1)
            gsm = parts[0] if parts and parts[0].startswith("GSM") else ""
            title = parts[1] if len(parts) > 1 else ""
            sample = gsm if gsm in sample_to_title else title_to_sample.get(title, "")
            if not sample:
                continue
            f = tar.extractfile(member)
            if f is None:
                continue
            series[sample] = parse_rnaseq_txt_bytes(f.read())
    if not series:
        raise ValueError("No GSE156451 sample expression files parsed from TAR.")
    df = pd.DataFrame(series)
    return df


def main() -> None:
    reg = load_registry()
    expr_reg = reg[reg["omic"].eq("transcriptome")]
    manifest = load_manifest()
    out_gene = PROCESSED / "gene_expression"
    out_rank = PROCESSED / "rank_normal_scores"
    out_gene.mkdir(parents=True, exist_ok=True)
    out_rank.mkdir(parents=True, exist_ok=True)

    dims = []
    maps = []
    vals = []
    gene_sets = {}

    for row in expr_reg.itertuples(index=False):
        acc = row.accession
        m = manifest[(manifest["accession"] == acc) & (manifest["primary_include"].str.lower().isin(["true", "1"]))]
        samples = m["sample_id"].tolist()
        if not samples:
            raise SystemExit(f"No primary samples for {acc}")

        if row.download_mode == "series_matrix":
            mat = read_series_expression(RAW_GEO / f"{acc}_series_matrix.txt.gz", use_samples=samples)
            annot = read_platform_annotation(RAW_PLATFORM / f"{row.platform}.annot.gz")
            gene, summ = collapse_probes_to_genes(mat, annot)
            maps.append({"accession": acc, "platform": row.platform, **summ})
        elif row.download_mode == "rnaseq_supplement":
            gene = read_gse156451_tar(RAW_GEO / f"{acc}_RAW.tar", m)
            gene = gene[[s for s in samples if s in gene.columns]]
            # Keep only simple gene symbols; duplicate rows were already removed per sample parser.
            gene.index = gene.index.astype(str).str.upper().str.strip()
            gene = gene[gene.index.str.match(r"^[A-Z0-9_.\-]+$")]
            gene = gene.groupby(level=0).median(numeric_only=True)
            maps.append({"accession": acc, "platform": row.platform, "n_input_probes": "", "n_annotated_unambiguous_probes": "", "n_gene_symbols": gene.shape[0]})
        else:
            continue

        if gene.shape[1] != len(samples):
            missing = sorted(set(samples) - set(gene.columns))
            raise SystemExit(f"{acc}: expected {len(samples)} primary samples but expression has {gene.shape[1]}; missing={missing[:10]}")

        gene = gene.replace([np.inf, -np.inf], np.nan)
        gene.to_pickle(out_gene / f"{acc}.pkl.gz", compression="gzip")
        gene_sets[acc] = set(gene.index[gene.notna().any(axis=1)])

        flat = gene.to_numpy(dtype=float)
        finite = flat[np.isfinite(flat)]
        vals.append({
            "accession": acc,
            "min": float(np.min(finite)) if finite.size else np.nan,
            "q01": float(np.quantile(finite, 0.01)) if finite.size else np.nan,
            "median": float(np.median(finite)) if finite.size else np.nan,
            "q99": float(np.quantile(finite, 0.99)) if finite.size else np.nan,
            "max": float(np.max(finite)) if finite.size else np.nan,
            "missing_fraction": float(np.isnan(flat).mean()),
        })
        dims.append({"accession": acc, "genes_before_intersection": gene.shape[0], "primary_samples": gene.shape[1]})

    common = set.intersection(*gene_sets.values())
    common_sorted = sorted(common)
    (RESULTS / "qc" / "gene_intersection.txt").write_text("\n".join(common_sorted) + "\n", encoding="utf-8")

    for row in expr_reg.itertuples(index=False):
        acc = row.accession
        path = out_gene / f"{acc}.pkl.gz"
        if not path.exists():
            continue
        gene = pd.read_pickle(path, compression="gzip")
        sub = gene.reindex(common_sorted)
        ranked = rank_normal_scores(sub)
        ranked.to_pickle(out_rank / f"{acc}.pkl.gz", compression="gzip")
        for d in dims:
            if d["accession"] == acc:
                d["genes_in_all_9_cohorts"] = len(common_sorted)
                d["ranked_missing_fraction"] = float(ranked.isna().to_numpy().mean())

    pd.DataFrame(dims).to_csv(RESULTS / "qc" / "expression_dimensions.csv", index=False)
    pd.DataFrame(maps).to_csv(RESULTS / "qc" / "platform_mapping_summary.csv", index=False)
    pd.DataFrame(vals).to_csv(RESULTS / "qc" / "value_summary.csv", index=False)

    if len(common_sorted) < 5000:
        raise SystemExit(f"Cross-platform intersection unexpectedly small ({len(common_sorted)} genes). Audit mappings before proceeding.")
    print(f"Preprocessing audit completed. Cross-platform intersection: {len(common_sorted)} genes.")


if __name__ == "__main__":
    main()
