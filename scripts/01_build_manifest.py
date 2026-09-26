from __future__ import annotations

import pandas as pd

from common import RAW_GEO, RESULTS, classify_sample, load_registry, make_primary_include, sample_metadata_frame


def num(x: str):
    return int(x) if str(x).strip() else None


def main() -> None:
    reg = load_registry()
    expr = reg[reg["omic"].eq("transcriptome")]
    all_frames = []
    summaries = []
    pair_rows = []
    failures = []

    for row in expr.itertuples(index=False):
        acc = row.accession
        path = RAW_GEO / f"{acc}_series_matrix.txt.gz"
        if not path.exists():
            raise SystemExit(f"Missing {path}. Run 00_download_geo.py first.")
        meta = sample_metadata_frame(path)
        classified = meta.apply(lambda r: classify_sample(acc, r), axis=1, result_type="expand")
        classified.columns = ["phenotype", "participant_id", "control_source"]
        df = pd.concat([meta, classified], axis=1)
        df.insert(0, "accession", acc)
        df.insert(1, "role", row.role)
        df["primary_include"] = make_primary_include(acc, df)
        df["case"] = df["phenotype"].eq("tumor").astype("Int64")
        all_frames.append(df)

        n_tumor = int((df["phenotype"] == "tumor").sum())
        n_normal = int(df["phenotype"].isin(["normal_adjacent", "normal_healthy", "normal_other"]).sum())
        n_ambig = int((df["phenotype"] == "ambiguous").sum())
        n_excl = int((df["phenotype"] == "exclude").sum())
        pri = df[df["primary_include"]]
        p_tumor = int((pri["phenotype"] == "tumor").sum())
        p_normal = int(pri["phenotype"].isin(["normal_adjacent", "normal_healthy", "normal_other"]).sum())
        summaries.append({
            "accession": acc,
            "role": row.role,
            "n_series_samples": len(df),
            "tumor_all": n_tumor,
            "normal_all": n_normal,
            "ambiguous": n_ambig,
            "excluded": n_excl,
            "primary_tumor": p_tumor,
            "primary_normal": p_normal,
            "expected_tumor_all": row.expected_tumor_all,
            "expected_normal_all": row.expected_normal_all,
            "expected_primary_tumor": row.primary_tumor,
            "expected_primary_normal": row.primary_normal,
        })

        for field, observed, expected in [
            ("tumor_all", n_tumor, num(row.expected_tumor_all)),
            ("normal_all", n_normal, num(row.expected_normal_all)),
            ("primary_tumor", p_tumor, num(row.primary_tumor)),
            ("primary_normal", p_normal, num(row.primary_normal)),
        ]:
            if expected is not None and observed != expected:
                failures.append({"accession": acc, "check": field, "observed": observed, "expected": expected})

        if row.paired_primary == "1":
            z = df[df["phenotype"].isin(["tumor", "normal_adjacent"])]
            counts = z.groupby("participant_id")["phenotype"].nunique()
            n_complete = int((counts == 2).sum())
            pair_rows.append({
                "accession": acc,
                "complete_pairs": n_complete,
                "primary_included_samples": int(df["primary_include"].sum()),
            })

    full = pd.concat(all_frames, ignore_index=True)
    full.to_csv(RESULTS / "manifests" / "sample_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(RESULTS / "qc" / "cohort_sample_counts.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(RESULTS / "qc" / "pairing_summary.csv", index=False)
    pd.DataFrame(failures).to_csv(RESULTS / "qc" / "inclusion_failures.csv", index=False)

    if failures:
        f = pd.DataFrame(failures)
        print(f.to_string(index=False))
        raise SystemExit("Frozen inclusion-count audit FAILED. Resolve metadata rules before preprocessing.")
    print("Frozen inclusion-count audit PASSED.")


if __name__ == "__main__":
    main()
