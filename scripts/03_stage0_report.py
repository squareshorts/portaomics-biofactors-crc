from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from common import CONFIG, REPORTS, RESULTS


def table_md(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    cols = [str(c) for c in df.columns]

    def esc(v):
        if pd.isna(v):
            return ""
        return str(v).replace("|", "\\|").replace("\n", " ")

    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"

    rows = []
    for row in df.itertuples(index=False, name=None):
        rows.append("| " + " | ".join(esc(v) for v in row) + " |")

    return "\n".join([header, sep] + rows)


def main() -> None:
    counts = pd.read_csv(RESULTS / "qc" / "cohort_sample_counts.csv")
    pairs = pd.read_csv(RESULTS / "qc" / "pairing_summary.csv") if (RESULTS / "qc" / "pairing_summary.csv").exists() else pd.DataFrame()
    dims = pd.read_csv(RESULTS / "qc" / "expression_dimensions.csv")
    maps = pd.read_csv(RESULTS / "qc" / "platform_mapping_summary.csv")
    vals = pd.read_csv(RESULTS / "qc" / "value_summary.csv")
    failures_path = RESULTS / "qc" / "inclusion_failures.csv"
    try:
        failures = pd.read_csv(failures_path) if failures_path.exists() else pd.DataFrame()
    except pd.errors.EmptyDataError:
        failures = pd.DataFrame()
    lock = json.loads((CONFIG / "analysis_lock.json").read_text(encoding="utf-8"))
    inter = (RESULTS / "qc" / "gene_intersection.txt").read_text(encoding="utf-8").splitlines()

    status = "PASS" if failures.empty else "FAIL"
    lines = [
        "# Stage-0 data-integrity and preprocessing audit",
        "",
        f"Protocol lock: {lock['lock_version']} ({lock['lock_date']})",
        "",
        f"Overall inclusion-count status: **{status}**",
        "",
        f"Genes measurable after harmonization in all nine transcriptomic cohorts: **{len(inter):,}**",
        "",
        "No disease-association statistic, differential-expression result, AUC, classifier, PAFS score, or external predictive performance was computed in this stage.",
        "",
        "## Sample-count audit",
        "",
        table_md(counts),
        "",
        "## Pairing audit",
        "",
        table_md(pairs),
        "",
        "## Expression dimensions",
        "",
        table_md(dims),
        "",
        "## Probe-to-gene mapping",
        "",
        table_md(maps),
        "",
        "## Value-range / missingness audit",
        "",
        table_md(vals),
        "",
        "## Frozen-analysis gate",
        "",
        "`classifier_analysis_authorized` remains **false**. Stage 1 must not be run until this audit has been reviewed and any technical discrepancies have been resolved without reference to downstream predictive performance.",
    ]
    if not failures.empty:
        lines += ["", "## Inclusion failures", "", table_md(failures)]
    out = REPORTS / "stage0_audit.md"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
