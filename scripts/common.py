from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
import tarfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"
RAW_GEO = ROOT / "data" / "raw" / "geo"
RAW_PLATFORM = ROOT / "data" / "raw" / "platform"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
REPORTS = ROOT / "reports"

for p in [RAW_GEO, RAW_PLATFORM, PROCESSED, RESULTS / "manifests", RESULTS / "qc", REPORTS]:
    p.mkdir(parents=True, exist_ok=True)


def load_registry() -> pd.DataFrame:
    return pd.read_csv(CONFIG / "datasets.csv", dtype=str).fillna("")


def geo_bucket(accession: str) -> str:
    m = re.fullmatch(r"([A-Z]+)(\d+)", accession.upper())
    if not m:
        raise ValueError(f"Invalid GEO accession: {accession}")
    prefix, digits = m.groups()
    n = int(digits)
    return f"{prefix}{n // 1000}nnn" if n >= 1000 else f"{prefix}nnn"


def series_matrix_url(accession: str) -> str:
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_bucket(accession)}/"
        f"{accession}/matrix/{accession}_series_matrix.txt.gz"
    )


def series_soft_url(accession: str) -> str:
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_bucket(accession)}/"
        f"{accession}/soft/{accession}_family.soft.gz"
    )


def series_supp_url(accession: str, filename: str) -> str:
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_bucket(accession)}/"
        f"{accession}/suppl/{filename}"
    )


def platform_annot_url(platform: str) -> str:
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{geo_bucket(platform)}/"
        f"{platform}/annot/{platform}.annot.gz"
    )


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def parse_series_matrix_header(path: Path) -> Dict[str, List[str]]:
    """Read !Sample_* header rows from a GEO Series Matrix file."""
    out: Dict[str, List[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))
            key = parts[0].lstrip("!")
            values = parts[1:]
            out.setdefault(key, []).append(values)
    return out


def sample_metadata_frame(path: Path) -> pd.DataFrame:
    h = parse_series_matrix_header(path)
    if "Sample_geo_accession" not in h:
        raise ValueError(f"No Sample_geo_accession in {path}")
    accessions = h["Sample_geo_accession"][0]
    n = len(accessions)
    rows = [{"sample_id": accessions[i]} for i in range(n)]
    for key, blocks in h.items():
        if key == "Sample_geo_accession":
            continue
        for bidx, values in enumerate(blocks, start=1):
            if len(values) != n:
                continue
            col = key if len(blocks) == 1 else f"{key}_{bidx}"
            for i, v in enumerate(values):
                rows[i][col] = v
    df = pd.DataFrame(rows)
    return df


def metadata_text(row: pd.Series) -> str:
    vals = []
    for k, v in row.items():
        if k == "sample_id" or pd.isna(v):
            continue
        vals.append(str(v))
    return " | ".join(vals).lower()


def extract_characteristic(row: pd.Series, key: str) -> str:
    target = key.lower().strip()
    for col, val in row.items():
        if not str(col).startswith("Sample_characteristics_ch1") or pd.isna(val):
            continue
        s = str(val)
        if ":" in s:
            k, v = s.split(":", 1)
            if k.lower().strip() == target:
                return v.strip()
    return ""


def classify_sample(accession: str, row: pd.Series) -> Tuple[str, str, str]:
    """Return phenotype, participant_id, control_source.

    phenotype is one of tumor, normal_adjacent, normal_healthy, normal_other, exclude, ambiguous.
    The rules are accession-specific and intentionally conservative.
    """
    title = str(row.get("Sample_title", ""))
    source = str(row.get("Sample_source_name_ch1", ""))
    text = metadata_text(row)


    if accession == "GSE41258":
        tissue = extract_characteristic(row, "tissue").lower().strip()
        pid = extract_characteristic(row, "patient id")

        if tissue == "primary tumor":
            return "tumor", pid, ""

        if tissue == "normal colon":
            return "normal_other", pid, "corresponding_normal"

        return "exclude", pid, ""

    if accession == "GSE39582":
        if "primary colorectal adenocarcinoma" in text:
            return "tumor", "", ""
        if "non-tumoral" in text or "non tumoral" in text or "normal mucosa" in text:
            return "normal_other", "", "non_tumoral_mucosa"
        return "ambiguous", "", ""

    if accession == "GSE9348":
        if "crc patient's tumor" in text or "crc patient" in text and "tumor" in text:
            return "tumor", "", ""
        if "healthy control" in text:
            return "normal_healthy", "", "healthy_control"
        return "ambiguous", "", ""

    if accession == "GSE23878":
        t = title.strip().upper()
        if re.match(r"^CC\d+", t):
            return "tumor", "", ""
        if re.match(r"^N\d+", t):
            return "normal_other", "", "non_cancerous_colorectal"
        if "colorectal cancer" in text or "crc" in text and "normal" not in text:
            return "tumor", "", ""
        if "normal" in text or "non-cancerous" in text:
            return "normal_other", "", "non_cancerous_colorectal"
        return "ambiguous", "", ""


    if accession == "GSE44861":
        tissue = extract_characteristic(row, "tissue").lower().strip()
        pid = extract_characteristic(row, "case_id")

        if tissue == "tumor":
            return "tumor", pid, ""

        if tissue in {
            "adjacent nontumor",
            "adjacent non-tumor",
            "adjacent noncancerous",
            "adjacent non-cancerous",
        }:
            return "normal_adjacent", pid, "adjacent_non_cancerous"

        return "ambiguous", pid, ""


    if accession == "GSE103512":
        cancer_type = extract_characteristic(
            row, "cancer type"
        ).upper().strip()

        normal_flag = extract_characteristic(
            row, "normal"
        ).lower().strip()

        pid = re.sub(
            r"_normal$",
            "",
            title.strip(),
            flags=re.I
        )

        if cancer_type != "CRC":
            return "exclude", pid, ""

        if normal_flag in {"yes", "y", "true", "1"}:
            return "normal_adjacent", pid, "matched_normal"

        if normal_flag in {"no", "n", "false", "0"}:
            return "tumor", pid, ""

        return "ambiguous", pid, ""

    if accession == "GSE156451":
        m = re.match(r"^([TN])(\d+)-RNA$", title.strip(), flags=re.I)
        if not m:
            return "ambiguous", "", ""
        typ, pid = m.groups()
        if typ.upper() == "T":
            return "tumor", pid, ""
        return "normal_adjacent", pid, "native_tissue"

    if accession == "GSE106582":
        tissue = extract_characteristic(row, "tissue").lower()
        pid = extract_characteristic(row, "patientid")
        if not pid:
            pid = re.split(r"[_-]", title.strip())[0]
        if tissue == "tumor" or "tissue tumor" in text:
            return "tumor", pid, ""
        if tissue == "mucosa" or "tissue mucosa" in text:
            return "normal_adjacent", pid, "adjacent_mucosa"
        return "ambiguous", pid, ""

    if accession == "GSE44076":
        m = re.search(r"from\s+([A-Za-z]\d+)\s+patient", title, flags=re.I)
        pid = m.group(1) if m else ""
        if title.lower().startswith("tumor sample"):
            return "tumor", pid, ""
        if title.lower().startswith("normal paired sample"):
            return "normal_adjacent", pid, "adjacent_mucosa"
        if "healthy" in title.lower() or "healthy" in text:
            return "normal_healthy", pid, "healthy_control"
        return "ambiguous", pid, ""

    if accession in {"GSE101764", "GSE131013"}:
        tissue = extract_characteristic(row, "tissue").lower()
        pid = extract_characteristic(row, "patientid")
        if not pid:
            m = re.search(r"from\s+([A-Za-z]\d+)\s+patient", title, flags=re.I)
            pid = m.group(1) if m else re.split(r"[_-]", title.strip())[0]
        if "tumor" in tissue or "tumor" in text:
            return "tumor", pid, ""
        if "healthy" in text:
            return "normal_healthy", pid, "healthy_control"
        if "mucosa" in tissue or "normal" in text:
            return "normal_adjacent", pid, "adjacent_mucosa"
        return "ambiguous", pid, ""

    return "ambiguous", "", ""


def make_primary_include(accession: str, df: pd.DataFrame) -> pd.Series:
    phen = df["phenotype"]
    if accession == "GSE106582":
        complete = (
            df[df["phenotype"].isin(["tumor", "normal_adjacent"])]
            .groupby("participant_id")["phenotype"]
            .nunique()
        )
        ids = set(complete[complete == 2].index)
        return df["participant_id"].isin(ids) & phen.isin(["tumor", "normal_adjacent"])
    if accession == "GSE44076":
        complete = (
            df[df["phenotype"].isin(["tumor", "normal_adjacent"])]
            .groupby("participant_id")["phenotype"]
            .nunique()
        )
        ids = set(complete[complete == 2].index)
        return df["participant_id"].isin(ids) & phen.isin(["tumor", "normal_adjacent"])
    if accession == "GSE156451":
        complete = (
            df[df["phenotype"].isin(["tumor", "normal_adjacent"])]
            .groupby("participant_id")["phenotype"]
            .nunique()
        )
        ids = set(complete[complete == 2].index)
        return df["participant_id"].isin(ids) & phen.isin(["tumor", "normal_adjacent"])
    return phen.isin(["tumor", "normal_adjacent", "normal_healthy", "normal_other"])


def read_series_expression(path: Path, use_samples: Optional[List[str]] = None) -> pd.DataFrame:
    """Read a GEO Series Matrix expression table as probes x samples."""
    skip = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            if line.startswith("!series_matrix_table_begin"):
                skip = i + 1
                break
        else:
            raise ValueError(f"No series_matrix_table_begin in {path}")
    df = pd.read_csv(path, sep="\t", compression="gzip", skiprows=skip, comment="!", dtype={0: str})
    first = df.columns[0]
    df = df.rename(columns={first: "probe_id"}).set_index("probe_id")
    df.columns = [str(c).strip('"') for c in df.columns]
    if use_samples is not None:
        keep = [x for x in use_samples if x in df.columns]
        df = df[keep]
    return df.apply(pd.to_numeric, errors="coerce")


def read_platform_annotation(path: Path) -> pd.DataFrame:
    """
    Read either an official GEO GPL .annot.gz file or our locally derived
    two-column GPL13667 annotation.

    GEO annotation files may contain non-tabular metadata before the real
    header, so locate the actual probe/gene-symbol header first.
    """

    def norm_col(x):
        return re.sub(
            r"\s+",
            " ",
            str(x).strip().strip('"').replace("_", " ").lower()
        )

    header_row = None
    detected_columns = None

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as f:

        for i, line in enumerate(f):
            raw = line.rstrip("\r\n")

            if not raw or raw.startswith("#"):
                continue

            fields = [
                x.strip().strip('"')
                for x in raw.split("\t")
            ]

            normalized = [norm_col(x) for x in fields]

            id_present = any(
                x in {
                    "id",
                    "id ref",
                    "probe id",
                    "probe set id",
                    "probeset id",
                    "reporter id",
                }
                for x in normalized
            )

            symbol_present = any(
                x in {
                    "gene symbol",
                    "symbol",
                    "gene symbols",
                }
                for x in normalized
            )

            if id_present and symbol_present:
                header_row = i
                detected_columns = fields
                break

    if header_row is None:
        raise ValueError(
            f"Could not locate real annotation-table header in {path}"
        )

    print(
        f"[annotation] {path.name}: "
        f"header line {header_row + 1}; "
        f"{len(detected_columns)} columns"
    )

    df = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        skiprows=header_row,
        dtype=str,
        low_memory=False,
        comment="#"
    )

    df.columns = [
        str(c).strip().strip('"')
        for c in df.columns
    ]

    normalized_map = {
        norm_col(c): c
        for c in df.columns
    }

    id_col = next(
        (
            normalized_map[x]
            for x in [
                "id",
                "id ref",
                "probe id",
                "probe set id",
                "probeset id",
                "reporter id",
            ]
            if x in normalized_map
        ),
        None
    )

    sym_col = next(
        (
            normalized_map[x]
            for x in [
                "gene symbol",
                "symbol",
                "gene symbols",
            ]
            if x in normalized_map
        ),
        None
    )

    if id_col is None or sym_col is None:
        raise ValueError(
            f"Could not identify probe/gene columns in {path}; "
            f"columns={list(df.columns)[:30]}"
        )

    out = (
        df[[id_col, sym_col]]
        .rename(
            columns={
                id_col: "probe_id",
                sym_col: "gene_symbol"
            }
        )
        .dropna()
    )

    out["probe_id"] = (
        out["probe_id"]
        .astype(str)
        .str.strip()
        .str.strip('"')
    )

    out["gene_symbol"] = (
        out["gene_symbol"]
        .astype(str)
        .str.strip()
        .str.strip('"')
    )

    out = out[
        (out["probe_id"] != "") &
        (out["gene_symbol"] != "")
    ]

    print(
        f"[annotation] {path.name}: "
        f"{len(out):,} probe rows loaded"
    )

    return out

def clean_single_gene_symbol(s: str) -> str:
    s = str(s).strip()
    if not s or s in {"---", "NA", "nan"}:
        return ""
    # GEO annotations use separators such as ///, //, ; for multi-mapping.
    if "///" in s or "//" in s or ";" in s or " | " in s:
        parts = [p.strip() for p in re.split(r"///|//|;|\s\|\s", s) if p.strip()]
        uniq = {p for p in parts if p not in {"---", "NA"}}
        if len(uniq) != 1:
            return ""
        s = next(iter(uniq))
    if not re.fullmatch(r"[A-Za-z0-9_.\-]+", s):
        return ""
    return s.upper()


def collapse_probes_to_genes(expr: pd.DataFrame, annot: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, int]]:
    a = annot.copy()
    a["gene_symbol_clean"] = a["gene_symbol"].map(clean_single_gene_symbol)
    a = a[a["gene_symbol_clean"] != ""].drop_duplicates("probe_id")
    common = expr.index.intersection(a["probe_id"])
    sub = expr.loc[common].copy()
    mapper = a.set_index("probe_id")["gene_symbol_clean"]
    sub["__gene__"] = mapper.reindex(common).values
    gene = sub.groupby("__gene__", sort=True).median(numeric_only=True)
    summary = {
        "n_input_probes": int(expr.shape[0]),
        "n_annotated_unambiguous_probes": int(len(common)),
        "n_gene_symbols": int(gene.shape[0]),
    }
    return gene, summary


def rank_normal_scores(gene_by_sample: pd.DataFrame) -> pd.DataFrame:
    """Independent within-sample rank -> N(0,1), genes x samples."""
    arr = gene_by_sample.to_numpy(dtype=float)
    out = np.full_like(arr, np.nan, dtype=float)
    for j in range(arr.shape[1]):
        x = arr[:, j]
        ok = np.isfinite(x)
        n = int(ok.sum())
        if n == 0:
            continue
        r = rankdata(x[ok], method="average")
        pct = (r - 0.5) / n
        out[ok, j] = norm.ppf(pct)
    return pd.DataFrame(out, index=gene_by_sample.index, columns=gene_by_sample.columns)


def parse_rnaseq_txt_bytes(data: bytes) -> pd.Series:
    """Parse one GSE156451 processed per-sample TXT.GZ into a gene->RPKM Series.

    The parser is defensive because GEO supplements sometimes vary in column labels.
    """
    raw = gzip.decompress(data)
    text = raw.decode("utf-8", errors="replace")
    frame = pd.read_csv(io.StringIO(text), sep="\t", comment="#")
    if frame.shape[1] < 2:
        frame = pd.read_csv(io.StringIO(text), sep=r"\s+", engine="python", comment="#")
    cols = {c.lower().strip(): c for c in frame.columns}
    gene_col = next((cols[k] for k in ["gene_name", "gene", "gene_symbol", "symbol", "geneid", "gene_id"] if k in cols), frame.columns[0])
    value_col = next((cols[k] for k in ["rpkm", "fpkm", "value", "expression"] if k in cols), None)
    if value_col is None:
        numeric_candidates = []
        for c in frame.columns[1:]:
            z = pd.to_numeric(frame[c], errors="coerce")
            numeric_candidates.append((z.notna().mean(), c))
        numeric_candidates.sort(reverse=True)
        if not numeric_candidates or numeric_candidates[0][0] < 0.8:
            raise ValueError(f"Could not identify expression column; columns={list(frame.columns)}")
        value_col = numeric_candidates[0][1]
    genes = frame[gene_col].astype(str).str.strip().str.upper()
    vals = pd.to_numeric(frame[value_col], errors="coerce")
    s = pd.Series(vals.values, index=genes).replace([np.inf, -np.inf], np.nan).dropna()
    s = s[~s.index.duplicated(keep=False)].copy()
    s = s[s.index.str.match(r"^[A-Z0-9_.\-]+$")]
    return s
