from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

from common import (
    RAW_GEO,
    RAW_PLATFORM,
    RESULTS,
    load_registry,
    platform_annot_url,
    series_matrix_url,
    series_supp_url,
    sha256_file,
)


def download(url: str, path: Path, timeout: int = 90) -> None:
    if path.exists() and path.stat().st_size > 0:
        print(f"[keep] {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    print(f"[get]  {url}")
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with tmp.open("wb") as f, tqdm(total=total, unit="B", unit_scale=True, leave=False) as bar:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    f.write(chunk)
                    bar.update(len(chunk))
    tmp.replace(path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    reg = load_registry()
    expr = reg[reg["omic"].eq("transcriptome")].copy()

    todo = []
    for row in expr.itertuples(index=False):
        acc = row.accession
        if row.download_mode == "series_matrix":
            todo.append((series_matrix_url(acc), RAW_GEO / f"{acc}_series_matrix.txt.gz"))
        elif row.download_mode == "rnaseq_supplement":
            # Matrix is useful for sample metadata; supplement contains RPKM values.
            todo.append((series_matrix_url(acc), RAW_GEO / f"{acc}_series_matrix.txt.gz"))
            todo.append((series_supp_url(acc, f"{acc}_RAW.tar"), RAW_GEO / f"{acc}_RAW.tar"))

    platforms = sorted(set(expr.loc[expr["download_mode"].eq("series_matrix"), "platform"]))
    for gpl in platforms:
        todo.append((platform_annot_url(gpl), RAW_PLATFORM / f"{gpl}.annot.gz"))

    if args.dry_run:
        for url, path in todo:
            print(f"{path.relative_to(path.parents[3]) if len(path.parents) > 3 else path}\t{url}")
        return

    for url, path in todo:
        try:
            download(url, path)
        except requests.HTTPError as e:
            raise SystemExit(f"Download failed for {url}: {e}")
        time.sleep(0.2)

    rows = []
    for _, path in todo:
        if path.exists():
            rows.append({"path": str(path.relative_to(path.parents[3])), "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    out = RESULTS / "qc" / "file_checksums.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
