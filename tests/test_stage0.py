from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from common import clean_single_gene_symbol, geo_bucket, rank_normal_scores  # noqa: E402


def test_geo_bucket():
    assert geo_bucket("GSE9348") == "GSE9nnn"
    assert geo_bucket("GSE103512") == "GSE103nnn"
    assert geo_bucket("GPL96") == "GPLnnn"
    assert geo_bucket("GPL10558") == "GPL10nnn"


def test_gene_cleaning():
    assert clean_single_gene_symbol("TP53") == "TP53"
    assert clean_single_gene_symbol("TP53 /// TP53") == "TP53"
    assert clean_single_gene_symbol("TP53 /// WRAP53") == ""
    assert clean_single_gene_symbol("---") == ""


def test_rank_normal_scores_independent_columns():
    x = pd.DataFrame({"s1": [1.0, 2.0, 3.0], "s2": [30.0, 10.0, 20.0]}, index=["A", "B", "C"])
    z = rank_normal_scores(x)
    assert z.loc["A", "s1"] < z.loc["B", "s1"] < z.loc["C", "s1"]
    assert z.loc["B", "s2"] < z.loc["C", "s2"] < z.loc["A", "s2"]
    assert np.allclose(z.mean(axis=0).values, 0, atol=1e-12)
