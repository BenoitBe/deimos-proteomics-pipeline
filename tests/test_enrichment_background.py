"""Univers de l'enrichissement (enrichment_background) : protéines quantifiées par
défaut pour gProfiler (domain_scope='custom_annotated'), STRING et les clusters.
Aucun appel réseau : l'API gProfiler et STRING sont simulées."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

go = pytest.importorskip("go_enrichment")

PARAMS = {"volcano_use_padj": False, "volcano_p_thresh": 0.05,
          "volcano_ratio_min": 1.5, "volcano_lfc_min": float(np.log2(1.5))}


def _comparison(n=40, n_sig=8):
    names = [f"P{i:03d};Q{i:03d}" for i in range(n)]
    p = np.r_[np.full(n_sig, 0.001), np.full(n - n_sig, 0.5)]
    return pd.DataFrame({"name": names, "C_vs_D_p.val": p, "C_vs_D_p.adj": p,
                         "C_vs_D_diff": np.r_[np.full(n_sig, 2.0), np.zeros(n - n_sig)]})


class _FakeResponse:
    status_code = 200

    def __init__(self, query):
        self._q = query

    def json(self):
        q0 = self._q[0]
        return {"meta": {"query_metadata": {"queries": {"query_1": list(self._q)}},
                         "genes_metadata": {"query": {"query_1": {
                             "mapping": {q0: ["ENSG1"]}, "ensgs": ["ENSG1"]}}}},
                "result": [{"source": "GO:BP", "native": "GO:0000001", "name": "term",
                            "p_value": 0.001, "significant": True, "description": "d",
                            "term_size": 10, "query_size": len(self._q),
                            "intersection_size": 1, "effective_domain_size": 30,
                            "precision": 0.1, "recall": 0.1, "query": "query_1",
                            "parents": [], "intersections": [["IDA"]]}]}


@pytest.fixture
def gprofiler_api(monkeypatch):
    gpm = pytest.importorskip("gprofiler.gprofiler")
    sent = []
    real = gpm.requests

    class _Requests:
        def __getattr__(self, k):
            return getattr(real, k)

        def post(self, url, json=None, **kw):
            sent.append(json)
            return _FakeResponse(json["query"])

    monkeypatch.setattr(gpm, "requests", _Requests())
    return gpm, sent


class TestRunGost:
    def test_custom_annotated_with_background(self, gprofiler_api):
        gpm, sent = gprofiler_api
        before = gpm.requests
        res = go.run_gost(["P001", "P002"], "hsapiens", background=["P003", "P001"])
        assert res is not None and len(res) == 1
        assert sent[-1]["domain_scope"] == "custom_annotated"
        assert set(sent[-1]["background"]) == {"P001", "P002", "P003"}  # requête ⊂ fond
        assert gpm.requests is before                                 # patch retiré

    def test_known_without_background(self, gprofiler_api):
        _, sent = gprofiler_api
        go.run_gost(["P001", "P002"], "hsapiens")
        assert sent[-1]["domain_scope"] == "known" and not sent[-1]["background"]


class TestUniverse:
    def _capture(self, monkeypatch):
        calls = []

        def fake_gost(query, organism, background=None):
            calls.append({"query": list(query), "background": background})
            return None

        monkeypatch.setattr(go, "run_gost", fake_gost)
        return calls

    def test_contrast_background_is_the_tested_proteins(self, monkeypatch, tmp_path):
        calls = self._capture(monkeypatch)
        df = _comparison()
        go.run_go_enrichment(df, ["C_vs_D"], PARAMS,
                             {"organism": "hsapiens", "backend": "gprofiler"}, str(tmp_path))
        assert calls and calls[0]["background"] == [f"P{i:03d}" for i in range(40)]
        calls.clear()
        go.run_go_enrichment(df, ["C_vs_D"], PARAMS,
                             {"organism": "hsapiens", "backend": "gprofiler",
                              "background_mode": "genome"}, str(tmp_path))
        assert calls[0]["background"] is None

    def test_orthologue_background(self, monkeypatch, tmp_path):
        calls = self._capture(monkeypatch)
        orth = {f"P{i:03d}": f"AT{i:03d}" for i in range(0, 40, 2)} | \
               {f"P{i:03d}": f"AT{i:03d}" for i in range(8)}
        go.run_go_enrichment(_comparison(), ["C_vs_D"], PARAMS,
                             {"organism": "athaliana", "backend": "gprofiler",
                              "ortholog_map": orth}, str(tmp_path))
        assert sorted(calls[0]["background"]) == sorted(set(orth.values()))

    def test_clusters_use_the_full_universe(self, monkeypatch, tmp_path):
        calls = self._capture(monkeypatch)
        clusters = [(f"P{i:03d};X", "Cluster_1") for i in range(8)]
        universe = [f"P{i:03d};Q" for i in range(40)]
        go.run_go_enrichment_clusters(clusters, PARAMS, {"organism": "hsapiens"},
                                      str(tmp_path), universe=universe)
        assert len(calls[0]["background"]) == 40          # pas seulement le cluster

    def test_string_clusters_resolve_the_universe(self, monkeypatch, tmp_path):
        seen = {}

        def fake_prepare(ids, taxid, caller, mode="genome"):
            seen["ids"], seen["mode"] = list(ids), mode
            return None, None

        monkeypatch.setattr(go, "_prepare_string", fake_prepare)
        clusters = [(f"P{i:03d}", "Cluster_1") for i in range(8)]
        go.run_go_enrichment_clusters(
            clusters, PARAMS, {"backend": "string", "species_taxid": 9606},
            str(tmp_path), universe=[f"P{i:03d}" for i in range(40)])
        assert len(seen["ids"]) == 40 and seen["mode"] == "quantified"

    def test_explicit_string_background_wins(self, monkeypatch, tmp_path):
        seen = {}

        def fake_prepare(ids, taxid, caller, mode="genome"):
            seen["mode"] = mode
            return None, None

        monkeypatch.setattr(go, "_prepare_string", fake_prepare)
        go.run_go_enrichment(_comparison(), ["C_vs_D"], PARAMS,
                             {"backend": "string", "species_taxid": 9606,
                              "background_mode": "quantified",
                              "string_background": "genome"}, str(tmp_path))
        assert seen["mode"] == "genome"


class TestConfig:
    def test_defaults_and_validation(self):
        import config
        assert config._DEFAULTS["enrichment_background"] == "quantified"
        assert config._DEFAULTS["string_background"] is None
        base = dict(config._DEFAULTS, tsv_path="x", design_path="y")
        assert not [e for e in config._validate(base) if "background" in e]
        assert any("enrichment_background" in e for e in
                   config._validate(dict(base, enrichment_background="proteome")))
        assert any("string_background" in e for e in
                   config._validate(dict(base, string_background="all")))
