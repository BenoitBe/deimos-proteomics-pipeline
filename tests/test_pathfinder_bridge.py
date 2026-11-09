"""Deimos -> Pathfinder: the _manifest contract and the Pathfinder launch step.

Fast tests (no Deimos run, no network). The round trip through Pathfinder's own
reader runs only when a Pathfinder clone is found (../Pathfinder or $PATHFINDER_HOME).
"""
import os
import sys
import textwrap
from pathlib import Path

import openpyxl
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pathfinder_bridge as pb                         # noqa: E402
from go_enrichment import export_go_sheets             # noqa: E402

# The five naming situations of INTEGRATION_PATHFINDER.md §2, where Pathfinder's
# name-based reconstruction found 1/6, 0/12, 3/3, 3/3 and 0/6 sheets.
_STRAINS = ["K12", "M18", "JIM8777", "CCHSS3", "ATCC25975", "NCTC8618",
            "SK126", "PS4", "HSISS4", "DSM20560", "YU10", "CIRMBP1243"]
NAME_SETS = {
    "ko_wt_stress": ["KO_vs_WT_stress", "KO_vs_KO_stress", "KO_vs_WT",
                     "WT_stress_vs_KO_stress", "WT_stress_vs_WT", "KO_stress_vs_WT"],
    "intra_strain": [f"Live_{s}_vs_Heat_killed_{s}" for s in _STRAINS],
    "long_prefix": ["X18_009IG01di6h_vs_X18_009IG01Ctrl",
                    "X18_009IG01di24h_vs_X18_009IG01Ctrl",
                    "X18_009IG01di24h_vs_X18_009IG01di6h"],
    "nested": ["EPHPS7_vs_HPS7", "IPHPS7_vs_HPS7", "EPHPS7_vs_IPHPS7"],
    "long_close": ["AlphaTreatment_day01_vs_BetaTreatmentX_day01",
                   "AlphaTreatment_day14_vs_BetaTreatmentX_day14",
                   "AlphaTreatment_day14_vs_AlphaTreatment_day01",
                   "BetaTreatmentX_day14_vs_BetaTreatmentX_day01",
                   "AlphaTreatment_day01_vs_BetaTreatmentX_day14",
                   "AlphaTreatment_day14_vs_BetaTreatmentX_day01"],
}

DESIGN = pd.DataFrame({"label": [f"S{i}" for i in range(6)],
                       "condition": ["A"] * 3 + ["B"] * 3,
                       "replicate": [1, 2, 3] * 2,
                       "batch": ["b1", "b2", "b1", "b2", "b1", "b2"]})

BASE_PARAMS = {"volcano_use_padj": False, "volcano_p_thresh": 0.05,
               "volcano_ratio_min": 1.5, "fdr_global": False,
               "impute_method": "qrilc", "seed": 42}


def _go_results(contrasts):
    """Same shape as run_go_enrichment(): one table per contrast, gProfiler columns,
    a distinct term count per contrast so a mis-assignment would show."""
    out = {}
    for k, c in enumerate(contrasts):
        n = 3 + k
        out[c] = {"table": pd.DataFrame({
            "source": ["GO:BP"] * n,
            "term_id": [f"GO:{1000000 + 100 * k + i}" for i in range(n)],
            "term_name": [f"term_{k}_{i}" for i in range(n)],
            "p_value": [1e-4] * n,
            "intersection": ["P1,P2,P3"] * n})}
    return out


def _workbook_with_go(path, contrasts, params=None):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Differential_Expression"
    ws.append(["Protein.Group"] + [f"{c}_{s}" for c in contrasts
                                   for s in ("diff", "p.val", "p.adj")])
    go_map = export_go_sheets(wb, _go_results(contrasts), None, lambda *a, **k: None)
    pb.write_manifest(wb, {**BASE_PARAMS, **(params or {})}, contrasts, go_map, DESIGN,
                      deimos_commit="test")
    wb.save(path)
    return go_map


def _manifest(path):
    df = pd.read_excel(path, sheet_name="_manifest", dtype=str).fillna("")
    out = {}
    for s, k, v in df[["section", "key", "value"]].itertuples(index=False):
        out.setdefault(s, {})[k] = v
    return out


# --------------------------------------------------------------------------- #
# Manifest                                                                     #
# --------------------------------------------------------------------------- #
class TestManifest:
    @pytest.mark.parametrize("name", sorted(NAME_SETS))
    def test_every_go_sheet_is_declared_with_its_exact_contrast(self, tmp_path, name):
        contrasts = NAME_SETS[name]
        path = tmp_path / f"{name}.xlsx"
        go_map = _workbook_with_go(path, contrasts)
        sheets = openpyxl.load_workbook(path, read_only=True).sheetnames
        go_tabs = [s for s in sheets if s.startswith("GO_")]
        assert len(go_tabs) == len(contrasts)
        declared = _manifest(path)["go_sheet"]
        assert set(declared) == set(go_tabs)
        for c, tab in go_map.items():
            assert declared[tab] == c

    def test_sheet_is_last_and_hidden(self, tmp_path):
        path = tmp_path / "wb.xlsx"
        _workbook_with_go(path, NAME_SETS["nested"])
        wb = openpyxl.load_workbook(path)
        assert wb.sheetnames[-1] == "_manifest"
        assert wb["_manifest"].sheet_state == "hidden"

    def test_rewriting_does_not_duplicate(self, tmp_path):
        wb = openpyxl.Workbook()
        pb.write_manifest(wb, BASE_PARAMS, ["A_vs_B"], {}, DESIGN, deimos_commit="x")
        pb.write_manifest(wb, BASE_PARAMS, ["A_vs_B"], {}, DESIGN, deimos_commit="x")
        assert wb.sheetnames.count("_manifest") == 1

    def test_core_facts(self):
        rows = pb.manifest_rows({**BASE_PARAMS, "fdr_global": True}, ["KO_vs_WT"], {},
                                DESIGN, deimos_commit="abc1234")
        m = {}
        for s, k, v in rows:
            m.setdefault(s, {})[k] = v
        assert m["schema"] == {"name": "deimos-results", "version": pb.SCHEMA_VERSION}
        assert m["run"]["deimos_commit"] == "abc1234"
        assert m["stats"]["sheet"] == "Differential_Expression"
        assert m["stats"]["p_adjust"] == "BH global"
        assert m["stats"]["volcano_rule"] == "p.value < 0.05 & ratio >= 1.5"
        assert m["matrix"] == {"sheet": "Log2_Impute", "imputed": "true",
                               "imputation": "QRILC, seed 42"}
        assert m["contrast"] == {"KO_vs_WT": "KO|WT"}
        assert "go" not in m                         # GO did not run
        assert m["design"]["factors"] == "batch"     # replicate is not a factor

    @pytest.mark.parametrize("params,expected", [
        ({"go_backend": "string", "go_species_taxid": "1304", "go_string_mode": "ortholog"},
         {"taxid": "1304", "string_mode": "ortholog"}),
        ({"go_backend": "string", "go_species_taxid": "auto"}, None),
        ({"go_backend": "gprofiler", "go_organism": "bnapus"}, {"gprofiler": "bnapus"}),
        ({"go_backend": "gprofiler", "go_organism": None}, None),
    ])
    def test_organism(self, params, expected):
        rows = pb.manifest_rows({**BASE_PARAMS, **params}, [], {}, DESIGN, deimos_commit="x")
        org = {k: v for s, k, v in rows if s == "organism"}
        assert org == (expected or {})

    def test_perseverance_declares_no_organism(self):
        """GO on the REFERENCE species: the measured organism is not that one."""
        params = {**BASE_PARAMS, "go_organism": "athaliana",
                  "perseverance_reference_organism": "athaliana"}
        rows = pb.manifest_rows(params, ["A_vs_B"], {"A_vs_B": "GO_AvB"}, DESIGN,
                                perseverance=True, deimos_commit="x")
        assert not [r for r in rows if r[0] == "organism"]
        assert ("go", "backend", "perseverance") in rows
        assert ("go", "reference_organism", "athaliana") in rows
        assert pb.pathfinder_org(params, perseverance=True) is None


# --------------------------------------------------------------------------- #
# Round trip through Pathfinder's own reader                                   #
# --------------------------------------------------------------------------- #
PF_DIR = pb.locate_pathfinder({})


@pytest.mark.skipif(PF_DIR is None, reason="no Pathfinder clone (../Pathfinder, PATHFINDER_HOME)")
class TestPathfinderReadsTheManifest:
    @pytest.fixture(autouse=True)
    def _io_layer(self, monkeypatch):
        monkeypatch.syspath_prepend(PF_DIR)
        import io_layer
        if not hasattr(io_layer, "read_manifest"):
            pytest.skip("this Pathfinder predates the manifest (P1)")
        self.io_layer = io_layer

    @pytest.mark.parametrize("name", sorted(NAME_SETS))
    def test_all_go_sheets_mapped_via_manifest(self, tmp_path, name):
        contrasts = NAME_SETS[name]
        path = tmp_path / f"{name}.xlsx"
        _workbook_with_go(path, contrasts)
        go = self.io_layer.load_go(str(path), contrasts, verbose=False)
        assert go.unmatched == [] and go.duplicates == {}
        assert set(go.mapping_how.values()) == {"manifest"}
        terms = go.long.groupby("contrast").size().to_dict()
        assert terms == {c: 3 + k for k, c in enumerate(contrasts)}


# --------------------------------------------------------------------------- #
# Launch step (fake Pathfinder: fast, no dependency)                           #
# --------------------------------------------------------------------------- #
def _fake_pathfinder(tmp_path, body):
    d = tmp_path / "FakePathfinder"
    d.mkdir()
    (d / "pathfinder.py").write_text(
        "import sys, os, json, time\n__version__ = '9.9.9'\n" + textwrap.dedent(body),
        encoding="utf-8")
    return d


WRITE_REPORT = """
    a = sys.argv[1:]
    out = a[a.index('--outdir') + 1]
    with open(os.path.join(out, 'argv.json'), 'w') as f:
        json.dump({'argv': a, 'stdin_tty': sys.stdin.isatty(),
                   'stdin_eof': sys.stdin.read() == ''}, f)
    print('Pathfinder ⚠ ε² —')             # non-ASCII on a redirected stdout
    open(os.path.join(out, 'pathfinder_report.html'), 'w').write('<html></html>')
"""


class TestRunPathfinderStep:
    def test_not_found_skips(self, tmp_path, capsys, monkeypatch):
        monkeypatch.delenv("PATHFINDER_HOME", raising=False)
        monkeypatch.setattr(pb, "_HERE", str(tmp_path / "nowhere" / "deimos"))
        r = pb.run_pathfinder_step(tmp_path / "x.xlsx", DESIGN,
                                   {"pathfinder_path": str(tmp_path / "missing")}, str(tmp_path))
        assert r is None
        assert "[SKIP]" in capsys.readouterr().out

    def test_success_and_arguments(self, tmp_path):
        pf = _fake_pathfinder(tmp_path, WRITE_REPORT)
        xlsx = tmp_path / "ProteomicAnalysis_Results.xlsx"
        xlsx.write_bytes(b"")
        params = {"pathfinder_path": str(pf), "go_backend": "string",
                  "go_species_taxid": "1304", "pathfinder_options": {"pthr": 0.01}}
        html = pb.run_pathfinder_step(xlsx, DESIGN, params, str(tmp_path), go_ran=False)
        assert html and html.endswith(os.path.join("pathfinder", "pathfinder_report.html"))
        import json
        info = json.loads((tmp_path / "pathfinder" / "argv.json").read_text())
        a = info["argv"]
        assert a[a.index("--go") + 1] == "none"           # never auto-detect a GO file
        assert a[a.index("--stats") + 1] == str(xlsx.resolve())
        assert "--no-online" in a and "--org" not in a    # Publications off by default
        assert a[a.index("--pthr") + 1] == "0.01"
        assert info["stdin_eof"] and not info["stdin_tty"]
        design = pd.read_csv(tmp_path / "pathfinder" / "design_pathfinder.csv")
        assert list(design.columns) == ["label", "condition", "batch"]
        log = (tmp_path / "pathfinder" / "pathfinder_log.txt").read_text(encoding="utf-8")
        assert "⚠ ε² —" in log

    def test_go_sheets_are_read_from_the_workbook(self, tmp_path):
        pf = _fake_pathfinder(tmp_path, WRITE_REPORT)
        xlsx = tmp_path / "wb.xlsx"
        xlsx.write_bytes(b"")
        pb.run_pathfinder_step(xlsx, DESIGN, {"pathfinder_path": str(pf)}, str(tmp_path),
                               go_ran=True)
        import json
        a = json.loads((tmp_path / "pathfinder" / "argv.json").read_text())["argv"]
        assert a[a.index("--go") + 1] == str(xlsx.resolve())

    @pytest.mark.parametrize("params,perseverance,org", [
        ({"go_backend": "string", "go_species_taxid": "1304"}, False, "1304"),
        ({"go_backend": "gprofiler", "go_organism": "bnapus"}, False, "bnapus"),
        ({"go_backend": "gprofiler", "go_organism": "athaliana"}, True, None),
        ({"go_backend": "gprofiler"}, False, None),
    ])
    def test_publications_species(self, tmp_path, params, perseverance, org):
        argv = pb.build_pathfinder_argv("/pf", "/x.xlsx", "/d.csv", "/o",
                                        {**params, "pathfinder_publications": True},
                                        perseverance=perseverance)
        if org:
            assert argv[argv.index("--org") + 1] == org and "--no-online" not in argv
        else:
            assert "--no-online" in argv and "--org" not in argv

    def test_failure_warns_and_returns_none(self, tmp_path, capsys):
        pf = _fake_pathfinder(tmp_path, "print('boom: no _p.val columns'); sys.exit(2)\n")
        r = pb.run_pathfinder_step(tmp_path / "wb.xlsx", DESIGN, {"pathfinder_path": str(pf)},
                                   str(tmp_path))
        out = capsys.readouterr().out
        assert r is None and "[WARN]" in out and "exit code 2" in out and "boom" in out

    def test_crash_without_report_is_a_failure(self, tmp_path, capsys):
        pf = _fake_pathfinder(tmp_path, "pass\n")           # exit 0, no report
        assert pb.run_pathfinder_step(tmp_path / "wb.xlsx", DESIGN,
                                      {"pathfinder_path": str(pf)}, str(tmp_path)) is None
        assert "[WARN]" in capsys.readouterr().out

    def test_timeout(self, tmp_path, capsys):
        pf = _fake_pathfinder(tmp_path, "time.sleep(30)\n")
        r = pb.run_pathfinder_step(tmp_path / "wb.xlsx", DESIGN,
                                   {"pathfinder_path": str(pf), "pathfinder_timeout": 1},
                                   str(tmp_path))
        out = capsys.readouterr().out
        assert r is None and "[WARN]" in out and "1 s" in out

    def test_never_raises(self, tmp_path, capsys):
        pf = _fake_pathfinder(tmp_path, WRITE_REPORT)
        r = pb.run_pathfinder_step(tmp_path / "wb.xlsx", "not a DataFrame",
                                   {"pathfinder_path": str(pf)}, str(tmp_path))
        assert r is None and "[WARN]" in capsys.readouterr().out

    def test_version_string(self, tmp_path):
        pf = _fake_pathfinder(tmp_path, "pass\n")
        assert pb.pathfinder_version(str(pf)).startswith("9.9.9")


class TestValidation:
    def test_defaults_are_valid(self):
        assert pb.validate_pathfinder_params({"pathfinder_options": {},
                                              "pathfinder_timeout": 1800}) == []

    @pytest.mark.parametrize("params", [
        {"pathfinder_options": "pthr=0.05"},
        {"pathfinder_options": {"pthr": 2}},
        {"pathfinder_options": {"p_thr": 0.05}},
        {"pathfinder_timeout": -1},
        {"pathfinder_timeout": "soon"},
    ])
    def test_invalid(self, params):
        assert pb.validate_pathfinder_params(params)


def test_reading_options_reach_pathfinder_13(tmp_path):
    """Pathfinder >= 1.3.0 porte la lecture biologique (ex-InSight) : options de
    pathfinder_options, sinon d'insight_options ; jamais passées à un Pathfinder 1.2."""
    import pathfinder_bridge as pb
    params = {"pathfinder_options": {"k": 5, "groups_from": "gsea", "prog_boot": 0},
              "insight_options": {"k": 7, "volcanos": ["A_vs_B"], "categorical": True}}
    argv = pb.build_pathfinder_argv("/pf", "/x.xlsx", "/d.csv", "/o", params, reading=True)
    s = " ".join(argv)
    assert "--k 5" in s and "--k 7" not in s and "--volcano A_vs_B" in s
    assert "--categorical" in s and "--groups-from gsea" in s and "--prog-boot 0" in s
    assert "--k" not in pb.build_pathfinder_argv("/pf", "/x.xlsx", "/d.csv", "/o", params)
    old = tmp_path / "pf12"
    old.mkdir()
    (old / "pathfinder.py").write_text("ap.add_argument('--jaccard')", encoding="utf-8")
    new = tmp_path / "pf13"
    new.mkdir()
    (new / "pathfinder.py").write_text("rd.add_argument('--no-reading')", encoding="utf-8")
    assert not pb.pathfinder_has_reading(str(old)) and pb.pathfinder_has_reading(str(new))
    assert pb.reading_argv({"pathfinder_options": {"no_reading": True, "k": 3}}) == ["--no-reading"]
    assert pb.validate_pathfinder_params({"pathfinder_options": {"groups_from": "x", "k": 1}})
    assert not pb.validate_pathfinder_params(params)
    rx = r"(?P<genotype>[^_]+)_T(?P<time>\d+)"
    assert pb.reading_argv({"pathfinder_options": {"factor_regex": rx}}) == ["--factor-regex", rx]
    assert not pb.validate_pathfinder_params({"pathfinder_options": {"factor_regex": rx}})
    assert pb.validate_pathfinder_params({"pathfinder_options": {"factor_regex": "(a"}})
    assert pb.validate_pathfinder_params({"pathfinder_options": {"factor_regex": "T(\\d+)"}})


def test_make_insight_is_absorbed_by_pathfinder_13(tmp_path, monkeypatch, capsys):
    """Un seul rapport HTML : avec un Pathfinder >= 1.3.0, make_insight lance Pathfinder
    (onglet Biological reading) et plus InSight ; sans lui, InSight reste à part."""
    import pathfinder_bridge as pb
    new = tmp_path / "pf13"
    new.mkdir()
    (new / "pathfinder.py").write_text("rd.add_argument('--no-reading')", encoding="utf-8")
    old = tmp_path / "pf12"
    old.mkdir()
    (old / "pathfinder.py").write_text("ap.add_argument('--jaccard')", encoding="utf-8")
    monkeypatch.setattr(pb, "locate_pathfinder", lambda params: str(new))
    p = {"make_insight": True, "make_pathfinder": False, "insight_options": {"k": 5}}
    assert pb.resolve_reading(p) == "pathfinder"
    assert p["make_pathfinder"] is True and p["make_insight"] is False
    assert "no separate insight_report.html" in capsys.readouterr().out
    assert "--k 5" in " ".join(pb.reading_argv(p))                 # options d'InSight reprises
    p = {"make_insight": True, "make_pathfinder": True}
    assert pb.resolve_reading(p) == "pathfinder" and not p["make_insight"]
    assert pb.resolve_reading({"make_insight": False}) is None
    monkeypatch.setattr(pb, "locate_pathfinder", lambda params: str(old))
    p = {"make_insight": True, "make_pathfinder": True}
    assert pb.resolve_reading(p) == "insight" and p["make_insight"] is True
    monkeypatch.setattr(pb, "locate_pathfinder", lambda params: None)
    p = {"make_insight": True}
    assert pb.resolve_reading(p) == "insight" and not p.get("make_pathfinder")
