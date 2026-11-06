"""Deimos -> InSight : localisation du clone, ligne de commande, validation des
clés YAML, et contrat du lancement (processus séparé, jamais bloquant).
Un faux clone InSight suffit : aucune dépendance au vrai dépôt."""
import os
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import insight_bridge as ib                            # noqa: E402


def _fake_clone(tmp, body):
    d = tmp / "InSight"
    (d / "insight").mkdir(parents=True)
    (d / "insight" / "__init__.py").write_text('__version__ = "9.9.9"\n')
    (d / "insight.py").write_text(textwrap.dedent(body))
    return d


OK_BODY = """
    import os, sys
    out = sys.argv[sys.argv.index("--out") + 1]
    os.makedirs(out, exist_ok=True)
    print("args:", " ".join(sys.argv[1:]))
    print("  [WARN] InSight: themes skipped (no GO tables)")
    open(os.path.join(out, "insight_report.html"), "w").write("<html></html>")
"""


def test_locate_and_version(tmp_path, monkeypatch):
    monkeypatch.delenv("INSIGHT_HOME", raising=False)
    assert ib.locate_insight({"insight_path": str(tmp_path / "nowhere")}) in (
        None, ib.locate_insight({}))
    d = _fake_clone(tmp_path, OK_BODY)
    assert ib.locate_insight({"insight_path": str(d)}) == str(d)
    monkeypatch.setenv("INSIGHT_HOME", str(d))
    assert ib.locate_insight({}) == str(d)
    assert ib.insight_version(str(d)) == "9.9.9"


def test_argv_from_options(tmp_path):
    p = {"insight_options": {"k": 5, "max_volcanos": 2, "volcanos": ["A_vs_B", "C_vs_D"],
                             "categorical": True, "order_regex": r"T(\d+)",
                             "title": "My project", "reference_fasta": "ref.fasta"}}
    a = ib.build_insight_argv("/x/InSight", str(tmp_path), str(tmp_path / "r.xlsx"), p)
    assert a[1].endswith("insight.py") and a[2] == str(tmp_path)
    assert a[a.index("--out") + 1] == os.path.join(str(tmp_path), "insight")
    assert a[a.index("--k") + 1] == "5" and a[a.index("--max-volcanos") + 1] == "2"
    assert a.count("--volcano") == 2 and "--categorical" in a
    assert a[a.index("--order-regex") + 1] == r"T(\d+)"
    assert a[a.index("--title") + 1] == "My project"
    assert os.path.isabs(a[a.index("--reference-fasta") + 1])
    # sans options : seulement le dossier, --out et --xlsx
    a = ib.build_insight_argv("/x/InSight", str(tmp_path), None, {})
    assert a[3:] == ["--out", os.path.join(str(tmp_path), "insight")]


def test_validation():
    assert ib.validate_insight_params({}) == []
    assert ib.validate_insight_params({"insight_options": {"k": 6, "volcanos": ["a"]}}) == []
    assert ib.validate_insight_params({"insight_options": "k=6"})
    assert ib.validate_insight_params({"insight_options": {"kk": 6}})
    assert ib.validate_insight_params({"insight_options": {"k": 1}})
    assert ib.validate_insight_params({"insight_options": {"volcanos": "A_vs_B"}})
    assert ib.validate_insight_params({"insight_timeout": -5})


def test_run_step_success(tmp_path, capsys):
    d = _fake_clone(tmp_path, OK_BODY)
    out = tmp_path / "proteomics_output"
    out.mkdir()
    html = ib.run_insight_step(str(out / "r.xlsx"), {"insight_path": str(d),
                               "insight_options": {"title": "T T"}}, str(out))
    assert html == str(out / "insight" / "insight_report.html")
    log = (out / "insight" / "insight_log.txt").read_text()
    assert "--title" in log and "--xlsx" in log
    o = capsys.readouterr().out
    assert "[OK] InSight 9.9.9" in o and "themes skipped" in o     # [WARN] relayés


def test_run_step_failure_never_raises(tmp_path, capsys):
    d = _fake_clone(tmp_path, "import sys\nprint('boom')\nsys.exit(3)\n")
    out = tmp_path / "o"
    out.mkdir()
    assert ib.run_insight_step(None, {"insight_path": str(d)}, str(out)) is None
    o = capsys.readouterr().out
    assert "exit code 3" in o and "boom" in o


def test_run_step_not_found(tmp_path, capsys, monkeypatch):
    monkeypatch.delenv("INSIGHT_HOME", raising=False)
    monkeypatch.setattr(ib, "insight_search_paths", lambda p: [str(tmp_path / "none")])
    assert ib.run_insight_step(None, {}, str(tmp_path)) is None
    assert "[SKIP] InSight not found" in capsys.readouterr().out


def test_deimos_keys_and_hook():
    import config
    for k in ("make_insight", "insight_path", "insight_options", "insight_timeout"):
        assert k in config._DEFAULTS, k
    src = (ROOT / "deimos.py").read_text(encoding="utf-8")
    assert "run_insight_step(output_name, params, out_dir)" in src
    # après Pathfinder, avant le nettoyage
    assert src.index("15bis. Pathfinder") < src.index("15ter. InSight") < src.index("16. Nettoyage")
