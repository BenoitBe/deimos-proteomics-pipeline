"""Configuration of the Pathfinder step: YAML keys, --pathfinder, interactive question."""
import builtins
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import config                                          # noqa: E402
import insight_bridge as ib                            # noqa: E402
import pathfinder_bridge as pb                         # noqa: E402


@pytest.fixture(autouse=True)
def _no_insight_clone(monkeypatch):
    """Un clone InSight voisin (../InSight) ajouterait sa propre question : ces
    tests ne portent que sur Pathfinder (InSight : tests/test_config_insight.py)."""
    monkeypatch.delenv("INSIGHT_HOME", raising=False)
    monkeypatch.setattr(ib, "insight_search_paths", lambda params: [])


def _resolve(monkeypatch, argv, answers=()):
    monkeypatch.setattr(sys, "argv", ["deimos.py", *argv])
    it = iter(answers)
    asked = []

    def fake_input(prompt=""):
        asked.append(prompt)
        return next(it)
    monkeypatch.setattr(builtins, "input", fake_input)
    params = config.resolve_config(ask_params_fn=lambda: {}, ask_go_params_fn=lambda: None,
                                   go_available=False, dash_available=True,
                                   pr_path="/nonexistent/report.pr_matrix.tsv")
    return params, asked


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    def write(extra=""):
        p = tmp_path / "cfg.yaml"
        p.write_text(f"out_dir: {tmp_path / 'out'}\nrun_perseverance: false\n{extra}",
                     encoding="utf-8")
        return str(p)
    return write


def test_defaults():
    d = config._DEFAULTS
    assert d["make_pathfinder"] is False and d["pathfinder_path"] is None
    assert d["pathfinder_publications"] is False and d["pathfinder_options"] == {}
    assert d["pathfinder_timeout"] == 1800


def test_yaml_off_by_default(monkeypatch, cfg):
    params, asked = _resolve(monkeypatch, ["--config", cfg()])
    assert params["make_pathfinder"] is False and asked == []


def test_cli_flag_wins_over_yaml(monkeypatch, cfg):
    params, _ = _resolve(monkeypatch, ["--config", cfg("make_pathfinder: false\n"), "--pathfinder"])
    assert params["make_pathfinder"] is True


def test_yaml_on(monkeypatch, cfg):
    params, _ = _resolve(monkeypatch, ["--config", cfg("make_pathfinder: true\n"
                                                         "pathfinder_options: {pthr: 0.01}\n")])
    assert params["make_pathfinder"] is True
    assert params["pathfinder_options"] == {"pthr": 0.01}


def test_invalid_options_stop_before_the_run(monkeypatch, cfg):
    with pytest.raises(SystemExit):
        _resolve(monkeypatch, ["--config", cfg("pathfinder_options: {p_thr: 0.05}\n")])


def test_saved_config_keeps_the_keys(tmp_path):
    p = tmp_path / "saved.yaml"
    config.save_config({**config._DEFAULTS, "make_pathfinder": True}, str(p))
    assert config.load_config(str(p))["make_pathfinder"] is True


class TestInteractive:
    def test_question_when_a_clone_is_found(self, tmp_path, monkeypatch):
        clone = tmp_path / "Pathfinder"
        clone.mkdir()
        (clone / "pathfinder.py").write_text("__version__ = '1.1.0'\n")
        monkeypatch.setenv("PATHFINDER_HOME", str(clone))
        monkeypatch.chdir(tmp_path)
        params, asked = _resolve(monkeypatch, [], answers=["n", "y", "y"])
        assert params["make_pathfinder"] is True
        assert any("Pathfinder" in a for a in asked)

    def test_default_answer_is_no(self, tmp_path, monkeypatch):
        clone = tmp_path / "Pathfinder"
        clone.mkdir()
        (clone / "pathfinder.py").write_text("")
        monkeypatch.setenv("PATHFINDER_HOME", str(clone))
        monkeypatch.chdir(tmp_path)
        params, _ = _resolve(monkeypatch, [], answers=["n", "y", ""])
        assert params["make_pathfinder"] is False

    def test_no_question_without_a_clone(self, tmp_path, monkeypatch):
        monkeypatch.delenv("PATHFINDER_HOME", raising=False)
        monkeypatch.setattr(pb, "_HERE", str(tmp_path / "isolated" / "deimos"))
        monkeypatch.chdir(tmp_path)
        params, asked = _resolve(monkeypatch, [], answers=["n", "y"])
        assert params["make_pathfinder"] is False
        assert not any("Pathfinder" in a for a in asked)
