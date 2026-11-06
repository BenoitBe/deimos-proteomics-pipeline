"""DIAMOND tourne en silence : Deimos affiche chaque étape du RBH et sa durée, et
utilise par défaut tous les cœurs moins un (perseverance_threads pour forcer)."""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import rbh                                             # noqa: E402


class _Ok:
    returncode, stderr = 0, ""


def test_progress_lines_and_threads(tmp_path, monkeypatch, capsys):
    calls = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        if "--out" in cmd:
            open(cmd[cmd.index("--out") + 1], "w").close()
        return _Ok()

    class FakePopen:                    # DIAMOND est lancé sous surveillance (Popen)
        def __init__(self, cmd, **kw):
            fake_run(cmd)
            self.returncode = 0

        def poll(self):
            return 0

        def wait(self, timeout=None):
            return 0

        def kill(self):
            pass
    monkeypatch.setattr(rbh.subprocess, "run", fake_run)
    monkeypatch.setattr(rbh.subprocess, "Popen", FakePopen)
    q, r = tmp_path / "q.fa", tmp_path / "r.fa"
    q.write_text(">a\nMA\n>b\nMA\n")
    r.write_text(">c\nMA\n")
    rbh.compute_rbh(str(q), str(r), str(tmp_path / "w"), threads=7, verbose=True)
    out = capsys.readouterr().out
    assert "2 séquences de la recherche <-> 1 séquences de référence" in out
    assert out.count("[DIAMOND]   ok") == 4
    assert [c[1] for c in calls] == ["makedb", "makedb", "blastp", "blastp"]
    assert calls[2][calls[2].index("--threads") + 1] == "7"
    rbh.compute_rbh(str(q), str(r), str(tmp_path / "w2"))          # défaut : silencieux
    assert "[DIAMOND]" not in capsys.readouterr().out


def test_thread_setting():
    import deimos
    assert deimos._diamond_threads({"perseverance_threads": 3}) == 3
    assert deimos._diamond_threads({}) == max(1, (os.cpu_count() or 4) - 1)
    assert deimos._diamond_threads({"perseverance_threads": "x"}) >= 1
