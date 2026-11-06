"""DIAMOND 2.2 se bloque après « sequences ... only contain DNA letters » : FASTA
contrôlés avant DIAMOND, processus tué si une erreur n'est pas suivie de sa fin,
ligne de vie, journal par étape. Test d'intégration si `diamond` est installé."""
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import rbh                                             # noqa: E402

PROT = "MKVLAAGIVGLLLAQSHAEDKTPWYRRQCNFMEHGT" * 4
DNA = "ATGCGTACGTAGCTAGCTAGCTAGCATCGATCGGCTANNACGT" * 4


def _fasta(path, seqs, prefix="s"):
    path.write_text("".join(f">{prefix}{i} desc\n{s}\n" for i, s in enumerate(seqs)))
    return str(path)


def test_check_protein_fasta(tmp_path):
    assert rbh.check_protein_fasta(_fasta(tmp_path / "p.fa", [PROT] * 5)) == []
    nt = rbh.check_protein_fasta(_fasta(tmp_path / "n.fa", [DNA] * 5))
    assert len(nt) == 1 and "nucléotidiques" in nt[0] and "PROTÉINES" in nt[0]
    # DIAMOND 2.2.8 se bloque si ses 10 premières séquences sont nucléotidiques
    mixed = tmp_path / "m.fa"
    mixed.write_text(_fasta_text([DNA] * 10) + _fasta_text([PROT] * 400, "p"))
    msg = rbh.check_protein_fasta(str(mixed))
    assert len(msg) == 1 and "10 premières séquences" in msg[0]
    nine = tmp_path / "m9.fa"                          # 9 puis protéines : accepté
    nine.write_text(_fasta_text([DNA] * 9) + _fasta_text([PROT] * 400, "p"))
    assert rbh.check_protein_fasta(str(nine)) == []
    empty = tmp_path / "e.fa"
    empty.write_text(f">a\n\n>b\n{PROT}\n>c\n")
    msg = rbh.check_protein_fasta(str(empty))
    assert len(msg) == 1 and "2 séquence(s) vide(s)" in msg[0] and "a, c" in msg[0]
    bom = tmp_path / "bom.fa"
    bom.write_bytes(b"\xef\xbb\xbf>a\n" + PROT.encode() + b"\n")
    assert "BOM" in rbh.check_protein_fasta(str(bom))[0]
    u16 = tmp_path / "u16.fa"
    u16.write_bytes(f">a\n{PROT}\n".encode("utf-16"))
    assert "UTF-16" in rbh.check_protein_fasta(str(u16))[0]
    cr = tmp_path / "cr.fa"
    cr.write_bytes(f">a\r{PROT}\r>b\r{PROT}\r".encode())
    assert "CR" in rbh.check_protein_fasta(str(cr))[0]
    crlf = tmp_path / "crlf.fa"
    crlf.write_bytes(f">a\r\n{PROT}\r\n>b\r\n{PROT}\r\n".encode())
    assert rbh.check_protein_fasta(str(crlf)) == []
    notfa = tmp_path / "x.fa"
    notfa.write_text("protein_id\tGO\n")
    assert "'>'" in rbh.check_protein_fasta(str(notfa))[0]


def test_truncated_uniprot_download_and_bad_characters(tmp_path):
    """Cas réel (référence souris, 06/11) : flux UniProt coupé à 40 000 séquences, le
    fichier se termine par le message d'erreur du serveur -> DIAMOND 2.2.8 bloqué."""
    cut = tmp_path / "ref.fa"
    cut.write_text(_fasta_text([PROT] * 12) + "MKVLA\n\n"
                   "Error encountered when streaming data. Please try again later.\n")
    msg = rbh.check_protein_fasta(str(cut))
    assert len(msg) == 1 and "UniProt interrompu" in msg[0] and "12 séquences" in msg[0]
    dot = tmp_path / "dot.fa"
    dot.write_text(_fasta_text([PROT] * 3) + ">x\nMKV.LA\n>y\nMK VL\n")
    msg = rbh.check_protein_fasta(str(dot))
    assert len(msg) == 1 and "2 ligne(s)" in msg[0] and "'.'" in msg[0] and "(x)" in msg[0]
    ok = tmp_path / "ok.fa"                         # '*' et '-' sont admis par DIAMOND
    ok.write_text(_fasta_text([PROT + "*", "MKV-LA" * 20]))
    assert rbh.check_protein_fasta(str(ok)) == []


def _fasta_text(seqs, prefix="s"):
    return "".join(f">{prefix}{i}\n{s}\n" for i, s in enumerate(seqs))


def test_nucleotide_reference_stops_before_diamond(tmp_path, monkeypatch):
    called = []
    monkeypatch.setattr(rbh.subprocess, "Popen", lambda *a, **k: called.append(a))
    q = _fasta(tmp_path / "q.fa", [PROT] * 3)
    r = _fasta(tmp_path / "r.fa", [DNA] * 3)
    with pytest.raises(RuntimeError, match="nucléotidiques"):
        rbh.compute_rbh(q, r, str(tmp_path / "w"), verbose=True)
    assert called == []


def test_error_then_hang_is_killed(tmp_path, monkeypatch):
    """Simule DIAMOND 2.2 : message d'erreur puis processus qui ne se termine pas."""
    monkeypatch.setattr(rbh, "ERROR_GRACE_S", 0.5)
    monkeypatch.setattr(rbh, "POLL_S", 0.1)
    code = ("import time; print('Loading sequences... Error: The sequences are expected "
            "to be proteins but only contain DNA letters.', flush=True); time.sleep(120)")
    log = tmp_path / "diamond_makedb_reference.log"
    t0 = time.time()
    with pytest.raises(RuntimeError) as e:
        rbh._run([sys.executable, "-c", code], "Construction de la base DIAMOND (r.fa)",
                 log_path=str(log))
    assert time.time() - t0 < 30
    msg = str(e.value)
    assert "DNA letters" in msg and "PROTÉINES" in msg and "arrêté" in msg
    assert str(log) in msg and log.read_text().startswith("# ")


def test_idle_process_is_killed(tmp_path, monkeypatch):
    """Blocage sans message (cas réel : DIAMOND 2.2 garde l'erreur dans son tampon)."""
    if rbh._cpu_seconds(subprocess.Popen([sys.executable, "-c", "pass"])) is None \
            and sys.platform == "darwin":
        pytest.skip("temps CPU indisponible")
    monkeypatch.setattr(rbh, "IDLE_KILL_S", 1.0)
    monkeypatch.setattr(rbh, "POLL_S", 0.1)
    code = "import time; print('Loading sequences... ', end='', flush=True); time.sleep(120)"
    t0 = time.time()
    with pytest.raises(RuntimeError) as e:
        rbh._run([sys.executable, "-c", code], "makedb r.fa",
                 log_path=str(tmp_path / "d.log"))
    assert time.time() - t0 < 30
    assert "ne consommait plus de CPU" in str(e.value) and "nucléotidiques" in str(e.value)
    # un processus qui calcule n'est pas tué
    busy = "import time\nt = time.time()\nwhile time.time() - t < 2.5: pass"
    rbh._run([sys.executable, "-c", busy], "blastp", log_path=str(tmp_path / "b.log"))


def test_nonzero_exit_and_heartbeat(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(rbh, "POLL_S", 0.1)
    code = "import sys, time; print('Processing query block 1', flush=True); time.sleep(1.2)"
    rbh._run([sys.executable, "-c", code], "blastp q -> r", verbose=True, heartbeat=0.3)
    out = capsys.readouterr().out
    assert "toujours en cours" in out and "Processing query block 1" in out
    assert "[DIAMOND]   ok" in out
    code = "import sys; print('Error: Error opening file x.dmnd', flush=True); sys.exit(1)"
    with pytest.raises(RuntimeError, match="verrouillé"):
        rbh._run([sys.executable, "-c", code], "makedb")


@pytest.mark.skipif(shutil.which("diamond") is None, reason="DIAMOND absent")
def test_real_diamond_rbh(tmp_path):
    import random
    rnd = random.Random(1)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    seqs = ["".join(rnd.choice(aa) for _ in range(300)) for _ in range(10)]
    mut = ["".join(c if rnd.random() > 0.1 else rnd.choice(aa) for c in s) for s in seqs]
    q = _fasta(tmp_path / "q.fa", seqs, "q")
    r = _fasta(tmp_path / "r.fa", mut, "r")
    res = rbh.compute_rbh(q, r, str(tmp_path / "w"), min_identity=25, threads=1)
    assert len(res) == 10 and (res["query_id"].str[1:] == res["reference_id"].str[1:]).all()
    assert (tmp_path / "w" / "diamond_blastp_rev.log").exists()
    out = subprocess.run(["diamond", "version"], capture_output=True, text=True).stdout
    assert "diamond" in out.lower()
