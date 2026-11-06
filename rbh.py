"""
rbh.py — Perseverance
-----------------
Calcule les Reciprocal Best Hits (RBH) entre un protéome requête (organisme
non-modèle) et un protéome de référence bien annoté, via DIAMOND (BLAST-
compatible, ~100-10000x plus rapide que BLAST classique, adapté à un
protéome entier plutôt qu'à quelques séquences).

Principe RBH : A (query) est le meilleur hit de B (reference) ET B est le
meilleur hit de A -> paire orthologue robuste (bien plus fiable qu'un simple
best-hit unidirectionnel, qui confond facilement orthologues et paralogues
proches — cf. le cas AT2G_KINASE1/AT4G_KINASE2 dans les données de test).

Robustesse (DIAMOND 2.2.x) : sur un FASTA nucléotidique, DIAMOND écrit
« The sequences are expected to be proteins but only contain DNA letters »
puis ne rend JAMAIS la main (makedb comme blastp). Chaque appel est donc
surveillé : FASTA contrôlés avant DIAMOND (check_protein_fasta), sortie de
DIAMOND dans un journal par étape (diamond_<étape>.log du dossier de travail),
processus tué s'il ne consomme plus de CPU pendant 3 min ou si une ligne
« Error: » n'est pas suivie de sa fin, et une ligne de vie toutes les 30 s en
mode verbeux.
"""
from __future__ import annotations
import subprocess
import tempfile
import time
import os
import pandas as pd

# Surveillance des appels DIAMOND
HEARTBEAT_S = 30.0       # ligne de vie (mode verbeux)
ERROR_GRACE_S = 5.0      # délai laissé à DIAMOND pour sortir après « Error: »
IDLE_KILL_S = 180.0      # DIAMOND sans aucun temps CPU pendant ce délai = bloqué
POLL_S = 0.5

# Contrôle des FASTA. Une séquence est nucléotidique si >= 90 % de ses lettres sont
# A/C/G/T/U/N (une protéine en contient ~25 %). DIAMOND 2.2.8 se bloque quand ses 10
# premières séquences sont toutes nucléotidiques (9 puis des protéines : accepté ;
# vérifié sur la version livrée avec Deimos) ; le fichier est aussi refusé si les
# 2 premiers millions de résidus sont nucléotidiques à >= 90 %.
NT_FRACTION_MAX = 0.90
NT_CHECK_FIRST_SEQS = 10
NT_CHECK_RESIDUES = 2_000_000
_NT = (b"A", b"C", b"G", b"T", b"U", b"N")


def clean_seqid(seqid: str) -> str:
    """Identifiant utilisable par gProfiler / l'Excel : DIAMOND rend le premier mot
    de l'en-tête FASTA tel quel, soit 'sp|Q9C5R8|NAME_ARATH' pour un FASTA UniProt,
    que gProfiler ne reconnaît pas -> accession seule ('Q9C5R8'). Les autres
    formats (AT1G01010.1, BnaA01g00010D, accession nue) sont rendus inchangés."""
    s = str(seqid)
    parts = s.split("|")
    if len(parts) >= 3 and parts[0] in ("sp", "tr") and parts[1]:
        return parts[1]
    return s


BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send",
                 "evalue", "bitscore", "qlen", "slen"]


def _fmt_duration(s: float) -> str:
    s = int(round(s))
    return f"{s // 60} min {s % 60:02d} s" if s >= 60 else f"{s} s"


def _diamond_error(text: str) -> str | None:
    """Première ligne de DIAMOND contenant « Error: », réduite à partir de ce mot."""
    for line in text.splitlines():
        i = line.find("Error:")
        if i >= 0:
            return line[i:].strip()
    return None


def _last_line(text: str) -> str:
    for line in reversed(text.splitlines()):
        if line.strip() and not line.startswith("# "):
            return line.strip()[:120]
    return ""


def _cpu_seconds(proc) -> float | None:
    """Temps CPU cumulé (user + system) du processus DIAMOND, None si indisponible.
    psutil s'il est installé, sinon /proc (Linux) ou GetProcessTimes (Windows)."""
    try:
        import psutil                                    # optionnel
        t = psutil.Process(proc.pid).cpu_times()
        return float(t.user + t.system)
    except ImportError:
        pass
    except Exception:
        return None
    try:
        if os.name == "nt":
            import ctypes
            from ctypes import wintypes
            ft = [wintypes.FILETIME() for _ in range(4)]
            ok = ctypes.windll.kernel32.GetProcessTimes(
                wintypes.HANDLE(int(proc._handle)), *[ctypes.byref(f) for f in ft])
            if not ok:
                return None
            sec = lambda f: ((f.dwHighDateTime << 32) + f.dwLowDateTime) / 1e7   # noqa: E731
            return sec(ft[2]) + sec(ft[3])               # noyau + utilisateur
        with open(f"/proc/{proc.pid}/stat") as fh:
            fields = fh.read().rsplit(")", 1)[1].split()
        return (int(fields[11]) + int(fields[12])) / os.sysconf("SC_CLK_TCK")
    except Exception:
        return None


def _error_hint(err: str | None) -> str:
    e = (err or "").lower()
    if "dna letters" in e:
        return ("FASTA nucléotidique (CDS, transcrits ou génome) : Perseverance attend des "
                "PROTÉINES (protéome UniProt, RefSeq protein .faa). DIAMOND 2.2 reste "
                "bloqué après ce message : il a été arrêté.")
    if "length 0" in e:
        return "Séquence vide dans un FASTA (en-tête sans séquence)."
    if "first line must begin" in e:
        return ("Le FASTA ne commence pas par '>' : encodage UTF-16 ou BOM UTF-8 (fichier "
                "écrit par PowerShell ou Excel ?) — le réenregistrer en UTF-8 sans BOM.")
    if "opening file" in e:
        return ("Fichier illisible ou verrouillé : un diamond.exe d'un run précédent tourne "
                "peut-être encore (gestionnaire des tâches).")
    return ""


def _run(cmd: list, description: str, verbose: bool = False,
         log_path: str | None = None, heartbeat: float | None = None) -> None:
    """Lance DIAMOND sous surveillance.

    Sortie de DIAMOND (stdout + stderr) écrite dans `log_path` (fichier temporaire
    supprimé après succès si absent), stdin fermé. Toutes les POLL_S secondes :
      - une ligne « Error: » non suivie de la fin du processus sous ERROR_GRACE_S
        -> processus tué (DIAMOND 2.2 se bloque après certaines erreurs) ;
      - aucun temps CPU consommé pendant IDLE_KILL_S -> processus tué (blocage) ;
      - en mode verbeux, une ligne de vie toutes les `heartbeat` s (HEARTBEAT_S par
        défaut) avec la durée et la dernière ligne de DIAMOND.
    Ctrl+C tue DIAMOND (plus de diamond.exe orphelin). En mode verbeux, une ligne
    avant et après chaque étape, avec la durée."""
    t0 = time.time()
    heartbeat = HEARTBEAT_S if heartbeat is None else heartbeat
    if verbose:
        print(f"    [DIAMOND] {description}...", flush=True)
    own_log = log_path is None
    if own_log:
        fd, log_path = tempfile.mkstemp(prefix="diamond_", suffix=".log")
        os.close(fd)
    with open(log_path, "wb") as fh:
        fh.write(("# " + " ".join(str(c) for c in cmd) + "\n").encode("utf-8"))
        fh.flush()
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=fh,
                                stderr=subprocess.STDOUT)
    text, pos, err, t_err, killed, last_beat = "", 0, None, None, False, t0
    idle, cpu_last, t_cpu = False, None, t0
    try:
        while True:
            done = proc.poll() is not None
            with open(log_path, "rb") as rf:
                rf.seek(pos)
                chunk = rf.read()
            pos += len(chunk)
            text += chunk.decode("utf-8", errors="replace")
            if done:
                break
            now = time.time()
            if err is None:
                err = _diamond_error(text)
                t_err = now if err else None
            elif now - t_err >= ERROR_GRACE_S:
                proc.kill()
                proc.wait()
                killed = True
                break
            # DIAMOND 2.2 bloqué (verrou entre threads) : 0 % CPU indéfiniment, et le
            # message d'erreur reste dans son tampon quand la sortie va dans un fichier.
            cpu = _cpu_seconds(proc)
            if cpu is not None:
                if cpu_last is None or cpu > cpu_last + 0.02:
                    cpu_last, t_cpu = cpu, now
                elif IDLE_KILL_S and now - t_cpu >= IDLE_KILL_S:
                    proc.kill()
                    proc.wait()
                    killed = idle = True
                    break
            if verbose and heartbeat and now - last_beat >= heartbeat:
                last_beat = now
                print(f"    [DIAMOND]   ... toujours en cours ({_fmt_duration(now - t0)})"
                      f" — {_last_line(text) or 'pas encore de sortie'}", flush=True)
            time.sleep(POLL_S)
    except KeyboardInterrupt:
        proc.kill()
        raise
    err = err or _diamond_error(text)
    if killed or proc.returncode != 0:
        tail = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("# ")]
        detail = err or "\n".join(tail[-3:])
        hint = _error_hint(err)
        msg = f"[ERREUR] {description} a échoué :\n{detail}"
        if idle:
            msg += (f"\nDIAMOND ne consommait plus de CPU depuis {_fmt_duration(IDLE_KILL_S)} "
                    f"(bloqué, dernière étape : « {_last_line(text)} ») : arrêté au bout de "
                    f"{_fmt_duration(time.time() - t0)}.")
            if not hint and "Loading sequences" in _last_line(text):
                hint = ("Blocage pendant la lecture du FASTA : vérifier qu'il contient des "
                        "protéines (DIAMOND 2.2 se bloque sur des séquences nucléotidiques).")
        elif killed:
            msg += f"\nDIAMOND ne s'est pas arrêté après l'erreur : arrêté au bout de {_fmt_duration(time.time() - t0)}."
        if hint:
            msg += f"\n{hint}"
        msg += f"\n(journal DIAMOND : {log_path})"
        raise RuntimeError(msg)
    if own_log:
        try:
            os.remove(log_path)
        except OSError:
            pass
    if verbose:
        print(f"    [DIAMOND]   ok ({time.time() - t0:.0f} s)", flush=True)


def check_protein_fasta(fasta_path: str) -> list[str]:
    """Contrôle d'un FASTA avant DIAMOND. Renvoie les problèmes bloquants (liste vide
    = utilisable) : encodage UTF-16 ou BOM UTF-8, fins de ligne CR seules (ancien
    Mac), première ligne sans '>', aucune séquence, séquences vides, contenu
    nucléotidique (DIAMOND 2.2 se bloque dessus). Lecture binaire en flux."""
    name = os.path.basename(fasta_path)
    try:
        fh = open(fasta_path, "rb")
    except OSError as e:
        return [f"{name} : illisible ({e})."]
    with fh:
        head = fh.read(1 << 20)
        if head.startswith(b"\xff\xfe") or head.startswith(b"\xfe\xff"):
            return [f"{name} : encodé en UTF-16 (fichier écrit par PowerShell 5 avec '>' "
                    f"ou Out-File ?) — le réenregistrer en UTF-8."]
        if head.startswith(b"\xef\xbb\xbf"):
            return [f"{name} : commence par un BOM UTF-8, que DIAMOND refuse — le "
                    f"réenregistrer en UTF-8 sans BOM."]
        if b"\r" in head and b"\n" not in head:
            return [f"{name} : fins de ligne CR seules (format Mac ancien) — convertir "
                    f"en LF ou CRLF."]
        fh.seek(0)
        n_seq = n_empty = 0
        cur_len, cur_id, empty_ids = None, "", []
        nt_all = let_all = nt_cur = let_cur = 0
        head_nt = []                     # 10 premières séquences : nucléotidique ?
        first = True

        def _close_seq():
            if n_seq and len(head_nt) < NT_CHECK_FIRST_SEQS and n_seq <= NT_CHECK_FIRST_SEQS:
                head_nt.append(let_cur > 0 and nt_cur / let_cur >= NT_FRACTION_MAX)

        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if first:
                first = False
                if not line.startswith(b">"):
                    return [f"{name} : la première ligne ne commence pas par '>' "
                            f"(pas un FASTA ?)."]
            if line.startswith(b">"):
                _close_seq()
                if cur_len == 0:
                    n_empty += 1
                    if len(empty_ids) < 3:
                        empty_ids.append(cur_id)
                n_seq += 1
                cur_len, nt_cur, let_cur = 0, 0, 0
                cur_id = line[1:].split(maxsplit=1)[0].decode("utf-8", "replace") \
                    if len(line) > 1 else "?"
                continue
            cur_len += len(line)
            if let_all < NT_CHECK_RESIDUES or n_seq <= NT_CHECK_FIRST_SEQS:
                up = line.upper()
                k = sum(up.count(c) for c in _NT)
                nt_cur += k
                let_cur += len(up)
                if let_all < NT_CHECK_RESIDUES:
                    nt_all += k
                    let_all += len(up)
        _close_seq()
        if cur_len == 0:
            n_empty += 1
            empty_ids.append(cur_id)
    problems = []
    if n_seq == 0:
        return [f"{name} : aucune séquence."]
    frac = nt_all / let_all if let_all else 0.0
    head_dna = bool(head_nt) and all(head_nt)
    if head_dna and frac < NT_FRACTION_MAX:
        problems.append(
            f"{name} : les {len(head_nt)} premières séquences sont nucléotidiques — "
            f"DIAMOND 2.2 se bloquerait sur ce fichier ; Perseverance attend des PROTÉINES.")
    elif frac >= NT_FRACTION_MAX:
        problems.append(
            f"{name} : séquences nucléotidiques ({100 * frac:.0f} % de A/C/G/T/U/N) — "
            f"Perseverance attend des PROTÉINES (protéome UniProt, RefSeq protein .faa), "
            f"pas des CDS ni des transcrits. DIAMOND 2.2 se bloquerait sur ce fichier.")
    if n_empty:
        problems.append(f"{name} : {n_empty} séquence(s) vide(s) (ex. {', '.join(empty_ids[:3])})"
                        f" — DIAMOND s'arrête dessus.")
    return problems


def _n_seqs(fasta_path: str) -> int:
    try:
        with open(fasta_path, encoding="utf-8", errors="replace") as fh:
            return sum(1 for line in fh if line.startswith(">"))
    except OSError:
        return 0


def build_diamond_db(fasta_path: str, db_path: str, verbose: bool = False,
                     log_path: str | None = None) -> None:
    _run(["diamond", "makedb", "--in", fasta_path, "-d", db_path],
         f"Construction de la base DIAMOND ({os.path.basename(fasta_path)})", verbose,
         log_path=log_path)


def run_blastp(query_fasta: str, db_path: str, out_tsv: str,
               evalue: float = 1e-5, threads: int = 4,
               verbose: bool = False, log_path: str | None = None) -> pd.DataFrame:
    # sans --quiet : la progression de DIAMOND (blocs traités) va dans le journal et
    # alimente la ligne de vie ; la console n'affiche que les lignes [DIAMOND].
    _run(["diamond", "blastp",
          "--query", query_fasta, "--db", db_path,
          "--out", out_tsv, "--outfmt", "6"] + BLAST_COLUMNS +
         ["--evalue", str(evalue), "--max-target-seqs", "1",
          "--threads", str(threads)],
         f"blastp {os.path.basename(query_fasta)} ({_n_seqs(query_fasta)} séquences) -> "
         f"{os.path.basename(db_path)}, {threads} threads", verbose, log_path=log_path)
    if os.path.getsize(out_tsv) == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)
    return pd.read_csv(out_tsv, sep="\t", names=BLAST_COLUMNS)


def _best_hit_per_query(df: pd.DataFrame) -> pd.DataFrame:
    """--max-target-seqs 1 garantit déjà 1 ligne par query, mais on trie par
    bitscore décroissant par sécurité si jamais plusieurs lignes existent."""
    if df.empty:
        return df
    return df.sort_values("bitscore", ascending=False).drop_duplicates("qseqid", keep="first")


def compute_rbh(query_fasta: str, reference_fasta: str, workdir: str,
                evalue: float = 1e-5, min_identity: float = 25.0,
                min_coverage: float = 50.0, threads: int = 4,
                verbose: bool = False) -> pd.DataFrame:
    """
    Pipeline RBH complet.

    Parameters
    ----------
    min_identity  : % identité minimum pour retenir un hit (avant même de
                    tester la réciprocité) — filtre les homologies trop
                    lointaines pour être fonctionnellement fiables.
    min_coverage  : % de couverture minimum (sur la plus courte des deux
                    séquences) — évite les hits locaux sur un seul domaine
                    partagé entre protéines par ailleurs non-orthologues.

    Returns
    -------
    DataFrame : query_id, reference_id, pident, coverage, evalue, bitscore,
                confidence (high/medium/low)
    """
    os.makedirs(workdir, exist_ok=True)
    db_query = os.path.join(workdir, "query_db")
    db_ref = os.path.join(workdir, "reference_db")
    if verbose:
        print(f"    [DIAMOND] RBH : {_n_seqs(query_fasta)} séquences de la recherche <-> "
              f"{_n_seqs(reference_fasta)} séquences de référence (4 étapes ; fichiers et "
              f"journaux diamond_*.log dans {workdir})", flush=True)
    problems = check_protein_fasta(query_fasta) + check_protein_fasta(reference_fasta)
    if problems:
        raise RuntimeError("[ERREUR] FASTA inutilisable(s) pour DIAMOND :\n  - "
                           + "\n  - ".join(problems))
    log = lambda step: os.path.join(workdir, f"diamond_{step}.log")      # noqa: E731
    build_diamond_db(query_fasta, db_query, verbose, log_path=log("makedb_query"))
    build_diamond_db(reference_fasta, db_ref, verbose, log_path=log("makedb_reference"))

    fwd = run_blastp(query_fasta, db_ref,
                     os.path.join(workdir, "fwd.tsv"), evalue, threads, verbose,
                     log_path=log("blastp_fwd"))
    rev = run_blastp(reference_fasta, db_query,
                     os.path.join(workdir, "rev.tsv"), evalue, threads, verbose,
                     log_path=log("blastp_rev"))

    fwd_best = _best_hit_per_query(fwd)
    rev_best = _best_hit_per_query(rev)

    if fwd_best.empty or rev_best.empty:
        return pd.DataFrame(columns=["query_id", "reference_id", "pident",
                                     "coverage", "evalue", "bitscore", "confidence"])

    # RBH : qseqid(fwd) -> sseqid(fwd) doit correspondre exactement a
    # qseqid(rev) -> sseqid(rev) en sens inverse
    rev_map = dict(zip(rev_best["qseqid"], rev_best["sseqid"]))
    rows = []
    for _, r in fwd_best.iterrows():
        q, ref = r["qseqid"], r["sseqid"]
        if rev_map.get(ref) == q:
            # longueur d'alignement (gaps compris) / plus courte séquence, plafonnée à 100
            coverage = min(100.0, 100.0 * r["length"] / min(r["qlen"], r["slen"]))
            if r["pident"] >= min_identity and coverage >= min_coverage:
                rows.append({
                    "query_id": q, "reference_id": ref,
                    "pident": round(r["pident"], 1),
                    "coverage": round(coverage, 1),
                    "evalue": r["evalue"], "bitscore": r["bitscore"],
                })

    result = pd.DataFrame(rows)
    if result.empty:
        return result

    def _confidence(row):
        if row["pident"] >= 50 and row["coverage"] >= 70:
            return "high"
        if row["pident"] >= 30 and row["coverage"] >= 50:
            return "medium"
        return "low"

    result["confidence"] = result.apply(_confidence, axis=1)
    # En-têtes UniProt (sp|ACC|NAME) -> accession, des deux côtés
    result["query_id"] = result["query_id"].map(clean_seqid)
    result["reference_id"] = result["reference_id"].map(clean_seqid)
    return result.sort_values("bitscore", ascending=False).reset_index(drop=True)
