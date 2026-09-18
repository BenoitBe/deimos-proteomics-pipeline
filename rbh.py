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
"""
from __future__ import annotations
import subprocess
import os
import pandas as pd


BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send",
                 "evalue", "bitscore", "qlen", "slen"]


def _run(cmd: list, description: str) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"[ERREUR] {description} a échoué :\n{result.stderr}")


def build_diamond_db(fasta_path: str, db_path: str) -> None:
    _run(["diamond", "makedb", "--in", fasta_path, "-d", db_path],
         f"Construction de la base DIAMOND ({fasta_path})")


def run_blastp(query_fasta: str, db_path: str, out_tsv: str,
               evalue: float = 1e-5, threads: int = 4) -> pd.DataFrame:
    _run(["diamond", "blastp",
          "--query", query_fasta, "--db", db_path,
          "--out", out_tsv, "--outfmt", "6"] + BLAST_COLUMNS +
         ["--evalue", str(evalue), "--max-target-seqs", "1",
          "--threads", str(threads), "--quiet"],
         f"Recherche DIAMOND blastp ({query_fasta} vs {db_path})")
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
                min_coverage: float = 50.0, threads: int = 4) -> pd.DataFrame:
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
    build_diamond_db(query_fasta, db_query)
    build_diamond_db(reference_fasta, db_ref)

    fwd = run_blastp(query_fasta, db_ref,
                     os.path.join(workdir, "fwd.tsv"), evalue, threads)
    rev = run_blastp(reference_fasta, db_query,
                     os.path.join(workdir, "rev.tsv"), evalue, threads)

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
            coverage = 100.0 * r["length"] / min(r["qlen"], r["slen"])
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
    return result.sort_values("bitscore", ascending=False).reset_index(drop=True)
