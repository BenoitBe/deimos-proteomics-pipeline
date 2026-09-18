"""
deimos_to_beagle.py — pont d'intégration Beagle <-> Deimos
=============================================================
Deimos exporte des IDENTIFIANTS de protéines (Protein.Group, souvent des
groupes composites "A0A076VEN9;A0A078I920") mais AUCUNE séquence. Beagle a
besoin de séquences (FASTA) pour calculer l'homologie par RBH. Ce module
fait le pont dans les deux sens :

  1. export_query_fasta() : à partir du FASTA1 (le FASTA de recherche
     DIA-NN/PEAKS déjà utilisé par Deimos) et de la liste des Protein.Group
     réellement quantifiés par Deimos, extrait UNIQUEMENT ces séquences —
     pas besoin de faire tourner Beagle sur un protéome entier alors que
     Deimos n'en a identifié qu'une fraction.

  2. remap_annotation_to_protein_groups() : la sortie de Beagle est indexée
     par accession individuelle (un header FASTA = une ligne). Pour un
     groupe composite Deimos "A;B", si Beagle a annoté A et/ou B, cette
     fonction fusionne (union des GO, confiance = la meilleure des deux)
     pour reconstituer une ligne par Protein.Group Deimos — le format que
     go_enrichment.py attend en aval.
"""
from __future__ import annotations
import re
import pandas as pd


def _parse_fasta_headers(fasta_path: str) -> dict:
    """
    Retourne {accession: header_complet} pour chaque séquence du FASTA.

    Gère deux conventions de header courantes :
      - UniProt pipe-delimited : >sp|P12345|NAME_SPECIES description...
                                  >tr|A0A076U3R7|A0A076U3R7_BRANA ...
        -> l'accession est le 2e champ entre pipes.
      - Header simple : >A0A076U3R7 description...
        -> l'accession est le 1er token.
    """
    headers = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                if header.startswith(("sp|", "tr|")) or "|" in header.split()[0]:
                    parts = header.split("|")
                    if len(parts) >= 2:
                        headers[parts[1]] = header
                        continue
                accession = header.split()[0]
                headers[accession] = header
    return headers


def _extract_sequences(fasta_path: str, wanted: set) -> dict:
    """Extrait les séquences (accession -> séquence) pour les accessions demandées."""
    seqs = {}
    current_acc = None
    current_seq = []

    def _flush():
        if current_acc is not None and current_acc in wanted:
            seqs[current_acc] = "".join(current_seq)

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                _flush()
                header = line[1:].strip()
                if header.startswith(("sp|", "tr|")):
                    parts = header.split("|")
                    current_acc = parts[1] if len(parts) >= 2 else header.split()[0]
                else:
                    current_acc = header.split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        _flush()
    return seqs


def export_query_fasta(search_fasta: str, protein_groups: list,
                       out_fasta: str, out_mapping_tsv: str) -> dict:
    """
    Parameters
    ----------
    search_fasta    : FASTA1 — le FASTA de recherche déjà utilisé par Deimos
                      (celui donné à DIA-NN/PEAKS)
    protein_groups  : liste des valeurs Protein.Group de Deimos (peut
                      contenir des groupes composites "A;B")
    out_fasta       : FASTA de sortie, une séquence par accession
                      INDIVIDUELLE — c'est ce fichier qui sert de --query à
                      beagle.py
    out_mapping_tsv : table accession -> protein_group (groupe composite
                      d'origine), nécessaire pour remap_annotation_to_protein_groups()

    Returns
    -------
    dict : n_groups, n_accessions_wanted, n_found, n_missing, missing (liste)
    """
    all_accessions = set()
    mapping_rows = []
    for pg in protein_groups:
        accs = [a.strip() for a in str(pg).split(";")]
        for a in accs:
            all_accessions.add(a)
            mapping_rows.append({"accession": a, "protein_group": pg})

    sequences = _extract_sequences(search_fasta, all_accessions)
    found = set(sequences.keys())
    missing = sorted(all_accessions - found)

    with open(out_fasta, "w") as f:
        for acc, seq in sequences.items():
            f.write(f">{acc}\n{seq}\n")

    pd.DataFrame(mapping_rows).to_csv(out_mapping_tsv, sep="\t", index=False)

    report = {
        "n_groups": len(protein_groups),
        "n_accessions_wanted": len(all_accessions),
        "n_found": len(found),
        "n_missing": len(missing),
        "missing": missing,
    }
    print(f"  [OK] {report['n_found']}/{report['n_accessions_wanted']} accessions "
          f"trouvées dans {search_fasta} -> {out_fasta}")
    if missing:
        print(f"  [WARN] {len(missing)} accession(s) absente(s) du FASTA1 "
              "(ex: entrées obsolètes, mismatch de version de base) : "
              f"{missing[:5]}{'...' if len(missing) > 5 else ''}")
    return report


def remap_annotation_to_protein_groups(beagle_annotation: pd.DataFrame,
                                       mapping_tsv: str) -> pd.DataFrame:
    """
    Remonte l'annotation Beagle (indexée par accession individuelle) vers
    la convention Protein.Group de Deimos (groupes composites "A;B").

    Règle de fusion pour un groupe composite dont plusieurs accessions sont
    annotées : union des GO/KEGG/Pfam, confiance = la meilleure des deux
    (high > medium > low), source='orthologue (RBH)' si au moins un membre
    provient d'un transfert (le detail par accession reste consultable dans
    le fichier Beagle brut si besoin d'audit).
    """
    mapping = pd.read_csv(mapping_tsv, sep="\t")
    merged = mapping.merge(beagle_annotation, left_on="accession",
                           right_on="protein_id", how="inner")
    if merged.empty:
        return pd.DataFrame(columns=["protein_id", "GO", "KEGG", "Pfam",
                                     "source", "confidence"])

    conf_rank = {"high": 3, "medium": 2, "low": 1}
    merged["_conf_rank"] = merged["confidence"].map(conf_rank).fillna(0)

    def _union_semicolon(series):
        vals = set()
        for v in series.dropna():
            vals.update(x.strip() for x in str(v).split(";") if x.strip())
        return ";".join(sorted(vals))

    grouped = merged.groupby("protein_group").agg(
        GO=("GO", _union_semicolon),
        KEGG=("KEGG", _union_semicolon),
        Pfam=("Pfam", _union_semicolon),
        source=("source", "first"),
        _best_conf_rank=("_conf_rank", "max"),
    ).reset_index()

    rank_to_conf = {v: k for k, v in conf_rank.items()}
    grouped["confidence"] = grouped["_best_conf_rank"].map(rank_to_conf)
    grouped = grouped.drop(columns=["_best_conf_rank"])

    # IMPORTANT : go_enrichment.py (isolate_significant, process_enrichment/
    # _zscore) nettoie systematiquement les groupes composites en ne gardant
    # que le PREMIER accession ("A;B" -> "A") pour tout matching. On aligne
    # la cle de sortie sur cette meme convention plutot que de garder la
    # chaine composite complete, sous peine de mismatch silencieux en aval.
    grouped["protein_id"] = grouped["protein_group"].str.split(";").str[0]
    grouped = grouped.drop(columns=["protein_group"])

    return grouped[["protein_id", "GO", "KEGG", "Pfam", "source", "confidence"]]
