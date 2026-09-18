"""
deimos_to_perseverance.py — pont d'intégration Perseverance <-> Deimos
=============================================================
Deimos exporte des IDENTIFIANTS de protéines (Protein.Group, souvent des
groupes composites "A0A076VEN9;A0A078I920") mais AUCUNE séquence. Perseverance a
besoin de séquences (FASTA) pour calculer l'homologie par RBH. Ce module
fait le pont dans les deux sens :

  1. export_query_fasta() : à partir du FASTA1 (le FASTA de recherche
     DIA-NN/PEAKS déjà utilisé par Deimos) et de la liste des Protein.Group
     réellement quantifiés par Deimos, extrait UNIQUEMENT ces séquences —
     pas besoin de faire tourner Perseverance sur un protéome entier alors que
     Deimos n'en a identifié qu'une fraction.

  2. remap_annotation_to_protein_groups() : la sortie de Perseverance est indexée
     par accession individuelle (un header FASTA = une ligne). Pour un
     groupe composite Deimos "A;B", si Perseverance a annoté A et/ou B, cette
     fonction fusionne (union des GO, confiance = la meilleure des deux)
     pour reconstituer une ligne par Protein.Group Deimos — le format que
     go_enrichment.py attend en aval.

  3. detect_fasta_roles() / build_ortholog_map() : pour le flux "RBH pur"
     (Perseverance sert uniquement à trouver l'orthologue, pas à transférer
     de GO — l'enrichissement final passe par le VRAI gProfiler interrogé
     sur l'espèce de référence). detect_fasta_roles() distingue automatiquement,
     parmi deux FASTA fournis, lequel correspond au FASTA de recherche DIA-NN
     de CE run (fort recouvrement avec Protein.Group) et lequel est le
     protéome de référence pour l'orthologie. build_ortholog_map() convertit
     la sortie RBH en dict {protein_id (Deimos) -> ortholog_id (référence)},
     directement utilisable pour traduire une liste de protéines
     significatives avant de la soumettre à gProfiler.
"""
from __future__ import annotations
import os
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
                      perseverance.py
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


def remap_annotation_to_protein_groups(perseverance_annotation: pd.DataFrame,
                                       mapping_tsv: str) -> pd.DataFrame:
    """
    Remonte l'annotation Perseverance (indexée par accession individuelle) vers
    la convention Protein.Group de Deimos (groupes composites "A;B").

    Règle de fusion pour un groupe composite dont plusieurs accessions sont
    annotées : union des GO/KEGG/Pfam, confiance = la meilleure des deux
    (high > medium > low), source='orthologue (RBH)' si au moins un membre
    provient d'un transfert (le detail par accession reste consultable dans
    le fichier Perseverance brut si besoin d'audit).
    """
    mapping = pd.read_csv(mapping_tsv, sep="\t")
    merged = mapping.merge(perseverance_annotation, left_on="accession",
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


def _fasta_accessions(fasta_path: str) -> set:
    """Liste des accessions (mêmes règles de parsing que _parse_fasta_headers)."""
    accs = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                if header.startswith(("sp|", "tr|")):
                    parts = header.split("|")
                    accs.add(parts[1] if len(parts) >= 2 else header.split()[0])
                else:
                    accs.add(header.split()[0])
    return accs


def detect_fasta_roles(fasta_a: str, fasta_b: str, protein_groups: list,
                       ambiguous_margin: float = 0.15) -> dict:
    """
    Détecte automatiquement, parmi deux FASTA, lequel est le FASTA de
    RECHERCHE (celui utilisé par DIA-NN/PEAKS pour CE run Deimos — fort
    recouvrement attendu avec Protein.Group) et lequel est le FASTA de
    RÉFÉRENCE pour l'orthologie (l'autre, ex: Arabidopsis).

    Principe : calcule le % d'accessions de chaque FASTA retrouvées dans
    les Protein.Group du run en cours. Le FASTA avec le recouvrement le
    plus élevé est désigné "search" (= FASTA1, --query de perseverance.py),
    l'autre "reference" (= FASTA2, --reference de perseverance.py).

    Parameters
    ----------
    ambiguous_margin : écart minimum entre les deux taux de recouvrement
                       pour trancher automatiquement. En-dessous, le
                       résultat est marqué ambigu (à l'appelant de
                       demander confirmation plutôt que de deviner).

    Returns
    -------
    dict : search_fasta, reference_fasta, overlap_a, overlap_b, ambiguous (bool)
    """
    wanted = set()
    for pg in protein_groups:
        wanted.update(a.strip() for a in str(pg).split(";"))

    acc_a = _fasta_accessions(fasta_a)
    acc_b = _fasta_accessions(fasta_b)

    overlap_a = len(acc_a & wanted) / len(wanted) if wanted else 0.0
    overlap_b = len(acc_b & wanted) / len(wanted) if wanted else 0.0

    # Ambigu seulement si les deux ont un recouvrement non-trivial et proche —
    # si l'un des deux est a 0 (ou quasi), il n'y a pas d'ambiguite possible
    # meme si l'ecart absolu est petit en valeur.
    if min(overlap_a, overlap_b) < 0.01:
        ambiguous = False
    else:
        ambiguous = abs(overlap_a - overlap_b) < ambiguous_margin

    if overlap_a >= overlap_b:
        search_fasta, reference_fasta = fasta_a, fasta_b
    else:
        search_fasta, reference_fasta = fasta_b, fasta_a

    print(f"  [AUTO] Recouvrement avec les Protein.Group du run : "
          f"{os.path.basename(fasta_a)}={overlap_a:.1%}, "
          f"{os.path.basename(fasta_b)}={overlap_b:.1%}")
    if ambiguous:
        print(f"  [WARN] Écart de recouvrement trop faible (<{ambiguous_margin:.0%}) "
              f"pour trancher automatiquement avec confiance.")

    return {
        "search_fasta": search_fasta, "reference_fasta": reference_fasta,
        "overlap_a": overlap_a, "overlap_b": overlap_b, "ambiguous": ambiguous,
    }


def build_ortholog_map(rbh_df: pd.DataFrame, mapping_tsv: str) -> dict:
    """
    Convertit la sortie RBH (indexée par accession individuelle) en dict
    {protein_id (convention Deimos, 1er accession d'un groupe) -> ortholog_id
    (accession de l'espèce de référence)} — pour traduire une liste de
    protéines significatives avant une requête gProfiler.

    Si plusieurs accessions d'un même groupe composite ont chacune un hit
    RBH (rare mais possible), on garde celui au meilleur bitscore — un seul
    orthologue par groupe, pas une liste, pour rester compatible avec
    l'usage simple "traduire un ID -> un ID".
    """
    mapping = pd.read_csv(mapping_tsv, sep="\t")
    merged = mapping.merge(rbh_df, left_on="accession", right_on="query_id", how="inner")
    if merged.empty:
        return {}

    merged["protein_id"] = merged["protein_group"].str.split(";").str[0]
    best = (merged.sort_values("bitscore", ascending=False)
           .drop_duplicates("protein_id", keep="first"))

    return dict(zip(best["protein_id"], best["reference_id"]))
