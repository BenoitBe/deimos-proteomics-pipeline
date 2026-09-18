"""
transfer.py — Perseverance
----------------------
Transfère les annotations fonctionnelles (GO, KEGG, Pfam) des orthologues
de référence identifiés par RBH (rbh.py) vers les protéines requête.

Principe de traçabilité : CHAQUE annotation transférée conserve son niveau
de confiance (high/medium/low, dérivé de %identité et %couverture RBH) —
jamais mélangée silencieusement avec une annotation directe. C'est la
différence essentielle avec une simple jointure de table : un consommateur
en aval (go_enrichment.py, Pathfinder Meta-GO) doit pouvoir choisir de ne
garder que les annotations 'high', ou de pondérer l'enrichissement par
confiance, plutôt que de traiter un transfert à 25% d'identité comme
équivalent à une annotation directe UniProt.
"""
from __future__ import annotations
import pandas as pd


def transfer_annotations(rbh_df: pd.DataFrame, ref_annot_path: str,
                         direct_annot_path: "str | None" = None
                         ) -> pd.DataFrame:
    """
    Construit la table finale d'annotation, une ligne par protéine query
    annotée (directement ou par transfert).

    Parameters
    ----------
    rbh_df             : sortie de rbh.compute_rbh()
    ref_annot_path     : TSV référence (protein_id, GO, KEGG, Pfam) — GO en
                         liste ';'-séparée
    direct_annot_path  : TSV optionnel d'annotations DIRECTES déjà connues
                         pour l'espèce query (ex: si gProfiler couvre déjà
                         partiellement l'espèce) — mêmes colonnes que
                         ref_annot_path mais indexées par protein_id query.
                         Ces annotations ont priorité sur le transfert par
                         orthologie et sont marquées source='direct'.

    Returns
    -------
    DataFrame : protein_id, GO, KEGG, Pfam, source (direct|orthologue),
                confidence (NA pour source=direct), ortholog_id (NA pour
                source=direct), pident, coverage
    """
    ref_annot = pd.read_csv(ref_annot_path, sep="\t").set_index("protein_id")

    merged = rbh_df.merge(ref_annot, left_on="reference_id", right_index=True,
                          how="left")
    merged["source"] = "orthologue (RBH)"
    merged = merged.rename(columns={"reference_id": "ortholog_id"})

    out_cols = ["query_id", "GO", "KEGG", "Pfam", "source", "confidence",
               "ortholog_id", "pident", "coverage"]
    result = merged[out_cols].rename(columns={"query_id": "protein_id"})

    if direct_annot_path:
        direct = pd.read_csv(direct_annot_path, sep="\t")
        direct = direct.rename(columns={"protein_id": "protein_id"})
        direct["source"] = "direct"
        direct["confidence"] = pd.NA
        direct["ortholog_id"] = pd.NA
        direct["pident"] = pd.NA
        direct["coverage"] = pd.NA
        direct = direct[out_cols]

        # Priorite au direct : on retire du transfert les proteines deja
        # couvertes directement, plutot que de dupliquer les lignes
        already_direct = set(direct["protein_id"])
        result = result[~result["protein_id"].isin(already_direct)]
        result = pd.concat([direct, result], ignore_index=True)

    return result.sort_values(["source", "protein_id"]).reset_index(drop=True)


def summarize_coverage(annotation_table: pd.DataFrame, all_query_ids: list) -> dict:
    """
    Résumé de couverture — combien de protéines de l'espèce query ont
    finalement une annotation, et via quelle voie/confiance.
    """
    annotated = set(annotation_table["protein_id"])
    n_total = len(all_query_ids)
    n_annotated = len(annotated)

    by_source = annotation_table["source"].value_counts().to_dict()
    by_confidence = (annotation_table.loc[annotation_table["source"] == "orthologue (RBH)",
                                          "confidence"]
                    .value_counts().to_dict())

    return {
        "n_total": n_total,
        "n_annotated": n_annotated,
        "pct_annotated": round(100 * n_annotated / n_total, 1) if n_total else 0.0,
        "n_unannotated": n_total - n_annotated,
        "by_source": by_source,
        "by_confidence_orthologue": by_confidence,
    }


def print_coverage_summary(summary: dict) -> None:
    print(f"  [STATS] Couverture d'annotation : {summary['n_annotated']}/{summary['n_total']} "
          f"protéines ({summary['pct_annotated']}%)")
    for source, n in summary["by_source"].items():
        print(f"      - {source}: {n}")
    if summary["by_confidence_orthologue"]:
        print("  [STATS] Répartition de confiance (transfert par orthologie) :")
        for conf, n in summary["by_confidence_orthologue"].items():
            print(f"      - {conf}: {n}")
    if summary["n_unannotated"] > 0:
        print(f"  [INFO] {summary['n_unannotated']} protéine(s) sans annotation "
              "(ni directe, ni orthologue RBH trouvé) — probablement "
              "spécifiques à l'espèce ou trop divergentes de la référence.")
