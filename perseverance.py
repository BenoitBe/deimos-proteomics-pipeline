#!/usr/bin/env python3
"""
Perseverance — Transfert d'annotation fonctionnelle par orthologie
================================================================
Pour un organisme non-modèle mal couvert par les bases GO/KEGG usuelles
(gProfiler, etc.) : trouve les orthologues dans un protéome de référence
bien annoté (ex: Arabidopsis thaliana pour les Brassicacées) par Reciprocal
Best Hit (DIAMOND), et transfère GO/KEGG/Pfam avec un niveau de confiance
explicite. Sortie directement injectable comme annotation custom dans
go_enrichment.py / topGO / clusterProfiler, à la place de (ou en complément
de) gProfiler.

Usage :
    python perseverance.py --query monorganisme.fasta \\
                      --reference arabidopsis.fasta \\
                      --ref-annotations arabidopsis_go.tsv \\
                      --out annotation_transferee.tsv \\
                      [--direct-annotations partiel_deja_connu.tsv] \\
                      [--min-identity 25] [--min-coverage 50] [--threads 4]
"""
import argparse
import sys

from rbh import compute_rbh
from transfer import transfer_annotations, summarize_coverage, print_coverage_summary


def fasta_ids(path: str) -> list:
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0].strip())
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", required=True, help="FASTA proteines de l'organisme non-modele")
    ap.add_argument("--reference", required=True, help="FASTA proteome de reference bien annote")
    ap.add_argument("--ref-annotations", required=True,
                    help="TSV: protein_id, GO, KEGG, Pfam (GO ';'-separe)")
    ap.add_argument("--direct-annotations", default=None,
                    help="TSV optionnel d'annotations directes deja connues pour le query "
                         "(memes colonnes) — priorite sur le transfert par orthologie")
    ap.add_argument("--out", required=True, help="Table d'annotation finale (TSV)")
    ap.add_argument("--rbh-out", default=None, help="Sauver aussi la table RBH brute (optionnel)")
    ap.add_argument("--workdir", default="perseverance_work", help="Dossier de travail DIAMOND")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--min-identity", type=float, default=25.0,
                    help="%% identite minimum pour retenir un RBH (defaut: 25)")
    ap.add_argument("--min-coverage", type=float, default=50.0,
                    help="%% couverture minimum pour retenir un RBH (defaut: 50)")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    print("[1/3] Calcul des Reciprocal Best Hits (DIAMOND)...")
    rbh_df = compute_rbh(args.query, args.reference, args.workdir,
                         evalue=args.evalue, min_identity=args.min_identity,
                         min_coverage=args.min_coverage, threads=args.threads)
    print(f"  [OK] {len(rbh_df)} orthologue(s) identifie(s) par RBH.")
    if args.rbh_out:
        rbh_df.to_csv(args.rbh_out, sep="\t", index=False)
        print(f"  [OK] Table RBH brute ecrite : {args.rbh_out}")

    print("\n[2/3] Transfert d'annotation...")
    annot = transfer_annotations(rbh_df, args.ref_annotations, args.direct_annotations)
    annot.to_csv(args.out, sep="\t", index=False)
    print(f"  [OK] Table d'annotation ecrite : {args.out}")

    print("\n[3/3] Résumé de couverture...")
    query_ids = fasta_ids(args.query)
    summary = summarize_coverage(annot, query_ids)
    print_coverage_summary(summary)

    if summary["pct_annotated"] < 30:
        print("\n  [WARN] Couverture < 30% : verifier que l'espece de reference "
              "est bien la plus proche disponible, ou assouplir --min-identity/--min-coverage.")


if __name__ == "__main__":
    main()
