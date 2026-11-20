"""Jeu de plasma synthétique à contamination connue, au format DIA-NN de Deimos.

Utilisé par tests/test_plasma_qc.py et par le job `smoke` de la CI :

    python tests/plasma_fixture.py <dossier> [outlier|confounded|clean]

Modèle : intensité = plasma + h x lysat érythrocytaire + q x plaquettes, puis
consommation des fibrinogènes (fraction c) et libération plaquettaire au caillot.
Un vrai effet biologique (30 protéines hors panels, |log2FC| = 1) sépare B de A.
« RBCTRK » est une protéine érythrocytaire HORS panel : seule la corrélation
intra-condition avec l'indice peut la signaler.

Scénarios (2 conditions x 6 réplicats) :
  outlier    : un échantillon de A fortement hémolysé, les autres propres ;
  confounded : tous les échantillons de B hémolysés (h variable d'un échantillon
               à l'autre), A propre — aucun échantillon aberrant, mais confusion ;
  clean      : aucune contamination au-delà du fond (contrôle négatif).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import plasma_qc  # noqa: E402

CLASSIC = ["ALB", "APOA1", "APOB", "TF", "HP", "A2M", "C3", "SERPINA1", "IGHG1", "IGKC",
           "APOA2", "HPX", "AHSG", "CP", "ORM1", "GC", "AMBP", "ITIH2", "C4B", "VTN",
           "CFB", "PLG", "APOH", "AGT", "SERPINA3", "LRG1", "AZGP1", "RBP4", "TTR", "APOE"]
EFFECT_N = 30


def _first(entries):
    return [e[0] for e in entries if not e[0].startswith("HB") or "-" not in e[0]]


def make_plasma_dataset(scenario: str = "outlier", seed: int = 0, n_rep: int = 6,
                        n_filler: int = 420):
    rng = np.random.default_rng(seed)
    ery = _first(plasma_qc.DEFAULT_PANELS["erythrocyte"])
    plt_ = _first(plasma_qc.DEFAULT_PANELS["platelet"])
    coag = _first(plasma_qc.DEFAULT_PANELS["coagulation"])
    genes = list(dict.fromkeys(CLASSIC + ery + plt_ + coag + ["RBCTRK"]))
    genes += [f"PLS{i:03d}" for i in range(n_filler)]
    n = len(genes)
    gi = {g: i for i, g in enumerate(genes)}

    # Plasma propre (log2) : albumine en tête, décroissance, panels bas
    base = rng.normal(19.5, 2.2, n)
    for k, g in enumerate(CLASSIC):
        base[gi[g]] = 34 - 0.25 * k
    for g in ("FGA", "FGB", "FGG"):
        base[gi[g]] = 30.5
    # ordre trié : itérer un set dépend du hash, randomisé à chaque processus
    for g in sorted(set(ery) | set(plt_) | {"RBCTRK"}):
        base[gi[g]] = rng.normal(17.5, 1.0)
    base[gi["HBB"]], base[gi["HBA1"]] = 21.0, 20.8

    # Profils des contaminants (log2, à l'échelle d'un lysat)
    rbc = np.full(n, -np.inf)
    for g in ery:
        rbc[gi[g]] = rng.normal(27, 1.2)
    rbc[gi["HBB"]], rbc[gi["HBA1"]], rbc[gi["CA1"]] = 35.0, 34.8, 30.0
    rbc[gi["RBCTRK"]] = 28.0
    plt_prof = np.full(n, -np.inf)
    for g in plt_:
        plt_prof[gi[g]] = rng.normal(27, 1.0)

    labels, conds = [], []
    for cond in ("A", "B"):
        for r in range(1, n_rep + 1):
            labels.append(f"Sple-{cond}-{r:02d}")
            conds.append(cond)
    S = len(labels)
    h = np.exp(rng.normal(np.log(2e-4), 0.3, S))      # hémolyse de fond
    q = np.exp(rng.normal(np.log(2e-4), 0.3, S))
    c = np.clip(rng.normal(0.05, 0.02, S), 0, 0.9)    # fibrinogène consommé
    conds_a = np.array(conds)
    if scenario == "outlier":
        h[0] = 1.5e-2
    elif scenario == "confounded":
        h[conds_a == "B"] = np.exp(rng.normal(np.log(3e-3), 0.5, (conds_a == "B").sum()))
    elif scenario == "clean":
        pass
    else:
        raise ValueError(scenario)

    effect_idx = rng.choice([gi[g] for g in genes if g.startswith("PLS")], EFFECT_N, replace=False)
    effect = np.zeros(n)
    effect[effect_idx] = rng.choice([-1.0, 1.0], EFFECT_N)

    lin = np.zeros((n, S))
    for s in range(S):
        b = base + (effect if conds[s] == "B" else 0.0) + rng.normal(0, 0.25, n)
        x = np.exp2(b) + h[s] * np.exp2(rbc) + q[s] * np.exp2(plt_prof)
        for g in ("FGA", "FGB", "FGG"):
            x[gi[g]] *= (1 - c[s])
        lin[:, s] = x
    # Non-détection sous ~2^16 (probabilité logistique)
    lg = np.log2(lin)
    p_det = 1 / (1 + np.exp(-(lg - 16.5) * 2.0))
    lin[rng.uniform(size=lin.shape) > p_det] = 0.0

    cols = [f"C:\\DIA-NN\\{lab}_Slot1-{i + 1:02d}_1_3898.d" for i, lab in enumerate(labels)]
    pg = pd.DataFrame({"Protein.Group": [f"P{50000 + i}" for i in range(n)],
                       "Protein.Ids": [f"P{50000 + i}" for i in range(n)],
                       "Protein.Names": [f"{g}_HUMAN" for g in genes],
                       "Genes": genes,
                       "First.Protein.Description": [f"{g} protein" for g in genes],
                       "Proteotypic": 1})
    pg = pd.concat([pg, pd.DataFrame(np.round(lin), columns=cols)], axis=1)
    design = pd.DataFrame({"label": labels, "condition": conds,
                           "replicate": [int(l[-2:]) for l in labels]})
    truth = {"h": h, "q": q, "c": c, "effect_genes": [genes[i] for i in effect_idx],
             "labels": labels, "conditions": conds, "genes": genes}
    return pg, design, truth


def write_dataset(out_dir: str, scenario: str = "outlier", seed: int = 0):
    os.makedirs(out_dir, exist_ok=True)
    pg, design, truth = make_plasma_dataset(scenario, seed)
    pg.to_csv(os.path.join(out_dir, "report.pg_matrix.tsv"), sep="\t", index=False)
    design.to_csv(os.path.join(out_dir, "ExperimentalDesign.csv"), sep=";", index=False)
    return truth


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else "plasma_example"
    sc = sys.argv[2] if len(sys.argv) > 2 else "outlier"
    write_dataset(d, sc)
    print(f"plasma example ({sc}) written to {d}/")
