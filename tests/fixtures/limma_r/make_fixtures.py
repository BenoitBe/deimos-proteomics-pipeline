"""Régénère les jeux de référence R de tests/test_limma_reference.py (développement
seulement ; les tests eux-mêmes n'ont besoin ni de R ni de réseau).

  python make_fixtures.py                      # données simulées (graines fixes)
  DEQMS_R=/chemin/DEqMS/R/DEqMS.R Rscript run_reference.R case_*   # sorties R

Références : R 4.3.3, limma 3.58.1, DEqMS 1.31.0 (fonction spectraCounteBayes
sourcée depuis le code du paquet). Données entièrement simulées.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..")))
from limma_ebayes import make_all_contrasts, make_design_matrix, make_paired_design_matrix  # noqa: E402


def save(name, E, X, C, counts=None):
    d = os.path.join(HERE, name)
    os.makedirs(d, exist_ok=True)
    np.savetxt(os.path.join(d, "expr.csv.gz"), E, delimiter=",", fmt="%.6f")
    np.savetxt(os.path.join(d, "design.csv"), X, delimiter=",", fmt="%.17g")
    np.savetxt(os.path.join(d, "contrasts.csv"), C, delimiter=",", fmt="%.17g")
    if counts is not None:
        np.savetxt(os.path.join(d, "counts.csv"), counts, delimiter=",", fmt="%d")


def sim(n, groups, seed, d0=4.0, s0=0.04, de=0.1, counts=False):
    rng = np.random.default_rng(seed)
    X, g = make_design_matrix(groups)
    C, _ = make_all_contrasts(g)
    k = None
    if counts:
        k = rng.geometric(0.15, n) + 1
        s2 = s0 * (0.3 + 1.5 / k) * d0 / rng.chisquare(d0, n)
    else:
        s2 = s0 * d0 / rng.chisquare(d0, n)
    mu = rng.normal(20, 2, n)
    eff = np.where(rng.random(n) < de, rng.choice([-1, 1], n) * rng.uniform(0.5, 2, n), 0)
    gi = np.argmax(X, axis=1)
    E = (mu[:, None] + eff[:, None] * (gi[None, :] == 0)
         + rng.normal(0, 1, (n, len(groups))) * np.sqrt(s2)[:, None])
    return E, X, C, k


if __name__ == "__main__":
    save("case_2g", *sim(300, ["A"] * 3 + ["B"] * 3, 11)[:3])
    save("case_4g", *sim(300, sum([[g] * 3 for g in "ABCD"], []), 12)[:3])
    E, X, C, k = sim(600, ["A"] * 4 + ["B"] * 4, 13, counts=True)
    save("case_deqms", E, X, C, k)
    rng = np.random.default_rng(14)
    cond, subj = ["T", "C"] * 6, [f"s{i // 2}" for i in range(12)]
    X, g, cols = make_paired_design_matrix(cond, subj)
    C, _ = make_all_contrasts(g, n_total_cols=len(cols))
    E = (rng.normal(20, 2, 200)[:, None] + rng.normal(0, 1, (200, 6))[:, [i // 2 for i in range(12)]]
         + rng.normal(0, 0.3, (200, 12)) + (rng.random(200) < 0.1)[:, None] * np.array([1.0, 0] * 6))
    save("case_paired", E, X, C)
    E, X, C, _ = sim(200, ["A"] * 3 + ["B"] * 3, 15)
    E[:5] = 18.0                                         # variances exactement nulles
    save("case_zero", E, X, C)
    rng = np.random.default_rng(16)
    X, g = make_design_matrix(["A"] * 3 + ["B"] * 3)
    C, _ = make_all_contrasts(g)
    r = rng.normal(0, 0.2, 6)
    r -= np.r_[np.full(3, r[:3].mean()), np.full(3, r[3:].mean())]
    E = (rng.normal(20, 2, 40)[:, None] + r[None, :]
         + (np.arange(40) < 5)[:, None] * np.r_[np.ones(3), np.zeros(3)])
    save("case_homog", E, X, C)                          # d0 = Inf
