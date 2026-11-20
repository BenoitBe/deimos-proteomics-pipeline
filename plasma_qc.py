#!/usr/bin/env python3
"""plasma_qc.py — qualité des échantillons de plasma / sérum (panels de Geyer et al. 2019).

Trois sources d'artefact dominent les études de plasma : la lyse des érythrocytes,
la contamination plaquettaire et la coagulation (prélèvement, délai, centrifugation).
Geyer et al. (EMBO Mol Med 2019, 11:e10427) en ont tiré trois panels de marqueurs ;
la moitié des études de biomarqueurs qu'ils ont relues proposaient des protéines de
ces panels comme candidats.

Ce module :
  1. calcule, par échantillon, un indice par panel (définitions de l'article) :
       érythrocytes, plaquettes = somme des intensités du panel / somme de toutes
                                  les protéines quantifiées ;
       coagulation              = somme de toutes les protéines / somme des
                                  fibrinogènes (FGA, FGB, FGG), consommés par le
                                  caillot — l'article n'utilise que les protéines
                                  qui diminuent ;
     les trois sont orientés de la même façon : plus haut = plus contaminé ;
  2. signale les échantillons hors norme : règle de l'article (log2 de l'indice
     > moyenne + 2 SD) et règle robuste (z modifié médiane/MAD > 3,5, Iglewicz &
     Hoaglin), qui ne se laisse pas masquer par l'échantillon contaminé lui-même ;
  3. teste si un indice diffère entre conditions (Kruskal-Wallis, epsilon^2 ; par
     contraste, Welch sur le log2) : une contamination confondue avec la condition
     ne produit AUCUN échantillon aberrant, mais des DEP artefactuelles ;
  4. pour chaque DEP : appartenance à un panel, corrélation intra-condition avec
     l'indice (moyennes de condition retirées des deux côtés : on mesure la
     co-variation d'échantillon à échantillon, pas l'effet de la condition), et
     sens du changement comparé à celui de l'indice entre les deux groupes.

Branché dans Deimos par la clé `plasma_qc: true` (feuilles Plasma_QC et
Plasma_QC_DEP du classeur). Utilisable aussi sur un run déjà terminé :

    python plasma_qc.py proteomics_output

Aucune dépendance au-delà de celles de Deimos (numpy, pandas, scipy, matplotlib,
openpyxl).

Provenance des panels : Table EV2 de l'article (30 protéines par panel). Le
fichier n'étant pas téléchargeable automatiquement, la liste a été transcrite
depuis la copie qu'en donne CAT-APP (The Hong Wang Lab, GitHub), puis les deux
exclusions énoncées dans le texte de l'article ont été appliquées (NIF3L1 retiré
du panel érythrocytaire, F13A1 du panel plaquettaire). À vérifier contre la
Table EV2 avant publication ; `plasma_qc_panels` permet de fournir sa propre liste.
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

PANEL_SOURCE = "Geyer et al. 2019, EMBO Mol Med 11:e10427 (Table EV2)"
PANEL_ORDER = ("erythrocyte", "platelet", "coagulation")

# Une entrée = un groupe protéique ; les synonymes d'une même entrée sont les
# gènes qui partagent ce groupe dans DIA-NN (ubiquitine, EIF5A...) ou un symbole
# renommé depuis 2019 (SDPR -> CAVIN2). Correspondance insensible à la casse.
# Hémoglobines, seuls symboles qui changent d'une espèce à l'autre (tous les autres
# membres des panels ont le même symbole chez la souris, le rat et le chien, à la
# casse près) : souris Hba-a1/a2, Hbb-bs/bt/b1/b2 ; rat Hba-a1/a2/a3, Hba1, Hbb,
# Hbb-b1 ; chien HBA, HBB (UniProt). Sans effet sur l'humain.
_ERY = (("HBA1", "HBA2", "HBA", "HBA-A1", "HBA-A2", "HBA-A3"),
        ("HBB", "HBB-BS", "HBB-BT", "HBB-B1", "HBB-B2"), ("CA1",), ("HBD",), ("PRDX2",),
        ("CA2",), ("CAT",), ("SLC4A1",), ("BLVRB",), ("SPTA1",), ("SPTB",), ("ANK1",),
        ("GAPDH",), ("SOD1",), ("PRDX6",), ("BPGM",), ("ACTB", "ACTG1"),
        ("SELENBP1",), ("EPB41",), ("LDHB",), ("EIF5A", "EIF5A2", "EIF5AL1"),
        ("PNP",), ("ALDOA",), ("TPI1",),
        ("UBC", "UBB", "RPS27A", "UBA52", "UBBP4"), ("HBE1",), ("HSPA8",),
        ("PGK1",), ("ALAD",))
_PLT = (("FLNA",), ("TLN1",), ("MYH9",), ("ACTB",), ("VCL",), ("ACTN1",),
        ("TPM4",), ("THBS1",), ("TUBB1",), ("YWHAZ",), ("GSN",), ("TUBA1B",),
        ("ITGA2B",), ("PFN1",), ("ITGB3",), ("TAGLN2",), ("FERMT3",),
        ("RAP1B",), ("PLEK",), ("PPBP",), ("GAPDH",), ("MMRN1",), ("MYL6",),
        ("CFL1",), ("PARVB",), ("SDPR", "CAVIN2"), ("TUBB4B",), ("TMSB4X",),
        ("PKM",))
_COAG = (("FGB",), ("FGG",), ("FGA",), ("F13A1",), ("SERPINC1",), ("PPBP",),
         ("F2",), ("GP1BA",), ("ECM1",), ("CLU",), ("DSP",), ("WDR1",),
         ("ATRN",), ("GP5",), ("SERPINA5",), ("C1RL",), ("MAN1A1",), ("F11",),
         ("KNG1",), ("BCHE",), ("GPLD1",), ("THBS1",), ("MASP1",), ("TENM4",),
         ("HSPG2",), ("APOC3",), ("CDH5",), ("LTF",), ("CST3",), ("PROC",),
         ("PF4", "PF4V1"))

DEFAULT_PANELS = {"erythrocyte": _ERY, "platelet": _PLT, "coagulation": _COAG}
# Protéines qui entrent dans l'indice de coagulation : celles que le caillot consomme
FIBRINOGEN = (("FGA",), ("FGB",), ("FGG",))

SD_K = 2.0               # règle de l'article : > moyenne + 2 SD
ROBUST_Z = 3.5           # z modifié (Iglewicz & Hoaglin 1993)
R_THRESHOLD = 0.5        # seuil de corrélation de l'article
MIN_PANEL_FOUND = 3      # en dessous, l'indice n'est pas calculé
MIN_PAIRS = 6            # paires minimales pour une corrélation intra-condition
ALPHA = 0.05

INDEX_LABEL = {"erythrocyte": "Erythrocyte lysis", "platelet": "Platelet contamination",
               "coagulation": "Coagulation (fibrinogen loss)"}


# ─────────────────────────── panels ─────────────────────────────────────────

def gene_tokens(genes) -> set:
    """Symboles d'une cellule Genes de DIA-NN ('UBB;UBC' -> {'UBB', 'UBC'})."""
    if not isinstance(genes, str) or not genes.strip():
        return set()
    return {t.strip().upper() for t in re.split(r"[;,|\s]+", genes) if t.strip()}


def _entries_to_sets(entries) -> list:
    return [frozenset(str(a).upper() for a in e) for e in entries]


def load_panels(path: "str | None" = None) -> dict:
    """Panels par défaut, ou fichier utilisateur (TSV/CSV, colonnes gene + panel ;
    plusieurs synonymes d'une entrée séparés par ';' dans la colonne gene)."""
    if not path:
        return {k: _entries_to_sets(v) for k, v in DEFAULT_PANELS.items()}
    sep = "," if str(path).lower().endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str)
    cols = {c.lower().strip(): c for c in df.columns}
    if "gene" not in cols or "panel" not in cols:
        raise ValueError(f"{path}: columns 'gene' and 'panel' are required "
                         f"(found {list(df.columns)})")
    out = {}
    for _, r in df.iterrows():
        p = str(r[cols["panel"]]).strip().lower()
        syn = frozenset(t.strip().upper() for t in str(r[cols["gene"]]).split(";") if t.strip())
        if syn:
            out.setdefault(p, []).append(syn)
    unknown = sorted(set(out) - set(PANEL_ORDER))
    if unknown:
        raise ValueError(f"{path}: unknown panel(s) {unknown}; expected {list(PANEL_ORDER)}")
    return out


def index_entries(panels: dict, p: str) -> list:
    """Entrées qui entrent dans la SOMME de l'indice p.
    Érythrocytes / plaquettes : membres propres au panel. Les protéines communes aux
    deux panels (ACTB, GAPDH dans la Table EV2) sont retirées des deux sommes : sinon
    une hémolysie fait monter l'indice plaquettaire, et inversement (écart assumé à
    l'article, qui les garde). Coagulation : fibrinogènes (FGA, FGB, FGG)."""
    if p == "coagulation":
        return _entries_to_sets(FIBRINOGEN)
    other = "platelet" if p == "erythrocyte" else "erythrocyte"
    other_alias = set().union(*panels.get(other, [])) if panels.get(other) else set()
    return [e for e in panels.get(p, []) if not (e & other_alias)]


def membership(genes, entries) -> np.ndarray:
    """Vecteur booléen : la protéine (groupe) appartient-elle à un des groupes du panel."""
    alias = set().union(*entries) if entries else set()
    return np.array([bool(gene_tokens(g) & alias) for g in genes], dtype=bool)


def panel_coverage(genes, panels: dict) -> pd.DataFrame:
    """Une ligne par entrée de panel : quantifiée ou non, et sous quel groupe."""
    toks = [gene_tokens(g) for g in genes]
    rows = []
    for p in PANEL_ORDER:
        for e in panels.get(p, []):
            hit = [str(genes[i]) for i, t in enumerate(toks) if t & e]
            rows.append({"panel": p, "entry": ";".join(sorted(e)),
                         "quantified": bool(hit), "protein_group_genes": " | ".join(hit)})
    return pd.DataFrame(rows)


# ─────────────────────────── ces échantillons sont-ils du plasma ? ───────────

# Protéines plasmatiques classiques (préfixes : sous-familles de souris et chaînes
# d'immunoglobulines). Sert à proposer l'étape, jamais à la lancer seule.
_PLASMA_EXACT = {"ALB", "APOA1", "APOA2", "APOB", "APOH", "APOE", "APOC3", "TF", "TRF",
                 "HP", "HPX", "A2M", "PZP", "MUG1", "C3", "FGA", "FGB", "FGG", "ORM1",
                 "ORM2", "AHSG", "CP", "GC", "AMBP", "PLG", "TTR", "CFB", "VTN", "RBP4",
                 "CLU", "KNG1", "HRG", "FN1", "AGT", "C4B", "C4A", "IGKC", "IGHM", "LRG1"}
_PLASMA_PREFIX = ("SERPINA1", "SERPINA3", "ITIH", "IGHG", "IGHA", "IGLC", "IGHV", "IGKV")


def _is_plasma_protein(genes) -> bool:
    toks = gene_tokens(genes)
    return bool(toks & _PLASMA_EXACT) or any(t.startswith(_PLASMA_PREFIX) for t in toks)


def plasma_signature(mat_log2: pd.DataFrame, genes, top: int = 20,
                     min_plasma_top: int = 8, min_alb_share: float = 0.05) -> dict:
    """Signature d'un protéome de plasma / sérum : l'albumine est la protéine la plus
    abondante (médiane des échantillons), porte au moins min_alb_share de l'intensité
    totale, et au moins min_plasma_top des `top` protéines les plus abondantes sont des
    protéines plasmatiques classiques. Le LCR passe aussi (albumine d'origine
    plasmatique) : le panel érythrocytaire y a un sens, la contamination sanguine du
    LCR étant le problème classique."""
    genes = list(genes)
    I = np.exp2(mat_log2.to_numpy(float))
    I[~np.isfinite(I)] = 0.0
    med = np.median(I, axis=1)
    tot = med.sum()
    order = np.argsort(-med)
    alb = [i for i, g in enumerate(genes) if "ALB" in gene_tokens(g)]
    alb_rank = int(np.where(order == alb[0])[0][0]) + 1 if alb else None
    alb_share = float(med[alb].sum() / tot) if alb and tot > 0 else 0.0
    n_top = int(sum(_is_plasma_protein(genes[i]) for i in order[:top]))
    looks = bool(alb_rank == 1 and alb_share >= min_alb_share and n_top >= min_plasma_top)
    ev = (f"albumin rank {alb_rank if alb_rank else 'n/a'}, {100 * alb_share:.0f}% of the "
          f"intensity; {n_top}/{top} most abundant proteins are classic plasma proteins")
    return {"looks_like": looks, "alb_rank": alb_rank, "alb_share": alb_share,
            "n_plasma_top": n_top, "evidence": ev}


# ─────────────────────────── indices ────────────────────────────────────────

def compute_indices(mat_log2: pd.DataFrame, genes, panels: dict) -> tuple:
    """Indices par échantillon. mat_log2 : protéines × échantillons, NaN = absent.
    Retourne (DataFrame échantillons × indices linéaires, dict d'infos)."""
    I = np.exp2(mat_log2.to_numpy(float))
    I[~np.isfinite(I)] = 0.0
    total = I.sum(axis=0)
    genes = list(genes)
    out = pd.DataFrame(index=[str(c) for c in mat_log2.columns])
    info = {}
    for p in ("erythrocyte", "platelet"):
        m = membership(genes, index_entries(panels, p))
        n_found = int(m.sum())
        info[p] = {"n_entries": len(panels.get(p, [])),
                   "n_found": int(membership(genes, panels.get(p, [])).sum()),
                   "n_in_index": n_found,
                   "n_shared_excluded": len(panels.get(p, [])) - len(index_entries(panels, p))}
        if n_found < MIN_PANEL_FOUND:
            out[p] = np.nan
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            out[p] = np.where(total > 0, I[m].sum(axis=0) / total, np.nan)
    m = membership(genes, _entries_to_sets(FIBRINOGEN))
    info["coagulation"] = {"n_entries": len(panels.get("coagulation", [])),
                           "n_found": int(membership(genes, panels.get("coagulation", [])).sum()),
                           "n_fibrinogen": int(m.sum())}
    fib = I[m].sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        out["coagulation"] = np.where(fib > 0, total / fib, np.nan) if m.any() else np.nan
    return out, info


def flag_samples(idx: pd.DataFrame, k_sd: float = SD_K, robust_z: float = ROBUST_Z) -> pd.DataFrame:
    """Statistiques par échantillon sur le log2 des indices, règle de l'article
    (moyenne + k SD, ddof = 1) et règle robuste (z modifié > robust_z), unilatérales."""
    res = pd.DataFrame(index=idx.index)
    for p in PANEL_ORDER:
        x = np.log2(idx[p].to_numpy(float)) if p in idx else np.full(len(idx), np.nan)
        ok = np.isfinite(x)
        res[f"{p}_index"] = idx[p] if p in idx else np.nan
        res[f"{p}_log2"] = x
        z = np.full(len(x), np.nan)
        rz = np.full(len(x), np.nan)
        if ok.sum() >= 3:
            mu, sd = x[ok].mean(), x[ok].std(ddof=1)
            if sd > 0:
                z[ok] = (x[ok] - mu) / sd
            med = np.median(x[ok])
            mad = 1.4826 * np.median(np.abs(x[ok] - med))
            if mad > 0:
                rz[ok] = (x[ok] - med) / mad
        res[f"{p}_z"] = z
        res[f"{p}_robust_z"] = rz
        res[f"{p}_flag_sd"] = np.nan_to_num(z, nan=-np.inf) > k_sd
        res[f"{p}_flag_robust"] = np.nan_to_num(rz, nan=-np.inf) > robust_z
    flags = [f"{p}_flag_{r}" for p in PANEL_ORDER for r in ("sd", "robust")]
    res["flagged"] = res[flags].any(axis=1)
    res["flag_reason"] = [
        ", ".join(f"{p} ({'+'.join(r for r in ('sd', 'robust') if res.at[s, f'{p}_flag_{r}'])})"
                  for p in PANEL_ORDER
                  if res.at[s, f"{p}_flag_sd"] or res.at[s, f"{p}_flag_robust"])
        for s in res.index]
    return res


def max_z_reachable(n: int) -> float:
    """Plus grand z (ddof=1) qu'un échantillon peut atteindre parmi n : (n-1)/sqrt(n).
    Pour n <= 5 il reste sous 2 — la règle moyenne + 2 SD ne peut alors rien signaler."""
    return (n - 1) / np.sqrt(n) if n > 1 else 0.0


# ─────────────────────────── confusion avec la condition ────────────────────

def _kruskal_eps2(values, groups) -> tuple:
    """Kruskal-Wallis + epsilon^2 (Tomczak & Tomczak 2014), même formule que Deimos."""
    from scipy.stats import kruskal
    v = np.asarray(values, float)
    g = np.asarray(groups, dtype=object)
    ok = np.isfinite(v)
    v, g = v[ok], g[ok]
    levels = [lv for lv in dict.fromkeys(g) if (g == lv).sum() >= 1]
    samples = [v[g == lv] for lv in levels]
    n = len(v)
    if len(samples) < 2 or n < 3 or all(s.size < 2 for s in samples):
        return np.nan, np.nan, np.nan
    try:
        H, p = kruskal(*samples)
    except ValueError:            # toutes valeurs identiques
        return np.nan, np.nan, np.nan
    eps2 = float(np.clip(H * (n + 1) / (n ** 2 - 1), 0.0, 1.0))
    return float(H), float(p), eps2


def confounding(flags: pd.DataFrame, conditions, contrast_names=None) -> pd.DataFrame:
    """Un indice diffère-t-il entre conditions ? Kruskal-Wallis global, puis par
    contraste A_vs_B : différence des moyennes du log2 (A - B) et test de Welch."""
    from scipy.stats import ttest_ind
    # Les noms de contrastes de Deimos remplacent '.' par '_' (make_all_contrasts)
    cond = np.asarray([str(c).replace(".", "_") for c in conditions], dtype=object)
    rows = []
    for p in PANEL_ORDER:
        x = flags[f"{p}_log2"].to_numpy(float)
        H, pk, e2 = _kruskal_eps2(x, cond)
        rows.append({"index": p, "comparison": "all conditions (Kruskal-Wallis)",
                     "delta_log2": np.nan, "p_value": pk, "epsilon2": e2})
        for c in (contrast_names or []):
            if "_vs_" not in c:
                continue
            a, b = c.split("_vs_", 1)
            xa = x[(cond == a) & np.isfinite(x)]
            xb = x[(cond == b) & np.isfinite(x)]
            d = xa.mean() - xb.mean() if xa.size and xb.size else np.nan
            pv = np.nan
            if xa.size >= 2 and xb.size >= 2 and (xa.std() > 0 or xb.std() > 0):
                pv = float(ttest_ind(xa, xb, equal_var=False).pvalue)
            rows.append({"index": p, "comparison": c, "delta_log2": d,
                         "p_value": pv, "epsilon2": np.nan})
    return pd.DataFrame(rows)


# ─────────────────────────── corrélation intra-condition ────────────────────

def within_condition_r(Y: np.ndarray, x: np.ndarray, groups) -> tuple:
    """Corrélation de Pearson entre chaque ligne de Y (protéines × échantillons, NaN
    autorisé) et x (échantillons), après retrait des moyennes de condition des deux
    côtés, calculées sur les échantillons où la protéine est observée. C'est la
    corrélation partielle à la condition : elle mesure la co-variation d'échantillon
    à échantillon, pas l'effet de la condition.
    Retourne (r, n_paires, ddl) ; ddl = n_paires - n_conditions_utilisées - 1."""
    Y = np.asarray(Y, float)
    x = np.asarray(x, float)
    g = np.asarray(groups, dtype=object)
    n_prot = Y.shape[0]
    Yr = np.full(Y.shape, np.nan)
    Xr = np.full(Y.shape, np.nan)
    n_groups_used = np.zeros(n_prot, int)
    for lv in dict.fromkeys(g):
        cols = np.where(g == lv)[0]
        Yg = Y[:, cols]
        Xg = np.broadcast_to(x[cols], Yg.shape)
        M = np.isfinite(Yg) & np.isfinite(Xg)
        c = M.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            ym = np.where(M, Yg, 0.0).sum(axis=1) / c
            xm = np.where(M, Xg, 0.0).sum(axis=1) / c
        use = M & (c >= 2)[:, None]
        Yr[:, cols] = np.where(use, Yg - ym[:, None], np.nan)
        Xr[:, cols] = np.where(use, Xg - xm[:, None], np.nan)
        n_groups_used += (c >= 2)
    M = np.isfinite(Yr) & np.isfinite(Xr)
    n = M.sum(axis=1)
    yr = np.where(M, Yr, 0.0)
    xr = np.where(M, Xr, 0.0)
    num = (yr * xr).sum(axis=1)
    den = np.sqrt((yr ** 2).sum(axis=1) * (xr ** 2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where((den > 0) & (n >= MIN_PAIRS), num / den, np.nan)
    df = n - n_groups_used - 1
    r = np.where(df >= 2, r, np.nan)
    return r, n, df


def r_pvalue(r, df, two_sided: bool = False) -> np.ndarray:
    """p-value d'une corrélation partielle (t à df ddl) ; unilatérale (r > 0) par défaut."""
    from scipy.stats import t as tdist
    r = np.asarray(r, float)
    df = np.asarray(df, float)
    out = np.full(r.shape, np.nan)
    ok = np.isfinite(r) & (df >= 1) & (np.abs(r) < 1)
    tt = r[ok] * np.sqrt(df[ok] / (1 - r[ok] ** 2))
    out[ok] = 2 * tdist.sf(np.abs(tt), df[ok]) if two_sided else tdist.sf(tt, df[ok])
    out[np.isfinite(r) & (np.abs(r) >= 1)] = 0.0
    return out


def protein_correlations(mat_obs: pd.DataFrame, genes, flags: pd.DataFrame,
                         conditions, panels: dict, mat_all: "pd.DataFrame | None" = None,
                         genes_all=None) -> pd.DataFrame:
    """Corrélation intra-condition de chaque protéine de mat_obs avec chaque indice.
    Pour une protéine qui appartient elle-même au panel, l'indice est recalculé
    SANS elle (sinon HBB, qui domine la somme érythrocytaire, corrèle à ~1 avec
    l'indice par construction)."""
    samples = [str(c) for c in mat_obs.columns]
    Y = mat_obs.to_numpy(float)
    genes = list(genes)
    out = pd.DataFrame(index=range(len(Y)))
    # Base des indices « sans la protéine » : intensités de la matrice d'origine
    base = mat_all if mat_all is not None else mat_obs
    gbase = list(genes_all) if genes_all is not None else genes
    Ib = np.exp2(base[samples].to_numpy(float))
    Ib[~np.isfinite(Ib)] = 0.0
    total = Ib.sum(axis=0)
    Iobs = np.exp2(Y)
    Iobs[~np.isfinite(Iobs)] = 0.0
    for p in PANEL_ORDER:
        x = flags.loc[samples, f"{p}_log2"].to_numpy(float)
        r, n, df = within_condition_r(Y, x, conditions)
        entries = index_entries(panels, p)
        mem = membership(genes, entries)
        if mem.any() and np.isfinite(x).any():
            num_b = Ib[membership(gbase, entries)].sum(axis=0)
            for i in np.where(mem)[0]:
                with np.errstate(divide="ignore", invalid="ignore"):
                    if p == "coagulation":
                        xi = np.log2(total / (num_b - Iobs[i]))
                    else:
                        xi = np.log2((num_b - Iobs[i]) / total)
                xi[~np.isfinite(xi)] = np.nan
                ri, ni, dfi = within_condition_r(Y[i:i + 1], xi, conditions)
                r[i], n[i], df[i] = ri[0], ni[0], dfi[0]
        out[f"r_within_{p}"] = r
        pr = r_pvalue(r, df, two_sided=(p == "coagulation"))
        out[f"p_r_{p}"] = pr
        # BH sur toutes les protéines testées, par indice : trois indices x des
        # centaines de protéines, sans correction un DEP sur sept « suivrait » un
        # indice par hasard (vérifié sur simulation)
        out[f"q_r_{p}"] = bh(pr)
        out[f"n_pairs_{p}"] = n
    return out


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg sur les valeurs finies ; NaN conservés."""
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = int(ok.sum())
    if m == 0:
        return q
    order = np.argsort(p[ok])
    ranked = p[ok][order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty(m)
    tmp[order] = np.clip(ranked, 0, 1)
    q[ok] = tmp
    return q


# ─────────────────────────── DEP ────────────────────────────────────────────

def dep_flags(df_results: pd.DataFrame, contrast_names, sig_masks: dict,
              genes, corr: pd.DataFrame, conf: pd.DataFrame, panels: dict,
              r_thr: float = R_THRESHOLD, alpha: float = ALPHA) -> pd.DataFrame:
    """Une ligne par DEP et par contraste. Niveaux :
      'high'  = indice qui diffère entre les deux groupes (Welch p < alpha) ET protéine
                du panel changeant dans le même sens, ou qui suit l'indice ;
      'check' = membre d'un panel, ou suit un indice (r intra-condition >= r_thr et
                q < alpha, BH sur toutes les protéines testées pour cet indice) ;
      ''      = rien à signaler."""
    genes = list(genes)
    mem = {p: membership(genes, panels.get(p, [])) for p in PANEL_ORDER}
    name_col = "name" if "name" in df_results.columns else df_results.columns[0]
    rows = []
    for c in contrast_names:
        mask = np.asarray(sig_masks.get(c, np.zeros(len(df_results), bool)), bool)
        if not mask.any():
            continue
        lfc = df_results[f"{c}_diff"].to_numpy(float) if f"{c}_diff" in df_results else \
            np.full(len(df_results), np.nan)
        pv = df_results[f"{c}_p.val"].to_numpy(float) if f"{c}_p.val" in df_results else \
            np.full(len(df_results), np.nan)
        cf = conf[conf["comparison"] == c].set_index("index")
        for i in np.where(mask)[0]:
            reasons, level = [], ""
            members = [p for p in PANEL_ORDER if mem[p][i]]
            for p in members:
                reasons.append(f"member of {p} panel")
            tracking = []
            for p in PANEL_ORDER:
                r = corr.at[i, f"r_within_{p}"]
                qr = corr.at[i, f"q_r_{p}"]
                strong = (abs(r) if p == "coagulation" else r) >= r_thr
                if np.isfinite(r) and strong and np.isfinite(qr) and qr < alpha:
                    tracking.append(p)
                    reasons.append(f"tracks {p} index (within-condition r = {r:.2f})")
            if members or tracking:
                level = "check"
            for p in [q for q in PANEL_ORDER if q in members or q in tracking]:
                if p not in cf.index:
                    continue
                d, pd_ = cf.at[p, "delta_log2"], cf.at[p, "p_value"]
                if not (np.isfinite(d) and np.isfinite(pd_) and pd_ < alpha):
                    continue
                if p == "coagulation":
                    same = p in tracking
                else:
                    same = (np.sign(lfc[i]) == np.sign(d)) and np.isfinite(lfc[i])
                if same:
                    level = "high"
                    reasons.append(f"{p} index differs between groups "
                                   f"(delta log2 = {d:+.2f}, p = {pd_:.2g}), same direction")
            rows.append({"contrast": c, "protein": df_results.iloc[i][name_col],
                         "Genes": genes[i], "log2FC": lfc[i], "p.value": pv[i],
                         "panel_member": ", ".join(members),
                         **{f"r_within_{p}": corr.at[i, f"r_within_{p}"] for p in PANEL_ORDER},
                         "flag": level, "reason": "; ".join(reasons)})
    return pd.DataFrame(rows, columns=["contrast", "protein", "Genes", "log2FC", "p.value",
                                       "panel_member",
                                       *[f"r_within_{p}" for p in PANEL_ORDER],
                                       "flag", "reason"])


# ─────────────────────────── figure ─────────────────────────────────────────

def plot_indices(flags: pd.DataFrame, conditions, conf: pd.DataFrame, out_png: str) -> "str | None":
    """Log2 des trois indices par condition ; seuils des deux règles ; échantillons signalés en rouge."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    cond = np.asarray(conditions, dtype=object)
    levels = list(dict.fromkeys(cond))
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, 3, figsize=(4.2 * 3, 4.4))
    for ax, p in zip(axes, PANEL_ORDER):
        x = flags[f"{p}_log2"].to_numpy(float)
        ok = np.isfinite(x)
        if not ok.any():
            ax.text(0.5, 0.5, "not computable\n(too few panel proteins)",
                    ha="center", va="center", transform=ax.transAxes, color="#7f7f7f")
            ax.set_title(INDEX_LABEL[p], fontsize=10)
            ax.set_xticks([])
            continue
        fl = flags[f"{p}_flag_sd"].to_numpy(bool) | flags[f"{p}_flag_robust"].to_numpy(bool)
        for j, lv in enumerate(levels):
            m = (cond == lv) & ok
            jit = rng.uniform(-0.18, 0.18, m.sum())
            ax.scatter(j + jit, x[m], s=26, c=np.where(fl[m], "#c0392b", "#5b7fa6"),
                       edgecolors="white", linewidths=0.5, zorder=3)
            for xx, yy, lab, f in zip(j + jit, x[m], flags.index[m], fl[m]):
                if f:
                    ax.annotate(str(lab), (xx, yy), fontsize=7, color="#c0392b",
                                xytext=(4, 2), textcoords="offset points")
        mu, sd = x[ok].mean(), x[ok].std(ddof=1) if ok.sum() > 1 else 0.0
        med = np.median(x[ok])
        mad = 1.4826 * np.median(np.abs(x[ok] - med))
        if sd > 0:
            ax.axhline(mu + SD_K * sd, ls="--", lw=0.9, color="#7f7f7f",
                       label=f"mean + {SD_K:g} SD")
        if mad > 0:
            ax.axhline(med + ROBUST_Z * mad, ls=":", lw=0.9, color="#7f7f7f",
                       label=f"median + {ROBUST_Z:g} MAD-SD")
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels(levels, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("log2 index (higher = more contaminated)", fontsize=8)
        row = conf[(conf["index"] == p) & conf["comparison"].str.startswith("all")]
        pk = row["p_value"].iloc[0] if len(row) else np.nan
        sub = f"Kruskal-Wallis p = {pk:.2g}" if np.isfinite(pk) else ""
        ax.set_title(f"{INDEX_LABEL[p]}\n{sub}", fontsize=10)
        ax.legend(fontsize=7, loc="best", frameon=False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle("Plasma sample quality — " + PANEL_SOURCE, fontsize=9, y=1.0)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_png


# ─────────────────────────── orchestration ──────────────────────────────────

def run_sample_qc(mat_log2: pd.DataFrame, genes, conditions, panels_path=None,
                  contrast_names=None, out_dir: "str | None" = None, verbose: bool = True,
                  matrix_label: str = "unfiltered matrix, before imputation") -> dict:
    """Étape 1 (avant les statistiques) : indices, échantillons signalés, confusion."""
    panels = load_panels(panels_path)
    idx, info = compute_indices(mat_log2, genes, panels)
    flags = flag_samples(idx)
    flags.insert(0, "condition", list(conditions))
    conf = confounding(flags, conditions, contrast_names)
    cover = panel_coverage(list(genes), panels)
    png = None
    if out_dir:
        png = plot_indices(flags, conditions, conf, os.path.join(out_dir, "plasma_qc_indices.png"))
    n = int(np.isfinite(flags["erythrocyte_log2"]).sum() or len(flags))
    res = {"panels": panels, "panels_path": panels_path, "indices": idx, "flags": flags,
           "confounding": conf, "coverage": cover, "info": info, "plot": png,
           "n_samples": len(flags), "max_z": max_z_reachable(n),
           "matrix_label": matrix_label, "out_dir": out_dir}
    res["confounded"] = [r["index"] for _, r in conf.iterrows()
                         if r["comparison"].startswith("all")
                         and np.isfinite(r["p_value"]) and r["p_value"] < ALPHA]
    if verbose:
        _print_sample_summary(res)
    return res


def run_dep_qc(res: dict, mat_obs: pd.DataFrame, genes_obs, df_results: pd.DataFrame,
               contrast_names, sig_masks: dict, conditions, mat_all=None, genes_all=None,
               verbose: bool = True) -> dict:
    """Étape 2 (après le différentiel) : confusion par contraste, corrélations et
    drapeaux sur les DEP."""
    conf = confounding(res["flags"], res["flags"]["condition"].tolist(), contrast_names)
    corr = protein_correlations(mat_obs, genes_obs, res["flags"], conditions, res["panels"],
                                mat_all=mat_all, genes_all=genes_all)
    deps = dep_flags(df_results, contrast_names, sig_masks, genes_obs, corr,
                     conf, res["panels"])
    res = {**res, "confounding": conf, "correlations": corr, "dep": deps}
    res["dep_summary"] = (deps.assign(n=1).pivot_table(index="contrast", columns="flag",
                                                         values="n", aggfunc="sum", fill_value=0)
                          if len(deps) else pd.DataFrame())
    if verbose and len(deps):
        for c, sub in deps.groupby("contrast", sort=False):
            n_hi = int((sub["flag"] == "high").sum())
            n_ck = int((sub["flag"] == "check").sum())
            n_mem = int((sub["panel_member"] != "").sum())
            if n_hi or n_ck:
                print(f"  [PLASMA-QC] {c}: {n_hi} DEP 'high', {n_ck} 'check' "
                      f"(of {len(sub)} DEP; {n_mem} panel member(s))")
            if n_mem >= 3 and n_mem >= 0.1 * len(sub):
                print(f"  [WARN] {c}: {n_mem}/{len(sub)} DEP belong to contamination panels "
                      f"— possible systematic bias (Geyer et al. 2019).")
    return res


def _print_sample_summary(res: dict) -> None:
    info, flags = res["info"], res["flags"]
    print(f"  [PLASMA-QC] Panels: erythrocyte {info['erythrocyte']['n_found']}/"
          f"{info['erythrocyte']['n_entries']}, platelet {info['platelet']['n_found']}/"
          f"{info['platelet']['n_entries']}, coagulation {info['coagulation']['n_found']}/"
          f"{info['coagulation']['n_entries']} quantified "
          f"(fibrinogen chains: {info['coagulation']['n_fibrinogen']}/3)")
    low = [p for p in ("erythrocyte", "platelet")
           if info[p]["n_found"] < max(MIN_PANEL_FOUND, info[p]["n_entries"] // 3)]
    if low:
        print(f"  [WARN] Few panel proteins matched ({', '.join(low)}). Matching uses gene "
              f"symbols (Genes column, from GN= in a UniProt FASTA): a FASTA without gene "
              f"names (NCBI RefSeq, unnamed TrEMBL entries) matches nothing. Check the "
              f"'Panel coverage' table, or give your own list (plasma_qc_panels).")
    if res["n_samples"] <= 5:
        print(f"  [INFO] {res['n_samples']} samples: the mean + 2 SD rule cannot flag any "
              f"sample (largest reachable z = {res['max_z']:.2f}); robust rule only.")
    fl = flags[flags["flagged"]]
    computable = [p for p in PANEL_ORDER if np.isfinite(flags[f"{p}_log2"]).any()]
    if not computable:
        print("  [WARN] No plasma index could be computed: no sample can be assessed.")
        return
    if len(fl):
        for s, r in fl.iterrows():
            print(f"  [WARN] Sample {s} ({r['condition']}): {r['flag_reason']}")
    else:
        print("  [OK] No sample flagged by the plasma quality panels.")
    for p in res["confounded"]:
        row = res["confounding"]
        r = row[(row["index"] == p) & row["comparison"].str.startswith("all")].iloc[0]
        print(f"  [WARN] {INDEX_LABEL[p]} index differs between conditions "
              f"(Kruskal-Wallis p = {r['p_value']:.2g}, eps2 = {r['epsilon2']:.2f}): "
              f"DEP of this panel may be artefacts — see sheet Plasma_QC_DEP.")


def methods_rows(res: dict) -> list:
    """Lignes (clé, valeur, méthode) pour la feuille Methods de Deimos."""
    info = res["info"]
    flagged = res["flags"].index[res["flags"]["flagged"]].tolist()
    rows = [
        ("Marker panels", PANEL_SOURCE if not res.get("panels_path")
         else f"user panels ({os.path.basename(str(res['panels_path']))})",
         "Erythrocyte lysis, platelet contamination, coagulation; quantified: erythrocyte "
         f"{info['erythrocyte']['n_found']}/{info['erythrocyte']['n_entries']}, platelet "
         f"{info['platelet']['n_found']}/{info['platelet']['n_entries']}, coagulation "
         f"{info['coagulation']['n_found']}/{info['coagulation']['n_entries']}"),
        ("Sample indices", "summed intensity ratios (Geyer et al. 2019)",
         "Erythrocyte and platelet: panel / all quantified proteins, without the proteins "
         "shared by the two panels (ACTB, GAPDH), so that hemolysis does not raise the "
         "platelet index (deviation from the article, which keeps them). Coagulation: all "
         "quantified proteins / fibrinogen chains (FGA, FGB, FGG), consumed by the clot. "
         "Higher = more contaminated for all three. Computed on the "
         + res.get("matrix_label", "unfiltered matrix, before imputation") + "."),
        ("Sample flags", f"log2 index > mean + {SD_K:g} SD, or modified z > {ROBUST_Z:g}",
         "Rule of the article (one-sided, SD with ddof = 1), plus a median/MAD rule that the "
         f"contaminated sample cannot mask (Iglewicz & Hoaglin). Flagged: "
         f"{', '.join(flagged) if flagged else 'none'}. Samples are flagged, never removed."),
        ("Confounding", "index vs condition",
         "Kruskal-Wallis across conditions (epsilon^2) and Welch test of the log2 index per "
         "contrast: a contamination confounded with the condition flags no sample but creates "
         "artefactual DEP. "
         + (f"Index differing between conditions (p < {ALPHA}): "
            f"{', '.join(res['confounded'])}." if res["confounded"] else
            f"No index differs between conditions (p >= {ALPHA}).")),
    ]
    if "dep" in res:
        d = res["dep"]
        rows.append(("DEP flags", f"{int((d['flag'] == 'high').sum())} high, "
                     f"{int((d['flag'] == 'check').sum())} check (sheet Plasma_QC_DEP)",
                     "Per DEP and contrast: panel membership; within-condition Pearson "
                     "correlation with each index (condition means removed from both sides; "
                     f"r >= {R_THRESHOLD:g} with BH q < {ALPHA} over all tested proteins, "
                     "one-sided, two-sided for coagulation; the index is recomputed without "
                     "the protein for panel "
                     "members); 'high' when the index also differs between the two groups "
                     "in the same direction."))
    return rows


# ─────────────────────────── écriture Excel ─────────────────────────────────

def _clean(v):
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else round(float(v), 6)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    return v


def _append_table(ws, df: pd.DataFrame, title: "str | None" = None) -> int:
    from openpyxl.styles import Font, PatternFill
    if title:
        ws.append([title])
        ws.cell(ws.max_row, 1).font = Font(bold=True, size=11)
    ws.append([str(c) for c in df.columns])
    hdr = ws.max_row
    for cell in ws[hdr]:
        cell.font = Font(bold=True)
        cell.fill = PatternFill("solid", fgColor="D9D9D9")
    for row in df.itertuples(index=False):
        ws.append([_clean(v) for v in row])
    ws.append([])
    return hdr


def write_sheets(wb, res: dict, ins_img_fn=None) -> list:
    """Feuilles Plasma_QC (échantillons, confusion, couverture des panels, figure) et
    Plasma_QC_DEP (si l'étape 2 a tourné). Retourne les noms écrits."""
    from openpyxl.styles import Font
    written = []
    ws = wb.create_sheet("Plasma_QC")
    ws.append(["Plasma sample quality — " + PANEL_SOURCE])
    ws.cell(1, 1).font = Font(bold=True, size=13)
    ws.append(["Higher index = more contaminated. Samples are flagged, never removed: "
               "decide on exclusion or a covariate before reading the DEP."])
    ws.append([])
    for k, v, m in methods_rows(res):
        ws.append([k, v, m])
        ws.cell(ws.max_row, 1).font = Font(bold=True)
    ws.append([])
    fl = res["flags"].reset_index().rename(columns={"index": "sample"})
    keep = ["sample", "condition", "flagged", "flag_reason"]
    for p in PANEL_ORDER:
        keep += [f"{p}_index", f"{p}_log2", f"{p}_z", f"{p}_robust_z"]
    _append_table(ws, fl[keep], "Samples")
    _append_table(ws, res["confounding"], "Index vs condition")
    _append_table(ws, res["coverage"], "Panel coverage (entries of the panels found in this run)")
    for col, w in zip("ABCD", (26, 22, 14, 60)):
        ws.column_dimensions[col].width = w
    if res.get("plot") and ins_img_fn is not None:
        ins_img_fn(ws, res["plot"], "F5", max_w_px=900, max_h_px=360)
    elif res.get("plot"):
        try:
            from openpyxl.drawing.image import Image as XLImage
            ws.add_image(XLImage(res["plot"]), "F5")
        except Exception:
            pass
    written.append("Plasma_QC")
    if "dep" in res:
        ws2 = wb.create_sheet("Plasma_QC_DEP")
        d = res["dep"]
        order = {"high": 0, "check": 1, "": 2}
        d = d.assign(_o=d["flag"].map(order)).sort_values(["contrast", "_o", "p.value"]) \
             .drop(columns="_o") if len(d) else d
        _append_table(ws2, d)
        ws2.column_dimensions["B"].width = 18
        ws2.column_dimensions["C"].width = 16
        ws2.column_dimensions["K"].width = 90
        written.append("Plasma_QC_DEP")
    return written


# ─────────────────────────── ligne de commande (run terminé) ────────────────

def _read_run(out_dir: str, workbook: "str | None"):
    mdir = os.path.join(out_dir, "matrices")
    mpath = os.path.join(mdir, "log2_filtered_NA.tsv")
    spath = os.path.join(mdir, "samples.tsv")
    if not (os.path.exists(mpath) and os.path.exists(spath)):
        raise FileNotFoundError(f"{mdir}: log2_filtered_NA.tsv and samples.tsv are required "
                                f"(Deimos export_matrices: true)")
    mat = pd.read_csv(mpath, sep="\t")
    smp = pd.read_csv(spath, sep="\t", dtype=str)
    labels = [str(s) for s in smp["label"]]
    missing = [s for s in labels if s not in mat.columns]
    if missing:
        raise ValueError(f"samples.tsv labels absent from the matrix: {missing[:5]}")
    wb_path = workbook or os.path.join(out_dir, "ProteomicAnalysis_Results.xlsx")
    df_res = None
    if os.path.exists(wb_path):
        try:
            df_res = pd.read_excel(wb_path, sheet_name="Differential_Expression")
        except Exception as e:
            print(f"  [WARN] {wb_path}: Differential_Expression not read ({e}) — sample QC only.")
    return mat, smp, labels, df_res


def _volcano_params(out_dir: str) -> dict:
    p = {"volcano_use_padj": False, "volcano_p_thresh": 0.05, "volcano_ratio_min": 1.5}
    cfg = os.path.join(out_dir, "last_config.yaml")
    if os.path.exists(cfg):
        try:
            import yaml
            with open(cfg, encoding="utf-8") as fh:
                raw = yaml.safe_load(fh) or {}
            for k in p:
                if k in raw and raw[k] is not None:
                    p[k] = raw[k]
        except Exception:
            pass
    p["volcano_lfc_min"] = float(np.log2(float(p["volcano_ratio_min"])))
    return p


def sig_masks_from_results(df_res: pd.DataFrame, params: dict) -> tuple:
    """Règle des volcanos de Deimos (volcano_status) : statistique < seuil ET |log2FC| >= seuil."""
    contrasts = [c[:-len("_diff")] for c in df_res.columns
                 if str(c).endswith("_diff") and f"{c[:-len('_diff')]}_p.val" in df_res.columns]
    key = "p.adj" if params["volcano_use_padj"] else "p.val"
    masks = {}
    for c in contrasts:
        lfc = df_res[f"{c}_diff"].to_numpy(float)
        ps = df_res[f"{c}_{key}"].to_numpy(float)
        masks[c] = (ps < params["volcano_p_thresh"]) & (np.abs(lfc) >= params["volcano_lfc_min"])
    return contrasts, masks


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Plasma / serum sample quality (Geyer et al. 2019 panels) on a finished "
                    "Deimos run: per-sample indices, confounding with the condition, DEP flags.")
    ap.add_argument("out_dir", help="Deimos output folder (with matrices/)")
    ap.add_argument("--workbook", default=None,
                    help="Deimos workbook (default: <out_dir>/ProteomicAnalysis_Results.xlsx)")
    ap.add_argument("--panels", default=None,
                    help="custom panels, TSV/CSV with columns gene, panel")
    ap.add_argument("--out", default=None,
                    help="output workbook (default: <out_dir>/plasma_qc/Plasma_QC.xlsx)")
    a = ap.parse_args(argv)
    mat, smp, labels, df_res = _read_run(a.out_dir, a.workbook)
    out_sub = os.path.join(a.out_dir, "plasma_qc")
    os.makedirs(out_sub, exist_ok=True)
    genes = mat["Genes"].fillna("").astype(str).tolist() if "Genes" in mat else [""] * len(mat)
    conds = smp["condition"].astype(str).tolist()
    contrasts, masks = ([], {})
    if df_res is not None:
        contrasts, masks = sig_masks_from_results(df_res, _volcano_params(a.out_dir))
    print("[PLASMA-QC] Sample indices (filtered matrix of the run)...")
    res = run_sample_qc(mat[labels], genes, conds, panels_path=a.panels,
                        contrast_names=contrasts, out_dir=out_sub,
                        matrix_label="filtered matrix of the run (matrices/log2_filtered_NA.tsv)")
    if df_res is not None and contrasts:
        # Lignes du classeur alignées sur la matrice par identifiant de protéine
        key = mat["protein"].astype(str)
        dfr = df_res.copy()
        dfr["_k"] = dfr["name"].astype(str) if "name" in dfr else key.values
        dfr = dfr.set_index("_k").reindex(key.values).reset_index(drop=True)
        dfr["name"] = key.values
        masks = {c: np.nan_to_num(m, nan=0).astype(bool) for c, m in
                 sig_masks_from_results(dfr, _volcano_params(a.out_dir))[1].items()}
        res = run_dep_qc(res, mat[labels], genes, dfr, contrasts, masks, conds)
    from openpyxl import Workbook
    wb = Workbook()
    wb.remove(wb.active)
    write_sheets(wb, res)
    out = a.out or os.path.join(out_sub, "Plasma_QC.xlsx")
    wb.save(out)
    print(f"[DONE] {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
