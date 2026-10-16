# ==============================================================================
# dpc_sensitivity.py — Sensibilité aux valeurs manquantes (en complément)
# Deimos — DIA Expression Integrated Multi-Omics Suite
# ==============================================================================
# Deux analyses AJOUTÉES à la statistique principale (imputation + limma/DEqMS),
# qui reste celle des volcanos, du GO, du score de robustesse et du dashboard :
#
#   1. Analyse différentielle SANS imputation, fondée sur la courbe de
#      probabilité de détection (DPC) :
#        P(détecté | y) = 1 / (1 + exp(-(beta0 + beta1 * y)))      y = log2 intensité
#      - DPC estimée par maximum de vraisemblance sous le modèle « normal
#        complet » (Li & Smyth 2023, Bioinformatics 39:btad200) ;
#      - valeurs protéiques complétées par maximum a posteriori, avec erreur
#        standard (logique de limpa::dpcQuantByRow, Li et al. 2025) ;
#      - limma pondéré : poids de précision vooma avec l'erreur standard comme
#        prédicteur et correction des ddl résiduels quand un groupe entier
#        n'est pas observé (logique de limpa::dpcDE).
#      Implémentation Python indépendante, écrite d'après les méthodes
#      publiées ; validée contre limma (R) sur la partie vooma/eBayes.
#
#   2. Détection différentielle (protéines on/off) : pour une protéine absente
#      (ou presque) d'une condition, probabilité d'observer autant de
#      non-détections si elle avait l'abondance de l'autre condition, d'après
#      la DPC (binomiale intégrée sur l'incertitude de cette abondance).
#      Inspiré de la « differential detection » de MS-DAP (Koopmans 2023),
#      mais avec une probabilité calibrée par la DPC plutôt qu'un z-score.
#
# Rien n'est modifié dans l'analyse principale : ce module ne fait que lire
# la matrice filtrée non imputée, la design matrix et les résultats.
# ==============================================================================

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from numpy.polynomial.hermite_e import hermegauss
from scipy import optimize, stats
from scipy.special import digamma, expit, polygamma

from limma_ebayes import _bh_correction, _inverse_trigamma

# Quadrature de Gauss-Hermite pour une espérance sous N(0, 1)
_GH_X, _GH_W = hermegauss(16)
_GH_W = _GH_W / np.sqrt(2.0 * np.pi)


def _log1m_expit(eta):
    """log(1 - expit(eta)), stable."""
    return -np.logaddexp(0.0, eta)


def _log_expit(eta):
    """log(expit(eta)), stable."""
    return -np.logaddexp(0.0, -eta)


# ==============================================================================
# 1. Complétion par inclinaison exponentielle (valeurs de départ, hyperparamètres)
# ==============================================================================

def exp_tilt_complete(y: np.ndarray, slope: float) -> np.ndarray:
    """Complète les NA par inclinaison exponentielle (Li & Smyth 2023).

    Si P(détection) croît comme exp(slope * y) aux faibles intensités, une valeur
    manquante suit approximativement N(mu_obs - slope * sigma², sigma²). On
    combine l'estimation par ligne (protéine) et par colonne (échantillon),
    pondérées par l'inverse de leur variance de prédiction. Sert uniquement de
    point de départ et à fixer les hyperparamètres.
    """
    y = np.asarray(y, dtype=float)
    out = y.copy()
    miss = np.isnan(y)
    if not miss.any():
        return out
    nobs_r = (~miss).sum(axis=1)
    nobs_c = (~miss).sum(axis=0)

    # Lignes
    m_r = np.nanmean(np.where(nobs_r[:, None] > 0, y, 0.0), axis=1)
    m_r[nobs_r == 0] = np.nan
    with np.errstate(invalid="ignore", divide="ignore"):
        s2_r = np.nansum((y - m_r[:, None]) ** 2, axis=1) / (nobs_r - 1)
    ok = nobs_r > 1
    s2_row = float(np.quantile(s2_r[ok], 0.9)) if ok.any() else np.nan
    df_row = max(int(nobs_r.max()) - 1, 1)
    mumis_r = m_r - slope * s2_row
    s2_row_f = max(s2_row, 1e-3) if np.isfinite(s2_row) else np.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        pv_r = (nobs_r + 1) / nobs_r * s2_row_f + 2.0 / df_row * (slope * s2_row_f) ** 2

    # Colonnes
    m_c = np.nanmean(y, axis=0)
    dfc = np.maximum(nobs_c - 1, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        s2_c = np.nansum((y - m_c[None, :]) ** 2, axis=0) / dfc
    s2_c = np.where(np.isfinite(s2_c), np.maximum(s2_c, 0.1), np.nan)
    small = dfc < 5
    if small.any() and np.isfinite(s2_c).any():
        s2_c[small] = np.nanquantile(s2_c, 0.75)
        dfc = np.where(small, dfc.max(), dfc)
    mumis_c = m_c - slope * s2_c
    with np.errstate(divide="ignore", invalid="ignore"):
        pv_c = (nobs_c + 1) / nobs_c * s2_c + 2.0 / np.maximum(dfc, 1) * (slope * s2_c) ** 2

    # Repli global si une ligne/colonne n'a pas d'estimation
    g_obs = y[~miss]
    g_m, g_s2 = g_obs.mean(), g_obs.var(ddof=1) if g_obs.size > 1 else 1.0
    g_mu = g_m - slope * g_s2
    g_pv = (g_obs.size + 1) / g_obs.size * g_s2 + 2.0 / max(g_obs.size - 1, 1) * (slope * g_s2) ** 2
    bad_r = ~np.isfinite(mumis_r) | ~np.isfinite(pv_r)
    mumis_r[bad_r], pv_r[bad_r] = g_mu, g_pv
    bad_c = ~np.isfinite(mumis_c) | ~np.isfinite(pv_c)
    mumis_c[bad_c], pv_c[bad_c] = g_mu, g_pv

    i, j = np.where(miss)
    w_r, w_c = 1.0 / pv_r[i], 1.0 / pv_c[j]
    out[i, j] = (mumis_r[i] * w_r + mumis_c[j] * w_c) / (w_r + w_c)
    return out


def estimate_dpc_intercept(y: np.ndarray, slope: float) -> float:
    """Ordonnée à l'origine de la DPC pour une pente fixée : celle qui redonne la
    proportion globale de valeurs observées sur la matrice complétée."""
    y = np.asarray(y, dtype=float)
    obs = ~np.isnan(y)
    n_obs = obs.sum()
    if abs(slope) < 1e-8:
        p = n_obs / obs.size
        return float(np.log(p / (1 - p)))
    yc = exp_tilt_complete(y, slope).ravel(order="F")
    isobs = obs.ravel(order="F").astype(float)
    if yc.size > 10000:
        # Agrégation en 300 intervalles d'intensité de même largeur (comme limpa) :
        # même estimation, calcul plus léger et identique à la référence R.
        lo, hi = yc.min(), yc.max()
        dx = hi - lo
        edges = np.linspace(lo, hi, 301)
        edges[0], edges[-1] = lo - dx / 1000, hi + dx / 1000
        b = np.clip(np.searchsorted(edges, yc, side="left") - 1, 0, 299)
        N = np.bincount(b, minlength=300).astype(float)
        ok = N > 0
        mean_y = np.bincount(b, weights=yc, minlength=300)[ok] / N[ok]
        n_obs_b = np.bincount(b, weights=isobs, minlength=300)[ok]
        N = N[ok]
    else:
        mean_y, n_obs_b, N = yc, isobs, np.ones_like(yc)
    target = float(n_obs_b.sum())

    def f(b0):
        return (N * expit(b0 + slope * mean_y)).sum() - target

    lo, hi = -200.0, 200.0
    return float(optimize.brentq(f, lo, hi))


# ==============================================================================
# 2. Estimation de la DPC — modèle normal complet (Li & Smyth 2023)
# ==============================================================================

def _row_map_mu_sigma(yi, b0, b1, mu0, sig0, prior):
    """MAP de (mu, log sigma) pour une ligne, DPC fixée (y complet ~ N(mu, sigma²))."""
    obs = ~np.isnan(yi)
    yo = yi[obs]
    nobs, nmis = obs.sum(), (~obs).sum()
    mm, msd, lm, lsd = prior

    def fun(par):
        mu, ls = par
        s = np.exp(ls)
        q = ((mu - mm) / msd) ** 2 + ((ls - lm) / lsd) ** 2
        g_mu = 2 * (mu - mm) / msd ** 2
        g_ls = 2 * (ls - lm) / lsd ** 2
        if nobs:
            r = yo - mu
            q += 2 * (nobs - 1) * ls + np.sum(r ** 2) / s ** 2
            g_mu += -2 * np.sum(r) / s ** 2
            g_ls += 2 * (nobs - 1) - 2 * np.sum(r ** 2) / s ** 2
        if nmis:
            eta = b0 + b1 * (mu + s * _GH_X)
            pm = expit(-eta)                       # P(manquant | y)
            G = np.sum(_GH_W * pm)
            dl = expit(eta) * pm                   # densité logistique
            q -= 2 * nmis * np.log(G)
            g_mu += 2 * nmis * b1 * np.sum(_GH_W * dl) / G
            g_ls += 2 * nmis * s * b1 * np.sum(_GH_W * _GH_X * dl) / G
        return q, np.array([g_mu, g_ls])

    res = optimize.minimize(fun, x0=[mu0, np.log(max(sig0, 1e-3))], jac=True,
                            method="L-BFGS-B")
    return res.x[0], np.exp(res.x[1])


def estimate_dpc(y: np.ndarray, slope_start: float = 0.7, iterations: int = 2,
                 subset: int = 2000) -> dict:
    """Estime (beta0, beta1) de la DPC par maximum de vraisemblance sous le
    modèle normal complet : y_ij ~ N(mu_i, sigma_i²) avant détection, et
    P(détecté | y) logistique. Alterne MAP de (mu_i, sigma_i) par ligne et MLE
    de la DPC (quadrature de Gauss-Hermite pour les valeurs manquantes).

    À appliquer de préférence sur la matrice log2 NON FILTRÉE : le filtrage sur
    le nombre de détections sélectionne les lignes bien détectées et biaise la
    DPC vers le haut."""
    y = np.asarray(y, dtype=float)
    y = y[(~np.isnan(y)).any(axis=1)]
    n_rows, n_samples = y.shape
    if n_rows < 3 or n_samples < 2:
        raise ValueError("too few rows or samples to estimate the DPC")
    if not np.isnan(y).any():
        raise ValueError("no missing value: the DPC is not identifiable")

    if n_rows > subset:                     # sous-échantillon systématique
        nmis = np.isnan(y).sum(axis=1)
        ave = np.nanmean(y, axis=1)
        o = np.lexsort((ave, nmis))
        idx = np.round(((np.arange(1, subset + 1) - 0.5) * n_rows / subset)).astype(int) - 1
        y = y[o[np.clip(idx, 0, n_rows - 1)]]
        n_rows = subset

    b1 = float(slope_start)
    b0 = estimate_dpc_intercept(y, b1)
    yc = exp_tilt_complete(y, b1)
    mu = yc.mean(axis=1)
    sig = np.sqrt(((yc - mu[:, None]) ** 2).sum(axis=1) / (n_samples - 1))
    logsig = np.log(sig + 0.001)
    prior = (mu.mean(), max(mu.std(ddof=1), 1e-3), logsig.mean(), max(logsig.std(ddof=1), 1e-3))

    obs = ~np.isnan(y)
    nmis = (~obs).sum(axis=1)
    y0 = np.where(obs, y, 0.0)
    history = [(b0, b1)]

    for _ in range(int(iterations)):
        for i in range(n_rows):
            mu[i], sig[i] = _row_map_mu_sigma(y[i], b0, b1, mu[i], sig[i], prior)

        z = mu[:, None] + sig[:, None] * _GH_X[None, :]          # (n_rows × 16)

        def negll(par):
            c0, t = par
            c1 = np.exp(t / 10.0)
            eta_o = c0 + c1 * y0
            q = np.sum(np.where(obs, _log_expit(eta_o), 0.0))
            d0 = np.sum(np.where(obs, expit(-eta_o), 0.0))
            d1 = np.sum(np.where(obs, y0 * expit(-eta_o), 0.0))
            eta_m = c0 + c1 * z
            pm = expit(-eta_m)
            G = pm @ _GH_W
            dl = expit(eta_m) * pm
            has = nmis > 0
            q += np.sum(nmis[has] * np.log(G[has]))
            d0 -= np.sum(nmis[has] * (dl[has] @ _GH_W) / G[has])
            d1 -= np.sum(nmis[has] * ((dl[has] * z[has]) @ _GH_W) / G[has])
            return -2 * q, -2 * np.array([d0, d1 * c1 / 10.0])

        res = optimize.minimize(negll, x0=[b0, 10.0 * np.log(b1)], jac=True, method="BFGS")
        b0, b1 = float(res.x[0]), float(np.exp(res.x[1] / 10.0))
        history.append((b0, b1))

    return {"beta0": b0, "beta1": b1, "n_rows": int(n_rows), "n_samples": int(n_samples),
            "history": history, "prop_observed": float(obs.mean())}


def detection_probability(y, beta0: float, beta1: float):
    return expit(beta0 + beta1 * np.asarray(y, dtype=float))


# ==============================================================================
# 3. Complétion MAP par protéine, avec erreurs standard (logique dpcQuantByRow)
# ==============================================================================

def dpc_complete(y: np.ndarray, beta0: float, beta1: float,
                 sigma_obs: float = 0.01, sd_quantile: float = 0.9,
                 max_iter: int = 100, tol: float = 1e-8) -> dict:
    """Valeurs protéiques complétées par maximum a posteriori et leur erreur
    standard, sans connaître le design (a priori échangeable entre échantillons).

    Modèle, pour la protéine i :
      m_i ~ N(mu0, s0²)                       (abondance moyenne)
      z_ij ~ N(m_i, tau²)                     (valeur dans l'échantillon j)
      y_ij | z_ij ~ N(z_ij, sigma_obs²)       si observée
      P(détectée | z_ij) = expit(beta0 + beta1 z_ij)
    mu0, s0 et tau viennent de la matrice complétée par inclinaison
    exponentielle. La log-postérieure est concave : Newton sur (m_i, z_i.), dont
    la hessienne en « flèche » s'inverse analytiquement (vectorisé)."""
    y = np.asarray(y, dtype=float)
    obs = ~np.isnan(y)
    n_rows, n = y.shape
    yc = exp_tilt_complete(y, beta1)
    rm = yc.mean(axis=1)
    mu0 = float(rm.mean())
    s0 = float(max(rm.std(ddof=1), 1e-3))
    rv = ((yc - rm[:, None]) ** 2).sum(axis=1) / max(n - 1, 1)
    tau = float(np.sqrt(max(np.quantile(rv, sd_quantile), 1e-6)))

    a_prior, it2, ie2 = 1.0 / s0 ** 2, 1.0 / tau ** 2, 1.0 / sigma_obs ** 2
    z = yc.copy()
    m = rm.copy()
    y0 = np.where(obs, y, 0.0)

    for _ in range(max_iter):
        p = expit(beta0 + beta1 * z)
        g_m = (m - mu0) * a_prior - ((z - m[:, None]) * it2).sum(axis=1)
        g_z = (z - m[:, None]) * it2 + np.where(obs, (z - y0) * ie2 - beta1 * (1 - p), beta1 * p)
        D = it2 + np.where(obs, ie2, 0.0) + beta1 ** 2 * p * (1 - p)
        c = -it2
        a = a_prior + n * it2
        S = a - c ** 2 * (1.0 / D).sum(axis=1)
        dm = np.clip((-g_m + c * (g_z / D).sum(axis=1)) / S, -2.0, 2.0)
        dz = np.clip((-g_z - c * dm[:, None]) / D, -2.0, 2.0)
        m += dm
        z += dz
        if max(np.abs(dm).max(), np.abs(dz).max()) < tol:
            break

    p = expit(beta0 + beta1 * z)
    D = it2 + np.where(obs, ie2, 0.0) + beta1 ** 2 * p * (1 - p)
    S = (a_prior + n * it2) - it2 ** 2 * (1.0 / D).sum(axis=1)
    var_z = 1.0 / D + it2 ** 2 / (D ** 2 * S[:, None])
    return {"E": z, "SE": np.sqrt(var_z), "imputed": ~obs, "row_mean": m,
            "prior_mean": mu0, "prior_sd": s0, "prior_logfc": tau}


# ==============================================================================
# 3bis. Quantification protéique À PARTIR DES PRÉCURSEURS (logique dpcQuant)
# ==============================================================================

def precursor_hyperparams(Y, pid, slope, sd_quantile=0.9):
    """Hyperparamètres de la quantification à partir des précurseurs : a priori
    sur l'abondance protéique (moyenne, écart-type, variation entre échantillons)
    et écart-type résiduel du modèle additif précurseur + échantillon, modéré par
    protéine (eBayes)."""
    Y = np.asarray(Y, float)
    n = Y.shape[1]
    yc = exp_tilt_complete(Y, slope)
    rm = yc.mean(axis=1)
    rv = ((yc - rm[:, None]) ** 2).sum(axis=1) / max(n - 1, 1)
    P = int(pid.max()) + 1
    K = np.bincount(pid, minlength=P).astype(float)
    y0 = yc - rm[:, None]
    cm = np.zeros((P, n))
    np.add.at(cm, pid, y0)
    cm /= K[:, None]
    res = y0 - cm[pid]
    ss = np.zeros(P)
    np.add.at(ss, pid, (res ** 2).sum(axis=1))
    df = (n - 1) * (K - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        s2 = ss / df
    d0, s0 = fit_fdist(s2, df)
    if np.isinf(d0):
        s2_post = np.full(P, s0)
    else:
        s2_post = np.where(df > 0, (df * np.nan_to_num(s2) + d0 * s0) / (df + d0), s0)
    return {"prior_mean": float(rm.mean()),
            "prior_sd": float(max(rm.std(ddof=1), 1.0)),
            "prior_logfc": float(max(np.sqrt(np.quantile(rv, sd_quantile)), 1.0)),
            "sigma": np.sqrt(s2_post), "start": yc}


def _quant_one_protein(Yp, ystart, sigma, b0, b1, m0, s0, tau, max_iter=60, tol=1e-7):
    """MAP des abondances protéiques (une par échantillon) et des effets
    précurseurs (somme nulle) pour UNE protéine ; erreurs standard de Laplace.

    Modèle : y_kj ~ N(z_j + b_k, sigma²) avant détection ;
    P(détecté) = expit(beta0 + beta1 * y) (non-détection intégrée sur sigma) ;
    a priori : moyenne(z) ~ N(m0, s0²), z_j - moyenne(z) ~ N(0, tau²)."""
    K, n = Yp.shape
    obs = ~np.isnan(Yp)
    mis = ~obs
    y0 = np.where(obs, Yp, 0.0)
    s2 = sigma ** 2
    z = ystart.mean(axis=0)
    bfull = ystart.mean(axis=1) - ystart.mean()
    b = bfull[:-1].copy() if K > 1 else np.zeros(0)
    has_mis = mis.any()
    npar = n + K - 1
    J = np.full((n, n), 1.0 / n)
    Hp = (np.eye(n) - J) / tau ** 2 + J / (n * s0 ** 2)

    def unpack(b):
        return np.append(b, -b.sum()) if K > 1 else np.zeros(1)

    def objective(z, b):
        mu = z[None, :] + unpack(b)[:, None]
        f = np.sum(np.where(obs, (y0 - mu) ** 2, 0.0)) / (2 * s2)
        if has_mis:
            eta = b0 + b1 * (mu[mis][:, None] + sigma * _GH_X[None, :])
            G = expit(-eta) @ _GH_W
            f -= np.sum(np.log(np.maximum(G, 1e-300)))
        zb = z.mean()
        f += (zb - m0) ** 2 / (2 * s0 ** 2) + np.sum((z - zb) ** 2) / (2 * tau ** 2)
        return f

    def derivs(z, b):
        mu = z[None, :] + unpack(b)[:, None]
        g = np.where(obs, -(y0 - mu) / s2, 0.0)
        w = np.where(obs, 1.0 / s2, 0.0)
        if has_mis:
            eta = b0 + b1 * (mu[mis][:, None] + sigma * _GH_X[None, :])
            p = expit(eta)
            pm = 1.0 - p
            G = np.maximum(pm @ _GH_W, 1e-300)
            G1 = -b1 * ((p * pm) @ _GH_W)
            G2 = -b1 ** 2 * ((p * pm * (1 - 2 * p)) @ _GH_W)
            g[mis] = -G1 / G
            w[mis] = (G1 ** 2 - G * G2) / G ** 2
        zb = z.mean()
        gz = g.sum(axis=0) + (zb - m0) / (n * s0 ** 2) + (z - zb) / tau ** 2
        H = np.zeros((npar, npar))
        H[:n, :n] = np.diag(w.sum(axis=0)) + Hp
        if K > 1:
            gb = g[:-1].sum(axis=1) - g[-1].sum()
            Hzb = (w[:-1] - w[-1][None, :]).T                        # (n × K-1)
            H[:n, n:] = Hzb
            H[n:, :n] = Hzb.T
            H[n:, n:] = np.diag(w[:-1].sum(axis=1)) + w[-1].sum()
            grad = np.concatenate([gz, gb])
        else:
            grad = gz
        return grad, H

    f = objective(z, b)
    for _ in range(max_iter):
        grad, H = derivs(z, b)
        try:
            step = -np.linalg.solve(H, grad)
        except np.linalg.LinAlgError:
            step = -np.linalg.lstsq(H, grad, rcond=None)[0]
        t = 1.0
        while True:
            zn, bn = z + t * step[:n], b + t * step[n:]
            fn = objective(zn, bn)
            if fn <= f + 1e-4 * t * grad @ step or t < 1e-6:
                break
            t *= 0.5
        z, b, f_old, f = zn, bn, f, fn
        if np.max(np.abs(t * step)) < tol or abs(f_old - f) < 1e-12:
            break
    grad, H = derivs(z, b)
    try:
        cov = np.linalg.inv(H)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(H)
    return z, np.sqrt(np.maximum(np.diag(cov)[:n], 0.0))


def dpc_quant_precursors(Y, protein_ids, beta0, beta1, verbose=False):
    """Abondances protéiques et erreurs standard à partir de la matrice des
    précurseurs (NaN = non détecté), sous la DPC. Retourne un dict avec
    'proteins' (ordre des protéines), 'E', 'SE', 'n_obs' (précurseurs observés
    par protéine × échantillon), 'n_prec'."""
    Y = np.asarray(Y, float)
    keep = (~np.isnan(Y)).any(axis=1)
    Y = Y[keep]
    ids = np.asarray(protein_ids)[keep]
    proteins, pid = np.unique(ids, return_inverse=True)
    h = precursor_hyperparams(Y, pid, beta1)
    order = np.argsort(pid, kind="mergesort")
    bounds = np.searchsorted(pid[order], np.arange(len(proteins) + 1))
    n = Y.shape[1]
    E = np.zeros((len(proteins), n))
    SE = np.zeros_like(E)
    nobs = np.zeros_like(E, dtype=int)
    for p in range(len(proteins)):
        rows = order[bounds[p]:bounds[p + 1]]
        Yp = Y[rows]
        E[p], SE[p] = _quant_one_protein(Yp, h["start"][rows], h["sigma"][p], beta0, beta1,
                                         h["prior_mean"], h["prior_sd"], h["prior_logfc"])
        nobs[p] = (~np.isnan(Yp)).sum(axis=0)
        if verbose and (p + 1) % 1000 == 0:
            print(f"    {p + 1}/{len(proteins)} proteins quantified", flush=True)
    return {"proteins": proteins, "E": E, "SE": SE, "n_obs": nobs,
            "n_prec": np.bincount(pid, minlength=len(proteins)), "hyper": h}


# ==============================================================================
# 4. limma pondéré : vooma avec prédicteur + correction des ddl (logique dpcDE)
# ==============================================================================

def choose_lowess_span(n: int, small_n: int = 50, min_span: float = 0.3,
                       power: float = 1 / 3) -> float:
    """Équivalent de limma::chooseLowessSpan."""
    return float(min(min_span + (1 - min_span) * (small_n / n) ** power, 1.0))


def _lowess_r(x, y, f, iterations=3, delta=None):
    """Lowess de Cleveland tel que stats::lowess de R (fenêtre de round(f n)
    points, noyau tricube, itérations robustes bisquare, raccourci « delta »).
    x, y triés par x croissant. Retourne les valeurs lissées."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    if n < 2:
        return y.copy()
    if delta is None:
        delta = 0.01 * (x[-1] - x[0])
    ns = max(min(int(f * n + 1e-7), n), 2)
    xrange_ = x[-1] - x[0]
    ys = np.zeros(n)
    rw = np.ones(n)

    def lowest(xs, nleft, nright, userw):
        h = max(xs - x[nleft], x[nright] - xs)
        h9, h1 = 0.999 * h, 0.001 * h
        # points j >= nleft dans |x - xs| <= h9, en s'arrêtant au premier x > xs hors fenêtre
        jj = np.arange(nleft, n)
        r = np.abs(x[jj] - xs)
        stop = np.where((r > h9) & (x[jj] > xs))[0]
        if len(stop):
            jj, r = jj[: stop[0]], r[: stop[0]]
        inside = r <= h9
        w = np.where(r <= h1, 1.0, (1 - (r / h) ** 3) ** 3) if h > 0 else np.ones_like(r)
        w = np.where(inside, w, 0.0)
        if userw:
            w = w * rw[jj]
        a = w.sum()
        if a <= 0:
            return None
        w = w / a
        if h > 0:
            a = (w * x[jj]).sum()
            b = xs - a
            c = (w * (x[jj] - a) ** 2).sum()
            if np.sqrt(c) > 0.001 * xrange_:
                b /= c
                w = w * (b * (x[jj] - a) + 1)
        return (w * y[jj]).sum()

    for it in range(1, iterations + 2):
        nleft, nright, last, i = 0, ns - 1, -1, 0
        while True:
            if nright < n - 1:
                d1, d2 = x[i] - x[nleft], x[nright + 1] - x[i]
                if d1 > d2:
                    nleft += 1
                    nright += 1
                    continue
            v = lowest(x[i], nleft, nright, it > 1)
            ys[i] = y[i] if v is None else v
            if last < i - 1:
                denom = x[i] - x[last]
                for jx in range(last + 1, i):
                    alpha = (x[jx] - x[last]) / denom
                    ys[jx] = alpha * ys[i] + (1 - alpha) * ys[last]
            last = i
            cut = x[last] + delta
            i = last + 1
            while i < n:
                if x[i] > cut:
                    break
                if x[i] == x[last]:
                    ys[i] = ys[last]
                    last = i
                i += 1
            i = max(last + 1, i - 1)
            if last >= n - 1:
                break
        res = y - ys
        sc = np.abs(res).sum() / n
        if it > iterations:
            break
        ar = np.abs(res)
        m1 = n // 2
        srt = np.sort(ar)
        cmad = 3.0 * (srt[m1] + srt[n - m1 - 1]) if n % 2 == 0 else 6.0 * srt[m1]
        if cmad < 1e-7 * sc:
            break
        c9, c1 = 0.999 * cmad, 0.001 * cmad
        rw = np.where(ar <= c1, 1.0, np.where(ar <= c9, (1 - (ar / cmad) ** 2) ** 2, 0.0))
    return ys


def _weighted_lowess_limma(x, y, w, span, iterations=4, npts=200):
    """limma::weightedLowess : la fenêtre contient une fraction 'span' du POIDS
    total (et non du nombre de points). x, y, w triés par x croissant."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    w = np.asarray(w, float)
    n = len(x)
    if npts >= n:
        delta = 0.0
    else:
        dx = np.sort(np.diff(x))
        cum = np.cumsum(dx)
        k = np.arange(npts)
        delta = float(np.min(cum[len(dx) - 1 - k] / (npts - k)))
    seeds = [0]
    lastp = 0
    for p in range(1, n - 1):
        if x[p] - x[lastp] > delta:
            seeds.append(p)
            lastp = p
    seeds.append(n - 1)
    total = w.sum()
    spanweight = total * span
    limits = []
    for cp in seeds:
        left = right = cp
        curw = w[cp]
        ende, ends = cp == n - 1, cp == 0
        mdist = 0.0
        while curw < spanweight and (not ende or not ends):
            if ende:
                left -= 1
                curw += w[left]
                ends = left == 0
                mdist = max(mdist, x[cp] - x[left])
            elif ends:
                right += 1
                curw += w[right]
                ende = right == n - 1
                mdist = max(mdist, x[right] - x[cp])
            else:
                ld, rd = x[cp] - x[left - 1], x[right + 1] - x[cp]
                if ld < rd:
                    left -= 1
                    curw += w[left]
                    ends = left == 0
                    mdist = max(mdist, ld)
                else:
                    right += 1
                    curw += w[right]
                    ende = right == n - 1
                    mdist = max(mdist, rd)
        while left > 0 and x[left] == x[left - 1]:
            left -= 1
        while right < n - 1 and x[right] == x[right + 1]:
            right += 1
        limits.append((left, right, mdist))
    subrange = (x[-1] - x[0]) / n
    rob = np.ones(n)
    fitted = np.zeros(n)
    TH = 1e-7

    def fit_at(cp, left, right, dist):
        sl = slice(left, right + 1)
        if dist < TH:
            ww = w[sl] * rob[sl]
            return (ww * y[sl]).sum() / ww.sum()
        ww = (1 - (np.abs(x[cp] - x[sl]) / dist) ** 3) ** 3 * w[sl] * rob[sl]
        aw = ww.sum()
        xm, ym = (ww * x[sl]).sum() / aw, (ww * y[sl]).sum() / aw
        var = (ww * (x[sl] - xm) ** 2).sum()
        if var < TH:
            return ym
        slope = (ww * (x[sl] - xm) * (y[sl] - ym)).sum() / var
        return slope * x[cp] + ym - slope * xm

    for _ in range(iterations):
        lastp = 0
        fitted[0] = fit_at(0, *limits[0])
        for s_i in range(1, len(seeds)):
            p = seeds[s_i]
            fitted[p] = fit_at(p, *limits[s_i])
            if p - lastp > 1:
                cur = x[p] - x[lastp]
                if cur > TH * subrange:
                    slope = (fitted[p] - fitted[lastp]) / cur
                    icpt = fitted[p] - slope * x[p]
                    fitted[lastp + 1:p] = slope * x[lastp + 1:p] + icpt
                else:
                    fitted[lastp + 1:p] = 0.5 * (fitted[p] + fitted[lastp])
            lastp = p
        ar = np.abs(y - fitted)
        resid_scale = ar.mean()
        o = np.argsort(ar, kind="mergesort")
        cw = np.cumsum(w[o])
        half = total / 2
        k = int(np.searchsorted(cw, half, side="left"))
        if k < n and cw[k] == half and k + 1 < n:
            cmad = 3 * (ar[o[k]] + ar[o[k + 1]])
        else:
            k = int(np.searchsorted(cw, half, side="right"))
            cmad = 6 * ar[o[min(k, n - 1)]]
        if cmad <= TH * resid_scale:
            break
        rob = np.where(ar < cmad, (1 - (ar / cmad) ** 2) ** 2, 0.0)
    return fitted


def _variance_trend(sx, sy, span, weights=None):
    """Tendance sqrt(sigma) ~ sx comme voomaLmFitWithImputation : lowess de R sans
    poids, limma::loessFit (weightedLowess, poids = ddl résiduels) avec poids.
    Retourne (x triés, y lissés) avec les ex-aequo moyennés (approxfun)."""
    o = np.argsort(sx, kind="mergesort")
    xs, ys_ = sx[o], sy[o]
    if weights is None:
        fit = _lowess_r(xs, ys_, span, iterations=3)
    else:
        w = np.clip(np.where(np.isfinite(weights[o]), weights[o], 0.0), 1e-5, 1e5)
        if w.max() - w.min() < 1e-15:
            fit = _lowess_r(xs, ys_, span, iterations=3)
        else:
            fit = _weighted_lowess_limma(xs, ys_, w, span, iterations=4, npts=200)
    ux, inv = np.unique(xs, return_inverse=True)
    uy = np.bincount(inv, weights=fit) / np.bincount(inv)
    return ux, uy


def _weighted_rowwise_fit(E, W, X):
    """Moindres carrés pondérés ligne par ligne. W = 0 exclut une cellule.
    Retourne coefficients, matrice (X'WX)^-1 par ligne, sigma, df résiduels."""
    n_rows, n = E.shape
    p = X.shape[1]
    XtWX = np.einsum("sp,ns,sq->npq", X, W, X)
    cov = np.linalg.pinv(XtWX)
    XtWy = np.einsum("sp,ns,ns->np", X, W, np.nan_to_num(E))
    coef = np.einsum("npq,nq->np", cov, XtWy)
    fitted = coef @ X.T
    res = np.where(W > 0, np.nan_to_num(E) - fitted, 0.0)
    # rang effectif par ligne (cellules de poids > 0)
    rank = np.array([np.linalg.matrix_rank(X[W[i] > 0]) if (W[i] > 0).any() else 0
                     for i in range(n_rows)]) if (W == 0).any() else \
        np.full(n_rows, np.linalg.matrix_rank(X))
    df = (W > 0).sum(axis=1) - rank
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma = np.sqrt((W * res ** 2).sum(axis=1) / df)
    sigma[df <= 0] = np.nan
    return coef, cov, sigma, df.astype(float), fitted


def fit_fdist(s2, df):
    """Hyperparamètres (d0, s0²) de la loi a priori des variances, ddl par ligne
    (limma::fitFDist sans tendance)."""
    s2 = np.asarray(s2, float)
    df = np.asarray(df, float)
    ok = np.isfinite(s2) & (s2 > 0) & (df > 0)
    if ok.sum() < 3:
        return np.inf, float(np.nanmedian(s2[ok])) if ok.any() else 1.0
    z = np.log(s2[ok])
    d = df[ok]
    e = z - digamma(d / 2) + np.log(d / 2)
    emean = e.mean()
    evar = ((e - emean) ** 2).sum() / (len(e) - 1) - polygamma(1, d / 2).mean()
    if evar > 0:
        d0 = 2.0 * _inverse_trigamma(evar)
        s0 = float(np.exp(emean + digamma(d0 / 2) - np.log(d0 / 2)))
    else:
        d0, s0 = np.inf, float(np.exp(emean))
    return float(d0), s0


def dpc_de(E, SE, imputed, design, contrast_mat, fdr_global=False, eps=1e-6):
    """Analyse différentielle sur les valeurs complétées : poids vooma prédits par
    l'abondance moyenne ET l'erreur standard de quantification, ddl résiduels
    corrigés pour les groupes entièrement non observés, puis eBayes."""
    E = np.asarray(E, float)
    n_rows, n = E.shape
    X = np.asarray(design, float)
    C = np.asarray(contrast_mat, float)
    P = np.log(np.asarray(SE, float) + eps)
    imp = np.asarray(imputed, bool)
    if imp.sum(axis=1).max() <= 1:
        imp = np.zeros_like(imp)

    A = E.mean(axis=1)
    ones = np.ones_like(E)
    coef0, cov0, sigma, df, fitted = _weighted_rowwise_fit(E, ones, X)

    # Cellules imputées dont la valeur ajustée ne dépend d'aucune valeur observée
    # (ex. groupe entièrement non observé) : elles ne portent pas de ddl résiduel.
    H = X @ np.linalg.pinv(X)                                  # matrice chapeau
    dep_obs = (np.abs(H) > 1e-8).astype(float) @ (~imp).T.astype(float)   # (n × n_rows)
    all_imp = imp & (dep_obs.T < 0.5)
    has_loss = all_imp.any(axis=1)
    if has_loss.any():
        Wna = np.where(all_imp, 0.0, 1.0)
        _, _, sig_na, df_na, _ = _weighted_rowwise_fit(E[has_loss], Wna[has_loss], X)
        no_df = ~(df_na > 0.5)
        sig_na[no_df] = np.nanmedian(sigma)
        sigma[has_loss], df[has_loss] = sig_na, df_na
        A_na = np.nanmean(np.where(all_imp, np.nan, E), axis=1)
        A[has_loss] = np.where(np.isfinite(A_na[has_loss]), A_na[has_loss], np.nanmedian(A))

    # Tendance de variance : sqrt(sigma) ~ abondance + log(erreur standard)
    sxc = P.mean(axis=1)
    if has_loss.any():
        pna = np.nanmean(np.where(all_imp, np.nan, P), axis=1)
        sxc[has_loss] = np.where(np.isfinite(pna[has_loss]), pna[has_loss], np.nanmedian(sxc))
    sy = np.sqrt(sigma)
    okt = np.isfinite(sy)
    M = np.column_stack([np.ones(n_rows), A, sxc])
    beta, *_ = np.linalg.lstsq(M[okt], sy[okt], rcond=None)
    sx = M @ beta
    mu_cells = beta[0] + beta[1] * fitted + beta[2] * P
    span = choose_lowess_span(n_rows)
    lx, ly = _variance_trend(sx[okt], sy[okt], span,
                             weights=df[okt] if has_loss.any() else None)
    ly = np.maximum(ly, 1e-6)
    W = 1.0 / np.interp(mu_cells, lx, ly) ** 4

    # Ajustement pondéré final ; sigma/ddl recalculés sans les cellules « tout imputé »
    coef, cov, sigma_w, df_w, _ = _weighted_rowwise_fit(E, W, X)
    if has_loss.any():
        W_na = np.where(all_imp, 0.0, W)
        _, _, sig_na, df_na, _ = _weighted_rowwise_fit(E[has_loss], W_na[has_loss], X)
        sig_na[~(df_na > 0.5)] = np.nanmedian(sigma_w)
        sigma_w[has_loss], df_w[has_loss] = sig_na, df_na

    # Contrastes (covariance exacte par ligne)
    coef_c = coef @ C
    su_c = np.sqrt(np.einsum("pk,npq,qk->nk", C, cov, C))

    # eBayes à ddl variables
    s2 = sigma_w ** 2
    d0, s0 = fit_fdist(s2, df_w)
    if np.isinf(d0):
        s2_post = np.full(n_rows, s0)
        df_total = np.full(n_rows, df_w[df_w > 0].sum())
    else:
        dfr = np.where(df_w > 0, df_w, 0.0)
        s2r = np.where(df_w > 0, s2, 0.0)
        s2_post = (dfr * s2r + d0 * s0) / (dfr + d0)
        df_total = np.minimum(dfr + d0, dfr.sum())
    t = coef_c / (su_c * np.sqrt(s2_post)[:, None])
    pval = 2 * stats.t.sf(np.abs(t), df=df_total[:, None])
    if fdr_global:
        padj = _bh_correction(pval.ravel()).reshape(pval.shape)
    else:
        padj = np.column_stack([_bh_correction(pval[:, k]) for k in range(pval.shape[1])])
    return {"logFC": coef_c, "t": t, "p_value": pval, "p_adj": padj,
            "df_prior": d0, "s2_prior": s0, "df_total": df_total, "weights": W,
            "sigma": sigma_w, "df_residual": df_w, "df_corrected": has_loss,
            "n_df_corrected": int(has_loss.sum()), "span": span}


# ==============================================================================
# 5. Détection différentielle (protéines on/off)
# ==============================================================================

def _prob_missing(mu, sigma, b0, b1):
    """P(non détectée) pour une valeur ~ N(mu, sigma²) (quadrature)."""
    mu = np.asarray(mu, float)[..., None]
    s = np.asarray(sigma, float)[..., None]
    return expit(-(b0 + b1 * (mu + s * _GH_X))) @ _GH_W


def _group_abundance_map(yo_sum, yo_n, nmis, sigma, b0, b1, mu_start, iters=50):
    """Abondance d'un groupe (normale, sigma connu) estimée par maximum de
    vraisemblance à partir de ses valeurs observées ET de ses non-détections
    (DPC). Vectorisé ; log-vraisemblance concave. Retourne (mu, erreur standard
    de Laplace) ; NaN pour un groupe sans aucune valeur observée (pas de maximum).
    """
    mu = np.asarray(mu_start, float).copy()
    ok = yo_n > 0
    mu[~ok] = np.nan

    def derivs(m, s, ysum, k, nm):
        eta = b0 + b1 * (m[:, None] + s[:, None] * _GH_X)
        pm = expit(-eta)
        pdl = expit(eta) * pm
        G = pm @ _GH_W
        G1 = -b1 * (pdl @ _GH_W)                                   # dG/dmu
        G2 = -b1 ** 2 * ((pdl * (1 - 2 * expit(eta))) @ _GH_W)     # d²G/dmu²
        g = (ysum - k * m) / s ** 2 + nm * G1 / G
        h = -k / s ** 2 + nm * (G2 / G - (G1 / G) ** 2)
        return g, np.minimum(h, -1e-10)

    m, s = mu[ok], sigma[ok]
    ysum, k, nm = yo_sum[ok], yo_n[ok], nmis[ok]
    for _ in range(iters):
        g, h = derivs(m, s, ysum, k, nm)
        step = np.clip(-g / h, -2.0, 2.0)
        m = m + step
        if np.abs(step).max() < 1e-8:
            break
    _, h = derivs(m, s, ysum, k, nm)
    mu[ok] = m
    se = np.full_like(mu, np.nan)
    se[ok] = np.sqrt(-1.0 / h)
    return mu, se


def squeeze_within_group_variance(y, groups):
    """Variance intra-condition des valeurs OBSERVÉES, modérée (eBayes)."""
    y = np.asarray(y, float)
    rss = np.zeros(y.shape[0])
    df = np.zeros(y.shape[0])
    for g in np.unique(groups):
        yy = y[:, groups == g]
        k = (~np.isnan(yy)).sum(axis=1)
        m = np.nanmean(np.where(k[:, None] > 0, yy, 0.0), axis=1)
        rss += np.nansum((yy - m[:, None]) ** 2, axis=1)
        df += np.maximum(k - 1, 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        s2 = rss / df
    d0, s0 = fit_fdist(s2, df)
    if np.isinf(d0):
        return np.full(len(s2), s0)
    return np.where(df > 0, (df * np.nan_to_num(s2) + d0 * s0) / (df + d0), s0)


def estimate_dropout_prior(y, groups, beta0, beta1, p_min=0.98, min_cells=50,
                           n_grid=60, n_iter=500):
    """Loi a priori du taux de non-détection SPORADIQUE pi (échecs d'identification
    indépendants de l'abondance, mais aussi protéines présentes chez certains
    individus seulement), propre à chaque protéine.

    Estimée sur les cellules protéine × condition où la protéine est assez abondante
    pour que la DPC prédise une détection quasi certaine (P >= p_min à la moyenne
    observée) : leurs non-détections ne s'expliquent pas par l'abondance. Loi
    discrète sur une grille de pi, estimée par maximum de vraisemblance non
    paramétrique (EM ; vraisemblance binomiale tronquée : au moins une détection,
    sinon pas de moyenne observée). Une loi Bêta ne suffit pas : sur données
    réelles, la masse est presque toute près de 0 avec une petite fraction de
    protéines à pi de 0,2-0,5 (présence/absence selon l'individu), et ce sont ces
    dernières qui font les faux positifs. Retourne None si trop peu de cellules."""
    y = np.asarray(y, float)
    groups = np.asarray(groups)
    y_hi = (np.log(p_min / (1 - p_min)) - beta0) / beta1
    ms, ns = [], []
    for g in np.unique(groups):
        yy = y[:, groups == g]
        k = (~np.isnan(yy)).sum(axis=1)
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mo = np.nanmean(yy, axis=1)
        sel = (k > 0) & (mo >= y_hi)
        ms.append(yy.shape[1] - k[sel])
        ns.append(np.full(sel.sum(), yy.shape[1]))
    m, n = np.concatenate(ms), np.concatenate(ns)
    if len(m) < min_cells or n.max() < 2:
        return None
    tab = pd.DataFrame({"n": n, "m": m}).value_counts().reset_index(name="w")
    grid = np.r_[0.0, np.geomspace(1e-4, 0.95, n_grid - 1)]
    tm, tn, tw = tab["m"].to_numpy(), tab["n"].to_numpy(), tab["w"].to_numpy(float)
    L = stats.binom.pmf(tm[:, None], tn[:, None], grid[None, :])
    L = L / np.maximum(1 - grid[None, :] ** tn[:, None], 1e-300)       # tronquée : m < n
    g = np.full(len(grid), 1.0 / len(grid))
    for _ in range(n_iter):
        post = L * g[None, :]
        post /= post.sum(axis=1, keepdims=True)
        g_new = (tw[:, None] * post).sum(axis=0) / tw.sum()
        if np.abs(g_new - g).max() < 1e-10:
            g = g_new
            break
        g = g_new
    return {"grid": grid, "weights": g, "mean": float(g @ grid),
            "frac_above_0.1": float(g[grid > 0.1].sum()), "n_cells": int(len(m)),
            "p_min": p_min}


def differential_detection(y, groups, pairs, beta0, beta1, fdr_global=False,
                           dropout="auto"):
    """Pour chaque contraste (num, den) et chaque protéine : probabilité des
    non-détections d'un groupe si la protéine y avait l'abondance de l'autre
    groupe. p bilatérale = min(1, 2 min(p_num_absent, p_den_absent)). BH sur
    toutes les protéines du contraste (ou global).

    Probabilité de non-détection d'un échantillon d'abondance y :
        pi + (1 - pi) * (1 - DPC(y))
    pi = taux de non-détection sporadique PROPRE À LA PROTÉINE, loi a priori
    estimée sur les données (estimate_dropout_prior), a posteriori mise à jour par
    les non-détections du groupe présent (toutes comptées comme sporadiques : choix
    conservateur). Sans ce terme (dropout=None), une absence
    loin au-dessus du seuil est jugée quasi impossible : sur données réelles, des
    protéines abondantes manquant dans 3-4 échantillons sur 8 sans aucune
    différence de groupe sortaient alors très significatives (30-60 faux positifs
    par comparaison nulle 8 contre 8).

    y       : matrice log2 filtrée, NON imputée (NaN = non détecté)
    groups  : étiquette de groupe par colonne
    pairs   : liste de (nom_contraste, groupe_num, groupe_den)
    dropout : "auto" (estimé), None (modèle DPC seul), ou un résultat de
              estimate_dropout_prior
    """
    y = np.asarray(y, float)
    groups = np.asarray(groups)
    n_rows = y.shape[0]
    sigma = np.sqrt(squeeze_within_group_variance(y, groups))
    prior = (estimate_dropout_prior(y, groups, beta0, beta1) if dropout == "auto"
             else dropout)

    # Abondance (MAP sous DPC) et erreur standard, par protéine × groupe
    ab, ab_se, det, ntot, mobs = {}, {}, {}, {}, {}
    for g in np.unique(groups):
        yy = y[:, groups == g]
        k = (~np.isnan(yy)).sum(axis=1)
        n = yy.shape[1]
        ysum = np.nansum(yy, axis=1)
        mean_obs = np.where(k > 0, ysum / np.maximum(k, 1), np.nan)
        start = np.where(k > 0, mean_obs, np.nanmin(y))
        mu, se = _group_abundance_map(ysum, k.astype(float), (n - k).astype(float),
                                      sigma, beta0, beta1, start)
        ab[g], ab_se[g], det[g], ntot[g], mobs[g] = mu, se, k, n, mean_obs

    # Rang d'abondance (percentile) des moyennes observées de chaque protéine
    all_means = np.nanmean(y, axis=1)

    out = {}
    pvals = []
    for cname, gn, gd in pairs:
        p_dir = {}
        for absent, present in ((gn, gd), (gd, gn)):
            n_a = ntot[absent]
            m_a = n_a - det[absent]                       # non-détections
            # abondance du groupe présent, incertitude intégrée (quadrature)
            mu_nodes = ab[present][:, None] + ab_se[present][:, None] * _GH_X[None, :]
            pmis = _prob_missing(mu_nodes, sigma[:, None], beta0, beta1)    # (n_rows × 16)
            if prior is None:
                tail = stats.binom.sf(m_a[:, None] - 1, n_a, np.clip(pmis, 1e-300, 1.0))
                p = tail @ _GH_W
            else:
                # a posteriori de pi : grille a priori × vraisemblance du groupe présent
                # (toutes ses non-détections comptées comme sporadiques : conservateur)
                grid, gw = prior["grid"], prior["weights"]
                m_p = ntot[present] - det[present]
                w = gw[None, :] * stats.binom.pmf(m_p[:, None], ntot[present], grid[None, :])
                w /= np.maximum(w.sum(axis=1, keepdims=True), 1e-300)        # (n_rows × K)
                keep_k = gw > 1e-8
                q = grid[None, None, keep_k] + (1 - grid[None, None, keep_k]) * pmis[:, :, None]
                tail = stats.binom.sf(m_a[:, None, None] - 1, n_a, np.clip(q, 1e-300, 1.0))
                p = (tail * w[:, None, keep_k]).sum(axis=2) @ _GH_W
            p = np.where(m_a == 0, 1.0, p)
            p = np.where(det[present] == 0, 1.0, p)       # rien d'observé de l'autre côté
            p_dir[absent] = np.clip(p, 0.0, 1.0)
        p_two = np.minimum(1.0, 2.0 * np.minimum(p_dir[gn], p_dir[gd]))
        # sens : le groupe le MOINS détecté est l'absent
        higher_in = np.where(p_dir[gd] <= p_dir[gn], gn, gd)
        out[cname] = {"p": p_two, "higher_in": higher_in,
                      "det_num": det[gn], "n_num": ntot[gn],
                      "det_den": det[gd], "n_den": ntot[gd],
                      "mean_num": mobs[gn], "mean_den": mobs[gd],
                      "ab_num": ab[gn], "ab_den": ab[gd]}
        pvals.append(p_two)

    P = np.column_stack(pvals) if pvals else np.zeros((n_rows, 0))
    if P.size:
        if fdr_global:
            Padj = _bh_correction(P.ravel()).reshape(P.shape)
        else:
            Padj = np.column_stack([_bh_correction(P[:, k]) for k in range(P.shape[1])])
        for k, (cname, _, _) in enumerate(pairs):
            out[cname]["p_adj"] = Padj[:, k]

    pct = pd.Series(all_means).rank(pct=True).to_numpy()
    return {"by_contrast": out, "sigma": sigma, "abundance_pct": pct, "dropout": prior}


def plot_dpc(y_all, beta0, beta1, out_path, title="Detection probability curve"):
    """Proportion de détection par protéine en fonction de son abondance estimée
    (complétion par inclinaison exponentielle), et courbe ajustée. Retourne le
    chemin du PNG, ou None si matplotlib échoue."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    y = np.asarray(y_all, float)
    y = y[(~np.isnan(y)).any(axis=1)]
    yc = exp_tilt_complete(y, beta1)
    x = yc.mean(axis=1)
    prop = (~np.isnan(y)).mean(axis=1)
    sd = np.median(yc.std(axis=1, ddof=1))
    grid = np.linspace(np.quantile(x, 0.005), np.quantile(x, 0.995), 200)
    expected = expit(beta0 + beta1 * (grid[:, None] + sd * _GH_X)) @ _GH_W
    edges = np.quantile(x, np.linspace(0, 1, 21))
    bins = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, 19)
    bx = np.array([x[bins == b].mean() for b in range(20) if (bins == b).any()])
    by = np.array([prop[bins == b].mean() for b in range(20) if (bins == b).any()])
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.scatter(x, prop, s=6, alpha=0.15, color="#7f8c8d", lw=0, label="proteins")
    ax.plot(bx, by, "o", color="#2c3e50", ms=5, label="binned proportion")
    ax.plot(grid, expected, color="#c0392b", lw=2,
            label=f"fitted DPC: logit p = {beta0:.2f} + {beta1:.2f} y")
    ax.axvline(-beta0 / beta1, color="#c0392b", ls=":", lw=1)
    ax.set_xlabel("log2 intensity (estimated, missing values included)")
    ax.set_ylabel("Proportion of samples detected")
    ax.set_ylim(-0.03, 1.03)
    ax.set_title(title)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# ==============================================================================
# 6. Orchestration pour Deimos : tableaux prêts pour l'Excel
# ==============================================================================

def _pairs_from_contrasts(contrast_mat, contrast_names, n_groups):
    """(nom, indice groupe num, indice groupe den) pour les contrastes par paires."""
    pairs = []
    for k, c in enumerate(contrast_names):
        col = contrast_mat[:n_groups, k]
        pos, neg = np.where(np.isclose(col, 1))[0], np.where(np.isclose(col, -1))[0]
        if len(pos) == 1 and len(neg) == 1 and np.isclose(np.abs(col).sum(), 2):
            pairs.append((c, int(pos[0]), int(neg[0])))
    return pairs


def _volcano_call(p, lfc, params):
    lfc_min = params.get("volcano_lfc_min", float(np.log2(params.get("volcano_ratio_min", 1.5))))
    return (np.asarray(p) < params.get("volcano_p_thresh", 0.05)) & \
           (np.abs(np.asarray(lfc)) >= lfc_min)


def run_missingness_sensitivity(mat_all, mat_filt, df_results, design_mat,
                                contrast_mat, contrast_names, n_groups, params,
                                do_dpc=True, do_detection=True, verbose=True,
                                precursors=None, primary=None):
    """Les deux analyses de sensibilité. Ne modifie aucun objet reçu.

    mat_all   : matrice log2 NON filtrée (lignes = protéines quantifiées), pour la DPC
    mat_filt  : matrice log2 filtrée NON imputée, mêmes lignes et ordre que df_results
    design_mat, contrast_mat, contrast_names : ceux de l'analyse principale
    n_groups  : nombre de colonnes « condition » en tête de design_mat
    precursors: None, ou {"Y": matrice log2 des précurseurs (colonnes dans l'ordre
                de mat_filt), "protein_ids": Protein.Group de chaque précurseur}.
                Si fournie, l'analyse sans imputation quantifie les protéines à
                partir des précurseurs (logique dpcQuant), comme limpa le recommande ;
                les protéines absentes du pr_matrix restent au niveau protéine.
    primary   : résultat de dpc_primary_analysis déjà calculé (statistique
                principale DPC) ; évite de refaire la quantification.
    """
    Yall = np.asarray(mat_all, float)
    Y = np.asarray(mat_filt, float)
    X = np.asarray(design_mat, float)
    if len(df_results) != Y.shape[0]:
        raise ValueError("df_results and the filtered matrix are not aligned")
    use_padj = bool(params.get("volcano_use_padj"))
    pkey = "p.adj" if use_padj else "p.val"
    fdr_global = bool(params.get("fdr_global"))
    out = {"pkey": pkey}

    dpc = primary["dpc"] if primary is not None else estimate_dpc(Yall)
    b0, b1 = dpc["beta0"], dpc["beta1"]
    out["dpc"] = dpc
    if verbose:
        lod50 = -b0 / b1
        print(f"  [DPC] logit P(detected) = {b0:.2f} + {b1:.2f} x log2 intensity "
              f"(50 % detection at log2 = {lod50:.1f}; {dpc['n_rows']} rows)")

    ids = [c for c in ("name", "Protein.Group", "Genes") if c in df_results.columns]
    base = df_results[ids].reset_index(drop=True)
    imputed_n = np.isnan(Y).sum(axis=1)

    # ---- 1. Analyse sans imputation -----------------------------------------
    # Niveau précurseur pour les protéines présentes dans le pr_matrix, niveau
    # protéine pour les autres (dpc_primary_analysis). Si l'analyse a déjà été
    # faite comme statistique principale, elle est réutilisée telle quelle.
    dpc_res = None
    if do_dpc:
        if primary is None:
            pg = (df_results["Protein.Group"].astype(str).to_numpy()
                  if "Protein.Group" in df_results.columns else None)
            primary = dpc_primary_analysis(Yall, Y, pg, X, contrast_mat, fdr_global=fdr_global,
                                           precursors=precursors if pg is not None else None,
                                           verbose=verbose)
        level = primary["level"]
        dpc_res = {k: primary[k] for k in ("logFC", "t", "p_value", "p_adj", "n_df_corrected",
                                           "df_prior", "s2_prior")}
        if primary.get("dpc_precursor"):
            out["dpc_precursor"] = primary["dpc_precursor"]
        out["dpc_level"] = level
        sheet = base.copy()
        sheet["n_imputed"] = imputed_n
        if level == "precursor":
            sheet["dpc_source"] = primary["source"]
            sheet["n_precursors"] = primary["n_precursors"]
        for k, c in enumerate(contrast_names):
            sheet[f"{c}_logFC_dpc"] = dpc_res["logFC"][:, k]
            sheet[f"{c}_p.val_dpc"] = dpc_res["p_value"][:, k]
            sheet[f"{c}_p.adj_dpc"] = dpc_res["p_adj"][:, k]
            prim = _volcano_call(df_results[f"{c}_{pkey}"], df_results[f"{c}_diff"], params)
            alt = _volcano_call(dpc_res["p_adj" if use_padj else "p_value"][:, k],
                                dpc_res["logFC"][:, k], params)
            na = np.isnan(dpc_res["p_value"][:, k])
            sheet[f"{c}_concordance"] = np.select(
                [na, prim & alt, prim & ~alt, ~prim & alt],
                ["not assessed", "both", "imputation only", "DPC only"],
                default="")
        out["dpc_sheet"] = sheet
        out["dpc_de"] = dpc_res
        out["dpc_prior"] = primary["prior"]
        if verbose:
            print(f"  [OK] Imputation-free analysis (DPC, {level} level): "
                  f"{dpc_res['n_df_corrected']} protein(s) with an unobserved group "
                  f"(residual df corrected)")

    # ---- 2. Détection différentielle ----------------------------------------
    if do_detection:
        groups = np.argmax(X[:, :n_groups], axis=1)
        pairs = _pairs_from_contrasts(contrast_mat, contrast_names, n_groups)
        dd = differential_detection(Y, groups, pairs, b0, b1, fdr_global=fdr_global)
        blocks = []
        for c, gi, gj in pairs:
            r = dd["by_contrast"][c]
            num, _, den = c.partition("_vs_")
            prim = _volcano_call(df_results[f"{c}_{pkey}"], df_results[f"{c}_diff"], params)
            kn, kd = r["det_num"], r["det_den"]
            n_n, n_d = r["n_num"], r["n_den"]
            # motif on/off : un groupe détecté (>= max(2, n-1)), l'autre dans moins
            # de la moitié de ses échantillons
            pres_n, pres_d = kn >= max(2, n_n - 1), kd >= max(2, n_d - 1)
            abs_n, abs_d = kn <= (n_n - 1) // 2, kd <= (n_d - 1) // 2
            onoff = (pres_n & abs_d) | (pres_d & abs_n)
            sig_det = r["p_adj"] < 0.05
            keep = np.where(onoff | sig_det)[0]
            if not len(keep):
                continue
            num_higher = r["higher_in"][keep] == gi
            absent_k = np.where(num_higher, kd[keep], kn[keep])
            pattern = np.where(~onoff[keep], "partial",
                               np.where(absent_k == 0, "on/off", "mostly absent"))
            verdict = np.select(
                [sig_det[keep] & prim[keep], sig_det[keep], prim[keep]],
                ["confirmed", "detection only (missed by the imputed analysis)",
                 "imputation only (absence compatible with the detection limit)"],
                default="not significant")
            b = base.iloc[keep].reset_index(drop=True)
            b["contrast"] = c
            b["pattern"] = pattern
            b["higher_in"] = np.where(num_higher, num, den)
            b["detected_num"] = [f"{k}/{n_n}" for k in kn[keep]]
            b["detected_den"] = [f"{k}/{n_d}" for k in kd[keep]]
            b["mean_intensity_present"] = np.where(num_higher, r["mean_num"][keep],
                                                   r["mean_den"][keep])
            b["abundance_rank_pct"] = np.round(100 * dd["abundance_pct"][keep], 1)
            b["p.val_detection"] = r["p"][keep]
            b["p.adj_detection"] = r["p_adj"][keep]
            b["logFC_imputed"] = df_results[f"{c}_diff"].to_numpy()[keep]
            b[f"{pkey}_imputed"] = df_results[f"{c}_{pkey}"].to_numpy()[keep]
            b["significant_imputed"] = np.where(prim[keep], "yes", "no")
            if dpc_res is not None:
                k = contrast_names.index(c)
                b["logFC_dpc"] = dpc_res["logFC"][keep, k]
                b["p.adj_dpc"] = dpc_res["p_adj"][keep, k]
            b["verdict"] = verdict
            blocks.append(b)
        out["detection_sheet"] = (pd.concat(blocks, ignore_index=True) if blocks
                                  else pd.DataFrame(columns=ids + ["contrast", "pattern",
                                                                   "p.adj_detection", "verdict"]))
        out["detection"] = dd
        out["dropout_prior"] = dd.get("dropout")
        if verbose:
            if dd.get("dropout"):
                print(f"  [DPC] sporadic dropout prior: mean {dd['dropout']['mean']:.3f}, "
                      f"{100 * dd['dropout']['frac_above_0.1']:.1f} % of proteins with "
                      f"rate > 0.1 ({dd['dropout']['n_cells']} abundant protein x "
                      f"condition cells)")
            d = out["detection_sheet"]
            n_sig = int((d["p.adj_detection"] < 0.05).sum()) if len(d) else 0
            print(f"  [OK] Differential detection: {len(d)} protein x contrast pair(s) "
                  f"reported, {n_sig} significant (p.adj < 0.05)")

    # ---- Résumé par contraste ----------------------------------------------
    summ = []
    for c in contrast_names:
        prim = _volcano_call(df_results[f"{c}_{pkey}"], df_results[f"{c}_diff"], params)
        rec = {"contrast": c, "count_sig_imputed": int(prim.sum())}
        if dpc_res is not None:
            conc = out["dpc_sheet"][f"{c}_concordance"]
            rec.update({"count_sig_dpc": int(conc.isin(["both", "DPC only"]).sum()),
                        "count_both": int((conc == "both").sum()),
                        "count_imputation_only": int((conc == "imputation only").sum()),
                        "count_dpc_only": int((conc == "DPC only").sum())})
            rob = f"Robustness_Score_{c}"
            if rob in df_results.columns:
                rs = df_results[rob].to_numpy()
                for lab, key in (("both", "robust_median_both"),
                                 ("imputation only", "robust_median_imputation_only")):
                    m = (conc == lab).to_numpy()
                    rec[key] = float(np.median(rs[m])) if m.any() else np.nan
        if do_detection and "detection_sheet" in out and len(out["detection_sheet"]):
            d = out["detection_sheet"][out["detection_sheet"]["contrast"] == c]
            rec.update({
                "count_onoff_pairs": int(d["pattern"].isin(["on/off", "mostly absent"]).sum()),
                "count_detection_sig": int((d["p.adj_detection"] < 0.05).sum()),
                "count_detection_only": int(d["verdict"].str.startswith("detection only").sum()),
                "count_imputation_only_onoff": int(d["verdict"].str.startswith("imputation only").sum())})
        summ.append(rec)
    out["summary"] = pd.DataFrame(summ)
    return out


# ==============================================================================
# 7. Statistique PRINCIPALE sans imputation (option primary_statistic: dpc)
# ==============================================================================
# Même chaîne que l'analyse de sensibilité ci-dessus (logique limpa : DPC ->
# quantification -> limma pondéré), organisée pour servir de statistique
# principale :
#   - niveau précurseur (dpcQuant) pour toutes les protéines présentes dans le
#     pr_matrix ; repli au niveau protéine (dpcQuantByRow) pour les autres, avec
#     un ajustement limma estimé sur toutes les protéines (hyperparamètres
#     stables), dont on ne garde que les lignes de repli ;
#   - correction BH commune à toutes les protéines (par contraste, ou globale) ;
#   - score de robustesse « leave-one-sample-out » (aucun tirage aléatoire ici,
#     le score par ré-imputation n'a donc pas d'équivalent direct).

def _bh_matrix(pval, fdr_global=False):
    """BH par contraste (colonnes) ou global, en ignorant les NaN."""
    pval = np.asarray(pval, float)
    out = np.full_like(pval, np.nan)
    if fdr_global:
        ok = np.isfinite(pval)
        if ok.any():
            out[ok] = _bh_correction(pval[ok])
        return out
    for k in range(pval.shape[1]):
        ok = np.isfinite(pval[:, k])
        if ok.any():
            out[ok, k] = _bh_correction(pval[ok, k])
    return out


def _fit_dpc_sources(blocks, n_prot, X, C, keep_cols=None):
    """dpc_de sur chaque source (précurseur / protéine), résultats replacés dans
    l'ordre des protéines. keep_cols : échantillons conservés (leave-one-out)."""
    n_c = C.shape[1]
    lfc = np.full((n_prot, n_c), np.nan)
    tt = np.full((n_prot, n_c), np.nan)
    pv = np.full((n_prot, n_c), np.nan)
    fits = {}
    for src, b in blocks.items():
        E, SE, imp, Xs = b["E"], b["SE"], b["imputed"], X
        if keep_cols is not None:
            E, SE, imp, Xs = E[:, keep_cols], SE[:, keep_cols], imp[:, keep_cols], X[keep_cols]
        r = dpc_de(E, SE, imp, Xs, C)
        take, rows = b["take"], b["rows"]
        lfc[rows], tt[rows], pv[rows] = r["logFC"][take], r["t"][take], r["p_value"][take]
        fits[src] = r
    return lfc, tt, pv, fits


def dpc_primary_analysis(mat_all, mat_filt, protein_groups, design_mat, contrast_mat,
                         fdr_global=False, precursors=None, verbose=True) -> dict:
    """Analyse différentielle sans imputation, pour servir de statistique principale.

    mat_all        : matrice log2 protéines NON filtrée (estimation de la DPC protéine)
    mat_filt       : matrice log2 filtrée non imputée (lignes = protéines analysées)
    protein_groups : Protein.Group de chaque ligne de mat_filt (lien avec le pr_matrix)
    precursors     : None ou {"Y": log2 précurseurs (colonnes = celles de mat_filt),
                     "protein_ids": Protein.Group de chaque précurseur}

    Retourne logFC, t, p_value, p_adj (protéines × contrastes), source par protéine
    ('precursor' / 'protein'), E / SE / imputed (protéines × échantillons) de la
    source retenue, paramètres des DPC et les blocs nécessaires au leave-one-out."""
    Yall = np.asarray(mat_all, float)
    Y = np.asarray(mat_filt, float)
    X = np.asarray(design_mat, float)
    C = np.asarray(contrast_mat, float)
    n_prot, n = Y.shape
    dpc = estimate_dpc(Yall)
    out = {"dpc": dpc, "dpc_precursor": None}
    blocks = {}
    has = np.zeros(n_prot, bool)
    n_prec = np.full(n_prot, np.nan)
    prior = None
    if precursors is not None and protein_groups is not None:
        Ypr = np.asarray(precursors["Y"], float)
        okr = (~np.isnan(Ypr)).any(axis=1)
        if okr.sum() >= 10:
            dpr = estimate_dpc(Ypr[okr])
            q = dpc_quant_precursors(Ypr[okr], np.asarray(precursors["protein_ids"])[okr],
                                     dpr["beta0"], dpr["beta1"])
            pos = pd.Series(np.arange(len(q["proteins"])), index=q["proteins"])
            ix = pos.reindex(np.asarray(protein_groups).astype(str))
            has = ix.notna().to_numpy()
            if has.sum() >= 3:
                ii = ix[has].astype(int).to_numpy()
                blocks["precursor"] = {"rows": np.where(has)[0], "take": np.arange(has.sum()),
                                       "E": q["E"][ii], "SE": q["SE"][ii],
                                       "imputed": q["n_obs"][ii] == 0}
                n_prec[has] = q["n_prec"][ii]
                prior = {k: q["hyper"][k] for k in ("prior_mean", "prior_sd", "prior_logfc")}
                out["dpc_precursor"] = {
                    "beta0": dpr["beta0"], "beta1": dpr["beta1"], "n_rows": dpr["n_rows"],
                    "n_precursors": int(okr.sum()), "n_proteins_quantified": int(has.sum()),
                    "n_proteins_missing": int((~has).sum())}
            else:
                has[:] = False
    if (~has).any():
        comp = dpc_complete(Y, dpc["beta0"], dpc["beta1"])
        blocks["protein"] = {"rows": np.where(~has)[0], "take": np.where(~has)[0],
                             "E": comp["E"], "SE": comp["SE"], "imputed": comp["imputed"]}
        if prior is None:
            prior = {k: comp[k] for k in ("prior_mean", "prior_sd", "prior_logfc")}
    lfc, tt, pv, fits = _fit_dpc_sources(blocks, n_prot, X, C)
    E = np.full((n_prot, n), np.nan)
    SE = np.full((n_prot, n), np.nan)
    imp = np.zeros((n_prot, n), bool)
    for b in blocks.values():
        E[b["rows"]], SE[b["rows"]], imp[b["rows"]] = (b["E"][b["take"]], b["SE"][b["take"]],
                                                        b["imputed"][b["take"]])
    main = fits.get("precursor") or fits["protein"]
    # protéines gardées dont les ddl résiduels ont été corrigés (groupe non observé)
    n_df = sum(int(r["df_corrected"][blocks[src]["take"]].sum()) for src, r in fits.items())
    out.update({
        "logFC": lfc, "t": tt, "p_value": pv, "p_adj": _bh_matrix(pv, fdr_global),
        "source": np.where(has, "precursor", "protein"), "n_precursors": n_prec,
        "E": E, "SE": SE, "imputed": imp, "level": "precursor" if "precursor" in blocks
        else "protein", "n_df_corrected": n_df, "df_prior": main["df_prior"],
        "s2_prior": main["s2_prior"], "prior": prior, "_blocks": blocks,
        "_design": X, "_contrast": C})
    if verbose:
        lod = -dpc["beta0"] / dpc["beta1"]
        if out["dpc_precursor"]:
            d = out["dpc_precursor"]
            print(f"  [DPC] precursor level: logit P(detected) = {d['beta0']:.2f} + "
                  f"{d['beta1']:.2f} x log2 intensity; {d['n_precursors']} precursors -> "
                  f"{d['n_proteins_quantified']} proteins"
                  + (f" ({d['n_proteins_missing']} at protein level, no precursor)"
                     if d["n_proteins_missing"] else ""))
        else:
            print(f"  [DPC] protein level: logit P(detected) = {dpc['beta0']:.2f} + "
                  f"{dpc['beta1']:.2f} x log2 intensity (50 % detection at log2 = {lod:.1f})")
    return out


def dpc_loo_robustness(primary, n_groups, params, fdr_global=False, n_workers=None,
                       verbose=True) -> dict:
    """Score de robustesse de la statistique DPC : pourcentage des ajustements
    « leave-one-sample-out » où la protéine passe encore les seuils des volcanos.

    Pour un contraste, seuls les échantillons des conditions comparées sont retirés
    tour à tour (un retrait ailleurs ne change que l'a priori de variance). Les
    abondances et erreurs standard de la quantification complète sont conservées :
    on mesure la dépendance de l'appel à un échantillon (injection aberrante,
    réplicat qui porte seul l'effet), sans re-quantifier. Retrait ignoré s'il rend
    le design non estimable (condition à un seul échantillon)."""
    blocks = primary["_blocks"]
    X, C = primary["_design"], primary["_contrast"]
    n, n_c = X.shape[0], C.shape[1]
    n_prot = primary["logFC"].shape[0]
    use_padj = bool(params.get("volcano_use_padj"))
    p_thr = float(params.get("volcano_p_thresh", 0.05))
    lfc_min = float(params.get("volcano_lfc_min",
                               np.log2(params.get("volcano_ratio_min", 1.5))))
    cond = np.argmax(X[:, :n_groups], axis=1)
    rel = [set(np.where(~np.isclose(C[:n_groups, k], 0))[0]) for k in range(n_c)]
    todo = [s for s in range(n) if any(cond[s] in r for r in rel)]
    full_rank = np.linalg.matrix_rank(X)

    def one(s):
        keep = np.delete(np.arange(n), s)
        if np.linalg.matrix_rank(X[keep]) < full_rank:
            return s, None
        lf, _, pv, _ = _fit_dpc_sources(blocks, n_prot, X, C, keep_cols=keep)
        pk = _bh_matrix(pv, fdr_global) if use_padj else pv
        return s, (pk < p_thr) & (np.abs(lf) >= lfc_min)

    success = np.zeros((n_prot, n_c))
    n_fits = np.zeros(n_c, int)
    results = []
    try:
        from concurrent.futures import ThreadPoolExecutor
        import os as _os
        nw = n_workers or min(4, _os.cpu_count() or 1)
        with ThreadPoolExecutor(max_workers=nw) as ex:
            results = list(ex.map(one, todo))
    except Exception:
        results = [one(s) for s in todo]
    skipped = 0
    for s, call in results:
        if call is None:
            skipped += 1
            continue
        for k in range(n_c):
            if cond[s] in rel[k]:
                success[:, k] += call[:, k]
                n_fits[k] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        score = np.where(n_fits > 0, np.round(100.0 * success / np.maximum(n_fits, 1)), np.nan)
    if verbose:
        print(f"  [OK] DPC robustness (leave-one-sample-out): {len(todo) - skipped} refits"
              + (f", {skipped} sample(s) skipped (design not estimable without them)"
                 if skipped else ""))
    return {"score": score, "n_fits": n_fits, "n_skipped": skipped}
