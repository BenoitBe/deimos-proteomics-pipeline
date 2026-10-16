# ==============================================================================
# limma_ebayes.py — Implémentation Python de l'eBayes de limma (Smyth 2004)
# Pipeline protéomique — Proteogen
# ==============================================================================
# Reproduit :
#   lmFit()         → régression OLS par protéine sur design matrix
#   makeContrasts() → combinaisons linéaires de coefficients
#   contrasts.fit() → recompute sur les contrastes
#   eBayes()        → modération empirique (emprunt de variance entre features)
#   topTable()      → statistiques finales avec p-values ajustées
# ==============================================================================

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma
from itertools import combinations


# ------------------------------------------------------------------------------
# 1. lmFit — régression OLS sur la design matrix
# ------------------------------------------------------------------------------
def lm_fit(expr_matrix: np.ndarray, design: np.ndarray):
    """
    Régression OLS protéine par protéine.

    Parameters
    ----------
    expr_matrix : np.ndarray (n_proteins × n_samples)
    design      : np.ndarray (n_samples × n_coefficients)

    Returns
    -------
    dict avec : coefficients, stdev_unscaled, sigma, df_residual
    """
    n_proteins, n_samples = expr_matrix.shape
    n_coef = design.shape[1]
    # ddl résiduels = n - rang(X), comme lmFit (QR) : un design de rang
    # incomplet ne doit pas gonfler les ddl.
    rank = int(np.linalg.matrix_rank(design))
    if rank < n_coef:
        print(f"  [WARN] Design matrix not of full rank ({rank} < {n_coef} columns): "
              "some coefficients are not estimable.")
    df_residual = n_samples - rank

    # Pseudo-inverse de la design matrix (stable même si colonnes corrélées)
    XtX_inv = np.linalg.pinv(design.T @ design)
    hat = XtX_inv @ design.T  # (n_coef × n_samples)

    coefficients = expr_matrix @ hat.T  # (n_proteins × n_coef)
    fitted = coefficients @ design.T    # (n_proteins × n_samples)
    residuals = expr_matrix - fitted    # (n_proteins × n_samples)

    # sigma² par protéine
    rss = np.sum(residuals ** 2, axis=1)
    sigma2 = rss / df_residual
    sigma = np.sqrt(np.maximum(sigma2, 0))

    # stdev_unscaled : écart-type non-scalé des coefficients (sans sigma)
    stdev_unscaled = np.sqrt(np.diag(XtX_inv))  # (n_coef,)
    stdev_unscaled_mat = np.tile(stdev_unscaled, (n_proteins, 1))  # broadcast

    return {
        "coefficients": coefficients,
        "stdev_unscaled": stdev_unscaled_mat,
        "sigma": sigma,
        "df_residual": df_residual,
        "design": design,
    }


# ------------------------------------------------------------------------------
# 2. contrasts_fit — projette le fit sur les contrastes
# ------------------------------------------------------------------------------
def contrasts_fit(fit: dict, contrast_matrix: np.ndarray):
    """
    Recalcule coefficients et stdev_unscaled pour chaque contraste.

    Parameters
    ----------
    fit             : dict retourné par lm_fit
    contrast_matrix : np.ndarray (n_coef × n_contrasts)

    Returns
    -------
    dict enrichi avec les statistiques par contraste
    """
    coef = fit["coefficients"]       # (n_proteins × n_coef)
    su   = fit["stdev_unscaled"]     # (n_proteins × n_coef)
    C    = contrast_matrix           # (n_coef × n_contrasts)

    # Nouveaux coefficients (logFC par contraste)
    new_coef = coef @ C              # (n_proteins × n_contrasts)

    # Variance non-scalée par contraste : diag(C^T * XtX_inv * C)
    XtX_inv = np.linalg.pinv(fit["design"].T @ fit["design"])
    stdev_contrasts = np.sqrt(np.diag(C.T @ XtX_inv @ C))  # (n_contrasts,)
    new_su = np.tile(stdev_contrasts, (coef.shape[0], 1))   # broadcast

    fit_c = fit.copy()
    fit_c["coefficients"] = new_coef
    fit_c["stdev_unscaled"] = new_su
    return fit_c


# ------------------------------------------------------------------------------
# 3. eBayes — modération empirique de Bayes (Smyth 2004)
# ------------------------------------------------------------------------------
def ebayes(fit: dict, fdr_global: bool = False):
    """
    Modération empirique de Bayes des variances résiduelles.
    Estime les hyperparamètres (d0, s0²) par fitting de la distribution
    a priori sur les sigma² observés.

    fdr_global : si True, la correction BH est appliquée sur l'ENSEMBLE des
        p-values (tous contrastes confondus, une seule famille de tests) au
        lieu d'une correction indépendante par contraste. Plus conservateur :
        contrôle le FDR au niveau de toute l'étude.

    Returns
    -------
    dict enrichi : t_stat, p_value, p_adj (BH), logFC, s2_post, df_total
    """
    sigma   = np.asarray(fit["sigma"], float)   # (n_proteins,)
    su      = fit["stdev_unscaled"]             # (n_proteins × n_contrasts)
    df_res  = float(fit["df_residual"])
    coef    = fit["coefficients"]               # (n_proteins × n_contrasts)

    n_proteins, n_contrasts = coef.shape

    # --- Hyperparamètres d0, s0² : limma::squeezeVar -> fitFDist (ddl égaux) ---
    # Reproduction exacte de limma 3.58 (Smyth 2004) :
    #   x = max(s², 1e-5 × médiane(s²))   (variances nulles décalées, pas exclues)
    #   e = log(x) - digamma(df/2) + log(df/2)
    #   evar = var(e) - trigamma(df/2)
    #   evar > 0  : d0 = 2 trigamma^-1(evar), s0² = exp(mean(e) + digamma(d0/2) - log(d0/2))
    #   evar <= 0 : d0 = Inf, s0² = mean(x)  (variance poolée, limma >= 2017)
    s2 = sigma ** 2
    ok = np.isfinite(s2) & (df_res > 1e-15)
    d0, s0_sq = _fit_fdist(s2[ok], df_res)

    # --- Variance postérieure modérée (sur les s² d'origine) ---
    if np.isinf(d0):
        s2_post = np.full_like(s2, s0_sq)
    else:
        s2_post = (df_res * s2 + d0 * s0_sq) / (df_res + d0)   # shrinkage vers s0²
    s2_post[~np.isfinite(s2)] = np.nan
    # ddl totaux plafonnés aux ddl poolés (eBayes : pmin(df.total, df.pooled))
    df_pooled = df_res * int(ok.sum())
    df_total = float(min(df_res + d0, df_pooled))

    # --- t-statistiques modérées et p-values ---
    with np.errstate(divide="ignore", invalid="ignore"):
        t_stat = coef / su / np.sqrt(s2_post)[:, None]
    p_value = 2 * stats.t.sf(np.abs(t_stat), df=df_total)
    if fdr_global:
        # Correction BH GLOBALE : toutes les p-values en une seule famille
        p_adj = _bh_correction(p_value.ravel()).reshape(p_value.shape)
    else:
        # Correction BH par contraste (une famille de tests par comparaison)
        p_adj = np.column_stack([_bh_correction(p_value[:, j])
                                 for j in range(n_contrasts)])

    fit_eb = fit.copy()
    fit_eb.update({
        "t_stat":   t_stat,
        "p_value":  p_value,
        "p_adj":    p_adj,
        "s2_post":  s2_post,
        "df_total": df_total,
        "d0":       d0,
        "s0_sq":    s0_sq,
    })
    return fit_eb


# ------------------------------------------------------------------------------
# 3bis. spectra_count_ebayes — modération DEqMS (Zhu et al. 2020)
# ------------------------------------------------------------------------------
def spectra_count_ebayes(fit: dict, pep_count: np.ndarray,
                         loess_span: float = 0.75, fdr_global: bool = False):
    """
    Reproduction de DEqMS::spectraCounteBayes (Zhu et al. 2020, DEqMS 1.31,
    fit.method = "loess") puis BH comme outputDEqMS.

    La variance a priori dépend du nombre de peptides quantifiés :
      1. logVAR = log(s²) ; tendance loess(logVAR ~ log2(count), span = 0.75),
         paramètres par défaut de R (degré 2, gaussien, surface interpolée) ;
      2. eg = logVAR - digamma(df/2) + log(df/2), egpred idem sur la tendance ;
         m = moyenne((eg - egpred)² - trigamma(df/2)) ;
      3. d0 = première valeur de la grille 0.1, 0.2, ... qui minimise
         |m - trigamma(d0/2)| (recherche séquentielle de DEqMS) ;
      4. s0²(count) = exp(egpred + digamma(d0/2) - log(d0/2)) ;
      5. variance postérieure (d0 s0² + df s²)/(d0 + df), t sur d0 + df ddl.

    La loess est celle de R (même code Fortran, via scikit-misc) quand le paquet
    est installé ; sinon une loess directe (approximation, voir _loess_direct).

    Parameters
    ----------
    fit       : dict retourné par lm_fit + contrasts_fit (avant ou après ebayes)
    pep_count : np.ndarray (n_proteins,) — comptage peptidique par protéine
                (déjà augmenté du pseudo-comptage de 1, comme dans la vignette
                LFQ de DEqMS) ; NaN ou < 1 = protéine exclue
    loess_span: span de la loess (0.75 dans DEqMS)

    Returns
    -------
    dict enrichi : t_stat, p_value, p_adj, s2_post, s2_prior, df_prior, df_total,
                   pep_count, count_valid, loess_fit, loess_engine
    """
    sigma   = np.asarray(fit["sigma"], float)
    su      = fit["stdev_unscaled"]           # (n_proteins × n_contrasts)
    df_res  = float(fit["df_residual"])
    coef    = fit["coefficients"]             # (n_proteins × n_contrasts)
    n_proteins, n_contrasts = coef.shape

    s2 = sigma ** 2
    pep_count = np.asarray(pep_count, dtype=float)
    valid = (s2 > 0) & np.isfinite(s2) & (pep_count >= 1) & np.isfinite(pep_count)
    if valid.sum() < 10:
        raise ValueError("too few proteins with a peptide count for DEqMS")

    log_var = np.full(n_proteins, np.nan)
    log_var[valid] = np.log(s2[valid])
    x = np.full(n_proteins, np.nan)
    x[valid] = np.log2(pep_count[valid])

    # --- 1. Tendance loess de log(s²) sur log2(comptage) ---
    y_pred = np.full(n_proteins, np.nan)
    y_pred[valid], engine = _loess_fit(x[valid], log_var[valid], span=loess_span)
    if not np.all(np.isfinite(y_pred[valid])):
        raise ValueError("loess fit failed (degenerate peptide counts)")

    # --- 2-3. Degrés de liberté a priori (grille séquentielle de DEqMS) ---
    shift = -digamma(df_res / 2.0) + np.log(df_res / 2.0)
    eg = log_var + shift
    egpred = y_pred + shift
    mean_myfct = float(np.nanmean((eg - egpred) ** 2 - polygamma(1, df_res / 2.0)))
    df_prior = _deqms_prior_df(mean_myfct, int(valid.sum()))

    # --- 4-5. Variance a priori par protéine et postérieure ---
    s2_prior = np.exp(egpred + digamma(df_prior / 2.0) - np.log(df_prior / 2.0))
    df_total = df_res + df_prior
    s2_post = (df_prior * s2_prior + df_res * s2) / df_total

    with np.errstate(divide="ignore", invalid="ignore"):
        t_stat = coef / (su * np.sqrt(s2_post)[:, None])
    p_value = 2 * stats.t.sf(np.abs(t_stat), df=df_total)
    t_stat[~valid, :] = np.nan
    p_value[~valid, :] = np.nan

    if fdr_global:
        p_adj = _bh_correction(p_value.ravel()).reshape(p_value.shape)
    else:
        p_adj = np.column_stack([_bh_correction(p_value[:, j])
                                 for j in range(n_contrasts)])

    fit_dq = fit.copy()
    fit_dq.update({
        "t_stat":   t_stat,
        "p_value":  p_value,
        "p_adj":    p_adj,
        "s2_post":  s2_post,
        "s2_prior": s2_prior,
        "df_prior": float(df_prior),
        "df_total": float(df_total),
        "pep_count": pep_count,
        "count_valid": valid,
        "loess_fit": y_pred,
        "loess_engine": engine,
    })
    return fit_dq


def _deqms_prior_df(mean_myfct: float, numgenes: int) -> float:
    """Recherche de d0 de DEqMS : grille i/10 (i = 1..10 x numgenes), arrêt dès que
    la fonction |mean_myfct - trigamma(d0/2)| remonte, minimum des valeurs vues
    (premier en cas d'égalité), comme la boucle R d'origine."""
    if not np.isfinite(mean_myfct):
        raise ValueError("DEqMS prior df: non-finite variance residuals")
    n_max = max(numgenes * 10, 3)
    best_i, prev2, prev1 = None, None, None
    vals = []
    # parcours par blocs pour rester rapide sans changer la règle d'arrêt
    start = 1
    while start <= n_max:
        stop = min(start + 4096, n_max + 1)
        d0 = np.arange(start, stop) / 10.0
        v = np.abs(mean_myfct - polygamma(1, d0 / 2.0))
        for k, val in enumerate(v):
            i = start + k
            vals.append(val)
            if i > 2 and vals[i - 3] < vals[i - 2]:      # priordf[i-2] < priordf[i-1]
                arr = np.asarray(vals)
                return float((int(np.argmin(arr)) + 1) / 10.0)
        start = stop
    arr = np.asarray(vals)
    return float((int(np.argmin(arr)) + 1) / 10.0)


def _loess_fit(x: np.ndarray, y: np.ndarray, span: float = 0.75):
    """Valeurs ajustées de loess(y ~ x, span) avec les réglages par défaut de R.
    Utilise scikit-misc (même implémentation Fortran que R : dloess de Cleveland,
    Grosse et Shyu) ; à défaut, une loess directe (même modèle, sans
    l'interpolation par arbre k-d : écarts faibles). Retourne (fitted, moteur)."""
    try:
        from skmisc.loess import loess as _skloess
    except ImportError:
        return _loess_direct(x, y, span=span), "direct loess (scikit-misc not installed)"
    try:
        m = _skloess(np.asarray(x, float), np.asarray(y, float), span=span, degree=2,
                     family="gaussian", surface="interpolate", statistics="approximate")
        m.fit()
        return np.asarray(m.outputs.fitted_values, float), "R loess (scikit-misc)"
    except Exception as e:
        # Comptages très dégénérés (la plupart des protéines au même comptage) :
        # R n'émet qu'un avertissement (pseudo-inverse), scikit-misc lève une
        # erreur. On passe à la loess directe, finie dans tous les cas.
        return (_loess_direct(x, y, span=span),
                "direct loess (R loess degenerate on these peptide counts)")


def _loess_direct(x: np.ndarray, y: np.ndarray, span: float = 0.75,
                  degree: int = 2) -> np.ndarray:
    """Loess locale quadratique pondérée tricube (surface « direct » de R) :
    voisinage des q = floor(n x span) points les plus proches."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    q = max(int(np.floor(n * span)), degree + 1)
    out = np.empty(n)
    ux, inv = np.unique(x, return_inverse=True)
    fit_u = np.empty(len(ux))
    for j, x0 in enumerate(ux):
        d = np.abs(x - x0)
        h = np.partition(d, q - 1)[q - 1]
        if h <= 0:
            h = max(np.min(d[d > 0]) if np.any(d > 0) else 1.0, 1e-12)
        h = h * 1.0000001
        w = np.clip(1 - (d / h) ** 3, 0, None) ** 3
        A = np.vander(x - x0, degree + 1, increasing=True)
        W = w[:, None]
        beta, *_ = np.linalg.lstsq(A * np.sqrt(W), y * np.sqrt(w), rcond=None)
        fit_u[j] = beta[0]
    out[:] = fit_u[inv]
    return out


# ------------------------------------------------------------------------------
# 4. top_table — extraction des résultats pour un contraste donné
# ------------------------------------------------------------------------------
def top_table(fit_eb: dict, coef_idx: int, protein_names=None) -> pd.DataFrame:
    """
    Retourne un DataFrame trié pour le contraste coef_idx.

    Parameters
    ----------
    fit_eb       : dict retourné par ebayes()
    coef_idx     : index du contraste (0-based)
    protein_names: list/array de noms (optionnel)

    Returns
    -------
    DataFrame : logFC, t, P.Value, adj.P.Val
    """
    n = fit_eb["coefficients"].shape[0]
    names = protein_names if protein_names is not None else np.arange(n)

    df = pd.DataFrame({
        "name":      names,
        "logFC":     fit_eb["coefficients"][:, coef_idx],
        "t":         fit_eb["t_stat"][:, coef_idx],
        "P.Value":   fit_eb["p_value"][:, coef_idx],
        "adj.P.Val": fit_eb["p_adj"][:, coef_idx],
    })
    return df


# ------------------------------------------------------------------------------
# 5. make_contrasts — génère la matrice de contrastes (toutes paires)
# ------------------------------------------------------------------------------
DEFAULT_CONTROL_SYNONYMS = [
    "ctl", "ctrl", "control", "wt", "mock", "vehicle", "veh",
    "untreated", "naive", "temoin", "témoin", "dmso",
]


def detect_control_condition(group_names: list, synonyms: "list | None" = None,
                             explicit: "str | None" = None) -> "str | None":
    """
    Détecte automatiquement une condition de type "contrôle" par correspondance
    EXACTE (insensible à la casse) contre une liste de synonymes courants — pas
    de correspondance partielle/substring, pour éviter un faux positif du type
    "CS8" match "ctl" sur un préfixe accidentel.

    Parameters
    ----------
    explicit  : nom exact d'une condition à forcer comme contrôle (prioritaire
                sur la détection automatique) — utile quand le nom réel ne
                matche aucun synonyme standard (ex: "CS8" dans un design où
                cette condition EST le contrôle mais nommée autrement).
    synonyms  : liste de synonymes à tester, défaut DEFAULT_CONTROL_SYNONYMS.

    Returns
    -------
    Le nom du groupe contrôle détecté (str), ou None si :
      - aucune correspondance trouvée
      - plusieurs conditions matchent (ambigu : on ne devine pas, on
        prévient l'appelant de trancher explicitement via `explicit`)
    """
    if explicit:
        if explicit in group_names:
            return explicit
        print(f"  [WARN] control_condition='{explicit}' introuvable parmi "
              f"les conditions {group_names} — ignoré.")
        return None

    syn = set(s.lower() for s in (synonyms or DEFAULT_CONTROL_SYNONYMS))
    matches = [g for g in group_names if g.lower() in syn]

    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        print(f"  [WARN] Plusieurs conditions ressemblent à un contrôle "
              f"{matches} — ambigu, aucune n'est choisie automatiquement. "
              f"Préciser control_condition dans la config pour trancher.")
    return None


def reorder_control_last(group_names: list, control_name: str) -> list:
    """
    Replace control_name en DERNIÈRE position de group_names.

    Dans make_all_contrasts (i<j -> nom "{g1}_vs_{g2}", g1=numérateur), une
    condition à l'index le plus haut est TOUJOURS dénominateur (g2) dans
    chacun de ses contrastes -> c'est la position en dernier qui garantit
    "exp_vs_ctl" plutôt que "ctl_vs_exp" pour toutes les paires impliquant
    le contrôle, pas la première position.
    """
    if control_name not in group_names:
        return group_names
    reordered = [g for g in group_names if g != control_name] + [control_name]
    return reordered


def make_all_contrasts(group_names: list, n_total_cols: "int | None" = None
                       ) -> tuple[np.ndarray, list[str]]:
    """
    Génère toutes les paires de groupes sous forme de matrice de contrastes.

    Parameters
    ----------
    group_names : noms des colonnes de CONDITION (celles sur lesquelles on
                  contraste).
    n_total_cols : si fourni (design apparié ~0 + condition + subject), la
                  matrice est paddée à n_total_cols lignes avec des zéros pour
                  les colonnes nuisance (sujet) placées APRÈS les conditions.
                  Les contrastes ne portent donc que sur les conditions.

    Returns
    -------
    contrast_matrix : np.ndarray (n_rows × n_contrasts)
                      n_rows = n_total_cols si fourni, sinon len(group_names)
    contrast_names  : list de str, ex: ["GroupA_vs_GroupB", ...]
    """
    n = len(group_names)
    n_rows = n_total_cols if n_total_cols is not None else n
    pairs = list(combinations(range(n), 2))
    contrast_matrix = np.zeros((n_rows, len(pairs)))
    contrast_names = []

    for k, (i, j) in enumerate(pairs):
        contrast_matrix[i, k] = 1
        contrast_matrix[j, k] = -1
        g1 = group_names[i].replace(".", "_")
        g2 = group_names[j].replace(".", "_")
        contrast_names.append(f"{g1}_vs_{g2}")

    return contrast_matrix, contrast_names


# ------------------------------------------------------------------------------
# 6. make_design_matrix — matrice design ~0 + condition
# ------------------------------------------------------------------------------
def make_design_matrix(conditions: list,
                       control_condition: "str | None" = None
                       ) -> tuple[np.ndarray, list[str]]:
    """
    Crée la design matrix indicator ~0 + condition (one-hot encoding).

    Parameters
    ----------
    control_condition : si fourni (déjà résolu par detect_control_condition),
                        placé en DERNIÈRE position de group_names — garantit
                        que ce groupe est systématiquement le dénominateur
                        ("exp_vs_ctl") dans tous les contrastes générés par
                        make_all_contrasts (cf. reorder_control_last).

    Returns
    -------
    design      : np.ndarray (n_samples × n_groups)
    group_names : list des noms de groupes (colonnes)
    """
    import re
    # Nettoyage similaire à make.names() de R
    clean = [re.sub(r"[^A-Za-z0-9_]", ".", c) for c in conditions]
    clean = ["X" + c if c[0].isdigit() else c for c in clean]

    # Ordre stable = ordre d'apparition
    seen = []
    for c in clean:
        if c not in seen:
            seen.append(c)
    group_names = seen

    if control_condition is not None:
        clean_control = re.sub(r"[^A-Za-z0-9_]", ".", control_condition)
        if clean_control[0].isdigit():
            clean_control = "X" + clean_control
        group_names = reorder_control_last(group_names, clean_control)

    n_samples = len(clean)
    n_groups = len(group_names)
    design = np.zeros((n_samples, n_groups))
    for i, c in enumerate(clean):
        j = group_names.index(c)
        design[i, j] = 1.0

    return design, group_names


def make_paired_design_matrix(conditions: list, subjects: list,
                              control_condition: "str | None" = None
                              ) -> tuple[np.ndarray, list[str], list[str]]:
    """
    Design matrix pour un plan APPARIÉ : ~0 + condition + subject.

    Le même sujet apparaît dans plusieurs conditions. On ajoute un intercept par
    sujet (effet fixe) pour absorber sa ligne de base individuelle : c'est
    l'équivalent d'un test t apparié généralisé, qui retire la variabilité
    inter-sujets et récupère la puissance statistique.

    Encodage : one-hot des conditions (comme le modèle non apparié) PUIS one-hot
    des sujets en retirant le premier sujet (référence) pour éviter la
    colinéarité avec l'intercept implicite des conditions. Les contrastes ne
    portent QUE sur les colonnes de condition (le sujet est un nuisance factor).

    Parameters
    ----------
    control_condition : cf. make_design_matrix — placé en dernière position
                        de group_names (dénominateur systématique).

    Returns
    -------
    design        : np.ndarray (n_samples × (n_cond + n_subject-1))
    group_names   : noms des colonnes de CONDITION (pour makeContrasts)
    all_col_names : noms de TOUTES les colonnes (condition + subject)
    """
    import re
    # --- Colonnes de condition (identique à make_design_matrix) ---
    clean_c = [re.sub(r"[^A-Za-z0-9_]", ".", c) for c in conditions]
    clean_c = ["X" + c if c[0].isdigit() else c for c in clean_c]
    group_names = []
    for c in clean_c:
        if c not in group_names:
            group_names.append(c)

    if control_condition is not None:
        clean_control = re.sub(r"[^A-Za-z0-9_]", ".", control_condition)
        if clean_control[0].isdigit():
            clean_control = "X" + clean_control
        group_names = reorder_control_last(group_names, clean_control)

    # --- Colonnes de sujet (one-hot, 1er sujet = référence, retiré) ---
    clean_s = [re.sub(r"[^A-Za-z0-9_]", ".", str(s)) for s in subjects]
    clean_s = ["S" + s if s[0].isdigit() else s for s in clean_s]
    subj_levels = []
    for s in clean_s:
        if s not in subj_levels:
            subj_levels.append(s)
    subj_ref = subj_levels[0]                    # référence (absorbée)
    subj_cols = [s for s in subj_levels if s != subj_ref]

    n_samples = len(clean_c)
    n_cond = len(group_names)
    n_subj = len(subj_cols)

    design = np.zeros((n_samples, n_cond + n_subj))
    for i in range(n_samples):
        design[i, group_names.index(clean_c[i])] = 1.0     # condition
        if clean_s[i] in subj_cols:                        # subject (non-ref)
            design[i, n_cond + subj_cols.index(clean_s[i])] = 1.0

    all_col_names = group_names + [f"subject_{s}" for s in subj_cols]
    return design, group_names, all_col_names


def add_covariate_columns(design: np.ndarray, all_col_names: list,
                          values: list, prefix: str
                          ) -> tuple[np.ndarray, list]:
    """
    Ajoute une covariable catégorielle nuisance (batch, plaque, ordre
    d'injection catégorisé, etc.) à une design matrix EXISTANTE — même
    principe que le terme 'subject' de make_paired_design_matrix : one-hot
    avec le premier niveau retiré (référence absorbée dans l'intercept
    implicite des conditions) pour préserver le rang plein de la matrice.

    Compose avec un design standard OU apparié : appeler APRÈS
    make_design_matrix / make_paired_design_matrix, puis repasser
    len(all_col_names_étendu) à make_all_contrasts(..., n_total_cols=...)
    pour que les contrastes ne portent que sur les conditions.

    Parameters
    ----------
    design        : design matrix existante (n_samples × n_cols)
    all_col_names : noms de TOUTES les colonnes existantes de `design`
    values        : valeurs de la covariable, une par échantillon (même
                    ordre que les lignes de `design`)
    prefix        : préfixe des nouveaux noms de colonnes (ex. "batch")

    Returns
    -------
    design_étendue, all_col_names_étendus
    """
    import re
    clean = [re.sub(r"[^A-Za-z0-9_]", ".", str(v)) for v in values]
    clean = ["X" + c if c[0].isdigit() else c for c in clean]

    levels = []
    for c in clean:
        if c not in levels:
            levels.append(c)

    if len(levels) < 2:
        print(f"  [WARN] Covariate '{prefix}' has a single level — ignored "
              "(nothing to control for).")
        return design, all_col_names

    ref = levels[0]                          # référence (absorbée)
    extra_levels = [l for l in levels if l != ref]
    n_samples = design.shape[0]
    extra = np.zeros((n_samples, len(extra_levels)))
    for i, c in enumerate(clean):
        if c != ref:
            extra[i, extra_levels.index(c)] = 1.0

    design_ext = np.hstack([design, extra])
    names_ext = list(all_col_names) + [f"{prefix}_{l}" for l in extra_levels]
    return design_ext, names_ext


# ------------------------------------------------------------------------------
# Helpers internes
# ------------------------------------------------------------------------------
def _inverse_trigamma(x: float, max_iter: int = 50) -> float:
    """
    Inverse de la fonction trigamma : retourne y tel que trigamma(y) = x.
    Algorithme de Newton sur 1/y (Smyth, limma::trigammaInverse).
    """
    x = float(x)
    if np.isnan(x):
        return np.nan
    if x > 1e7:
        return 1.0 / np.sqrt(x)
    if x < 1e-6:
        return 1.0 / x
    # Initialisation (approximation de départ de Smyth)
    y = 0.5 + 1.0 / x
    for _ in range(max_iter):
        tri = float(polygamma(1, y))          # trigamma(y)
        tetra = float(polygamma(2, y))        # tetragamma(y) = trigamma'(y)
        # Newton sur g(y)=trigamma(y)-x, mais reparamétré pour stabilité :
        # delta = tri * (1 - tri/x) / tetra
        if tetra == 0:
            break
        delta = tri * (1.0 - tri / x) / tetra
        y += delta
        if -delta / y < 1e-8:
            break
    return max(y, 1e-6)


def _fit_fdist(s2: np.ndarray, df1: float) -> "tuple[float, float]":
    """limma::fitFDist (sans covariable, ddl égaux) : retourne (df2, scale) =
    (d0, s0²). Variances nulles décalées à 1e-5 x médiane comme limma ; si plus
    de la moitié sont nulles, la médiane est remplacée par 1 (limma avertit que
    eBayes n'est alors pas fiable)."""
    x = np.asarray(s2, float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n <= 1:
        return np.nan, np.nan
    x = np.maximum(x, 0.0)
    m = float(np.median(x))
    if m == 0:
        print("  [WARN] More than half of residual variances are exactly zero: "
              "eBayes unreliable")
        m = 1.0
    x = np.maximum(x, 1e-5 * m)
    e = np.log(x) - digamma(df1 / 2.0) + np.log(df1 / 2.0)
    emean = float(np.mean(e))
    evar = float(np.sum((e - emean) ** 2) / (n - 1)) - float(polygamma(1, df1 / 2.0))
    if evar > 0:
        d0 = 2.0 * _inverse_trigamma(evar)
        s0 = float(np.exp(emean + digamma(d0 / 2.0) - np.log(d0 / 2.0)))
    else:
        d0 = np.inf
        s0 = float(np.mean(x))
    return float(d0), s0


def _bh_correction(p_values: np.ndarray) -> np.ndarray:
    """Correction Benjamini-Hochberg (FDR), identique à p.adjust(method = "BH") :
    les NaN restent NaN et ne comptent pas dans le nombre de tests."""
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    n = int(ok.sum())
    if n == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)[::-1]                       # décroissant
    i = np.arange(n, 0, -1)
    adj = np.minimum.accumulate(n / i * pv[order])
    res = np.empty(n)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out
