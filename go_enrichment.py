# ==============================================================================
# go_enrichment.py — Enrichissement GO/gProfiler optionnel (portage EnrichGO.R)
# Version Python — Proteogen
# ==============================================================================
# Portage fidèle de 260416-EnrichGO.R :
#   - gProfiler (gost) via le package officiel gprofiler-official
#   - sources : GO:BP, GO:CC, GO:MF, REAC, KEGG
#   - correction FDR, exclusion des termes racines, term_size 5–1000
#   - z-score = moyenne des LFC des protéines du terme
#   - graphiques : Manhattan, lollipop facetté, dotplot facetté, chord (GO:BP top5)
#   - export Excel : un onglet par contraste (données + plots)
#
# SEUILS : identiques à ceux des volcanos du pipeline principal (params).
# ESPÈCE : demandée à l'utilisateur (ex: hsapiens, rnorvegicus, mmusculus...).
#
# ENTIÈREMENT NON BLOQUANT : toute erreur (réseau, espèce inconnue, pas assez
# de protéines, API indisponible) est capturée et n'interrompt jamais le
# pipeline principal.
# ==============================================================================

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


def _strip_common_affix(conditions):
    """Retire le préfixe ET le suffixe communs à TOUTES les conditions, pour ne
    garder que la partie discriminante. Identique à la fonction du dashboard
    (build_dashboardv7.py) — gardée synchronisée pour un nommage d'onglets et un
    appariement cohérents de part et d'autre.

    Ex: ['X18_009IG01Ctrl','X18_009IG01di6h','X18_009IG01di24h']
        -> {'...Ctrl':'Ctrl', '...di6h':'di6h', '...di24h':'di24h'}

    Garanties : jamais de label vide, labels uniques (suffixe #i si collision).
    """
    conds = list(dict.fromkeys(conditions))
    if len(conds) <= 1:
        return {c: (c[:10] if len(c) > 10 else c) for c in conds}

    pre = os.path.commonprefix(conds)
    suf = os.path.commonprefix([c[::-1] for c in conds])[::-1]

    def _trim(c):
        core = c
        if pre and core.startswith(pre):
            core = core[len(pre):]
        if suf and core.endswith(suf) and len(core) > len(suf):
            core = core[:len(core) - len(suf)]
        core = core.strip(" _-.")
        return core if core else c

    trimmed = {c: _trim(c) for c in conds}
    vals = list(trimmed.values())
    if any(not v for v in vals) or len(set(vals)) < len(vals):
        def _trim_pre_only(c):
            core = c[len(pre):] if (pre and c.startswith(pre)) else c
            core = core.strip(" _-.")
            return core if core else c
        trimmed2 = {c: _trim_pre_only(c) for c in conds}
        if len(set(trimmed2.values())) == len(trimmed2):
            trimmed = trimmed2

    out, seen = {}, {}
    for c, lbl in trimmed.items():
        s = lbl if len(lbl) <= 8 else lbl[:8]
        if s in seen.values():
            i = 2
            while f"{s[:6]}#{i}" in seen.values():
                i += 1
            s = f"{s[:6]}#{i}"
        seen[c] = s
        out[c] = s
    return out

# Sources d'enrichissement (comme le script R)
GO_SOURCES = ["GO:BP", "GO:CC", "GO:MF", "REAC", "KEGG"]

# Termes racines GO à exclure (molecular_function, biological_process, cellular_component)
GO_ROOT_TERMS = {"GO:0003674", "GO:0008150", "GO:0005575"}

# Palette par source pour le Manhattan plot
SOURCE_COLORS = {
    "GO:BP": "#FF7F0E", "GO:CC": "#2CA02C", "GO:MF": "#1F77B4",
    "REAC":  "#9467BD", "KEGG":  "#D62728", "TF": "#8C564B",
    "MIRNA": "#E377C2", "HPA": "#7F7F7F", "CORUM": "#BCBD22", "HP": "#17BECF",
}


# ==============================================================================
# 0. QUESTION INTERACTIVE (optionnelle)
# ==============================================================================

def ask_go_params() -> dict | None:
    """
    Demande si l'utilisateur veut lancer l'enrichissement GO et avec quelle espèce.
    Les SEUILS (p, ratio) ne sont PAS redemandés : ils proviennent des volcanos.
    Retourne None si refusé → le pipeline saute proprement cette étape.
    """
    print("\n" + "="*60)
    print("  GO ENRICHMENT / gProfiler (optional)")
    print("="*60)
    print("  [!] Requires internet access to g:Profiler (biit.cs.ut.ee).")
    print("      The thresholds used are those defined for the volcanos.")
    print("\n  Example gProfiler species codes:")
    print("     hsapiens (human)     | mmusculus (mouse)   | rnorvegicus (rat)")
    print("     drerio (zebrafish)   | scerevisiae (yeast) | ecoli / efaecalis")

    rep = input("\n  Run GO enrichment? (y/N) -> ").strip().lower()
    if rep not in ("o", "oui", "y", "yes"):
        print("  [SKIP] GO enrichment skipped.\n")
        return None

    organism = input("  gProfiler species code (e.g. rnorvegicus) -> ").strip()
    if not organism:
        print("  [WARN] No species provided — step skipped.\n")
        return None

    return {"organism": organism}


# ==============================================================================
# 1. APPEL gProfiler (gost)
# ==============================================================================

def run_gost(query: list[str], organism: str) -> pd.DataFrame | None:
    """
    Appelle gProfiler (équivalent de gost() en R).
    Retourne un DataFrame de résultats ou None en cas d'échec.

    Paramètres (équivalents du gost() R, adaptés à l'API Python) :
      user_threshold=0.05, no_iea=False (= exclude_iea=FALSE),
      significance_threshold_method='fdr' (= correction_method='fdr'),
      domain_scope='known', no_evidences=False (= evcodes=TRUE).
    """
    try:
        from gprofiler import GProfiler
    except ImportError:
        print("    [WARN] Package 'gprofiler-official' not installed "
              "(pip install gprofiler-official). GO step skipped.")
        return None

    try:
        gp = GProfiler(return_dataframe=True)
        res = gp.profile(
            organism=organism,
            query=query,
            sources=GO_SOURCES,
            user_threshold=0.05,
            no_evidences=False,       # evcodes=TRUE → récupère les intersections
            no_iea=False,             # exclude_iea=FALSE
            domain_scope="known",
            significance_threshold_method="fdr",   # correction_method='fdr'
        )
        if res is None or len(res) == 0:
            return None
        return res
    except Exception as e:
        print(f"    [WARN] gProfiler failed ({type(e).__name__}: {str(e)[:80]})")
        return None


# ==============================================================================
# 2. ISOLATION DES PROTÉINES SIGNIFICATIVES PAR CONTRASTE
# ==============================================================================

def isolate_significant(df_comparaison: pd.DataFrame, contrast: str,
                        params: dict) -> tuple[list[str], pd.DataFrame]:
    """
    Réplique l'isolation du script R (idx_sig) :
    - filtre p < seuil ET |LFC| > lfc_min   (seuils = volcanos)
    - tri par p croissante
    - IDs nettoyés (1er ID avant ';'), dédupliqués
    Retourne (liste protéines, df_local des stats du contraste).

    NB : le script R utilise p.val brute (col_p = _p.val). On respecte ici le
    choix volcano de l'utilisateur (p.val ou p.adj) pour rester cohérent avec
    le pipeline principal.
    """
    suffixe_p = "_p.adj" if params["volcano_use_padj"] else "_p.val"
    col_p    = f"{contrast}{suffixe_p}"
    col_diff = f"{contrast}_diff"

    if col_p not in df_comparaison.columns or col_diff not in df_comparaison.columns:
        return [], pd.DataFrame()

    id_col = "name" if "name" in df_comparaison.columns else df_comparaison.columns[0]
    seuil_p = params["volcano_p_thresh"]
    lfc_min = params["volcano_lfc_min"]

    df_local = df_comparaison[[id_col, col_p, col_diff]].copy()
    df_local.columns = ["name", "p_val", "diff"]
    df_local = df_local.dropna(subset=["p_val", "diff"])
    df_local = df_local[(df_local["p_val"] < seuil_p) &
                        (df_local["diff"].abs() > lfc_min)]
    df_local = df_local.sort_values("p_val")

    # Nettoyage IDs (1er avant ';'), déduplication en gardant l'ordre
    df_local["name_clean"] = (df_local["name"].astype(str)
                              .str.split(";").str[0].str.strip())
    prots_clean = df_local["name_clean"].drop_duplicates().tolist()

    return prots_clean, df_local


# ==============================================================================
# 3. POST-TRAITEMENT DU TABLEAU D'ENRICHISSEMENT (z-score, filtres)
# ==============================================================================

def process_enrichment(gost_df: pd.DataFrame, df_local: pd.DataFrame
                        ) -> pd.DataFrame:
    """
    Réplique le calcul tab_enrich du script R :
    - exclusion des termes racines GO
    - filtre term_size 5–1000
    - gene_ratio = intersection_size / query_size
    - z_score = moyenne des LFC des protéines du terme
    - tri par (source, p_value)
    """
    df = gost_df.copy()

    # Harmonisation des noms de colonnes gProfiler Python → noms du script R
    rename = {
        "native": "term_id", "name": "term_name", "p_value": "p_value",
        "term_size": "term_size", "query_size": "query_size",
        "intersection_size": "intersection_size", "source": "source",
        "intersections": "intersection",   # liste des gènes (evcodes)
    }
    for old, new in rename.items():
        if old in df.columns and old != new:
            df = df.rename(columns={old: new})

    # Exclusion termes racines
    if "term_id" in df.columns:
        df = df[~df["term_id"].isin(GO_ROOT_TERMS)]

    # Filtre term_size
    if "term_size" in df.columns:
        df = df[(df["term_size"] >= 5) & (df["term_size"] <= 1000)]

    if df.empty:
        return df

    # gene_ratio
    if "intersection_size" in df.columns and "query_size" in df.columns:
        df["gene_ratio"] = df["intersection_size"] / df["query_size"].replace(0, np.nan)

    # Map nom → LFC pour le z-score (df_local peut être vide ou sans name_clean)
    if "name_clean" in df_local.columns and "diff" in df_local.columns:
        lfc_map = dict(zip(df_local["name_clean"], df_local["diff"]))
    else:
        lfc_map = {}

    def _zscore(intersection):
        # gProfiler renvoie soit une liste, soit une string séparée par ','
        if isinstance(intersection, (list, tuple, np.ndarray)):
            prots = [str(p) for p in intersection]
        elif isinstance(intersection, str):
            prots = [p.strip() for p in intersection.split(",")]
        else:
            return 0.0
        prots = [p.split(";")[0] for p in prots]
        vals = [lfc_map[p] for p in prots if p in lfc_map]
        return float(np.mean(vals)) if vals else 0.0

    if "intersection" in df.columns:
        df["z_score"] = df["intersection"].apply(_zscore)
    else:
        df["intersection"] = ""
        df["z_score"] = 0.0

    # --- Garantir la présence de TOUTES les colonnes utilisées en aval ---
    # (les plots et l'export y accèdent directement ; on évite tout KeyError
    #  si gProfiler n'a pas renvoyé une colonne attendue).
    defaults = {
        "term_id": "", "term_name": "", "source": "GO:BP",
        "p_value": 1.0, "term_size": 0, "query_size": 0,
        "intersection_size": 0, "gene_ratio": 0.0, "z_score": 0.0,
        "intersection": "",
    }
    for col, val in defaults.items():
        if col not in df.columns:
            df[col] = val

    sort_cols = [c for c in ["source", "p_value"] if c in df.columns]
    df = df.sort_values(sort_cols).reset_index(drop=True)
    return df


# ==============================================================================
# 4. GRAPHIQUES
# ==============================================================================

def plot_manhattan(tab: pd.DataFrame, contrast: str, out_dir: str) -> str | None:
    """
    Manhattan plot façon gostplot : -log10(p) par terme, groupé par source.
    Les 2 termes les plus significatifs par source sont annotés.
    """
    if tab.empty or "source" not in tab.columns:
        return None

    df = tab.copy()
    df["logp"] = -np.log10(df["p_value"].clip(lower=1e-300))
    sources = [s for s in GO_SOURCES if s in df["source"].unique()]
    sources += [s for s in df["source"].unique() if s not in sources]

    total_pts = len(df)
    # Largeur adaptative : ni trop étalé (peu de points) ni trop serré
    fig_w = float(np.clip(total_pts * 0.45 + 3, 7, 16))
    fig, ax = plt.subplots(figsize=(fig_w, 6))
    x_offset = 0
    xticks, xlabels = [], []
    texts = []

    for src in sources:
        sub = df[df["source"] == src].sort_values("p_value")
        n = len(sub)
        if n == 0:
            continue
        xs = np.arange(x_offset, x_offset + n)
        color = SOURCE_COLORS.get(src, "#666666")
        ax.scatter(xs, sub["logp"].values, c=color, s=28, alpha=0.75,
                   label=src, edgecolor="none")
        xticks.append(x_offset + n / 2)
        xlabels.append(src)
        # Top 2 par source à annoter (texte court, repoussé ensuite)
        for rank in range(min(2, n)):
            nm = str(sub["term_name"].values[rank])
            nm = nm if len(nm) <= 32 else nm[:30] + "…"
            texts.append(ax.text(xs[rank], sub["logp"].values[rank], nm,
                                 fontsize=6.5))
        x_offset += n + max(1, n // 8)

    # Anti-chevauchement des labels (adjustText si dispo)
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, only_move={"text": "y"},
                    arrowprops=dict(arrowstyle="-", color="grey", lw=0.3))
    except ImportError:
        pass

    ax.axhline(-np.log10(0.05), ls="--", color="grey", lw=0.8)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=10, fontweight="bold")
    ax.set_ylabel("-log10(p-value FDR)")
    ax.set_title(f"Manhattan GO/gProfiler — {contrast.replace('_vs_', ' vs ')}")
    # Marges : éviter que points/labels touchent les bords
    ax.set_xlim(-1, x_offset)
    ymax = df["logp"].max()
    ax.set_ylim(0, ymax * 1.15 + 0.5)
    ax.legend(fontsize=8, loc="upper right", ncol=max(1, len(sources)))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    f = os.path.join(out_dir, f"go_manhattan_{contrast}.png")
    fig.savefig(f, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return f


def _top_terms_per_source(tab: pd.DataFrame, n_per_source: int = 8) -> pd.DataFrame:
    """Top n termes par source, avec term_name tronqué pour l'affichage."""
    df = tab.copy()
    # Garantir les colonnes de base AVANT tout traitement
    if "source" not in df.columns:
        df["source"] = "GO"
    if "p_value" not in df.columns:
        df["p_value"] = 1.0
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce").fillna(1.0)

    # Top n par source via boucle (évite les pièges de groupby.apply selon
    # la version de pandas, qui pouvaient faire disparaître la colonne 'source')
    parts = []
    for src in df["source"].dropna().unique():
        sub = df[df["source"] == src].nsmallest(n_per_source, "p_value")
        parts.append(sub)
    top = (pd.concat(parts, ignore_index=True) if parts
           else df.head(0).copy())

    # Garantir toutes les colonnes utilisées par les plots
    top["minus_log10_p"] = -np.log10(top["p_value"].clip(lower=1e-300))
    if "z_score" not in top.columns:
        top["z_score"] = 0.0
    if "gene_ratio" not in top.columns:
        top["gene_ratio"] = 0.0
    if "intersection_size" not in top.columns:
        top["intersection_size"] = 1
    if "term_name" in top.columns:
        top["term_short"] = top["term_name"].astype(str).apply(
            lambda s: s if len(s) <= 45 else s[:43] + "…")
    else:
        top["term_short"] = [f"term_{i}" for i in range(len(top))]
    return top


def plot_lollipop_faceted(tab: pd.DataFrame, contrast: str, out_dir: str) -> str | None:
    """Lollipop facetté par source, couleur = z-score (LFC moyen)."""
    if tab.empty:
        return None
    top = _top_terms_per_source(tab, 8)
    sources = [s for s in GO_SOURCES if s in top["source"].unique()]
    sources += [s for s in top["source"].unique() if s not in sources]
    n_src = len(sources)
    if n_src == 0:
        return None

    n_src = len(sources)
    if n_src == 0:
        return None

    # Layout vertical (1 colonne) : évite tout chevauchement des labels longs.
    fig, axes = plt.subplots(n_src, 1,
                             figsize=(9, max(2.4, 2.2 * n_src)),
                             squeeze=False)
    axes = axes.flatten()

    zmax = max(abs(top["z_score"]).max(), 1e-6)
    for idx, src in enumerate(sources):
        ax = axes[idx]
        sub = top[top["source"] == src].sort_values("minus_log10_p")
        ypos = np.arange(len(sub))
        ax.hlines(ypos, 0, sub["minus_log10_p"], color="grey", lw=1)
        sc = ax.scatter(sub["minus_log10_p"], ypos, c=sub["z_score"],
                        cmap="RdBu_r", vmin=-zmax, vmax=zmax, s=90,
                        edgecolor="black", lw=0.5, zorder=3)
        ax.set_yticks(ypos)
        ax.set_yticklabels(sub["term_short"], fontsize=8)
        ax.set_xlabel("-log10(p)", fontsize=9)
        ax.set_title(src, fontweight="bold", fontsize=11, loc="left")
        ax.set_xlim(0, sub["minus_log10_p"].max() * 1.15 + 0.3)
        ax.margins(y=0.12)
        ax.spines[["top", "right"]].set_visible(False)

    for idx in range(n_src, len(axes)):
        axes[idx].set_visible(False)

    sm = plt.cm.ScalarMappable(cmap="RdBu_r",
                               norm=plt.Normalize(vmin=-zmax, vmax=zmax))
    cbar = fig.colorbar(sm, ax=axes.tolist(), fraction=0.025, pad=0.04)
    cbar.set_label("Z-score (LFC moyen)", fontsize=9)
    fig.suptitle(f"Lollipop GO — {contrast.replace('_vs_', ' vs ')}",
                 fontsize=13, fontweight="bold", y=1.0)
    fig.subplots_adjust(left=0.32, right=0.86, hspace=0.55, top=0.90)
    f = os.path.join(out_dir, f"go_lollipop_{contrast}.png")
    fig.savefig(f, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return f


def plot_dotplot_faceted(tab: pd.DataFrame, contrast: str, out_dir: str) -> str | None:
    """Dotplot facetté : gene_ratio en x, taille = intersection_size, couleur = p."""
    if tab.empty:
        return None
    top = _top_terms_per_source(tab, 8)
    sources = [s for s in GO_SOURCES if s in top["source"].unique()]
    sources += [s for s in top["source"].unique() if s not in sources]
    n_src = len(sources)
    if n_src == 0:
        return None

    n_src = len(sources)
    if n_src == 0:
        return None

    # Layout vertical (1 colonne) : les noms de termes GO sont longs, empiler
    # les sources évite tout chevauchement horizontal entre facettes.
    fig, axes = plt.subplots(n_src, 1,
                             figsize=(9, max(2.4, 2.2 * n_src)),
                             squeeze=False)
    axes = axes.flatten()

    # Échelle de couleur commune à toutes les facettes
    all_logp = -np.log10(top["p_value"].clip(lower=1e-300))
    vmin, vmax = float(all_logp.min()), float(all_logp.max())
    sc = None
    for idx, src in enumerate(sources):
        ax = axes[idx]
        sub = top[top["source"] == src].sort_values("gene_ratio")
        ypos = np.arange(len(sub))
        sizes = sub["intersection_size"].values * 28 if "intersection_size" in sub.columns else 60
        sc = ax.scatter(sub["gene_ratio"], ypos, s=sizes,
                        c=-np.log10(sub["p_value"].clip(lower=1e-300)),
                        cmap="plasma", vmin=vmin, vmax=vmax,
                        edgecolor="black", lw=0.4, zorder=3)
        ax.set_yticks(ypos)
        ax.set_yticklabels(sub["term_short"], fontsize=8)
        ax.set_xlabel("Gene Ratio", fontsize=9)
        ax.set_title(src, fontweight="bold", fontsize=11, loc="left")
        ax.margins(x=0.18, y=0.12)
        ax.grid(axis="x", ls=":", color="grey", alpha=0.4)
        ax.spines[["top", "right"]].set_visible(False)

    for idx in range(n_src, len(axes)):
        axes[idx].set_visible(False)

    # Colorbar unique partagée à droite, bien détachée des panneaux
    if sc is not None:
        cbar = fig.colorbar(sc, ax=axes.tolist(), fraction=0.025, pad=0.04)
        cbar.set_label("-log10(p)", fontsize=9)
    # Légende des tailles (Count) sur la 1re facette
    if "intersection_size" in top.columns:
        for cnt in sorted(top["intersection_size"].unique())[:5]:
            axes[0].scatter([], [], s=cnt * 28, c="grey",
                            edgecolor="black", lw=0.4, label=str(int(cnt)))
        axes[0].legend(title="Count", loc="upper left",
                       bbox_to_anchor=(1.01, 1.0), fontsize=7,
                       title_fontsize=8, framealpha=0.9)

    fig.suptitle(f"Dotplot GO — {contrast.replace('_vs_', ' vs ')}",
                 fontsize=13, fontweight="bold", y=1.0)
    fig.subplots_adjust(left=0.32, right=0.86, hspace=0.55, top=0.90)
    f = os.path.join(out_dir, f"go_dotplot_{contrast}.png")
    fig.savefig(f, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return f


def plot_chord(tab: pd.DataFrame, df_local: pd.DataFrame,
               contrast: str, out_dir: str) -> str | None:
    """
    Chord plot (GO:BP top 5) version matplotlib : matrice protéines × termes,
    protéines ordonnées par LFC, barre LFC à gauche (GOChord-like).
    """
    if tab.empty:
        return None
    bp = tab[tab["source"] == "GO:BP"].nsmallest(5, "p_value")
    if bp.empty:
        return None

    # Récupérer les protéines de chaque terme
    def _prots(intersection):
        if isinstance(intersection, (list, tuple, np.ndarray)):
            return [str(p).split(";")[0] for p in intersection]
        elif isinstance(intersection, str):
            return [p.strip().split(";")[0] for p in intersection.split(",")]
        return []

    term_prots = {row["term_name"]: _prots(row["intersection"])
                  for _, row in bp.iterrows() if "intersection" in bp.columns}
    all_prots = sorted(set().union(*term_prots.values())) if term_prots else []
    lfc_map = dict(zip(df_local["name_clean"], df_local["diff"]))
    all_prots = [p for p in all_prots if p in lfc_map]
    if len(all_prots) <= 2:
        return None

    terms = list(term_prots.keys())
    # Matrice : NaN si gène absent du terme, sinon LFC du gène (pour la couleur)
    mat = np.full((len(all_prots), len(terms)), np.nan)
    for j, t in enumerate(terms):
        for p in term_prots[t]:
            if p in all_prots:
                mat[all_prots.index(p), j] = lfc_map[p]

    lfc_vals = np.array([lfc_map[p] for p in all_prots])
    order = np.argsort(lfc_vals)
    mat = mat[order]
    prots_ordered = [all_prots[i] for i in order]
    lfc_ordered = lfc_vals[order]

    fig, (ax_lfc, ax_mat) = plt.subplots(
        1, 2, figsize=(11, max(6, len(all_prots) * 0.26)),
        gridspec_kw={"width_ratios": [0.16, 1]})

    vmax = max(np.abs(lfc_ordered).max(), 0.5) if len(lfc_ordered) else 1
    # Barre LFC à gauche
    ax_lfc.barh(range(len(prots_ordered)), lfc_ordered,
                color=plt.cm.RdBu_r((lfc_ordered / (2 * vmax)) + 0.5),
                edgecolor="grey", linewidth=0.3)
    ax_lfc.set_yticks(range(len(prots_ordered)))
    ax_lfc.set_yticklabels(prots_ordered, fontsize=6)
    ax_lfc.set_xlabel("logFC", fontsize=8)
    ax_lfc.axvline(0, color="black", lw=0.5)
    ax_lfc.set_ylim(-0.5, len(prots_ordered) - 0.5)
    ax_lfc.invert_yaxis()
    ax_lfc.spines[["top", "right"]].set_visible(False)

    # Matrice colorée par LFC (cases vides = gène hors du terme)
    masked = np.ma.masked_invalid(mat)
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color="#f0f0f0")   # gris clair pour les cases vides
    im = ax_mat.imshow(masked, aspect="auto", cmap=cmap,
                       vmin=-vmax, vmax=vmax, interpolation="nearest")
    # Grille fine pour séparer les cases
    ax_mat.set_xticks(np.arange(-0.5, len(terms), 1), minor=True)
    ax_mat.set_yticks(np.arange(-0.5, len(prots_ordered), 1), minor=True)
    ax_mat.grid(which="minor", color="white", linewidth=1.2)
    ax_mat.set_xticks(range(len(terms)))
    wrapped = ["\n".join([t[k:k+22] for k in range(0, min(len(t), 44), 22)])
               for t in terms]
    ax_mat.set_xticklabels(wrapped, rotation=40, ha="right", fontsize=7)
    ax_mat.set_yticks([])
    ax_mat.set_title(f"Genes <-> GO:BP terms (top 5) — "
                     f"{contrast.replace('_vs_', ' vs ')}\n"
                     f"color = gene log2FC", fontsize=9)
    cbar = fig.colorbar(im, ax=ax_mat, fraction=0.025, pad=0.02)
    cbar.set_label("log2FC", fontsize=8)

    fig.tight_layout()
    f = os.path.join(out_dir, f"go_chord_{contrast}.png")
    fig.savefig(f, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return f


# ==============================================================================
# 5. POINT D'ENTRÉE PRINCIPAL
# ==============================================================================

def run_go_enrichment(df_comparaison: pd.DataFrame, contrast_names: list[str],
                      params: dict, go_params: dict, out_dir: str) -> dict | None:
    """
    Lance l'enrichissement GO sur tous les contrastes.
    Entièrement protégé : retourne None si tout échoue, sans jamais lever
    d'exception qui interromprait le pipeline principal.

    Seuils de significativité = ceux des volcanos (params).
    Espèce = go_params['organism'].

    Retourne {contrast: {"table": df, "manhattan":.., "lollipop":.., "dotplot":.., "chord":..}}
    """
    if go_params is None:
        return None

    try:
        organism = go_params["organism"]
        print(f"\n[GO] GO/gProfiler enrichment — species '{organism}'")
        print(f"   Seuils (volcanos) : "
              f"{'p.adj' if params['volcano_use_padj'] else 'p.val'} < "
              f"{params['volcano_p_thresh']} | ratio ≥ {params['volcano_ratio_min']}")

        results = {}
        for contrast in contrast_names:
            try:
                prots, df_local = isolate_significant(df_comparaison, contrast, params)

                if len(prots) <= 5:
                    print(f"  [SKIP] {contrast}: too few proteins ({len(prots)}), skipped.")
                    continue

                print(f"  [GO] {contrast}: {len(prots)} proteins -> gProfiler...")
                gost_df = run_gost(prots, organism)
                if gost_df is None or gost_df.empty:
                    print(f"     Aucun enrichissement significatif.")
                    continue

                tab = process_enrichment(gost_df, df_local)
                if tab.empty:
                    print(f"     No term after filtering (term_size, roots).")
                    continue

                # Chaque plot est isolé : un échec graphique ne doit pas faire
                # perdre le tableau d'enrichissement (l'essentiel pour l'Excel).
                def _safe_plot(fn, *a):
                    try:
                        return fn(*a)
                    except Exception as ex:
                        import traceback as _tb
                        detail = f"clé {ex}" if isinstance(ex, KeyError) else str(ex)[:80]
                        print(f"        (plot {fn.__name__} skipped: "
                              f"{type(ex).__name__} — {detail})")
                        line = _tb.format_exc().strip().splitlines()[-2:]
                        for l in line:
                            print(f"          {l.strip()}")
                        return None
                f_lol = _safe_plot(plot_lollipop_faceted, tab, contrast, out_dir)
                f_dot = _safe_plot(plot_dotplot_faceted, tab, contrast, out_dir)

                # Nettoyer colonnes complexes (listes) avant export Excel
                tab_export = tab.copy()
                if "intersection" in tab_export.columns:
                    tab_export["intersection"] = tab_export["intersection"].apply(
                        lambda x: ",".join(map(str, x)) if isinstance(x, (list, tuple, np.ndarray))
                        else str(x))
                # Retirer colonnes d'affichage internes
                tab_export = tab_export.drop(columns=["term_short", "minus_log10_p"],
                                             errors="ignore")

                results[contrast] = {
                    "table": tab_export,
                    "lollipop": f_lol, "dotplot": f_dot,
                }
                n_bp = (tab["source"] == "GO:BP").sum()
                print(f"     [OK] {len(tab)} enriched terms ({n_bp} GO:BP).")

            except Exception as e:
                import traceback
                detail = str(e)
                if isinstance(e, KeyError):
                    detail = f"missing column/key: {e}"
                print(f"     [WARN] Error on {contrast}: {type(e).__name__} "
                      f"({detail}) — skipped.")
                # Trace courte pour diagnostic (1ère ligne utile)
                tb = traceback.format_exc().strip().splitlines()
                for line in tb[-3:]:
                    print(f"        {line.strip()}")
                continue

        if not results:
            print("  [INFO] No significant GO enrichment across all contrasts.")
            return None

        return results

    except Exception as e:
        print(f"  [WARN] GO enrichment interrupted ({type(e).__name__}) — "
              f"the main pipeline is not affected.")
        return None


def run_go_enrichment_clusters(cluster_mapping, params: dict, go_params: dict,
                               out_dir: str) -> dict | None:
    """
    Enrichissement GO par cluster de la heatmap ANOVA.

    Contrairement à run_go_enrichment (qui part des DEP par contraste), on prend
    ici les protéines regroupées par profil d'expression (clusters k-means de la
    heatmap). Background = protéome connu de l'espèce (domain_scope='known').

    Parameters
    ----------
    cluster_mapping : liste de tuples (protein_name, "Cluster_N") — telle que
                      produite par run_anova_heatmaps.
    params, go_params, out_dir : identiques à run_go_enrichment.

    Retourne {cluster_id: {"table": df, "lollipop": path, "dotplot": path}}.
    Entièrement protégé : n'interrompt jamais le pipeline.

    Note z-score : les protéines d'un cluster proviennent de contrastes variés,
    il n'existe pas de LFC unique par protéine au niveau cluster. Le z_score des
    termes est donc neutralisé (0) ici — l'information portée par le cluster est
    le profil d'expression, pas un sens de régulation par contraste.
    """
    if go_params is None:
        return None

    try:
        organism = go_params["organism"]
        # Regrouper les protéines par cluster
        from collections import defaultdict
        clusters = defaultdict(list)
        for prot, cid in cluster_mapping:
            if cid:  # ignorer les non-assignés ("")
                clusters[cid].append(str(prot).split(";")[0].strip())

        if not clusters:
            return None

        print(f"\n[GO] Cluster GO enrichment — species '{organism}' "
              f"({len(clusters)} clusters)")

        results = {}
        for cid in sorted(clusters.keys()):
            prots = list(dict.fromkeys(clusters[cid]))  # dédup, ordre conservé
            try:
                if len(prots) <= 5:
                    print(f"  [SKIP] {cid}: too few proteins ({len(prots)}), skipped.")
                    continue

                print(f"  [GO] {cid}: {len(prots)} proteins -> gProfiler...")
                gost_df = run_gost(prots, organism)
                if gost_df is None or gost_df.empty:
                    print(f"     No significant enrichment.")
                    continue

                # df_local vide → process_enrichment neutralise le z-score (0.0)
                tab = process_enrichment(gost_df, pd.DataFrame())
                if tab.empty:
                    print(f"     No term after filtering (term_size, roots).")
                    continue

                def _safe_plot(fn, *a):
                    try:
                        return fn(*a)
                    except Exception as ex:
                        print(f"        (plot {fn.__name__} skipped: "
                              f"{type(ex).__name__})")
                        return None
                # On réutilise les plots existants ; l'étiquette = cid
                f_lol = _safe_plot(plot_lollipop_faceted, tab, cid, out_dir)
                f_dot = _safe_plot(plot_dotplot_faceted, tab, cid, out_dir)

                tab_export = tab.copy()
                if "intersection" in tab_export.columns:
                    tab_export["intersection"] = tab_export["intersection"].apply(
                        lambda x: ",".join(map(str, x))
                        if isinstance(x, (list, tuple, np.ndarray)) else str(x))
                tab_export = tab_export.drop(columns=["term_short", "minus_log10_p"],
                                             errors="ignore")

                results[cid] = {"table": tab_export,
                                "lollipop": f_lol, "dotplot": f_dot}
                n_bp = (tab["source"] == "GO:BP").sum()
                print(f"     [OK] {len(tab)} enriched terms ({n_bp} GO:BP).")

            except Exception as e:
                detail = f"missing column/key: {e}" if isinstance(e, KeyError) else str(e)[:80]
                print(f"     [WARN] Error on {cid}: {type(e).__name__} ({detail}) — skipped.")
                continue

        if not results:
            print("  [INFO] No significant GO enrichment across all clusters.")
            return None
        return results

    except Exception as e:
        print(f"  [WARN] Cluster GO enrichment interrupted ({type(e).__name__}) — "
              f"the main pipeline is not affected.")
        return None


def export_go_cluster_sheets(wb, go_cluster_results: dict, write_df_fn, ins_img_fn):
    """Ajoute un onglet Cluster_Enrich_N par cluster au classeur Excel.

    IMPORTANT — nommage : on n'utilise PAS le préfixe 'GO_' ici. Le dashboard
    (build_dashboardv7.py) capte toute feuille dont le nom contient 'go' pour
    construire sa heatmap GO cross-contrastes, puis apparie de force chaque
    feuille à un contraste (sans seuil de rejet). Des onglets de CLUSTER (qui ne
    correspondent à aucun contraste) y seraient appariés à tort et casseraient la
    heatmap. Le préfixe 'Cluster_Enrich_' échappe à tous les filtres du dashboard
    ('go' in name / startswith('GO_')), donc ces onglets restent visibles dans
    l'Excel sans perturber le dashboard.
    """
    if not go_cluster_results:
        return
    for cid, data in go_cluster_results.items():
        # cid = "Cluster_1" → onglet "Cluster_Enrich_1"
        n = cid.replace("Cluster_", "")
        sheet_name = f"Cluster_Enrich_{n}"[:31]
        ws = wb.create_sheet(sheet_name)
        df = data["table"]
        start_col = 21
        for j, col in enumerate(df.columns):
            ws.cell(row=1, column=start_col + j, value=str(col))
        for i, row in enumerate(df.itertuples(index=False), start=2):
            for j, val in enumerate(row):
                if isinstance(val, float) and np.isnan(val):
                    val = None
                elif isinstance(val, (list, tuple, np.ndarray)):
                    val = ",".join(map(str, val))
                ws.cell(row=i, column=start_col + j, value=val)
        ins_img_fn(ws, data.get("lollipop"), "A2")
        ins_img_fn(ws, data.get("dotplot"),  "A60")


# ==============================================================================
# 6. EXPORT EXCEL (onglets GO ajoutés au classeur principal)
# ==============================================================================

def _make_go_sheet_names(contrasts):
    """Construit des noms d'onglets GO courts, lisibles et GARANTIS uniques
    (<= 31 car., contrainte Excel).

    Problème résolu : avec de longs préfixes communs aux conditions
    (ex: 'X18_009IG01Ctrl_vs_X18_009IG01di6h'), une troncature naïve à 31 car.
    garde le préfixe commun et produit des noms IDENTIQUES pour deux contrastes
    différents -> Excel corrompt /xl/workbook.xml ('Réparations...').

    Stratégie : retirer le préfixe/suffixe commun aux conditions (partie
    discriminante), recomposer 'GO_<a>v<b>', tronquer, et dédupliquer avec un
    suffixe numérique si nécessaire.

    Retourne {contrast: sheet_name}.
    """
    # Conditions = les 2 côtés de chaque contraste
    all_conds = []
    for c in contrasts:
        all_conds.extend(c.split("_vs_"))
    cond_label = _strip_common_affix(all_conds)

    names = {}
    used = set()
    for c in contrasts:
        parts = c.split("_vs_")
        a = cond_label.get(parts[0], parts[0])[:12] if parts else c[:12]
        b = cond_label.get(parts[1], parts[1])[:12] if len(parts) > 1 else ""
        base = f"GO_{a}v{b}" if b else f"GO_{a}"
        base = base[:31]
        name = base
        i = 2
        # Dédup anti-collision (réserve 3 car. pour le suffixe #NN)
        while name in used:
            suffix = f"#{i}"
            name = (base[:31 - len(suffix)]) + suffix
            i += 1
        used.add(name)
        names[c] = name
    return names


def export_go_sheets(wb, go_results: dict, write_df_fn, ins_img_fn):
    """
    Ajoute un onglet par contraste au classeur Excel existant.
    Mise en page proche du script R : données en colonne U (21),
    plots empilés à gauche (Manhattan, lollipop, dotplot, chord).
    """
    if not go_results:
        return

    sheet_names = _make_go_sheet_names(list(go_results.keys()))

    for contrast, data in go_results.items():
        sheet_name = sheet_names[contrast]
        ws = wb.create_sheet(sheet_name)

        # Données décalées en colonne 21 (U) comme le script R
        df = data["table"]
        start_col = 21
        for j, col in enumerate(df.columns):
            ws.cell(row=1, column=start_col + j, value=str(col))
        for i, row in enumerate(df.itertuples(index=False), start=2):
            for j, val in enumerate(row):
                if isinstance(val, float) and np.isnan(val):
                    val = None
                elif isinstance(val, (list, tuple, np.ndarray)):
                    val = ",".join(map(str, val))
                ws.cell(row=i, column=start_col + j, value=val)

        # Plots empilés à gauche : lollipop (z-score) puis dotplot (gene ratio)
        # Plots empilés à gauche : lollipop (z-score) puis dotplot (gene ratio).
        # Layout vertical = images plus hautes → espacement vertical large.
        ins_img_fn(ws, data.get("lollipop"),  "A2")
        ins_img_fn(ws, data.get("dotplot"),   "A60")
