# ==============================================================================
# gsea_enrichment.py — GSEA rank-based (pré-ranked) optionnel pour Deimos
# Version Python — Proteogen
# ==============================================================================
# Complément de l'ORA (go_enrichment.py). Là où l'ORA part de la LISTE des DEP
# (coupée au seuil), le GSEA classe TOUTES les protéines par un score signé et
# détecte les gene sets concentrés en haut (up) ou en bas (down) du classement —
# captant les effets coordonnés que l'ORA rate (beaucoup de protéines qui
# bougent un peu dans le même sens, sans passer le seuil).
#
# Métrique de classement : -log10(p) * sign(LFC)
#   -> sépare mieux les extrêmes que le LFC seul (combine sens ET significativité)
#
# Moteur : gseapy.prerank (algorithme Broad, permutations). L'algorithme tourne
# en local ; seul le téléchargement des gene sets (Enrichr) nécessite le réseau.
#
# ENTIÈREMENT NON BLOQUANT : toute erreur (gseapy absent, réseau, mapping vide,
# pas assez de protéines) est capturée et n'interrompt jamais le pipeline.
# ==============================================================================

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# Gene sets Enrichr par défaut (couvrent GO + pathways, cohérent avec l'ORA).
# Modifiable via go_params['gsea_gene_sets'] si besoin.
DEFAULT_GENE_SETS = [
    "GO_Biological_Process_2021",
    "GO_Molecular_Function_2021",
    "GO_Cellular_Component_2021",
    "KEGG_2021_Human",
    "Reactome_2022",
]


def _ranking_metric(df_comparaison, contrast, params):
    """Construit le vecteur de ranking signé pour un contraste.

    Métrique = -log10(p) * sign(LFC), sur TOUTES les protéines (pas seulement
    les DEP). Utilise le choix volcano (p.val ou p.adj) pour rester cohérent
    avec le pipeline. Mapping accession -> symbole de gène via la colonne 'Genes'.

    Retourne une pd.Series {gene_symbol: score} triée décroissante, ou None.
    """
    suffixe_p = "_p.adj" if params.get("volcano_use_padj") else "_p.val"
    col_p    = f"{contrast}{suffixe_p}"
    col_diff = f"{contrast}_diff"
    if col_p not in df_comparaison.columns or col_diff not in df_comparaison.columns:
        return None

    gene_col = "Genes" if "Genes" in df_comparaison.columns else None
    if gene_col is None:
        return None

    df = df_comparaison[[gene_col, col_p, col_diff]].copy()
    df.columns = ["gene", "p", "lfc"]
    df = df.dropna(subset=["p", "lfc"])
    # Symbole propre : premier gène avant ';', majuscule, sans espaces
    df["gene"] = (df["gene"].astype(str).str.split(";").str[0]
                  .str.strip().str.upper())
    df = df[df["gene"].notna() & (df["gene"] != "") & (df["gene"] != "NAN")]
    if df.empty:
        return None

    # Métrique signée : -log10(p) * sign(LFC)
    p_clip = df["p"].clip(lower=1e-300)
    df["score"] = -np.log10(p_clip) * np.sign(df["lfc"])
    # En cas de doublons de symboles : garder le score de plus grande amplitude
    df["absscore"] = df["score"].abs()
    df = df.sort_values("absscore", ascending=False).drop_duplicates("gene")
    rnk = df.set_index("gene")["score"].sort_values(ascending=False)
    # Retirer les scores nuls (p=1 exactement) qui n'apportent rien
    rnk = rnk[rnk != 0]
    return rnk if len(rnk) >= 15 else None


def _plot_gsea_nes(res2d, contrast, out_dir, top_n=15):
    """Barplot horizontal des top voies par |NES| (rouge=up, bleu=down)."""
    try:
        df = res2d.copy()
        # Colonnes gseapy : 'Term', 'NES', 'FDR q-val'
        df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
        df["FDR"] = pd.to_numeric(df.get("FDR q-val", np.nan), errors="coerce")
        df = df.dropna(subset=["NES"])
        if df.empty:
            return None
        df["absNES"] = df["NES"].abs()
        df = df.sort_values("absNES", ascending=False).head(top_n)
        df = df.sort_values("NES")
        colors = ["#3498db" if v < 0 else "#e74c3c" for v in df["NES"]]
        labels = [t[:45] + ("…" if len(t) > 45 else "") for t in df["Term"]]
        fig, ax = plt.subplots(figsize=(9, max(3, 0.42 * len(df))))
        ax.barh(range(len(df)), df["NES"], color=colors, edgecolor="black",
                linewidth=0.4)
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_xlabel("NES (Normalized Enrichment Score)")
        ax.set_title(f"GSEA — {contrast}\n(red = enriched in up / blue = in down)")
        fig.tight_layout()
        f = os.path.join(out_dir, f"gsea_nes_{contrast}.png")
        fig.savefig(f, dpi=150)
        plt.close(fig)
        return f
    except Exception:
        return None


def run_gsea(df_comparaison, contrast_names, params, go_params, out_dir):
    """
    GSEA pré-ranked par contraste.

    Retourne {contrast: {"table": df, "nes_plot": path}} ou None.
    Entièrement protégé : n'interrompt jamais le pipeline.
    """
    if go_params is None or not go_params.get("run_gsea", False):
        return None

    try:
        import gseapy
    except ImportError:
        print("  [WARN] GSEA skipped: gseapy not installed "
              "(pip install gseapy).")
        return None

    gene_sets = go_params.get("gsea_gene_sets", DEFAULT_GENE_SETS)
    fdr_max   = go_params.get("gsea_fdr", 0.25)
    n_perm    = go_params.get("gsea_permutations", 1000)

    print(f"\n[GSEA] Pre-ranked GSEA — metric: -log10(p) x sign(LFC)")
    results = {}
    for contrast in contrast_names:
        try:
            rnk = _ranking_metric(df_comparaison, contrast, params)
            if rnk is None:
                print(f"  [SKIP] {contrast}: not enough ranked genes.")
                continue

            print(f"  [GSEA] {contrast}: {len(rnk)} ranked genes -> prerank...")
            pre = gseapy.prerank(
                rnk=rnk, gene_sets=gene_sets,
                min_size=10, max_size=1000,
                permutation_num=n_perm, seed=42,
                no_plot=True, verbose=False, threads=4,
                outdir=None)

            res = pre.res2d.copy()
            # Filtrer par FDR et trier
            res["FDR q-val"] = pd.to_numeric(res["FDR q-val"], errors="coerce")
            res["NES"] = pd.to_numeric(res["NES"], errors="coerce")
            res_sig = res[res["FDR q-val"] <= fdr_max].copy()
            res_sig = res_sig.sort_values("NES", ascending=False)

            if res_sig.empty:
                print(f"     No pathway below FDR {fdr_max}.")
                continue

            f_nes = _plot_gsea_nes(res_sig, contrast, out_dir)
            # Nettoyer les colonnes pour l'export (listes -> texte)
            for col in res_sig.columns:
                if res_sig[col].apply(lambda x: isinstance(x, (list, tuple))).any():
                    res_sig[col] = res_sig[col].apply(
                        lambda x: ";".join(map(str, x))
                        if isinstance(x, (list, tuple)) else x)

            results[contrast] = {"table": res_sig, "nes_plot": f_nes}
            n_up = (res_sig["NES"] > 0).sum()
            n_dn = (res_sig["NES"] < 0).sum()
            print(f"     [OK] {len(res_sig)} pathways (FDR<={fdr_max}): "
                  f"{n_up} up, {n_dn} down.")

        except Exception as e:
            print(f"     [WARN] GSEA error on {contrast}: "
                  f"{type(e).__name__} ({str(e)[:70]}) — skipped.")
            continue

    if not results:
        print("  [INFO] No significant GSEA result across contrasts.")
        return None
    return results


def export_gsea_sheets(wb, gsea_results, ins_img_fn):
    """Ajoute un onglet GSEA_<contrast> par contraste au classeur Excel.

    Nommage sans 'GO' ni 'go' pour ne pas être capté par les filtres GO du
    dashboard (heatmap/réseau GO cross-contrastes).
    """
    if not gsea_results:
        return
    for contrast, data in gsea_results.items():
        # Nom court unique, sans collision avec les onglets GO
        base = f"GSEA_{contrast.replace('_vs_', 'v').replace('_', '')}"
        sheet_name = base[:31]
        ws = wb.create_sheet(sheet_name)
        df = data["table"]
        start_col = 16  # laisse de la place au plot NES à gauche
        for j, col in enumerate(df.columns):
            ws.cell(row=1, column=start_col + j, value=str(col))
        for i, row in enumerate(df.itertuples(index=False), start=2):
            for j, val in enumerate(row):
                if isinstance(val, float) and np.isnan(val):
                    val = None
                elif isinstance(val, (list, tuple, np.ndarray)):
                    val = ";".join(map(str, val))
                ws.cell(row=i, column=start_col + j, value=val)
        ins_img_fn(ws, data.get("nes_plot"), "A2")
