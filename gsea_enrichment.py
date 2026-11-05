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
# Jeux humains : défaut de gsea_organism: human (comportement d'origine).
DEFAULT_GENE_SETS = [
    "GO_Biological_Process_2021",
    "GO_Molecular_Function_2021",
    "GO_Cellular_Component_2021",
    "KEGG_2021_Human",
    "Reactome_2022",
]

# gsea_organism: mouse — collections MSigDB souris (symboles murins). Pas de KEGG
# dans MSigDB souris : Reactome + WikiPathways pour les voies.
DEFAULT_GENE_SETS_MOUSE = [
    "msigdb:mh.all",
    "msigdb:m2.cp.reactome",
    "msigdb:m2.cp.wikipathways",
    "msigdb:m5.go.bp",
    "msigdb:m5.go.mf",
    "msigdb:m5.go.cc",
]
DEFAULT_MSIGDB_VERSION = {"human": "2025.1.Hs", "mouse": "2025.1.Mm"}
GSEA_ORGANISMS = ("human", "mouse")


def gsea_settings(params, go_params=None) -> dict:
    """Réglages GSEA lus dans la configuration (clés gsea_*), avec repli sur
    go_params (ancien passage par go_params['gsea_gene_sets'])."""
    gp = go_params or {}
    org = str(params.get("gsea_organism") or gp.get("gsea_organism") or "human").lower()
    if org not in GSEA_ORGANISMS:
        org = "human"
    sets = params.get("gsea_gene_sets") or gp.get("gsea_gene_sets")
    if isinstance(sets, str):
        sets = [sets]
    if not sets:
        sets = list(DEFAULT_GENE_SETS_MOUSE if org == "mouse" else DEFAULT_GENE_SETS)
    return {
        "organism": org,
        "gene_sets": list(sets),
        "msigdb_version": str(params.get("gsea_msigdb_version")
                              or DEFAULT_MSIGDB_VERSION[org]),
        "fdr": float(params.get("gsea_fdr", gp.get("gsea_fdr", 0.25))),
        "permutations": int(params.get("gsea_permutations",
                                       gp.get("gsea_permutations", 1000))),
        "rank_p": str(params.get("gsea_rank_p", "pvalue")).lower(),
    }


def read_gmt(path: str) -> dict:
    """Fichier GMT (nom <tab> description <tab> gènes...) -> {nom: [gènes]}."""
    out = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            row = line.rstrip("\n\r").split("\t")
            if len(row) >= 3 and row[0]:
                out[row[0]] = [g for g in row[2:] if g]
    return out


def load_gene_sets(specs, organism: str = "human",
                   msigdb_version: "str | None" = None, verbose: bool = True) -> dict:
    """Charge les jeux de gènes UNE fois pour tous les contrastes.

    specs : liste de
      - chemins de fichiers .gmt locaux (version figée, sans réseau) ;
      - "msigdb:<collection>" (ex. msigdb:m5.go.bp, msigdb:mh.all), version
        msigdb_version (ex. 2025.1.Mm) ;
      - noms de bibliothèques Enrichr (ex. GO_Biological_Process_2021).
    Retourne {"<source>__<terme>": [SYMBOLES EN MAJUSCULES]} : même convention de
    nommage que gseapy avec plusieurs bibliothèques ; symboles en majuscules comme
    le classement (_ranking_metric), donc Actb (souris) = ACTB (humain).
    Une source introuvable est signalée et ignorée."""
    import gseapy
    version = msigdb_version or DEFAULT_MSIGDB_VERSION.get(organism, "2025.1.Hs")
    sets = {}
    for spec in specs:
        spec = str(spec).strip()
        try:
            if os.path.isfile(spec):
                src = os.path.splitext(os.path.basename(spec))[0]
                lib = read_gmt(spec)
            elif spec.lower().startswith("msigdb:"):
                cat = spec.split(":", 1)[1]
                src = cat
                lib = gseapy.Msigdb.get_gmt(category=cat, dbver=version)
                if not lib:
                    raise ValueError(f"MSigDB {cat} v{version} not found")
            else:
                src = spec
                lib = gseapy.get_library(name=spec, organism=organism.capitalize())
        except Exception as e:                                        # noqa: BLE001
            if verbose:
                print(f"  [WARN] GSEA gene sets '{spec}' unavailable "
                      f"({type(e).__name__}: {str(e)[:80]}) — skipped. Offline: give "
                      f"a local .gmt file in gsea_gene_sets.")
            continue
        n0 = len(sets)
        for term, genes in lib.items():
            g = list(dict.fromkeys(str(x).strip().upper() for x in genes if str(x).strip()))
            if g:
                sets[f"{src}__{term}"] = g
        if verbose:
            print(f"     {src}: {len(sets) - n0} gene sets")
    return sets


def _ranking_metric(df_comparaison, contrast, params):
    """Construit le vecteur de ranking signé pour un contraste.

    Métrique = -log10(p) * sign(LFC), sur TOUTES les protéines (pas seulement
    les DEP). Mapping accession -> symbole de gène via la colonne 'Genes'.
    p : p brute par défaut (gsea_rank_p: pvalue). La p.adj BH forme des paliers
    (beaucoup de protéines à la même valeur) qui écrasent l'ordre du milieu du
    classement, là où la GSEA détecte les effets coordonnés ; gsea_rank_p: volcano
    reprend l'ancien comportement (p.adj si volcano_use_padj).

    Retourne une pd.Series {gene_symbol: score} triée décroissante, ou None.
    """
    if str(params.get("gsea_rank_p", "pvalue")).lower() == "volcano":
        suffixe_p = "_p.adj" if params.get("volcano_use_padj") else "_p.val"
    else:
        suffixe_p = "_p.val"
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

    cfg       = gsea_settings(params, go_params)
    fdr_max   = cfg["fdr"]
    n_perm    = cfg["permutations"]
    p_lbl = ("p.adj" if cfg["rank_p"] == "volcano" and params.get("volcano_use_padj")
             else "p.value")

    print(f"\n[GSEA] Pre-ranked GSEA — metric: -log10({p_lbl}) x sign(LFC) — "
          f"{cfg['organism']} gene sets"
          + (f" (MSigDB {cfg['msigdb_version']})"
             if any(str(s).lower().startswith("msigdb:") for s in cfg["gene_sets"]) else ""))
    gene_sets = load_gene_sets(cfg["gene_sets"], cfg["organism"], cfg["msigdb_version"])
    if not gene_sets:
        print("  [WARN] GSEA skipped: no gene set could be loaded.")
        return None
    universe = {g for v in gene_sets.values() for g in v}
    results = {}
    coverage_done = False
    for contrast in contrast_names:
        try:
            rnk = _ranking_metric(df_comparaison, contrast, params)
            if rnk is None:
                print(f"  [SKIP] {contrast}: not enough ranked genes.")
                continue
            if not coverage_done:
                n_in = int(rnk.index.isin(list(universe)).sum())
                print(f"  [GSEA] {n_in}/{len(rnk)} ranked genes found in the gene sets "
                      f"({100 * n_in / len(rnk):.0f} %)"
                      + ("  [WARN] low (expected only for a small custom collection): "
                         "check the Genes column / gsea_organism"
                         if n_in < 0.4 * len(rnk) else ""))
                coverage_done = True

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
