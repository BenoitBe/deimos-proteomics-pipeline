#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_dashboard.py — Générateur du dashboard protéomique HTML (module 4 du skill
ProteomicDataProcess), porté de la logique du skill vers Python.

PRINCIPE FONDAMENTAL
--------------------
Toutes les valeurs statistiques (p.val, p.adj, LFC, Pi, z-score GO, q-value,
clusters, HubScore...) proviennent EXCLUSIVEMENT des fichiers Excel produits par
les scripts R. Le script ne RECALCULE jamais ces statistiques : il lit, extrait,
embarque.

Seules exceptions explicitement prévues par le skill : les projections de
visualisation PCA et Scatter, car les feuilles `PCA` / `Scatter_Plots` du rapport
Excel ne contiennent que des images (pas de coordonnées). Ces projections sont
dérivées géométriquement des colonnes d'intensité log2 LFQ (aucune statistique
n'est recalculée).

ENTRÉES
-------
  rapportstatistique.xlsx        (obligatoire)  -> sheets Differential_Expression, Zscore_Heatmap,
                                                   ANOVA_Results, UpSet_Intersections
  rapportenrichissementGO.xlsx   (optionnel)    -> 1 sheet GO_* par contraste + images
  WGCNA_Emotional.xlsx           (optionnel)    -> sheet WGCNA_Results

SORTIE
------
  proteogen_dashboard.html  (template + bloc de données injecté)

USAGE
-----
  python build_dashboard.py                         # auto-détection dans /mnt/user-data/uploads
  python build_dashboard.py --uploads ./data        # auto-détection dans un dossier
  python build_dashboard.py --stats a.xlsx --go b.xlsx --wgcna c.xlsx --out dash.html
"""

import os
import re
import io
import sys
import json
import base64
import argparse
import openpyxl #type: ignore
import numpy as np
import pandas as pd

# ──────────────────────────────────────────────────────────────────────────────
# Constantes de seuils (identiques au skill / dashboard-python.md)
# ──────────────────────────────────────────────────────────────────────────────
P_THRESH = 0.05
LFC_THRESH = 0.26
NS_VOLCANO_MAX = 1200   # NS échantillonnés dans le volcano (nuage dense)
NS_SCATTER_MAX = 2000   # NS échantillonnés dans le scatter (nuage dense comme l'Excel)
NS_CV          = 400    # Points sous-échantillonnés pour le CV plot
SEED = 42

META = {'Protein.Group', 'Protein.Names', 'Genes', 'First.Protein.Description',
        'N.Sequences', 'N.Proteotypic.Sequences', 'name',
        'imputed', 'num_NAs', 'Peptide_Count_min', 'significant'}
NON_SAMPLE = {'name', 'Protein.Group', 'Protein.Names', 'Genes',
              'First.Protein.Description', 'N.Sequences', 'N.Proteotypic.Sequences',
              'ID', 'imputed', 'num_NAs', 'significant', 'Peptide_Count_min'}

IMG_LABELS = ['lollipop', 'dotplot']  # ordre d'insertion dans les sheets GO


# ──────────────────────────────────────────────────────────────────────────────
# Détection des fichiers (par contenu, conformément au skill)
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# VALIDATION DES FICHIERS D'ENTRÉE
# ──────────────────────────────────────────────────────────────────────────────

def validate_inputs(stats_path, go_path, wgcna_path, umap_path):
    """Vérifie la présence des sheets et colonnes obligatoires.
    Retourne une liste de warnings (non bloquants) et d'erreurs (bloquantes).
    """
    warnings = []
    errors   = []

    # ── Stats (obligatoire) ───────────────────────────────────────────────────
    try:
        xl = pd.ExcelFile(stats_path)
        sheets = set(xl.sheet_names)

        required_sheets = ['Differential_Expression', 'Zscore_Heatmap', 'ANOVA_Results']
        for s in required_sheets:
            if s not in sheets:
                errors.append(f"stats : sheet '{s}' manquante")

        if 'Differential_Expression' in sheets:
            df = pd.read_excel(stats_path, sheet_name='Differential_Expression', nrows=2)
            has_pval = any(str(c).endswith('_p.val') for c in df.columns)
            if not has_pval:
                errors.append("stats/Differential_Expression : aucune colonne *_p.val détectée "
                               "(contrastes manquants)")
            if 'name' not in df.columns:
                warnings.append("stats/Differential_Expression : colonne 'name' absente — "
                                 "les accessions ne seront pas disponibles")

        if 'UpSet_Intersections' not in sheets:
            warnings.append("stats : sheet 'UpSet_Intersections' absente — "
                             "UpSet plot désactivé")
        if 'raw_data' not in sheets:
            warnings.append("stats : sheet 'raw_data' absente — "
                             "figure Missing Values désactivée")

    except Exception as e:
        errors.append(f"stats : impossible de lire le fichier ({e})")

    # ── GO (optionnel) ────────────────────────────────────────────────────────
    if go_path:
        try:
            xl = pd.ExcelFile(go_path)
            go_sheets = [s for s in xl.sheet_names
                         if s.startswith('GO_') or s.lower().startswith('go')]
            if not go_sheets:
                warnings.append("GO : aucune sheet GO_* détectée")
            else:
                df = pd.read_excel(go_path, sheet_name=go_sheets[0], nrows=2)
                df = df[[c for c in df.columns if not str(c).startswith('Unnamed')]]
                for col in ['term_name', 'p_value', 'z_score']:
                    if col not in df.columns:
                        warnings.append(f"GO/{go_sheets[0]} : colonne '{col}' absente")
        except Exception as e:
            warnings.append(f"GO : impossible de lire le fichier ({e})")

    # ── WGCNA (optionnel) ─────────────────────────────────────────────────────
    if wgcna_path:
        try:
            xl = pd.ExcelFile(wgcna_path)
            wgcna_sheets = [s for s in xl.sheet_names if 'wgcna' in s.lower()]
            if not wgcna_sheets:
                warnings.append("WGCNA : aucune sheet contenant 'wgcna' détectée")
            else:
                df = pd.read_excel(wgcna_path, sheet_name=wgcna_sheets[0], nrows=2)
                for col in ['Gene_Name', 'Module_Color', 'HubScore']:
                    if col not in df.columns:
                        warnings.append(f"WGCNA/{wgcna_sheets[0]} : "
                                         f"colonne '{col}' absente")
        except Exception as e:
            warnings.append(f"WGCNA : impossible de lire le fichier ({e})")

    # ── UMAP (optionnel) ──────────────────────────────────────────────────────
    if umap_path:
        try:
            xl = pd.ExcelFile(umap_path)
            exact = [s for s in xl.sheet_names if s.strip().lower() == 'umap']
            umap_sheets = exact or [s for s in xl.sheet_names
                                    if 'umap' in s.lower() and 'pca' not in s.lower()]
            if not umap_sheets:
                warnings.append("UMAP : aucune sheet contenant 'umap' détectée")
            else:
                df = pd.read_excel(umap_path, sheet_name=umap_sheets[0], nrows=2)
                has_umap1 = any(str(c).lower() in ('umap1','umap_1','dim1') for c in df.columns)
                has_umap2 = any(str(c).lower() in ('umap2','umap_2','dim2') for c in df.columns)
                if not has_umap1 or not has_umap2:
                    warnings.append("UMAP : colonnes UMAP1/UMAP2 non trouvées")
        except Exception as e:
            warnings.append(f"UMAP : impossible de lire le fichier ({e})")

    return warnings, errors

def detect_files(uploads):
    """Reconnaît stats / GO / WGCNA / UMAP par leur contenu, pas par leur nom.

    Un même fichier peut remplir plusieurs rôles (ex : stats + WGCNA dans le même
    Excel). L'ordre de priorité pour l'attribution unique reste :
      1. Stats   (Differential_Expression ou Zscore présent)
      2. WGCNA   (sheet contenant 'wgcna')
      3. GO      (sheet commençant par GO_)
      4. UMAP    (sheet contenant 'umap')
    Si un fichier contient à la fois stats et WGCNA, il est assigné aux DEUX.
    """
    found = {'stats': None, 'go': None, 'wgcna': None, 'umap': None}
    if not os.path.isdir(uploads):
        return found
    for fn in sorted(os.listdir(uploads)):
        if not fn.lower().endswith('.xlsx'):
            continue
        path = os.path.join(uploads, fn)
        try:
            xl = pd.ExcelFile(path)
        except Exception:
            continue
        sheets = set(xl.sheet_names)

        is_stats = ('Differential_Expression' in sheets or
                    any('zscore' in s.lower() for s in sheets))
        is_wgcna = any('wgcna' in s.lower() for s in sheets)
        is_go    = any(s.startswith('GO_') or s.lower().startswith('go')
                       for s in sheets)
        is_umap  = any('umap' in s.lower() for s in sheets)

        # Assigner chaque rôle indépendamment (premier fichier qui convient gagne)
        if is_stats and found['stats'] is None:
            found['stats'] = path
        if is_wgcna and found['wgcna'] is None:
            found['wgcna'] = path
        if is_go and found['go'] is None:
            found['go'] = path
        if is_umap and found['umap'] is None:
            found['umap'] = path

    return found



# ──────────────────────────────────────────────────────────────────────────────
# MODE OFFLINE — embarquer Chart.js et xlsx.js dans le HTML
# ──────────────────────────────────────────────────────────────────────────────

def embed_libs(html, script_dir):
    """Remplace les balises <script src="cdnjs..."> par des scripts inline embarqués.
    Cherche chart.umd.js et xlsx.full.min.js dans script_dir (répertoire du script).
    Si introuvables, conserve les balises CDN (fallback).
    """
    import re as _re

    lib_files = {
        'chart.js':     ['chart.umd.js', 'chart.umd.min.js', 'chart.min.js'],
        'xlsx':         ['xlsx.full.min.js', 'xlsx.min.js'],
    }
    cdn_patterns = {
        'chart.js': _re.compile(r'<script[^>]+cdnjs[^>]+[Cc]hart[^>]+></script>'),
        'xlsx':     _re.compile(r'<script[^>]+cdnjs[^>]+xlsx[^>]+></script>'),
    }

    for lib_key, filenames in lib_files.items():
        content = None
        for fn in filenames:
            path = os.path.join(script_dir, fn)
            if os.path.isfile(path):
                with open(path, encoding='utf-8', errors='replace') as f:
                    content = f.read()
                size_kb = len(content.encode('utf-8')) // 1024
                print(f"  embed_libs : {fn} embarqué ({size_kb} KB).")
                break

        if content is None:
            print(f"  [warn] embed_libs : {lib_key} non trouvé dans {script_dir} — CDN conservé.")
            continue

        # Remplacer la balise CDN par un script inline
        pattern = cdn_patterns[lib_key]
        # Sécurité : seul '</script>' (avec >) peut fermer la balise prématurément
        # On remplace uniquement la séquence exacte '</script>' (insensible casse, avec >)
        import re as _re2
        safe_content = _re2.sub(
            r'<(/\s*script\s*>)',
            lambda m: '<\\' + m.group(1),
            content,
            flags=_re2.IGNORECASE
        )
        inline = '<script>' + safe_content + '</script>'
        if pattern.search(html):
            html = pattern.sub(lambda _: inline, html, count=1)
        else:
            print(f"  [warn] embed_libs : balise CDN {lib_key} introuvable dans le HTML.")

    return html

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────
def norm(s):
    """Normalise un nom de condition pour le matching (Polyculture.1 == Polyculture_1
    == Polyculture-1). Retire points, underscores ET tirets, casse ignorée."""
    return str(s).lower().replace('.', '').replace('_', '').replace('-', '')


def find_cond_cols(cond_name, cond_map):
    """Colonnes d'une condition à partir d'un nom issu d'un contraste (noms normalisés).
    cond_map est {condition: [colonnes]}."""
    n = norm(cond_name)
    for cond, cols in cond_map.items():
        if norm(cond) == n:
            return list(cols)
    return []


def safe_str(v, limit=None):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ''
    s = str(v)
    return s[:limit] if limit else s


# ──────────────────────────────────────────────────────────────────────────────
# Lecture du fichier statistique
# ──────────────────────────────────────────────────────────────────────────────
def read_stats(path):
    df_comp = pd.read_excel(path, sheet_name='Differential_Expression')

    # Détection automatique des contrastes
    contrasts = [c[:-len('_p.val')] for c in df_comp.columns if str(c).endswith('_p.val')]
    if not contrasts:
        raise ValueError("Aucun contraste détecté (colonnes *_p.val absentes de 'Differential_Expression').")

    df_z = pd.read_excel(path, sheet_name='Zscore_Heatmap')
    # sample_cols = colonnes numériques d'intensité (hors métadonnées)
    sample_cols = [c for c in df_z.columns
                   if c not in META and pd.api.types.is_numeric_dtype(df_z[c])]
    # Fallback : si df_z vide ou pas de colonnes numériques,
    # déduire sample_cols depuis la sheet Comparaison
    if not sample_cols:
        stat_suffixes = ('_p.val','_p.adj','_diff','_Pi_Score','_Robustness_Score',
                         'Pi_Score_','Robustness_Score_')
        sample_cols = [c for c in df_comp.columns
                       if c not in NON_SAMPLE
                       and not any(c.endswith(s) or c.startswith(s)
                                   for s in stat_suffixes)
                       and pd.api.types.is_numeric_dtype(df_comp[c])]

    df_anova = pd.read_excel(path, sheet_name='ANOVA_Results')

    # Intersections : header en ligne 2
    try:
        df_inter = pd.read_excel(path, sheet_name='UpSet_Intersections', header=1)
        df_inter = df_inter[[c for c in df_inter.columns
                             if not str(c).startswith('Unnamed')]]
        if 'Protein.Group' in df_inter.columns:
            df_inter = df_inter.dropna(subset=['Protein.Group'])
    except Exception:
        df_inter = pd.DataFrame()

    return df_comp, df_z, df_anova, df_inter, contrasts, sample_cols


# ──────────────────────────────────────────────────────────────────────────────
# Conditions / cond_map
# ──────────────────────────────────────────────────────────────────────────────
def build_conditions(sample_cols):
    cond_map_inv = {col: col.rsplit('_', 1)[0] for col in sample_cols}
    conditions = list(dict.fromkeys(cond_map_inv.values()))
    cond_map = {cond: [c for c in sample_cols if cond_map_inv[c] == cond]
                for cond in conditions}
    return conditions, cond_map


# ──────────────────────────────────────────────────────────────────────────────
# VOLCANO
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# PROTEOME_DATA — structure unifiée remplaçant VOLCANO + SCATTER_DATA + MA_DATA
# ──────────────────────────────────────────────────────────────────────────────

def build_proteome_data(df_comp, contrasts, cond_map):
    """Structure unifiée par contraste : un seul point par protéine contenant
    tous les champs nécessaires aux volcanos, scatter et MA plots.
    Champs : name, g (gène), desc, d (LFC), p, padj, r, pi, s, x (moy condA), y (moy condB)
    Champs dérivés calculés en JS : lp=-log10(p), a=(x+y)/2
    """
    result = {}
    for c in contrasts:
        parts = c.split('_vs_')
        # Utiliser find_cond_cols pour gérer les différences . vs _ dans les noms
        cols_a = [col for col in find_cond_cols(parts[0], cond_map) if col in df_comp.columns]
        cols_b = [col for col in find_cond_cols(parts[1] if len(parts) > 1 else '', cond_map) if col in df_comp.columns]

        pcol  = c + '_p.val'
        dcol  = c + '_diff'
        pacol = c + '_p.adj'
        rcol  = 'Robustness_Score_' + c
        picol = 'Pi_Score_' + c

        if pcol not in df_comp.columns or dcol not in df_comp.columns:
            continue

        mean_a = df_comp[cols_a].fillna(0).mean(axis=1) if cols_a else pd.Series([0.0]*len(df_comp))
        mean_b = df_comp[cols_b].fillna(0).mean(axis=1) if cols_b else pd.Series([0.0]*len(df_comp))

        pts = []
        for _, row in df_comp.iterrows():
            pv = row.get(pcol)
            dv = row.get(dcol)
            if pd.isna(pv) or pd.isna(dv):
                continue
            pv = float(pv); dv = float(dv)
            padj = float(row[pacol]) if pacol in df_comp.columns and pd.notna(row.get(pacol)) else None
            r    = int(row[rcol])    if rcol  in df_comp.columns and pd.notna(row.get(rcol))  else 0
            pi   = float(row[picol]) if picol in df_comp.columns and pd.notna(row.get(picol)) else None
            s    = 'N'
            if pv < P_THRESH and dv >  LFC_THRESH: s = 'U'
            if pv < P_THRESH and dv < -LFC_THRESH: s = 'D'
            xa = round(float(mean_a[row.name]), 3)
            xb = round(float(mean_b[row.name]), 3)
            pt = {
                'name': safe_str(row.get('name')),
                'g':    safe_str(row.get('Genes')),
                'desc': safe_str(row.get('First.Protein.Description'), 70),
                'd':    round(dv, 4),
                'p':    round(pv, 6),
                'padj': round(padj, 6) if padj is not None else None,
                'r':    r,
                'pi':   round(pi, 4) if pi is not None else None,
                'lp':   round(float(-np.log10(max(pv, 1e-10))), 4),
                's':    s,
                'x':    xa,
                'y':    xb,
            }
            pts.append(pt)

        # Sous-échantillonner les NS pour alléger
        import random as _rnd
        _rnd.seed(SEED)
        sig = [p for p in pts if p['s'] != 'N']
        ns  = [p for p in pts if p['s'] == 'N']
        ns_sampled = _rnd.sample(ns, min(NS_VOLCANO_MAX, len(ns)))
        result[c] = sig + ns_sampled

    return result

def build_volcano(df_comp, contrasts):
    volcano = {}
    for c in contrasts:
        cols = ['name', 'Genes', 'First.Protein.Description',
                c + '_diff', c + '_p.val', c + '_p.adj',
                'Robustness_Score_' + c, 'Pi_Score_' + c]
        cols = [col for col in cols if col in df_comp.columns]
        sub = df_comp[cols].copy()
        ren = {c + '_diff': 'd', c + '_p.val': 'p', c + '_p.adj': 'padj',
               'Robustness_Score_' + c: 'r', 'Pi_Score_' + c: 'pi',
               'First.Protein.Description': 'desc'}
        sub = sub.rename(columns=ren)
        sub = sub.dropna(subset=['p', 'd'])
        for need, default in (('padj', np.nan), ('r', 0), ('pi', np.nan), ('desc', '')):
            if need not in sub.columns:
                sub[need] = default
        sub['lp'] = np.round(-np.log10(sub['p'].clip(lower=1e-10)), 4)
        sub['d'] = sub['d'].round(4)
        sub['p'] = sub['p'].round(6)
        sub['padj'] = sub['padj'].round(6)
        sub['pi'] = sub['pi'].round(4)
        sub['r'] = sub['r'].fillna(0).round().astype(int)
        sub['desc'] = sub['desc'].fillna('').astype(str).str.slice(0, 70)
        sub['s'] = 'N'
        sub.loc[(sub['p'] < P_THRESH) & (sub['d'] > LFC_THRESH), 's'] = 'U'
        sub.loc[(sub['p'] < P_THRESH) & (sub['d'] < -LFC_THRESH), 's'] = 'D'
        sig = sub[sub['s'] != 'N']
        ns_all = sub[sub['s'] == 'N']
        ns = ns_all.sample(min(NS_VOLCANO_MAX, len(ns_all)), random_state=SEED) \
            if len(ns_all) else ns_all
        out = pd.concat([sig, ns])[
            ['name', 'Genes', 'desc', 'd', 'p', 'padj', 'r', 'pi', 'lp', 's']
        ].where(pd.notna, None)
        volcano[c] = json.loads(out.to_json(orient='records'))
    return volcano


# ──────────────────────────────────────────────────────────────────────────────
# STATS
# ──────────────────────────────────────────────────────────────────────────────
def build_stats(df_comp, df_anova, contrasts):
    sig_anova = int((df_anova.get('significant') == True).sum()) if 'significant' in df_anova else 0
    clusters = {}
    if 'Cluster_ID' in df_anova.columns:
        clusters = {str(k): int(v)
                    for k, v in df_anova['Cluster_ID'].value_counts().dropna().items()}
    summary = {}
    dep_union = set()  # protéines DEP dans au moins un contraste
    for c in contrasts:
        pcol, dcol = c + '_p.val', c + '_diff'
        if pcol in df_comp.columns and dcol in df_comp.columns:
            up_mask = (df_comp[pcol] < P_THRESH) & (df_comp[dcol] > LFC_THRESH)
            dn_mask = (df_comp[pcol] < P_THRESH) & (df_comp[dcol] < -LFC_THRESH)
            up = int(up_mask.sum())
            dn = int(dn_mask.sum())
            summary[c] = {'up': up, 'down': dn}
            # Protéines significatives (up ou down) pour ce contraste
            sig_idx = df_comp.index[up_mask | dn_mask]
            dep_union.update(sig_idx)
    total_dep = len(dep_union)
    return {
        'total_proteins': int(len(df_comp)),
        'sig_anova': sig_anova,
        'total_dep': total_dep,
        'clusters': clusters,
        'contrasts_summary': summary,
    }


# ──────────────────────────────────────────────────────────────────────────────
# HEATMAP (protéines ANOVA significatives, z-scores groupés)
# ──────────────────────────────────────────────────────────────────────────────
def build_heatmap(df_z, df_anova, sample_cols):
    # Cas 2 conditions : df_z peut être vide ou sans colonnes d'échantillons
    if df_z.empty or not any(c in df_z.columns for c in sample_cols):
        return []
    keep = [c for c in ['Protein.Group', 'Cluster_ID', 'qvalue', 'significant']
            if c in df_anova.columns]
    df_z_m = df_z.merge(df_anova[keep], on='Protein.Group', how='left') \
        if 'Protein.Group' in df_z.columns and 'Protein.Group' in keep else df_z.copy()
    if 'significant' in df_z_m.columns:
        df_z_m = df_z_m[df_z_m['significant'] == True]
    sort_cols = [c for c in ['Cluster_ID', 'qvalue'] if c in df_z_m.columns]
    if sort_cols:
        df_z_m = df_z_m.sort_values(sort_cols)

    rows = []
    for _, row in df_z_m.iterrows():
        # Ignorer les lignes dont les colonnes sample ne sont pas numériques
        try:
            vals = [round(float(row[c]), 3) if c in row.index and pd.notna(row[c]) else 0.0
                    for c in sample_cols]
        except (ValueError, TypeError):
            continue
        gene = row['Genes'] if 'Genes' in row and pd.notna(row['Genes']) else row.get('Protein.Group')
        rows.append({
            'g': safe_str(gene),
            'cl': safe_str(row.get('Cluster_ID')) or 'NA',
            'q': round(float(row['qvalue']), 5) if pd.notna(row.get('qvalue')) else 1.0,
            'pg': safe_str(row.get('Protein.Group')),
            'v': vals,
        })
    return rows


# ──────────────────────────────────────────────────────────────────────────────
# TOP 20 Pi
# ──────────────────────────────────────────────────────────────────────────────
def build_top20_pi(df_comp, contrasts):
    all_pi = {}
    for c in contrasts:
        pi_col = 'Pi_Score_' + c
        if pi_col not in df_comp.columns:
            continue
        cols = [col for col in ['name', 'Genes', 'First.Protein.Description', pi_col]
                if col in df_comp.columns]
        for _, row in df_comp[cols].dropna(subset=[pi_col]).iterrows():
            k = row.get('name')
            pi = float(row[pi_col])
            if k not in all_pi or all_pi[k]['pi'] < pi:
                all_pi[k] = {
                    'name': safe_str(k),
                    'gene': safe_str(row.get('Genes')),
                    'desc': safe_str(row.get('First.Protein.Description'), 70),
                    'pi': pi,            # non arrondi (conforme à l'exemple)
                }
    return sorted(all_pi.values(), key=lambda x: x['pi'], reverse=True)[:20]


# ──────────────────────────────────────────────────────────────────────────────
# INTERSECTIONS + UPSET
# ──────────────────────────────────────────────────────────────────────────────
def short_label(part):
    """Label court d'une condition, SANS contexte (rétrocompat / fallback).
    Conserve le mapping historique Monoculture/Polyculture, mais ne tronque plus
    aveuglément au début (source du bug X180/X180 sur préfixes communs) : on
    privilégie la fin distinctive si le nom est long.
    """
    s = part.replace('Polyculture_', 'P').replace('Polyculture.', 'P').replace('Monoculture', 'M')
    s = ''.join(s.split('_'))
    if len(s) <= 6:
        return s
    return s[-6:]   # garde la fin (souvent distinctive) plutôt que le début


def _strip_common_affix(conditions):
    """À partir de la liste complète des conditions, retire le préfixe ET le
    suffixe communs à TOUTES, pour ne garder que la partie discriminante.

    Ex: ['18-009IG01Ctrl','18-009IG01di6h','18-009IG01di24h']
        -> {'18-009IG01Ctrl':'Ctrl', '18-009IG01di6h':'di6h', '18-009IG01di24h':'di24h'}

    Garanties :
    - jamais de label vide (repli sur le nom complet tronqué)
    - unicité des labels (suffixe numérique #i ajouté en cas de collision)
    """
    conds = list(dict.fromkeys(conditions))  # uniques, ordre conservé
    if len(conds) <= 1:
        return {c: (c[:10] if len(c) > 10 else c) for c in conds}

    # Préfixe commun (sur caractères)
    pre = os.path.commonprefix(conds)
    # Suffixe commun (via reverse)
    suf = os.path.commonprefix([c[::-1] for c in conds])[::-1]

    # On ne coupe le préfixe/suffixe que s'il reste quelque chose de non vide
    # pour TOUTES les conditions (sinon on annulerait la coupe).
    def _trim(c):
        core = c
        if pre and core.startswith(pre):
            core = core[len(pre):]
        if suf and core.endswith(suf) and len(core) > len(suf):
            core = core[:len(core) - len(suf)]
        core = core.strip(' _-.')
        return core if core else c   # jamais vide

    trimmed = {c: _trim(c) for c in conds}

    # Si la coupe a rendu des labels vides ou tous identiques, on annule la coupe
    vals = list(trimmed.values())
    if any(not v for v in vals) or len(set(vals)) < len(vals):
        # Tentative : ne couper que le préfixe (garder le suffixe distinctif)
        def _trim_pre_only(c):
            core = c[len(pre):] if (pre and c.startswith(pre)) else c
            core = core.strip(' _-.')
            return core if core else c
        trimmed2 = {c: _trim_pre_only(c) for c in conds}
        if len(set(trimmed2.values())) == len(trimmed2):
            trimmed = trimmed2

    # Raccourcir les labels trop longs (garde le début de la partie distinctive)
    out = {}
    seen = {}
    for c, lbl in trimmed.items():
        s = lbl if len(lbl) <= 8 else lbl[:8]
        # Garantir l'unicité
        if s in seen.values():
            i = 2
            while f"{s[:6]}#{i}" in seen.values():
                i += 1
            s = f"{s[:6]}#{i}"
        seen[c] = s
        out[c] = s
    return out


def build_intersections(df_inter, contrasts):
    inter_rows = []
    if df_inter.empty:
        return inter_rows, {'sets': [], 'full_labels': {}, 'set_sizes': {}, 'intersections': []}

    # Table INTERSECTIONS (renommage contrastes -> c1..cN)
    cmap = {c: 'c%d' % (i + 1) for i, c in enumerate(contrasts)}
    for _, row in df_inter.iterrows():
        rec = {
            'pg': safe_str(row.get('Protein.Group')),
            'g': safe_str(row.get('Genes')),
            'desc': safe_str(row.get('First.Protein.Description'), 70),
            'nb': int(row['Nb_Conditions']) if pd.notna(row.get('Nb_Conditions')) else 0,
        }
        for c in contrasts:
            v = row.get(c)
            rec[cmap[c]] = 'X' if (isinstance(v, str) and v.strip() == 'X') else ''
        inter_rows.append(rec)

    # UPSET : labels courts par contraste
    # On collecte d'abord TOUTES les conditions (les 2 côtés de chaque contraste)
    # pour retirer le préfixe/suffixe commun à l'ensemble (évite le bug X180/X180
    # quand les conditions partagent un long préfixe de code projet/date/lignée).
    all_conds = []
    for c in contrasts:
        parts = c.split('_vs_')
        all_conds.extend(parts)
    cond_label = _strip_common_affix(all_conds)

    contrast_short = {}
    for c in contrasts:
        parts = c.split('_vs_')
        a = cond_label.get(parts[0], short_label(parts[0])) if len(parts) > 0 else c[:3]
        b = cond_label.get(parts[1], short_label(parts[1])) if len(parts) > 1 else ''
        contrast_short[c] = (a + '/' + b) if b else a

    sets = list(contrast_short.values())
    full_labels = {v: k for k, v in contrast_short.items()}
    set_sizes = {contrast_short[c]: int((df_inter[c] == 'X').sum())
                 for c in contrasts if c in df_inter.columns}

    # Signatures binaires exactes
    present = [c for c in contrasts if c in df_inter.columns]
    bin_df = pd.DataFrame({c: (df_inter[c] == 'X').astype(int) for c in present})
    bin_df['sig'] = bin_df[present].apply(tuple, axis=1)
    counts = bin_df.groupby('sig').size().reset_index(name='count') \
        .sort_values('count', ascending=False)

    intersections = []
    for _, row in counts.iterrows():
        active = [contrast_short[present[i]] for i, v in enumerate(row['sig']) if v == 1]
        if active:
            intersections.append({'sets': active, 'count': int(row['count'])})

    upset = {
        'sets': sets, 'full_labels': full_labels,
        'set_sizes': set_sizes, 'intersections': intersections[:30],
    }
    return inter_rows, upset


# ──────────────────────────────────────────────────────────────────────────────
# PCA (projection de visualisation — colonnes d'intensité log2 LFQ, toutes protéines)
# ──────────────────────────────────────────────────────────────────────────────
def build_pca(df_comp, sample_cols, cond_map_inv):
    # Colonnes d'intensité présentes dans Comparaison (priorité aux noms SAMPLE_COLS)
    pca_cols = [c for c in sample_cols if c in df_comp.columns]
    if len(pca_cols) < 2:
        # repli : détection des colonnes d'intensité dans Comparaison
        stat_like = set()
        for col in df_comp.columns:
            for suf in ('_diff', '_p.val', '_p.adj'):
                if str(col).endswith(suf):
                    stat_like.add(col)
            if str(col).startswith(('Robustness_Score_', 'Pi_Score_')):
                stat_like.add(col)
        pca_cols = [c for c in df_comp.columns if c not in NON_SAMPLE and c not in stat_like]
    if len(pca_cols) < 2:
        return {'points': [], 'var_exp': []}

    mat = df_comp[pca_cols].fillna(0).values.T          # samples × proteins
    mat_c = mat - mat.mean(axis=0)                       # centrage par protéine
    U, S, Vt = np.linalg.svd(mat_c, full_matrices=False)
    scores = U * S
    var_exp = (S ** 2 / np.sum(S ** 2) * 100).round(2)

    points = []
    for i, samp in enumerate(pca_cols):
        cond = cond_map_inv.get(samp, samp.rsplit('_', 1)[0])
        points.append({
            'sample': samp, 'cond': cond,
            'pc1': round(float(scores[i, 0]), 3),
            'pc2': round(float(scores[i, 1]), 3) if scores.shape[1] > 1 else 0.0,
            'pc3': round(float(scores[i, 2]), 3) if scores.shape[1] > 2 else 0.0,
        })
    return {'points': points, 'var_exp': var_exp[:3].tolist()}


# ──────────────────────────────────────────────────────────────────────────────
# SCATTER (intensités log2 LFQ moyennes par condition + r de Pearson)
# ──────────────────────────────────────────────────────────────────────────────
def build_scatter(df_comp, contrasts, cond_map):
    scatter = {}
    for c in contrasts:
        parts = c.split('_vs_')
        if len(parts) != 2:
            continue
        cols_a = [col for col in find_cond_cols(parts[0], cond_map) if col in df_comp.columns]
        cols_b = [col for col in find_cond_cols(parts[1], cond_map) if col in df_comp.columns]
        if not cols_a or not cols_b:
            continue
        mean_a = df_comp[cols_a].fillna(0).mean(axis=1)
        mean_b = df_comp[cols_b].fillna(0).mean(axis=1)
        pval = df_comp[c + '_p.val'].fillna(1) if c + '_p.val' in df_comp.columns else pd.Series([1.0] * len(df_comp))
        diff = df_comp[c + '_diff'].fillna(0) if c + '_diff' in df_comp.columns else pd.Series([0.0] * len(df_comp))
        padj = df_comp[c + '_p.adj'] if c + '_p.adj' in df_comp.columns else pd.Series([np.nan] * len(df_comp))

        status = pd.Series(['N'] * len(df_comp), index=df_comp.index)
        status[(pval < P_THRESH) & (diff > LFC_THRESH)] = 'U'
        status[(pval < P_THRESH) & (diff < -LFC_THRESH)] = 'D'

        df_sc = pd.DataFrame({
            'x': mean_a.round(3), 'y': mean_b.round(3), 's': status,
            'g': df_comp['Genes'].fillna('') if 'Genes' in df_comp else '',
            'n': df_comp['name'] if 'name' in df_comp else '',
            'd': (df_comp['First.Protein.Description'].fillna('').astype(str).str.slice(0, 70)
                  if 'First.Protein.Description' in df_comp else ''),
            'lfc': diff.round(3), 'pv': pval.round(6), 'pa': padj.round(6),
        })
        # r de Pearson sur les intensités positives
        mask = (mean_a > 0) & (mean_b > 0)
        r = float(np.corrcoef(mean_a[mask], mean_b[mask])[0, 1]) if mask.sum() > 1 else 0.0

        sig = df_sc[df_sc['s'] != 'N']
        ns_all = df_sc[df_sc['s'] == 'N']
        ns = ns_all.sample(min(NS_SCATTER_MAX, len(ns_all)), random_state=SEED) if len(ns_all) else ns_all
        pts = pd.concat([sig, ns]).where(pd.notna, None)
        scatter[c] = {
            'pts': json.loads(pts.to_json(orient='records')),
            'r': round(r, 4),
            'label_a': parts[0].replace('_', ' '),
            'label_b': parts[1].replace('_', ' '),
        }
    return scatter


# ──────────────────────────────────────────────────────────────────────────────
# GO — données, images, chord
# ──────────────────────────────────────────────────────────────────────────────
def map_go_sheet(sheet, contrasts):
    """Apparie une sheet GO_* au contraste correspondant.

    Les onglets sont nommés 'GO_<contrast>' avec '_vs_' -> 'v' et tronqués à 31
    caractères (limite Excel). L'ancienne version comparait des caractères triés
    (similarité floue), ce qui échoue dès que les contrastes partagent un long
    préfixe commun : plusieurs onglets se mappaient au même contraste et un
    contraste restait orphelin (colonne vide dans la heatmap).

    Nouvelle stratégie : on encode chaque contraste comme le pipeline encode le
    nom d'onglet (remplacement '_vs_' -> 'v', suppression de '_', troncature à
    28 caractères = 31 - len('GO_')), puis on cherche une correspondance EXACTE
    sur ce préfixe. Repli sur la plus longue sous-chaîne commune seulement si
    aucune correspondance exacte n'est trouvée.
    """
    s = sheet[3:] if sheet.startswith('GO_') else sheet  # retire 'GO_'

    # Les noms d'onglets sont désormais construits à partir des conditions
    # STRIPPÉES de leur préfixe/suffixe commun (cf. _make_go_sheet_names dans
    # go_enrichment.py). On reproduit ici le même stripping pour matcher.
    all_conds = []
    for c in contrasts:
        all_conds.extend(c.split('_vs_'))
    cond_label = _strip_common_affix(all_conds)

    def _encode_stripped(contrast):
        parts = contrast.split('_vs_')
        a = cond_label.get(parts[0], parts[0])[:12] if parts else contrast[:12]
        b = cond_label.get(parts[1], parts[1])[:12] if len(parts) > 1 else ""
        return f"{a}v{b}" if b else a

    def _encode_legacy(contrast):
        # ancien schéma (conditions complètes) — rétrocompat fichiers anciens
        return contrast.replace('_vs_', 'v').replace('_', '')

    # 1) Correspondance exacte sur le schéma strippé (nominal)
    for c in contrasts:
        enc = _encode_stripped(c)
        if enc[:len(s)] == s or enc == s:
            return c

    # 2) Correspondance exacte sur l'ancien schéma (anciens classeurs)
    for c in contrasts:
        for enc in (_encode_legacy(c),
                    c.replace('_vs_', 'v')):   # ancien format gardant les '_'
            if enc[:len(s)] == s or enc == s:
                return c

    # 3) Repli : meilleur préfixe commun (strippé), avec seuil de fiabilité
    def _lcp(a, b):
        n = 0
        for x, y in zip(a, b):
            if x == y:
                n += 1
            else:
                break
        return n
    best, best_len = None, -1
    for c in contrasts:
        n = max(_lcp(_encode_stripped(c), s), _lcp(_encode_legacy(c), s))
        if n > best_len:
            best_len, best = n, c
    return best if best_len >= max(3, len(s) // 2) else None


def read_go_table(path, sheet):
    """Lit une sheet GO_* : la table R commence colonne 21 -> on retire les colonnes vides."""
    df = pd.read_excel(path, sheet_name=sheet)
    df = df[[c for c in df.columns if not str(c).startswith('Unnamed')]]
    if 'term_name' in df.columns:
        df = df.dropna(subset=['term_name'])
    return df


def build_go(path, contrasts, df_comp):
    go_data, chord_data = {}, {}
    if not path:
        return go_data, chord_data
    xl = pd.ExcelFile(path)
    go_sheets = [s for s in xl.sheet_names if s.startswith('GO_') or s.lower().startswith('go')]
    # mapping sheet -> contraste
    sheet_map = {s: map_go_sheet(s, contrasts) for s in go_sheets}

    # lookup accession -> (gene, lfc par contraste) pour chord
    name_to_gene = {}
    if 'name' in df_comp.columns:
        gcol = df_comp['Genes'] if 'Genes' in df_comp else pd.Series([''] * len(df_comp))
        name_to_gene = dict(zip(df_comp['name'].astype(str), gcol.fillna('').astype(str)))

    for sheet in go_sheets:
        c = sheet_map.get(sheet)
        if not c:
            continue
        df = read_go_table(path, sheet)
        keep = [k for k in ['source', 'term_name', 'term_id', 'p_value',
                            'intersection_size', 'query_size', 'z_score', 'gene_ratio']
                if k in df.columns]
        recs = []
        for _, r in df[keep].iterrows():
            rec = {
                'source': safe_str(r.get('source')),
                'term_name': safe_str(r.get('term_name')),
                'term_id': safe_str(r.get('term_id')),
                'p_value': round(float(r['p_value']), 5) if pd.notna(r.get('p_value')) else None,
                'intersection_size': int(r['intersection_size']) if pd.notna(r.get('intersection_size')) else 0,
                'query_size': int(r['query_size']) if pd.notna(r.get('query_size')) else 0,
                'z_score': round(float(r['z_score']), 5) if pd.notna(r.get('z_score')) else 0.0,
                'gene_ratio': round(float(r['gene_ratio']), 5) if pd.notna(r.get('gene_ratio')) else 0.0,
            }
            recs.append(rec)
        go_data[c] = recs

        # CHORD : top termes GO:BP avec colonne intersection (accessions) -> LFC
        if 'intersection' in df.columns:
            diff_col = c + '_diff'
            lfc_lookup = {}
            if diff_col in df_comp.columns and 'name' in df_comp.columns:
                lfc_lookup = dict(zip(df_comp['name'].astype(str),
                                      df_comp[diff_col]))
            bp = df[df['source'] == 'GO:BP'] if 'source' in df.columns else df
            bp = bp.dropna(subset=['intersection'])
            if 'p_value' in bp.columns:
                bp = bp.sort_values('p_value').head(10)
            else:
                bp = bp.head(10)
            chord_terms = []
            for _, r in bp.iterrows():
                accs = [a.strip() for a in str(r['intersection']).split(',') if a.strip()]
                genes = {}
                for a in accs:
                    lfc = lfc_lookup.get(a)
                    if lfc is None or pd.isna(lfc):
                        continue
                    label = name_to_gene.get(a) or a
                    genes[label] = round(float(lfc), 3)
                if genes:
                    chord_terms.append({
                        'term': safe_str(r.get('term_name')),
                        'p': round(float(r['p_value']), 5) if pd.notna(r.get('p_value')) else None,
                        'z': round(float(r['z_score']), 3) if pd.notna(r.get('z_score')) else 0.0,
                        'genes': genes,
                    })
            if chord_terms:
                chord_data[c] = chord_terms
    return go_data, chord_data



# ──────────────────────────────────────────────────────────────────────────────
# GO BUBBLE PLOT — données structurées pour bubble plot + table interactive
# ──────────────────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────────────────
# CORRÉLATION MODULE-TRAIT WGCNA
# ──────────────────────────────────────────────────────────────────────────────

def build_module_trait(wgcna_path, df_comp, sample_cols, cond_map):
    """Calcule la matrice de corrélation Pearson eigengène × condition.

    Eigengène = moyenne des z-scores par réplique pour chaque module.
    Retourne : {module: {corr: {cond: r}, pval: {cond: p}, n: int}}
    """
    if not wgcna_path or not os.path.isfile(wgcna_path):
        return {}
    try:
        from scipy import stats as _stats
    except ImportError:
        print("  [warn] build_module_trait : scipy manquant.")
        return {}
    try:
        df_wgcna = pd.read_excel(wgcna_path)
        df_wgcna = df_wgcna[df_wgcna['Module_Color'].str.lower() != 'grey'].copy()
    except Exception as e:
        print(f"  [warn] build_module_trait : {e}")
        return {}

    df_merged = df_wgcna.merge(
        df_comp[['name'] + [c for c in sample_cols if c in df_comp.columns]]
              .rename(columns={'name': 'Gene_Name'}),
        on='Gene_Name', how='inner'
    )
    if df_merged.empty:
        return {}

    conditions = list(dict.fromkeys(cond_map.values()))
    modules    = sorted(df_merged['Module_Color'].unique())
    results    = {}

    for mod in modules:
        sub = df_merged[df_merged['Module_Color'] == mod][
            [c for c in sample_cols if c in df_merged.columns]
        ].fillna(0).values

        if len(sub) < 3:
            continue

        # Eigengène = moyenne des z-scores par protéine puis par réplique
        means = sub.mean(axis=1, keepdims=True)
        stds  = sub.std(axis=1,  keepdims=True) + 1e-8
        eigengene = ((sub - means) / stds).mean(axis=0)

        corrs, pvals = {}, {}
        for cond in conditions:
            trait = np.array(
                [1.0 if cond_map.get(c) == cond else 0.0
                 for c in sample_cols if c in df_merged.columns]
            )
            if trait.sum() == 0:
                continue
            r, p = _stats.pearsonr(eigengene, trait)
            corrs[cond] = round(float(r), 4)
            pvals[cond] = round(float(p), 5)

        results[str(mod)] = {'corr': corrs, 'pval': pvals, 'n': int(len(sub))}

    print(f"  build_module_trait : {len(results)} modules × {len(conditions)} conditions.")
    return results


# ──────────────────────────────────────────────────────────────────────────────
# WGCNA HUB MAP — lookup accession → {mod, hub} pour le volcano
# ──────────────────────────────────────────────────────────────────────────────
HUB_SCORE_THRESHOLD = 0.85  # Top ~15% des protéines par module

def build_wgcna_hub_map(wgcna_path, hub_threshold=HUB_SCORE_THRESHOLD):
    """Retourne un dict {accession: {mod, hub}} pour les hubs au-dessus du seuil."""
    if not wgcna_path or not os.path.isfile(wgcna_path):
        return {}
    try:
        df = pd.read_excel(wgcna_path)
        df = df[df['Module_Color'].str.lower() != 'grey']
        df = df[df['HubScore'] >= hub_threshold]
        hub_map = {}
        for _, r in df.iterrows():
            acc = safe_str(r['Gene_Name'])
            if acc:
                hub_map[acc] = {
                    'mod': safe_str(r['Module_Color']),
                    'hub': round(float(r['HubScore']), 4),
                }
        print(f"  build_wgcna_hub_map : {len(hub_map)} hubs (seuil={hub_threshold}).")
        return hub_map
    except Exception as e:
        print(f"  [warn] build_wgcna_hub_map : {e}")
        return {}

def build_go_bubble(path, contrasts):
    """Lit toutes les sheets GO et retourne un dict prêt pour le bubble plot JS.
    Structure : {contraste: [{term_name, term_id, source, p_value, neg_log10_p,
                              z_score, gene_ratio, intersection_size, intersection}]}
    """
    if not path or not os.path.isfile(path):
        return {}
    try:
        xl = pd.ExcelFile(path)
    except Exception:
        return {}

    bubble = {}
    go_sheets = [s for s in xl.sheet_names if 'go' in s.lower() or s.startswith('GO_')]

    for sheet in go_sheets:
        c = map_go_sheet(sheet, contrasts)
        if not c:
            continue
        try:
            df = pd.read_excel(path, sheet_name=sheet)
            df = df[[col for col in df.columns if not str(col).startswith('Unnamed')]]
        except Exception:
            continue

        keep = ['term_name','term_id','source','p_value','z_score',
                'gene_ratio','intersection_size','intersection']
        keep = [k for k in keep if k in df.columns]
        df = df[keep].dropna(subset=['term_name','p_value'])

        recs = []
        for _, r in df.iterrows():
            p = float(r['p_value'])
            recs.append({
                'term':   safe_str(r.get('term_name'))[:80],
                'tid':    safe_str(r.get('term_id')),
                'src':    safe_str(r.get('source','?')),
                'p':      round(p, 6),
                'lp':     round(float(-np.log10(max(p, 1e-10))), 3),
                'z':      round(float(r.get('z_score', 0) or 0), 4),
                'gr':     round(float(r.get('gene_ratio', 0) or 0), 4),
                'is':     int(r.get('intersection_size', 0) or 0),
                'genes':  safe_str(r.get('intersection', ''))[:500],
            })
        if recs:
            bubble[c] = recs

    return bubble


def build_go_heatmap(path, contrasts):
    """Construit la matrice z_score termes × contrastes pour la heatmap GO.
    Structure : {
      terms   : [{tid, name, src, row: [z|null × n_contrasts]}],
      contrasts: [str],
      z_min, z_max
    }
    Termes triés par source (BP/MF/CC) puis par variance décroissante.
    """
    if not path or not os.path.isfile(path):
        return {'terms': [], 'contrasts': contrasts, 'labels': contrasts, 'z_min': -1, 'z_max': 1}
    try:
        xl = pd.ExcelFile(path)
    except Exception:
        return {'terms': [], 'contrasts': contrasts, 'labels': contrasts, 'z_min': -1, 'z_max': 1}

    go_sheets = [s for s in xl.sheet_names if 'go' in s.lower() or s.startswith('GO_')]
    term_dict = {}  # tid -> {name, src, data: {contrast: z}}

    for sheet in go_sheets:
        c = map_go_sheet(sheet, contrasts)
        if not c:
            continue
        try:
            df = pd.read_excel(path, sheet_name=sheet)
            df = df[[col for col in df.columns if not str(col).startswith('Unnamed')]]
        except Exception:
            continue
        if 'term_id' not in df.columns or 'z_score' not in df.columns:
            continue
        for _, r in df.dropna(subset=['term_id', 'z_score']).iterrows():
            tid = safe_str(r['term_id'])
            if tid not in term_dict:
                term_dict[tid] = {
                    'name': safe_str(r.get('term_name', tid))[:60],
                    'src':  safe_str(r.get('source', '?')),
                    'data': {}
                }
            term_dict[tid]['data'][c] = round(float(r['z_score']), 4)

    if not term_dict:
        return {'terms': [], 'contrasts': contrasts, 'labels': contrasts, 'z_min': -1, 'z_max': 1}

    # Trier : source BP > MF > CC > autre, puis variance décroissante
    src_order = {'GO:BP': 0, 'GO:MF': 1, 'GO:CC': 2}
    def sort_key(item):
        tid, v = item
        zvals = list(v['data'].values())
        variance = float(np.var(zvals)) if len(zvals) > 1 else 0.0
        return (src_order.get(v['src'], 3), -variance)

    sorted_terms = sorted(term_dict.items(), key=sort_key)

    # Construire les lignes de la matrice
    all_z = []
    terms_out = []
    for tid, v in sorted_terms:
        row = [round(v['data'].get(c), 4) if c in v['data'] else None
               for c in contrasts]
        zvals = [z for z in row if z is not None]
        all_z.extend(zvals)
        terms_out.append({
            'tid':  tid,
            'name': v['name'],
            'src':  v['src'],
            'row':  row,
            'n':    len(zvals),   # nombre de contrastes où le terme apparaît
        })

    z_abs_max = max(abs(z) for z in all_z) if all_z else 1.0

    # Labels courts lisibles pour l'axe colonnes : on retire le préfixe/suffixe
    # commun aux conditions (évite 'X18009IG01Ctrl/X18009IG01di6h' illisible).
    all_conds = []
    for c in contrasts:
        all_conds.extend(c.split('_vs_'))
    cond_label = _strip_common_affix(all_conds)
    contrast_labels = []
    for c in contrasts:
        parts = c.split('_vs_')
        a = cond_label.get(parts[0], parts[0][:8]) if parts else c
        b = cond_label.get(parts[1], parts[1][:8]) if len(parts) > 1 else ''
        contrast_labels.append(f"{a}/{b}" if b else a)

    return {
        'terms':     terms_out,
        'contrasts': contrasts,          # noms complets (tooltips / mapping)
        'labels':    contrast_labels,    # labels courts lisibles (affichage)
        'z_min':     round(-z_abs_max, 3),
        'z_max':     round(z_abs_max, 3),
    }


# ──────────────────────────────────────────────────────────────────────────────
# RÉSEAU GO HIÉRARCHIQUE — graphe DAG termes enrichis avec relations parent/enfant
# ──────────────────────────────────────────────────────────────────────────────

def build_go_network(path, contrasts):
    """Construit le graphe hiérarchique GO pour chaque contraste.
    Structure : {contraste: {nodes: [{id,name,src,p,z,is,level}], edges: [{from,to}]}}
    Level = profondeur dans le DAG (0 = racine des termes enrichis).
    """
    if not path or not os.path.isfile(path):
        return {}
    try:
        xl = pd.ExcelFile(path)
    except Exception:
        return {}

    go_sheets = [s for s in xl.sheet_names if 'go' in s.lower() or s.startswith('GO_')]
    result = {}

    for sheet in go_sheets:
        c = map_go_sheet(sheet, contrasts)
        if not c:
            continue
        try:
            df = pd.read_excel(path, sheet_name=sheet)
            df = df[[col for col in df.columns if not str(col).startswith('Unnamed')]]
        except Exception:
            continue
        if 'term_id' not in df.columns or 'parents' not in df.columns:
            continue
        df = df.dropna(subset=['term_id'])

        # Index des termes enrichis
        term_ids = set(df['term_id'].astype(str))

        # Nœuds
        nodes = {}
        for _, row in df.iterrows():
            tid = safe_str(row['term_id'])
            nodes[tid] = {
                'id':   tid,
                'name': safe_str(row.get('term_name', tid))[:60],
                'src':  safe_str(row.get('source', '?')),
                'p':    round(float(row['p_value']), 6) if pd.notna(row.get('p_value')) else 1.0,
                'z':    round(float(row['z_score']), 4)  if pd.notna(row.get('z_score'))  else 0.0,
                'is':   int(row['intersection_size'])    if pd.notna(row.get('intersection_size')) else 0,
                'level': 0,  # sera calculé
            }

        # Arêtes (seulement intra-enrichis)
        edges = []
        children = {tid: [] for tid in nodes}
        parents_of = {tid: [] for tid in nodes}
        for _, row in df.iterrows():
            child_id = safe_str(row['term_id'])
            parents_str = safe_str(row.get('parents', ''))
            plist = [p.strip() for p in parents_str.strip('[]\'\" ').split(',')
                     if p.strip().startswith('GO:') and p.strip() in term_ids]
            for p in plist:
                edges.append({'from': p, 'to': child_id})
                children[p].append(child_id)
                parents_of[child_id].append(p)

        # Calcul des niveaux (BFS depuis les racines)
        roots = [tid for tid in nodes if not parents_of[tid]]
        from collections import deque
        q = deque([(r, 0) for r in roots])
        visited = set()
        while q:
            tid, lvl = q.popleft()
            if tid in visited:
                continue
            visited.add(tid)
            nodes[tid]['level'] = lvl
            for ch in children.get(tid, []):
                q.append((ch, lvl + 1))

        result[c] = {
            'nodes': list(nodes.values()),
            'edges': edges,
        }

    return result

def build_go_images(path, contrasts):
    """Extrait les images embarquées dans chaque sheet GO (ordre : manhattan, lollipop, dotplot, chord)."""
    go_images = {}
    if not path:
        return go_images
    try:
        from PIL import Image
    except Exception as e:
        print("  [!] openpyxl/Pillow indisponibles, images GO ignorées (%s)" % e)
        return go_images

    wb = openpyxl.load_workbook(path)
    for sheet in wb.sheetnames:
        if not (sheet.startswith('GO_') or sheet.lower().startswith('go')):
            continue
        c = map_go_sheet(sheet, contrasts)
        if not c:
            continue
        ws = wb[sheet]
        imgs = getattr(ws, '_images', [])
        bucket = {}
        for i, img in enumerate(imgs[:4]):
            try:
                data = img._data()
                pil = Image.open(io.BytesIO(data)).convert('RGB')
                if pil.width > 1200:
                    pil = pil.resize((1200, int(pil.height * 1200 / pil.width)), Image.LANCZOS)
                buf = io.BytesIO()
                pil.save(buf, format='JPEG', quality=85)
                b64 = 'data:image/jpeg;base64,' + base64.b64encode(buf.getvalue()).decode()
                bucket[IMG_LABELS[i]] = b64
            except Exception:
                continue
        if bucket:
            go_images[c] = bucket
    return go_images


# ──────────────────────────────────────────────────────────────────────────────
# WGCNA
# ──────────────────────────────────────────────────────────────────────────────
def build_wgcna(path):
    rows, top5, mods = [], {}, []
    if not path:
        return rows, top5, mods
    xl = pd.ExcelFile(path)
    sheet = next((s for s in xl.sheet_names if 'wgcna' in s.lower()), xl.sheet_names[0])
    df = pd.read_excel(path, sheet_name=sheet)

    gcol = 'Gene_Name' if 'Gene_Name' in df.columns else ('name' if 'name' in df.columns else None)
    mcol = 'Module_Color' if 'Module_Color' in df.columns else None
    hcol = 'HubScore' if 'HubScore' in df.columns else None
    dcol = next((c for c in ['First.Protein.Description', 'Protein.Names'] if c in df.columns), None)
    if not (gcol and mcol and hcol):
        return rows, top5, mods

    df = df[df[mcol].astype(str).str.lower() != 'grey']      # exclure module grey
    df = df.sort_values([mcol, hcol], ascending=[True, False])

    for _, r in df.iterrows():
        rows.append({
            'g': safe_str(r.get(gcol)),
            'mod': safe_str(r.get(mcol)),
            'hub': round(float(r[hcol]), 4) if pd.notna(r.get(hcol)) else 0.0,
            'desc': safe_str(r.get(dcol)) if dcol else '',
        })
    mods = sorted(df[mcol].dropna().astype(str).unique().tolist())
    for m in mods:
        top5[m] = [rec for rec in rows if rec['mod'] == m][:5]
    return rows, top5, mods


# ──────────────────────────────────────────────────────────────────────────────
# UMAP — lecture des coordonnées pré-calculées par R (sheet UMAP_Output)
# ──────────────────────────────────────────────────────────────────────────────

def build_umap(path):
    """Lit la sheet UMAP_Output (colonnes UMAP1, UMAP2, Sample, Condition).
    Retourne un dict prêt pour injection JS :
      points  : [{u1, u2, sample, cond}]
      conditions : liste ordonnée des conditions uniques
    Robuste : tolère les noms de colonnes avec casse variable et les espaces.
    """
    if not path or not os.path.isfile(path):
        return {'points': [], 'conditions': []}
    try:
        xl = pd.ExcelFile(path)
        # Prioriser une sheet nommée EXACTEMENT 'UMAP' (PCA_UMAP ne contient
        # que des images), sinon toute sheet contenant 'umap'.
        exact = [s for s in xl.sheet_names if s.strip().lower() == 'umap']
        if exact:
            sheet = exact[0]
        else:
            sheet = next(
                (s for s in xl.sheet_names
                 if 'umap' in s.lower() and 'pca' not in s.lower()),
                None
            )
        if sheet is None:
            return {'points': [], 'conditions': []}
        df = pd.read_excel(path, sheet_name=sheet)
        # Normaliser les noms de colonnes
        df.columns = [str(c).strip() for c in df.columns]
        col_map = {c.lower(): c for c in df.columns}
        u1_col   = col_map.get('umap1') or col_map.get('umap_1') or col_map.get('dim1')
        u2_col   = col_map.get('umap2') or col_map.get('umap_2') or col_map.get('dim2')
        samp_col = (col_map.get('sample') or col_map.get('sample_id')
                    or col_map.get('label') or col_map.get('id'))
        cond_col = col_map.get('condition') or col_map.get('group') or col_map.get('cond')
        if not (u1_col and u2_col):
            print("  [warn] build_umap : colonnes UMAP1/UMAP2 introuvables dans '%s'." % sheet)
            return {'points': [], 'conditions': []}
        df = df.dropna(subset=[u1_col, u2_col])
        conditions = list(dict.fromkeys(
            df[cond_col].astype(str).tolist() if cond_col else ['Unknown'] * len(df)
        ))
        points = []
        for _, row in df.iterrows():
            points.append({
                'u1':    round(float(row[u1_col]), 5),
                'u2':    round(float(row[u2_col]), 5),
                'sample': safe_str(row[samp_col]) if samp_col else '',
                'cond':   safe_str(row[cond_col]) if cond_col else 'Unknown',
            })
        print("  build_umap : %d points, %d conditions." % (len(points), len(conditions)))
        return {'points': points, 'conditions': conditions}
    except Exception as e:
        print("  [warn] build_umap : erreur lecture (%s)." % e)
        return {'points': [], 'conditions': []}


# ──────────────────────────────────────────────────────────────────────────────
# WGCNA NETWORK — graphe étoile hub-centric (sans matrice TOM)
# ──────────────────────────────────────────────────────────────────────────────
WGCNA_SPOKE_N = 20   # top-N protéines par module (hors hub)

def build_wgcna_network(path):
    """Construit un graphe étoile par module WGCNA depuis la sheet WGCNA_Results.

    Structure retournée : {module_color: {hub, nodes, n_total, h_max, h_min}}
      hub   : {g, hub, desc}          — protéine avec le HubScore le plus élevé
      nodes : [{g, hub, desc, w}]     — top WGCNA_SPOKE_N, w = épaisseur arête normalisée
      n_total : int                   — total protéines dans le module (y.c. hors top-N)
    Les edges sont implicites (hub → chaque node). Pas de matrice TOM disponible.
    """
    if not path or not os.path.isfile(path):
        return {}
    try:
        xl = pd.ExcelFile(path)
        sheet = next((s for s in xl.sheet_names if 'wgcna' in s.lower()), xl.sheet_names[0])
        df = pd.read_excel(path, sheet_name=sheet)

        gcol  = 'Gene_Name'      if 'Gene_Name'      in df.columns else 'name'
        mcol  = 'Module_Color'   if 'Module_Color'   in df.columns else None
        hcol  = 'HubScore'       if 'HubScore'       in df.columns else None
        dcol  = next((c for c in ['First.Protein.Description', 'Protein.Names']
                      if c in df.columns), None)
        if not (gcol in df.columns and mcol and hcol):
            print("  [warn] build_wgcna_network : colonnes manquantes.")
            return {}

        df = df[df[mcol].astype(str).str.lower() != 'grey'].copy()
        network = {}
        for mod in sorted(df[mcol].dropna().unique()):
            grp = df[df[mcol] == mod].sort_values(hcol, ascending=False)
            if len(grp) < 2:
                continue
            h_max   = float(grp[hcol].max())
            h_min   = float(grp[hcol].min())
            h_range = h_max - h_min if h_max > h_min else 1.0

            hub_row = grp.iloc[0]
            spokes  = grp.iloc[1:WGCNA_SPOKE_N + 1]

            network[str(mod)] = {
                'hub': {
                    'g':    safe_str(hub_row[gcol]),
                    'hub':  round(float(hub_row[hcol]), 4),
                    'desc': safe_str(hub_row[dcol])[:60] if dcol else '',
                },
                'nodes': [
                    {
                        'g':    safe_str(r[gcol]),
                        'hub':  round(float(r[hcol]), 4),
                        'desc': safe_str(r[dcol])[:60] if dcol else '',
                        'w':    round(0.2 + 0.8 * (float(r[hcol]) - h_min) / h_range, 3),
                    }
                    for _, r in spokes.iterrows()
                ],
                'n_total': int(len(grp)),
                'h_max':   round(h_max, 4),
                'h_min':   round(h_min, 4),
            }
        print("  build_wgcna_network : %d modules." % len(network))
        return network
    except Exception as e:
        print("  [warn] build_wgcna_network : erreur (%s)." % e)
        return {}


# ──────────────────────────────────────────────────────────────────────────────
# VIOLIN PLOT — distributions log2 LFQ post-imputation par échantillon
# ──────────────────────────────────────────────────────────────────────────────
VIOLIN_N_BINS = 60   # résolution de l'histogramme (KDE approchée côté JS)

def build_violin(df_comp, sample_cols, cond_map):
    """Retourne pour chaque échantillon : histogramme à VIOLIN_N_BINS bins,
    quartiles [q0,q1,q2,q3,q4] et condition d'appartenance.
    La plage globale (vmin/vmax) est calculée sur l'ensemble des échantillons
    pour que les bins soient comparables entre violins.
    """
    all_vals = df_comp[sample_cols].values.flatten()
    all_vals = all_vals[~np.isnan(all_vals)]
    vmin, vmax = float(all_vals.min()), float(all_vals.max())
    bins = np.linspace(vmin, vmax, VIOLIN_N_BINS + 1)
    bin_centers = ((bins[:-1] + bins[1:]) / 2).round(3).tolist()

    violin = {}
    for samp in sample_cols:
        vals = df_comp[samp].dropna().values
        hist, _ = np.histogram(vals, bins=bins)
        q = np.percentile(vals, [0, 25, 50, 75, 100]).round(3).tolist()
        violin[samp] = {
            'hist': hist.tolist(),
            'q': q,                        # [min, Q1, median, Q3, max]
            'cond': cond_map.get(samp, samp.rsplit('_', 1)[0]),
            'n': int(len(vals)),
        }
    return {
        'samples': violin,
        'bin_centers': bin_centers,
        'vmin': round(vmin, 3),
        'vmax': round(vmax, 3),
    }


# ──────────────────────────────────────────────────────────────────────────────
# MISSING VALUES — matrice présence/absence depuis raw_data (avant imputation)
# ──────────────────────────────────────────────────────────────────────────────
MV_MAX_PROT = 600    # nombre max de protéines affichées (triées par NA rate desc)

def build_missing_values(stats_path, sample_cols, cond_map):
    """Lit la sheet raw_data pour obtenir les valeurs brutes avant imputation.
    Retourne :
      - matrix : liste de {pg, na_vec} où na_vec = [0|1] par échantillon brut
      - lfq_cols : noms de colonnes LFQ brutes (ordre identique à na_vec)
      - lfq_cond : condition de chaque colonne LFQ (mapping par position)
      - na_rates : taux NA par colonne LFQ (pour la légende)
      - summary : stats globales
    """
    try:
        df_brut = pd.read_excel(stats_path, sheet_name='raw_data')
    except Exception as e:
        print("  [warn] build_missing_values : sheet raw_data inaccessible (%s)" % e)
        return {'matrix': [], 'lfq_cols': [], 'lfq_cond': [], 'na_rates': [], 'summary': {}}

    lfq_cols = [c for c in df_brut.columns if 'LFQ' in str(c)]
    if not lfq_cols:
        return {'matrix': [], 'lfq_cols': [], 'lfq_cond': [], 'na_rates': [], 'summary': {}}

    # Matching colonnes LFQ brutes → conditions via sample_cols post-imputation
    # Stratégie : normaliser les deux noms et chercher la meilleure correspondance
    def _norm_col(s):
        return re.sub(r'[^a-z0-9]', '', str(s).lower())

    norm_to_cond = {_norm_col(c): cond_map.get(c, c.rsplit('_', 1)[0])
                    for c in sample_cols}

    # Mapping LFQ brutes → conditions par position
    # Les n colonnes LFQ de raw_data correspondent dans l'ordre aux n sample_cols
    # post-imputation (même design, même ordre garanti par le script R)
    if len(lfq_cols) == len(sample_cols):
        lfq_cond = [cond_map.get(sample_cols[i], 'Unknown') for i in range(len(lfq_cols))]
    else:
        # Fallback : mapper par index modulo si longueurs différentes (réplicats techniques A/B)
        # raw_data peut avoir 2x plus de colonnes (A+B) que les sample_cols post-imputation
        cond_list = [cond_map.get(c, 'Unknown') for c in sample_cols]
        ratio = len(lfq_cols) // len(sample_cols) if len(sample_cols) > 0 else 1
        if ratio >= 1:
            lfq_cond = []
            for i in range(len(lfq_cols)):
                idx = i // max(ratio, 1)
                lfq_cond.append(cond_list[idx] if idx < len(cond_list) else 'Unknown')
        else:
            lfq_cond = ['Unknown'] * len(lfq_cols)

    # Matrice présence/absence (1 = manquant)
    na_matrix = df_brut[lfq_cols].isna().astype(int)
    na_per_prot = na_matrix.sum(axis=1)

    # Trier par taux NA décroissant, garder les MV_MAX_PROT premières avec NA > 0
    # + un échantillon aléatoire de protéines sans NA pour contexte
    has_na_idx = na_per_prot[na_per_prot > 0].sort_values(ascending=False).index
    no_na_idx  = na_per_prot[na_per_prot == 0].sample(
        min(100, (na_per_prot == 0).sum()), random_state=SEED).index
    keep_idx   = list(has_na_idx[:MV_MAX_PROT]) + list(no_na_idx)
    df_keep    = df_brut.loc[keep_idx]
    na_keep    = na_matrix.loc[keep_idx]

    matrix = []
    for idx in keep_idx:
        pg  = safe_str(df_brut.at[idx, 'Protein.Group']) if 'Protein.Group' in df_brut else ''
        gene = safe_str(df_brut.at[idx, 'Genes']) if 'Genes' in df_brut else ''
        vec = na_keep.loc[idx].tolist()    # [0|1] × n_lfq_cols
        matrix.append({'pg': pg, 'g': gene, 'v': vec})

    na_rates = na_matrix[lfq_cols].mean().round(4).tolist()

    summary = {
        'total_prot': int(len(df_brut)),
        'prot_with_na': int((na_per_prot > 0).sum()),
        'global_na_rate': round(float(na_matrix.values.mean()), 4),
        'n_lfq_cols': len(lfq_cols),
    }

    return {
        'matrix': matrix,
        'lfq_cols': [c.replace('LFQ.intensity.', '') for c in lfq_cols],  # labels courts
        'lfq_cond': lfq_cond,
        'na_rates': na_rates,
        'summary': summary,
    }


# ──────────────────────────────────────────────────────────────────────────────
# MA PLOT — A = intensité moyenne log2 LFQ, M = LFC (diff)
# ──────────────────────────────────────────────────────────────────────────────
NS_MA_MAX = 600     # NS sous-échantillonnés (même logique que volcano)

def build_ma_plot(df_comp, contrasts, cond_map):
    """Pour chaque contraste, construit les points du MA plot :
      A = (mean_log2_A + mean_log2_B) / 2   (axe X : intensité moyenne)
      M = _diff                               (axe Y : LFC, déjà en log2)
    Statut UP/DOWN/NS dérivé de p.val et diff (identique au volcano).
    Sous-échantillonnage NS pour alléger le HTML.
    """
    # cond_map ici est cond_map_inv : {col -> condition}
    cond_to_cols = {}
    for col, cond in cond_map.items():
        cond_to_cols.setdefault(cond, []).append(col)

    def _find_cols(name):
        n = name.lower().replace('_', '').replace('.', '')
        for cond, cols in cond_to_cols.items():
            if cond.lower().replace('_', '').replace('.', '') == n:
                return [c for c in cols if c in df_comp.columns]
        return []

    ma = {}
    for c in contrasts:
        parts = c.split('_vs_')
        if len(parts) != 2:
            continue
        cols_a = _find_cols(parts[0])
        cols_b = _find_cols(parts[1])
        if not cols_a or not cols_b:
            continue

        mean_a = df_comp[cols_a].fillna(0).mean(axis=1)
        mean_b = df_comp[cols_b].fillna(0).mean(axis=1)
        A = ((mean_a + mean_b) / 2).round(3)

        diff_col = c + '_diff'
        pval_col = c + '_p.val'
        padj_col = c + '_p.adj'
        if diff_col not in df_comp.columns or pval_col not in df_comp.columns:
            continue

        M    = df_comp[diff_col].round(3)
        pval = df_comp[pval_col].fillna(1)
        padj = df_comp[padj_col].fillna(1) if padj_col in df_comp.columns else pval
        gene = df_comp['Genes'].fillna('') if 'Genes' in df_comp.columns else pd.Series([''] * len(df_comp))
        name = df_comp['name'] if 'name' in df_comp.columns else pd.Series([''] * len(df_comp))

        status = pd.Series(['N'] * len(df_comp), index=df_comp.index)
        status[(pval < P_THRESH) & (M > LFC_THRESH)]  = 'U'
        status[(pval < P_THRESH) & (M < -LFC_THRESH)] = 'D'

        df_ma = pd.DataFrame({
            'a': A, 'm': M, 's': status,
            'g': gene.values, 'n': name.values,
            'p': pval.round(6).values, 'pa': padj.round(6).values,
        })
        df_ma = df_ma.dropna(subset=['a', 'm'])

        sig   = df_ma[df_ma['s'] != 'N']
        ns_all = df_ma[df_ma['s'] == 'N']
        ns     = ns_all.sample(min(NS_MA_MAX, len(ns_all)), random_state=SEED) if len(ns_all) else ns_all
        pts    = pd.concat([sig, ns]).where(pd.notna, None)

        ma[c] = {
            'pts': json.loads(pts.to_json(orient='records')),
            'label_a': parts[0].replace('_', ' '),
            'label_b': parts[1].replace('_', ' '),
            'lfc_thresh': LFC_THRESH,
        }
    return ma


# ──────────────────────────────────────────────────────────────────────────────
# SCATTER INTER-RÉPLICATS — matrice de corrélations Pearson par paires intra-condition
# ──────────────────────────────────────────────────────────────────────────────
NS_REP_MAX = 120    # NS sous-échantillonnés par paire scatter (allégé : évite le gel navigateur)

def build_replicate_scatter(df_comp, cond_map):
    """Pour chaque condition, calcule toutes les paires de réplicats :
      - coordonnées (x, y) sous-échantillonnées
      - r de Pearson sur toutes les valeurs non-nulles
    Structure : {condition: {pairs: [{r1, r2, r, pts:[{x,y}]}], reps: [...]}}
    """
    from itertools import combinations as _comb

    cond_to_cols = {}
    for col, cond in cond_map.items():
        if col in df_comp.columns:
            cond_to_cols.setdefault(cond, []).append(col)

    rep_data = {}
    for cond, reps in cond_to_cols.items():
        if len(reps) < 2:
            continue
        pairs = []
        for r1, r2 in _comb(reps, 2):
            sub = df_comp[[r1, r2]].dropna()
            sub = sub[(sub[r1] > 0) & (sub[r2] > 0)]
            if len(sub) < 2:
                continue
            r = float(np.corrcoef(sub[r1], sub[r2])[0, 1])
            # sous-échantillonnage
            if len(sub) > NS_REP_MAX:
                sub = sub.sample(NS_REP_MAX, random_state=SEED)
            pts = [{'x': round(float(row[r1]), 3), 'y': round(float(row[r2]), 3)}
                   for _, row in sub.iterrows()]
            pairs.append({
                'r1': r1, 'r2': r2,
                'r':  round(r, 4),
                'n':  int(len(df_comp[[r1, r2]].dropna())),
                'pts': pts,
            })
        if pairs:
            rep_data[cond] = {
                'reps': reps,
                'pairs': pairs,
            }
    return rep_data


# ──────────────────────────────────────────────────────────────────────────────
# RANKED ABUNDANCE — intensité LFQ moyenne triée par rang décroissant
# ──────────────────────────────────────────────────────────────────────────────

def build_ranked_abundance(df_comp, sample_cols):
    """Ranked protein abundance — intensité LFQ moyenne triée par rang décroissant.
    Sous-échantillonne à N_RANKED points pour alléger le JSON (~40× plus léger).
    Structure conservée : {points:[{r,v,g?,n?}], vmin, vmax}
    """
    N_RANKED = 600   # suffisant pour une courbe lisse (~9945 → 600)
    intens = df_comp[sample_cols].fillna(0)
    mean_intens = intens.mean(axis=1)

    df_r = pd.DataFrame({
        'v': mean_intens,
        'g': df_comp.get('Genes', pd.Series([''] * len(df_comp))).fillna(''),
        'n': df_comp.get('name',  pd.Series([''] * len(df_comp))).fillna(''),
    })
    df_r = df_r[df_r['v'] > 0].sort_values('v', ascending=False).reset_index(drop=True)
    N_total = len(df_r)

    # Sous-échantillonner en conservant les top-5 et bottom-5
    keep_idx = set(range(5)) | set(range(N_total-5, N_total))
    step = max(1, N_total // N_RANKED)
    keep_idx |= set(range(0, N_total, step))
    keep_idx = sorted(keep_idx)

    points = []
    for i in keep_idx:
        row = df_r.iloc[i]
        pt = {'r': i + 1, 'v': round(float(row['v']), 3)}
        # Métadonnées seulement pour les extrêmes
        if i < 5 or i >= N_total - 5:
            pt['g'] = str(row['g'])[:12]
            pt['n'] = str(row['n'])
        points.append(pt)

    vmin = float(df_r['v'].min()) if len(df_r) else 0.0
    vmax = float(df_r['v'].max()) if len(df_r) else 1.0
    return {
        'points': points,
        'vmin': round(vmin, 3),
        'vmax': round(vmax, 3),
        'n_total': N_total,
    }

def build_cv_plot(df_comp, sample_cols, cond_map):
    """Pour chaque condition : scatter intensité log2 LFQ moyenne × CV% intra-condition.
    Sous-échantillonné à NS_CV points. Inclut la médiane de CV pour annotation.
    """
    cond_to_cols = {}
    for col, cond in cond_map.items():
        if col in df_comp.columns:
            cond_to_cols.setdefault(cond, []).append(col)

    cv_data = {}
    for cond, cols in cond_to_cols.items():
        if len(cols) < 2:
            continue
        sub = df_comp[cols].replace(0, np.nan)
        mean_c = sub.mean(axis=1)
        std_c  = sub.std(axis=1)
        cv     = (std_c / mean_c * 100).round(2)

        df_cv = pd.DataFrame({
            'x': mean_c.round(3),
            'y': cv,
            'g': df_comp['Genes'].fillna('') if 'Genes' in df_comp.columns else '',
        }).dropna()
        df_cv = df_cv[df_cv['x'] > 0]

        median_cv = round(float(df_cv['y'].median()), 2)

        if len(df_cv) > NS_CV:
            df_cv = df_cv.sample(NS_CV, random_state=SEED)

        cv_data[cond] = {
            'pts':    [{'x': float(r.x), 'y': float(r.y), 'g': str(r.g)}
                       for _, r in df_cv.iterrows()],
            'median': median_cv,
            'n':      int(len(df_cv)),
        }
    return cv_data


# ──────────────────────────────────────────────────────────────────────────────
# VERSIONING — panneau À propos avec métadonnées de génération
# ──────────────────────────────────────────────────────────────────────────────

def build_about(stats_path, go_path, wgcna_path, umap_path,
                contrasts, sample_cols, conditions, n_proteins,
                p_thresh, lfc_thresh, df_comp=None, chosen_params=None):
    """Construit le dict de métadonnées injecté dans le panneau À propos."""
    import hashlib as _hl
    from datetime import datetime as _dt

    def _file_meta(path):
        if not path or not os.path.isfile(path):
            return None
        size_kb = os.path.getsize(path) // 1024
        with open(path, 'rb') as f:
            md5 = _hl.md5(f.read()).hexdigest()[:8]
        return {'name': os.path.basename(path), 'size_kb': size_kb, 'md5': md5}

    # Inférer les seuils réels depuis les données (Robustness_Score)
    inferred_p = p_thresh
    inferred_lfc = lfc_thresh
    if df_comp is not None:
        p_vals_sig, lfc_vals_sig = [], []
        for c in contrasts:
            rcol = 'Robustness_Score_' + c
            pcol, dcol = c + '_p.val', c + '_diff'
            if rcol in df_comp.columns and pcol in df_comp.columns and dcol in df_comp.columns:
                sig = df_comp[df_comp[rcol] > 0]
                if len(sig) > 0:
                    p_vals_sig.append(float(sig[pcol].max()))
                    lfc_vals_sig.append(float(sig[dcol].abs().min()))
        if p_vals_sig:
            # Arrondir au seuil standard le plus proche (0.05, 0.01, 0.001)
            raw_p = float(np.median(p_vals_sig))
            for std in [0.001, 0.01, 0.05, 0.1]:
                if raw_p <= std * 10:
                    inferred_p = std
                    break
        if lfc_vals_sig:
            inferred_lfc = round(float(np.median(lfc_vals_sig)), 3)
    inferred_ratio = round(2 ** inferred_lfc, 3)

    # Si les seuils choisis dans le pipeline sont fournis, ils priment sur
    # l'inférence (qui reste un fallback).
    if chosen_params:
        v_padj  = chosen_params.get("volcano_use_padj", False)
        v_p     = chosen_params.get("volcano_p_thresh", inferred_p)
        v_ratio = chosen_params.get("volcano_ratio_min", inferred_ratio)
        v_lfc   = chosen_params.get("volcano_lfc_min", inferred_lfc)
        a_padj  = chosen_params.get("anova_use_padj", False)
        a_p     = chosen_params.get("anova_p_thresh", inferred_p)
        params_block = {
            'p_thresh':       v_p,
            'lfc_thresh':     v_lfc,
            'ratio_thresh':   v_ratio,
            'volcano_p_type': 'p.adj' if v_padj else 'p.value',
            'volcano_p':      v_p,
            'volcano_ratio':  v_ratio,
            'anova_p_type':   'p.adj' if a_padj else 'p.value',
            'anova_p':        a_p,
            'n_contrasts':    len(contrasts),
            'n_samples':      len(sample_cols),
            'n_conditions':   len(conditions),
            'n_proteins':     n_proteins,
        }
    else:
        params_block = {
            'p_thresh':     inferred_p,
            'lfc_thresh':   inferred_lfc,
            'ratio_thresh': inferred_ratio,
            'n_contrasts':  len(contrasts),
            'n_samples':    len(sample_cols),
            'n_conditions': len(conditions),
            'n_proteins':   n_proteins,
        }
    return {
        'generated_at': _dt.now().strftime('%Y-%m-%d %H:%M'),
        'files': {
            'stats': _file_meta(stats_path),
            'go':    _file_meta(go_path),
            'wgcna': _file_meta(wgcna_path),
            'umap':  _file_meta(umap_path),
        },
        'params': params_block,
        'contrasts': contrasts,
    }


# ──────────────────────────────────────────────────────────────────────────────
# Patch JS — Volcano comparatif multi-contrastes
# ──────────────────────────────────────────────────────────────────────────────

DRAW_VOLCANO_COMP_JS = """\
/* ── VOLCANO COMPARATIF ── */
(function(){

var _volcCompInit = false;

window.initVolcanoComp = function(){
  if(_volcCompInit) return;
  _volcCompInit = true;

  var volcPage = document.getElementById('volcano');
  if(!volcPage) return;

  var panel = document.createElement('div');
  panel.className = 'panel';
  panel.style.marginTop = '14px';
  panel.innerHTML =
    '<div class="ph">' +
    '<h2>Volcano comparatif (multi-contrastes)</h2>' +
    '<button id="volcCompExportBtn" class="exp-btn" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button>' +
    '</div>' +
    '<div id="volcCompChecks" style="display:flex;flex-wrap:wrap;gap:8px;margin-bottom:12px;"></div>' +
    '<canvas id="volcCompCanvas" style="display:block;width:100%;"></canvas>' +
    '<div id="volcCompTip" style="position:fixed;display:none;background:var(--bg2);border:1px solid var(--bd);' +
    'border-radius:6px;padding:6px 10px;font-size:11px;color:var(--tx);pointer-events:none;z-index:999;max-width:280px;"></div>';
  volcPage.appendChild(panel);
  // Attacher l'export après création pour éviter les problèmes d'apostrophes inline
  var exportBtn = document.getElementById('volcCompExportBtn');
  if(exportBtn) exportBtn.onclick = function(){ exportCanvas('volcCompCanvas','proteogen_volcano_comparatif'); };

  var COMP_COLORS = ['#388bfd','#3fb950','#d29922','#bc8cff','#f85149','#39d353'];
  var colorMap = {};
  CONTRASTS.forEach(function(c,i){ colorMap[c] = COMP_COLORS[i%COMP_COLORS.length]; });

  var checksDiv = document.getElementById('volcCompChecks');
  CONTRASTS.forEach(function(c, i){
    var col = COMP_COLORS[i % COMP_COLORS.length];
    var label = document.createElement('label');
    label.style.cssText = 'display:inline-flex;align-items:center;gap:5px;cursor:pointer;' +
      'font-size:11px;color:var(--tx2);padding:3px 8px;border-radius:4px;' +
      'border:1px solid var(--bd);background:var(--bg2);';
    label.innerHTML =
      '<input type="checkbox" value="'+c+'" checked style="accent-color:'+col+';cursor:pointer;">' +
      '<span style="color:'+col+';font-weight:600;">'+
      c.replace(/_vs_/,' vs ').replace(/_/g,' ')+'</span>';
    label.querySelector('input').addEventListener('change', drawVolcanoComp);
    checksDiv.appendChild(label);
  });

  drawVolcanoComp();
};

window.drawVolcanoComp = function(){
  var canvas = document.getElementById('volcCompCanvas');
  if(!canvas) return;
  var dpr = devicePixelRatio || 1;
  var lW=60, rW=110, tH=36, bH=44, PW=660, PH=380;
  var W = lW+PW+rW, H = tH+PH+bH;
  canvas.width = W*dpr; canvas.height = H*dpr;
  canvas.style.width = W+'px'; canvas.style.height = H+'px';
  var ctx = canvas.getContext('2d');
  ctx.scale(dpr, dpr);
  ctx.clearRect(0,0,W,H);

  var checks = document.querySelectorAll('#volcCompChecks input:checked');
  var active = [];
  checks.forEach(function(cb){ active.push(cb.value); });
  if(!active.length) return;

  var COMP_COLORS = ['#388bfd','#3fb950','#d29922','#bc8cff','#f85149','#39d353'];
  var colorMap = {};
  CONTRASTS.forEach(function(c,i){ colorMap[c] = COMP_COLORS[i%COMP_COLORS.length]; });

  // Plages globales
  var allD = [], allLP = [];
  active.forEach(function(c){
    (VOLCANO[c]||[]).forEach(function(p){
      if(p.d!=null) allD.push(p.d);
      if(p.lp!=null) allLP.push(p.lp);
    });
  });
  if(!allD.length) return;

  var dMin = Math.min.apply(null,allD)*1.05;
  var dMax = Math.max.apply(null,allD)*1.05;
  var lpMax = Math.max.apply(null,allLP)*1.1 || 4;

  function toX(d){ return lW + (d-dMin)/(dMax-dMin)*PW; }
  function toY(lp){ return tH + PH*(1 - lp/lpMax); }

  // Grille
  ctx.strokeStyle='#21262d'; ctx.lineWidth=1;
  [-4,-3,-2,-1,0,1,2,3,4,5].forEach(function(v){
    if(v<dMin||v>dMax) return;
    var x=toX(v);
    ctx.beginPath(); ctx.moveTo(x,tH); ctx.lineTo(x,tH+PH); ctx.stroke();
    ctx.fillStyle='#8b949e'; ctx.font='8px sans-serif'; ctx.textAlign='center';
    ctx.fillText(v, x, tH+PH+14);
  });
  [1,2,3,4,5,6,7].forEach(function(v){
    if(v>lpMax) return;
    var y=toY(v);
    ctx.beginPath(); ctx.moveTo(lW,y); ctx.lineTo(lW+PW,y); ctx.stroke();
    ctx.fillStyle='#8b949e'; ctx.font='8px sans-serif'; ctx.textAlign='right';
    ctx.fillText(v, lW-4, y+3);
  });

  // Ligne LFC=0 + seuil p
  ctx.strokeStyle='#8b949e'; ctx.lineWidth=1.5;
  ctx.beginPath(); ctx.moveTo(toX(0),tH); ctx.lineTo(toX(0),tH+PH); ctx.stroke();
  ctx.setLineDash([5,4]); ctx.strokeStyle='#6e7681'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(lW,toY(-Math.log10(0.05))); ctx.lineTo(lW+PW,toY(-Math.log10(0.05))); ctx.stroke();
  ctx.setLineDash([]);

  // NS d'abord (gris transparent)
  active.forEach(function(c){
    var col = colorMap[c];
    (VOLCANO[c]||[]).forEach(function(p){
      if(p.s!=='N'||p.d==null||p.lp==null) return;
      ctx.beginPath(); ctx.arc(toX(p.d), toY(p.lp), 1.5, 0, 2*Math.PI);
      ctx.fillStyle = col+'28'; ctx.fill();
    });
  });
  // Significatifs par-dessus
  active.forEach(function(c){
    var col = colorMap[c];
    (VOLCANO[c]||[]).forEach(function(p){
      if(p.s==='N'||p.d==null||p.lp==null) return;
      ctx.beginPath(); ctx.arc(toX(p.d), toY(p.lp), 3.5, 0, 2*Math.PI);
      ctx.fillStyle = col+'cc'; ctx.fill();
      ctx.strokeStyle = col; ctx.lineWidth=0.8; ctx.stroke();
    });
  });

  // Titre
  ctx.fillStyle='#c9d1d9'; ctx.font='bold 11px sans-serif'; ctx.textAlign='center';
  ctx.fillText('Volcano comparatif — '+active.length+' contraste(s) actif(s)', lW+PW/2, 22);

  // Labels axes
  ctx.fillStyle='#8b949e'; ctx.font='9px sans-serif'; ctx.textAlign='center';
  ctx.fillText('log\\u2082 Fold Change', lW+PW/2, tH+PH+30);
  ctx.save(); ctx.translate(12, tH+PH/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('-log\\u2081\\u2080(p-value)', 0, 0); ctx.restore();

  // Légende compacte à droite
  var lgX = lW+PW+8, lgY = tH+8;
  active.forEach(function(c, i){
    var col = colorMap[c];
    var lbl = c.replace(/_vs_/,' vs ').replace(/_/g,' ');
    var words = lbl.split(' vs ');
    ctx.fillStyle=col; ctx.beginPath();
    ctx.arc(lgX+5, lgY+i*30+5, 5, 0, 2*Math.PI); ctx.fill();
    ctx.fillStyle='#c9d1d9'; ctx.font='bold 8px sans-serif'; ctx.textAlign='left';
    ctx.fillText(words[0]||lbl, lgX+13, lgY+i*30+4);
    if(words[1]){
      ctx.fillStyle='#8b949e'; ctx.font='8px sans-serif';
      ctx.fillText('vs '+words[1], lgX+13, lgY+i*30+14);
    }
  });

  // Tooltip sur les significatifs
  canvas._vcData = active.map(function(c){ return {c:c, col:colorMap[c], pts:VOLCANO[c]||[]}; });
  var tip = document.getElementById('volcCompTip');
  canvas.onmousemove = function(e){
    var rect = canvas.getBoundingClientRect();
    var mx = (e.clientX-rect.left)*(W/rect.width);
    var my = (e.clientY-rect.top)*(H/rect.height);
    var best=null, bestD=1e9;
    canvas._vcData.forEach(function(cd){
      cd.pts.forEach(function(p){
        if(p.s==='N'||p.d==null) return;
        var dx=toX(p.d)-mx, dy=toY(p.lp)-my;
        var dist=Math.sqrt(dx*dx+dy*dy);
        if(dist<bestD){bestD=dist;best={p:p,c:cd.c,col:cd.col};}
      });
    });
    if(tip&&best&&bestD<14){
      tip.style.display='block';
      tip.style.left=(e.clientX+12)+'px'; tip.style.top=(e.clientY-28)+'px';
      tip.innerHTML='<b>'+(best.p.Genes||best.p.name||'—')+'</b>'+
        '<br><span style="color:'+best.col+';">'+best.c.replace(/_vs_/,' vs ').replace(/_/g,' ')+'</span>'+
        '<br>LFC='+best.p.d+'&nbsp;&nbsp;-log10p='+best.p.lp;
    } else if(tip){ tip.style.display='none'; }
  };
  canvas.onmouseleave=function(){ if(tip) tip.style.display='none'; };
};

})();\
"""




GO_HEATMAP_JS = """\
/* ── GO HEATMAP CROSS-CONTRASTES ── */
function drawGOHeatmap(){
  if(!GO_HEATMAP||!GO_HEATMAP.terms||!GO_HEATMAP.terms.length)return;
  const canvas=document.getElementById('goHmCanvas');
  if(!canvas)return;

  // Filtre source
  const srcSel=document.getElementById('goHmSrc');
  const srcFilter=srcSel?srcSel.value:'all';
  let terms=GO_HEATMAP.terms.slice();
  if(srcFilter!=='all') terms=terms.filter(t=>t.src===srcFilter);
  if(!terms.length)return;

  const contrasts=GO_HEATMAP.contrasts;
  const nR=terms.length, nC=contrasts.length;
  const dpr=devicePixelRatio||1;
  const cellW=70, cellH=10, lW=260, topH=50, lgH=30;
  const W=lW+nC*cellW+20, H=topH+nR*cellH+lgH+20;
  canvas.width=W*dpr; canvas.height=H*dpr;
  canvas.style.width=W+'px'; canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.clearRect(0,0,W,H);

  const zAbs=Math.max(Math.abs(GO_HEATMAP.z_min),Math.abs(GO_HEATMAP.z_max),0.1);

  // Palette Blue->White->Red (z_score)
  function zColor(z){
    if(z==null) return '#1c2128';
    const t=Math.max(-1,Math.min(1,z/zAbs));
    if(t>=0){
      const f=t;
      return 'rgb('+Math.round(255)+','+Math.round(255*(1-f))+','+Math.round(255*(1-f))+')';
    } else {
      const f=-t;
      return 'rgb('+Math.round(255*(1-f))+','+Math.round(255*(1-f))+',255)';
    }
  }

  const SRC_COL={'GO:BP':'#388bfd','GO:MF':'#3fb950','GO:CC':'#d29922'};

  // En-têtes colonnes — labels courts pré-calculés côté Python (préfixe commun
  // retiré). Repli sur le nom complet si 'labels' absent (rétrocompat).
  const colLabels = (GO_HEATMAP.labels && GO_HEATMAP.labels.length===nC)
                    ? GO_HEATMAP.labels : contrasts;
  contrasts.forEach((c,ci)=>{
    const x=lW+ci*cellW+cellW/2;
    ctx.fillStyle='#c9d1d9';ctx.font='bold 9px sans-serif';ctx.textAlign='center';
    ctx.fillText(colLabels[ci], x, topH-6);
  });

  // Titre
  ctx.fillStyle='#c9d1d9';ctx.font='bold 11px sans-serif';ctx.textAlign='center';
  ctx.fillText('Heatmap GO cross-contrastes — z-score ('+terms.length+' termes)',
    lW+nC*cellW/2, 16);

  // Séparateurs de source
  let prevSrc=null;
  terms.forEach((t,ri)=>{
    const y=topH+ri*cellH;
    if(t.src!==prevSrc&&prevSrc!==null){
      ctx.strokeStyle='#30363d';ctx.lineWidth=1.5;
      ctx.beginPath();ctx.moveTo(0,y);ctx.lineTo(W-10,y);ctx.stroke();
    }
    prevSrc=t.src;

    // Barre source à gauche
    const col=SRC_COL[t.src]||'#888';
    ctx.fillStyle=col;
    ctx.fillRect(0,y+1,4,cellH-2);

    // Nom du terme
    ctx.fillStyle=t.n>1?'#c9d1d9':'#8b949e';
    ctx.font=(t.n>1?'bold ':'')+'9px sans-serif';
    ctx.textAlign='right';
    ctx.fillText(t.name.slice(0,36),lW-8,y+cellH-3);

    // Cellules
    t.row.forEach((z,ci)=>{
      ctx.fillStyle=zColor(z);
      ctx.fillRect(lW+ci*cellW,y,cellW-1,cellH-1);
      // Valeur si présente
      if(z!=null){
        ctx.fillStyle=Math.abs(z)>0.4?'#ffffff':'#333333';
        ctx.font='7px sans-serif';ctx.textAlign='center';
        ctx.fillText(z.toFixed(2),lW+ci*cellW+cellW/2,y+cellH-3);
      }
    });
  });

  // Légende z_score
  const lgX=lW, lgY=topH+nR*cellH+8;
  const lgW=nC*cellW;
  const grad=ctx.createLinearGradient(lgX,lgY,lgX+lgW,lgY);
  grad.addColorStop(0,'rgb(0,0,255)');
  grad.addColorStop(0.5,'rgb(255,255,255)');
  grad.addColorStop(1,'rgb(255,0,0)');
  ctx.fillStyle=grad;
  ctx.fillRect(lgX,lgY,lgW,10);
  ctx.strokeStyle='#30363d';ctx.lineWidth=0.5;
  ctx.strokeRect(lgX,lgY,lgW,10);
  [[0,(-zAbs).toFixed(2)],[0.5,'0'],[1,zAbs.toFixed(2)]].forEach(([f,v])=>{
    const x=lgX+f*lgW;
    ctx.fillStyle='#8b949e';ctx.font='8px sans-serif';ctx.textAlign='center';
    ctx.fillText(v,x,lgY+22);
  });
  ctx.fillStyle='#8b949e';ctx.font='8px sans-serif';ctx.textAlign='left';
  ctx.fillText('z-score',lgX,lgY+22);

  // Légende sources
  const slegX=lW+lgW+10;
  Object.entries(SRC_COL).forEach(([src,col],i)=>{
    ctx.fillStyle=col;ctx.fillRect(slegX,lgY+i*12,8,8);
    ctx.fillStyle='#8b949e';ctx.font='8px sans-serif';ctx.textAlign='left';
    ctx.fillText(src.replace('GO:',''),slegX+11,lgY+i*12+8);
  });

  // Tooltip
  canvas._goHmTerms=terms;canvas._goHmCtrasts=contrasts;
  canvas._goHmLW=lW;canvas._goHmTopH=topH;
  canvas._goHmCellW=cellW;canvas._goHmCellH=cellH;
  const tip=document.getElementById('goHmTip');
  canvas.onmousemove=function(e){
    const rect=canvas.getBoundingClientRect();
    const mx=(e.clientX-rect.left)*(W/rect.width);
    const my=(e.clientY-rect.top)*(H/rect.height);
    const ri=Math.floor((my-topH)/cellH);
    const ci=Math.floor((mx-lW)/cellW);
    if(ri>=0&&ri<terms.length&&ci>=0&&ci<contrasts.length){
      const t=terms[ri], z=t.row[ci];
      if(tip){
        tip.style.display='block';
        tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY-28)+'px';
        tip.innerHTML='<b>'+t.name+'</b><br>'+t.src+' · '+t.tid+
          '<br>'+contrasts[ci].replace(/_vs_/,' vs ').replace(/_/g,' ')+
          '<br>z-score: '+(z!=null?z:'n/a');
      }
    } else if(tip){tip.style.display='none';}
  };
  canvas.onmouseleave=function(){if(tip)tip.style.display='none';};
}\
"""



# ──────────────────────────────────────────────────────────────────────────────
# Patch JS — Réseau GO hiérarchique (DAG canvas natif)
# ──────────────────────────────────────────────────────────────────────────────

GO_NETWORK_JS = """\
/* ── RÉSEAU GO HIÉRARCHIQUE ── */
(function(){

window.drawGONetwork = function(){
  if(!GO_NETWORK||!Object.keys(GO_NETWORK).length)return;
  var sel=document.getElementById('goNetContrast');
  var c=sel?sel.value:Object.keys(GO_NETWORK)[0];
  var data=GO_NETWORK[c];
  if(!data||!data.nodes||!data.nodes.length)return;

  var canvas=document.getElementById('goNetCanvas');
  if(!canvas)return;
  var dpr=devicePixelRatio||1;

  var nodes=data.nodes, edges=data.edges;
  var SRC_COL={'GO:BP':'#388bfd','GO:MF':'#3fb950','GO:CC':'#d29922'};

  // Layout hiérarchique : regrouper par level
  var levels={};
  nodes.forEach(function(n){
    if(!levels[n.level]) levels[n.level]=[];
    levels[n.level].push(n);
  });
  var maxLevel=Math.max.apply(null,nodes.map(function(n){return n.level;}));
  var maxPerLevel=Math.max.apply(null,Object.values(levels).map(function(l){return l.length;}));

  var nodeR=function(n){return 5+Math.min(n.is,50)/50*12;};
  var colW=Math.max(120,900/Math.max(maxLevel+1,1));
  var rowH=Math.max(40,600/Math.max(maxPerLevel,1));
  var PAD=40;
  var W=PAD*2+(maxLevel+1)*colW;
  var H=PAD*2+maxPerLevel*rowH;

  canvas.width=W*dpr;canvas.height=H*dpr;
  canvas.style.width=W+'px';canvas.style.height=H+'px';
  var ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);ctx.clearRect(0,0,W,H);

  // Position des nœuds
  var pos={};
  Object.keys(levels).forEach(function(lvl){
    var arr=levels[lvl];
    arr.forEach(function(n,i){
      pos[n.id]={
        x: PAD+parseInt(lvl)*colW+colW/2,
        y: PAD+(i-(arr.length-1)/2)*rowH+H/2
      };
    });
  });

  // Arêtes
  ctx.strokeStyle='#30363d';ctx.lineWidth=1.5;
  edges.forEach(function(e){
    var a=pos[e.from],b=pos[e.to];
    if(!a||!b)return;
    ctx.beginPath();
    // Courbe de Bézier pour éviter le chevauchement
    var mx=(a.x+b.x)/2;
    ctx.moveTo(a.x,a.y);
    ctx.bezierCurveTo(mx,a.y,mx,b.y,b.x,b.y);
    ctx.stroke();
    // Flèche
    var angle=Math.atan2(b.y-a.y,b.x-a.x);
    var ar=nodeR({is:0})+3;
    ctx.beginPath();
    ctx.moveTo(b.x-ar*Math.cos(angle),b.y-ar*Math.sin(angle));
    ctx.lineTo(b.x-ar*Math.cos(angle)-6*Math.cos(angle-0.4),
               b.y-ar*Math.sin(angle)-6*Math.sin(angle-0.4));
    ctx.lineTo(b.x-ar*Math.cos(angle)-6*Math.cos(angle+0.4),
               b.y-ar*Math.sin(angle)-6*Math.sin(angle+0.4));
    ctx.closePath();ctx.fillStyle='#30363d';ctx.fill();
  });

  // Nœuds
  var maxZ=Math.max.apply(null,nodes.map(function(n){return Math.abs(n.z);})) || 0.5;
  nodes.forEach(function(n){
    var p=pos[n.id];if(!p)return;
    var r=nodeR(n);
    var col=SRC_COL[n.src]||'#8b949e';
    // Couleur de fond = z_score (rouge=up, bleu=down)
    var t=Math.max(-1,Math.min(1,n.z/maxZ));
    var fill;
    if(t>=0){fill='rgb('+Math.round(255)+','+Math.round(255*(1-t))+','+Math.round(255*(1-t))+')';}
    else{var f=-t;fill='rgb('+Math.round(255*(1-f))+','+Math.round(255*(1-f))+',255)';}
    // Cercle rempli
    ctx.beginPath();ctx.arc(p.x,p.y,r,0,2*Math.PI);
    ctx.fillStyle=fill+'cc';ctx.fill();
    ctx.strokeStyle=col;ctx.lineWidth=2;ctx.stroke();
    // Label court
    var lbl=n.name.length>18?n.name.slice(0,16)+'\u2026':n.name;
    ctx.fillStyle='#c9d1d9';ctx.font='9px sans-serif';ctx.textAlign='center';
    ctx.fillText(lbl,p.x,p.y+r+12);
    // p-value sous le label
    ctx.fillStyle='#8b949e';ctx.font='8px sans-serif';
    ctx.fillText('p='+n.p.toExponential(1),p.x,p.y+r+22);
  });

  // Titre
  ctx.fillStyle='#c9d1d9';ctx.font='bold 11px sans-serif';ctx.textAlign='center';
  ctx.fillText('Réseau GO hiérarchique — '+c.replace(/_vs_/,' vs ').replace(/_/g,' '),W/2,20);

  // Légende
  var lgX=10,lgY=H-70;
  [['GO:BP','#388bfd'],['GO:MF','#3fb950'],['GO:CC','#d29922']].forEach(function(kv,i){
    ctx.beginPath();ctx.arc(lgX+6,lgY+i*18,6,0,2*Math.PI);
    ctx.strokeStyle=kv[1];ctx.lineWidth=2;ctx.stroke();
    ctx.fillStyle='#8b949e';ctx.font='9px sans-serif';ctx.textAlign='left';
    ctx.fillText(kv[0],lgX+15,lgY+i*18+3);
  });
  ctx.fillStyle='#8b949e';ctx.font='8px sans-serif';ctx.textAlign='left';
  ctx.fillText('Taille \u221d intersection  |  Rouge=up  Bleu=down',lgX,lgY+54);

  // Tooltip
  canvas._goNetNodes=nodes;canvas._goNetPos=pos;canvas._goNetNodeR=nodeR;
  var tip=document.getElementById('goNetTip');
  canvas.onmousemove=function(e){
    var rect=canvas.getBoundingClientRect();
    var mx=(e.clientX-rect.left)*(W/rect.width);
    var my=(e.clientY-rect.top)*(H/rect.height);
    var best=null,bestD=1e9;
    nodes.forEach(function(n){
      var p=pos[n.id];if(!p)return;
      var dx=p.x-mx,dy=p.y-my;
      var d=Math.sqrt(dx*dx+dy*dy);
      if(d<bestD){bestD=d;best=n;}
    });
    if(tip&&best&&bestD<nodeR(best)+10){
      tip.style.display='block';
      tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY-28)+'px';
      tip.innerHTML='<b>'+best.name+'</b><br>'+best.src+' · '+best.id+
        '<br>p='+best.p.toExponential(2)+'  z='+best.z.toFixed(3)+
        '<br>n='+best.is+' gènes  niveau='+best.level;
    } else if(tip){tip.style.display='none';}
  };
  canvas.onmouseleave=function(){if(tip)tip.style.display='none';};
};

window.initGONetwork=function(){
  var sel=document.getElementById('goNetContrast');
  if(!sel||sel.options.length)return;
  Object.keys(GO_NETWORK).forEach(function(c){
    var o=document.createElement('option');o.value=c;
    o.textContent=c.replace(/_vs_/,' vs ').replace(/_/g,' ');sel.appendChild(o);
  });
  sel.onchange=drawGONetwork;
  drawGONetwork();
};

})();\
"""

GO_NETWORK_HTML = (
    '<div class="panel">\n'
    '  <div class="ph"><h2>R\u00e9seau GO hi\u00e9rarchique</h2>'
    '<div style="display:flex;gap:8px;align-items:center;">'
    '<select id="goNetContrast" style="background:var(--bg2);border:1px solid var(--bd);'
    'color:var(--tx);border-radius:4px;padding:3px 8px;font-size:11px;"></select>'
    '<button id="goNetExportBtn" class="exp-btn" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button>'
    '</div></div>\n'
    '  <div style="overflow:auto;">'
    '<canvas id="goNetCanvas" style="display:block;"></canvas></div>\n'
    '  <div id="goNetTip" style="position:fixed;display:none;background:var(--bg2);'
    'border:1px solid var(--bd);border-radius:6px;padding:6px 10px;font-size:11px;'
    'color:var(--tx);pointer-events:none;z-index:999;max-width:300px;"></div>\n'
    '</div>\n'
)


def patch_go_network_js(html):
    """Injecte drawGONetwork, initGONetwork + panel dans #gopage."""
    import re as _re

    # JS avant </script>
    pos = html.rfind('</script>')
    if pos == -1:
        print("  [warn] patch_go_network_js : </script> introuvable.")
        return html
    html = html[:pos] + "\n\n" + GO_NETWORK_JS + "\n" + html[pos:]

    # HTML : insérer avant le panel Bubble Plot
    target = '<div class="panel">\n  <div class="ph">\n    <h2>GO Enrichment'
    if target in html:
        html = html.replace(target, GO_NETWORK_HTML + target, 1)
        print("  patch_go_network_js : réseau GO injecté dans #gopage.")
    else:
        # Fallback : insérer avant <!-- HEATMAP -->
        import re as _re2
        m = _re2.search(r'(</div>)\s*\n\s*<!-- HEATMAP -->', html)
        if m:
            html = html[:m.start(1)] + '\n' + GO_NETWORK_HTML + html[m.start(1):]
            print("  patch_go_network_js : réseau GO injecté (fallback).")
        else:
            print("  [warn] patch_go_network_js : point d'insertion introuvable.")

    # Brancher initGONetwork + export dans showPage gopage via re.sub sur le monkey-patch
    html = _re.sub(
        r"(if\(typeof initGOBubble === 'function'\) initGOBubble\(\);)",
        r"\1\n      if(typeof initGONetwork === 'function') initGONetwork();"
        r"\n      var nb=document.getElementById('goNetExportBtn');"
        r"\n      if(nb) nb.onclick=function(){exportCanvas('goNetCanvas','proteogen_GO_network');};",
        html, count=1
    )
    return html



def patch_go_heatmap_js(html):
    """Injecte drawGOHeatmap() avant </script> et ajoute le panel dans #gopage."""
    import re as _re

    # JS
    pos = html.rfind('</script>')
    if pos == -1:
        print("  [warn] patch_go_heatmap_js : </script> introuvable.")
        return html
    html = html[:pos] + "\n\n" + GO_HEATMAP_JS + "\n" + html[pos:]

    # HTML panel dans #gopage — insérer entre bubble plot et heatmap terme
    # Le panel GO bubble est déjà là, on ajoute le panel heatmap juste après
    HM_HTML = (
        '<div class="panel">\n'
        '  <div class="ph"><h2>Heatmap GO cross-contrastes</h2>'
        '<div style="display:flex;gap:8px;align-items:center;">'
        '<select id="goHmSrc" style="background:var(--bg2);border:1px solid var(--bd);'
        'color:var(--tx);border-radius:4px;padding:3px 8px;font-size:11px;">'
        '<option value="all">Toutes sources</option>'
        '<option value="GO:BP">BP</option>'
        '<option value="GO:MF">MF</option>'
        '<option value="GO:CC">CC</option>'
        '</select>'
        '<button id="goHmExportBtn" class="exp-btn" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button>'
        '</div></div>\n'
        '  <div style="overflow-y:auto;max-height:580px;">'
        '<canvas id="goHmCanvas" style="display:block;"></canvas>'
        '</div>\n'
        '  <div id="goHmTip" style="position:fixed;display:none;background:var(--bg2);'
        'border:1px solid var(--bd);border-radius:6px;padding:6px 10px;font-size:11px;'
        'color:var(--tx);pointer-events:none;z-index:999;max-width:300px;"></div>\n'
        '</div>\n'
    )

    # Insérer la heatmap GO juste après l'ouverture de #gopage
    # (avant goSourcePanels et le bubble plot)
    m = _re.search(r'(id=["\x27]gopage["\x27][^>]*>)', html)
    if m:
        insert_pos = m.end()
        html = html[:insert_pos] + '\n' + HM_HTML + html[insert_pos:]
        print("  patch_go_heatmap_js : heatmap GO cross-contrastes injectée.")

    # Les boutons export/filtre goHmSrc sont branchés via initGOBubble dans QC_INIT_SCRIPT
    # Plus besoin de monkey-patch ici
    return html





# ──────────────────────────────────────────────────────────────────────────────
# Patch JS — Corrélation module-trait WGCNA
# ──────────────────────────────────────────────────────────────────────────────

MODULE_TRAIT_JS = """\
/* ── MODULE-TRAIT WGCNA ── */
function drawModuleTrait(){
  if(!MODULE_TRAIT||!Object.keys(MODULE_TRAIT).length)return;
  const canvas=document.getElementById('modulTraitCanvas');
  if(!canvas)return;

  const modules=Object.keys(MODULE_TRAIT);
  const conditions=Object.keys(MODULE_TRAIT[modules[0]].corr);
  const nR=modules.length, nC=conditions.length;
  const dpr=devicePixelRatio||1;
  const cellW=110, cellH=50, lW=90, topH=80, botH=20;
  const W=lW+nC*cellW+20, H=topH+nR*cellH+botH;
  canvas.width=W*dpr; canvas.height=H*dpr;
  canvas.style.width=W+'px'; canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.clearRect(0,0,W,H);

  const MOD_COLORS={
    turquoise:'#1abc9c',blue:'#3498db',brown:'#a0522d',
    yellow:'#f1c40f',green:'#2ecc71',red:'#e74c3c',black:'#95a5a6',
    pink:'#e91e63',purple:'#9b59b6'
  };

  // Palette corrélation : Blue(-1) → White(0) → Red(+1)
  function corrColor(r){
    const t=Math.max(-1,Math.min(1,r));
    if(t>=0){
      const f=t;
      return 'rgb('+Math.round(255)+','+Math.round(255*(1-f))+','+Math.round(255*(1-f))+')';
    } else {
      const f=-t;
      return 'rgb('+Math.round(255*(1-f))+','+Math.round(255*(1-f))+',255)';
    }
  }

  // Titre
  ctx.fillStyle='#c9d1d9';ctx.font='bold 11px sans-serif';ctx.textAlign='center';
  ctx.fillText('Corrélation module–trait (eigengène × condition)',lW+nC*cellW/2,18);

  // En-têtes colonnes (conditions)
  conditions.forEach((cond,ci)=>{
    const x=lW+ci*cellW+cellW/2;
    const col=COND_COLORS[cond]||'#8b949e';
    ctx.fillStyle=col;ctx.font='bold 10px sans-serif';ctx.textAlign='center';
    const lbl=cond.replace(/Polyculture[._]?/i,'Poly').replace(/Monoculture/i,'Mono');
    // Wrap sur 2 lignes si nécessaire
    ctx.fillText(lbl,x,topH-28);
  });

  // Lignes = modules
  modules.forEach((mod,ri)=>{
    const y=topH+ri*cellH;
    const d=MODULE_TRAIT[mod];
    const modCol=MOD_COLORS[mod]||'#888';

    // Barre couleur module à gauche
    ctx.fillStyle=modCol;
    ctx.fillRect(0,y+2,6,cellH-4);

    // Label module
    ctx.fillStyle='#c9d1d9';ctx.font='bold 9px sans-serif';ctx.textAlign='right';
    ctx.fillText(mod,lW-10,y+cellH/2+4);

    // Cellules
    conditions.forEach((cond,ci)=>{
      const r=d.corr[cond]||0;
      const p=d.pval[cond]||1;
      const x=lW+ci*cellW;

      // Fond coloré
      ctx.fillStyle=corrColor(r);
      ctx.fillRect(x+2,y+2,cellW-4,cellH-4);

      // Valeur r
      const textColor=Math.abs(r)>0.5?'#ffffff':'#333333';
      ctx.fillStyle=textColor;ctx.font='bold 11px sans-serif';ctx.textAlign='center';
      ctx.fillText(r.toFixed(3),x+cellW/2,y+cellH/2);

      // p-value (plus petit, en dessous)
      ctx.font='8px sans-serif';
      const pStr=p<0.001?'p<0.001':'p='+p.toFixed(3);
      const sig=p<0.05?'*':'';
      ctx.fillText(pStr+sig,x+cellW/2,y+cellH/2+13);

      // Bordure
      ctx.strokeStyle='#21262d';ctx.lineWidth=0.5;
      ctx.strokeRect(x+2,y+2,cellW-4,cellH-4);
    });

    // N protéines
    ctx.fillStyle='#6e7681';ctx.font='8px sans-serif';ctx.textAlign='left';
    ctx.fillText('n='+d.n,4,y+cellH/2+14);
  });

  // Légende
  const lgX=lW+nC*cellW+4, lgY=topH+10;
  const lgW=16, lgH2=nR*cellH-20;
  const grad=ctx.createLinearGradient(lgX,lgY,lgX,lgY+lgH2);
  grad.addColorStop(0,'rgb(255,0,0)');
  grad.addColorStop(0.5,'rgb(255,255,255)');
  grad.addColorStop(1,'rgb(0,0,255)');
  ctx.fillStyle=grad;
  ctx.fillRect(lgX,lgY,lgW,lgH2);
  ctx.strokeStyle='#30363d';ctx.lineWidth=0.5;
  ctx.strokeRect(lgX,lgY,lgW,lgH2);
  [[0,'+1'],[0.5,'0'],[1,'-1']].forEach(([f,v])=>{
    const y2=lgY+f*lgH2;
    ctx.fillStyle='#8b949e';ctx.font='8px sans-serif';ctx.textAlign='left';
    ctx.fillText(v,lgX+lgW+3,y2+3);
  });
  ctx.fillStyle='#8b949e';ctx.font='8px sans-serif';ctx.textAlign='center';
  ctx.save();ctx.translate(lgX+lgW/2,lgY+lgH2/2);ctx.rotate(-Math.PI/2);
  ctx.fillText('r Pearson',0,0);ctx.restore();

  // Note * p<0.05
  ctx.fillStyle='#6e7681';ctx.font='9px sans-serif';ctx.textAlign='left';
  ctx.fillText('* p < 0.05',lW,topH+nR*cellH+14);

  // Tooltip
  canvas._mtModules=modules;canvas._mtConditions=conditions;
  canvas._mtLW=lW;canvas._mtTopH=topH;
  canvas._mtCellW=cellW;canvas._mtCellH=cellH;
  const tip=document.getElementById('modulTraitTip');
  canvas.onmousemove=function(e){
    const rect=canvas.getBoundingClientRect();
    const mx=(e.clientX-rect.left)*(W/rect.width);
    const my=(e.clientY-rect.top)*(H/rect.height);
    const ri=Math.floor((my-topH)/cellH);
    const ci=Math.floor((mx-lW)/cellW);
    if(ri>=0&&ri<modules.length&&ci>=0&&ci<conditions.length){
      const mod=modules[ri],cond=conditions[ci];
      const d=MODULE_TRAIT[mod];
      if(tip){
        tip.style.display='block';
        tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY-28)+'px';
        tip.innerHTML='<b>'+mod+'</b> × <b>'+cond+'</b><br>'+
          'r = '+d.corr[cond].toFixed(4)+'<br>'+
          'p = '+d.pval[cond]+'<br>'+
          'n = '+d.n+' protéines';
      }
    } else if(tip){tip.style.display='none';}
  };
  canvas.onmouseleave=function(){if(tip)tip.style.display='none';};
}\
"""

MODULE_TRAIT_HTML = (
    '<div class="panel">\n'
    '  <div class="ph"><h2>Corr\u00e9lation module\u2013trait</h2>'
    '<button id="mtExportBtn" class="exp-btn" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button>'
    '</div>\n'
    '  <canvas id="modulTraitCanvas" style="display:block;width:100%;"></canvas>\n'
    '  <div id="modulTraitTip" style="position:fixed;display:none;background:var(--bg2);'
    'border:1px solid var(--bd);border-radius:6px;padding:6px 10px;font-size:11px;'
    'color:var(--tx);pointer-events:none;z-index:999;"></div>\n'
    '</div>\n'
)


def patch_module_trait_js(html):
    """Injecte drawModuleTrait() et son panel dans #wgcna."""
    import re as _re

    # JS
    pos = html.rfind('</script>')
    if pos == -1:
        print("  [warn] patch_module_trait_js : </script> introuvable.")
        return html
    html = html[:pos] + "\n\n" + MODULE_TRAIT_JS + "\n" + html[pos:]

    # HTML : insérer avant le network WGCNA dans #wgcna
    idx_wg = html.find('id="wgcna"')
    if idx_wg < 0:
        print("  [warn] patch_module_trait_js : #wgcna introuvable.")
        return html
    # Trouver la fin de #wgcna (avant <script>)
    chunk = html[idx_wg:]
    m = _re.search(r'(</div>\s*\n)\s*\n\s*<script', chunk)
    if m:
        insert_pos = idx_wg + m.start(1)
        html = html[:insert_pos] + '\n' + MODULE_TRAIT_HTML + html[insert_pos:]
        print("  patch_module_trait_js : corrélation module-trait injectée dans #wgcna.")
    else:
        print("  [warn] patch_module_trait_js : fermeture #wgcna introuvable.")

    # Brancher drawModuleTrait dans monkey-patch showPage('wgcna')
    html = _re.sub(
        r"(if\(id === 'wgcna' && typeof initWGCNANetwork === 'function'\) initWGCNANetwork\(\);)",
        r"\1\n    if(id === 'wgcna' && typeof drawModuleTrait === 'function') { "
        r"drawModuleTrait(); "
        r"var b=document.getElementById('mtExportBtn'); "
        r"if(b)b.onclick=function(){exportCanvas('modulTraitCanvas','proteogen_WGCNA_module_trait');}; }",
        html, count=1
    )
    return html





# ──────────────────────────────────────────────────────────────────────────────
# Patch JS — Hubs WGCNA sur les volcanos (individuel + comparatif)
# ──────────────────────────────────────────────────────────────────────────────

WGCNA_HUB_VOLCANO_JS = """\
/* ── WGCNA HUBS SUR VOLCANO ── */
(function(){

// Couleurs modules WGCNA (identiques à drawModuleTrait)
var MOD_COLORS = {
  turquoise:'#1abc9c', blue:'#3498db', brown:'#a0522d',
  yellow:'#f1c40f', green:'#2ecc71', red:'#e74c3c',
  black:'#95a5a6', pink:'#e91e63', purple:'#9b59b6'
};

// Monkey-patch drawVolcano pour afficher les anneaux hubs
var _origDrawVolcano = window.drawVolcano;
window.drawVolcano = function(){
  if(_origDrawVolcano) _origDrawVolcano.apply(this, arguments);
  if(!WGCNA_HUB_MAP || !Object.keys(WGCNA_HUB_MAP).length) return;

  var canvas = document.getElementById('volcCanvas');
  if(!canvas) return;
  var ctx = canvas.getContext('2d');
  var dpr = devicePixelRatio || 1;
  var W = canvas.width/dpr, H = canvas.height/dpr;

  // Récupérer le contraste actif
  var sel = document.getElementById('volcContrast');
  var c = sel ? sel.value : (CONTRASTS && CONTRASTS[0]);
  if(!c || !VOLCANO[c]) return;
  var pts = VOLCANO[c];

  // Récupérer les axes depuis le canvas (re-calculer toX/toY)
  var lW=60, rW=20, tH=36, bH=44, PW=W-lW-rW, PH=H-tH-bH;
  var allD = pts.map(function(p){return p.d;}).filter(function(v){return v!=null;});
  var allLP = pts.map(function(p){return p.lp;}).filter(function(v){return v!=null;});
  if(!allD.length) return;
  var dMin=Math.min.apply(null,allD), dMax=Math.max.apply(null,allD);
  var pad=(dMax-dMin)*0.05; dMin-=pad; dMax+=pad;
  var lpMax=Math.max.apply(null,allLP)*1.1||4;
  function toX(d){return lW+(d-dMin)/(dMax-dMin)*PW;}
  function toY(lp){return tH+PH*(1-lp/lpMax);}

  // Dessiner les anneaux hubs
  pts.forEach(function(p){
    if(!p.name && !p.Genes) return;
    var hub = WGCNA_HUB_MAP[p.name];
    if(!hub) return;
    if(p.d==null || p.lp==null) return;
    var col = MOD_COLORS[hub.mod] || '#fff';
    var x = toX(p.d), y = toY(p.lp);
    var r = 2.5 + hub.hub * 3; // rayon proportionnel au hub score
    ctx.beginPath(); ctx.arc(x, y, r, 0, 2*Math.PI);
    ctx.strokeStyle = col;
    ctx.lineWidth = 1.5;
    ctx.stroke();
  });

  // Légende hubs
  var modsPresent = {};
  pts.forEach(function(p){
    var hub = WGCNA_HUB_MAP[p.name];
    if(hub) modsPresent[hub.mod] = true;
  });
  var mods = Object.keys(modsPresent);
  if(!mods.length) return;
  var lgX = lW+PW+4, lgY = tH+10;
  ctx.fillStyle='#8b949e'; ctx.font='bold 8px sans-serif'; ctx.textAlign='left';
  ctx.fillText('Hubs WGCNA', lgX, lgY);
  mods.forEach(function(mod, i){
    var col = MOD_COLORS[mod] || '#fff';
    var y2 = lgY + 14 + i*14;
    ctx.beginPath(); ctx.arc(lgX+5, y2, 5, 0, 2*Math.PI);
    ctx.strokeStyle=col; ctx.lineWidth=1.5; ctx.stroke();
    ctx.fillStyle='#c9d1d9'; ctx.font='8px sans-serif'; ctx.textAlign='left';
    ctx.fillText(mod, lgX+13, y2+3);
  });
};

// Monkey-patch drawVolcanoComp pour afficher les anneaux hubs
var _origDrawVolcanoComp = window.drawVolcanoComp;
window.drawVolcanoComp = function(){
  if(_origDrawVolcanoComp) _origDrawVolcanoComp.apply(this, arguments);
  if(!WGCNA_HUB_MAP || !Object.keys(WGCNA_HUB_MAP).length) return;

  var canvas = document.getElementById('volcCompCanvas');
  if(!canvas) return;
  var ctx = canvas.getContext('2d');
  var dpr = devicePixelRatio||1;
  var W=canvas.width/dpr, H=canvas.height/dpr;

  var checks = document.querySelectorAll('#volcCompChecks input:checked');
  var active=[];
  checks.forEach(function(cb){active.push(cb.value);});
  if(!active.length) return;

  // Recalculer les axes (copie de drawVolcanoComp)
  var lW=60,rW=110,tH=36,bH=44,PW=W-lW-rW,PH=H-tH-bH;
  var allD=[],allLP=[];
  active.forEach(function(c){
    (VOLCANO[c]||[]).forEach(function(p){
      if(p.d!=null) allD.push(p.d);
      if(p.lp!=null) allLP.push(p.lp);
    });
  });
  if(!allD.length) return;
  var dMin=Math.min.apply(null,allD)*1.05, dMax=Math.max.apply(null,allD)*1.05;
  var lpMax=Math.max.apply(null,allLP)*1.1||4;
  function toX(d){return lW+(d-dMin)/(dMax-dMin)*PW;}
  function toY(lp){return tH+PH*(1-lp/lpMax);}

  // Anneaux hubs sur tous les contrastes actifs
  var modsPresent={};
  active.forEach(function(c){
    (VOLCANO[c]||[]).forEach(function(p){
      var hub=WGCNA_HUB_MAP[p.name];
      if(!hub||p.d==null||p.lp==null) return;
      modsPresent[hub.mod]=true;
      var col=MOD_COLORS[hub.mod]||'#fff';
      var x=toX(p.d),y=toY(p.lp);
      var r=2.5+hub.hub*3;
      ctx.beginPath();ctx.arc(x,y,r,0,2*Math.PI);
      ctx.strokeStyle=col;ctx.lineWidth=1.5;ctx.stroke();
    });
  });
};

})();\
"""


def patch_wgcna_hub_volcano_js(html):
    """Injecte les monkey-patches volcano + volcanoComp pour afficher les hubs WGCNA."""
    pos = html.rfind('</script>')
    if pos == -1:
        print("  [warn] patch_wgcna_hub_volcano_js : </script> introuvable.")
        return html
    html = html[:pos] + "\n\n" + WGCNA_HUB_VOLCANO_JS + "\n" + html[pos:]
    print("  patch_wgcna_hub_volcano_js : hubs WGCNA sur volcanos injectés.")
    return html



def patch_volcano_comp_js(html):
    """Injecte le volcano comparatif et branche initVolcanoComp dans showPage."""
    pos = html.rfind('</script>')
    if pos == -1:
        print("  [warn] patch_volcano_comp_js : </script> introuvable.")
        return html
    html = html[:pos] + "\n\n" + DRAW_VOLCANO_COMP_JS + "\n" + html[pos:]

    # Brancher dans monkey-patch showPage — chercher la ligne gopage
    import re as _re
    html = _re.sub(
        r"(if\(id === 'gopage' && typeof initGOBubble === 'function'\) initGOBubble\(\);)",
        r"\1\n    if(id === 'volcano' && typeof initVolcanoComp === 'function') initVolcanoComp();",
        html, count=1
    )
    print("  patch_volcano_comp_js : volcano comparatif injecté.")
    return html


# ──────────────────────────────────────────────────────────────────────────────
# Patch JS — Ranked abundance + CV plot + panneau À propos
# ──────────────────────────────────────────────────────────────────────────────

DRAW_RANKED_CV_JS = """\
/* ── RANKED ABUNDANCE ── */
function drawRankedAbundance(){
  const canvas=document.getElementById('rankedCanvas');
  if(!canvas||!RANKED_DATA||!RANKED_DATA.points||!RANKED_DATA.points.length)return;
  const pts=RANKED_DATA.points;
  const N=RANKED_DATA.n_total||pts.length;  // N total réel pour les axes
  const dpr=devicePixelRatio||1;
  const lW=60,rW=20,tH=36,bH=44,PW=600,PH=300;
  const W=lW+PW+rW,H=tH+PH+bH;
  canvas.width=W*dpr;canvas.height=H*dpr;
  canvas.style.width=W+'px';canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.clearRect(0,0,W,H);

  const vMin=RANKED_DATA.vmin,vMax=RANKED_DATA.vmax;
  function toX(r){return lW+(r-1)/(N-1)*PW;}
  function toY(v){return tH+PH*(1-(v-vMin)/(vMax-vMin));}

  // Grille
  ctx.strokeStyle='#21262d';ctx.lineWidth=1;
  for(let v=Math.ceil(vMin);v<=vMax;v+=2){
    const y=toY(v);
    ctx.beginPath();ctx.moveTo(lW,y);ctx.lineTo(lW+PW,y);ctx.stroke();
    ctx.fillStyle='#8b949e';ctx.font='8px sans-serif';ctx.textAlign='right';
    ctx.fillText(v,lW-4,y+3);
  }
  [0.1,0.25,0.5,0.75,0.9].forEach(f=>{
    const x=lW+f*PW;
    ctx.beginPath();ctx.moveTo(x,tH);ctx.lineTo(x,tH+PH);ctx.stroke();
    ctx.fillStyle='#8b949e';ctx.font='8px sans-serif';ctx.textAlign='center';
    ctx.fillText(Math.round(f*N),x,tH+PH+12);
  });

  // Titre et axes
  ctx.fillStyle='#c9d1d9';ctx.font='bold 11px sans-serif';ctx.textAlign='center';
  ctx.fillText('Ranked protein abundance ('+N+' protéines)',lW+PW/2,18);
  ctx.fillStyle='#8b949e';ctx.font='9px sans-serif';
  ctx.fillText('Rang (intensité décroissante)',lW+PW/2,tH+PH+28);
  ctx.fillStyle='#8b949e';ctx.save();ctx.translate(12,tH+PH/2);ctx.rotate(-Math.PI/2);
  ctx.fillText('log\\u2082 LFQ (moyenne)',0,0);ctx.restore();

  // Courbe dégradée — pts sous-échantillonnés avec rang réel dans p.r
  for(let i=0;i<pts.length-1;i++){
    const p=pts[i],pn=pts[i+1];
    const f=1-p.r/(N-1);
    const r=Math.round(f*220+35*(1-f));
    const b=Math.round((1-f)*220+35*f);
    ctx.strokeStyle='rgb('+r+',80,'+b+')';
    ctx.lineWidth=1.5;
    ctx.beginPath();ctx.moveTo(toX(p.r),toY(p.v));
    ctx.lineTo(toX(pn.r),toY(pn.v));ctx.stroke();
  }

  // Annotations top 5 et bottom 5
  const annot=[...pts.slice(0,5),...pts.slice(-5)];
  annot.forEach((p,i)=>{
    const x=toX(p.r),y=toY(p.v);
    const isTop=i<5;
    ctx.fillStyle=isTop?'#c0392b':'#2980b9';
    ctx.beginPath();ctx.arc(x,y,4,0,2*Math.PI);ctx.fill();
    if(p.g){
      ctx.fillStyle='#333';ctx.font='bold 8px sans-serif';
      ctx.textAlign=isTop?'left':'right';
      ctx.fillText(p.g.slice(0,12),x+(isTop?6:-6),y+(isTop?-6:10));
    }
  });

  // Tooltip
  const tip=document.getElementById('rankedTip');
  canvas.onmousemove=function(e){
    const rect=canvas.getBoundingClientRect();
    const mx=(e.clientX-rect.left)*(W/rect.width);
    const rank=Math.round((mx-lW)/PW*(N-1))+1;
    const idx=Math.max(0,Math.min(rank-1,pts.length-1));
    const p=pts[idx];
    if(tip&&mx>=lW&&mx<=lW+PW){
      tip.style.display='block';
      tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY-28)+'px';
      tip.innerHTML='<b>'+(p.g||p.n||'—')+'</b><br>Rang '+p.r+
        '<br>'+p.v.toFixed(3)+' log\\u2082 LFQ';
    }
  };
  canvas.onmouseleave=function(){if(tip)tip.style.display='none';};
}

/* ── CV PLOT ── */
function drawCVPlot(){
  const canvas=document.getElementById('cvCanvas');
  if(!canvas||!CV_DATA||!Object.keys(CV_DATA).length)return;
  const sel=document.getElementById('cvCond');
  const cond=sel?sel.value:Object.keys(CV_DATA)[0];
  const d=CV_DATA[cond];
  if(!d||!d.pts)return;

  const dpr=devicePixelRatio||1;
  const lW=55,rW=20,tH=36,bH=44,PW=500,PH=300;
  const W=lW+PW+rW,H=tH+PH+bH;
  canvas.width=W*dpr;canvas.height=H*dpr;
  canvas.style.width=W+'px';canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.clearRect(0,0,W,H);

  const xMin=Math.min(...d.pts.map(p=>p.x));
  const xMax=Math.max(...d.pts.map(p=>p.x));
  const yMax=Math.min(Math.max(...d.pts.map(p=>p.y)),100);
  function toX(v){return lW+(v-xMin)/(xMax-xMin)*PW;}
  function toY(v){return tH+PH*(1-Math.min(v,yMax)/yMax);}

  // Grille
  ctx.strokeStyle='#21262d';ctx.lineWidth=1;
  [0,20,40,60,80,100].forEach(cv=>{
    if(cv>yMax)return;
    const y=toY(cv);
    ctx.beginPath();ctx.moveTo(lW,y);ctx.lineTo(lW+PW,y);ctx.stroke();
    ctx.fillStyle='#8b949e';ctx.font='8px sans-serif';ctx.textAlign='right';
    ctx.fillText(cv+'%',lW-4,y+3);
  });

  // Ligne médiane CV
  const medY=toY(d.median);
  ctx.strokeStyle='#e67e22';ctx.lineWidth=1.5;ctx.setLineDash([5,3]);
  ctx.beginPath();ctx.moveTo(lW,medY);ctx.lineTo(lW+PW,medY);ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle='#e67e22';ctx.font='bold 9px sans-serif';ctx.textAlign='left';
  ctx.fillText('médiane '+d.median+'%',lW+4,medY-4);

  // Points
  const col=COND_COLORS[cond]||'#4a90d9';
  d.pts.forEach(p=>{
    ctx.beginPath();ctx.arc(toX(p.x),toY(p.y),2.5,0,2*Math.PI);
    ctx.fillStyle=col+'99';ctx.fill();
  });

  // Titre + axes
  ctx.fillStyle='#c9d1d9';ctx.font='bold 11px sans-serif';ctx.textAlign='center';
  ctx.fillText('CV intra-condition — '+cond+' ('+d.n+' prot.)',lW+PW/2,18);
  ctx.fillStyle='#8b949e';ctx.font='9px sans-serif';
  ctx.fillText('Intensité log\\u2082 LFQ moyenne',lW+PW/2,tH+PH+28);
  ctx.fillStyle='#8b949e';ctx.save();ctx.translate(12,tH+PH/2);ctx.rotate(-Math.PI/2);
  ctx.fillText('CV (%)',0,0);ctx.restore();

  // Tooltip
  const tip=document.getElementById('cvTip');
  canvas.onmousemove=function(e){
    const rect=canvas.getBoundingClientRect();
    const mx=(e.clientX-rect.left)*(W/rect.width);
    const my=(e.clientY-rect.top)*(H/rect.height);
    let best=null,bestD=1e9;
    d.pts.forEach(p=>{
      const dx=toX(p.x)-mx,dy=toY(p.y)-my;
      const dist=Math.sqrt(dx*dx+dy*dy);
      if(dist<bestD){bestD=dist;best=p;}
    });
    if(tip&&best&&bestD<14){
      tip.style.display='block';
      tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY-28)+'px';
      tip.innerHTML='<b>'+(best.g||'—')+'</b><br>Intensité: '+best.x+
        '<br>CV: '+best.y.toFixed(1)+'%';
    } else if(tip){tip.style.display='none';}
  };
  canvas.onmouseleave=function(){if(tip)tip.style.display='none';};
}
function initCVPlot(){
  const sel=document.getElementById('cvCond');
  if(!sel||sel.options.length)return;
  Object.keys(CV_DATA).forEach(c=>{
    const o=document.createElement('option');o.value=c;o.textContent=c;sel.appendChild(o);
  });
  sel.onchange=drawCVPlot;
  drawCVPlot();
}

/* ── PANNEAU À PROPOS ── */
function renderAbout(){
  const el=document.getElementById('aboutPanel');
  if(!el||!ABOUT_DATA)return;
  const d=ABOUT_DATA;
  let html='<div style="display:grid;grid-template-columns:1fr 1fr;gap:14px;">';
  // Génération
  html+='<div class="panel"><div class="ph"><h2>Génération</h2></div>';
  html+='<p style="font-size:11px;color:var(--tx2);">Date : <b style="color:var(--tx);">'+d.generated_at+'</b></p>';
  // Fichiers sources (dédupliqués : stats/go/umap peuvent pointer le même .xlsx)
  html+='<table style="width:100%;border-collapse:collapse;font-size:10px;margin-top:8px;">';
  html+='<tr style="color:var(--tx2);border-bottom:1px solid var(--bd);">';
  html+='<th style="text-align:left;padding:3px 6px;">Fichier</th>';
  html+='<th style="text-align:right;padding:3px 6px;">Taille</th>';
  html+='<th style="text-align:right;padding:3px 6px;">MD5</th></tr>';
  var _seen={};
  Object.entries(d.files).forEach(([k,f])=>{
    if(!f)return;
    var key=f.name+'|'+f.md5;
    if(_seen[key])return;        // évite le doublon (même fichier)
    _seen[key]=true;
    html+='<tr style="border-bottom:1px solid var(--bd);">';
    html+='<td style="padding:3px 6px;color:var(--tx);">'+f.name+'</td>';
    html+='<td style="padding:3px 6px;text-align:right;color:var(--tx2);">'+f.size_kb+' KB</td>';
    html+='<td style="padding:3px 6px;text-align:right;font-family:monospace;color:var(--bl);">'+f.md5+'</td></tr>';
  });
  html+='</table></div>';
  // Paramètres
  html+='<div class="panel"><div class="ph"><h2>Param\\u00e8tres d\\u0027analyse</h2></div>';
  const p=d.params;
  const rows=[
    ['Prot\u00e9ines totales', p.n_proteins],
    ['\u00c9chantillons', p.n_samples],
    ['Conditions', p.n_conditions],
    ['Contrastes', p.n_contrasts],
  ];
  // Seuils choisis (volcanos + ANOVA/heatmap), si disponibles
  if(p.volcano_p!==undefined){
    rows.push(['Seuil Volcanos', p.volcano_p_type+' < '+p.volcano_p+'  |  ratio \u2265 '+p.volcano_ratio]);
  }
  if(p.anova_p!==undefined){
    rows.push(['Seuil ANOVA/Heatmap', p.anova_p_type+' < '+p.anova_p]);
  }
  html+='<table style="width:100%;border-collapse:collapse;font-size:11px;">';
  rows.forEach(([k,v])=>{
    html+='<tr style="border-bottom:1px solid var(--bd);">';
    html+='<td style="padding:4px 6px;color:var(--tx2);">'+k+'</td>';
    html+='<td style="padding:4px 6px;font-weight:700;color:var(--tx);">'+v+'</td></tr>';
  });
  html+='</table>';
  html+='<div style="margin-top:10px;font-size:10px;color:var(--tx2);">Contrastes : ';
  html+=d.contrasts.map(c=>'<span style="display:inline-block;background:var(--bg);border:1px solid var(--bd);border-radius:4px;padding:1px 6px;margin:2px;font-size:9px;">'+c.replace(/_vs_/,' vs ').replace(/_/g,' ')+'</span>').join('');
  html+='</div></div></div>';
  el.innerHTML=html;
}\
"""


INIT_GUARD_JS = """\
/* ── GUARD INIT + réimplémentation renderGO en canvas natif ── */
(function(){

  window.renderGO = function(contrast){
    var c=contrast||(typeof currentGOContrast!=='undefined'?currentGOContrast:null)||(CONTRASTS&&CONTRASTS[0]);
    // Mettre à jour currentGOContrast pour éviter que l'ancien initGO réinitialise
    if(typeof currentGOContrast !== 'undefined') currentGOContrast = c;
    if(!c||!GO_DATA||!GO_DATA[c])return;
    var terms=(typeof goInvertActive!=='undefined'&&goInvertActive)
      ? GO_DATA[c].map(function(t){return Object.assign({},t,{z_score:-(t.z_score||0)});})
      : GO_DATA[c];
    var srcColors={'GO:BP':'#388bfd','GO:MF':'#3fb950','GO:CC':'#d29922'};

    // Panels par source
    var panelsDiv=document.getElementById('goSourcePanels');
    if(panelsDiv){
      panelsDiv.innerHTML='';
      ['GO:BP','GO:MF','GO:CC'].forEach(function(src){
        var st=terms.filter(function(t){return t.source===src;})
                    .sort(function(a,b){return a.p_value-b.p_value;}).slice(0,10);
        if(!st.length)return;
        var col=srcColors[src]||'#888';
        var maxZ=Math.max.apply(null,st.map(function(t){return Math.abs(t.z_score);}));
        var rows=st.map(function(t){
          var pct=Math.abs(t.z_score)/maxZ*100;
          var bc=t.z_score>=0?'#e74c3c':'#3498db';
          var nm=t.term_name.length>30?t.term_name.slice(0,28)+'\u2026':t.term_name;
          var ps=t.p_value<0.001?t.p_value.toExponential(1):'\u00b1'+t.p_value.toFixed(2);
          return '<tr style="border-bottom:1px solid var(--bd);">'+
            '<td style="padding:3px 6px;font-size:10px;color:var(--tx2);max-width:180px;">'+nm+'</td>'+
            '<td style="padding:3px 6px;width:120px;"><div style="height:10px;width:'+pct+'%;background:'+bc+';border-radius:2px;min-width:2px;"></div></td>'+
            '<td style="padding:3px 6px;font-size:10px;color:'+col+';text-align:right;white-space:nowrap;">'+ps+'</td></tr>';
        }).join('');
        var p=document.createElement('div');p.className='panel';p.style.flex='1';
        p.innerHTML='<div class="ph"><h2 style="color:'+col+';">'+src+'</h2></div>'+
          '<table style="width:100%;border-collapse:collapse;"><thead><tr style="border-bottom:1px solid var(--bd);">'+
          '<th style="padding:3px 6px;text-align:left;font-size:9px;color:var(--tx2);">Terme</th>'+
          '<th style="padding:3px 6px;font-size:9px;color:var(--tx2);">z-score</th>'+
          '<th style="padding:3px 6px;text-align:right;font-size:9px;color:var(--tx2);">p-value</th>'+
          '</tr></thead><tbody>'+rows+'</tbody></table>';
        panelsDiv.appendChild(p);
      });
    }

    // Bar chart Z-score
    var barDiv=document.getElementById('goZChartReplace');
    if(!barDiv){
      barDiv=document.createElement('div');barDiv.id='goZChartReplace';
      barDiv.className='panel';barDiv.style.marginBottom='14px';
      if(panelsDiv&&panelsDiv.parentNode) panelsDiv.parentNode.insertBefore(barDiv,panelsDiv.nextSibling);
    }
    // Vider complètement avant de reconstruire (changement de contraste)
    barDiv.innerHTML='';
    var sorted=terms.slice().sort(function(a,b){return b.z_score-a.z_score;});
    var N=Math.min(sorted.length,20);
    var dpr=devicePixelRatio||1;
    var lW=280,rW=60,tH=44,barH=28,PW=460,bH=40;
    var W=lW+PW+rW,H=tH+N*barH+bH;
    var cv=document.createElement('canvas');cv.id='goZCanvasNative';
    cv.width=W*dpr;cv.height=H*dpr;cv.style.cssText='display:block;width:100%;max-width:'+W+'px;';
    var ctx=cv.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,W,H);
    var maxZ=Math.max.apply(null,sorted.map(function(t){return Math.abs(t.z_score);})) || 0.5;
    var scale=PW*0.44/maxZ;
    var x0=lW+PW/2;
    // Grille verticale
    ctx.strokeStyle='#21262d';ctx.lineWidth=1;
    var step=maxZ>0.4?0.2:0.1;
    for(var gv=-maxZ;gv<=maxZ+0.001;gv+=step){
      gv=Math.round(gv*100)/100;
      var gx=x0+gv*scale;
      ctx.beginPath();ctx.moveTo(gx,tH);ctx.lineTo(gx,tH+N*barH);ctx.stroke();
      ctx.fillStyle='#8b949e';ctx.font='9px sans-serif';ctx.textAlign='center';
      ctx.fillText(gv.toFixed(2),gx,tH+N*barH+16);
    }
    // Ligne zéro
    ctx.strokeStyle='#6e7681';ctx.lineWidth=2;
    ctx.beginPath();ctx.moveTo(x0,tH-6);ctx.lineTo(x0,tH+N*barH+6);ctx.stroke();
    // Titre
    ctx.fillStyle='#c9d1d9';ctx.font='bold 12px sans-serif';ctx.textAlign='center';
    ctx.fillText('Z-score par terme  —  rouge = up · bleu = down',lW+PW/2,24);
    // Barres et labels
    sorted.slice(0,N).forEach(function(t,i){
      var y=tH+i*barH;
      var z=t.z_score;
      var bc=z>=0?'#e74c3c':'#3498db';
      var bW=Math.abs(z)*scale;
      var bX=z>=0?x0:x0-bW;
      // Fond léger de la ligne
      if(i%2===0){ctx.fillStyle='rgba(255,255,255,0.02)';ctx.fillRect(0,y,W,barH);}
      // Barre
      ctx.fillStyle=bc+'dd';ctx.fillRect(bX,y+4,bW,barH-8);
      // Label terme (à gauche)
      var nm=t.term_name.length>38?t.term_name.slice(0,36)+'\u2026':t.term_name;
      ctx.fillStyle='#c9d1d9';ctx.font='10px sans-serif';ctx.textAlign='right';
      ctx.fillText(nm,lW-8,y+barH/2+3);
      // Badge source (à droite)
      var sc=srcColors[t.source]||'#888';
      ctx.fillStyle=sc+'33';
      ctx.beginPath();
      ctx.roundRect ? ctx.roundRect(lW+PW+4,y+6,46,barH-12,3) :
                      ctx.rect(lW+PW+4,y+6,46,barH-12);
      ctx.fill();
      ctx.fillStyle=sc;ctx.font='bold 8px sans-serif';ctx.textAlign='center';
      ctx.fillText(t.source.replace('GO:',''),lW+PW+27,y+barH/2+3);
      // Valeur z au bout de la barre
      ctx.fillStyle=bc;ctx.font='bold 9px sans-serif';
      ctx.textAlign=z>=0?'left':'right';
      var valX=z>=0?bX+bW+3:bX-3;
      ctx.fillText(z.toFixed(3),valX,y+barH/2+3);
    });
    ctx.fillStyle='#8b949e';ctx.font='10px sans-serif';ctx.textAlign='center';
    ctx.fillText('Z-score (positif = up-regulated)',lW+PW/2,tH+N*barH+30);
    barDiv.innerHTML='<div class="ph"><h2>Z-score par terme \u2014 rouge=up \u00b7 bleu=down</h2>'+
      '<button id="goZExportBtn" class="exp-btn" style="font-size:9px;padding:2px 7px;">\u2b07 PNG</button></div>';
    barDiv.appendChild(cv);
    var eb=document.getElementById('goZExportBtn');
    if(eb) eb.onclick=function(){exportCanvas('goZCanvasNative','proteogen_GO_zscore');};
  };

  var _goInitDone = false;
  window.initGO=function(){
    // Recréer les boutons si goBtns est vide (une seule fois)
    var wrap=document.getElementById('goBtns');
    var needsInit = wrap && wrap.children.length === 0;
    if(needsInit){
      CONTRASTS.forEach(function(c,i){
        var btn=document.createElement('button');
        btn.className='cbtn'+(i===0?' active':'');
        var parts=c.split('_vs_');
        btn.textContent=parts[0].replace(/_/g,' ')+' vs '+
                        (parts[1]||'').replace(/_/g,' ');
        btn.dataset.c=c;
        btn.onclick=function(){
          document.querySelectorAll('#goBtns .cbtn').forEach(function(b){b.classList.remove('active');});
          btn.classList.add('active');
          renderGO(c);
          // Synchroniser bubble plot et réseau GO
          var bblSel=document.getElementById('goBubbleContrast');
          if(bblSel){bblSel.value=c;if(typeof drawGOBubble==='function')drawGOBubble();
                     if(typeof renderGOTable==='function')renderGOTable();}
          var netSel=document.getElementById('goNetContrast');
          if(netSel){netSel.value=c;if(typeof drawGONetwork==='function')drawGONetwork();}
        };
        wrap.appendChild(btn);
      });
    }
    // Ne rappeler renderGO qu'au premier init
    if(!_goInitDone){
      _goInitDone = true;
      if(typeof currentGOContrast!=='undefined'&&currentGOContrast) renderGO(currentGOContrast);
      else renderGO(CONTRASTS&&CONTRASTS[0]);
    }
  };

  // Re-exécuter les initialisations critiques
  ['initDepGrid','initProtTable','initOverview'].forEach(function(fn){
    try{ if(typeof window[fn]==='function') window[fn](); }catch(e){ console.warn('Init guard:',fn,e.message); }
  });

  // Forcer protContrast si vide
  var sel=document.getElementById('protContrast');
  if(sel&&sel.options.length===0&&typeof CONTRASTS!=='undefined'){
    CONTRASTS.forEach(function(c){var o=document.createElement('option');o.value=c;
      o.textContent=c.replace(/_vs_/g,' vs ').replace(/_/g,' ');sel.appendChild(o);});
    if(typeof renderProtTable==='function') renderProtTable();
  }
})();
"""

def patch_ranked_cv_about_js(html):
    """Injecte drawRankedAbundance(), drawCVPlot()/initCVPlot(), renderAbout(),
    renderQCReport() et la recherche globale avant la dernière </script>.
    """
    QC_REPORT_JS = """\
/* ── QC REPORT ── */
function renderQCReport(){
  const el=document.getElementById('qcReportPanel');
  if(!el||!QC_REPORT)return;
  const d=QC_REPORT;
  const S={'OK':'#3fb950','WARN':'#d29922','FAIL':'#f85149'};
  const SB={'OK':'#3fb95022','WARN':'#d2992222','FAIL':'#f8514922'};
  let html='';
  // Bandeau résumé global
  const sm=d.summary;
  const gs=sm.global_status;
  html+='<div style="display:flex;align-items:center;gap:12px;padding:10px 14px;'+
        'background:'+SB[gs]+';border:1px solid '+S[gs]+';border-radius:6px;margin-bottom:14px;">'+
        '<span style="font-size:20px;">'+(gs==='OK'?'✓':gs==='WARN'?'⚠':'✗')+'</span>'+
        '<div><b style="color:'+S[gs]+';font-size:13px;">Statut global : '+gs+'</b>'+
        '<div style="font-size:11px;color:var(--tx2);margin-top:2px;">'+
        sm.n_ok+' OK &nbsp;'+sm.n_warn+' WARN &nbsp;'+sm.n_fail+' FAIL</div></div></div>';
  // Sections
  html+='<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:10px;">';
  d.sections.forEach(function(sec){
    html+='<div class="panel"><div class="ph"><h2>'+sec.title+'</h2></div>';
    html+='<table style="width:100%;border-collapse:collapse;font-size:11px;">';
    sec.metrics.forEach(function(m){
      html+='<tr style="border-bottom:1px solid var(--bd);">'+
        '<td style="padding:5px 6px;color:var(--tx2);">'+m.label+'</td>'+
        '<td style="padding:5px 6px;font-weight:700;color:var(--tx);">'+m.value+
        (m.unit?' <span style="font-weight:400;color:var(--tx2);font-size:10px;">'+m.unit+'</span>':'')+
        '</td>'+
        '<td style="padding:5px 6px;text-align:center;">'+
        '<span style="display:inline-block;padding:1px 8px;border-radius:10px;font-size:10px;font-weight:700;'+
        'background:'+SB[m.status]+';color:'+S[m.status]+';">'+m.status+'</span></td>'+
        '<td style="padding:5px 6px;font-size:10px;color:var(--tx2);">'+
        (m.detail||'')+'</td></tr>';
    });
    html+='</table></div>';
  });
  html+='</div>';
  el.innerHTML=html;
}\
"""

    GLOBAL_SEARCH_JS = """\
/* ── RECHERCHE GLOBALE ── */
(function(){
  var _gsInput=null;
  function buildGlobalSearch(){
    if(document.getElementById('globalSearchInput'))return;
    // Créer la barre de recherche dans le header/nav
    var nav=document.querySelector('nav')||document.querySelector('header');
    if(!nav)return;
    var wrapper=document.createElement('div');
    wrapper.style.cssText='display:inline-flex;align-items:center;gap:6px;margin-left:16px;position:relative;';
    wrapper.innerHTML=
      '<input id="globalSearchInput" type="text" placeholder="Rechercher un gène / accession..." '+
      'style="background:var(--bg2);border:1px solid var(--bd);color:var(--tx);border-radius:20px;'+
      'padding:4px 12px;font-size:11px;width:220px;outline:none;" autocomplete="off">'+
      '<div id="globalSearchResults" style="display:none;position:absolute;top:28px;left:0;'+
      'background:var(--bg2);border:1px solid var(--bd);border-radius:6px;width:380px;'+
      'max-height:320px;overflow-y:auto;z-index:9999;box-shadow:0 4px 16px rgba(0,0,0,.4);"></div>';
    nav.appendChild(wrapper);
    _gsInput=document.getElementById('globalSearchInput');
    _gsInput.addEventListener('input',debounce(runGlobalSearch,250));
    _gsInput.addEventListener('keydown',function(e){
      if(e.key==='Escape'){closeGlobalSearch();}
    });
    document.addEventListener('click',function(e){
      if(!wrapper.contains(e.target))closeGlobalSearch();
    });
  }

  function debounce(fn,ms){var t;return function(){clearTimeout(t);t=setTimeout(fn,ms);};}

  function closeGlobalSearch(){
    var r=document.getElementById('globalSearchResults');
    if(r){r.style.display='none';r.innerHTML='';}
  }

  function runGlobalSearch(){
    var q=(_gsInput.value||'').toLowerCase().trim();
    var res=document.getElementById('globalSearchResults');
    if(!res)return;
    if(q.length<2){res.style.display='none';res.innerHTML='';return;}

    var hits=[];

    // Volcano / Comparaison (VOLCANO contient tous les contrastes)
    Object.entries(VOLCANO||{}).forEach(function(e){
      var contrast=e[0],pts=e[1];
      pts.forEach(function(p){
        if((p.Genes||'').toLowerCase().includes(q)||(p.name||'').toLowerCase().includes(q)){
          hits.push({section:'Volcano',contrast:contrast,label:(p.Genes||p.name),
            detail:'LFC='+p.d+' p='+p.p, page:'volcano',
            action:function(){showPage('volcano');}});
        }
      });
    });

    // Heatmap (WGCNA_ROWS)
    (WGCNA_ROWS||[]).forEach(function(r){
      if((r.g||'').toLowerCase().includes(q)){
        hits.push({section:'WGCNA',contrast:r.mod,label:r.g,
          detail:'HubScore='+r.hub, page:'wgcna',
          action:function(){showPage('wgcna');}});
      }
    });

    // Intersections UpSet
    (INTERSECTIONS||[]).forEach(function(r){
      if((r.g||'').toLowerCase().includes(q)||(r.pg||'').toLowerCase().includes(q)){
        hits.push({section:'Intersections',contrast:'',label:(r.g||r.pg),
          detail:'Nb contrastes='+r.nb, page:'intersect',
          action:function(){showPage('intersect');}});
      }
    });

    // GO terms
    Object.entries(GO_DATA||{}).forEach(function(e){
      var contrast=e[0],terms=e[1];
      terms.forEach(function(t){
        if((t.term_name||'').toLowerCase().includes(q)){
          hits.push({section:'GO',contrast:contrast,label:t.term_name,
            detail:'p='+t.p_value, page:'gopage',
            action:function(){showPage('gopage');}});
        }
      });
    });

    // Dédupliquer et limiter à 30
    var seen=new Set();
    hits=hits.filter(function(h){
      var k=h.section+'|'+h.label;
      if(seen.has(k))return false;
      seen.add(k);return true;
    }).slice(0,30);

    if(!hits.length){
      res.innerHTML='<div style="padding:10px 14px;color:var(--tx2);font-size:11px;">Aucun résultat pour "'+q+'"</div>';
      res.style.display='block';
      return;
    }

    var COL={'Volcano':'#388bfd','WGCNA':'#3fb950','Intersections':'#bc8cff','GO':'#d29922'};
    var html='<div style="padding:6px 10px;font-size:10px;color:var(--tx2);border-bottom:1px solid var(--bd);">'+
      hits.length+' résultat'+(hits.length>1?'s':'')+' pour "'+q+'"</div>';
    hits.forEach(function(h, i){
      var col=COL[h.section]||'#888';
      html+='<div class="gs-hit" data-idx="'+i+'" style="padding:7px 12px;cursor:pointer;border-bottom:1px solid var(--bd);'+
        'display:flex;gap:10px;align-items:center;">'+
        '<span style="min-width:80px;font-size:9px;font-weight:700;padding:2px 6px;border-radius:10px;'+
        'background:'+col+'22;color:'+col+'">'+h.section+'</span>'+
        '<div><div style="font-size:11px;font-weight:700;color:var(--tx);">'+h.label+'</div>'+
        '<div style="font-size:10px;color:var(--tx2);">'+(h.contrast?h.contrast+' — ':'')+h.detail+'</div></div>'+
        '</div>';
    });
    res.innerHTML=html;
    res.style.display='block';

    // Attacher les clics ET hover après injection
    res.querySelectorAll('.gs-hit').forEach(function(el,i){
      el.addEventListener('mouseover',function(){el.style.background='var(--bg)';});
      el.addEventListener('mouseout', function(){el.style.background='';});
      el.addEventListener('click',function(){
        hits[i].action();
        closeGlobalSearch();
        _gsInput.value='';
      });
    });
  }

  document.addEventListener('DOMContentLoaded',buildGlobalSearch);
  // Fallback si DOMContentLoaded déjà passé
  if(document.readyState==='complete'||document.readyState==='interactive'){
    setTimeout(buildGlobalSearch,100);
  }
})();\
"""

    pos = html.rfind('</script>')
    if pos == -1:
        print("  [warn] patch_ranked_cv_about_js : </script> introuvable.")
        return html
    inject = "\n\n" + DRAW_RANKED_CV_JS + "\n\n" + QC_REPORT_JS + "\n\n" + GLOBAL_SEARCH_JS + "\n\n" + INIT_GUARD_JS + "\n"
    html = html[:pos] + inject + html[pos:]
    print("  patch_ranked_cv_about_js : ranked + CV + about + QC report + recherche globale injectés.")
    return html


# Bloc HTML pour les nouvelles figures dans l'onglet QC + panneau about dans overview
RANKED_CV_HTML = """\
<div class="row2">
  <div class="panel">
    <div class="ph"><h2>Ranked protein abundance</h2><button class="exp-btn" onclick="exportCanvas('rankedCanvas','proteogen_ranked')" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button></div>
    <canvas id="rankedCanvas" style="display:block;width:100%;"></canvas>
    <div id="rankedTip" style="position:fixed;display:none;background:var(--bg2);border:1px solid var(--bd);border-radius:6px;padding:6px 10px;font-size:11px;color:var(--tx);pointer-events:none;z-index:999;"></div>
  </div>
  <div class="panel">
    <div class="ph"><h2>CV intra-condition</h2><div style="display:flex;gap:6px;align-items:center;"><select id="cvCond" style="background:var(--bg2);border:1px solid var(--bd);color:var(--tx);border-radius:4px;padding:3px 8px;font-size:11px;"></select><button class="exp-btn" onclick="exportCanvas('cvCanvas','proteogen_CV_'+document.getElementById('cvCond').value)" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button></div></div>
    <canvas id="cvCanvas" style="display:block;width:100%;"></canvas>
    <div id="cvTip" style="position:fixed;display:none;background:var(--bg2);border:1px solid var(--bd);border-radius:6px;padding:6px 10px;font-size:11px;color:var(--tx);pointer-events:none;z-index:999;"></div>
  </div>
</div>\
"""

ABOUT_HTML = """\
<div class="panel" style="margin-bottom:14px;">
  <div class="ph"><h2>&#8505; À propos de ce dashboard</h2></div>
  <div id="aboutPanel"></div>
</div>\
"""

# Palette de couleurs conditions : correspond aux couleurs pastel du TIFF R/pheatmap.
# L'ordre de fallback suit l'index de CONDITIONS ; le matching par nom prend la priorité.
COND_COLORS_JS = """\
// COND_COLORS déclaré sur window pour portée globale garantie
// (remplace const qui serait limité au script en cas de crash partiel)
window.COND_COLORS = window.COND_COLORS || {};
var _COND_PALETTE=['#f5e99a','#e8928a','#c3b8e8','#8ecfc9','#f0b87a','#a8d8a8'];
(CONDITIONS||[]).forEach((c,i)=>{
  if(c.match(/Monoculture/i))            window.COND_COLORS[c]='#8ecfc9';
  else if(c.match(/Polyculture[._]?1/i)) window.COND_COLORS[c]='#f5e99a';
  else if(c.match(/Polyculture[._]?2/i)) window.COND_COLORS[c]='#c3b8e8';
  else if(c.match(/Polyculture[._]?3/i)) window.COND_COLORS[c]='#e8928a';
  else window.COND_COLORS[c]=_COND_PALETTE[i%_COND_PALETTE.length];
});\
"""

# Nouvelle fonction drawHeatmap :
#   • Palette Blue→White→Red, échelle ±4 (identique au TIFF)
#   • En-têtes de condition en blocs colorés pleine hauteur
#   • Sidebar cluster avec numéro en texte blanc
#   • Séparateurs horizontaux entre clusters
#   • Légende Z-score dégradée à droite
#   • Fond blanc, texte sombre (rendu proche du TIFF R)
DRAW_HEATMAP_JS = """\
function buildCondLegend(){
  const el=document.getElementById('condLegend');
  if(!el)return;
  el.innerHTML='';
  CONDITIONS.forEach(cond=>{
    const col=COND_COLORS[cond]||'#888';
    const span=document.createElement('span');
    span.style.cssText='display:inline-flex;align-items:center;gap:4px;font-size:10px;';
    span.innerHTML='<span style="display:inline-block;width:10px;height:10px;border-radius:2px;background:'+col+';"></span>'+
      '<span>'+cond.replace(/Polyculture[._]?/i,'Poly').replace(/Monoculture/i,'Mono')+'</span>';
    el.appendChild(span);
  });
}
function drawHeatmap(){
  const clF=document.getElementById('hmCluster').value;
  const nMax=parseInt(document.getElementById('hmN').value);
  const srt=document.getElementById('hmSort').value;
  const search=(document.getElementById('hmSearch').value||'').toLowerCase().trim();
  let data=[...HEATMAP];
  if(clF!=='all')data=data.filter(r=>r.cl===clF);
  if(search)data=data.filter(r=>r.g.toLowerCase().includes(search));
  data=srt==='q'?data.sort((a,b)=>a.q-b.q):data.sort((a,b)=>a.cl.localeCompare(b.cl)||a.q-b.q);
  data=data.slice(0,nMax);
  document.getElementById('hmInfo').textContent=data.length+' protéines';
  // Guard : aucune protéine significative (ex: 2 conditions, pas d'ANOVA)
  if(!data.length){
    const c=document.getElementById('hmCanvas');
    if(c){const ctx=c.getContext('2d');c.width=400;c.height=80;
      ctx.fillStyle='#8b949e';ctx.font='13px sans-serif';ctx.textAlign='center';
      ctx.fillText('Aucune protéine significative à afficher',200,45);}
    return;
  }

  const nR=data.length,nC=SAMPLE_COLS.length;
  const cW=22,cH=8,lW=120,lH=130,clW=18,lgW=80;
  const W=clW+lW+nC*cW+lgW+30,H=lH+nR*cH+60;
  const canvas=document.getElementById('hmCanvas');
  const dpr=devicePixelRatio||1;
  canvas.width=W*dpr;canvas.height=H*dpr;
  canvas.style.width=W+'px';canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.fillStyle='#ffffff';ctx.fillRect(0,0,W,H);

  // Palette Blue->White->Red, echelle -4..+4 (identique TIFF R)
  function zColor(v){
    const clamped=Math.max(-4,Math.min(4,typeof v==='number'?v:0));
    if(clamped>=0){
      const f=clamped/4;
      return'rgb(255,'+Math.round(255*(1-f))+','+Math.round(255*(1-f))+')';
    } else {
      const f=-clamped/4;
      return'rgb('+Math.round(255*(1-f))+','+Math.round(255*(1-f))+',255)';
    }
  }

  // Sidebar cluster — généré dynamiquement selon les clusters présents
  // (le nombre de clusters dépend du paramètre choisi dans le pipeline).
  const _clPalette=['#5B8FF9','#F6BD16','#5AD8A6','#E8684A','#9270CA',
                    '#FF9D4D','#269A99','#FF99C3','#6DC8EC','#FF6B6B'];
  const CLC={'NA':'#aaaaaa'};
  const CLText={'NA':'NA'};
  (function(){
    var seen={};
    HEATMAP.forEach(function(r){var c=r.cl; if(c&&c!=='NA')seen[c]=true;});
    var keys=Object.keys(seen).sort();
    keys.forEach(function(k,i){
      CLC[k]=_clPalette[i%_clPalette.length];
      // 'Cluster_4' -> 'Cl.4'
      var m=String(k).match(/(\\d+)$/);
      CLText[k]='Cl.'+(m?m[1]:(i+1));
    });
  })();

  // Titre
  const totalN=HEATMAP.length;
  ctx.fillStyle='#111111';ctx.font='bold 11px sans-serif';ctx.textAlign='center';
  ctx.fillText(
    'Analyse par Clusters - Protéines Significatives (ANOVA, n='+totalN+')',
    (clW+lW+nC*cW)/2, 18
  );

  // En-têtes de condition : blocs colorés pleins (style TIFF)
  let xCursor=clW+lW;
  CONDITIONS.forEach(cond=>{
    const cols=COND_MAP[cond]||[];
    const condW=cols.length*cW;
    const col=COND_COLORS[cond]||'#cccccc';
    ctx.fillStyle=col;
    ctx.fillRect(xCursor,lH-28,condW,22);
    ctx.fillStyle='#222222';ctx.font='bold 9px sans-serif';ctx.textAlign='center';
    ctx.fillText(
      cond.replace(/Polyculture[._]?/i,'Poly').replace(/Monoculture/i,'Mono'),
      xCursor+condW/2, lH-12
    );
    // Séparateur vertical entre groupes
    if(xCursor>clW+lW){
      ctx.strokeStyle='rgba(255,255,255,0.7)';ctx.lineWidth=1.5;
      ctx.beginPath();ctx.moveTo(xCursor,lH-28);ctx.lineTo(xCursor,lH+nR*cH);ctx.stroke();
    }
    xCursor+=condW;
  });

  // Labels colonnes (numéro de réplique, inclinés)
  SAMPLE_COLS.forEach((col,i)=>{
    ctx.save();
    ctx.translate(clW+lW+i*cW+cW/2,lH-32);
    ctx.rotate(-Math.PI/3);
    ctx.fillStyle='#444444';ctx.font='7px sans-serif';ctx.textAlign='right';
    ctx.fillText(col.replace(/^.*_(\\d+)$/,'$1'),0,0);
    ctx.restore();
  });

  // Cellules Z-score
  data.forEach((row,ri)=>{
    const y=lH+ri*cH;
    // Barre latérale cluster
    ctx.fillStyle=CLC[row.cl]||'#888888';
    ctx.fillRect(0,y,clW-2,cH-0.5);
    ctx.fillStyle='#ffffff';ctx.font='bold 6px sans-serif';ctx.textAlign='center';
    ctx.fillText(CLText[row.cl]||'',clW/2,y+cH-1.5);
    // Nom du gène
    ctx.fillStyle='#222222';ctx.font='9px sans-serif';ctx.textAlign='right';
    ctx.fillText(row.g.slice(0,16),clW+lW-5,y+cH-1.5);
    // Valeurs Z-score
    row.v.forEach((val,ci)=>{
      ctx.fillStyle=zColor(val);
      ctx.fillRect(clW+lW+ci*cW,y,cW-0.5,cH-0.5);
    });
  });

  // Séparateurs horizontaux entre clusters
  let prevCl=null;
  data.forEach((row,ri)=>{
    if(prevCl!==null&&row.cl!==prevCl){
      const y=lH+ri*cH;
      ctx.strokeStyle='#ffffff';ctx.lineWidth=2;
      ctx.beginPath();ctx.moveTo(0,y);ctx.lineTo(clW+lW+nC*cW,y);ctx.stroke();
    }
    prevCl=row.cl;
  });

  // Légende Z-score (dégradé vertical, droite)
  const lgX=clW+lW+nC*cW+10;
  const lgBarH=120,lgY=lH;
  const grad=ctx.createLinearGradient(lgX,lgY,lgX,lgY+lgBarH);
  grad.addColorStop(0,'rgb(255,0,0)');
  grad.addColorStop(0.5,'rgb(255,255,255)');
  grad.addColorStop(1,'rgb(0,0,255)');
  ctx.fillStyle=grad;
  ctx.fillRect(lgX,lgY,14,lgBarH);
  ctx.strokeStyle='#aaaaaa';ctx.lineWidth=0.5;
  ctx.strokeRect(lgX,lgY,14,lgBarH);
  ctx.fillStyle='#333333';ctx.font='8px sans-serif';ctx.textAlign='left';
  [[0,4],[0.25,2],[0.5,0],[0.75,-2],[1,-4]].forEach(([f,v])=>{
    const ty=lgY+f*lgBarH;
    ctx.fillText(v,lgX+17,ty+3);
  });
  ctx.fillStyle='#333333';ctx.font='bold 9px sans-serif';ctx.textAlign='center';
  ctx.save();ctx.translate(lgX+7,lgY+lgBarH/2);ctx.rotate(-Math.PI/2);
  ctx.fillText('Z-score',0,0);ctx.restore();
}\
"""


# ──────────────────────────────────────────────────────────────────────────────
# QC REPORT — métriques de qualité calculées côté Python
# ──────────────────────────────────────────────────────────────────────────────
QC_NA_WARN  = 0.20
QC_NA_FAIL  = 0.40
QC_R_WARN   = 0.95
QC_R_FAIL   = 0.90
QC_CV_WARN  = 20.0
QC_IMP_WARN = 0.30

def _qc_status(val, warn_thr, fail_thr, higher_is_bad=True):
    if higher_is_bad:
        if val >= fail_thr: return "FAIL"
        if val >= warn_thr: return "WARN"
    else:
        if val <= fail_thr: return "FAIL"
        if val <= warn_thr: return "WARN"
    return "OK"


def build_qc_report(stats_path, df_comp, sample_cols, cond_map):
    """Calcule les métriques QC : NA rates, corrélations, CV, imputation."""
    from itertools import combinations as _comb

    cond_to_cols = {}
    for col, cond in cond_map.items():
        if col in df_comp.columns:
            cond_to_cols.setdefault(cond, []).append(col)

    sections = []
    all_statuses = []

    # ── Section 1 : Valeurs manquantes ───────────────────────────────────────
    metrics_na = []
    try:
        df_brut = pd.read_excel(stats_path, sheet_name="raw_data")
        lfq_cols = [c for c in df_brut.columns if "LFQ" in str(c)]
        if lfq_cols:
            na_rates = df_brut[lfq_cols].isna().mean()
            mean_na = float(na_rates.mean())
            max_na  = float(na_rates.max())
            max_col = str(na_rates.idxmax()).replace("LFQ.intensity.", "")
            s_mean = _qc_status(mean_na, QC_NA_WARN, QC_NA_FAIL)
            s_max  = _qc_status(max_na,  QC_NA_WARN, QC_NA_FAIL)
            all_statuses += [s_mean, s_max]
            metrics_na.append({"label": "Taux NA moyen (global)",
                                "value": round(mean_na*100,1), "unit": "%",
                                "status": s_mean, "detail": ""})
            metrics_na.append({"label": "Taux NA max (echantillon)",
                                "value": round(max_na*100,1), "unit": "%",
                                "status": s_max, "detail": max_col})
            n_quant = [(col, int((df_comp[col]>0).sum()))
                       for col in sample_cols if col in df_comp.columns]
            min_n = min(n for _, n in n_quant)
            max_n = max(n for _, n in n_quant)
            min_col = next(col for col, n in n_quant if n == min_n)
            s_q = "OK" if min_n > len(df_comp)*0.7 else "WARN"
            all_statuses.append(s_q)
            metrics_na.append({"label": "Proteines quantifiees min/max",
                                "value": str(min_n)+" / "+str(max_n), "unit": "prot.",
                                "status": s_q, "detail": "min: "+min_col})
    except Exception as e:
        metrics_na.append({"label": "raw_data", "value": "N/A", "unit": "",
                            "status": "WARN", "detail": str(e)[:60]})
    sections.append({"title": "Valeurs manquantes", "metrics": metrics_na})

    # ── Section 2 : Corrélations inter-réplicats ──────────────────────────────
    metrics_r = []
    for cond, cols in cond_to_cols.items():
        if len(cols) < 2:
            continue
        rs = []
        for c1, c2 in _comb(cols, 2):
            sub = df_comp[[c1,c2]].dropna()
            sub = sub[(sub[c1]>0)&(sub[c2]>0)]
            if len(sub) > 1:
                rs.append(float(np.corrcoef(sub[c1], sub[c2])[0,1]))
        if not rs:
            continue
        min_r  = round(min(rs), 4)
        mean_r = round(float(np.mean(rs)), 4)
        s = _qc_status(min_r, QC_R_WARN, QC_R_FAIL, higher_is_bad=False)
        all_statuses.append(s)
        n_w = sum(1 for r in rs if r < QC_R_WARN)
        metrics_r.append({"label": cond,
                           "value": str(mean_r)+" (min "+str(min_r)+")",
                           "unit": "r Pearson", "status": s,
                           "detail": str(n_w)+" paires <"+str(QC_R_WARN) if n_w else ""})
    sections.append({"title": "Correlations inter-replicats", "metrics": metrics_r})

    # ── Section 3 : CV ────────────────────────────────────────────────────────
    metrics_cv = []
    for cond, cols in cond_to_cols.items():
        if len(cols) < 2:
            continue
        sub = df_comp[cols].replace(0, np.nan)
        cv  = (sub.std(axis=1)/sub.mean(axis=1)*100).dropna()
        med = round(float(cv.median()),1)
        pct = round(float((cv>QC_CV_WARN).mean()*100),1)
        s = _qc_status(pct, 5.0, 15.0)
        all_statuses.append(s)
        metrics_cv.append({"label": cond, "value": med, "unit": "% CV median",
                            "status": s,
                            "detail": str(pct)+"% prot. CV>"+str(QC_CV_WARN)+"%"})
    sections.append({"title": "Coefficient de variation (intra-condition)", "metrics": metrics_cv})

    # ── Section 4 : Imputation ────────────────────────────────────────────────
    metrics_imp = []
    if "imputed" in df_comp.columns:
        # La colonne 'imputed' peut être VRAI/FAUX (texte, format R), booléenne
        # ou 0/1. On normalise en booléen avant de compter.
        imp_col = df_comp["imputed"]
        if imp_col.dtype == object:
            n_imp = int(imp_col.astype(str).str.upper().isin(
                ["VRAI", "TRUE", "1", "OUI", "YES"]).sum())
        else:
            n_imp = int(pd.to_numeric(imp_col, errors="coerce").fillna(0).astype(bool).sum())
        n_tot = int(len(df_comp))
        rate  = round(n_imp/n_tot, 4) if n_tot else 0.0
        s = _qc_status(rate, QC_IMP_WARN, 0.50)
        all_statuses.append(s)
        metrics_imp.append({"label": "Proteines avec imputation",
                             "value": round(rate*100,1), "unit": "%",
                             "status": s, "detail": str(n_imp)+"/"+str(n_tot)})
    if "num_NAs" in df_comp.columns:
        mna = round(float(df_comp["num_NAs"].mean()),2)
        s2 = "OK" if mna < 2 else "WARN"
        all_statuses.append(s2)
        metrics_imp.append({"label": "NA moyen par proteine",
                             "value": mna, "unit": "NA/prot.", "status": s2, "detail": ""})
    if metrics_imp:
        sections.append({"title": "Imputation", "metrics": metrics_imp})

    n_ok   = all_statuses.count("OK")
    n_warn = all_statuses.count("WARN")
    n_fail = all_statuses.count("FAIL")
    return {"sections": sections,
            "summary": {"n_ok": n_ok, "n_warn": n_warn, "n_fail": n_fail,
                        "global_status": "FAIL" if n_fail else ("WARN" if n_warn else "OK")}}



def patch_heatmap_js(html):
    """Remplace dans le HTML final :
    1. COND_COLORS : dans le script qui contient 'const COND_COLORS='
       (script applicatif original, PAS le bloc données JSON).
    2. drawHeatmap() : dans le script qui contient 'function drawHeatmap'.
    Les deux peuvent être dans des scripts différents — on les cible séparément.
    """
    import re as _re

    script_re = _re.compile(r'(<script[^>]*>)(.*?)(</script>)', _re.DOTALL)
    scripts = list(script_re.finditer(html))

    # ── 1. Patch COND_COLORS — script contenant 'const COND_COLORS=' ──────────
    cc_match = None
    for m in scripts:
        body = m.group(2)
        if 'const COND_COLORS=' in body and 'const CONTRASTS=' not in body:
            cc_match = m
            break

    if cc_match:
        body = cc_match.group(2)
        # Injecter le guard Chart.js + neutralisation initGO/renderGO au début du script
        # initGO/renderGO créent des Chart sur goZChart/goSrcChart supprimés → crash
        CHART_GUARD = (
            '/* Guard: Chart.js + COND_COLORS global + neutralisation initGO */\n'
            'if(typeof Chart==="undefined"){\n'
            '  function Chart(){};\n'
            '  Chart.defaults={color:"#8b949e",borderColor:"rgba(139,148,158,0.15)",font:{size:11}};\n'
            '  Chart.register=function(){};\n'
            '}\n'
            'if(Chart.defaults&&!Chart.defaults.font){Chart.defaults.font={size:11};}\n'
            'window.COND_COLORS=window.COND_COLORS||{};\n'
            '// initGO/renderGO utilisent goZChart/goSrcChart supprimés du DOM\n'
            '// On les neutralise avant que le bloc INIT les appelle\n'
            'window.initGO=function(){};\n'
            'window.renderGO=function(){};\n'
        )
        body = CHART_GUARD + body
        p1 = _re.compile(
            r'const COND_COLORS=\{\};.*?CONDITIONS\.forEach.*?\}\s*\);',
            _re.DOTALL
        )
        if p1.search(body):
            body = p1.sub(lambda _: COND_COLORS_JS, body, count=1)
            html = html[:cc_match.start(2)] + body + html[cc_match.end(2):]
            # Re-scanner les scripts après modification
            scripts = list(script_re.finditer(html))
            print("  patch_heatmap_js : COND_COLORS patché.")
        else:
            print("  [warn] patch_heatmap_js : pattern COND_COLORS non trouvé.")
    else:
        print("  [warn] patch_heatmap_js : script COND_COLORS introuvable.")

    # ── 2. Patch drawHeatmap — neutralisation + injection à la fin du script ──
    # Le regex de fin de fonction est fragile (accolades dans les lambdas/strings).
    # Stratégie : renommer l'ancienne fonction _origDrawHeatmap, injecter la nouvelle.
    hm_match = None
    for m in scripts:
        body = m.group(2)
        if 'function drawHeatmap' in body and 'const CONTRASTS=' not in body:
            hm_match = m
            break

    if hm_match:
        body = hm_match.group(2)
        # Renommer l'ancienne fonction pour la neutraliser
        body = body.replace('function drawHeatmap(){', 'function _origDrawHeatmap_unused(){', 1)
        # Injecter la nouvelle drawHeatmap à la fin du script
        body = body + '\n\n' + DRAW_HEATMAP_JS + '\n'
        html = html[:hm_match.start(2)] + body + html[hm_match.end(2):]
        print("  patch_heatmap_js : drawHeatmap() patché (neutralisation + injection).")
    else:
        print("  [warn] patch_heatmap_js : script drawHeatmap introuvable.")

    return html




# ──────────────────────────────────────────────────────────────────────────────
# Patch JS — Violin plot et Missing Values (onglet QC)
# ──────────────────────────────────────────────────────────────────────────────

# drawViolinPlot : distributions log2 LFQ post-imputation
#   • KDE approchée par histogramme lissé (fenêtre gaussienne simple, côté JS)
#   • Whiskers = min/max, boîte = IQR, trait médiane rouge
#   • Groupés et colorés par condition (palette pastel TIFF)
DRAW_VIOLIN_JS = """\
function drawViolinPlot(){
  if(!VIOLIN_DATA||!VIOLIN_DATA.samples)return;
  const canvas=document.getElementById('violinCanvas');
  if(!canvas)return;
  const samples=SAMPLE_COLS;
  const binCenters=VIOLIN_DATA.bin_centers;
  const nS=samples.length;
  const dpr=devicePixelRatio||1;
  const lW=110,bW=28,bPad=3,topH=26,plotH=200,botH=50,lgH=20;
  const W=lW+nS*(bW+bPad)+20, H=topH+plotH+botH+lgH;
  canvas.width=W*dpr; canvas.height=H*dpr;
  canvas.style.width=W+'px'; canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.clearRect(0,0,W,H);

  const yMin=VIOLIN_DATA.vmin, yMax=VIOLIN_DATA.vmax;
  function toY(v){ return topH+plotH*(1-(v-yMin)/(yMax-yMin)); }

  // Grille et axe Y
  ctx.strokeStyle='#21262d'; ctx.lineWidth=1;
  for(let v=Math.ceil(yMin);v<=yMax;v+=2){
    const y=toY(v);
    ctx.beginPath(); ctx.moveTo(lW,y); ctx.lineTo(W-10,y); ctx.stroke();
    ctx.fillStyle='#8b949e'; ctx.font='8px sans-serif'; ctx.textAlign='right';
    ctx.fillText(v, lW-4, y+3);
  }
  ctx.save(); ctx.translate(12, topH+plotH/2); ctx.rotate(-Math.PI/2);
  ctx.fillStyle='#444'; ctx.font='bold 9px sans-serif'; ctx.textAlign='center';
  ctx.fillText('log\u2082 LFQ intensity', 0, 0); ctx.restore();

  // Titre
  ctx.fillStyle='#111'; ctx.font='bold 11px sans-serif'; ctx.textAlign='center';
  ctx.fillText('Distribution des intensités LFQ par échantillon (post-imputation)',
    lW+(W-lW)/2, 18);

  // Lissage gaussien simple (sigma=1 bin)
  function smooth(hist){
    const k=[0.25,0.5,0.25];
    return hist.map((_,i)=>k.reduce((s,w,ki)=>s+w*(hist[i+ki-1]||0),0));
  }

  // Fond coloré par groupe de condition
  let xCur=lW;
  CONDITIONS.forEach(cond=>{
    const cols=(COND_MAP[cond]||[]);
    const gW=cols.length*(bW+bPad);
    ctx.fillStyle=COND_COLORS[cond]||'#ccc';
    ctx.globalAlpha=0.15; ctx.fillRect(xCur,topH,gW,plotH); ctx.globalAlpha=1;
    if(xCur>lW){
      ctx.strokeStyle='#ddd'; ctx.lineWidth=1;
      ctx.beginPath(); ctx.moveTo(xCur,topH); ctx.lineTo(xCur,topH+plotH); ctx.stroke();
    }
    ctx.fillStyle='#c9d1d9'; ctx.font='bold 9px sans-serif'; ctx.textAlign='center';
    ctx.fillText(
      cond.replace(/Polyculture[._]?/i,'Poly').replace(/Monoculture/i,'Mono'),
      xCur+gW/2, topH+plotH+botH+14
    );
    xCur+=gW;
  });

  // Violins
  samples.forEach((samp,si)=>{
    const d=VIOLIN_DATA.samples[samp];
    if(!d)return;
    const xC=lW+(si+0.5)*(bW+bPad);
    const col=COND_COLORS[d.cond]||'#cccccc';
    const hist=smooth(d.hist);
    const hMax=Math.max(...hist)||1;
    const halfW=(bW/2-2);
    // Contour violin
    ctx.beginPath();
    binCenters.forEach((bc,bi)=>{
      const hw=hist[bi]/hMax*halfW;
      const y=toY(bc);
      bi===0?ctx.moveTo(xC-hw,y):ctx.lineTo(xC-hw,y);
    });
    [...binCenters].reverse().forEach((bc,bi)=>{
      const bi2=binCenters.length-1-bi;
      ctx.lineTo(xC+hist[bi2]/hMax*halfW, toY(bc));
    });
    ctx.closePath();
    ctx.fillStyle=col+'bb'; ctx.fill();
    ctx.strokeStyle=col; ctx.lineWidth=1; ctx.stroke();
    // Boîte IQR + médiane + whiskers
    const [q0,q1,q2,q3,q4]=d.q;
    const yQ1=toY(q1),yQ3=toY(q3),yMed=toY(q2);
    const bxW=8;
    ctx.fillStyle='rgba(255,255,255,0.85)'; ctx.fillRect(xC-bxW/2,yQ3,bxW,yQ1-yQ3);
    ctx.strokeStyle='#444'; ctx.lineWidth=1.2; ctx.strokeRect(xC-bxW/2,yQ3,bxW,yQ1-yQ3);
    ctx.strokeStyle='#cc0000'; ctx.lineWidth=2;
    ctx.beginPath(); ctx.moveTo(xC-bxW/2,yMed); ctx.lineTo(xC+bxW/2,yMed); ctx.stroke();
    ctx.strokeStyle='#555'; ctx.lineWidth=1;
    ctx.beginPath();
    ctx.moveTo(xC,yQ1); ctx.lineTo(xC,toY(q0));
    ctx.moveTo(xC-3,toY(q0)); ctx.lineTo(xC+3,toY(q0));
    ctx.moveTo(xC,yQ3); ctx.lineTo(xC,toY(q4));
    ctx.moveTo(xC-3,toY(q4)); ctx.lineTo(xC+3,toY(q4));
    ctx.stroke();
    // Numéro réplique
    ctx.save();
    ctx.translate(xC,topH+plotH+8); ctx.rotate(-Math.PI/4);
    ctx.fillStyle='#8b949e'; ctx.font='8px sans-serif'; ctx.textAlign='right';
    ctx.fillText(samp.replace(/^.*_(\\d+)$/,'$1'),0,0);
    ctx.restore();
  });

  // Légende
  const lgY=topH+plotH+botH+2;
  [['— Médiane','#cc0000'],['□ IQR','#8b949e'],['| Min/Max','#8b949e']].forEach(([lbl,col],i)=>{
    ctx.fillStyle=col; ctx.font='8px sans-serif'; ctx.textAlign='left';
    ctx.fillText(lbl, lW+i*110, lgY);
  });
}\
"""

# drawMissingValues : vue compacte — bar chart taux NA par échantillon
#   • Barre horizontale par échantillon, colorée par condition
#   • Taux NA global en titre
#   • Hauteur fixe ~220px, largeur adaptée au nombre d'échantillons
DRAW_MISSING_JS = """\
function drawMissingValues(){
  const canvas=document.getElementById('missingCanvas');
  if(!canvas)return;
  if(!MISSING_DATA||!MISSING_DATA.matrix||!MISSING_DATA.matrix.length){
    canvas.width=400; canvas.height=80;
    const ctx=canvas.getContext('2d');
    ctx.fillStyle='#fff'; ctx.fillRect(0,0,400,80);
    ctx.fillStyle='#888'; ctx.font='13px sans-serif'; ctx.textAlign='center';
    ctx.fillText('Données brutes (raw_data) non disponibles',200,45);
    return;
  }
  const lfqCols=MISSING_DATA.lfq_cols;
  const lfqCond=MISSING_DATA.lfq_cond;
  const naRates=MISSING_DATA.na_rates;
  const summ=MISSING_DATA.summary;
  const nC=lfqCols.length;

  // Mise en page compacte
  const dpr=devicePixelRatio||1;
  const lW=130, rW=90, tH=44, bH=20;
  const BAR_H=12, BAR_PAD=2;
  const PW=340;
  const W=lW+PW+rW, H=tH+(BAR_H+BAR_PAD)*nC+bH;
  canvas.width=W*dpr; canvas.height=H*dpr;
  canvas.style.width=W+'px'; canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.clearRect(0,0,W,H);

  // Titre
  ctx.fillStyle='#c9d1d9'; ctx.font='bold 11px sans-serif'; ctx.textAlign='center';
  ctx.fillText('Valeurs manquantes par échantillon (avant imputation)', lW+PW/2, 16);
  ctx.fillStyle='#8b949e'; ctx.font='9px sans-serif';
  ctx.fillText(
    summ.prot_with_na+' protéines avec ≥1 NA / '+summ.total_prot+
    ' — taux global : '+Math.round(summ.global_na_rate*100)+'%',
    lW+PW/2, 30
  );

  // Axe X (0–100%)
  ctx.strokeStyle='#21262d'; ctx.lineWidth=1;
  [0,25,50,75,100].forEach(pct=>{
    const x=lW+pct/100*PW;
    ctx.beginPath(); ctx.moveTo(x,tH); ctx.lineTo(x,tH+(BAR_H+BAR_PAD)*nC); ctx.stroke();
    ctx.fillStyle='#8b949e'; ctx.font='8px sans-serif'; ctx.textAlign='center';
    ctx.fillText(pct+'%', x, tH+(BAR_H+BAR_PAD)*nC+12);
  });

  // Barres
  naRates.forEach((r,ci)=>{
    const y=tH+ci*(BAR_H+BAR_PAD);
    const col=COND_COLORS[lfqCond[ci]]||'#aaa';
    const pct=r||0;
    const bW=Math.round(pct*PW);

    // Fond gris clair
    ctx.fillStyle='#21262d';
    ctx.fillRect(lW, y, PW, BAR_H);

    // Barre colorée
    if(bW>0){
      ctx.fillStyle=col+'cc';
      ctx.fillRect(lW, y, bW, BAR_H);
    }

    // Label échantillon (gauche)
    ctx.fillStyle='#8b949e'; ctx.font='9px sans-serif'; ctx.textAlign='right';
    const lbl=lfqCols[ci].replace(/^LFQ[.]intensity[.]/,'').slice(-18);
    ctx.fillText(lbl, lW-5, y+BAR_H-4);

    // Valeur % (droite)
    ctx.fillStyle='#c9d1d9'; ctx.font='bold 9px sans-serif'; ctx.textAlign='left';
    ctx.fillText(Math.round(pct*100)+'%', lW+PW+5, y+BAR_H-4);
  });

  // Légende conditions
  const seen=new Set();
  let lgY=tH+4;
  lfqCond.forEach((cond,ci)=>{
    if(seen.has(cond))return; seen.add(cond);
    const col=COND_COLORS[cond]||'#aaa';
    ctx.fillStyle=col;
    ctx.fillRect(lW+PW+50, lgY, 10, 10);
    ctx.fillStyle='#333'; ctx.font='8px sans-serif'; ctx.textAlign='left';
    ctx.fillText(cond.replace(/Polyculture[._]?/i,'Poly').replace(/Monoculture/i,'Mono'),
      lW+PW+63, lgY+9);
    lgY+=16;
  });
}\
"""



# drawMAPlot : A = intensité moyenne log2 LFQ, M = LFC
#   • Même code couleur que le volcano (U=vert, D=rouge, NS=gris)
#   • Ligne M=0 et seuils LFC_THRESH en pointillés
#   • Sélecteur de contraste, tooltip interactif, export PNG
DRAW_MA_JS = """\
function drawMAPlot(){
  if(!MA_DATA||!Object.keys(MA_DATA).length)return;
  const sel=document.getElementById('maContrast');
  const contrast=sel?sel.value:Object.keys(MA_DATA)[0];
  const d=MA_DATA[contrast];
  if(!d||!d.pts)return;

  const canvas=document.getElementById('maCanvas');
  if(!canvas)return;
  const dpr=devicePixelRatio||1;
  const lW=60,rW=20,tH=40,bH=50,PW=700,PH=380;
  const W=lW+PW+rW, H=tH+PH+bH;
  canvas.width=W*dpr; canvas.height=H*dpr;
  canvas.style.width=W+'px'; canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.clearRect(0,0,W,H);

  // Titre
  ctx.fillStyle='#c9d1d9'; ctx.font='bold 11px sans-serif'; ctx.textAlign='center';
  ctx.fillText('MA plot — '+contrast.replace(/_vs_/,' vs ').replace(/_/g,' '), lW+PW/2, 18);
  ctx.fillStyle='#8b949e'; ctx.font='9px sans-serif';
  const up=d.pts.filter(p=>p.s==='U').length;
  const dn=d.pts.filter(p=>p.s==='D').length;
  ctx.fillText('UP: '+up+'  DOWN: '+dn+'  (p<0.05, |LFC|>'+d.lfc_thresh+')', lW+PW/2, 32);

  // Plages
  const pts=d.pts;
  const aVals=pts.map(p=>p.a).filter(v=>v!=null);
  const mVals=pts.map(p=>p.m).filter(v=>v!=null);
  const aMin=Math.min(...aVals), aMax=Math.max(...aVals);
  const mAbs=Math.max(Math.abs(Math.min(...mVals)), Math.abs(Math.max(...mVals)))*1.1;
  function toX(a){ return lW+((a-aMin)/(aMax-aMin))*PW; }
  function toY(m){ return tH+PH/2 - (m/mAbs)*(PH/2); }

  // Grille
  ctx.strokeStyle='#21262d'; ctx.lineWidth=1;
  for(let v=Math.ceil(aMin);v<=aMax;v+=2){
    ctx.beginPath(); ctx.moveTo(toX(v),tH); ctx.lineTo(toX(v),tH+PH); ctx.stroke();
    ctx.fillStyle='#8b949e'; ctx.font='8px sans-serif'; ctx.textAlign='center';
    ctx.fillText(v, toX(v), tH+PH+14);
  }
  const mStep=mAbs>4?2:1;
  for(let v=-Math.ceil(mAbs);v<=Math.ceil(mAbs);v+=mStep){
    ctx.beginPath(); ctx.moveTo(lW,toY(v)); ctx.lineTo(lW+PW,toY(v)); ctx.stroke();
    ctx.fillStyle='#8b949e'; ctx.font='8px sans-serif'; ctx.textAlign='right';
    ctx.fillText(v, lW-4, toY(v)+3);
  }

  // Ligne M=0 (noire), seuils LFC (pointillés orange)
  ctx.strokeStyle='#333'; ctx.lineWidth=1.5;
  ctx.beginPath(); ctx.moveTo(lW,toY(0)); ctx.lineTo(lW+PW,toY(0)); ctx.stroke();
  ctx.setLineDash([5,4]); ctx.strokeStyle='#e67e22'; ctx.lineWidth=1;
  [d.lfc_thresh, -d.lfc_thresh].forEach(v=>{
    ctx.beginPath(); ctx.moveTo(lW,toY(v)); ctx.lineTo(lW+PW,toY(v)); ctx.stroke();
  });
  ctx.setLineDash([]);

  // Points
  const COL={U:'#27ae60',D:'#e74c3c',N:'#aaaaaa'};
  const ALPHA={U:'dd',D:'dd',N:'55'};
  pts.forEach(p=>{
    if(p.a==null||p.m==null)return;
    const x=toX(p.a), y=toY(p.m);
    ctx.fillStyle=(COL[p.s]||'#aaa')+(ALPHA[p.s]||'55');
    ctx.beginPath(); ctx.arc(x,y,p.s==='N'?2:3.5,0,2*Math.PI); ctx.fill();
  });

  // Labels axes
  ctx.fillStyle='#8b949e'; ctx.font='bold 9px sans-serif'; ctx.textAlign='center';
  ctx.fillText('A = ('+d.label_a+' + '+d.label_b+') / 2  [log\u2082 LFQ]', lW+PW/2, tH+PH+30);
  ctx.fillStyle='#8b949e'; ctx.save(); ctx.translate(13, tH+PH/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('M = log\u2082 Fold Change', 0, 0); ctx.restore();

  // Légende
  [['UP','#27ae60'],['DOWN','#e74c3c'],['NS','#8b949e']].forEach(([l,c],i)=>{
    ctx.fillStyle=c; ctx.beginPath();
    ctx.arc(lW+PW-80+i*55, tH+10, 4, 0, 2*Math.PI); ctx.fill();
    ctx.fillStyle='#c9d1d9'; ctx.font='8px sans-serif'; ctx.textAlign='left';
    ctx.fillText(l, lW+PW-74+i*55, tH+13);
  });

  // Tooltip interactif
  canvas._maPts=pts; canvas._maToX=toX; canvas._maToY=toY;
  canvas._maLW=lW; canvas._maTH=tH;
  canvas.onmousemove=function(e){
    const rect=canvas.getBoundingClientRect();
    const mx=(e.clientX-rect.left)*(W/rect.width);
    const my=(e.clientY-rect.top)*(H/rect.height);
    let best=null, bestD=1e9;
    pts.forEach(p=>{
      if(p.a==null)return;
      const dx=toX(p.a)-mx, dy=toY(p.m)-my;
      const dist=Math.sqrt(dx*dx+dy*dy);
      if(dist<bestD){bestD=dist;best=p;}
    });
    const tip=document.getElementById('maTip');
    if(tip&&best&&bestD<18){
      tip.style.display='block';
      tip.style.left=(e.clientX+12)+'px'; tip.style.top=(e.clientY-28)+'px';
      tip.innerHTML='<b>'+(best.g||best.n||'—')+'</b><br>'
        +'A='+best.a+'  M='+best.m+'<br>'
        +'p='+( best.p!=null?best.p.toExponential(2):'—')
        +'  padj='+(best.pa!=null?best.pa.toExponential(2):'—');
    } else if(tip){ tip.style.display='none'; }
  };
  canvas.onmouseleave=function(){
    const tip=document.getElementById('maTip');
    if(tip)tip.style.display='none';
  };
}
function initMAPlot(){
  const sel=document.getElementById('maContrast');
  if(!sel||sel.options.length)return;
  Object.keys(MA_DATA).forEach(c=>{
    const o=document.createElement('option'); o.value=c;
    o.textContent=c.replace(/_vs_/,' vs ').replace(/_/g,' '); sel.appendChild(o);
  });
  sel.onchange=drawMAPlot;
  drawMAPlot();
}\
"""

# drawRepScatter : matrice de mini-panels scatter inter-réplicats par condition
#   • Un panel par paire de réplicats (grille n×n)
#   • r de Pearson affiché dans chaque panel
#   • Couleur des points par condition (palette TIFF)
DRAW_REP_SCATTER_JS = """\
function drawRepScatter(){
  if(!REP_SCATTER_DATA||!Object.keys(REP_SCATTER_DATA).length)return;
  const sel=document.getElementById('repScatterCond');
  const cond=sel?sel.value:Object.keys(REP_SCATTER_DATA)[0];
  const d=REP_SCATTER_DATA[cond];
  if(!d||!d.pairs||!d.pairs.length)return;

  const reps=d.reps;
  const n=reps.length;
  // Matrice triangulaire supérieure : n*(n-1)/2 panels
  // On dispose les panels en grille (n-1) colonnes x (n-1) lignes (diagonale = nom réplique)
  const cellSize=110, pad=4, labelH=20;
  const W=n*cellSize+pad, H=n*cellSize+pad+labelH;
  const canvas=document.getElementById('repScatterCanvas');
  if(!canvas)return;
  const dpr=devicePixelRatio||1;
  canvas.width=W*dpr; canvas.height=H*dpr;
  canvas.style.width=W+'px'; canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.clearRect(0,0,W,H);

  // Titre
  ctx.fillStyle='#c9d1d9'; ctx.font='bold 10px sans-serif'; ctx.textAlign='center';
  ctx.fillText('Corrélations inter-réplicats — '+cond, W/2, 13);

  const col=COND_COLORS[cond]||'#4a90d9';

  // Plage globale pour axes cohérents
  let gMin=Infinity, gMax=-Infinity;
  d.pairs.forEach(p=>p.pts.forEach(pt=>{
    gMin=Math.min(gMin,pt.x,pt.y); gMax=Math.max(gMax,pt.x,pt.y);
  }));
  const gRange=gMax-gMin||1;

  // Lookup paire -> données
  const pairMap={};
  d.pairs.forEach(p=>{ pairMap[p.r1+'|'+p.r2]=p; pairMap[p.r2+'|'+p.r1]=p; });

  // Dessin matrice n×n
  for(let ri=0;ri<n;ri++){
    for(let ci=0;ci<n;ci++){
      const xOff=ci*cellSize+pad, yOff=ri*cellSize+pad+labelH;
      const cs=cellSize-pad;

      if(ri===ci){
        // Diagonale : nom du réplique
        ctx.fillStyle='#21262d'; ctx.fillRect(xOff,yOff,cs,cs);
        ctx.fillStyle='#c9d1d9'; ctx.font='bold 9px sans-serif'; ctx.textAlign='center';
        const label=reps[ri].replace(/^.*_(\\d+)$/,'Rep $1');
        ctx.fillText(label, xOff+cs/2, yOff+cs/2-4);
        ctx.font='8px sans-serif'; ctx.fillStyle='#8b949e';
        ctx.fillText(reps[ri].split('_').slice(0,-1).join(' ').replace(/Polyculture\\./,'Poly').replace('Monoculture','Mono'),
          xOff+cs/2, yOff+cs/2+8);
      } else {
        // Panel scatter
        ctx.fillStyle='#0d1117'; ctx.fillRect(xOff,yOff,cs,cs);
        ctx.strokeStyle='#21262d'; ctx.lineWidth=0.5;
        ctx.strokeRect(xOff,yOff,cs,cs);

        const pair=pairMap[reps[ci]+'|'+reps[ri]];
        if(!pair){continue;}

        function toXp(v){ return xOff+4+((v-gMin)/gRange)*(cs-8); }
        function toYp(v){ return yOff+cs-4-((v-gMin)/gRange)*(cs-8); }

        // Diagonale de référence (ligne y=x)
        ctx.strokeStyle='#30363d'; ctx.lineWidth=0.8;
        ctx.beginPath(); ctx.moveTo(toXp(gMin),toYp(gMin));
        ctx.lineTo(toXp(gMax),toYp(gMax)); ctx.stroke();

        // Points (garde anti-gel : on plafonne le nb de points dessinés par
        // cellule, au cas où les données ne seraient pas sous-échantillonnées)
        ctx.fillStyle=col+'99';
        var _pp=pair.pts;
        var _step=_pp.length>200?Math.ceil(_pp.length/200):1;
        for(var _k=0;_k<_pp.length;_k+=_step){
          var pt=_pp[_k];
          ctx.beginPath(); ctx.arc(toXp(pt.x),toYp(pt.y),1.8,0,2*Math.PI); ctx.fill();
        }

        // r de Pearson
        const rVal=pair.r;
        const rColor=rVal>0.99?'#27ae60':rVal>0.95?'#f39c12':'#e74c3c';
        ctx.fillStyle=rColor; ctx.font='bold 9px sans-serif'; ctx.textAlign='center';
        ctx.fillText('r='+rVal.toFixed(4), xOff+cs/2, yOff+cs-4);
        ctx.fillStyle='#6e7681'; ctx.font='7px sans-serif';
        ctx.fillText('n='+pair.n, xOff+cs/2, yOff+10);
      }
    }
  }
}
function initRepScatter(){
  const sel=document.getElementById('repScatterCond');
  if(!sel||sel.options.length)return;
  Object.keys(REP_SCATTER_DATA).forEach(c=>{
    const o=document.createElement('option'); o.value=c; o.textContent=c; sel.appendChild(o);
  });
  sel.onchange=drawRepScatter;
  drawRepScatter();
}\
"""


# ──────────────────────────────────────────────────────────────────────────────
# Patch JS — UMAP interactif
# ──────────────────────────────────────────────────────────────────────────────

DRAW_UMAP_JS = """\
function drawUMAP(){
  if(!UMAP_DATA||!UMAP_DATA.points||!UMAP_DATA.points.length){
    const cv=document.getElementById('umapCanvas');
    if(cv){
      cv.width=400;cv.height=80;
      const ctx=cv.getContext('2d');
      ctx.fillStyle='#fff';ctx.fillRect(0,0,400,80);
      ctx.fillStyle='#888';ctx.font='13px sans-serif';ctx.textAlign='center';
      ctx.fillText('Donnees UMAP non disponibles',200,45);
    }
    return;
  }
  const pts=UMAP_DATA.points;
  const conds=UMAP_DATA.conditions;
  // Palette HSL : teintes regulierement espacees
  const UMAP_COLORS={};
  conds.forEach((c,i)=>{
    const h=Math.round(i*360/conds.length);
    UMAP_COLORS[c]='hsl('+h+',65%,52%)';
  });
  const canvas=document.getElementById('umapCanvas');
  if(!canvas)return;
  const dpr=devicePixelRatio||1;
  const lW=50,rW=200,tH=36,bH=36,PW=480,PH=420;
  const W=lW+PW+rW,H=tH+PH+bH;
  canvas.width=W*dpr;canvas.height=H*dpr;
  canvas.style.width=W+'px';canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.clearRect(0,0,W,H);
  // Plages avec marges 7%
  const u1s=pts.map(p=>p.u1),u2s=pts.map(p=>p.u2);
  const u1Min=Math.min(...u1s),u1Max=Math.max(...u1s);
  const u2Min=Math.min(...u2s),u2Max=Math.max(...u2s);
  const u1R=(u1Max-u1Min)||1,u2R=(u2Max-u2Min)||1;
  const mg=0.07;
  function toX(v){return lW+(v-u1Min+u1R*mg)/(u1R*(1+2*mg))*PW;}
  function toY(v){return tH+PH-(v-u2Min+u2R*mg)/(u2R*(1+2*mg))*PH;}
  // Grille legere
  ctx.strokeStyle='#21262d';ctx.lineWidth=1;
  [0.25,0.5,0.75].forEach(f=>{
    const xv=u1Min+f*u1R,yv=u2Min+f*u2R;
    ctx.beginPath();ctx.moveTo(toX(xv),tH);ctx.lineTo(toX(xv),tH+PH);ctx.stroke();
    ctx.beginPath();ctx.moveTo(lW,toY(yv));ctx.lineTo(lW+PW,toY(yv));ctx.stroke();
  });
  // Axes zero
  ctx.strokeStyle='#30363d';ctx.lineWidth=1.5;
  if(u1Min<0&&u1Max>0){ctx.beginPath();ctx.moveTo(toX(0),tH);ctx.lineTo(toX(0),tH+PH);ctx.stroke();}
  if(u2Min<0&&u2Max>0){ctx.beginPath();ctx.moveTo(lW,toY(0));ctx.lineTo(lW+PW,toY(0));ctx.stroke();}
  // Labels axes
  ctx.fillStyle='#8b949e';ctx.font='9px sans-serif';ctx.textAlign='center';
  ctx.fillText('UMAP1',lW+PW/2,tH+PH+22);
  ctx.save();ctx.translate(14,tH+PH/2);ctx.rotate(-Math.PI/2);
  ctx.fillStyle='#8b949e';ctx.fillText('UMAP2',0,0);ctx.restore();
  // Titre
  ctx.fillStyle='#c9d1d9';ctx.font='bold 11px sans-serif';ctx.textAlign='center';
  ctx.fillText('UMAP - '+pts.length+' echantillons, '+conds.length+' conditions',lW+PW/2,20);
  // Ellipses de confiance (centroide +/- 1 sigma)
  conds.forEach(cond=>{
    const cp=pts.filter(p=>p.cond===cond);
    if(cp.length<3)return;
    const mu1=cp.reduce((s,p)=>s+p.u1,0)/cp.length;
    const mu2=cp.reduce((s,p)=>s+p.u2,0)/cp.length;
    const s1=Math.sqrt(cp.reduce((s,p)=>s+(p.u1-mu1)**2,0)/cp.length)||0.01;
    const s2=Math.sqrt(cp.reduce((s,p)=>s+(p.u2-mu2)**2,0)/cp.length)||0.01;
    const col=UMAP_COLORS[cond]||'#888';
    ctx.save();
    ctx.strokeStyle=col;ctx.lineWidth=1.2;ctx.setLineDash([4,3]);ctx.globalAlpha=0.5;
    ctx.beginPath();
    ctx.ellipse(toX(mu1),toY(mu2),
      s1/(u1R*(1+2*mg))*PW,s2/(u2R*(1+2*mg))*PH,0,0,2*Math.PI);
    ctx.stroke();
    ctx.setLineDash([]);ctx.restore();
  });
  // Points
  pts.forEach(p=>{
    const col=UMAP_COLORS[p.cond]||'#888';
    const x=toX(p.u1),y=toY(p.u2);
    ctx.beginPath();ctx.arc(x,y,6,0,2*Math.PI);
    ctx.fillStyle=col+'dd';ctx.fill();
    ctx.strokeStyle=col;ctx.lineWidth=1;ctx.stroke();
    ctx.fillStyle='#8b949e';ctx.font='7px sans-serif';ctx.textAlign='center';
    ctx.fillText(p.sample.replace(/^.*_(\\d+)$/,'$1'),x,y-8);
  });
  // Legende
  const lgX=lW+PW+12,lgY=tH+10;
  ctx.fillStyle='#8b949e';ctx.font='bold 9px sans-serif';ctx.textAlign='left';
  ctx.fillText('Conditions',lgX,lgY);
  conds.forEach((cond,i)=>{
    const col=UMAP_COLORS[cond]||'#888';
    const y=lgY+18+i*22;
    ctx.fillStyle=col;ctx.beginPath();ctx.arc(lgX+7,y-4,6,0,2*Math.PI);ctx.fill();
    ctx.fillStyle='#c9d1d9';ctx.font='9px sans-serif';ctx.textAlign='left';
    ctx.fillText(cond.slice(0,24),lgX+18,y);
  });
  // Tooltip
  canvas._umapPts=pts;canvas._umapToX=toX;canvas._umapToY=toY;
  const tip=document.getElementById('umapTip');
  canvas.onmousemove=function(e){
    const rect=canvas.getBoundingClientRect();
    const mx=(e.clientX-rect.left)*(W/rect.width);
    const my=(e.clientY-rect.top)*(H/rect.height);
    let best=null,bestD=1e9;
    pts.forEach(p=>{
      const dx=toX(p.u1)-mx,dy=toY(p.u2)-my;
      const d=Math.sqrt(dx*dx+dy*dy);
      if(d<bestD){bestD=d;best=p;}
    });
    if(tip&&best&&bestD<14){
      tip.style.display='block';
      tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY-28)+'px';
      tip.innerHTML='<b>'+best.sample+'</b><br><span style="color:#888">'+best.cond+'</span>'+
        '<br><span style="font-size:9px">UMAP1='+best.u1.toFixed(3)+'  UMAP2='+best.u2.toFixed(3)+'</span>';
    } else if(tip){tip.style.display='none';}
  };
  canvas.onmouseleave=function(){if(tip)tip.style.display='none';};
}\
"""


def patch_umap_js(html):
    """Injecte drawUMAP() dans le HTML final avant </script>.
    Non-bloquant si la balise est absente.
    """
    pos = html.rfind('</script>')
    if pos == -1:
        print("  [warn] patch_umap_js : balise </script> introuvable.")
        return html
    html = html[:pos] + "\n\n" + DRAW_UMAP_JS + "\n" + html[pos:]
    print("  patch_umap_js : drawUMAP() injectee (%d chars)." % len(DRAW_UMAP_JS))
    return html


def patch_qc_figures(html):
    """Injecte drawViolinPlot(), drawMissingValues(), drawMAPlot() et drawRepScatter()
    dans le HTML final. Insertion juste avant la dernière balise </script>.
    Non-bloquant si la balise n'est pas trouvée.
    """
    new_funcs = (
        "\n\n" + DRAW_VIOLIN_JS +
        "\n\n" + DRAW_MISSING_JS +
        "\n\n" + DRAW_MA_JS +
        "\n\n" + DRAW_REP_SCATTER_JS +
        "\n"
    )
    pos = html.rfind('</script>')
    if pos == -1:
        print("  [warn] patch_qc_figures : balise </script> introuvable, figures QC non injectées.")
        return html
    html = html[:pos] + new_funcs + html[pos:]
    # Corriger l'affichage du r de Pearson (obj.r peut être NaN depuis SCATTER_DATA calculé en JS)
    html = html.replace(
        "ctx.fillText('r = '+(obj.r||'').toFixed(4),",
        "ctx.fillText('r = '+((typeof obj.r==='number'&&!isNaN(obj.r))?obj.r.toFixed(4):'n/a'),"
    )
    # Fix overlap titre/comptages dans drawVolcano — déplacer UP/DOWN à gauche
    html = html.replace(
        "ctx.fillStyle='rgba(63,185,80,.9)';ctx.fillText('↑ '+su,pad.l+pw-52,pad.t+16);",
        "ctx.fillStyle='rgba(63,185,80,.9)';ctx.fillText('↑ '+su,pad.l+6,pad.t+16);"
    ).replace(
        "ctx.fillStyle='rgba(248,81,73,.9)';ctx.fillText('↓ '+sd,pad.l+pw-52,pad.t+32);",
        "ctx.fillStyle='rgba(248,81,73,.9)';ctx.fillText('↓ '+sd,pad.l+6,pad.t+32);"
    )
    # Ajouter filtre LFC min dans la table DEP
    html = html.replace(
        '<div class="radio-grp">\n      <input type="radio" name="ptypeProt" id="pp_pval" value="p" checked>',
        '<span class="clbl">|LFC| ≥</span>'
        '<input type="range" id="protLFC" min="0" max="3" step="0.1" value="0.26" style="width:80px;" oninput="document.getElementById(\'protLFCVal\').textContent=parseFloat(this.value).toFixed(2);renderProtTable()">'
        '<span class="sv" id="protLFCVal">0.26</span>'
        '<div class="radio-grp">\n      <input type="radio" name="ptypeProt" id="pp_pval" value="p" checked>'
    )
    # Adapter renderProtTable pour utiliser le filtre LFC
    html = html.replace(
        "const fc=0.26;\n  let pts=(VOLCANO[c]||[]).filter(p=>{",
        "const fc=parseFloat((document.getElementById('protLFC')||{value:'0.26'}).value)||0.26;\n  let pts=(VOLCANO[c]||[]).filter(p=>{"
    )
    # Ajouter tri par colonne dans la table DEP
    html = html.replace(
        '<table class="pt"><thead><tr>\n      <th>Protéine</th><th>Gène</th><th>Description</th><th>Statut</th>\n      <th>LFC</th><th>p.val</th><th>p.adj</th><th>Pi score</th><th>Robustesse</th>\n    </tr></thead>',
        '<table class="pt"><thead><tr>'
        '<th>Protéine</th><th>Gène</th><th>Description</th><th>Statut</th>'
        '<th style="cursor:pointer;user-select:none;" onclick="sortProtTable(\'lfc\')" id="th_lfc">LFC ↕</th>'
        '<th style="cursor:pointer;user-select:none;" onclick="sortProtTable(\'p\')"   id="th_p">p.val ↕</th>'
        '<th style="cursor:pointer;user-select:none;" onclick="sortProtTable(\'padj\')" id="th_padj">p.adj ↕</th>'
        '<th style="cursor:pointer;user-select:none;" onclick="sortProtTable(\'pi\')"  id="th_pi">Pi score ↕</th>'
        '<th style="cursor:pointer;user-select:none;" onclick="sortProtTable(\'rob\')" id="th_rob">Robustesse ↕</th>'
        '</tr></thead>'
    )
    # Injecter la logique de tri dans renderProtTable
    html = html.replace(
        ".sort((a,b)=>parseFloat(a.p)-parseFloat(b.p));",
        ".sort((a,b)=>{\n"
        "    var col=window._protSortCol||'p', asc=window._protSortAsc!==false;\n"
        "    var va,vb;\n"
        "    if(col==='lfc'){va=Math.abs(parseFloat(a.d));vb=Math.abs(parseFloat(b.d));}\n"
        "    else if(col==='p'){va=parseFloat(a.p);vb=parseFloat(b.p);}\n"
        "    else if(col==='padj'){va=parseFloat(a.padj);vb=parseFloat(b.padj);}\n"
        "    else if(col==='pi'){va=parseFloat(a.pi)||0;vb=parseFloat(b.pi)||0;}\n"
        "    else if(col==='rob'){va=parseInt(a.r)||0;vb=parseInt(b.r)||0;}\n"
        "    else{va=parseFloat(a.p);vb=parseFloat(b.p);}\n"
        "    return asc?(va-vb):(vb-va);\n"
        "  });"
    )
    # Injecter la fonction sortProtTable
    html = html.replace(
        "function renderProtTable(){",
        "window._protSortCol='p'; window._protSortAsc=true;\n"
        "window.sortProtTable=function(col){\n"
        "  if(window._protSortCol===col){ window._protSortAsc=!window._protSortAsc; }\n"
        "  else { window._protSortCol=col; window._protSortAsc=(col==='p'||col==='padj'); }\n"
        "  // Mettre à jour les flèches dans les en-têtes\n"
        "  ['lfc','p','padj','pi','rob'].forEach(function(c){\n"
        "    var th=document.getElementById('th_'+c); if(!th)return;\n"
        "    var name={'lfc':'LFC','p':'p.val','padj':'p.adj','pi':'Pi score','rob':'Robustesse'}[c];\n"
        "    th.textContent=name+(window._protSortCol===c?(window._protSortAsc?' ↑':' ↓'):' ↕');\n"
        "  });\n"
        "  renderProtTable();\n"
        "};\n"
        "function renderProtTable(){"
    )
    print("  patch_qc_figures : violin + missing + MA plot + scatter réplicats injectés (%d chars)" % len(new_funcs))
    return html


# ──────────────────────────────────────────────────────────────────────────────
# Patch JS — UpSet plot interactif (clic sur barre → tableau protéines)
# ──────────────────────────────────────────────────────────────────────────────

# Stratégie :
#   • Réécriture complète de initUpset() dans le template
#   • Chaque barre mémorise sa signature exacte (ensemble de contrastes actifs)
#   • Clic → filtre INTERSECTIONS sur la signature → peuple #upsetDetail
#   • Barre sélectionnée surlignée, redraw partiel via une fonction dédiée
#   • Hover tooltip conservé sur les barres non sélectionnées
DRAW_UPSET_JS = """\
/* ── UPSET INTERACTIF (patch) ── */
(function(){
// Remplace initUpset dès que le DOM est prêt ; si la fonction d'origine a déjà
// été appelée (upsetDrawn=true) on repart d'un état propre.
const _origInit = typeof initUpset === 'function' ? initUpset : null;

window.initUpset = function(){
  // reset au cas où le template aurait déjà exécuté l'ancienne version
  const canvas = document.getElementById('upsetCanvas');
  if(!canvas) return;

  const U = UPSET_DATA;
  const SETS = U.sets, INTERS = U.intersections,
        SET_SIZES = U.set_sizes, FULL = U.full_labels;
  const nS = SETS.length, nI = INTERS.length;

  // Guard : UpSet nécessite ≥2 contrastes
  if(nS < 2 || !nI){
    const dpr = devicePixelRatio||1;
    canvas.width=500*dpr; canvas.height=80*dpr;
    canvas.style.width='500px'; canvas.style.height='80px';
    const ctx=canvas.getContext('2d'); ctx.scale(dpr,dpr);
    ctx.clearRect(0,0,500,80);
    ctx.fillStyle='#8b949e'; ctx.font='13px sans-serif'; ctx.textAlign='center';
    ctx.fillText(nS<2
      ? 'UpSet nécessite au moins 2 contrastes'
      : 'Aucune intersection à afficher', 250, 45);
    return;
  }
  const BAR_H=180, DOT_H=22, COL_W=46, LBL_W=52, SET_W=80, PAD_T=24, PAD_B=10;
  const MATRIX_H = nS*DOT_H + PAD_B;
  const TOTAL_H  = PAD_T + BAR_H + MATRIX_H + 42;
  const TOTAL_W  = SET_W + LBL_W + nI*COL_W + 20;
  const dpr = devicePixelRatio || 1;
  canvas.width  = TOTAL_W*dpr; canvas.height = TOTAL_H*dpr;
  canvas.style.width  = TOTAL_W + 'px';
  canvas.style.height = TOTAL_H + 'px';

  const ctx = canvas.getContext('2d');
  ctx.scale(dpr, dpr);

  const maxCnt = Math.max(...INTERS.map(x=>x.count), 1);
  const maxSz  = Math.max(...Object.values(SET_SIZES), 1);
  const colX   = i => SET_W + LBL_W + i*COL_W + COL_W/2;
  const rowY   = s => PAD_T + BAR_H + s*DOT_H + DOT_H/2;
  const barBot = PAD_T + BAR_H;

  // ── Lookup signature → protéines (depuis INTERSECTIONS) ──
  // INTERSECTIONS utilise les clés c1..cN dans l'ordre de CONTRASTS
  // UPSET_DATA.full_labels mappe label court -> nom contraste complet
  // On reconstruit : pour chaque protéine, un Set de labels courts actifs
  const shortLabels = SETS;  // ["M/P1","M/P2",...]

  // Pré-indexer : signature (labels actifs triés, join '|') -> [protéines]
  const sigIndex = {};
  INTERSECTIONS.forEach(row => {
    const active = shortLabels.filter((_,si) => {
      // La colonne dans INTERSECTIONS est 'c1','c2',...
      const key = 'c' + (si+1);
      return row[key] === 'X';
    });
    const sig = active.slice().sort().join('|');
    if(!sigIndex[sig]) sigIndex[sig] = [];
    sigIndex[sig].push(row);
  });

  let selectedSig = null;  // signature actuellement sélectionnée

  // ── Dessin ──
  function draw(highlightIdx) {
    ctx.fillStyle = '#161b22'; ctx.fillRect(0, 0, TOTAL_W, TOTAL_H);

    // Barres d'intersection
    INTERS.forEach((inter, i) => {
      const x  = colX(i);
      const bh = Math.round(inter.count / maxCnt * (BAR_H - 8));
      const by = barBot - bh;
      const isHL = (i === highlightIdx);

      // Fond surligné si sélectionné
      if(isHL){
        ctx.fillStyle = 'rgba(56,139,253,0.12)';
        ctx.fillRect(x - COL_W/2 + 1, PAD_T, COL_W - 2, BAR_H + MATRIX_H);
      }

      ctx.fillStyle = isHL ? '#60a5fa' : 'rgba(56,139,253,.75)';
      ctx.beginPath();
      ctx.roundRect(x - COL_W/2 + 6, by, COL_W - 12, bh, 3);
      ctx.fill();

      ctx.fillStyle = isHL ? '#ffffff' : '#c9d1d9';
      ctx.font = 'bold 9px sans-serif'; ctx.textAlign = 'center';
      ctx.fillText(inter.count, x, by - 3);
    });

    // Graduations axe Y
    ctx.fillStyle = '#6e7681'; ctx.font = '9px sans-serif'; ctx.textAlign = 'right';
    [0,.25,.5,.75,1].forEach(f => {
      const val = Math.round(f * maxCnt);
      const y   = barBot - Math.round(f * (BAR_H - 8));
      ctx.fillText(val, SET_W + LBL_W - 4, y + 3);
      ctx.strokeStyle = 'rgba(240,246,252,.04)'; ctx.lineWidth = .5;
      ctx.beginPath(); ctx.moveTo(SET_W+LBL_W, y); ctx.lineTo(TOTAL_W-10, y); ctx.stroke();
    });

    // Bandes alternées matrice
    SETS.forEach((_, si) => {
      ctx.fillStyle = si%2===0 ? 'rgba(255,255,255,.02)' : 'transparent';
      ctx.fillRect(SET_W+LBL_W, rowY(si)-DOT_H/2, nI*COL_W, DOT_H);
    });

    // Labels sets (axe gauche)
    SETS.forEach((s, si) => {
      ctx.fillStyle = '#c9d1d9'; ctx.font = '11px sans-serif'; ctx.textAlign = 'right';
      ctx.fillText(s, SET_W+LBL_W-8, rowY(si)+4);
    });

    // Barres set sizes (gauche)
    SETS.forEach((s, si) => {
      const y  = rowY(si);
      const bw = Math.round((SET_SIZES[s]||0) / maxSz * (SET_W-8));
      ctx.fillStyle = 'rgba(188,140,255,.6)';
      ctx.beginPath();
      ctx.roundRect(SET_W-4-bw, y-DOT_H/2+3, bw, DOT_H-6, 2);
      ctx.fill();
      ctx.fillStyle = '#8b949e'; ctx.font = '9px sans-serif'; ctx.textAlign = 'right';
      ctx.fillText(SET_SIZES[s]||0, SET_W-8-bw, y+4);
    });

    // Points + lignes de connexion matrice
    INTERS.forEach((inter, i) => {
      const x = colX(i);
      const isHL = (i === highlightIdx);
      const actIdx = inter.sets
        .map(s => SETS.indexOf(s))
        .filter(idx => idx >= 0)
        .sort((a,b) => a-b);

      if(actIdx.length > 1){
        ctx.strokeStyle = isHL ? '#60a5fa' : 'rgba(56,139,253,.85)';
        ctx.lineWidth = isHL ? 2.5 : 2;
        ctx.beginPath();
        ctx.moveTo(x, rowY(actIdx[0]));
        ctx.lineTo(x, rowY(actIdx[actIdx.length-1]));
        ctx.stroke();
      }
      SETS.forEach((s, si) => {
        const isA = inter.sets.includes(s);
        ctx.beginPath(); ctx.arc(x, rowY(si), isA ? 6 : 4, 0, Math.PI*2);
        ctx.fillStyle = isA ? (isHL ? '#60a5fa' : '#388bfd') : '#2d333b';
        ctx.fill();
        if(isA){
          ctx.strokeStyle = isHL ? 'rgba(96,165,250,.5)' : 'rgba(56,139,253,.35)';
          ctx.lineWidth = 1.5; ctx.stroke();
        }
      });
    });

    // Instruction clic (première fois)
    if(selectedSig === null){
      ctx.fillStyle = 'rgba(139,148,158,0.7)'; ctx.font = '9px sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText('↑ Cliquez sur une barre pour voir les protéines',
        SET_W + LBL_W + nI*COL_W/2, TOTAL_H - 6);
    }
  }

  draw(null);

  // ── Recherche de l'index de barre le plus proche du pointeur ──
  function barAt(mx) {
    for(let i = 0; i < nI; i++){
      if(Math.abs(mx - colX(i)) < COL_W/2) return i;
    }
    return -1;
  }

  // ── Tooltip hover ──
  const tt = document.getElementById('tooltip');
  canvas.addEventListener('mousemove', e => {
    const rect = canvas.getBoundingClientRect();
    const mx = (e.clientX - rect.left) * (TOTAL_W / rect.width);
    const my = (e.clientY - rect.top)  * (TOTAL_H / rect.height);
    const idx = (my < barBot) ? barAt(mx) : -1;
    canvas.style.cursor = idx >= 0 ? 'pointer' : 'default';
    if(idx >= 0 && tt){
      const inter = INTERS[idx];
      tt.innerHTML = '<div class="tg">' + inter.count + ' protéines</div>' +
        '<div class="tl" style="margin-top:3px;">' +
        inter.sets.map(s => FULL[s]||s).join('<br>') + '</div>' +
        '<div style="font-size:9px;margin-top:4px;color:#8b949e;">Cliquez pour filtrer</div>';
      tt.style.display = 'block';
      tt.style.left = (e.clientX + 14) + 'px';
      tt.style.top  = (e.clientY - 10) + 'px';
    } else if(tt) {
      tt.style.display = 'none';
    }
  });
  canvas.addEventListener('mouseleave', () => {
    if(tt) tt.style.display = 'none';
    canvas.style.cursor = 'default';
  });

  // ── Clic : filtre le tableau de protéines ──
  canvas.addEventListener('click', e => {
    const rect = canvas.getBoundingClientRect();
    const mx = (e.clientX - rect.left) * (TOTAL_W / rect.width);
    const my = (e.clientY - rect.top)  * (TOTAL_H / rect.height);
    if(my >= barBot) return;  // clic dans la matrice de points : ignoré
    const idx = barAt(mx);
    if(idx < 0) return;

    const inter = INTERS[idx];
    const sig   = inter.sets.slice().sort().join('|');

    if(selectedSig === sig){
      // Deuxième clic sur la même barre : désélection
      selectedSig = null;
      draw(null);
      renderUpsetDetail(null);
      return;
    }
    selectedSig = sig;
    draw(idx);
    renderUpsetDetail(inter, sigIndex[sig] || []);
  });
};

// ── Panneau de détail (tableau protéines de l'intersection) ──
window.renderUpsetDetail = function(inter, rows){
  let panel = document.getElementById('upsetDetail');
  if(!panel){
    // Créer le panneau s'il n'existe pas encore dans le template
    panel = document.createElement('div');
    panel.id = 'upsetDetail';
    panel.style.cssText =
      'margin-top:14px;padding:10px 14px;background:#0d1117;border:1px solid #30363d;' +
      'border-radius:6px;font-size:11px;';
    const canvas = document.getElementById('upsetCanvas');
    if(canvas && canvas.parentNode) canvas.parentNode.appendChild(panel);
  }

  if(!inter || !rows || !rows.length){
    panel.innerHTML = '';
    panel.style.display = 'none';
    return;
  }

  const U = UPSET_DATA;
  const FULL = U.full_labels;
  const title = inter.sets.map(s => FULL[s]||s)
    .map(s => s.replace(/_vs_/g,' vs ').replace(/_/g,' '))
    .join(' ∩ ');

  // En-tête
  let html = '<div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px;">' +
    '<span style="font-weight:700;color:#c9d1d9;">' + inter.count + ' protéines — ' + title + '</span>' +
    '<div style="display:flex;gap:6px;">' +
    '<button onclick="exportUpsetCSV()" style="background:none;border:1px solid #388bfd;color:#388bfd;border-radius:4px;padding:2px 8px;cursor:pointer;font-size:10px;">⬇ CSV</button>' +
    '<button onclick="renderUpsetDetail(null)" style="background:none;border:1px solid #30363d;color:#8b949e;border-radius:4px;padding:2px 8px;cursor:pointer;font-size:10px;">✕ Fermer</button>' +
    '</div></div>';

  // Champ de recherche dans le panneau
  html += '<input id="upsetDetailSearch" type="text" placeholder="Filtrer gène / description…" ' +
    'oninput="filterUpsetDetail()" ' +
    'style="width:100%;box-sizing:border-box;background:#161b22;border:1px solid #30363d;' +
    'color:#c9d1d9;border-radius:4px;padding:4px 8px;font-size:10px;margin-bottom:8px;">';

  // Tableau
  html += '<div style="overflow-y:auto;max-height:320px;">' +
    '<table style="width:100%;border-collapse:collapse;" id="upsetDetailTable">' +
    '<thead><tr style="position:sticky;top:0;background:#161b22;">' +
    '<th style="text-align:left;padding:4px 6px;color:#8b949e;font-weight:600;border-bottom:1px solid #21262d;">Gène</th>' +
    '<th style="text-align:left;padding:4px 6px;color:#8b949e;font-weight:600;border-bottom:1px solid #21262d;">Accession</th>' +
    '<th style="text-align:left;padding:4px 6px;color:#8b949e;font-weight:600;border-bottom:1px solid #21262d;">Description</th>' +
    '<th style="text-align:center;padding:4px 6px;color:#8b949e;font-weight:600;border-bottom:1px solid #21262d;">Nb contrastes</th>' +
    '</tr></thead><tbody id="upsetDetailBody"></tbody></table></div>';

  panel.innerHTML = html;
  panel.style.display = 'block';

  // Stocker les lignes et la signature pour le filtre et l'export
  panel._rows = rows;
  panel._sig  = inter.sets.slice().sort().join('_vs_');
  _renderUpsetRows(rows);
};

function _renderUpsetRows(rows){
  const tbody = document.getElementById('upsetDetailBody');
  if(!tbody) return;
  tbody.innerHTML = '';
  rows.slice(0, 500).forEach(r => {
    const tr = document.createElement('tr');
    tr.style.borderBottom = '1px solid #21262d';
    tr.innerHTML =
      '<td style="padding:3px 6px;font-weight:600;color:#c9d1d9;">' + (r.g||'—') + '</td>' +
      '<td style="padding:3px 6px;font-family:monospace;font-size:10px;color:#8b949e;">' + (r.pg||'') + '</td>' +
      '<td style="padding:3px 6px;color:#8b949e;font-size:10px;" title="' + (r.desc||'') + '">' +
        ((r.desc||'').slice(0,55)||(r.desc?'…':'—')) + '</td>' +
      '<td style="text-align:center;padding:3px 6px;">' +
        '<span style="font-weight:700;color:' + (r.nb>=4?'#3fb950':r.nb>=2?'#d29922':'#8b949e') + '">' +
        r.nb + '</span></td>';
    tbody.appendChild(tr);
  });
  const info = document.getElementById('upsetDetailInfo');
  if(info) info.textContent = rows.length + ' protéines' + (rows.length>500?' (500 affichées)':'');
}

window.filterUpsetDetail = function(){
  const panel = document.getElementById('upsetDetail');
  if(!panel || !panel._rows) return;
  const q = (document.getElementById('upsetDetailSearch').value||'').toLowerCase().trim();
  const filtered = q
    ? panel._rows.filter(r =>
        (r.g||'').toLowerCase().includes(q) ||
        (r.pg||'').toLowerCase().includes(q) ||
        (r.desc||'').toLowerCase().includes(q))
    : panel._rows;
  _renderUpsetRows(filtered);
};

window.exportUpsetCSV = function(){
  const panel = document.getElementById('upsetDetail');
  if(!panel || !panel._rows || !panel._rows.length) return;
  const q = (document.getElementById('upsetDetailSearch') || {}).value || '';
  const qn = q.toLowerCase().trim();
  const rows = qn
    ? panel._rows.filter(r =>
        (r.g||'').toLowerCase().includes(qn) ||
        (r.pg||'').toLowerCase().includes(qn) ||
        (r.desc||'').toLowerCase().includes(qn))
    : panel._rows;
  // En-tête CSV
  const lines = ['Accession,Gene,Description,Nb_Contrastes'];
  rows.forEach(r => {
    const esc = v => '"' + String(v||'').replace(/"/g,'""') + '"';
    lines.push([esc(r.pg), esc(r.g), esc(r.desc), r.nb].join(','));
  });
  const blob = new Blob([lines.join('\\n')], {type:'text/csv;charset=utf-8;'});
  const url  = URL.createObjectURL(blob);
  const a    = document.createElement('a');
  a.href = url;
  // Nom de fichier : intersection active
  const sig = panel._sig || 'intersection';
  a.download = 'proteogen_upset_' + sig.replace(/[^a-z0-9]/gi,'_') + '.csv';
  a.click();
  URL.revokeObjectURL(url);
};

})();
/* ── FIN UPSET INTERACTIF ── */\
"""


def patch_upset_js(html):
    """Remplace la fonction initUpset() du template par la version interactive
    (clic barre → tableau protéines filtré). Injection juste avant </script>.
    Non-bloquant si la balise n'est pas trouvée.
    """
    import re as _re

    # Neutraliser l'ancienne initUpset dans le template (évite le double-init)
    # On renomme _origInitUpset pour garder une référence sans l'exécuter
    pattern = _re.compile(
        r'(?<!\w)function\s+initUpset\s*\(\s*\)\s*\{',
        _re.DOTALL
    )
    if pattern.search(html):
        html = pattern.sub('function _origInitUpset(){', html, count=1)
        print("  patch_upset_js : ancienne initUpset() neutralisée.")
    else:
        print("  [warn] patch_upset_js : initUpset() non trouvée dans le template, injection directe.")

    # Injecter la nouvelle version avant la dernière </script>
    pos = html.rfind('</script>')
    if pos == -1:
        print("  [warn] patch_upset_js : balise </script> introuvable.")
        return html
    html = html[:pos] + "\n\n" + DRAW_UPSET_JS + "\n" + html[pos:]
    print("  patch_upset_js : UpSet interactif injecté (%d chars)." % len(DRAW_UPSET_JS))
    return html


# ──────────────────────────────────────────────────────────────────────────────
# Patch JS — WGCNA Network (graphe étoile hub-centric, canvas)
# ──────────────────────────────────────────────────────────────────────────────

# drawWGCNANetwork :
#   • Hub central (grand cercle coloré par module) + spokes (petits cercles)
#   • Arêtes d'épaisseur proportionnelle au HubScore
#   • Sélecteur de module — redessine instantanément
#   • Tooltip au survol : gène, HubScore, description
#   • Avertissement "edges simulés — TOM non disponible" affiché dans le canvas
DRAW_WGCNA_NETWORK_JS = """\
function drawWGCNANetwork(){
  if(!WGCNA_NET||!Object.keys(WGCNA_NET).length)return;
  const sel=document.getElementById('wgcnaNetMod');
  const mod=sel?sel.value:Object.keys(WGCNA_NET)[0];
  const d=WGCNA_NET[mod];
  if(!d)return;

  const canvas=document.getElementById('wgcnaNetCanvas');
  if(!canvas)return;
  const dpr=devicePixelRatio||1;
  const W=680, H=580;
  canvas.width=W*dpr; canvas.height=H*dpr;
  canvas.style.width=W+'px'; canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.fillStyle='#0d1117'; ctx.fillRect(0,0,W,H);

  // Couleur du module (CSS named color si reconnue, sinon fallback)
  const MOD_COLORS={
    turquoise:'#1abc9c',blue:'#3498db',brown:'#a0522d',
    yellow:'#f1c40f',green:'#2ecc71',red:'#e74c3c',black:'#95a5a6',
    pink:'#e91e63',purple:'#9b59b6',orange:'#e67e22',grey:'#7f8c8d'
  };
  const modCol=MOD_COLORS[mod]||('#'+Math.abs(mod.split('').reduce(
    (h,c)=>Math.imul(31,h)+c.charCodeAt(0)|0,0)).toString(16).padStart(6,'0').slice(-6));

  const cx=W/2, cy=H/2;
  const R_HUB=32, R_NODE=10;
  const ORBIT=Math.min(cx,cy)-R_NODE-24;

  // Positions des spokes sur un cercle
  const nodes=d.nodes;
  const nN=nodes.length;
  const positions=nodes.map((_,i)=>{
    const a=-Math.PI/2 + (2*Math.PI*i/nN);
    return {x:cx+ORBIT*Math.cos(a), y:cy+ORBIT*Math.sin(a)};
  });

  // Arêtes hub → spokes (épaisseur ∝ HubScore normalisé)
  nodes.forEach((node,i)=>{
    const p=positions[i];
    ctx.strokeStyle=modCol+'88';
    ctx.lineWidth=Math.max(0.5, node.w*3.5);
    ctx.beginPath(); ctx.moveTo(cx,cy); ctx.lineTo(p.x,p.y); ctx.stroke();
  });

  // Nœuds spokes
  nodes.forEach((node,i)=>{
    const p=positions[i];
    const r=R_NODE*0.6 + R_NODE*0.4*node.w;
    ctx.beginPath(); ctx.arc(p.x,p.y,r,0,2*Math.PI);
    ctx.fillStyle=modCol+'cc'; ctx.fill();
    ctx.strokeStyle=modCol; ctx.lineWidth=1; ctx.stroke();
    // Label court
    ctx.fillStyle='#c9d1d9'; ctx.font='7px sans-serif'; ctx.textAlign='center';
    const lbl=(node.g||'').slice(0,8);
    ctx.fillText(lbl, p.x, p.y+r+10);
  });

  // Hub central
  ctx.beginPath(); ctx.arc(cx,cy,R_HUB,0,2*Math.PI);
  const grad=ctx.createRadialGradient(cx-4,cy-4,4,cx,cy,R_HUB);
  grad.addColorStop(0,'#ffffff33');
  grad.addColorStop(1,modCol);
  ctx.fillStyle=grad; ctx.fill();
  ctx.strokeStyle='#ffffff44'; ctx.lineWidth=2; ctx.stroke();
  // Label hub — gène + HubScore
  ctx.fillStyle='#ffffff'; ctx.font='bold 9px sans-serif'; ctx.textAlign='center';
  ctx.fillText((d.hub.g||'').slice(0,10), cx, cy-5);
  ctx.font='8px sans-serif'; ctx.fillStyle='#ffffffcc';
  ctx.fillText('hub '+d.hub.hub.toFixed(3), cx, cy+7);
  // Description hub (sous le cercle, wrappée sur 2 lignes max)
  if(d.hub.desc){
    ctx.fillStyle='#c9d1d9'; ctx.font='8px sans-serif'; ctx.textAlign='center';
    const words=d.hub.desc.split(' ');
    let line='', lines=[];
    words.forEach(w=>{
      const test=line?line+' '+w:w;
      if(ctx.measureText(test).width>140&&line){lines.push(line);line=w;}
      else line=test;
    });
    if(line)lines.push(line);
    lines.slice(0,2).forEach((l,i)=>ctx.fillText(l, cx, cy+R_HUB+14+i*13));
  }

  // Titre et infos
  ctx.fillStyle='#c9d1d9'; ctx.font='bold 12px sans-serif'; ctx.textAlign='center';
  ctx.fillText('Module '+mod+' — '+d.n_total+' protéines', W/2, 22);

  // Avertissement TOM
  ctx.fillStyle='#6e7681'; ctx.font='9px sans-serif'; ctx.textAlign='center';
  ctx.fillText('⚠ Edges simulés par HubScore (matrice TOM non exportée)', W/2, H-8);

  // Légende HubScore
  ctx.fillStyle='#8b949e'; ctx.font='8px sans-serif'; ctx.textAlign='left';
  ctx.fillText('Épaisseur arête ∝ HubScore  ('+d.h_min.toFixed(2)+' – '+d.h_max.toFixed(2)+')',
    10, H-8);

  // Stocker positions pour tooltip
  canvas._wgcnaPos=positions; canvas._wgcnaNodes=nodes;
  canvas._wgcnaHub={x:cx,y:cy,r:R_HUB,data:d.hub};
  canvas._wgcnaR=R_NODE;

  const tip=document.getElementById('wgcnaNetTip');
  canvas.onmousemove=function(e){
    const rect=canvas.getBoundingClientRect();
    const mx=(e.clientX-rect.left)*(W/rect.width);
    const my=(e.clientY-rect.top)*(H/rect.height);
    // Hub ?
    const h=canvas._wgcnaHub;
    const dh=Math.sqrt((mx-h.x)**2+(my-h.y)**2);
    if(dh<=h.r+4){
      if(tip){
        tip.style.display='block';
        tip.style.left=(e.clientX+12)+'px'; tip.style.top=(e.clientY-28)+'px';
        tip.innerHTML='<b>'+h.data.g+'</b> <span style="color:#f1c40f">HUB</span><br>'
          +'HubScore: '+h.data.hub+'<br>'
          +'<span style="color:#8b949e;font-size:9px;">'+h.data.desc+'</span>';
      }
      canvas.style.cursor='pointer';
      return;
    }
    // Spokes
    let best=null, bestD=1e9;
    canvas._wgcnaPos.forEach((p,i)=>{
      const dist=Math.sqrt((mx-p.x)**2+(my-p.y)**2);
      if(dist<bestD){bestD=dist;best=i;}
    });
    const rN=canvas._wgcnaR;
    if(best!==null && bestD<=rN+6){
      const node=canvas._wgcnaNodes[best];
      if(tip){
        tip.style.display='block';
        tip.style.left=(e.clientX+12)+'px'; tip.style.top=(e.clientY-28)+'px';
        tip.innerHTML='<b>'+node.g+'</b><br>'
          +'HubScore: '+node.hub+'<br>'
          +'<span style="color:#8b949e;font-size:9px;">'+node.desc+'</span>';
      }
      canvas.style.cursor='pointer';
    } else {
      if(tip)tip.style.display='none';
      canvas.style.cursor='default';
    }
  };
  canvas.onmouseleave=function(){
    if(tip)tip.style.display='none';
    canvas.style.cursor='default';
  };
}
function initWGCNANetwork(){
  const sel=document.getElementById('wgcnaNetMod');
  if(!sel||sel.options.length)return;
  Object.keys(WGCNA_NET).forEach(mod=>{
    const o=document.createElement('option'); o.value=mod;
    o.textContent=mod.charAt(0).toUpperCase()+mod.slice(1)+
      ' ('+WGCNA_NET[mod].n_total+' prot.)';
    sel.appendChild(o);
  });
  sel.onchange=drawWGCNANetwork;
  drawWGCNANetwork();
}\
"""


WGCNA_TOP5_PATCH_JS = """\
/* ── WGCNA TOP5 : description complète + export XLSX ── */
(function(){
  // Surcharge renderWGCNATable pour ajouter description complète + bouton export
  var _origRender = window.renderWGCNATable;
  window.renderWGCNATable = function(){
    if(_origRender) _origRender();
    // Ajouter/mettre à jour le bouton export WGCNA si absent
    if(!document.getElementById('wgcnaExportBtn')){
      var countDiv = document.getElementById('wgcnaCount');
      if(countDiv){
        var btn = document.createElement('button');
        btn.id = 'wgcnaExportBtn';
        btn.className = 'exp-btn';
        btn.textContent = '\\u2b07 Export XLSX';
        btn.style.cssText = 'margin-left:12px;padding:3px 10px;font-size:10px;';
        btn.onclick = exportWGCNATable;
        countDiv.parentNode && countDiv.parentNode.insertBefore(btn, countDiv.nextSibling);
      }
    }
  };

  // Export XLSX de la table WGCNA filtrée (gène, module, HubScore, description)
  window.exportWGCNATable = function(){
    var mod = document.getElementById('wgcnaMod').value;
    var search = (document.getElementById('wgcnaSearch').value||'').toLowerCase().trim();
    var rows = WGCNA_ROWS.slice();
    if(mod !== 'all') rows = rows.filter(function(r){ return r.mod === mod; });
    if(search) rows = rows.filter(function(r){
      return r.g.toLowerCase().includes(search) || (r.desc||'').toLowerCase().includes(search);
    });
    rows.sort(function(a,b){ return b.hub - a.hub; });
    var data = rows.map(function(r){
      return {
        'Gene':        r.g,
        'Module':      r.mod,
        'HubScore':    r.hub,
        'Description': r.desc || ''
      };
    });
    var fn = 'proteogen_WGCNA_' + (mod === 'all' ? 'tous_modules' : mod) + '.xlsx';
    if(typeof makeWb === 'function' && typeof xls === 'function'){
      xls(makeWb([{name: 'WGCNA_' + mod.slice(0,25), data: data}]), fn);
    } else {
      // Fallback CSV si XLSX non disponible
      var lines = ['Gene,Module,HubScore,Description'];
      data.forEach(function(r){
        lines.push([r.Gene, r.Module, r.HubScore,
          '"' + (r.Description||'').replace(/"/g,'""') + '"'].join(','));
      });
      var blob = new Blob([lines.join('\\n')], {type:'text/csv'});
      var a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = fn.replace('.xlsx','.csv');
      a.click();
    }
  };

  // Patch renderWGCNATable pour afficher la description complète dans le titre de colonne
  // et étendre la troncature à 80 chars au lieu de 55
  var _origRender2 = window.renderWGCNATable;
  window.renderWGCNATable = function(){
    var mod = document.getElementById('wgcnaMod').value;
    var search = (document.getElementById('wgcnaSearch').value||'').toLowerCase().trim();
    var rows = WGCNA_ROWS.slice();
    if(mod !== 'all') rows = rows.filter(function(r){ return r.mod === mod; });
    if(search) rows = rows.filter(function(r){
      return r.g.toLowerCase().includes(search) || (r.desc||'').toLowerCase().includes(search);
    });
    rows.sort(function(a,b){ return b.hub - a.hub; });
    var modColors = {blue:'#388bfd',turquoise:'#00bcd4',green:'#3fb950',yellow:'#d29922',
                     red:'#f85149',brown:'#cd853f',black:'#8b949e'};
    var tbody = document.getElementById('wgcnaBody');
    if(!tbody) return;
    tbody.innerHTML = '';
    rows.slice(0,300).forEach(function(r){
      var tr = document.createElement('tr');
      var desc = r.desc || '';
      var shortDesc = desc.length > 80 ? desc.slice(0,78) + '\\u2026' : desc;
      tr.innerHTML =
        '<td style="color:#c9d1d9;font-weight:600;">' + r.g + '</td>' +
        '<td><span style="display:inline-block;padding:2px 8px;border-radius:10px;font-size:10px;' +
          'background:' + (modColors[r.mod]||'#888') + '22;color:' + (modColors[r.mod]||'#888') +
          ';font-weight:600;">' + r.mod + '</span></td>' +
        '<td style="color:var(--bl);font-weight:600;">' + r.hub.toFixed(4) + '</td>' +
        '<td style="font-size:11px;color:var(--tx2);" title="' + desc.replace(/"/g,"&quot;") + '">' +
          shortDesc + '</td>';
      tbody.appendChild(tr);
    });
    var countEl = document.getElementById('wgcnaCount');
    if(countEl) countEl.textContent = rows.length + ' protéines (max 300 affichées)';
    // Bouton export
    if(!document.getElementById('wgcnaExportBtn')){
      var btn = document.createElement('button');
      btn.id = 'wgcnaExportBtn';
      btn.className = 'exp-btn';
      btn.textContent = '\\u2b07 Export XLSX';
      btn.style.cssText = 'margin-left:12px;padding:3px 10px;font-size:10px;';
      btn.onclick = exportWGCNATable;
      if(countEl) countEl.parentNode && countEl.parentNode.insertBefore(btn, countEl.nextSibling);
    }
  };
})();\
"""




# ──────────────────────────────────────────────────────────────────────────────
# Patch JS — GO Bubble plot + Table filtrée
# ──────────────────────────────────────────────────────────────────────────────

GO_BUBBLE_JS = """\
/* ── GO BUBBLE PLOT + TABLE ── */
function drawGOBubble(){
  if(!GO_BUBBLE||!Object.keys(GO_BUBBLE).length)return;
  const sel=document.getElementById('goBubbleContrast');
  const srcSel=document.getElementById('goBubbleSrc');
  const contrast=sel?sel.value:Object.keys(GO_BUBBLE)[0];
  const srcFilter=srcSel?srcSel.value:'all';
  let data=(GO_BUBBLE[contrast]||[]);
  if(srcFilter!=='all') data=data.filter(d=>d.src===srcFilter);
  data=data.slice().sort((a,b)=>b.lp-a.lp).slice(0,40);
  if(!data.length)return;

  const canvas=document.getElementById('goBubbleCanvas');
  if(!canvas)return;
  const dpr=devicePixelRatio||1;
  const lW=200,rW=80,tH=36,bH=50,PW=520,PH=Math.max(300,data.length*22);
  const W=lW+PW+rW,H=tH+PH+bH;
  canvas.width=W*dpr;canvas.height=H*dpr;
  canvas.style.width=W+'px';canvas.style.height=H+'px';
  const ctx=canvas.getContext('2d');
  ctx.scale(dpr,dpr);
  ctx.clearRect(0,0,W,H);

  // Plages
  const lpMax=Math.max(...data.map(d=>d.lp))*1.1||1;
  const isMax=Math.max(...data.map(d=>d.is))||1;
  const SRC_COL={'GO:BP':'#388bfd','GO:MF':'#3fb950','GO:CC':'#d29922','?':'#8b949e'};
  const R_MIN=5,R_MAX=20;

  function toX(lp){return lW+lp/lpMax*PW;}
  function toY(i){return tH+i*22+11;}
  function toR(is){return Math.min(R_MIN+(is/isMax)*(R_MAX-R_MIN), R_MAX);}

  // Titre
  ctx.fillStyle='#c9d1d9';ctx.font='bold 11px sans-serif';ctx.textAlign='center';
  ctx.fillText('GO Enrichment — '+contrast.replace(/_vs_/,' vs ').replace(/_/g,' '),
    lW+PW/2, 20);

  // Grille X
  ctx.strokeStyle='#21262d';ctx.lineWidth=1;
  [0,1,2,3,4].forEach(v=>{
    if(v>lpMax)return;
    const x=toX(v);
    ctx.beginPath();ctx.moveTo(x,tH);ctx.lineTo(x,tH+PH);ctx.stroke();
    ctx.fillStyle='#8b949e';ctx.font='8px sans-serif';ctx.textAlign='center';
    ctx.fillText(v,x,tH+PH+14);
  });
  ctx.fillStyle='#8b949e';ctx.font='9px sans-serif';ctx.textAlign='center';
  ctx.fillText('-log\\u2081\\u2080(p-value)',lW+PW/2,tH+PH+30);

  // Ligne p=0.05
  const x05=toX(-Math.log10(0.05));
  ctx.strokeStyle='#f85149';ctx.lineWidth=1;ctx.setLineDash([4,3]);
  ctx.beginPath();ctx.moveTo(x05,tH);ctx.lineTo(x05,tH+PH);ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle='#f85149';ctx.font='8px sans-serif';ctx.textAlign='left';
  ctx.fillText('p=0.05',x05+3,tH+10);

  // Points et labels
  data.forEach((d,i)=>{
    const x=toX(d.lp),y=toY(i);
    const r=toR(d.is);
    const col=SRC_COL[d.src]||'#8b949e';
    // z_score : opacité/couleur du remplissage
    const alpha=d.z>=0?'cc':'66';
    ctx.beginPath();ctx.arc(x,y,r,0,2*Math.PI);
    ctx.fillStyle=col+alpha;ctx.fill();
    ctx.strokeStyle=col;ctx.lineWidth=1;ctx.stroke();
    // Label terme
    ctx.fillStyle='#c9d1d9';ctx.font='9px sans-serif';ctx.textAlign='right';
    const lbl=d.term.length>32?d.term.slice(0,30)+'\\u2026':d.term;
    ctx.fillText(lbl,lW-8,y+3);
    // Source badge
    ctx.fillStyle=col;ctx.font='bold 7px sans-serif';ctx.textAlign='left';
    ctx.fillText(d.src.replace('GO:',''),lW+PW+4,y+3);
  });

  // Légende taille (valeurs dérivées du max réel pour rester cohérentes)
  const lgX=lW+PW+30,lgY=tH+60;
  ctx.fillStyle='#8b949e';ctx.font='bold 8px sans-serif';ctx.textAlign='center';
  ctx.fillText('Taille',lgX,lgY-14);
  var _legVals=[Math.max(1,Math.round(isMax)),
                Math.max(1,Math.round(isMax/2)),
                Math.max(1,Math.round(isMax/4))];
  _legVals=_legVals.filter((v,i,a)=>a.indexOf(v)===i);  // dédoublonne
  var _accY=lgY;
  _legVals.forEach(function(is){
    var r=toR(is);
    var y2=_accY+r;
    ctx.beginPath();ctx.arc(lgX,y2,r,0,2*Math.PI);
    ctx.fillStyle='#388bfd44';ctx.fill();
    ctx.strokeStyle='#388bfd';ctx.lineWidth=1;ctx.stroke();
    ctx.fillStyle='#8b949e';ctx.font='7px sans-serif';ctx.textAlign='left';
    ctx.fillText(is,lgX+r+5,y2+3);
    _accY=y2+r+10;
  });

  // Stocker pour tooltip
  canvas._goData=data;canvas._goToX=toX;canvas._goToY=toY;canvas._goToR=toR;
  const tip=document.getElementById('goBubbleTip');
  canvas.onmousemove=function(e){
    const rect=canvas.getBoundingClientRect();
    const mx=(e.clientX-rect.left)*(W/rect.width);
    const my=(e.clientY-rect.top)*(H/rect.height);
    let best=null,bestD=1e9;
    data.forEach((d,i)=>{
      const dx=toX(d.lp)-mx,dy=toY(i)-my;
      const dist=Math.sqrt(dx*dx+dy*dy);
      if(dist<bestD){bestD=dist;best={d,i};}
    });
    if(tip&&best&&bestD<toR(best.d.is)+6){
      tip.style.display='block';
      tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY-28)+'px';
      const d=best.d;
      tip.innerHTML='<b>'+d.term+'</b><br>'+d.src+' · '+d.tid+'<br>'+
        'p='+d.p.toExponential(2)+'  z='+d.z+'  n='+d.is;
    } else if(tip){tip.style.display='none';}
  };
  canvas.onmouseleave=function(){if(tip)tip.style.display='none';};
  // Clic → filtrer table
  canvas.onclick=function(e){
    const rect=canvas.getBoundingClientRect();
    const mx=(e.clientX-rect.left)*(W/rect.width);
    const my=(e.clientY-rect.top)*(H/rect.height);
    data.forEach((d,i)=>{
      const dx=toX(d.lp)-mx,dy=toY(i)-my;
      if(Math.sqrt(dx*dx+dy*dy)<toR(d.is)+6){
        var si=document.getElementById('goTableSearch');
        if(si){si.value=d.term;filterGOTable();}
      }
    });
  };
}

function initGOBubble(){
  const sel=document.getElementById('goBubbleContrast');
  if(!sel||sel.options.length)return;
  Object.keys(GO_BUBBLE).forEach(c=>{
    const o=document.createElement('option');o.value=c;
    o.textContent=c.replace(/_vs_/,' vs ').replace(/_/g,' ');sel.appendChild(o);
  });
  sel.onchange=function(){drawGOBubble();renderGOTable();};
  const srcSel=document.getElementById('goBubbleSrc');
  if(srcSel){srcSel.onchange=function(){drawGOBubble();renderGOTable();};}
  drawGOBubble();
  renderGOTable();
}

function renderGOTable(){
  if(!GO_BUBBLE||!Object.keys(GO_BUBBLE).length)return;
  const sel=document.getElementById('goBubbleContrast');
  const srcSel=document.getElementById('goBubbleSrc');
  const contrast=sel?sel.value:Object.keys(GO_BUBBLE)[0];
  const srcFilter=srcSel?srcSel.value:'all';
  const q=(document.getElementById('goTableSearch').value||'').toLowerCase().trim();
  let data=(GO_BUBBLE[contrast]||[]).slice();
  if(srcFilter!=='all') data=data.filter(d=>d.src===srcFilter);
  if(q) data=data.filter(d=>d.term.toLowerCase().includes(q)||d.tid.toLowerCase().includes(q));
  data.sort((a,b)=>a.p-b.p);
  const SRC_COL={'GO:BP':'#388bfd','GO:MF':'#3fb950','GO:CC':'#d29922'};
  const tbody=document.getElementById('goTableBody');
  if(!tbody)return;
  tbody.innerHTML='';
  data.slice(0,200).forEach(d=>{
    const tr=document.createElement('tr');
    tr.style.borderBottom='1px solid var(--bd)';
    tr.style.cursor='pointer';
    const col=SRC_COL[d.src]||'#888';
    tr.innerHTML=
      '<td style="padding:4px 6px;font-size:10px;"><span style="display:inline-block;'+
        'padding:1px 6px;border-radius:10px;background:'+col+'22;color:'+col+
        ';font-weight:700;font-size:9px;">'+d.src.replace('GO:','')+'</span></td>'+
      '<td style="padding:4px 6px;font-size:11px;color:var(--tx);">'+d.term+'</td>'+
      '<td style="padding:4px 6px;font-family:monospace;font-size:9px;color:var(--tx2);">'+d.tid+'</td>'+
      '<td style="padding:4px 6px;text-align:right;font-size:10px;color:var(--tx);">'+d.p.toExponential(2)+'</td>'+
      '<td style="padding:4px 6px;text-align:right;font-size:10px;color:'+(d.z>=0?'#3fb950':'#388bfd')+';font-weight:700;">'+d.z.toFixed(3)+'</td>'+
      '<td style="padding:4px 6px;text-align:center;font-size:10px;color:var(--tx);">'+d.is+'</td>'+
      '<td style="padding:4px 6px;text-align:right;font-size:10px;color:var(--tx2);">'+(d.gr*100).toFixed(1)+'%</td>';
    tr.onclick=function(){
      var si=document.getElementById('goTableSearch');
      if(si){si.value=d.term;filterGOTable();}
    };
    tbody.appendChild(tr);
  });
  const cnt=document.getElementById('goTableCount');
  if(cnt)cnt.textContent=data.length+' termes';
}

window.filterGOTable=function(){renderGOTable();drawGOBubble();};

function exportGOBubble(){
  exportCanvas('goBubbleCanvas','proteogen_GO_bubble_'+
    (document.getElementById('goBubbleContrast').value||''));
}
function exportGOTable(){
  const sel=document.getElementById('goBubbleContrast');
  const c=sel?sel.value:Object.keys(GO_BUBBLE)[0];
  const data=(GO_BUBBLE[c]||[]).slice().sort((a,b)=>a.p-b.p);
  const lines=['Source,Term,Term_ID,p_value,z_score,intersection_size,gene_ratio,genes'];
  data.forEach(d=>{
    const esc=v=>'"'+String(v||'').replace(/"/g,'""')+'"';
    lines.push([d.src,esc(d.term),d.tid,d.p,d.z,d.is,d.gr,esc(d.genes)].join(','));
  });
  const blob=new Blob([lines.join('\\n')],{type:'text/csv'});
  const a=document.createElement('a');
  a.href=URL.createObjectURL(blob);
  a.download='proteogen_GO_'+c+'.csv';
  a.click();
}\
"""


def patch_go_bubble_js(html):
    """Injecte drawGOBubble, renderGOTable et les helpers GO avant </script>."""
    pos = html.rfind('</script>')
    if pos == -1:
        print("  [warn] patch_go_bubble_js : </script> introuvable.")
        return html
    html = html[:pos] + "\n\n" + GO_BUBBLE_JS + "\n" + html[pos:]
    print("  patch_go_bubble_js : GO bubble + table injectés.")
    return html


# HTML à injecter dans #gopage (remplace la section images statiques)
GO_BUBBLE_HTML = """\
<!-- Suppression des images GO statiques, remplacement par bubble plot + table -->
<div class="panel">
  <div class="ph">
    <h2>GO Enrichment — Bubble Plot</h2>
    <div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap;">
      <select id="goBubbleContrast" style="background:var(--bg2);border:1px solid var(--bd);color:var(--tx);border-radius:4px;padding:3px 8px;font-size:11px;"></select>
      <select id="goBubbleSrc" style="background:var(--bg2);border:1px solid var(--bd);color:var(--tx);border-radius:4px;padding:3px 8px;font-size:11px;">
        <option value="all">Toutes sources</option>
        <option value="GO:BP">BP</option>
        <option value="GO:MF">MF</option>
        <option value="GO:CC">CC</option>
      </select>
      <button class="exp-btn" onclick="exportGOBubble()" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button>
    </div>
  </div>
  <canvas id="goBubbleCanvas" style="display:block;width:100%;"></canvas>
  <div id="goBubbleTip" style="position:fixed;display:none;background:var(--bg2);border:1px solid var(--bd);border-radius:6px;padding:6px 10px;font-size:11px;color:var(--tx);pointer-events:none;z-index:999;max-width:320px;"></div>
</div>

<div class="panel">
  <div class="ph">
    <h2>Termes GO enrichis</h2>
    <div style="display:flex;gap:8px;align-items:center;">
      <input id="goTableSearch" type="text" placeholder="Filtrer terme / ID..." oninput="filterGOTable()"
        style="background:var(--bg2);border:1px solid var(--bd);color:var(--tx);border-radius:4px;padding:3px 8px;font-size:11px;width:200px;">
      <span id="goTableCount" style="font-size:10px;color:var(--tx2);"></span>
      <button class="exp-btn" onclick="exportGOTable()" style="padding:3px 9px;font-size:10px;">&#11015; CSV</button>
    </div>
  </div>
  <div style="overflow-y:auto;max-height:400px;">
    <table style="width:100%;border-collapse:collapse;">
      <thead><tr style="position:sticky;top:0;background:var(--bg2);">
        <th style="padding:4px 6px;text-align:left;font-size:10px;color:var(--tx2);border-bottom:1px solid var(--bd);font-weight:600;">Source</th>
        <th style="padding:4px 6px;text-align:left;font-size:10px;color:var(--tx2);border-bottom:1px solid var(--bd);font-weight:600;">Terme</th>
        <th style="padding:4px 6px;text-align:left;font-size:10px;color:var(--tx2);border-bottom:1px solid var(--bd);font-weight:600;">ID</th>
        <th style="padding:4px 6px;text-align:right;font-size:10px;color:var(--tx2);border-bottom:1px solid var(--bd);font-weight:600;">p-value</th>
        <th style="padding:4px 6px;text-align:right;font-size:10px;color:var(--tx2);border-bottom:1px solid var(--bd);font-weight:600;">z-score</th>
        <th style="padding:4px 6px;text-align:center;font-size:10px;color:var(--tx2);border-bottom:1px solid var(--bd);font-weight:600;">N gènes</th>
        <th style="padding:4px 6px;text-align:right;font-size:10px;color:var(--tx2);border-bottom:1px solid var(--bd);font-weight:600;">Gene ratio</th>
      </tr></thead>
      <tbody id="goTableBody"></tbody>
    </table>
  </div>
</div>\
"""

def patch_wgcna_network_js(html):
    """Injecte drawWGCNANetwork(), initWGCNANetwork() et WGCNA_TOP5_PATCH_JS avant </script>."""
    pos = html.rfind('</script>')
    if pos == -1:
        print("  [warn] patch_wgcna_network_js : </script> introuvable.")
        return html
    inject = "\n\n" + DRAW_WGCNA_NETWORK_JS + "\n\n" + WGCNA_TOP5_PATCH_JS + "\n"
    html = html[:pos] + inject + html[pos:]
    print("  patch_wgcna_network_js : network WGCNA + top5 patch injectés (%d chars)." % len(inject))
    return html


# ──────────────────────────────────────────────────────────────────────────────
# Patch HTML — injection des canvas/UI de l'onglet QC dans le template
# ──────────────────────────────────────────────────────────────────────────────

# Bloc HTML complet de l'onglet QC : inséré dans le premier <div class="tab-content">
# du dashboard. Les canvas sont liés aux fonctions JS injectées par patch_qc_figures()
# et patch_umap_js(). Les appels JS se font via onclick/onload.
QC_TAB_HTML = (
    '<div class="panel">\n'
    '  <div class="ph"><h2>&#9989; Rapport de Contr\u00f4le Qualit\u00e9</h2></div>\n'
    '  <div id="qcReportPanel"></div>\n'
    '</div>\n'
) + RANKED_CV_HTML + """
<div class="panel">
  <div class="ph"><h2>UMAP</h2><button class="exp-btn" onclick="exportCanvas('umapCanvas','proteogen_UMAP')" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button></div>
  <canvas id="umapCanvas" style="display:block;width:100%;"></canvas>
  <div id="umapTip" style="position:fixed;display:none;background:var(--bg2);border:1px solid var(--bd);border-radius:6px;padding:6px 10px;font-size:11px;color:var(--tx);pointer-events:none;z-index:999;"></div>
</div>
<div class="row2">
  <div class="panel">
    <div class="ph"><h2>Distribution LFQ (post-imputation)</h2><button class="exp-btn" onclick="exportCanvas('violinCanvas','proteogen_violin')" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button></div>
    <canvas id="violinCanvas" style="display:block;width:100%;"></canvas>
  </div>
  <div class="panel">
    <div class="ph"><h2>Valeurs manquantes (avant imputation)</h2><button class="exp-btn" onclick="exportCanvas('missingCanvas','proteogen_missing')" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button></div>
    <canvas id="missingCanvas" style="display:block;width:100%;"></canvas>
  </div>
</div>
<div class="row2">
  <div class="panel">
    <div class="ph"><h2>MA Plot</h2><div style="display:flex;gap:6px;align-items:center;"><select id="maContrast" style="background:var(--bg2);border:1px solid var(--bd);color:var(--tx);border-radius:4px;padding:3px 8px;font-size:11px;"></select><button class="exp-btn" onclick="exportCanvas('maCanvas','proteogen_MA_'+document.getElementById('maContrast').value)" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button></div></div>
    <canvas id="maCanvas" style="display:block;width:100%;"></canvas>
    <div id="maTip" style="position:fixed;display:none;background:var(--bg2);border:1px solid var(--bd);border-radius:6px;padding:6px 10px;font-size:11px;color:var(--tx);pointer-events:none;z-index:999;"></div>
  </div>
  <div class="panel">
    <div class="ph"><h2>Corr. inter-r\u00e9plicats</h2><div style="display:flex;gap:6px;align-items:center;"><select id="repScatterCond" style="background:var(--bg2);border:1px solid var(--bd);color:var(--tx);border-radius:4px;padding:3px 8px;font-size:11px;"></select><button class="exp-btn" onclick="exportCanvas('repScatterCanvas','proteogen_repl_'+document.getElementById('repScatterCond').value)" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button></div></div>
    <canvas id="repScatterCanvas" style="display:block;width:100%;"></canvas>
  </div>
</div>
"""

WGCNA_NET_HTML = """\
<div class="panel">
  <div class="ph"><h2>Network de co-expression (hub-centric)</h2><div style="display:flex;gap:6px;align-items:center;"><select id="wgcnaNetMod" style="background:var(--bg2);border:1px solid var(--bd);color:var(--tx);border-radius:4px;padding:3px 8px;font-size:11px;"></select><button class="exp-btn" onclick="exportCanvas('wgcnaNetCanvas','proteogen_WGCNA_'+document.getElementById('wgcnaNetMod').value)" style="padding:3px 9px;font-size:10px;">&#11015; PNG</button></div></div>
  <canvas id="wgcnaNetCanvas" style="display:block;width:100%;"></canvas>
  <div id="wgcnaNetTip" style="position:fixed;display:none;background:var(--bg2);border:1px solid var(--bd);border-radius:6px;padding:6px 10px;font-size:11px;color:var(--tx);pointer-events:none;z-index:999;"></div>
  <p style="font-size:9px;color:var(--tx2);margin-top:8px;">&#9888; Edges simul\u00e9s par HubScore \u2014 matrice TOM non export\u00e9e</p>
</div>
"""


# Script d'initialisation des figures QC (appelé au clic sur l'onglet QC)
QC_INIT_SCRIPT = """\
<script>
(function(){
  var _qcInit = false;
  function onQcTab(){
    if(_qcInit) return;
    _qcInit = true;
    // requestAnimationFrame garantit que le layout est calculé (offsetWidth > 0)
    // avant de dessiner — critique sous Firefox (display:none -> block)
    // Chaque rendu est isolé : une fonction qui plante ne doit pas emporter
    // tout l'onglet (sinon clic = page figée sans message).
    function _safe(name, fn){
      try { if(typeof fn === 'function') fn(); }
      catch(e){ console.error('[QC] '+name+' a échoué :', e); }
    }
    requestAnimationFrame(function(){
      _safe('drawUMAP', window.drawUMAP);
      _safe('drawViolinPlot', window.drawViolinPlot);
      _safe('initMAPlot', window.initMAPlot);
      _safe('initRepScatter', window.initRepScatter);
      _safe('drawMissingValues', window.drawMissingValues);
      _safe('drawRankedAbundance', window.drawRankedAbundance);
      _safe('initCVPlot', window.initCVPlot);
      _safe('renderQCReport', window.renderQCReport);
      requestAnimationFrame(function(){
        _safe('drawPCA', window.drawPCA);
        _safe('drawScatter', window.drawScatter);
      });
    });
  }
  // Monkey-patch showPage — exécuté immédiatement (en fin de body)
  var _orig = window.showPage;
  window.showPage = function(id){
    if(_orig) _orig(id);
    if(id === 'qcpage') onQcTab();
    if(id === 'wgcna' && typeof initWGCNANetwork === 'function') initWGCNANetwork();
    if(id === 'overview' && typeof renderAbout === 'function') renderAbout();
    if(id === 'gopage'){
      if(typeof initGO === 'function') initGO();
      if(typeof initGOBubble === 'function') initGOBubble();
      if(typeof drawGOHeatmap === 'function') requestAnimationFrame(function(){
        drawGOHeatmap();
        var b=document.getElementById('goHmExportBtn');
        if(b) b.onclick=function(){exportCanvas('goHmCanvas','proteogen_GO_heatmap');};
        var s=document.getElementById('goHmSrc');
        if(s) s.onchange=drawGOHeatmap;
        // Brancher renderGO sur le sélecteur de contrastes GO
        var goBtns=document.querySelectorAll('#goBtns button');
        goBtns.forEach(function(btn){
          btn.addEventListener('click',function(){
            if(typeof renderGO==='function') setTimeout(function(){
              var c=typeof currentGOContrast!=='undefined'?currentGOContrast:null;
              renderGO(c);
            },50);
          });
        });
      });
    }
    if(id === 'volcano' && typeof initVolcanoComp === 'function') initVolcanoComp();
  };
  // Appel initial si overview est déjà active au chargement
  if(typeof renderAbout === 'function') renderAbout();
})();
</script>\
"""




def patch_gopage_html(html):
    """Supprime les images GO statiques et le bloc goFigImg/goFigsPanel,
    puis injecte le bubble plot + table dans #gopage.
    """
    import re as _re

    # Supprimer goFigsPanel s'il reste
    p1 = _re.compile(r'<div[^>]*id=[\x22\x27]goFigsPanel[\x22\x27][^>]*>.*?</div>\s*\n\s*</div>', _re.DOTALL)
    if p1.search(html):
        html = p1.sub('', html, count=1)
        print("  patch_gopage_html : goFigsPanel supprimé.")

    # Supprimer le bloc goFigNoData + goFigImg (résidu si pas dans goFigsPanel)
    p2 = _re.compile(r'<div[^>]*id=[\x22\x27]goFigNoData[\x22\x27].*?</div>\s*\n\s*<div[^>]*>\s*\n\s*<img[^>]*id=[\x22\x27]goFigImg[\x22\x27].*?</div>\s*\n\s*</div>', _re.DOTALL)
    if p2.search(html):
        html = p2.sub('', html, count=1)
        print("  patch_gopage_html : bloc goFigImg supprimé.")

    # Supprimer les radio buttons goFigTabs
    p3 = _re.compile(r'<div[^>]*id=[\x22\x27]goFigTabs[\x22\x27].*?</div>', _re.DOTALL)
    if p3.search(html):
        html = p3.sub('', html, count=1)

    # Neutraliser les références JS aux éléments GO supprimés du DOM
    for stub in [
        ("const img=document.getElementById('goFigImg');",
         "const img=document.getElementById('goFigImg')||{style:{},src:'',alt:''};"),
        ("const noData=document.getElementById('goFigNoData');",
         "const noData=document.getElementById('goFigNoData')||{style:{}};"),
        ("document.querySelectorAll('input[name=\"goFig\"]').forEach(r=>r.onchange=",
         "//removed: document.querySelectorAll('input[name=\"goFig\"]').forEach(r=>r.onchange="),
    ]:
        html = html.replace(stub[0], stub[1])

    # Supprimer les panels Chart.js natifs (goZChart + goSrcChart)
    # remplacés par le bubble plot et la heatmap cross-contrastes
    for panel_id in ['goZChart', 'goSrcChart']:
        p_chart = _re.compile(
            r'<div class=["\x27]panel["\x27][^>]*>(?:(?!</div>).)*' + panel_id + r'.*?</div>\s*</div>',
            _re.DOTALL
        )
        if p_chart.search(html):
            html = p_chart.sub('', html, count=1)
            print(f"  patch_gopage_html : panel {panel_id} supprimé.")

    # Injecter GO_BUBBLE_HTML avant <!-- HEATMAP -->
    m = _re.search(r'(</div>)\s*\n\s*<!-- HEATMAP -->', html)
    if m:
        insert_pos = m.start(1)
        html = html[:insert_pos] + '\n' + GO_BUBBLE_HTML + '\n' + html[insert_pos:]
        print("  patch_gopage_html : GO bubble + table injectés dans #gopage.")
    else:
        print("  [warn] patch_gopage_html : fermeture gopage introuvable.")

    return html


def patch_qc_tab_html(html):
    """Injecte les canvas QC dans #qcpage et le network WGCNA dans #wgcna.

    Pour #qcpage : le template a une structure imbriquée (row2 > panel > ph > ...).
    On cherche la fermeture du dernier panel natif (Scatter) et on insère après.

    Pour #wgcna : on cherche la fermeture du premier panel natif et on insère après.

    Idempotent sur umapCanvas.
    """
    import re as _re

    if 'id="umapCanvas"' in html:
        print("  patch_qc_tab_html : déjà injecté, skip.")
        if 'onQcTab' not in html:
            pos = html.find('<!--PROTEOGEN_BODY_END-->')
            if pos == -1:
                real_body = _re.search(r'</body\s*>\s*</html\s*>', html, _re.IGNORECASE)
                pos = real_body.start() if real_body else -1
            if pos != -1:
                html = html[:pos] + QC_INIT_SCRIPT + '\n' + html[pos:]
        return html

    # ── 0. Injection du panneau About dans #overview ───────────────────────────
    if 'id="aboutPanel"' not in html:
        idx_ov = html.find('id="overview"')
        if idx_ov > 0:
            # Insérer après le div#kpiGrid (premier panel de l'overview)
            idx_kpi = html.find('id="kpiGrid"', idx_ov)
            if idx_kpi > 0:
                end_kpi = html.find('</div>', idx_kpi) + 6
                html = html[:end_kpi] + '\n' + ABOUT_HTML + html[end_kpi:]
                print('  patch_qc_tab_html : panneau About injecté dans #overview.')

    # ── 1. Injection dans #qcpage ─────────────────────────────────────────────
    # La jonction entre qcpage et la page suivante est : </div>\n\n<!-- WGCNA -->\n
    # On insère juste avant cette jonction
    m_qc = _re.search(r'(id=["\']qcpage["\'].*?)(</div>\s*\n\s*<!--[^>]*-->\s*\n\s*<div[^>]+id=["\']wgcna["\'])',
                      html, _re.DOTALL)
    if m_qc:
        insert_pos = m_qc.start(2)
        html = html[:insert_pos] + '\n' + QC_TAB_HTML + '\n' + html[insert_pos:]
        print("  patch_qc_tab_html : figures QC injectées dans #qcpage.")
    else:
        # fallback : avant le marqueur de body
        pos = html.find('<!--PROTEOGEN_BODY_END-->')
        if pos == -1:
            real_body = _re.search(r'</body\s*>\s*</html\s*>', html, _re.IGNORECASE)
            pos = real_body.start() if real_body else -1
        if pos != -1:
            html = html[:pos] + QC_TAB_HTML + '\n' + html[pos:]
        print("  patch_qc_tab_html : figures QC injectées avant </body> (fallback).")

    # ── 2. Injection dans #wgcna ──────────────────────────────────────────────
    if 'id="wgcnaNetCanvas"' not in html:
        idx_wg = html.find('id="wgcna"')
        if idx_wg == -1:
            idx_wg = html.find("id='wgcna'")
        if idx_wg > 0:
            # La page wgcna est la dernière — elle se ferme juste avant <script>
            chunk = html[idx_wg:]
            m_end = _re.search(r'(</div>\s*\n)\s*\n\s*<script', chunk)
            if m_end:
                insert_pos = idx_wg + m_end.start(1)
                html = html[:insert_pos] + '\n' + WGCNA_NET_HTML + '\n' + html[insert_pos:]
                print("  patch_qc_tab_html : network WGCNA injecté dans #wgcna.")
            else:
                print("  [warn] patch_qc_tab_html : fermeture de #wgcna introuvable.")

    # ── 3. Script d'init — injecté tel quel, wgcna déjà géré dedans ─────────
    pos_body = html.find('<!--PROTEOGEN_BODY_END-->')
    # Fallback si le marqueur n'est pas là : chercher </body> suivi de </html>
    if pos_body == -1:
        real_body = _re.search(r'</body\s*>\s*</html\s*>', html, _re.IGNORECASE)
        pos_body = real_body.start() if real_body else -1
    if pos_body != -1:
        html = html[:pos_body] + QC_INIT_SCRIPT + '\n' + html[pos_body:]
        print("  patch_qc_tab_html : QC_INIT_SCRIPT injecté avant </body>.")

    return html


# ──────────────────────────────────────────────────────────────────────────────
# Injection + écriture
# ──────────────────────────────────────────────────────────────────────────────
def inject(template_html, data):
    # Règle du skill : concaténation, JAMAIS de f-string pour injecter du JSON
    block = (
        "const CONTRASTS="     + json.dumps(data['contrasts'], ensure_ascii=False)   + ";\n"
        "const PROTEOME_DATA=" + json.dumps(data['proteome_data'], ensure_ascii=False) + ";\n"
        "const GO_DATA="       + json.dumps(data['go_data'], ensure_ascii=False)     + ";\n"
        "const HEATMAP="       + json.dumps(data['heatmap'], ensure_ascii=False)     + ";\n"
        "const SAMPLE_COLS="   + json.dumps(data['sample_cols'], ensure_ascii=False) + ";\n"
        "const CONDITIONS="    + json.dumps(data['conditions'], ensure_ascii=False)  + ";\n"
        "const COND_MAP="      + json.dumps(data['cond_map'], ensure_ascii=False)    + ";\n"
        "const STATS="         + json.dumps(data['stats'], ensure_ascii=False)       + ";\n"
        "const INTERSECTIONS=" + json.dumps(data['intersections'], ensure_ascii=False) + ";\n"
        "const UPSET_DATA="    + json.dumps(data['upset'], ensure_ascii=False)       + ";\n"
        "const TOP20_PI="      + json.dumps(data['top20_pi'], ensure_ascii=False)    + ";\n"
        "const WGCNA_ROWS="    + json.dumps(data['wgcna_rows'], ensure_ascii=False)  + ";\n"
        "const WGCNA_TOP5="    + json.dumps(data['wgcna_top5'], ensure_ascii=False)  + ";\n"
        "const PCA_DATA="      + json.dumps(data['pca'], ensure_ascii=False)         + ";\n"
        "const WGCNA_HUB_MAP="  + json.dumps(data['wgcna_hub_map'], ensure_ascii=False) + ";\n"
        "const GO_NETWORK="    + json.dumps(data['go_network'], ensure_ascii=False)  + ";\n"
        "const GO_HEATMAP="    + json.dumps(data['go_heatmap'], ensure_ascii=False)  + ";\n"
        "const GO_BUBBLE="     + json.dumps(data['go_bubble'], ensure_ascii=False)   + ";\n"
        "const GO_IMAGES="     + json.dumps({}, ensure_ascii=False)                 + ";\n"  # images GO supprimées
        "const CHORD_DATA="    + json.dumps(data['chord'], ensure_ascii=False)       + ";\n"
        "const WGCNA_MODS="       + json.dumps(data['wgcna_mods'], ensure_ascii=False)    + ";\n"
        "const MODULE_TRAIT="      + json.dumps(data['module_trait'], ensure_ascii=False)  + ";\n"
        "const VIOLIN_DATA="      + json.dumps(data['violin'], ensure_ascii=False)         + ";\n"
        "const MISSING_DATA="     + json.dumps(data['missing'], ensure_ascii=False)        + ";\n"
        "const WGCNA_NET="        + json.dumps(data['wgcna_net'], ensure_ascii=False)        + ";\n"
        "const UMAP_DATA="        + json.dumps(data['umap'], ensure_ascii=False)           + ";\n"
        "const REP_SCATTER_DATA=" + json.dumps(data['rep_scatter'], ensure_ascii=False)   + ";\n"
        "const RANKED_DATA="      + json.dumps(data['ranked'], ensure_ascii=False)        + ";\n"
        "const CV_DATA="          + json.dumps(data['cv_plot'], ensure_ascii=False)       + ";\n"
        "const ABOUT_DATA="       + json.dumps(data['about'], ensure_ascii=False)         + ";\n"
        "const QC_REPORT="        + json.dumps(data['qc_report'], ensure_ascii=False)     + ";\n"
        # PROTEOME_DATA contient lp précalculé — VOLCANO et SCATTER_DATA sont des vues dérivées.
        # Les points dessinés sont SOUS-ÉCHANTILLONNÉS (max ~1500/contraste) pour
        # éviter de figer le navigateur au rendu de l'onglet QC/PCA.
        "var _NS_MAX=2500;"
        "function _subsample(arr){if(arr.length<=_NS_MAX)return arr;"
        "var step=Math.ceil(arr.length/_NS_MAX),out=[];"
        "for(var i=0;i<arr.length;i+=step)out.push(arr[i]);return out;}\n"
        "const MA_DATA=(function(){"
        "var r={};Object.keys(PROTEOME_DATA).forEach(function(c){"
        "var pts=_subsample(PROTEOME_DATA[c]).map(function(p){"
        "return {a:Math.round((p.x+p.y)/2*1000)/1000,m:p.d,s:p.s,g:p.g,n:p.name,p:p.p,pa:p.padj};});"
        "r[c]={pts:pts,lfc_thresh:0.26,"
        "label_a:c.split('_vs_')[0].replace(/_/g,' '),"
        "label_b:(c.split('_vs_')[1]||'').replace(/_/g,' ')};"
        "});return r;})();\n"
        "const SCATTER_DATA=(function(){"
        "var r={};Object.keys(PROTEOME_DATA).forEach(function(c){"
        "var pts=PROTEOME_DATA[c];var xa=pts.map(function(p){return p.x;});"
        "var xb=pts.map(function(p){return p.y;});"
        "var mask=pts.map(function(p){return p.x>0&&p.y>0;});"
        "var xa2=xa.filter(function(_,i){return mask[i];});"
        "var xb2=xb.filter(function(_,i){return mask[i];});"
        "var n=xa2.length,sx=0,sy=0,sxy=0,sx2=0,sy2=0;"
        "for(var i=0;i<n;i++){sx+=xa2[i];sy+=xb2[i];sxy+=xa2[i]*xb2[i];sx2+=xa2[i]*xa2[i];sy2+=xb2[i]*xb2[i];}"
        "var rr=n>1?(sxy-sx*sy/n)/Math.sqrt((sx2-sx*sx/n)*(sy2-sy*sy/n)):0;"
        # r calculé sur TOUS les points ; seuls les points DESSINÉS sont échantillonnés
        "var mapped=_subsample(pts).map(function(p){"
        "return {x:p.x,y:p.y,s:p.s,g:p.g,n:p.name,d:p.desc,lfc:p.d,pv:p.p,pa:p.padj};});"
        "r[c]={pts:mapped,r:Math.round(rr*10000)/10000,"
        "label_a:c.split('_vs_')[0].replace(/_/g,' '),"
        "label_b:(c.split('_vs_')[1]||'').replace(/_/g,' ')};"
        "});return r;})();\n"
        # VOLCANO = PROTEOME_DATA (même structure, lp déjà présent)
        "const VOLCANO=PROTEOME_DATA;\n"
    )
    if '__DATA_INJECTION__' not in template_html:
        raise ValueError("Placeholder __DATA_INJECTION__ absent du template.")
    # Sécurité : empêcher une chaîne '</script>' présente dans une valeur (desc, terme...)
    # de fermer prématurément la balise <script>. '<\/' reste valide en JS.
    block = block.replace('</', '<\\/')
    return template_html.replace('__DATA_INJECTION__', block)


def main():
    ap = argparse.ArgumentParser(description="Génère le dashboard protéomique HTML depuis les XLSX des scripts R.")
    ap.add_argument('--uploads', default='/Data/Dashboard/upload', help="Dossier d'auto-détection des XLSX")
    ap.add_argument('--stats', help="Chemin rapportstatistique.xlsx (obligatoire si pas d'auto-détection)")
    ap.add_argument('--go', help="Chemin rapportenrichissementGO.xlsx (optionnel)")
    ap.add_argument('--wgcna', help="Chemin WGCNA_Emotional.xlsx (optionnel)")
    ap.add_argument('--umap',  help="Chemin umap.xlsx avec sheet UMAP_Output (optionnel)")
    ap.add_argument('--template', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dashboard_template.html'))
    ap.add_argument('--out', default='/Data/Dashboard/outputs/proteogen_dashboard.html')
    ap.add_argument('--params-json', dest='params_json', default=None,
                    help="JSON des seuils choisis dans le pipeline (optionnel)")
    args = ap.parse_args()

    # Lecture des seuils choisis (si fournis par le pipeline)
    chosen_params = None
    if getattr(args, 'params_json', None) and os.path.isfile(args.params_json):
        try:
            with open(args.params_json, encoding='utf-8') as _pf:
                chosen_params = json.load(_pf)
        except Exception:
            chosen_params = None

    # Résolution des entrées
    stats = args.stats
    go    = args.go
    wgcna = args.wgcna
    umap  = args.umap
    if not stats:
        det   = detect_files(args.uploads)
        stats = stats or det['stats']
        go    = go    or det['go']
        wgcna = wgcna or det['wgcna']
        umap  = umap  or det['umap']

    if not stats or not os.path.isfile(stats):
        sys.exit("ERREUR : fichier statistique (rapportstatistique.xlsx) introuvable. "
                 "Utilisez --stats ou placez-le dans --uploads.")

    print("Fichiers :")
    print("  stats :", stats)
    print("  GO    :", go    or "(absent — module GO désactivé)")
    print("  WGCNA :", wgcna or "(absent — module WGCNA désactivé)")
    print("  UMAP  :", umap  or "(absent — module UMAP désactivé)")

    # Validation des fichiers d'entrée
    _warns, _errs = validate_inputs(stats, go, wgcna, umap)
    for w in _warns:
        print(f"  [WARN]   {w}")
    for e in _errs:
        print(f"  [ERREUR] {e}")
    if _errs:
        sys.exit(f"ERREUR : {len(_errs)} problème(s) bloquant(s) — corrigez les fichiers avant de relancer.")

    # 1. Lecture des stats
    df_comp, df_z, df_anova, df_inter, contrasts, sample_cols = read_stats(stats)
    conditions, cond_map = build_conditions(sample_cols)
    cond_map_inv = {col: cond for cond, cols in cond_map.items() for col in cols}
    print("  -> %d contrastes, %d échantillons, %d conditions, %d protéines"
          % (len(contrasts), len(sample_cols), len(conditions), len(df_comp)))

    # 2. Construction des données (lecture pure pour les stats)
    data = {
        'contrasts': contrasts,
        'sample_cols': sample_cols,
        'conditions': conditions,
        'cond_map': cond_map,
        'volcano': build_volcano(df_comp, contrasts),
        'stats': build_stats(df_comp, df_anova, contrasts),
        'heatmap': build_heatmap(df_z, df_anova, sample_cols),
        'top20_pi': build_top20_pi(df_comp, contrasts),
    }
    inter_rows, upset = build_intersections(df_inter, contrasts)
    data['intersections'] = inter_rows
    data['upset'] = upset
    # Projections de visualisation (dérivées des intensités, pas de recalcul de stats)
    data['pca'] = build_pca(df_comp, sample_cols, cond_map_inv)
    data['proteome_data'] = build_proteome_data(df_comp, contrasts, cond_map)
    # Modules optionnels
    data['go_data'], data['chord'] = build_go(go, contrasts, df_comp)
    data['go_images'] = build_go_images(go, contrasts)
    data['go_bubble']  = build_go_bubble(go, contrasts)
    data['go_heatmap']  = build_go_heatmap(go, contrasts)
    data['go_network']  = build_go_network(go, contrasts)
    data['wgcna_rows'], data['wgcna_top5'], data['wgcna_mods'] = build_wgcna(wgcna)
    data['wgcna_hub_map'] = build_wgcna_hub_map(wgcna)
    data['module_trait'] = build_module_trait(wgcna, df_comp, sample_cols, cond_map_inv)
    data['wgcna_net'] = build_wgcna_network(wgcna)
    # QC figures
    data['umap']        = build_umap(umap)
    data['violin']      = build_violin(df_comp, sample_cols, cond_map_inv)
    data['missing']     = build_missing_values(stats, sample_cols, cond_map_inv)
    # Nouvelles figures analytiques
    # Nouvelles figures analytiques (NS_REP_MAX=150 → allégé vs 400)
    data['rep_scatter'] = build_replicate_scatter(df_comp, cond_map_inv)
    data['ranked']  = build_ranked_abundance(df_comp, sample_cols)
    data['cv_plot'] = build_cv_plot(df_comp, sample_cols, cond_map_inv)
    data['qc_report'] = build_qc_report(stats, df_comp, sample_cols, cond_map_inv)
    data['about']   = build_about(stats, go, wgcna, umap,
                          contrasts, sample_cols, conditions,
                          len(df_comp), P_THRESH, LFC_THRESH, df_comp,
                          chosen_params=chosen_params)

    print("  -> proteome:%d contrastes | heatmap:%d prot | GO:%d | chord:%d | images:%d | wgcna:%d prot"
          % (len(data['proteome_data']), len(data['heatmap']), len(data['go_data']),
             len(data['chord']), len(data['go_images']), len(data['wgcna_rows'])))

    # 3. Injection dans le template
    template_html = open(args.template, encoding='utf-8').read()
    # Marquer le vrai </body> avec un token unique AVANT d'embarquer les libs
    # (xlsx.js contient var Wm="</body></html>" qui pollue les rfind ultérieurs)
    BODY_MARKER = '<!--PROTEOGEN_BODY_END-->'
    template_html = template_html.replace('</body>', BODY_MARKER + '</body>', 1)
    # Embarquer Chart.js et xlsx en mode offline
    script_dir = os.path.dirname(os.path.abspath(__file__))
    template_html = embed_libs(template_html, script_dir)
    html = inject(template_html, data)

    # 3b. Patch JS heatmap (palette Blue-White-Red, couleurs conditions, légende Z-score)
    html = patch_heatmap_js(html)
    # 3c. Injection des figures QC (violin plot + missing values + MA plot + scatter réplicats)
    html = patch_qc_figures(html)
    # 3d. UpSet interactif (clic barre → tableau protéines filtré)
    html = patch_upset_js(html)
    # 3e. UMAP interactif (scatter + ellipses de confiance + tooltip)
    html = patch_umap_js(html)
    # 3f. WGCNA network (graphe étoile hub-centric)
    html = patch_wgcna_network_js(html)
    # 3f. GO bubble plot + table (JS + HTML)
    html = patch_go_bubble_js(html)
    html = patch_gopage_html(html)
    # 3g. Ranked abundance + CV plot + À propos (JS)
    html = patch_ranked_cv_about_js(html)
    # 3h. Injection HTML de l'onglet QC (canvas + sélecteurs + script d'init)
    html = patch_qc_tab_html(html)
    # 3i. Patches post-monkey-patch (nécessitent que showPage soit déjà patché)
    html = patch_module_trait_js(html)
    html = patch_volcano_comp_js(html)
    html = patch_go_network_js(html)
    html = patch_go_heatmap_js(html)

    # 4. Vérifications avant livraison (checklist du skill)
    assert 'const PROTEOME_DATA=' in html
    assert 'const VOLCANO=PROTEOME_DATA' in html
    assert '"padj":' in html or len(contrasts) <= 1
    assert '"pi":' in html or len(contrasts) == 0
    assert 'TOP20_PI=' in html
    assert 'drawHeatmap' in html
    assert 'renderUpsetDetail' in html
    assert 'UMAP_DATA=' in html
    assert 'QC_REPORT=' in html
    size_kb = len(html.encode('utf-8')) / 1024
    print("Taille HTML : %.0f KB" % size_kb)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    open(args.out, 'w', encoding='utf-8').write(html)
    print("OK ->", args.out)
    return args.out


if __name__ == '__main__':
    main()
