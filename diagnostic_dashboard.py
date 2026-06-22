# -*- coding: utf-8 -*-
"""
diagnostic_dashboard.py — Checks what the dashboard will find in the Excel
workbook for UpSet, PCA and QC. Run from the folder containing proteomics_output/.

    python diagnostic_dashboard.py
"""
import os
import sys
import pandas as pd

XLSX = os.path.join("proteomics_output", "ProteomicAnalysis_Results.xlsx")
if len(sys.argv) > 1:
    XLSX = sys.argv[1]

print(f"=== Dashboard diagnostic on: {XLSX} ===\n")
if not os.path.exists(XLSX):
    print("[ERROR] File not found. Run this script from the pipeline's folder.")
    sys.exit(1)

xl = pd.ExcelFile(XLSX)
print("Sheets present:", xl.sheet_names, "\n")

# ---- 1. Zscore_Heatmap (main source of sample_cols) ----
META = {"name", "Protein.Group", "Protein.Names", "Genes",
        "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences",
        "Cluster_ID", "Peptide_Count_min"}
print("--- 1. Zscore_Heatmap (-> sample_cols for PCA/QC) ---")
try:
    dfz = pd.read_excel(XLSX, sheet_name="Zscore_Heatmap")
    num_cols = [c for c in dfz.columns
                if c not in META and pd.api.types.is_numeric_dtype(dfz[c])]
    print(f"  Rows: {len(dfz)}")
    print(f"  Numeric columns (= sample_cols): {len(num_cols)}")
    print(f"  Examples: {num_cols[:6]}")
    if len(dfz) == 0:
        print("  [WARN] Zscore_Heatmap is EMPTY (0 significant ANOVA protein) -> PCA/QC will")
        print("         fall back to the Differential_Expression sheet (injected intensities).")
    if num_cols:
        # Check the Condition_Replicate format
        bad = [c for c in num_cols if "_" not in str(c)]
        if bad:
            print(f"  [WARN] Columns WITHOUT underscore (condition not derivable): {bad[:5]}")
        else:
            conds = sorted(set(str(c).rsplit("_", 1)[0] for c in num_cols))
            print(f"  [OK] Conditions inferred: {conds}")
except Exception as e:
    print(f"  [ERROR] {e}")

# ---- 2. Differential_Expression (fallback sample_cols + scatter + PCA) ----
print("\n--- 2. Differential_Expression (intensities injected?) ---")
try:
    dfc = pd.read_excel(XLSX, sheet_name="Differential_Expression")
    stat_suf = ("_p.val", "_p.adj", "_diff", "_DEqMS")
    stat_pre = ("Pi_Score_", "Robustness_Score_")
    intens = [c for c in dfc.columns if c not in META
              and not any(str(c).endswith(s) for s in stat_suf)
              and not any(str(c).startswith(p) for p in stat_pre)
              and pd.api.types.is_numeric_dtype(dfc[c])]
    print(f"  Per-sample intensity columns detected: {len(intens)}")
    print(f"  Examples: {intens[:6]}")
    if len(intens) < 2:
        print("  [WARN] NO per-sample intensities in Differential_Expression.")
        print("         -> PCA impossible if Zscore is also empty. (raw_data injection failed?)")
    else:
        print("  [OK] Intensities present -> PCA/scatter possible.")
except Exception as e:
    print(f"  [ERROR] {e}")

# ---- 3. UpSet_Intersections ----
print("\n--- 3. UpSet_Intersections ---")
try:
    # test both reads (header=0 and header=1)
    for h in (0, 1):
        d = pd.read_excel(XLSX, sheet_name="UpSet_Intersections", header=h)
        cols = [c for c in d.columns if not str(c).startswith("Unnamed")]
        has_nb = "Nb_Conditions" in d.columns
        has_pg = "Protein.Group" in d.columns
        print(f"  header={h}: {len(d)} rows | Nb_Conditions={has_nb} | "
              f"Protein.Group={has_pg} | first columns={list(d.columns)[:4]}")
    print("  (The dashboard reads this sheet with header=1: that row must")
    print("   show Nb_Conditions=True and Protein.Group=True.)")
except Exception as e:
    print(f"  [ERROR] {e}")

print("\n=== End of diagnostic ===")
print("Copy-paste everything above for diagnostics.")
