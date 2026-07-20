# -*- coding: utf-8 -*-
"""
check_paired.py — Vérifie AVANT de lancer Deimos si ton design sera
détecté comme apparié (paired design).

Usage : python check_paired.py ExperimentalDesign.csv
"""
import sys
import pandas as pd

if len(sys.argv) < 2:
    print("Usage: python check_paired.py <ExperimentalDesign.csv>")
    sys.exit(1)

path = sys.argv[1]
# Le design Deimos est séparé par ';'
try:
    design = pd.read_csv(path, sep=";")
except Exception as e:
    print(f"[ERROR] Cannot read {path}: {e}")
    sys.exit(1)

print("=" * 60)
print("  CHECK PAIRED DESIGN")
print("=" * 60)
print(f"\nColumns found: {list(design.columns)}")
print(f"Samples: {len(design)}")

# 1. Colonne subject présente ?
subj_col = None
for cand in ("subject", "pair", "patient", "individual"):
    if cand in design.columns:
        subj_col = cand
        break

if subj_col is None:
    print("\n[FAIL] No 'subject'/'pair'/'patient'/'individual' column found.")
    print("       -> Deimos will use the STANDARD (non-paired) model.")
    print("       Fix: add a column named exactly 'subject' to your CSV.")
    sys.exit(0)

print(f"\n[OK] Subject column detected: '{subj_col}'")

if "condition" not in design.columns:
    print("[FAIL] No 'condition' column — cannot check pairing.")
    sys.exit(1)

subjects = design[subj_col].astype(str)
conditions = design["condition"].astype(str)

# Normalisation IDENTIQUE à Deimos : retire espaces + suffixe '.0'
subjects_norm = subjects.str.strip().str.replace(r"\.0$", "", regex=True)
n_before = subjects.str.strip().nunique()
n_after = subjects_norm.nunique()
if n_after < n_before:
    print(f"\n[INFO] After normalization (spaces + '.0' removed): "
          f"{n_before} -> {n_after} unique subjects.")
    print("       Deimos applies this automatically — pairing will be correct.")
subjects = subjects_norm

# 2. Cohérence d'annotation résiduelle
stripped = subjects.str.strip()
if (stripped != subjects).any():
    print("[WARN] Some subject values have leading/trailing spaces — "
          "they will be treated as DIFFERENT subjects!")
lower_map = subjects.groupby(subjects.str.lower().str.strip()).nunique()
ambiguous = lower_map[lower_map > 1]
if len(ambiguous):
    print(f"[WARN] Case/spacing variants detected (seen as different subjects): "
          f"{list(ambiguous.index)}")

# 3. Conditions de déclenchement de Deimos
pairing = pd.DataFrame({"s": subjects, "c": conditions})
n_cond_per_subj = pairing.groupby("s")["c"].nunique()
n_subj = pairing["s"].nunique()
n_samples = len(subjects)

cond1 = (n_cond_per_subj >= 2).any()
cond2 = n_subj >= 2
cond3 = n_subj < n_samples

print(f"\n  Subjects: {n_subj}")
print(f"  Conditions covered per subject:")
for s, k in n_cond_per_subj.items():
    flag = "  <-- covers >=2 conditions" if k >= 2 else ""
    print(f"      {s}: {k} condition(s){flag}")

print(f"\n  Trigger conditions:")
print(f"    [{'OK' if cond1 else 'FAIL'}] at least one subject in >=2 conditions")
print(f"    [{'OK' if cond2 else 'FAIL'}] >=2 subjects total ({n_subj})")
print(f"    [{'OK' if cond3 else 'FAIL'}] fewer subjects than samples "
      f"({n_subj} < {n_samples})")

print()
if cond1 and cond2 and cond3:
    n_terms = n_subj - 1
    print(f"  ==> PAIRED design WILL be used.")
    print(f"      Model: ~0 + condition + subject ({n_terms} subject terms absorbed)")
    print(f"      You should see this line in the console at runtime:")
    print(f"      [MODEL] PAIRED design (~0 + condition + subject): "
          f"{n_subj} subjects, {n_terms} subject terms absorbed.")
else:
    print(f"  ==> STANDARD design will be used (NOT paired).")
    print(f"      Check the failed condition(s) above.")
print("=" * 60)
