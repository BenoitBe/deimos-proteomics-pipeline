# -*- coding: utf-8 -*-
"""
rebuild_dashboard.py — Regenerates proteogen_dashboard.html from a finished Deimos run,
without re-running the pipeline (e.g. after a dashboard update).

    python rebuild_dashboard.py                      # ./proteomics_output
    python rebuild_dashboard.py path/to/proteomics_output [--out other.html]

Thresholds come from the run's last_config.yaml (same rule as the Excel); the samples
and conditions from the PCA_Scores sheet (the samples actually analysed), else from the
design CSV named in the config.
"""
import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from config import _DEFAULTS, _resolve, load_config   # noqa: E402
from dashboard_integration import run_dashboard        # noqa: E402


def _design(xlsx, params, out_dir):
    try:
        sc = pd.read_excel(xlsx, sheet_name="PCA_Scores")
        if {"label", "condition"} <= set(sc.columns):
            return sc[["label", "condition"]].astype(str), "PCA_Scores sheet"
    except Exception:                                                 # noqa: BLE001
        pass
    for cand in (params.get("design_path"), os.path.join(out_dir, "..", "ExperimentalDesign.csv")):
        if cand and os.path.exists(cand):
            d = pd.read_csv(cand, sep=None, engine="python")
            return d[["label", "condition"]].astype(str), cand
    raise SystemExit("[ERROR] No design found (PCA_Scores sheet or design CSV).")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("out_dir", nargs="?", default="proteomics_output")
    ap.add_argument("--xlsx", help="results workbook (default: ProteomicAnalysis_Results.xlsx)")
    ap.add_argument("--out", help="output HTML (default: <out_dir>/proteogen_dashboard.html)")
    ap.add_argument("--config", help="YAML of the run, if last_config.yaml is missing (runs "
                    "launched with --config before 2026-11-06 did not write it)")
    a = ap.parse_args(argv)
    xlsx = a.xlsx or os.path.join(a.out_dir, "ProteomicAnalysis_Results.xlsx")
    if not os.path.exists(xlsx):
        raise SystemExit(f"[ERROR] Workbook not found: {xlsx}")
    cfg = a.config or os.path.join(a.out_dir, "last_config.yaml")
    if os.path.exists(cfg):
        params, src = load_config(cfg), cfg
    else:
        params, src = _resolve(dict(_DEFAULTS)), "Deimos defaults (no last_config.yaml)"
    design, dsrc = _design(xlsx, params, a.out_dir)
    print(f"[dashboard] workbook : {xlsx}\n[dashboard] thresholds: {src}\n"
          f"[dashboard] samples  : {len(design)} ({dsrc})")
    out = run_dashboard(xlsx, design, a.out or os.path.join(a.out_dir, "proteogen_dashboard.html"),
                        params=params)
    print(f"[OK] {out}" if out else "[WARN] Dashboard not generated (see messages above).")
    return 0 if out else 1


if __name__ == "__main__":
    sys.exit(main())
