# -*- coding: utf-8 -*-
"""
dashboard_integration.py — Bridge between deimos.py and build_dashboardv7.py

The dashboard (build_dashboardv7.py) expects:
  - a 'Differential_Expression' sheet with columns *_p.val, *_diff, *_p.adj,
    Pi_Score_*, Robustness_Score_*
  - a 'Zscore_Heatmap' sheet whose sample columns are named
    'Condition_Replicate' (condition inferred via rsplit('_', 1)[0])
  - GO_* sheets (detected by name) — compatible with our GO export
  - a SEPARATE WGCNA file whose 1st sheet carries Gene_Name/Module_Color/HubScore

Our pipeline produces a single workbook in which:
  - sample columns are the raw labels (e.g. 'Sple.S01'), without
    condition information encoded in the name;
  - WGCNA is a sheet ('WGCNA_Results') rather than a separate file.

This module prepares, in a temporary directory, "dashboard-ready" files
WITHOUT modifying the main workbook or the dashboard, then triggers generation.
"""

import os
import re
import sys
import shutil
import tempfile
import importlib.util
import pandas as pd
import openpyxl


def _load_dashboard_module(dashboard_path: str):
    """Imports build_dashboardv7.py as a module."""
    spec = importlib.util.spec_from_file_location("build_dashboard", dashboard_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _rename_sample_columns(df: pd.DataFrame, label_to_condrep: dict) -> pd.DataFrame:
    """Renames sample columns from 'label' to 'Condition_Replicate'."""
    return df.rename(columns={k: v for k, v in label_to_condrep.items()
                              if k in df.columns})


def _build_label_mapping(design: pd.DataFrame) -> dict:
    """
    Builds the label → 'Condition_Replicate' mapping expected by the dashboard.
    E.g. label='Sple.S01', condition='Ctrl' → 'Ctrl_1'
    Replicates are numbered by order of appearance within each condition.
    """
    mapping = {}
    counters = {}
    for _, row in design.iterrows():
        label = str(row["label"])
        cond  = str(row["condition"])
        # Sanitize the condition so it doesn't contain stray '_' characters
        cond_clean = re.sub(r"_+", "-", cond)
        counters[cond_clean] = counters.get(cond_clean, 0) + 1
        mapping[label] = f"{cond_clean}_{counters[cond_clean]}"
    return mapping


def prepare_dashboard_files(main_xlsx: str, design: pd.DataFrame,
                            work_dir: str) -> dict:
    """
    Creates the dashboard-ready files in work_dir:
      - stats.xlsx  : copy of the main workbook with sample columns renamed
                      to Condition_Replicate (Differential_Expression + Zscore sheets)
                      — also carries the GO_* sheets and raw_data / UpSet / ANOVA.
      - wgcna.xlsx  : dedicated file whose 1st sheet = WGCNA_Results.
    Returns {'stats':..., 'go':..., 'wgcna':..., 'umap':...} (paths or None).
    """
    label_map = _build_label_mapping(design)

    # ---- 1. Load the main workbook ----
    xl = pd.ExcelFile(main_xlsx)
    sheets = xl.sheet_names

    stats_path = os.path.join(work_dir, "stats.xlsx")

    # Retrieve per-sample intensities from raw_data (all proteins).
    # The dashboard expects these columns INSIDE the Differential_Expression
    # sheet (PCA, scatter, violin, CV...). Our pipeline doesn't put them
    # there → we inject them here.
    intensity_df = None
    if "raw_data" in sheets:
        try:
            df_brut = pd.read_excel(main_xlsx, sheet_name="raw_data")
            lfq_cols = [c for c in df_brut.columns if "LFQ" in str(c)]
            # Build the 'LFQ.intensity.<label>' → 'Condition_Rep' mapping
            lfq_rename = {}
            for c in lfq_cols:
                # extract the raw label (after the LFQ.intensity. prefix)
                lbl = re.sub(r"^LFQ\.intensity\.", "", str(c))
                if lbl in label_map:
                    lfq_rename[c] = label_map[lbl]
            key_col = "name" if "name" in df_brut.columns else "Protein.Group"
            if lfq_rename and key_col in df_brut.columns:
                intensity_df = df_brut[[key_col] + list(lfq_rename.keys())].copy()
                intensity_df = intensity_df.rename(columns=lfq_rename)
                intensity_df = intensity_df.rename(columns={key_col: "name"})
                # IMPORTANT: the dashboard (PCA, scatter, CV, violin) expects
                # intensities on a log2 scale, as the original R pipeline did.
                # The LFQ values in raw_data are RAW (up to several
                # millions) → without log2, the PCA produces huge values
                # (pc1 ~ -6e6) and the scatter rendering freezes the browser.
                import numpy as _np
                val_cols = [c for c in intensity_df.columns if c != "name"]
                for c in val_cols:
                    vals = pd.to_numeric(intensity_df[c], errors="coerce")
                    # log2(x) for x>0, NaN otherwise (missing/zero values)
                    intensity_df[c] = _np.where(vals > 0, _np.log2(vals), _np.nan)
        except Exception as e:
            print(f"  [dashboard] raw_data intensities not injected ({e}).")

    # Genes/description metadata to enrich the UpSet tab (the dashboard
    # expects 'Genes' and 'First.Protein.Description' columns there).
    meta_map = None
    try:
        df_brut_meta = pd.read_excel(main_xlsx, sheet_name="raw_data")
        keep = [c for c in ["name", "Protein.Group", "Genes",
                            "First.Protein.Description"] if c in df_brut_meta.columns]
        if "Protein.Group" in keep:
            meta_map = df_brut_meta[keep].drop_duplicates("Protein.Group")
    except Exception:
        meta_map = None

    # Rewrite a workbook renaming sample columns in Differential_Expression
    # and Zscore_Heatmap. Other sheets are copied through unchanged.
    with pd.ExcelWriter(stats_path, engine="openpyxl") as writer:
        for s in sheets:
            try:
                df = pd.read_excel(main_xlsx, sheet_name=s)   # header=0 everywhere
            except Exception:
                continue

            if s == "Differential_Expression":
                df = _rename_sample_columns(df, label_map)
                # Inject per-sample intensities (merge on 'name')
                if intensity_df is not None and "name" in df.columns:
                    new_cols = [c for c in intensity_df.columns
                                if c == "name" or c not in df.columns]
                    df = df.merge(intensity_df[new_cols], on="name", how="left")
            elif s == "Zscore_Heatmap":
                df = _rename_sample_columns(df, label_map)
            elif s == "UpSet_Intersections":
                # Enrich with Genes + description if missing (the dashboard reads them)
                if meta_map is not None and "Protein.Group" in df.columns:
                    add = [c for c in ["Genes", "First.Protein.Description"]
                           if c not in df.columns and c in meta_map.columns]
                    if add:
                        df = df.merge(meta_map[["Protein.Group"] + add],
                                      on="Protein.Group", how="left")
                # The dashboard re-reads this sheet with header=1: row 1 must
                # physically be a title, row 2 the column names, then the
                # data. We write it manually, without pandas' auto header.
                ws_up = writer.book.create_sheet(s[:31])
                ws_up.append(["Significant intersections table"])
                ws_up.append(list(df.columns))               # actual header (row 2)
                for row in df.itertuples(index=False):
                    ws_up.append(list(row))
                continue

            # truncate sheet name to 31 chars (Excel limit)
            df.to_excel(writer, sheet_name=s[:31], index=False)

    # ---- 2. Detect module presence ----
    has_go    = any(s.startswith("GO_") for s in sheets)
    has_wgcna = "WGCNA_Results" in sheets or any("wgcna" in s.lower() for s in sheets)

    # ---- 3. Dedicated WGCNA file (1st sheet = results) ----
    wgcna_path = None
    if has_wgcna:
        wgcna_sheet = "WGCNA_Results" if "WGCNA_Results" in sheets else \
            next(s for s in sheets if "wgcna" in s.lower())
        try:
            df_wgcna = pd.read_excel(main_xlsx, sheet_name=wgcna_sheet)
            # The dashboard expects Gene_Name / Module_Color / HubScore on the 1st sheet
            if {"Gene_Name", "Module_Color", "HubScore"}.issubset(df_wgcna.columns):
                wgcna_path = os.path.join(work_dir, "wgcna.xlsx")
                with pd.ExcelWriter(wgcna_path, engine="openpyxl") as w:
                    df_wgcna.to_excel(w, sheet_name="WGCNA_Results", index=False)
            else:
                print("  [dashboard] WGCNA sheet missing expected columns — WGCNA disabled.")
        except Exception as e:
            print(f"  [dashboard] WGCNA read failed ({e}) — WGCNA disabled.")

    # ---- 4. GO: reuse the stats workbook (includes GO_* sheets) ----
    go_path = stats_path if has_go else None

    # ---- 5. UMAP: the pipeline now exports a "UMAP" sheet with the
    # per-sample coordinates. If present, pass it to the dashboard.
    umap_path = None
    try:
        xl_check = pd.ExcelFile(stats_path)
        if any(s.lower() == "umap" for s in xl_check.sheet_names):
            umap_path = stats_path
    except Exception:
        umap_path = None

    return {"stats": stats_path, "go": go_path,
            "wgcna": wgcna_path, "umap": umap_path}


def run_dashboard(main_xlsx: str, design: pd.DataFrame, out_html: str,
                  dashboard_script: str = None, template: str = None,
                  params: dict = None) -> str | None:
    """
    Generates the HTML dashboard from the pipeline's Excel workbook.
    Fully guarded: returns None and never interrupts the pipeline on
    failure.

    Parameters
    ----------
    main_xlsx        : path to the ProteomicAnalysis_Results.xlsx workbook
    design           : design DataFrame (columns label, condition)
    out_html         : output path of the HTML dashboard
    dashboard_script : path to build_dashboardv7.py (default: alongside this module)
    template         : path to dashboard_template.html (default: alongside the script)

    Returns
    -------
    Path to the generated HTML, or None on failure.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    dashboard_script = dashboard_script or os.path.join(here, "build_dashboardv7.py")
    template = template or os.path.join(here, "dashboard_template.html")

    if not os.path.isfile(dashboard_script):
        print(f"  [dashboard] script not found ({dashboard_script}) — dashboard skipped.")
        return None
    if not os.path.isfile(template):
        print(f"  [dashboard] template not found ({template}) — dashboard skipped.")
        return None

    try:
        work_dir = tempfile.mkdtemp(prefix="dashboard_")
        files = prepare_dashboard_files(main_xlsx, design, work_dir)

        # Write the chosen thresholds to a JSON file the dashboard will read
        # to display/apply the ACTUAL parameters (instead of inferring them).
        if params is not None:
            try:
                import json as _json
                pj = os.path.join(work_dir, "dashboard_params.json")
                with open(pj, "w", encoding="utf-8") as _f:
                    _json.dump({
                        "volcano_use_padj":  bool(params.get("volcano_use_padj", False)),
                        "volcano_p_thresh":  float(params.get("volcano_p_thresh", 0.05)),
                        "volcano_ratio_min": float(params.get("volcano_ratio_min", 1.5)),
                        "volcano_lfc_min":   float(params.get("volcano_lfc_min", 0.585)),
                        "anova_use_padj":    bool(params.get("anova_use_padj", False)),
                        "anova_p_thresh":    float(params.get("anova_p_thresh", 0.05)),
                        "n_heatmap_clusters": int(params.get("n_heatmap_clusters", 3)),
                    }, _f)
            except Exception:
                pass

        mod = _load_dashboard_module(dashboard_script)

        # Build an argparse.Namespace equivalent to the dashboard's CLI call
        import argparse
        args = argparse.Namespace(
            uploads=work_dir,
            stats=files["stats"],
            go=files["go"],
            wgcna=files["wgcna"],
            umap=files["umap"],
            template=template,
            out=out_html,
            params_json=os.path.join(work_dir, "dashboard_params.json")
                if params is not None else None,
        )

        # The dashboard's main() parses sys.argv; we work around this by
        # calling its logic directly if it accepts args, otherwise we
        # simulate sys.argv.
        result = _invoke_dashboard_main(mod, args)
        # Cleanup
        shutil.rmtree(work_dir, ignore_errors=True)
        return result

    except Exception as e:
        print(f"  [dashboard] generation failed ({type(e).__name__}: {str(e)[:100]}) "
              f"— pipeline unaffected.")
        return None


def _invoke_dashboard_main(mod, args):
    """
    Calls the dashboard's main(). build_dashboardv7.main() parses sys.argv,
    so we rebuild sys.argv for the duration of the call.
    """
    argv_backup = sys.argv[:]
    try:
        sys.argv = ["build_dashboardv7.py",
                    "--stats", args.stats,
                    "--template", args.template,
                    "--out", args.out]
        if args.go:
            sys.argv += ["--go", args.go]
        if args.wgcna:
            sys.argv += ["--wgcna", args.wgcna]
        if args.umap:
            sys.argv += ["--umap", args.umap]
        if getattr(args, "params_json", None):
            sys.argv += ["--params-json", args.params_json]
        mod.main()
        return args.out if os.path.isfile(args.out) else None
    finally:
        sys.argv = argv_backup
