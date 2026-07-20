# Changelog — Deimos

All notable changes to the Deimos proteomics pipeline.
Repository: https://github.com/BenoitBe/deimos-proteomics-pipeline

---

## [Unreleased] — Major update

This release bundles new statistical capabilities, several correctness fixes,
and dashboard improvements accumulated since the previous public version.

### Added

- **Paired design support (repeated measures).** When the experimental design
  contains a `subject` column (aliases: `pair`, `patient`, `individual`), Deimos
  automatically switches to a paired model `~0 + condition + subject`. Each
  subject gets its own fixed-effect intercept, absorbing inter-individual
  baseline variability (generalized paired t-test). Triggered only when a subject
  spans ≥2 conditions; otherwise the standard `~0 + condition` model is used.
  Subject IDs are normalized (trailing spaces and `.0` float suffixes removed),
  and a warning is emitted for empty `subject` values (e.g. CSV separator typos).
  `check_paired.py` lets you verify the design *before* running.

- **Rank-based GSEA (`gsea_enrichment.py`).** Complements the existing ORA. Ranks
  ALL proteins by `-log10(p) × sign(LFC)` and detects coordinated pathway shifts
  even below the DEP threshold, using `gseapy.prerank` (Broad algorithm). Gene
  sets from Enrichr (GO / KEGG / Reactome). Optional (asked at launch), fully
  non-blocking: the pipeline runs normally if `gseapy` is absent or offline.
  Output: `GSEA_<contrast>` sheets + NES bar plots.

- **GO enrichment per ANOVA heatmap cluster.** Each co-expression cluster from
  the ANOVA heatmap is functionally annotated (ORA), exported as
  `Cluster_Enrich_N` sheets (named to avoid the dashboard's GO filters).

- **RLE plot in QC.** Relative Log Expression plot added to the QC sheet, to
  document residual normalization quality after DIA-NN (diagnostic only, no
  re-normalization). Title reports the worst per-sample median deviation.

- **Utility scripts:** `check_paired.py` (pre-flight paired-design check),
  `check_gsea.py` (GSEA integration diagnostic), `regen_dashboard.py` (rebuild
  the HTML dashboard from an existing Excel workbook without rerunning analysis).

### Fixed

- **Sample-to-column matching (critical).** Replaced substring matching with
  exact matching. Previously, numeric-prefix labels collided (`Sple-2` matched
  `Sple-28`/`Sple-29`), which could silently assign a sample to the wrong run.
  **Users should re-check past analyses that had prefix-related sample labels.**

- **Robustness score under DEqMS.** The robustness loop now recomputes DEqMS
  (not limma) when DEqMS is the primary statistic, so the score reflects the same
  statistic shown in the `_p.val` / `_p.adj` columns.

- **Scatter plots (dashboard).** Empty/deformed scatter fixed: intensity columns
  were not matched to contrasts for digit-prefixed conditions (missing `X`
  prefix). Axes are now square/symmetric, colors match the Excel scatter
  (Up=orange, Down=cyan), non-significant points drawn underneath, ±ratio lines
  added.

- **GO cross-contrast heatmap & network (dashboard).** (1) Empty rendering fixed:
  the term canvas exceeded the browser height limit (~32767 px) with thousands of
  terms — now capped at 150 terms (variance-ranked, per-source quota) and 150
  network nodes. (2) `REAC` (Reactome) and `KEGG` added to the source selectors
  and to the per-source panels. (3) Contrast→sheet mapping made deterministic
  (was confusing contrasts sharing a long common prefix).

- **Excel workbook corruption.** GO sheet names sharing a long common prefix
  truncated to identical 31-char names, corrupting `/xl/workbook.xml`. Names are
  now shortened (common affix stripped) and guaranteed unique.

- **UpSet plot labels.** Hard-coded `X180/X180` labels replaced by a generic
  common-prefix stripper (e.g. `Ctrl/di6h`).

- **GSEA parameter not applied.** The `run_gsea` flag chosen at launch was
  dropped in config resolution and never reached the pipeline; now transmitted
  end-to-end and persisted in `last_config.yaml`.

- **Duplicate prompt at launch.** The p-value question could appear twice due to
  a residual empty line in the input buffer (common on PowerShell). Empty inputs
  are now absorbed silently.

- **WGCNA `grey` module.** The unassigned-protein `grey` module is now shown in
  the module-size plot (hatched, labeled "unassigned", placed last) instead of
  being silently dropped — its size is a clustering-quality indicator.

### Changed

- **Full FR → EN translation** of the user-facing pipeline (prompts, progress
  messages, figure titles, Excel sheet names and headers, dashboard UI).
- **Dashboard header** homogenized with the sister suite (Phobos): icon +
  "Deimos" + expanded acronym + "Proteogen · Université de Caen".

### Dependencies

- **New:** `gseapy` (rank-based GSEA). See `requirements.txt`.
