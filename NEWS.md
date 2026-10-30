# OmniSeq (development version)

## Shiny app

* `omni_launch_app()` now covers the full workflow instead of QC + DE only:
  custom contrasts, signature masking, enrichment (ORA / GSEA / ssGSEA /
  cross-contrast), TF activity, STRING, deconvolution, WGCNA and module
  preservation, sPLS-DA, TE analysis, gene + TE multi-omics (sensor scores,
  DIABLO, mediation), and export (ZIP, HTML report, `.rds`).
* The app is split into one Shiny module per tab (`inst/shiny-app/R/`),
  shares a single app-wide active contrast, captures package messages into
  a session log, and offers a simulated demo dataset.
* Fixed: the app now attaches OmniSeq itself, so launching it via
  `OmniSeq::omni_launch_app()` without `library(OmniSeq)` no longer fails.
* New `omni_set_contrast()` switches the active contrast without
  recomputing anything.
* App support for the new reporting functions: interactive volcano
  (`plot_omni_volcano_interactive()`), searchable results table
  (`omni_gene_search_table()`), an "ORA Up/Down" panel
  (`run_omni_ora_updown()`, `plot_omni_enrichment_interactive()`,
  `plot_omni_enrichment_comparison()`), and a "Compare contrasts" tab
  (`plot_omni_lfc_scatter_interactive()`, `plot_omni_interaction_scatter()`,
  `plot_omni_interaction_scatter_static()`, `plot_omni_interaction_dumbbell()`).
  The heatmap and top-features plots now honour the app's |log2FC| threshold.
* Sensor scoring in the app accepts custom TE family lists, and lists the
  most abundant families in the paired TE data, for annotations whose family
  names differ from RepeatMasker's (e.g. Telescope's `HML2` for HERVK).
* The app still starts when plotly is not installed; interactive panels
  then show which package to install.

## Bug fixes (found by running every exported function on public data)

See `dev/real-data-tests/` for the scripts, datasets and results.

* `omni_build_te_coords()` failed with `could not find function "%>%"`
  unless dplyr/magrittr was attached; it no longer uses the pipe.
* `export_omniseq_report()` failed to render whenever `outdir` was a
  relative path (including the `"./omniseq_output"` default) and still
  returned the HTML path. The path is now made absolute, and a render
  failure gives a warning and returns `NULL`.
* `run_omni_string()` / `plot_omni_string_network()`: genes that STRING maps
  to several proteins (e.g. CT45A paralogs) produced duplicate vertices and
  a "Duplicate vertex names" error. Nodes are now one per gene and edges one
  per gene pair; objects saved before the fix also plot.
* `plot_omni_string_network()`, `plot_omni_tf_network()` and
  `plot_omni_wgcna_hub_network()` errored on PDF devices ("invalid font
  type") because `ggraph::theme_graph()` defaults to Arial Narrow; they now
  use the standard sans family.
* `run_omni_sensor_score()` now keeps every metadata column in the score
  table. `plot_omni_sensor_irgs_corr()` previously failed when the grouping
  column was not one of a fixed set of names, and its regression line is
  now a single fit per panel (matching the reported R and p) instead of one
  per point shape.
* `plot_omni_sc_profile()` always errored after drawing (it rounded a data
  frame containing the genotype labels).
* `plot_omni_sc_genotype()` / `plot_omni_sc_two_layer()` labelled the
  reference group `NA` when the genotype column was a factor, and matched
  Wilcoxon groups by metadata row order rather than by sample name.
* `run_omni_sc_deconv(method = "bisque")` always used BisqueRNA's
  same-individuals mode, which fails with an independent reference, and
  then returned the object unchanged. Overlap mode is now used only when
  bulk and reference samples overlap, and an explicit `"bisque"` run that
  produces nothing is an error.
* `plot_omni_deconv_heatmap()` (default `aggregate_main = TRUE`) dropped every
  cell type outside its immune map, so a tissue reference (e.g. islet
  alpha/beta cells via MuSiC) left too few rows and pheatmap failed to
  cluster. Unmapped cell types are now kept as they are.
* `plot_omni_volcano_interactive()` and `plot_omni_lfc_scatter_interactive()`
  ignored their colour schemes (plotly default palette); points now use the
  documented Up/Down and category colours.
* msigdbr calls use the current `collection` / `subcollection` arguments
  throughout (the minimum version is now msigdbr 10.0.0, which the existing
  `collection =` calls already required).
* `plot_omni_tf_genotype()` explains that `geno_contrasts` must be a named
  list instead of failing with "Could not extract per-genotype TF scores".
* Missing-package errors for MuSiC, immunedeconv, xCell and BisqueRNA now
  give their GitHub install command instead of `BiocManager::install()`.
  `omni_install_deps()` installs immunedeconv, MuSiC and BisqueRNA from
  GitHub, adds the interactive/app packages, and no longer reinstalls
  GitHub packages on every run.
* Packages used via `::` but missing from DESCRIPTION (tidyr, Biobase,
  BisqueRNA, BiocManager, remotes) are now listed in Suggests.

# OmniSeq 0.1.0

Initial release.

## Core object & workflow

* `init_omniseq()` builds a validated `OmniSeq` S3 object from a raw count
  matrix and sample metadata, auto-detecting gene vs. transposable element
  (TE) data from row name format.
* `params` / `model_info` / `results_cache` kept as separate, documented
  slots rather than one flat list; accessor functions (`omni_results()`,
  `omni_model()`, `omni_vst()`, `omni_counts()`, `omni_metadata()`,
  `omni_contrasts()`, `omni_qc()`) insulate downstream code from internal
  slot changes.
* VST is computed once per object and cached in `model_info$vst` — never
  recomputed by downstream QC, plotting, or enrichment functions.
* `interactive()` guards on every `menu()`/`readline()` call so all
  functions run cleanly in non-interactive/scripted sessions.

## Differential expression

* `run_omni_de()` / `run_omni_contrast()` support `apeglm` and `ashr` LFC
  shrinkage, with automatic fallback when a package is unavailable.
* The DESeq2 `stat` column is preserved through shrinkage on both engines,
  so `run_omni_tf()`'s Wald-statistic-based ULM input stays correct instead
  of silently falling back to log2FC.
* Rebuilt models (when a contrast uses a different metadata column) are
  cached per column in `model_info$cached_column_dds`, so multiple
  contrasts sharing a column reuse one fit instead of re-running `DESeq()`
  for each.
* `run_omni_mask()` masks a gene signature (built-in or custom) post-hoc
  without disturbing the original cached contrast.

## Enrichment, TF activity & networks

* `run_omni_enrichment()` (ORA), `run_omni_gsea()`, and `run_omni_ssgsea()`
  cover over-representation, gene set enrichment, and per-sample pathway
  scoring against MSigDB.
* `run_omni_tf()` infers transcription factor activity via `decoupleR` in
  two modes: contrast-based (Wald t-statistics, one score per TF per
  contrast) or sample-based (VST-restricted to one condition group, never
  mixing groups).
* Regulons are fetched once and cached in `params$tf_regulons`, avoiding
  repeat downloads across contrasts.
* `run_omni_string()` supports a persistent local cache directory to avoid
  re-downloading STRING database files every session.
* `run_omni_wgcna()` / `run_omni_wgcna_multi()` fix the `stats::cor` /
  `WGCNA::cor` masking conflict that otherwise crashes `blockwiseModules()`.

## Transposable element analysis

* TE locus IDs (`<id>:<Class>:<Family>`) are recognized automatically; TE
  identifiers are always preserved through the protein-coding gene filter.
* `run_omni_te_great()` runs GREAT proximity enrichment locally via
  `rGREAT::great()` (rGREAT >= 1.99.0), replacing the deprecated web-server
  API.
* `run_omni_sensor_score()` scores innate immune sensor activation
  (cGAS-STING, RIG-I/MDA5, TLR7/TLR8) from TE family expression.
* Built-in gene signatures for post-DE masking: `omniseq_isgs_28`,
  `omniseq_isgs_extended`, `omniseq_isgs_schoggins`, `omniseq_ribosomal`,
  `omniseq_mitochondrial`, `omniseq_immunoglobulin`.

## Multi-omics integration

* `run_omni_splsda()` / `run_omni_diablo()` provide sPLS-DA and DIABLO
  (gene + TE block) integration via mixOmics, with automatic component and
  `keepX` tuning through leave-one-out cross-validation.
* `run_omni_mediation()` tests causal mediation (e.g. mutation type → TE
  sensor activation → interferon response) with bootstrapped confidence
  intervals.
* `run_omni_sc_deconv()` auto-branches between MuSiC (raw scRNA-seq
  reference) and signature scoring / BisqueRNA (pseudo-bulk reference) for
  cell type deconvolution.

## Export

* `export_omniseq()` and `export_omniseq_report()` save all tables, plots,
  and a self-contained HTML report to a structured output directory, with
  masked-contrast outputs kept distinct from unmasked ones by filename.
