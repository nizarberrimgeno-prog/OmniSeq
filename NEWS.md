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
