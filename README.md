# OmniSeq

[![R-CMD-check](https://github.com/nizarberrimgeno-prog/OmniSeq/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nizarberrimgeno-prog/OmniSeq/actions/workflows/R-CMD-check.yaml)

Multi-omics RNA-seq analysis with transposable element (TE) integration.

OmniSeq wraps DESeq2, clusterProfiler, GSVA, WGCNA, mixOmics, decoupleR, and
immunedeconv into a single coherent interface built around one S3 object.
It takes a raw count matrix and sample metadata through QC, differential
expression, pathway enrichment, network/TF analysis, and multi-omics
integration — with built-in support for transposable element quantification
and innate immune sensor scoring.

## Installation

OmniSeq is not on CRAN/Bioconductor yet. Install from source:

```r
# from the package directory
devtools::install()
```

Most of OmniSeq's functionality comes from optional dependencies (Bioconductor
and GitHub packages) that aren't required just to load the package. After
installing OmniSeq, run:

```r
library(OmniSeq)
omni_check_deps()    # see what's installed / missing
omni_install_deps()  # install everything (CRAN + Bioconductor + GitHub)
```

You can install a subset instead of everything:

```r
omni_install_deps(modules = c("core", "enrichment"))
```

## Quickstart

```r
library(OmniSeq)

# 1. Build the object from a raw count matrix + sample metadata
obj <- init_omniseq(
  count_matrix  = my_counts,      # features (rows) x samples (cols), raw integers
  metadata      = my_meta,        # data.frame, rownames = sample IDs
  condition_col = "Genotype",     # main biological condition
  covariates    = "Batch",        # optional covariates in the model
  species       = "Homo sapiens",
  org_db        = "org.Hs.eg.db"
)

# 2. QC
obj <- run_omni_qc(obj)
plot_omni_qc_basics(obj, plot_type = "library_size")
plot_omni_qc_advanced(obj, plot_type = "pca_compare")

# 3. Differential expression (DESeq2 under the hood, apeglm shrinkage by default)
obj <- run_omni_de(obj)
res <- omni_results(obj)
head(as.data.frame(res))

plot_omni_volcano(obj, lfc_thresh = 1, top_n = 15)
plot_omni_pca(obj)
plot_omni_heatmap(obj, ntop = 50)

# 4. Enrichment
obj <- run_omni_enrichment(obj)   # ORA
obj <- run_omni_gsea(obj)         # GSEA

# 5. Export everything (plots + tables) to disk
export_omniseq(obj, outdir = "./omniseq_output")
```

`print(obj)` / `summary(obj)` give a quick status readout at any point —
what's been run, which contrasts are cached, and which is active.

### Full pipeline in one call

For a standard analysis with no manual steps, `run_omniseq_auto()` runs QC,
DE, visualizations, enrichment (ORA/GSEA/ssGSEA), TF activity, deconvolution,
STRING network, and export in sequence — skipping (not failing on) any step
whose dependency isn't installed:

```r
obj <- run_omniseq_auto(
  count_matrix  = my_counts,
  metadata      = my_meta,
  condition_col = "Genotype",
  outdir        = "./omniseq_output"
)
```

`run_omniseq_guided()` runs the same pipeline with interactive menu prompts
(design setup, coefficient selection, etc.) for exploratory console use.

## Interactive app

Every module is also available through a point-and-click Shiny app:

```r
OmniSeq::omni_launch_app()
```

| Tab | What it runs |
|---|---|
| Data | Upload counts + metadata (CSV/TSV), a simulated demo dataset, or reload a saved object; input validation; design setup |
| QC | `run_omni_qc()` and the QC plots |
| Differential expression | `run_omni_de()`, custom contrasts (`run_omni_contrast()`), signature masking (`run_omni_mask()`); volcano, MA, heatmap, PCA, top features, results table |
| Enrichment | ORA, GSEA, ssGSEA, cross-contrast comparison |
| Biology | TF activity, STRING network, immune deconvolution |
| Multivariate | WGCNA (incl. per-group networks + module preservation), sPLS-DA |
| TE analysis | TE family analysis, TE plots, proximity + GREAT |
| Multi-omics | Pair a gene and a TE object: sensor scores, TE-gene correlation, WGCNA TE traits, DIABLO, mediation |
| Export | ZIP of all tables/figures, HTML report, `.rds` object, session log, dependency status |

The active contrast is chosen once in the sidebar and drives every tab.
Tabs whose optional packages are missing stay visible and list what to
install. The app contains no analysis code of its own — it only calls the
exported package functions.

## The OmniSeq object

Everything hangs off one S3 object with a fixed structure:

| Slot | Contents |
|---|---|
| `$data` | The `SummarizedExperiment` (raw counts + metadata) |
| `$params` | Condition, covariates, design formula, species, data type |
| `$model_info` | Fitted `dds`, cached VST matrix, active contrast |
| `$results_cache` | DE results per contrast, keyed by contrast name |
| `$qc`, `$enrichment`, `$advanced_gsea`, ... | Populated as you run each step |

Use the accessor functions (`omni_results()`, `omni_model()`, `omni_vst()`,
`omni_counts()`, `omni_metadata()`, `omni_contrasts()`, `omni_qc()`) rather
than reaching into the object directly — they're what insulate your analysis
scripts from any internal slot changes.

## Built-in gene signatures

A handful of curated signatures ship with the package for post-DE masking
(e.g. checking whether a signal survives once interferon response genes are
excluded):

```r
omni_list_signatures()
obj <- run_omni_mask(obj, signature = "isgs_28")
```

Available: `omniseq_isgs_28`, `omniseq_isgs_extended`, `omniseq_ribosomal`,
`omniseq_mitochondrial`, `omniseq_immunoglobulin`.

## Transposable element analysis

If `init_omniseq()` detects TE locus IDs (format `<number>:<Class>:<Family>`,
e.g. `1000014:ERV:MSTB`) instead of gene symbols, the object is tagged
`data_type = "TE"` and TE-specific functions (`run_omni_te_family()`,
`run_omni_te_great()`, `run_omni_sensor_score()`, TE/gene correlation plots)
become relevant. A TE object and a gene object from the same samples can be
cross-referenced — e.g. `export_omniseq(te_obj, gene_obj = gene_obj)` for
TE-gene correlation plots, or `run_omni_diablo()` for multi-omics
integration.

## Status

Early-stage (v0.1.0), no CRAN/Bioconductor release yet, no automated test
suite. Treat exported functions as stable enough for active analysis work,
but expect breaking changes before a 1.0 release.

## License

MIT — see [LICENSE](LICENSE).
