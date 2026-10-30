# Real-data tests of every exported function

These scripts run every exported OmniSeq function on public RNA-seq data
(bulk gene + TE counts from the same samples, and bulk + single-cell data
from matching tissue). They are developer tooling: excluded from the package
build (`.Rbuildignore`) and too slow/heavy for CI.

## Datasets

| Dataset | What it is | Used for |
|---|---|---|
| **GSE182119** gene counts | 4 ovarian cancer lines (A2780, Hey, Kuramochi, TykNu) × Mock / ITF-2357 (HDACi) / 5-azacytidine (DNMTi) / AZA+ITF, 16 samples. NCBI-computed raw counts. | QC, DE, contrasts, masking, enrichment, TF, STRING, WGCNA, sPLS-DA, export, pipelines |
| **GSE182119** TE counts | Per-sample Telescope locus reports for the same 16 runs (11,448 L1/HERV loci kept), renamed to `<id>:<Class>:<Family>` | TE module, sensor scores, mediation, DIABLO, WGCNA TE traits |
| Telescope `retro.hg38.v1` GTF | Coordinates for those loci | `omni_build_te_coords()`, proximity, GREAT |
| **GSE50244** islet bulk | 77 donors with HbA1c (66 normal, 11 T2D; MuSiC tutorial copy) | single-cell reference deconvolution |
| **Xin et al. 2016** islet scRNA-seq | 1,492 cells, healthy + T2D donors (alpha, beta, delta, gamma) | `run_omni_sc_deconv()`, `run_omni_disease_state()`, sc plots |
| CollecTRI | TF regulons from the OmniPath web API | TF activity (see notes) |

DNMT inhibition de-represses TEs and induces a "viral mimicry" interferon
response, which is what the sensor/mediation modules model. The 2×2 drug
design blocked by cell line also gives a genuine AZA × ITF interaction term
for the interaction plots.

## Running

```bash
Rscript dev/real-data-tests/download_test_data.R   # ~110 MB into workspace/data/raw
Rscript dev/real-data-tests/prepare_test_data.R    # builds workspace/data/processed
Rscript dev/real-data-tests/test_all_functions.R   # ~40 min; results + plots in workspace/test_results
Rscript dev/real-data-tests/test_followup.R        # re-checks fixed functions
Rscript dev/real-data-tests/test_app_real_data.R   # drives the Shiny app modules
```

Each check records PASS/FAIL, runtime, warnings and output; every plot is
saved as PNG and also rendered to a PDF device. Export checks verify that the
expected files were written, not just that no error was raised.

## Results (R 4.3.3, Bioconductor 3.18)

180 checks covering 112 of the 112 exports: **174 pass, 4 fail, 2 skipped**.
`results/function_results*.csv` has every check.

- **Fail (4):** `run_omni_deconvolution()` (immunedeconv and MuSiC methods)
  and `run_omni_sc_deconv(method = "music" / "auto")`. immunedeconv and MuSiC
  are GitHub-only and could not be installed in the test container. The
  functions stop with a clear install hint.
  `plot_omni_deconv_heatmap()` / `plot_omni_deconv_corr()` were tested with
  BisqueRNA proportions instead.
- **Skipped (2):** `omni_install_deps()` (installs system-wide) and
  `omni_launch_app()` (covered by the app tests).
- **`run_omniseq_guided()`** is menu-driven. Only its non-interactive guard
  was checked.

### Bugs found and fixed

See `NEWS.md` (development version, "Bug fixes") for the full list. In
short: `omni_build_te_coords()` (`%>%` not imported); `export_omniseq_report()`
(broken with a relative `outdir`, including the default, while reporting
success); STRING duplicate vertices; graph plots failing on PDF devices;
sensor-score metadata and per-shape regression lines; `plot_omni_sc_profile()`
always erroring; `NA` reference group in the sc genotype plots; Bisque
same-individual mode and silent no-op; deconvolution heatmap dropping
non-immune cell types; interactive plot colours ignored; outdated msigdbr
arguments; misleading install hints and installer gaps.

### Notes on the data (not bugs)

- AZA alone yields 139 DEGs and ITF alone 0 across these four heterogeneous
  lines (n = 4 per group), while AZA+ITF yields 1,729. With weak contrasts,
  ashr shrinks most log2FCs to ~0, so GSEA on the shrunken fold changes finds
  nothing. Ranking by the Wald statistic (stored as `stat`) would be more
  robust for low-power contrasts.
- The AZA × ITF interaction has no significant genes with one replicate per
  cell line and condition, so the interaction plots show the null case.
- Telescope names HERV-K as `HML2` and has no Alu/SINE loci, so the default
  (RepeatMasker-named) sensor families barely match. The tests pass
  Telescope family names explicitly, and the Shiny app now accepts custom
  sensor families.

### Environment workarounds (container-specific)

- **OmnipathR** cannot build its organism table because ensembl.org returns
  403 to the container. CollecTRI was fetched from omnipathdb.org and placed
  in `obj$params$tf_regulons`, which `run_omni_tf()` reuses.
- **Bioconductor 3.18** has version conflicts with current CRAN releases.
  Pinned from the CRAN archive: ggplot2 3.5.2, ggiraph 0.8.13, tidytree
  0.4.6, BisqueRNA 1.0.5 (archived). Debian's `r-cran-dt` ships a broken
  noUiSlider, which breaks DT column filters in the browser, so DT was
  installed from CRAN.
