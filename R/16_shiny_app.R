# =============================================================================
# OmniSeq — 16_shiny_app.R
#
# Launcher for the bundled Shiny app (inst/shiny-app/). The app itself
# contains no analysis logic — it only orchestrates the existing
# init_omniseq() / run_omni_qc() / run_omni_de() / plot_omni_*() functions.
# =============================================================================

#' Launch the OmniSeq Interactive App
#'
#' Starts a local Shiny app that runs the OmniSeq workflow without writing
#' any R code: data upload (or a demo dataset / saved object), QC,
#' differential expression with custom contrasts and signature masking,
#' enrichment (ORA, GSEA, ssGSEA), TF activity, STRING, deconvolution,
#' WGCNA, sPLS-DA, TE analysis, gene + TE multi-omics (sensor scores,
#' DIABLO, mediation), and export of tables, figures and an HTML report.
#'
#' Modules whose optional dependencies are missing stay visible but show
#' which packages to install (see [omni_check_deps()]).
#'
#' @param ... Additional arguments passed to [shiny::runApp()] (e.g.
#'   `port`, `launch.browser`).
#'
#' @return Called for its side effect (starts the Shiny app). Does not return.
#' @export
#'
#' @examples
#' \dontrun{
#' omni_launch_app()
#' }
omni_launch_app <- function(...) {
  .require_pkg("shiny")
  .require_pkg("bslib")
  .require_pkg("DT")

  app_dir <- system.file("shiny-app", package = "OmniSeq")
  if (app_dir == "")
    stop("Could not find the OmniSeq Shiny app. Try re-installing the package.",
         call. = FALSE)

  shiny::runApp(app_dir, ...)
}
