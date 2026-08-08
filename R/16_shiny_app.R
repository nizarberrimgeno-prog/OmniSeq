# =============================================================================
# OmniSeq — 16_shiny_app.R
#
# Launcher for the bundled Shiny app (inst/shiny-app/). The app itself
# contains no analysis logic — it only orchestrates the existing
# init_omniseq() / run_omni_qc() / run_omni_de() / plot_omni_*() functions.
# =============================================================================

#' Launch the OmniSeq Interactive App
#'
#' Starts a local Shiny app for running the core OmniSeq workflow —
#' upload counts and metadata, run QC, run differential expression, and
#' browse the results — without writing any R code.
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
