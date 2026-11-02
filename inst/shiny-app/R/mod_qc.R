# =============================================================================
# QC module — run_omni_qc() and the QC plot family.
# =============================================================================

mod_qc_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 300,
      shiny::checkboxInput(ns("calc_distances"),
                           "Compute sample distances", TRUE),
      shiny::actionButton(ns("run"), "Run QC", class = "btn-primary w-100"),
      shiny::hr(),
      shiny::helpText("QC uses a blind VST of the raw counts. The DE tab",
                      "fits the model used by every downstream module.")
    ),
    bslib::layout_column_wrap(
      width = 1 / 2,
      plot_card("Library size", ns("libsize"), "380px"),
      plot_card("Detected features", ns("detected"), "380px"),
      plot_card("Sample distances", ns("distance"), "480px"),
      plot_card("Mean-variance", ns("mean_var"), "480px")
    ),
    plot_card("PCA: raw vs VST", ns("pca"), "480px")
  )
}

mod_qc_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(input$run, {
      obj <- need_obj(state)
      out <- run_step(state, "QC", {
        run_omni_qc(obj, calc_distances = input$calc_distances)
      })
      if (!is.null(out)) state$obj <- out
    })

    qc_obj <- shiny::reactive({
      obj <- need_obj(state)
      shiny::validate(shiny::need(!is.null(omni_qc(obj)),
                                  "Click 'Run QC' to compute QC metrics."))
      obj
    })

    output$libsize <- render_smart(
      safe_plot(plot_omni_qc_basics(qc_obj(), plot_type = "library_size")))
    output$detected <- render_smart(
      safe_plot(plot_omni_qc_basics(qc_obj(), plot_type = "detected_genes")))
    output$distance <- render_smart(
      safe_plot(plot_omni_qc_advanced(qc_obj(), plot_type = "distance")))
    output$mean_var <- render_smart(
      safe_plot(plot_omni_qc_advanced(qc_obj(), plot_type = "mean_var")))
    output$pca <- render_smart(
      safe_plot(plot_omni_qc_advanced(qc_obj(), plot_type = "pca_compare")))
  })
}
