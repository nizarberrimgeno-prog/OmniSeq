# =============================================================================
# OmniSeq Shiny App
#
# No analysis logic lives here. Every computation is delegated to the
# OmniSeq package functions; the modules in R/ (sourced automatically by
# shiny::runApp()) only wire them to reactive UI inputs and outputs.
#
# Shared state passed to every module:
#   state$obj      — the main OmniSeq object
#   state$partner  — optional paired object (gene + TE) for multi-omics
#   state$log      — messages emitted by package functions
# =============================================================================

suppressPackageStartupMessages(library(OmniSeq))

options(shiny.maxRequestSize = 1024^3)  # count matrices can be large

ui <- bslib::page_navbar(
  title = "OmniSeq",
  id    = "main_nav",
  theme = bslib::bs_theme(version = 5, primary = "#2166AC"),
  fillable = FALSE,

  sidebar = bslib::sidebar(
    width = 270,
    open  = "desktop",
    shiny::h6("Analysis status"),
    shiny::uiOutput("status"),
    shiny::hr(),
    shiny::selectInput("active_contrast", "Active contrast", choices = NULL),
    shiny::helpText("DE plots, enrichment, TF and network modules all",
                    "use the active contrast."),
    shiny::hr(),
    shiny::actionButton("reset", "Start over", class = "btn-sm btn-outline-danger")
  ),

  bslib::nav_panel("Data",             mod_data_ui("data")),
  bslib::nav_panel("QC",               mod_qc_ui("qc")),
  bslib::nav_panel("Differential expression", mod_de_ui("de")),
  bslib::nav_panel("Enrichment",       mod_enrichment_ui("enrich")),
  bslib::nav_panel("Biology",          mod_biology_ui("bio")),
  bslib::nav_panel("Multivariate",     mod_multivariate_ui("multi")),
  bslib::nav_panel("TE analysis",      mod_te_ui("te")),
  bslib::nav_panel("Multi-omics",      mod_multiomics_ui("omics")),
  bslib::nav_panel("Export",           mod_export_ui("export")),
  bslib::nav_spacer(),
  bslib::nav_item(shiny::tags$small(class = "text-muted",
                                    paste0("v", utils::packageVersion("OmniSeq"))))
)

server <- function(input, output, session) {

  state <- shiny::reactiveValues(obj = NULL, partner = NULL,
                                 log = character(0))

  mod_data_server("data", state)
  mod_qc_server("qc", state)
  mod_de_server("de", state)
  mod_enrichment_server("enrich", state)
  mod_biology_server("bio", state)
  mod_multivariate_server("multi", state)
  mod_te_server("te", state)
  mod_multiomics_server("omics", state)
  mod_export_server("export", state)

  # ---- Active contrast (app-wide) -----------------------------------------
  shiny::observeEvent(state$obj, {
    obj <- state$obj
    cs  <- if (is.null(obj)) character(0) else omni_contrasts(obj)
    shiny::updateSelectInput(session, "active_contrast", choices = cs,
                             selected = if (length(cs)) obj$active_contrast)
  }, ignoreNULL = FALSE)

  shiny::observeEvent(input$active_contrast, {
    obj <- state$obj
    shiny::req(obj, nzchar(input$active_contrast))
    if (!identical(obj$active_contrast, input$active_contrast))
      state$obj <- omni_set_contrast(obj, input$active_contrast)
  })

  # ---- Status checklist ----------------------------------------------------
  output$status <- shiny::renderUI({
    obj <- state$obj
    step <- function(done, label) {
      shiny::div(
        class = if (done) "text-success" else "text-muted",
        shiny::HTML(if (done) "&#10003; " else "&#9675; "), label
      )
    }
    if (is.null(obj))
      return(shiny::div(class = "text-muted small",
                        "No data loaded. Start on the Data tab."))
    has <- function(x) !is.null(x) && length(x) > 0
    ac  <- obj$active_contrast
    shiny::div(
      class = "small",
      shiny::div(sprintf("%s x %s %s", format(nrow(obj$data), big.mark = ","),
                         ncol(obj$data),
                         if (obj$params$data_type == "TE") "TE loci" else "genes")),
      shiny::div(class = "text-muted mb-2", deparse(obj$params$formula)),
      step(TRUE, "Object created"),
      step(has(omni_qc(obj)), "QC"),
      step(has(omni_model(obj)), "DESeq2 model"),
      step(has(omni_contrasts(obj)),
           sprintf("Contrasts (%d)", length(omni_contrasts(obj)))),
      step(has(obj$enrichment_cache[[ac]]), "ORA"),
      step(has(obj$gsea_cache[[ac]]), "GSEA"),
      step(has(obj$advanced_ssgsea), "ssGSEA"),
      step(has(obj$tf_results) || has(obj$tf_contrast_cache), "TF activity"),
      step(has(obj$string_network), "STRING"),
      step(has(obj$deconv), "Deconvolution"),
      step(has(obj$wgcna), "WGCNA"),
      step(has(obj$splsda), "sPLS-DA"),
      if (obj$params$data_type == "TE") step(has(obj$te_family), "TE families"),
      step(!is.null(state$partner), "Paired object")
    )
  })

  shiny::observeEvent(input$reset, {
    shiny::showModal(shiny::modalDialog(
      "Discard the current object and all results?",
      footer = shiny::tagList(
        shiny::modalButton("Cancel"),
        shiny::actionButton("reset_ok", "Discard", class = "btn-danger")
      )
    ))
  })
  shiny::observeEvent(input$reset_ok, {
    state$obj <- NULL
    state$partner <- NULL
    state$log <- character(0)
    shiny::removeModal()
    bslib::nav_select("main_nav", "Data")
  })
}

shiny::shinyApp(ui, server)
