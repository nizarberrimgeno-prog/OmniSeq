# =============================================================================
# Export module — full export_omniseq() bundle (ZIP), HTML report, the
# OmniSeq object itself (.rds), the session log, and dependency status.
# =============================================================================

mod_export_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_column_wrap(
    width = 1 / 2,
    bslib::card(
      bslib::card_header("Download results"),
      shiny::p("Everything is generated from the current object, including",
               "every cached contrast, enrichment and module result."),
      shiny::textInput(ns("prefix"), "File prefix (optional)"),
      shiny::fluidRow(
        shiny::column(4, shiny::numericInput(ns("width"), "Plot width (in)", 10,
                                             4, 30)),
        shiny::column(4, shiny::numericInput(ns("height"), "Plot height (in)", 8,
                                             4, 30)),
        shiny::column(4, shiny::numericInput(ns("dpi"), "DPI", 300, 72, 600))
      ),
      shiny::downloadButton(ns("dl_zip"), "All tables + figures (ZIP)",
                            class = "btn-primary mb-2"),
      shiny::textInput(ns("report_title"), "Report title",
                       "OmniSeq Analysis Report"),
      shiny::textInput(ns("report_author"), "Author"),
      shiny::downloadButton(ns("dl_report"), "HTML report",
                            class = "btn-outline-primary mb-2"),
      shiny::downloadButton(ns("dl_rds"), "OmniSeq object (.rds)",
                            class = "btn-outline-secondary mb-2"),
      shiny::helpText("Reload a saved object later from the Data tab",
                      "to continue where you left off.")
    ),
    bslib::card(
      bslib::card_header("Session log"),
      shiny::tags$style(shiny::HTML(sprintf(
        "#%s { max-height: 520px; overflow-y: auto; font-size: 0.8rem; }",
        ns("log")))),
      shiny::verbatimTextOutput(ns("log")),
      shiny::downloadButton(ns("dl_log"), "Download log", class = "btn-sm")
    ),
    bslib::card(
      bslib::card_header("Package dependencies"),
      shiny::p("Optional modules need extra packages. Install the missing ones",
               "with", shiny::code("omni_install_deps()"), "and restart the app."),
      DT::DTOutput(ns("deps"))
    )
  )
}

mod_export_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {

    gene_partner <- shiny::reactive({
      p <- state$partner
      if (!is.null(p) && identical(p$params$data_type, "gene") &&
          identical(state$obj$params$data_type, "TE")) p else NULL
    })

    output$dl_zip <- shiny::downloadHandler(
      filename = function() paste0("omniseq_results_", Sys.Date(), ".zip"),
      content  = function(file) {
        obj <- state$obj
        shiny::req(obj)
        outdir <- file.path(tempfile("omniseq_export_"), "omniseq_output")
        dir.create(outdir, recursive = TRUE)
        done <- run_step(state, "Export results", {
          export_omniseq(obj, outdir = outdir, prefix = .nz(input$prefix),
                         gene_obj = gene_partner(), width = input$width,
                         height = input$height, dpi = input$dpi,
                         save_rds = TRUE)
          TRUE
        })
        shiny::req(done)
        old <- setwd(dirname(outdir))
        on.exit(setwd(old), add = TRUE)
        utils::zip(file, files = list.files("omniseq_output", recursive = TRUE,
                                            full.names = TRUE),
                   flags = "-rq")
      }
    )

    output$dl_report <- shiny::downloadHandler(
      filename = function() paste0("omniseq_report_", Sys.Date(), ".html"),
      content  = function(file) {
        obj <- state$obj
        shiny::req(obj)
        outdir <- tempfile("omniseq_report_")
        dir.create(outdir)
        done <- run_step(state, "HTML report", {
          export_omniseq_report(obj, outdir = outdir,
                                report_file = "omniseq_report.html",
                                title = input$report_title,
                                author = input$report_author,
                                gene_obj = gene_partner(), open = FALSE)
          TRUE
        })
        shiny::req(done)
        html <- list.files(outdir, pattern = "\\.html$", recursive = TRUE,
                           full.names = TRUE)
        shiny::req(length(html) > 0)
        file.copy(html[1], file)
      }
    )

    output$dl_rds <- shiny::downloadHandler(
      filename = function() paste0("omniseq_object_", Sys.Date(), ".rds"),
      content  = function(file) {
        shiny::req(state$obj)
        saveRDS(state$obj, file)
      }
    )

    output$log <- shiny::renderText({
      log <- state$log
      if (length(log) == 0) return("Nothing run yet.")
      paste(utils::tail(log, 300), collapse = "\n")
    })
    output$dl_log <- shiny::downloadHandler(
      filename = function() paste0("omniseq_log_", Sys.Date(), ".txt"),
      content  = function(file) writeLines(state$log, file)
    )

    output$deps <- DT::renderDT({
      deps <- NULL
      utils::capture.output(deps <- omni_check_deps())
      DT::datatable(deps, rownames = FALSE,
                    options = list(pageLength = 40, dom = "ft")) |>
        DT::formatStyle("Installed",
                        color = DT::styleEqual(c(TRUE, FALSE),
                                               c("#1a7f37", "#cf222e")))
    })
  })
}
