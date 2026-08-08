# =============================================================================
# OmniSeq Shiny App — v1 (core workflow: upload -> QC -> DE -> results)
#
# No analysis logic lives here. Every computation is delegated to the
# OmniSeq package functions (init_omniseq, run_omni_qc, run_omni_de,
# plot_omni_*); this file only wires them to reactive UI inputs/outputs.
#
# Launched via OmniSeq::omni_launch_app(), which loads the package first —
# so init_omniseq()/run_omni_qc()/etc. are called unqualified here.
# =============================================================================

ui <- bslib::page_sidebar(
  title = "OmniSeq",
  theme = bslib::bs_theme(version = 5, primary = "#2166AC"),

  sidebar = bslib::sidebar(
    width = 340,

    shiny::h5("1. Upload data"),
    shiny::fileInput("counts_file", "Count matrix (CSV)",
                     accept = ".csv"),
    shiny::helpText("Rows = features, columns = samples. First column = feature IDs."),
    shiny::fileInput("meta_file", "Sample metadata (CSV)",
                     accept = ".csv"),
    shiny::helpText("Rows = samples. First column = sample IDs, must match count matrix columns."),

    shiny::hr(),
    shiny::h5("2. Design"),
    shiny::selectInput("condition_col", "Condition column", choices = NULL),
    shiny::selectInput("covariates", "Covariates (optional)",
                       choices = NULL, multiple = TRUE),
    shiny::actionButton("run_qc", "Initialize & Run QC",
                        class = "btn-primary w-100"),

    shiny::hr(),
    shiny::h5("3. Differential Expression"),
    shiny::selectInput("shrinkage", "LFC shrinkage",
                       choices = c("apeglm", "ashr", "none"),
                       selected = "apeglm"),
    shiny::actionButton("run_de", "Run DE", class = "btn-primary w-100"),

    shiny::hr(),
    shiny::downloadButton("download_results", "Download results (CSV)",
                          class = "w-100")
  ),

  bslib::navset_tab(
    bslib::nav_panel(
      "Overview",
      shiny::verbatimTextOutput("overview_text")
    ),
    bslib::nav_panel(
      "QC",
      bslib::layout_column_wrap(
        width = 1 / 2,
        shiny::plotOutput("qc_libsize"),
        shiny::plotOutput("qc_detected")
      ),
      shiny::plotOutput("qc_pca")
    ),
    bslib::nav_panel(
      "Differential Expression",
      bslib::layout_column_wrap(
        width = 1 / 2,
        shiny::plotOutput("de_volcano"),
        shiny::plotOutput("de_ma")
      ),
      DT::DTOutput("de_table")
    )
  )
)

server <- function(input, output, session) {

  metadata_rv <- shiny::reactive({
    shiny::req(input$meta_file)
    utils::read.csv(input$meta_file$datapath, row.names = 1,
                    check.names = FALSE)
  })

  shiny::observeEvent(metadata_rv(), {
    cols <- colnames(metadata_rv())
    shiny::updateSelectInput(session, "condition_col", choices = cols)
    shiny::updateSelectInput(session, "covariates", choices = cols)
  })

  obj_rv <- shiny::reactiveVal(NULL)

  shiny::observeEvent(input$run_qc, {
    shiny::req(input$counts_file, input$meta_file, input$condition_col)

    tryCatch({
      counts <- as.matrix(
        utils::read.csv(input$counts_file$datapath, row.names = 1,
                        check.names = FALSE)
      )
      meta <- metadata_rv()

      obj <- init_omniseq(
        count_matrix  = counts,
        metadata      = meta,
        condition_col = input$condition_col,
        covariates    = input$covariates
      )
      obj <- run_omni_qc(obj)
      obj_rv(obj)

      shiny::showNotification("QC complete.", type = "message")
    }, error = function(e) {
      shiny::showNotification(paste("Error:", e$message), type = "error",
                              duration = NULL)
    })
  })

  shiny::observeEvent(input$run_de, {
    shiny::req(obj_rv())

    tryCatch({
      obj <- run_omni_de(obj_rv(), shrinkage = input$shrinkage)
      obj_rv(obj)
      shiny::showNotification("Differential expression complete.",
                              type = "message")
    }, error = function(e) {
      shiny::showNotification(paste("Error:", e$message), type = "error",
                              duration = NULL)
    })
  })

  output$overview_text <- shiny::renderPrint({
    shiny::req(obj_rv())
    print(obj_rv())
  })

  output$qc_libsize <- shiny::renderPlot({
    shiny::req(obj_rv())
    plot_omni_qc_basics(obj_rv(), plot_type = "library_size")
  })

  output$qc_detected <- shiny::renderPlot({
    shiny::req(obj_rv())
    plot_omni_qc_basics(obj_rv(), plot_type = "detected_genes")
  })

  output$qc_pca <- shiny::renderPlot({
    shiny::req(obj_rv())
    plot_omni_qc_advanced(obj_rv(), plot_type = "pca_compare")
  })

  output$de_volcano <- shiny::renderPlot({
    shiny::req(obj_rv())
    shiny::req(length(omni_contrasts(obj_rv())) > 0)
    plot_omni_volcano(obj_rv())
  })

  output$de_ma <- shiny::renderPlot({
    shiny::req(obj_rv())
    shiny::req(length(omni_contrasts(obj_rv())) > 0)
    plot_omni_ma(obj_rv())
  })

  output$de_table <- DT::renderDT({
    shiny::req(obj_rv())
    shiny::req(length(omni_contrasts(obj_rv())) > 0)
    res <- as.data.frame(omni_results(obj_rv()))
    res <- cbind(Feature = rownames(res), res)
    DT::datatable(res, rownames = FALSE,
                  options = list(pageLength = 15, scrollX = TRUE))
  })

  output$download_results <- shiny::downloadHandler(
    filename = function() "omniseq_results.csv",
    content = function(file) {
      shiny::req(obj_rv())
      res <- as.data.frame(omni_results(obj_rv()))
      res <- cbind(Feature = rownames(res), res)
      utils::write.csv(res, file, row.names = FALSE)
    }
  )
}

shiny::shinyApp(ui, server)
