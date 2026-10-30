# =============================================================================
# Data module — upload counts + metadata (or demo data / a saved object),
# set the design, and build the OmniSeq object with init_omniseq().
# =============================================================================

species_choices <- c("Homo sapiens" = "org.Hs.eg.db",
                     "Mus musculus" = "org.Mm.eg.db")

mod_data_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::radioButtons(ns("source"), "Data source",
                          c("Upload files" = "upload",
                            "Demo dataset" = "demo",
                            "Saved OmniSeq object (.rds)" = "rds")),

      shiny::conditionalPanel(
        sprintf("input['%s'] == 'upload'", ns("source")),
        shiny::fileInput(ns("counts_file"), "Count matrix (CSV / TSV)",
                         accept = c(".csv", ".tsv", ".txt")),
        shiny::helpText("Rows = features (genes or TE loci), columns = samples.",
                        "First column = feature IDs."),
        shiny::fileInput(ns("meta_file"), "Sample metadata (CSV / TSV)",
                         accept = c(".csv", ".tsv", ".txt")),
        shiny::helpText("Rows = samples. First column = sample IDs matching",
                        "the count matrix columns."),
        shiny::checkboxInput(ns("round_counts"),
                             "Round non-integer counts (salmon / RSEM)", FALSE)
      ),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'demo'", ns("source")),
        shiny::helpText("Simulated 8-sample dataset (Control vs Treated,",
                        "two batches). Uses real human gene symbols when",
                        "org.Hs.eg.db is installed.")
      ),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'rds'", ns("source")),
        shiny::fileInput(ns("rds_file"), "OmniSeq object (.rds)",
                         accept = ".rds"),
        shiny::helpText("An object saved from this app or with saveRDS().",
                        "All cached results are restored."),
        shiny::actionButton(ns("load_rds"), "Load object",
                            class = "btn-primary w-100")
      ),

      shiny::conditionalPanel(
        sprintf("input['%s'] != 'rds'", ns("source")),
        shiny::hr(),
        shiny::h6("Design"),
        shiny::selectInput(ns("condition_col"), "Condition column",
                           choices = NULL),
        shiny::selectInput(ns("ref_level"), "Reference level (baseline)",
                           choices = NULL),
        shiny::selectInput(ns("covariates"), "Covariates (optional)",
                           choices = NULL, multiple = TRUE),
        shiny::selectInput(ns("species"), "Species",
                           choices = names(species_choices)),
        shiny::checkboxInput(ns("filter_coding"),
                             "Keep protein-coding gene symbols only", FALSE),
        shiny::textAreaInput(ns("exclude_genes"),
                             "Exclude features (one per line, optional)",
                             rows = 2),
        shiny::actionButton(ns("init"), "Create OmniSeq object",
                            class = "btn-primary w-100")
      )
    ),

    bslib::navset_card_tab(
      bslib::nav_panel("Status",
                       shiny::uiOutput(ns("problems")),
                       shiny::verbatimTextOutput(ns("summary"))),
      bslib::nav_panel("Metadata", DT::DTOutput(ns("meta_preview"))),
      bslib::nav_panel("Counts (first 200 rows)",
                       DT::DTOutput(ns("counts_preview")))
    )
  )
}

mod_data_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- Raw inputs ---------------------------------------------------------
    demo <- shiny::reactive({
      shiny::req(input$source == "demo")
      demo_data()
    })

    counts_raw <- shiny::reactive({
      if (input$source == "demo") return(as.data.frame(demo()$counts))
      shiny::req(input$counts_file)
      tryCatch(read_table_upload(input$counts_file$datapath,
                                 input$counts_file$name),
               error = function(e) {
                 shiny::showNotification(conditionMessage(e), type = "error")
                 NULL
               })
    })

    meta_raw <- shiny::reactive({
      if (input$source == "demo") return(demo()$meta)
      shiny::req(input$meta_file)
      tryCatch(read_table_upload(input$meta_file$datapath,
                                 input$meta_file$name),
               error = function(e) {
                 shiny::showNotification(conditionMessage(e), type = "error")
                 NULL
               })
    })

    checked <- shiny::reactive({
      counts <- counts_raw(); meta <- meta_raw()
      shiny::req(counts, meta)
      if (isTRUE(input$round_counts))
        counts[] <- lapply(counts, function(x) if (is.numeric(x)) round(x) else x)
      check_inputs(counts, meta)
    })

    # ---- Design selectors ---------------------------------------------------
    shiny::observeEvent(meta_raw(), {
      cols <- colnames(meta_raw())
      default <- if ("Condition" %in% cols) "Condition" else cols[1]
      shiny::updateSelectInput(session, "condition_col", choices = cols,
                               selected = default)
      shiny::updateSelectInput(session, "covariates",
                               choices = setdiff(cols, default),
                               selected = if (input$source == "demo") "Batch")
    })

    shiny::observeEvent(input$condition_col, {
      meta <- meta_raw()
      shiny::req(meta, input$condition_col %in% colnames(meta))
      lv <- unique(as.character(stats::na.omit(meta[[input$condition_col]])))
      if (is.factor(meta[[input$condition_col]]))
        lv <- levels(meta[[input$condition_col]])
      shiny::updateSelectInput(session, "ref_level", choices = lv,
                               selected = lv[1])
      shiny::updateSelectInput(session, "covariates",
                               choices = setdiff(colnames(meta),
                                                 input$condition_col),
                               selected = setdiff(input$covariates,
                                                  input$condition_col))
    })

    # ---- Build object -------------------------------------------------------
    shiny::observeEvent(input$init, {
      chk <- checked()
      if (length(chk$problems) > 0) {
        shiny::showNotification("Fix the input problems listed on the Status tab.",
                                type = "error")
        return()
      }
      meta <- meta_raw()[colnames(chk$counts), , drop = FALSE]
      excl <- trimws(strsplit(input$exclude_genes %||% "", "[\n,;]")[[1]])
      excl <- excl[nzchar(excl)]

      obj <- run_step(state, "Create OmniSeq object", {
        init_omniseq(
          count_matrix  = chk$counts,
          metadata      = meta,
          condition_col = input$condition_col,
          covariates    = .nz(input$covariates),
          exclude_genes = .nz(excl),
          filter_coding = input$filter_coding,
          species       = input$species,
          org_db        = species_choices[[input$species]],
          ref_level     = .nz(input$ref_level)
        )
      })
      if (!is.null(obj)) state$obj <- obj
    })

    shiny::observeEvent(input$load_rds, {
      shiny::req(input$rds_file)
      obj <- run_step(state, "Load OmniSeq object", {
        x <- readRDS(input$rds_file$datapath)
        validate_omniseq(x)
      })
      if (!is.null(obj)) state$obj <- obj
    })

    # ---- Outputs ------------------------------------------------------------
    output$problems <- shiny::renderUI({
      if (input$source == "rds") return(NULL)
      chk <- checked()
      if (length(chk$problems) == 0)
        return(shiny::div(class = "alert alert-success",
                          sprintf("Inputs look valid: %d features x %d samples.",
                                  nrow(chk$counts), ncol(chk$counts))))
      shiny::div(class = "alert alert-danger",
                 shiny::tags$ul(lapply(chk$problems, shiny::tags$li)))
    })

    output$summary <- shiny::renderPrint({
      obj <- need_obj(state)
      print(obj)
    })

    output$meta_preview <- DT::renderDT({
      meta <- if (!is.null(state$obj)) omni_metadata(state$obj) else meta_raw()
      shiny::req(meta)
      DT::datatable(meta, options = list(pageLength = 25, scrollX = TRUE))
    })

    output$counts_preview <- DT::renderDT({
      counts <- if (!is.null(state$obj)) omni_counts(state$obj) else counts_raw()
      shiny::req(counts)
      DT::datatable(utils::head(as.data.frame(counts), 200),
                    options = list(pageLength = 15, scrollX = TRUE))
    })
  })
}
