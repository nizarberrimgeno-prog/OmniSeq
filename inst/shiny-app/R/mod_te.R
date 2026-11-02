# =============================================================================
# TE module — transposable element analysis on a TE-type OmniSeq object:
# family-level enrichment, TE plots, and genomic proximity + GREAT.
# =============================================================================

mod_te_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::uiOutput(ns("te_status")),
      bslib::accordion(
        open = "Family analysis",
        bslib::accordion_panel(
          "Family analysis",
          shiny::sliderInput(ns("lfc"), "|log2FC| >", 0, 4, 1, step = 0.25),
          shiny::selectInput(ns("pval"), "padj <", c(0.001, 0.01, 0.05, 0.1),
                             selected = 0.05),
          shiny::selectInput(ns("direction"), "Direction",
                             c("Both" = "both", "Up" = "up", "Down" = "down")),
          shiny::actionButton(ns("run_family"), "Run TE family analysis",
                              class = "btn-primary w-100")
        ),
        bslib::accordion_panel(
          "Proximity + GREAT",
          deps_notice(c("GenomicRanges", "rGREAT")),
          shiny::fileInput(ns("coords"), "TE coordinates (.rds or .csv)",
                           accept = c(".rds", ".csv")),
          shiny::helpText("Output of omni_build_te_coords(): columns",
                          "Locus_ID, Chr, Start, End, Strand."),
          shiny::selectInput(ns("prox_dir"), "Direction",
                             c("Up" = "up", "Both" = "both", "Down" = "down")),
          shiny::selectInput(ns("genome"), "Genome",
                             c("hg38", "hg19", "mm10", "mm39")),
          shiny::actionButton(ns("run_great"), "Run proximity + GREAT",
                              class = "btn-outline-primary w-100")
        ),
        bslib::accordion_panel(
          "Plot settings",
          shiny::sliderInput(ns("n_fam"), "Families shown", 5, 50, 20),
          shiny::sliderInput(ns("top_n"), "Volcano labels", 0, 30, 10),
          shiny::selectInput(ns("label_col"), "Volcano labels by",
                             c("Family" = "family", "Locus" = "locus")),
          shiny::selectInput(ns("color_by"), "Family plot colour",
                             c("Class" = "class", "log2FC" = "lfc")),
          shiny::numericInput(ns("min_loci"), "Density: min loci per class",
                              50, 1, 1000)
        )
      )
    ),
    bslib::navset_card_tab(
      bslib::nav_panel("TE volcano",
                       plot_slot(ns("volcano"), height = "580px")),
      bslib::nav_panel("Families",
                       plot_slot(ns("family"), height = "620px")),
      bslib::nav_panel("Class composition",
                       plot_slot(ns("pie"), height = "520px")),
      bslib::nav_panel("log2FC density",
                       plot_slot(ns("density"), height = "520px")),
      bslib::nav_panel("Family heatmap",
                       plot_slot(ns("heatmap"), height = "700px")),
      bslib::nav_panel("GREAT",
                       plot_slot(ns("great"), height = "600px")),
      bslib::nav_panel("Family table", DT::DTOutput(ns("table")),
                       shiny::downloadButton(ns("dl_table"), "Download CSV",
                                             class = "btn-sm"))
    )
  )
}

mod_te_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {

    is_te <- shiny::reactive(!is.null(state$obj) &&
                               identical(state$obj$params$data_type, "TE"))

    output$te_status <- shiny::renderUI({
      if (is.null(state$obj) || is_te()) return(NULL)
      shiny::div(class = "alert alert-info small",
                 "The loaded object holds gene-level data. TE analysis needs",
                 "a count matrix with TE locus IDs (e.g.",
                 shiny::code("1000014:ERV:MSTB"), ").")
    })

    te_contrast <- shiny::reactive({
      obj <- need_contrast(state)
      shiny::validate(shiny::need(is_te(), "Load a TE-level count matrix to use this tab."))
      obj
    })

    shiny::observeEvent(input$run_family, {
      shiny::req(state$obj)
      out <- run_step(state, "TE family analysis", {
        run_omni_te_family(te_contrast(),
                           lfc_thresh  = input$lfc,
                           pval_thresh = as.numeric(input$pval),
                           direction   = input$direction)
      })
      if (!is.null(out)) state$obj <- out
    })

    shiny::observeEvent(input$run_great, {
      shiny::req(state$obj, input$coords)
      out <- run_step(state, "TE proximity + GREAT", {
        coords <- if (grepl("\\.rds$", input$coords$name, ignore.case = TRUE))
          readRDS(input$coords$datapath) else
            utils::read.csv(input$coords$datapath, stringsAsFactors = FALSE)
        obj <- run_omni_te_proximity(te_contrast(), te_coords = coords,
                                     lfc_thresh  = input$lfc,
                                     pval_thresh = as.numeric(input$pval),
                                     direction   = input$prox_dir)
        run_omni_te_great(obj, genome = input$genome)
      })
      if (!is.null(out)) state$obj <- out
    })

    fam_obj <- shiny::reactive({
      obj <- te_contrast()
      shiny::validate(shiny::need(obj$te_family,
                                  "Click 'Run TE family analysis'."))
      obj
    })

    output$volcano <- render_smart(safe_plot(
      plot_omni_te_volcano(te_contrast(), lfc_thresh = input$lfc,
                           pval_thresh = as.numeric(input$pval),
                           top_n = input$top_n, label_col = input$label_col)))
    output$family <- render_smart(safe_plot(
      plot_omni_te_family(fam_obj(), n_families = input$n_fam,
                          color_by = input$color_by)))
    output$pie <- render_smart(safe_plot(plot_omni_te_class_pie(fam_obj())))
    output$density <- render_smart(safe_plot(
      plot_omni_te_lfc_density(te_contrast(), min_loci = input$min_loci)))
    output$heatmap <- render_smart(safe_plot(
      plot_omni_te_heatmap(te_contrast(), n_families = input$n_fam)))
    output$great <- render_smart({
      obj <- te_contrast()
      shiny::validate(shiny::need(obj$te_great,
                                  "Upload TE coordinates and run proximity + GREAT."))
      safe_plot(plot_omni_te_great(obj, show_n = input$n_fam))
    })

    fam_df <- shiny::reactive({
      as.data.frame(fam_obj()$te_family$by_family)
    })
    output$table <- DT::renderDT(omni_dt(fam_df()))
    output$dl_table <- csv_download("TE_families.csv", fam_df)
  })
}
