# =============================================================================
# Biology module — TF activity (decoupleR), STRING PPI network, and bulk
# cell-type deconvolution (immunedeconv).
# =============================================================================

mod_biology_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::navset_card_tab(
    bslib::nav_panel(
      "TF activity",
      deps_notice(c("decoupleR", "OmnipathR")),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 300,
          shiny::radioButtons(ns("tf_mode"), "Mode",
                              c("Per contrast (Wald stat)" = "contrast",
                                "Per sample (one group)"   = "samples")),
          shiny::conditionalPanel(
            sprintf("input['%s'] == 'samples'", ns("tf_mode")),
            shiny::selectInput(ns("tf_group"), "Sample group", choices = NULL)
          ),
          shiny::selectInput(ns("tf_db"), "Regulon database",
                             c("CollecTRI" = "collectri",
                               "DoRothEA" = "dorothea", "Both" = "both")),
          shiny::selectInput(ns("tf_method"), "Method",
                             c("ULM" = "ulm", "MLM" = "mlm", "WSUM" = "wsum")),
          shiny::actionButton(ns("run_tf"), "Run TF inference",
                              class = "btn-primary w-100"),
          shiny::hr(),
          shiny::sliderInput(ns("tf_n"), "TFs shown", 10, 60, 30),
          shiny::sliderInput(ns("tf_min_reg"), "Min regulon size", 0, 50, 15),
          shiny::sliderInput(ns("tf_targets"), "Network: targets per TF",
                             1, 15, 5)
        ),
        bslib::navset_underline(
          bslib::nav_panel("Activity",
                           shiny::plotOutput(ns("tf_activity"), height = "650px")),
          bslib::nav_panel("TF-target network",
                           shiny::plotOutput(ns("tf_network"), height = "700px"))
        )
      )
    ),

    bslib::nav_panel(
      "STRING network",
      deps_notice(c("STRINGdb", "igraph", "ggraph")),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 300,
          shiny::sliderInput(ns("string_n"), "Top DEGs queried", 10, 300, 50,
                             step = 10),
          shiny::selectInput(ns("string_pval"), "padj <",
                             c(0.001, 0.01, 0.05, 0.1), selected = 0.05),
          shiny::sliderInput(ns("string_score"), "Min interaction score",
                             150, 900, 400, step = 50),
          shiny::actionButton(ns("run_string"), "Query STRING",
                              class = "btn-primary w-100"),
          shiny::helpText("First run downloads ~70 MB of STRING data into",
                          shiny::code("~/.omniseq_cache/string"), "."),
          shiny::hr(),
          shiny::sliderInput(ns("string_degree"), "Min node degree", 1, 10, 1),
          shiny::checkboxInput(ns("string_hubs"), "Label hubs", TRUE)
        ),
        shiny::plotOutput(ns("string_plot"), height = "700px")
      )
    ),

    bslib::nav_panel(
      "Deconvolution",
      deps_notice(c("immunedeconv")),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 300,
          shiny::selectInput(ns("deconv_method"), "Method",
                             c("xCell" = "xcell", "EPIC" = "epic",
                               "quanTIseq" = "quantiseq", "TIMER" = "timer",
                               "CIBERSORT" = "cibersort")),
          shiny::actionButton(ns("run_deconv"), "Run deconvolution",
                              class = "btn-primary w-100"),
          shiny::hr(),
          shiny::checkboxInput(ns("deconv_aggregate"),
                               "Aggregate to main cell types", TRUE),
          shiny::textInput(ns("deconv_feature"),
                           "Correlate with gene / pathway",
                           placeholder = "e.g. IFI27 or HALLMARK_...")
        ),
        bslib::navset_underline(
          bslib::nav_panel("Heatmap",
                           shiny::plotOutput(ns("deconv_heat"), height = "650px")),
          bslib::nav_panel("Correlation",
                           shiny::plotOutput(ns("deconv_corr"), height = "600px")),
          bslib::nav_panel("Table", DT::DTOutput(ns("deconv_table")),
                           shiny::downloadButton(ns("dl_deconv"), "Download CSV",
                                                 class = "btn-sm"))
        )
      )
    )
  )
}

mod_biology_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(state$obj, {
      obj <- state$obj
      lv  <- levels(as.factor(omni_metadata(obj)[[obj$params$condition]]))
      shiny::updateSelectInput(session, "tf_group", choices = lv,
                               selected = .nz(shiny::isolate(input$tf_group)) %||%
                                 lv[length(lv)])
    })

    # ---- TF ----------------------------------------------------------------
    shiny::observeEvent(input$run_tf, {
      shiny::req(state$obj)
      out <- run_step(state, "TF activity", {
        obj <- if (input$tf_mode == "contrast") need_contrast(state) else
          need_model(state)
        run_omni_tf(obj,
                    mode         = input$tf_mode,
                    sample_group = if (input$tf_mode == "samples") input$tf_group,
                    database     = input$tf_db,
                    method       = input$tf_method)
      })
      if (!is.null(out)) state$obj <- out
    })

    tf_obj <- shiny::reactive({
      obj <- need_obj(state)
      shiny::validate(shiny::need(
        !is.null(obj$tf_results) || !is.null(obj$tf_contrast_cache),
        "Click 'Run TF inference' (requires DE results)."))
      obj
    })

    tf_group <- shiny::reactive(if (input$tf_mode == "samples") input$tf_group)

    output$tf_activity <- shiny::renderPlot(safe_plot(
      plot_omni_tf_activity(tf_obj(), n = input$tf_n, group = tf_group(),
                            min_regulon = input$tf_min_reg)))
    output$tf_network <- shiny::renderPlot(safe_plot(
      plot_omni_tf_network(tf_obj(), n_tfs = min(input$tf_n, 40),
                           n_targets = input$tf_targets,
                           plot_type = "static", group = tf_group(),
                           min_regulon = input$tf_min_reg)))

    # ---- STRING ------------------------------------------------------------
    shiny::observeEvent(input$run_string, {
      shiny::req(state$obj)
      species_id <- if (identical(state$obj$params$species, "Mus musculus"))
        10090 else 9606
      out <- run_step(state, "STRING network", {
        run_omni_string(need_contrast(state),
                        n            = input$string_n,
                        pval_thresh  = as.numeric(input$string_pval),
                        string_score = input$string_score,
                        species_id   = species_id)
      })
      if (!is.null(out)) state$obj <- out
    })

    output$string_plot <- shiny::renderPlot({
      obj <- need_obj(state)
      shiny::validate(shiny::need(obj$string_network,
                                  "Click 'Query STRING' (requires DE results)."))
      safe_plot(plot_omni_string_network(obj, plot_type = "static",
                                         min_degree = input$string_degree,
                                         label_hubs = input$string_hubs))
    })

    # ---- Deconvolution -----------------------------------------------------
    shiny::observeEvent(input$run_deconv, {
      shiny::req(state$obj)
      out <- run_step(state, "Deconvolution", {
        run_omni_deconvolution(need_model(state), method = "immunedeconv",
                               deconv_method = input$deconv_method)
      })
      if (!is.null(out)) state$obj <- out
    })

    deconv_obj <- shiny::reactive({
      obj <- need_obj(state)
      shiny::validate(shiny::need(obj$deconv,
                                  "Click 'Run deconvolution' (requires a fitted model)."))
      obj
    })

    output$deconv_heat <- shiny::renderPlot(safe_plot(
      plot_omni_deconv_heatmap(deconv_obj(),
                               aggregate_main = input$deconv_aggregate)))
    output$deconv_corr <- shiny::renderPlot({
      shiny::validate(shiny::need(nzchar(input$deconv_feature),
                                  "Enter a gene or pathway name."))
      safe_plot(plot_omni_deconv_corr(deconv_obj(),
                                      feature = trimws(input$deconv_feature)))
    })

    deconv_df <- shiny::reactive({
      d <- deconv_obj()$deconv
      d <- as.data.frame(d)
      if (!"cell_type" %in% names(d))
        d <- cbind(cell_type = rownames(d), d, row.names = NULL)
      d
    })
    output$deconv_table <- DT::renderDT(omni_dt(deconv_df()))
    output$dl_deconv <- csv_download("deconvolution.csv", deconv_df)
  })
}
