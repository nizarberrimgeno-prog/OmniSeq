# =============================================================================
# Enrichment module — ORA (run_omni_enrichment), GSEA (run_omni_gsea),
# ssGSEA (run_omni_ssgsea) and the cross-contrast comparison plots.
# =============================================================================

msigdb_categories <- c("Hallmark (H)"            = "H",
                       "Curated (C2)"            = "C2",
                       "Ontology (C5)"           = "C5",
                       "Immunologic (C7)"        = "C7",
                       "Cell type (C8)"          = "C8",
                       "Regulatory targets (C3)" = "C3")

msigdb_inputs <- function(ns, prefix) {
  shiny::tagList(
    shiny::selectInput(ns(paste0(prefix, "_category")), "MSigDB collection",
                       msigdb_categories),
    shiny::textInput(ns(paste0(prefix, "_subcategory")),
                     "Subcollection (optional)",
                     placeholder = "e.g. CP:REACTOME, GO:BP")
  )
}

mod_enrichment_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::navset_card_tab(
    bslib::nav_panel(
      "ORA (GO)",
      deps_notice(c("clusterProfiler", "enrichplot")),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 300,
          shiny::selectInput(ns("ora_ont"), "GO ontology",
                             c("Biological process" = "BP",
                               "Molecular function" = "MF",
                               "Cellular component" = "CC")),
          shiny::selectInput(ns("ora_direction"), "Genes",
                             c("Up + down" = "both", "Up only" = "up",
                               "Down only" = "down")),
          shiny::sliderInput(ns("ora_lfc"), "|log2FC| >", 0, 3, 0, step = 0.25),
          shiny::selectInput(ns("ora_pval"), "padj <",
                             c(0.001, 0.01, 0.05, 0.1), selected = 0.05),
          shiny::checkboxInput(ns("ora_simplify"), "Simplify redundant terms",
                               TRUE),
          shiny::actionButton(ns("run_ora"), "Run ORA",
                              class = "btn-primary w-100"),
          shiny::hr(),
          shiny::selectInput(ns("ora_plot"), "Plot",
                             c("Dot plot" = "dot", "Gene-concept network" = "cnet",
                               "Enrichment map" = "emap")),
          shiny::sliderInput(ns("ora_n"), "Terms shown", 5, 40, 15)
        ),
        plot_slot(ns("ora_plot_out"), height = "620px"),
        DT::DTOutput(ns("ora_table")),
        shiny::downloadButton(ns("dl_ora"), "Download CSV", class = "btn-sm")
      )
    ),

    bslib::nav_panel(
      "ORA Up/Down",
      deps_notice(c("clusterProfiler", "msigdbr", "plotly", "htmltools")),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 300,
          shiny::checkboxGroupInput(ns("ud_db"), "Databases",
                                    c("GO-BP", "Hallmark"),
                                    selected = c("GO-BP", "Hallmark")),
          shiny::sliderInput(ns("ud_lfc"), "|log2FC| >", 0, 3, 1, step = 0.25),
          shiny::selectInput(ns("ud_pval"), "padj <",
                             c(0.001, 0.01, 0.05, 0.1), selected = 0.05),
          shiny::actionButton(ns("run_ud"), "Run ORA (active contrast)",
                              class = "btn-primary w-100"),
          shiny::helpText("Up- and down-regulated genes are tested as two",
                          "clusters in one call (run_omni_ora_updown).",
                          "Run it on several contrasts to compare them."),
          shiny::hr(),
          shiny::selectInput(ns("ud_show_db"), "Show database",
                             c("GO-BP", "Hallmark")),
          shiny::sliderInput(ns("ud_n"), "Terms per direction", 5, 40, 20),
          shiny::selectInput(ns("ud_cmp"), "Compare contrasts",
                             choices = NULL, multiple = TRUE)
        ),
        bslib::navset_underline(
          bslib::nav_panel("Active contrast (interactive)",
                           plotly_output(ns("ud_plot"), height = "650px"),
                           DT::DTOutput(ns("ud_table"))),
          bslib::nav_panel("Across contrasts",
                           plot_slot(ns("ud_compare"), height = "750px"))
        )
      )
    ),

    bslib::nav_panel(
      "GSEA",
      deps_notice(c("clusterProfiler", "enrichplot", "msigdbr")),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 300,
          msigdb_inputs(ns, "gsea"),
          shiny::selectInput(ns("gsea_pval"), "p-value cutoff",
                             c(0.01, 0.05, 0.1, 0.25), selected = 0.05),
          shiny::actionButton(ns("run_gsea"), "Run GSEA",
                              class = "btn-primary w-100"),
          shiny::hr(),
          shiny::selectInput(ns("gsea_plot"), "Plot",
                             c("Waterfall" = "waterfall", "Ridge" = "ridge",
                               "Bar" = "bar", "Leading-edge heatmap" = "leading",
                               "Pathway log2FC" = "pathway_lfc")),
          shiny::sliderInput(ns("gsea_n"), "Pathways shown", 5, 40, 15),
          shiny::selectInput(ns("gsea_pathway"), "Pathway (pathway log2FC)",
                             choices = NULL)
        ),
        plot_slot(ns("gsea_plot_out"), height = "620px"),
        DT::DTOutput(ns("gsea_table")),
        shiny::downloadButton(ns("dl_gsea"), "Download CSV", class = "btn-sm")
      )
    ),

    bslib::nav_panel(
      "ssGSEA",
      deps_notice(c("GSVA", "msigdbr", "pheatmap")),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 300,
          msigdb_inputs(ns, "ss"),
          shiny::actionButton(ns("run_ssgsea"), "Run ssGSEA",
                              class = "btn-primary w-100"),
          shiny::hr(),
          shiny::sliderInput(ns("ss_n"), "Pathways shown", 10, 60, 30),
          shiny::checkboxInput(ns("ss_isolate"),
                               "Only samples in the active contrast", TRUE)
        ),
        plot_slot(ns("ss_plot"), height = "720px"),
        shiny::downloadButton(ns("dl_ss"), "Download scores (CSV)",
                              class = "btn-sm")
      )
    ),

    bslib::nav_panel(
      "GO across contrasts",
      deps_notice(c("clusterProfiler", "enrichplot")),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 300,
          shiny::selectInput(ns("cmp_contrasts"), "Contrasts",
                             choices = NULL, multiple = TRUE),
          shiny::sliderInput(ns("cmp_lfc"), "|log2FC| >", 0, 3, 1, step = 0.25),
          shiny::selectInput(ns("cmp_pval"), "padj <",
                             c(0.001, 0.01, 0.05, 0.1), selected = 0.05),
          shiny::sliderInput(ns("cmp_n"), "Terms per contrast", 3, 30, 10),
          shiny::actionButton(ns("run_cmp"), "Compare",
                              class = "btn-primary w-100"),
          shiny::helpText("GO-BP over-representation across contrasts",
                          "(clusterProfiler::compareCluster).")
        ),
        plot_slot(ns("cmp_plot"), height = "720px")
      )
    )
  )
}

mod_enrichment_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- ORA ---------------------------------------------------------------
    shiny::observeEvent(input$run_ora, {
      shiny::req(state$obj)
      out <- run_step(state, "GO over-representation", {
        run_omni_enrichment(need_contrast(state),
                            lfc_thresh  = input$ora_lfc,
                            pval_thresh = as.numeric(input$ora_pval),
                            ont         = input$ora_ont,
                            simplify    = input$ora_simplify,
                            direction   = input$ora_direction)
      })
      if (!is.null(out)) state$obj <- out
    })

    ora_res <- shiny::reactive({
      obj <- need_contrast(state)
      res <- obj$enrichment_cache[[obj$active_contrast]]
      shiny::validate(shiny::need(res, paste0(
        "No ORA results for '", obj$active_contrast, "'. Click 'Run ORA'.")))
      res
    })

    # Dot plot: the package's interactive enrichment plot (hover shows the
    # term's genes). Network views (cnet / emap) are ggraph drawings and stay
    # images.
    output$ora_plot_out <- render_smart({
      res <- ora_res()
      if (identical(input$ora_plot, "dot") && pkgs_ok(c("plotly", "htmltools")))
        safe_widget(plot_omni_enrichment_interactive(
          res, n_top = input$ora_n, width = NULL, height = NULL,
          title = paste("GO", input$ora_ont, "--", state$obj$active_contrast)))[[1]]
      else
        safe_plot(plot_omni_ora(state$obj, plot_type = input$ora_plot,
                                show_n = input$ora_n))
    })
    output$ora_table <- DT::renderDT({
      df <- as_result_df(ora_res())
      shiny::validate(shiny::need(df, "No significant terms."))
      omni_dt(df, 10)
    })
    output$dl_ora <- csv_download("ORA_results.csv",
                                  function() as_result_df(ora_res()))

    # ---- ORA Up/Down ---------------------------------------------------------
    # run_omni_ora_updown() returns its results rather than storing them in
    # the object, so the app keeps them per contrast in state$ora_updown.
    shiny::observeEvent(input$run_ud, {
      obj <- need_contrast(state)
      shiny::req(length(input$ud_db) > 0)
      key <- obj$active_contrast
      out <- run_step(state, paste("ORA Up/Down:", key), {
        run_omni_ora_updown(obj,
                            pval_thresh = as.numeric(input$ud_pval),
                            lfc_thresh  = input$ud_lfc,
                            databases   = input$ud_db)
      })
      if (!is.null(out)) {
        if (length(out) == 0)
          shiny::showNotification("No significant terms in any database.",
                                  type = "warning")
        ud <- state$ora_updown
        ud[[key]] <- out
        state$ora_updown <- ud
      }
    })

    shiny::observeEvent(state$ora_updown, {
      done <- names(state$ora_updown)
      shiny::updateSelectInput(session, "ud_cmp", choices = done,
                               selected = .nz(shiny::isolate(input$ud_cmp)) %||% done)
    })

    # plot_omni_enrichment_interactive() returns tagList(plotly, datatable).
    # Rendered through renderUI, the datatable's filter JS (noUiSlider) is not
    # loaded and the resulting error also stops plotly drawing, so each
    # widget goes through its own output instead.
    ud_widgets <- shiny::reactive({
      obj <- need_contrast(state)
      res <- state$ora_updown[[obj$active_contrast]]
      shiny::validate(shiny::need(!is.null(res), paste0(
        "No Up/Down ORA for '", obj$active_contrast, "'. Click 'Run ORA'.")))
      shiny::validate(shiny::need(!is.null(res[[input$ud_show_db]]), paste0(
        input$ud_show_db, ": no significant terms (or not run) for this contrast.")))
      safe_widget(plot_omni_enrichment_interactive(
        res[[input$ud_show_db]], n_top = input$ud_n, width = NULL, height = 650,
        title = paste(input$ud_show_db, "--", obj$active_contrast)))
    })

    output$ud_plot  <- render_plotly(ud_widgets()[[1]])
    output$ud_table <- DT::renderDT({
      # The plot above already shows why nothing is available.
      w <- tryCatch(ud_widgets(), error = function(e) NULL)
      shiny::req(length(w) > 1)
      w[[2]]
    })

    output$ud_compare <- render_smart({
      keys <- input$ud_cmp
      shiny::validate(shiny::need(length(keys) > 0,
                                  "Run Up/Down ORA on one or more contrasts first."))
      res <- lapply(state$ora_updown[keys], function(x) x[[input$ud_show_db]])
      res <- res[!vapply(res, is.null, logical(1))]
      shiny::validate(shiny::need(length(res) > 0, paste0(
        "No ", input$ud_show_db, " terms in the selected contrasts.")))
      safe_plot(plot_omni_enrichment_comparison(
        res, n_top = min(input$ud_n, 15),
        title = paste(input$ud_show_db, "enrichment across contrasts")))
    })

    # ---- GSEA --------------------------------------------------------------
    shiny::observeEvent(input$run_gsea, {
      shiny::req(state$obj)
      out <- run_step(state, "GSEA", {
        run_omni_gsea(need_contrast(state),
                      category    = input$gsea_category,
                      subcategory = .nz(trimws(input$gsea_subcategory)),
                      pval_thresh = as.numeric(input$gsea_pval))
      })
      if (!is.null(out)) state$obj <- out
    })

    gsea_res <- shiny::reactive({
      obj <- need_contrast(state)
      res <- obj$gsea_cache[[obj$active_contrast]]
      shiny::validate(shiny::need(res, paste0(
        "No GSEA results for '", obj$active_contrast, "'. Click 'Run GSEA'.")))
      res
    })

    shiny::observe({
      df <- as_result_df(tryCatch(gsea_res(), error = function(e) NULL))
      ids <- if (is.null(df)) character(0) else df$ID
      shiny::updateSelectInput(session, "gsea_pathway", choices = ids,
                               selected = shiny::isolate(input$gsea_pathway))
    })

    output$gsea_plot_out <- render_smart({
      gsea_res()
      obj <- state$obj
      switch(input$gsea_plot,
             leading     = safe_plot(plot_omni_leading_edge(
               obj, n_pathways = min(input$gsea_n, 10))),
             pathway_lfc = safe_plot(plot_omni_pathway_lfc(
               obj, pathway_name = .nz(input$gsea_pathway))),
             safe_plot(plot_omni_gsea(obj, plot_type = input$gsea_plot,
                                      show_n = input$gsea_n)))
    })
    output$gsea_table <- DT::renderDT({
      df <- as_result_df(gsea_res())
      shiny::validate(shiny::need(df, "No significant gene sets."))
      df$core_enrichment <- NULL
      omni_dt(df, 10)
    })
    output$dl_gsea <- csv_download("GSEA_results.csv",
                                   function() as_result_df(gsea_res()))

    # ---- ssGSEA ------------------------------------------------------------
    shiny::observeEvent(input$run_ssgsea, {
      shiny::req(state$obj)
      out <- run_step(state, "ssGSEA", {
        run_omni_ssgsea(need_model(state),
                        category    = input$ss_category,
                        subcategory = .nz(trimws(input$ss_subcategory)))
      })
      if (!is.null(out)) state$obj <- out
    })

    ss_scores <- shiny::reactive({
      obj <- need_obj(state)
      shiny::validate(shiny::need(obj$advanced_ssgsea,
                                  "Click 'Run ssGSEA' (requires a fitted model)."))
      obj$advanced_ssgsea
    })

    output$ss_plot <- render_smart({
      ss_scores()
      safe_plot(plot_omni_ssgsea_heatmap(state$obj, show_n = input$ss_n,
                                         isolate_contrast = input$ss_isolate))
    })
    output$dl_ss <- csv_download("ssGSEA_scores.csv", function() {
      m <- ss_scores()
      cbind(Pathway = rownames(m), as.data.frame(m), row.names = NULL)
    })

    # ---- Compare contrasts --------------------------------------------------
    shiny::observeEvent(state$obj, {
      cs <- omni_contrasts(state$obj)
      shiny::updateSelectInput(session, "cmp_contrasts", choices = cs,
                               selected = .nz(shiny::isolate(input$cmp_contrasts)) %||%
                                 cs)
    })

    cmp_args <- shiny::eventReactive(input$run_cmp, {
      shiny::req(length(input$cmp_contrasts) >= 1)
      list(contrasts = input$cmp_contrasts, lfc = input$cmp_lfc,
           pval = as.numeric(input$cmp_pval), n = input$cmp_n)
    })

    output$cmp_plot <- render_smart({
      a   <- cmp_args()
      obj <- need_contrast(state)
      shiny::withProgress(message = "Comparing contrasts", value = 0.5,
        safe_plot(plot_omni_multicontrast(obj, contrasts = a$contrasts,
                                          lfc_thresh = a$lfc,
                                          pval_thresh = a$pval,
                                          show_n = a$n)))
    })
  })
}
