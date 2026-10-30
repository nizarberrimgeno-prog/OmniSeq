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
        shiny::plotOutput(ns("ora_plot_out"), height = "620px"),
        DT::DTOutput(ns("ora_table")),
        shiny::downloadButton(ns("dl_ora"), "Download CSV", class = "btn-sm")
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
        shiny::plotOutput(ns("gsea_plot_out"), height = "620px"),
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
        shiny::plotOutput(ns("ss_plot"), height = "720px"),
        shiny::downloadButton(ns("dl_ss"), "Download scores (CSV)",
                              class = "btn-sm")
      )
    ),

    bslib::nav_panel(
      "Compare contrasts",
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
        shiny::plotOutput(ns("cmp_plot"), height = "720px")
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

    output$ora_plot_out <- shiny::renderPlot({
      ora_res()
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

    output$gsea_plot_out <- shiny::renderPlot({
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

    output$ss_plot <- shiny::renderPlot({
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

    output$cmp_plot <- shiny::renderPlot({
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
