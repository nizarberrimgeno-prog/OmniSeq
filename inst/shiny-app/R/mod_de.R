# =============================================================================
# Differential expression module — run_omni_de(), run_omni_contrast(),
# run_omni_mask(), and the DE visualisation family. All plots follow the
# active contrast chosen in the app-wide sidebar.
# =============================================================================

signature_choices <- c("ISGs (28)"            = "isgs_28",
                       "ISGs (extended)"      = "isgs_extended",
                       "Ribosomal"            = "ribosomal",
                       "Mitochondrial"        = "mitochondrial",
                       "Immunoglobulin"       = "immunoglobulin")

mod_de_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 330,
      bslib::accordion(
        open = c("Fit model", "Add contrast"),
        bslib::accordion_panel(
          "Fit model",
          shiny::selectInput(ns("shrinkage"), "LFC shrinkage",
                             c("apeglm", "ashr", "none")),
          shiny::actionButton(ns("run_de"), "Run DESeq2",
                              class = "btn-primary w-100"),
          shiny::helpText("Fits the design formula and extracts the last",
                          "model coefficient (condition vs reference).")
        ),
        bslib::accordion_panel(
          "Add contrast",
          shiny::selectInput(ns("contrast_col"), "Metadata column",
                             choices = NULL),
          shiny::selectInput(ns("group1"), "Group 1 (numerator)",
                             choices = NULL),
          shiny::selectInput(ns("group2"), "Group 2 (baseline)",
                             choices = NULL),
          shiny::selectInput(ns("contrast_shrinkage"), "LFC shrinkage",
                             c("ashr", "apeglm", "none")),
          shiny::actionButton(ns("run_contrast"), "Extract contrast",
                              class = "btn-outline-primary w-100"),
          shiny::helpText("A different column than the main condition",
                          "refits the model once and caches it.")
        ),
        bslib::accordion_panel(
          "Mask a signature",
          shiny::selectInput(ns("mask_sig"), "Built-in signature",
                             c("(none)" = "", signature_choices)),
          shiny::textAreaInput(ns("mask_genes"),
                               "Extra features (one per line)", rows = 3),
          shiny::textInput(ns("mask_label"), "Label (optional)"),
          shiny::actionButton(ns("run_mask"), "Mask active contrast",
                              class = "btn-outline-primary w-100"),
          shiny::helpText("Stores a masked copy of the active contrast;",
                          "the original stays available.")
        ),
        bslib::accordion_panel(
          "Plot settings",
          shiny::sliderInput(ns("lfc_thresh"), "|log2FC| threshold",
                             0, 5, 1, step = 0.25),
          shiny::selectInput(ns("pval_thresh"), "Adjusted p-value threshold",
                             c(0.001, 0.01, 0.05, 0.1), selected = 0.05),
          shiny::sliderInput(ns("top_n"), "Labelled features", 0, 50, 15),
          shiny::sliderInput(ns("heat_n"), "Heatmap: top features",
                             10, 200, 50, step = 10),
          shiny::textAreaInput(ns("heat_genes"),
                               "Heatmap: custom features (optional)", rows = 2),
          shiny::selectInput(ns("color_by"), "Colour samples by",
                             choices = NULL),
          shiny::selectInput(ns("shape_by"), "Shape samples by",
                             choices = NULL)
        )
      )
    ),

    shiny::uiOutput(ns("summary_boxes")),
    bslib::navset_card_tab(
      bslib::nav_panel(
        "Volcano",
        shiny::helpText("Hover a point for gene, log2FC and padj; drag to zoom,",
                        "double-click to reset. The PNG download is the labelled",
                        "static version."),
        plot_slot(ns("volcano"), height = "600px"),
        shiny::downloadButton(ns("dl_volcano"), "PNG", class = "btn-sm")
      ),
      bslib::nav_panel(
        "MA",
        plot_slot(ns("ma"), height = "560px"),
        shiny::downloadButton(ns("dl_ma"), "PNG", class = "btn-sm")
      ),
      bslib::nav_panel(
        "Heatmap",
        plot_slot(ns("heatmap"), height = "700px"),
        shiny::downloadButton(ns("dl_heatmap"), "PNG", class = "btn-sm")
      ),
      bslib::nav_panel(
        "PCA",
        plot_slot(ns("pca"), height = "560px"),
        shiny::downloadButton(ns("dl_pca"), "PNG", class = "btn-sm")
      ),
      bslib::nav_panel(
        "Top features",
        plot_slot(ns("top_genes"), height = "650px")
      ),
      bslib::nav_panel(
        "p-value histogram",
        plot_slot(ns("phist"), height = "450px")
      ),
      bslib::nav_panel(
        "Results table",
        shiny::checkboxInput(ns("sig_only"),
                             "Only features below the padj threshold", FALSE),
        DT::DTOutput(ns("table")),
        shiny::downloadButton(ns("dl_table"), "Download CSV", class = "btn-sm")
      )
    )
  )
}

mod_de_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {

    pval <- shiny::reactive(as.numeric(input$pval_thresh))

    # ---- Selector population -----------------------------------------------
    shiny::observeEvent(state$obj, {
      obj  <- state$obj
      cols <- group_columns(obj)
      shiny::updateSelectInput(session, "contrast_col", choices = cols,
                               selected = .nz(shiny::isolate(input$contrast_col)) %||%
                                 obj$params$condition)
      meta_cols <- names(omni_metadata(obj))
      shiny::updateSelectInput(session, "color_by",
                               choices = c("(condition)" = "", meta_cols),
                               selected = shiny::isolate(input$color_by))
      shiny::updateSelectInput(session, "shape_by",
                               choices = c("(none)" = "", meta_cols),
                               selected = shiny::isolate(input$shape_by))
    }, ignoreNULL = TRUE)

    shiny::observeEvent(list(input$contrast_col, state$obj), {
      obj <- state$obj
      shiny::req(obj, input$contrast_col)
      x  <- omni_metadata(obj)[[input$contrast_col]]
      lv <- if (is.factor(x)) levels(droplevels(x)) else
        sort(unique(as.character(stats::na.omit(x))))
      shiny::updateSelectInput(session, "group1", choices = lv,
                               selected = lv[min(2, length(lv))])
      shiny::updateSelectInput(session, "group2", choices = lv,
                               selected = lv[1])
    })

    # ---- Actions -----------------------------------------------------------
    shiny::observeEvent(input$run_de, {
      shiny::req(state$obj)
      out <- run_step(state, "Differential expression", {
        run_omni_de(state$obj, shrinkage = input$shrinkage)
      })
      if (!is.null(out)) state$obj <- out
    })

    shiny::observeEvent(input$run_contrast, {
      shiny::req(state$obj)
      if (identical(input$group1, input$group2)) {
        shiny::showNotification("Pick two different groups.", type = "error")
        return()
      }
      if (is.null(omni_model(state$obj)) &&
          identical(input$contrast_col, state$obj$params$condition)) {
        shiny::showNotification("Run DESeq2 (Fit model) first.", type = "error")
        return()
      }
      out <- run_step(state, paste("Contrast", input$group1, "vs", input$group2), {
        run_omni_contrast(state$obj,
                          column    = input$contrast_col,
                          group1    = input$group1,
                          group2    = input$group2,
                          shrinkage = input$contrast_shrinkage)
      })
      if (!is.null(out)) state$obj <- out
    })

    shiny::observeEvent(input$run_mask, {
      shiny::req(state$obj, length(omni_contrasts(state$obj)) > 0)
      genes <- trimws(strsplit(input$mask_genes %||% "", "[\n,;]")[[1]])
      genes <- genes[nzchar(genes)]
      if (!nzchar(input$mask_sig) && length(genes) == 0) {
        shiny::showNotification("Choose a signature or list features to mask.",
                                type = "error")
        return()
      }
      out <- run_step(state, "Mask signature", {
        run_omni_mask(state$obj,
                      features   = .nz(genes),
                      signature  = .nz(input$mask_sig),
                      mask_label = .nz(input$mask_label))
      })
      if (!is.null(out)) state$obj <- out
    })

    # ---- Outputs -----------------------------------------------------------
    res_df <- shiny::reactive({
      obj <- need_contrast(state)
      de_table(obj)
    })

    output$summary_boxes <- shiny::renderUI({
      df <- res_df()
      sig  <- !is.na(df$padj) & df$padj < pval()
      up   <- sum(sig & df$log2FoldChange >  input$lfc_thresh)
      down <- sum(sig & df$log2FoldChange < -input$lfc_thresh)
      # Values sized so they stay on one line in a quarter-width box.
      num <- function(x) shiny::span(style = "font-size: 1.6rem; white-space: nowrap;",
                                     format(x, big.mark = ","))
      bslib::layout_column_wrap(
        width = "180px", fill = FALSE,
        bslib::value_box("Active contrast",
                         shiny::span(style = "font-size: 1.05rem; word-break: break-word;",
                                     state$obj$active_contrast),
                         theme = "primary"),
        bslib::value_box("Tested", num(sum(!is.na(df$padj)))),
        bslib::value_box("Up", num(up), theme = "danger"),
        bslib::value_box("Down", num(down), theme = "info")
      )
    })

    volcano_fn <- function() plot_omni_volcano(need_contrast(state),
                                               lfc_thresh  = input$lfc_thresh,
                                               pval_thresh = pval(),
                                               top_n       = input$top_n)
    ma_fn <- function() plot_omni_ma(need_contrast(state),
                                     pval_thresh = pval(),
                                     top_n       = input$top_n)
    heatmap_fn <- function() {
      genes <- trimws(strsplit(input$heat_genes %||% "", "[\n,;]")[[1]])
      plot_omni_heatmap(need_contrast(state),
                        ntop         = input$heat_n,
                        pval_thresh  = pval(),
                        lfc_thresh   = input$lfc_thresh,
                        custom_genes = .nz(genes[nzchar(genes)]))
    }
    pca_fn <- function() plot_omni_pca(need_model(state),
                                       color_by = .nz(input$color_by),
                                       shape_by = .nz(input$shape_by))

    # Interactive volcano (hover = gene); falls back to the static,
    # labelled plot when plotly is not installed.
    output$volcano <- render_smart({
      obj <- need_contrast(state)
      if (pkgs_ok(c("plotly", "htmltools")))
        safe_widget(plot_omni_volcano_interactive(obj,
                                                  pval_thresh = pval(),
                                                  lfc_thresh  = input$lfc_thresh,
                                                  width = NULL, height = NULL))
      else safe_plot(volcano_fn())
    })
    output$ma        <- render_smart(safe_plot(ma_fn()))
    output$heatmap   <- render_smart(safe_plot(heatmap_fn()))
    output$pca       <- render_smart(safe_plot(pca_fn()))
    output$top_genes <- render_smart(safe_plot(
      plot_omni_top_genes(need_contrast(state),
                          n           = 12,
                          pval_thresh = pval(),
                          lfc_thresh  = input$lfc_thresh,
                          color_by    = .nz(input$color_by),
                          shape_by    = .nz(input$shape_by))))
    output$phist <- render_smart(safe_plot(
      plot_omni_pvalue_hist(need_contrast(state))))

    table_df <- shiny::reactive({
      df <- res_df()
      if (isTRUE(input$sig_only))
        df <- df[!is.na(df$padj) & df$padj < pval(), , drop = FALSE]
      df[order(df$padj), , drop = FALSE]
    })

    output$table <- DT::renderDT(safe_widget(
      omni_gene_search_table(need_contrast(state),
                             pval_thresh = if (isTRUE(input$sig_only)) pval(),
                             page_length = 15)))

    safe_name <- function() gsub("[^A-Za-z0-9_.-]+", "_",
                                 state$obj$active_contrast %||% "contrast")
    output$dl_table <- shiny::downloadHandler(
      filename = function() paste0("DE_", safe_name(), ".csv"),
      content  = function(file) utils::write.csv(table_df(), file,
                                                 row.names = FALSE)
    )
    dl_png <- function(prefix, fn) shiny::downloadHandler(
      filename = function() paste0(prefix, "_", safe_name(), ".png"),
      content  = function(file) {
        grDevices::png(file, width = 10, height = 8, units = "in", res = 300)
        on.exit(grDevices::dev.off())
        show_plot(fn())
      }
    )
    output$dl_volcano <- dl_png("volcano", volcano_fn)
    output$dl_ma      <- dl_png("MA", ma_fn)
    output$dl_heatmap <- dl_png("heatmap", heatmap_fn)
    output$dl_pca     <- dl_png("PCA", pca_fn)
  })
}
