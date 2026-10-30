# =============================================================================
# Multi-omics module — pairs the main object with a second OmniSeq object
# from the same samples (gene + TE) for innate-immune sensor scoring,
# TE-gene correlation, WGCNA TE traits, DIABLO and causal mediation.
# Results are written to the gene-level object, as the package does.
# =============================================================================

mod_multiomics_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 330,
      bslib::accordion(
        open = "Paired object",
        bslib::accordion_panel(
          "Paired object",
          shiny::radioButtons(ns("pair_src"), NULL,
                              c("Count matrix (same samples)" = "counts",
                                "Saved OmniSeq object (.rds)" = "rds")),
          shiny::conditionalPanel(
            sprintf("input['%s'] == 'counts'", ns("pair_src")),
            shiny::fileInput(ns("pair_counts"), "Count matrix (CSV / TSV)",
                             accept = c(".csv", ".tsv", ".txt")),
            shiny::helpText("Built with the main object's metadata and design,",
                            "then DESeq2 is fitted automatically.")
          ),
          shiny::conditionalPanel(
            sprintf("input['%s'] == 'rds'", ns("pair_src")),
            shiny::fileInput(ns("pair_rds"), "OmniSeq object (.rds)",
                             accept = ".rds")
          ),
          shiny::actionButton(ns("load_pair"), "Load paired object",
                              class = "btn-primary w-100"),
          shiny::verbatimTextOutput(ns("pair_info"))
        ),
        bslib::accordion_panel(
          "Sensor scores",
          shiny::selectInput(ns("sens_col"), "Group column", choices = NULL),
          shiny::selectInput(ns("sens_excl"), "Exclude groups", choices = NULL,
                             multiple = TRUE),
          shiny::numericInput(ns("sens_min"), "Min loci per family", 5, 1, 100),
          shiny::helpText("Sensor TE families (comma-separated). Leave blank",
                          "for the package defaults (RepeatMasker names);",
                          "set them when your TE annotation names families",
                          "differently, e.g. Telescope's HML2 for HERVK."),
          shiny::textInput(ns("fam_cgas"), "cGAS-STING families"),
          shiny::textInput(ns("fam_rigi"), "RIG-I / MDA5 families"),
          shiny::textInput(ns("fam_tlr7"), "TLR7 / TLR8 families"),
          shiny::uiOutput(ns("fam_hint")),
          shiny::actionButton(ns("run_sens"), "Score sensors",
                              class = "btn-outline-primary w-100")
        ),
        bslib::accordion_panel(
          "DIABLO",
          deps_notice("mixOmics"),
          shiny::selectInput(ns("dia_col"), "Group column", choices = NULL),
          shiny::selectInput(ns("dia_excl"), "Exclude groups", choices = NULL,
                             multiple = TRUE),
          shiny::numericInput(ns("dia_ncomp"), "Components", 2, 1, 5),
          shiny::actionButton(ns("run_dia"), "Run DIABLO",
                              class = "btn-outline-primary w-100")
        ),
        bslib::accordion_panel(
          "Mediation",
          deps_notice("mediation"),
          shiny::selectInput(ns("med_mediator"), "Mediator",
                             c("cGAS_STING", "RIG_I_MDA5", "TLR7_TLR8")),
          shiny::selectInput(ns("med_treat"), "Treatment column", choices = NULL),
          shiny::selectInput(ns("med_g0"), "Treatment = 0", choices = NULL),
          shiny::selectInput(ns("med_g1"), "Treatment = 1", choices = NULL),
          shiny::selectInput(ns("med_outcome"), "Outcome", "IRG_S_Score"),
          shiny::numericInput(ns("med_boot"), "Bootstrap simulations", 500,
                              100, 5000, 100),
          shiny::checkboxInput(ns("med_all"), "Test all three sensors", FALSE),
          shiny::actionButton(ns("run_med"), "Run mediation",
                              class = "btn-outline-primary w-100")
        ),
        bslib::accordion_panel(
          "TE-gene correlation",
          shiny::textAreaInput(ns("corr_genes"), "Marker genes (one per line)",
                               "IFI27\nIFI44L\nISG15\nMX1\nRSAD2\nSIGLEC1",
                               rows = 4),
          shiny::selectInput(ns("corr_class"), "TE class",
                             c("ALL", "LINE", "SINE", "LTR", "ERV", "DNA")),
          shiny::textInput(ns("corr_name"), "Score name", "Marker Score")
        ),
        bslib::accordion_panel(
          "WGCNA TE traits",
          deps_notice("WGCNA"),
          shiny::selectInput(ns("wt_key"), "WGCNA network (gene object)",
                             choices = NULL),
          shiny::numericInput(ns("wt_nfam"), "Top TE families", 50, 5, 200),
          shiny::actionButton(ns("run_wt"), "Correlate modules with TEs",
                              class = "btn-outline-primary w-100")
        )
      )
    ),
    bslib::navset_card_tab(
      bslib::nav_panel("Sensor scores",
                       shiny::plotOutput(ns("sens_plot"), height = "520px"),
                       shiny::downloadButton(ns("dl_sens"), "Download CSV",
                                             class = "btn-sm")),
      bslib::nav_panel("Sensor heatmap",
                       shiny::plotOutput(ns("sens_heat"), height = "600px")),
      bslib::nav_panel(
        "Sensor vs IRGs",
        shiny::fluidRow(
          shiny::column(6, shiny::selectInput(ns("irg_color"), "Colour by",
                                              choices = NULL)),
          shiny::column(6, shiny::selectInput(ns("irg_shape"), "Shape by",
                                              choices = NULL))
        ),
        shiny::plotOutput(ns("sens_irg"), height = "560px")
      ),
      bslib::nav_panel("TE-gene correlation",
                       shiny::plotOutput(ns("corr_plot"), height = "600px")),
      bslib::nav_panel("WGCNA TE traits",
                       shiny::plotOutput(ns("wt_plot"), height = "650px")),
      bslib::nav_panel("DIABLO",
                       shiny::plotOutput(ns("dia_plot"), height = "560px"),
                       shiny::verbatimTextOutput(ns("dia_info"))),
      bslib::nav_panel("Mediation",
                       shiny::plotOutput(ns("med_plot"), height = "480px"),
                       shiny::plotOutput(ns("med_path"), height = "420px"))
    )
  )
}

mod_multiomics_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {

    # gene/TE roles: whichever object is TE-level is the TE partner.
    pair <- shiny::reactive({
      main <- need_obj(state)
      shiny::validate(shiny::need(state$partner,
                                  "Load a paired object (gene + TE) first."))
      main_is_te <- identical(main$params$data_type, "TE")
      part_is_te <- identical(state$partner$params$data_type, "TE")
      shiny::validate(shiny::need(main_is_te != part_is_te,
        "Need one gene-level and one TE-level object; both are the same type."))
      if (main_is_te)
        list(gene = state$partner, te = main, gene_is_main = FALSE)
      else
        list(gene = main, te = state$partner, gene_is_main = TRUE)
    })

    set_gene <- function(gene_obj) {
      if (pair()$gene_is_main) state$obj <- gene_obj else state$partner <- gene_obj
    }

    # ---- Paired object -----------------------------------------------------
    shiny::observeEvent(input$load_pair, {
      main <- state$obj
      shiny::req(main)
      if (is.null(omni_model(main))) {
        shiny::showNotification("Run DE on the main object first.", type = "error")
        return()
      }
      out <- run_step(state, "Load paired object", {
        if (input$pair_src == "rds") {
          shiny::req(input$pair_rds)
          validate_omniseq(readRDS(input$pair_rds$datapath))
        } else {
          shiny::req(input$pair_counts)
          counts <- read_table_upload(input$pair_counts$datapath,
                                      input$pair_counts$name)
          meta   <- omni_metadata(main)
          common <- intersect(colnames(counts), rownames(meta))
          if (length(common) < 3)
            stop("Fewer than 3 samples shared with the main object.",
                 call. = FALSE)
          obj <- init_omniseq(as.matrix(counts[, common]),
                              metadata      = meta[common, , drop = FALSE],
                              condition_col = main$params$condition,
                              covariates    = .nz(main$params$covariates),
                              species       = main$params$species,
                              org_db        = main$params$org_db)
          run_omni_de(obj, shrinkage = "none")
        }
      })
      if (!is.null(out)) state$partner <- out
    })

    output$pair_info <- shiny::renderPrint({
      shiny::req(state$partner)
      p <- state$partner
      cat("Paired:", p$params$data_type, "|", nrow(p$data), "features x",
          ncol(p$data), "samples\n")
    })

    shiny::observeEvent(list(state$obj, state$partner), {
      shiny::req(state$obj)
      cols <- group_columns(state$obj)
      cond <- state$obj$params$condition
      for (sel in c("sens_col", "dia_col", "med_treat"))
        shiny::updateSelectInput(session, sel, choices = cols,
                                 selected = .nz(shiny::isolate(input[[sel]])) %||% cond)
      meta_cols <- names(omni_metadata(state$obj))
      shiny::updateSelectInput(session, "irg_color", choices = meta_cols,
                               selected = .nz(shiny::isolate(input$irg_color)) %||% cond)
      shiny::updateSelectInput(session, "irg_shape",
                               choices = c("(none)" = "", meta_cols),
                               selected = shiny::isolate(input$irg_shape))
      gene <- if (identical(state$obj$params$data_type, "TE")) state$partner else
        state$obj
      shiny::updateSelectInput(session, "wt_key", choices = names(gene$wgcna))
    })

    lv <- function(col) {
      shiny::req(state$obj, col)
      sort(unique(as.character(stats::na.omit(omni_metadata(state$obj)[[col]]))))
    }
    shiny::observeEvent(input$sens_col,
      shiny::updateSelectInput(session, "sens_excl", choices = lv(input$sens_col)))
    shiny::observeEvent(input$dia_col,
      shiny::updateSelectInput(session, "dia_excl", choices = lv(input$dia_col)))
    shiny::observeEvent(input$med_treat, {
      g <- lv(input$med_treat)
      shiny::updateSelectInput(session, "med_g0", choices = g, selected = g[1])
      shiny::updateSelectInput(session, "med_g1", choices = g,
                               selected = g[min(2, length(g))])
    })

    # ---- Sensor scores -----------------------------------------------------
    shiny::observeEvent(input$run_sens, {
      p <- pair()
      out <- run_step(state, "Sensor scoring", {
        run_omni_sensor_score(p$gene, te_obj = p$te,
                              group_col      = input$sens_col,
                              exclude_groups = .nz(input$sens_excl),
                              min_loci       = as.integer(input$sens_min),
                              cgas_families  = split_families(input$fam_cgas),
                              rigi_families  = split_families(input$fam_rigi),
                              tlr7_families  = split_families(input$fam_tlr7))
      })
      if (!is.null(out)) set_gene(out)
    })

    split_families <- function(x) {
      f <- trimws(strsplit(x %||% "", "[,;\n]")[[1]])
      .nz(f[nzchar(f)])
    }

    # Most abundant TE families in the paired TE object, to help fill in the
    # sensor family fields with names that actually occur in the data.
    output$fam_hint <- shiny::renderUI({
      te <- tryCatch(pair()$te, error = function(e) NULL)
      if (is.null(te)) return(NULL)
      fam <- sort(table(sub("^[^:]+:[^:]+:", "", rownames(te$data))),
                  decreasing = TRUE)
      shiny::helpText(shiny::strong("Families in your TE data: "),
                      paste(utils::head(names(fam), 25), collapse = ", "))
    })

    sens_obj <- shiny::reactive({
      g <- pair()$gene
      shiny::validate(shiny::need(g$sensor_scores, "Click 'Score sensors'."))
      g
    })

    shiny::observeEvent(sens_obj(), {
      sc <- sens_obj()$sensor_scores$scores
      num <- names(sc)[vapply(sc, is.numeric, logical(1))]
      shiny::updateSelectInput(session, "med_mediator", choices = num,
                               selected = shiny::isolate(input$med_mediator))
      shiny::updateSelectInput(session, "med_outcome", choices = num,
                               selected = if ("IRG_S_Score" %in% num) "IRG_S_Score")
    })

    output$sens_plot <- shiny::renderPlot(safe_plot(
      plot_omni_sensor_scores(sens_obj(), group_col = input$sens_col)))
    output$sens_heat <- shiny::renderPlot(safe_plot(
      plot_omni_sensor_heatmap(sens_obj())))
    output$sens_irg <- shiny::renderPlot(safe_plot(
      plot_omni_sensor_irgs_corr(sens_obj(), color_col = input$irg_color,
                                 shape_col = .nz(input$irg_shape))))
    output$dl_sens <- csv_download("sensor_scores.csv",
                                   function() sens_obj()$sensor_scores$scores)

    # ---- TE-gene correlation ----------------------------------------------
    output$corr_plot <- shiny::renderPlot({
      p <- pair()
      genes <- trimws(strsplit(input$corr_genes %||% "", "[\n,;]")[[1]])
      genes <- genes[nzchar(genes)]
      shiny::validate(shiny::need(length(genes) > 0, "Enter marker genes."))
      safe_plot(plot_omni_te_gene_corr(p$gene, p$te, marker_genes = genes,
                                       te_class = input$corr_class,
                                       score_name = input$corr_name))
    })

    # ---- WGCNA TE traits ---------------------------------------------------
    shiny::observeEvent(input$run_wt, {
      p <- pair()
      shiny::req(nzchar(input$wt_key %||% ""))
      out <- run_step(state, "WGCNA TE traits", {
        run_omni_wgcna_te_traits(p$gene, p$te, wgcna_key = input$wt_key,
                                 n_top_families = as.integer(input$wt_nfam))
      })
      if (!is.null(out)) set_gene(out)
    })

    output$wt_plot <- shiny::renderPlot({
      g <- pair()$gene
      shiny::validate(shiny::need(
        nzchar(input$wt_key %||% "") && !is.null(g$wgcna[[input$wt_key]]$te_trait_cor),
        "Build a WGCNA network on the gene object (Multivariate tab), then click 'Correlate modules with TEs'."))
      safe_plot(plot_omni_wgcna_te_trait_heatmap(g, key = input$wt_key))
    })

    # ---- DIABLO ------------------------------------------------------------
    shiny::observeEvent(input$run_dia, {
      p <- pair()
      out <- run_step(state, "DIABLO", {
        run_omni_diablo(p$gene, p$te,
                        group_col      = input$dia_col,
                        exclude_groups = .nz(input$dia_excl),
                        ncomp          = as.integer(input$dia_ncomp))
      })
      if (!is.null(out)) set_gene(out)
    })

    dia_entry <- shiny::reactive({
      g <- pair()$gene
      shiny::validate(shiny::need(length(g$diablo) > 0, "Click 'Run DIABLO'."))
      g$diablo[[length(g$diablo)]]
    })
    output$dia_plot <- shiny::renderPlot({
      e <- dia_entry()
      shiny::validate(shiny::need(requireNamespace("mixOmics", quietly = TRUE),
                                  "mixOmics is not installed."))
      mixOmics::plotIndiv(e$model, ind.names = FALSE, legend = TRUE,
                          ellipse = TRUE, title = "DIABLO sample scores")
    })
    output$dia_info <- shiny::renderPrint({
      e <- dia_entry()
      cat("Model:", names(pair()$gene$diablo)[length(pair()$gene$diablo)], "\n")
      utils::str(e[setdiff(names(e), "model")], max.level = 1, give.attr = FALSE)
    })

    # ---- Mediation ---------------------------------------------------------
    shiny::observeEvent(input$run_med, {
      g <- sens_obj()
      if (identical(input$med_g0, input$med_g1)) {
        shiny::showNotification("Pick two different treatment groups.",
                                type = "error")
        return()
      }
      scores <- g$sensor_scores$scores
      tv  <- stats::setNames(c(input$med_g0, input$med_g1), c("0", "1"))
      excl <- setdiff(unique(as.character(scores[[input$med_treat]])),
                      c(input$med_g0, input$med_g1))
      out <- run_step(state, "Mediation", {
        if (isTRUE(input$med_all))
          run_omni_mediation_multi(g, sensor_scores = scores,
                                   treatment_col  = input$med_treat,
                                   treatment_vals = tv,
                                   exclude_groups = .nz(excl),
                                   outcome        = input$med_outcome,
                                   boot_sims      = as.integer(input$med_boot))
        else
          run_omni_mediation(g, sensor_scores = scores,
                             mediator       = input$med_mediator,
                             treatment_col  = input$med_treat,
                             treatment_vals = tv,
                             exclude_groups = .nz(excl),
                             outcome        = input$med_outcome,
                             boot_sims      = as.integer(input$med_boot))
      })
      if (!is.null(out)) set_gene(out)
    })

    med_obj <- shiny::reactive({
      g <- pair()$gene
      shiny::validate(shiny::need(length(g$mediation) > 0,
                                  "Score sensors, then click 'Run mediation'."))
      g
    })
    output$med_plot <- shiny::renderPlot(safe_plot(plot_omni_mediation(med_obj())))
    output$med_path <- shiny::renderPlot(safe_plot(
      plot_omni_mediation_path(med_obj())))
  })
}
