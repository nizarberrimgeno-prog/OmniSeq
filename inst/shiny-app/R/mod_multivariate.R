# =============================================================================
# Multivariate module — WGCNA co-expression networks (single, per-condition,
# module preservation) and sPLS-DA (mixOmics).
# =============================================================================

mod_multivariate_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::navset_card_tab(
    bslib::nav_panel(
      "WGCNA",
      deps_notice(c("WGCNA", "pheatmap")),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 300,
          shiny::textInput(ns("w_key"), "Network name", "full"),
          shiny::numericInput(ns("w_ntop"), "Top variable genes", 5000,
                              min = 500, max = 20000, step = 500),
          shiny::selectInput(ns("w_type"), "Network type",
                             c("signed", "unsigned")),
          shiny::numericInput(ns("w_minsize"), "Min module size", 30, 5, 200),
          shiny::sliderInput(ns("w_merge"), "Merge cut height",
                             0.05, 0.5, 0.25, step = 0.05),
          shiny::checkboxInput(ns("w_masked"),
                               "Use masked gene set (active masked contrast)",
                               FALSE),
          shiny::actionButton(ns("run_wgcna"), "Build network",
                              class = "btn-primary w-100"),
          shiny::hr(),
          shiny::selectInput(ns("w_show"), "Show network", choices = NULL),
          shiny::selectInput(ns("w_traits"), "Traits (metadata columns)",
                             choices = NULL, multiple = TRUE),
          shiny::sliderInput(ns("w_hubs"), "Hub genes per module", 1, 15, 5)
        ),
        bslib::navset_underline(
          bslib::nav_panel("Dendrogram",
                           plot_slot(ns("w_dendro"), height = "600px")),
          bslib::nav_panel("Module-trait",
                           plot_slot(ns("w_trait"), height = "650px")),
          bslib::nav_panel("Hub network",
                           plot_slot(ns("w_hub"), height = "700px")),
          bslib::nav_panel("Module membership",
                           DT::DTOutput(ns("w_table")),
                           shiny::downloadButton(ns("dl_w"), "Download CSV",
                                                 class = "btn-sm"))
        )
      )
    ),

    bslib::nav_panel(
      "Module preservation",
      deps_notice(c("WGCNA")),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 300,
          shiny::selectInput(ns("mp_col"), "Group column", choices = NULL),
          shiny::selectInput(ns("mp_exclude"), "Exclude groups",
                             choices = NULL, multiple = TRUE),
          shiny::numericInput(ns("mp_ntop"), "Top variable genes", 5000,
                              500, 20000, 500),
          shiny::actionButton(ns("run_multi"), "1. Per-group networks",
                              class = "btn-primary w-100"),
          shiny::hr(),
          shiny::selectInput(ns("mp_ref"), "Reference group", choices = NULL),
          shiny::selectInput(ns("mp_test"), "Test group", choices = NULL),
          shiny::numericInput(ns("mp_perm"), "Permutations", 200, 50, 2000, 50),
          shiny::actionButton(ns("run_pres"), "2. Test preservation",
                              class = "btn-outline-primary w-100"),
          shiny::helpText("Zsummary > 10 = strongly preserved;",
                          "< 2 = condition-specific.")
        ),
        plot_slot(ns("mp_plot"), height = "600px")
      )
    ),

    bslib::nav_panel(
      "sPLS-DA",
      deps_notice(c("mixOmics")),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 300,
          shiny::selectInput(ns("s_col"), "Group column", choices = NULL),
          shiny::selectInput(ns("s_exclude"), "Exclude groups",
                             choices = NULL, multiple = TRUE),
          shiny::numericInput(ns("s_ncomp"), "Components (blank = tune)", NA,
                              1, 5),
          shiny::checkboxInput(ns("s_masked"), "Use masked gene set", FALSE),
          shiny::actionButton(ns("run_splsda"), "Fit sPLS-DA",
                              class = "btn-primary w-100"),
          shiny::helpText("Component number and keepX are tuned by",
                          "leave-one-out CV when left blank (slow on",
                          "large datasets)."),
          shiny::hr(),
          shiny::selectInput(ns("s_key"), "Model", choices = NULL),
          shiny::numericInput(ns("s_comp"), "Component", 1, 1, 5),
          shiny::sliderInput(ns("s_nfeat"), "Features shown", 5, 50, 20)
        ),
        bslib::navset_underline(
          bslib::nav_panel("Scores",
                           plot_slot(ns("s_scores"), height = "600px")),
          bslib::nav_panel("Loadings",
                           plot_slot(ns("s_load"), height = "600px")),
          bslib::nav_panel("VIP",
                           plot_slot(ns("s_vip"), height = "600px")),
          bslib::nav_panel("ROC",
                           plot_slot(ns("s_roc"), height = "600px")),
          bslib::nav_panel("Signature", DT::DTOutput(ns("s_sig")),
                           shiny::downloadButton(ns("dl_sig"), "Download CSV",
                                                 class = "btn-sm"))
        )
      )
    )
  )
}

mod_multivariate_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- Selectors ---------------------------------------------------------
    shiny::observeEvent(state$obj, {
      obj  <- state$obj
      cols <- group_columns(obj)
      cond <- obj$params$condition
      for (sel in c("mp_col", "s_col"))
        shiny::updateSelectInput(session, sel, choices = cols,
                                 selected = .nz(shiny::isolate(input[[sel]])) %||% cond)
      shiny::updateSelectInput(session, "w_traits",
                               choices = names(omni_metadata(obj)),
                               selected = shiny::isolate(input$w_traits))
      shiny::updateSelectInput(session, "w_show", choices = names(obj$wgcna),
                               selected = shiny::isolate(input$w_show))
      shiny::updateSelectInput(session, "s_key", choices = names(obj$splsda),
                               selected = shiny::isolate(input$s_key))
    })

    levels_of <- function(col) {
      shiny::req(state$obj, col)
      sort(unique(as.character(stats::na.omit(omni_metadata(state$obj)[[col]]))))
    }
    shiny::observeEvent(input$mp_col, {
      lv <- levels_of(input$mp_col)
      shiny::updateSelectInput(session, "mp_exclude", choices = lv)
      shiny::updateSelectInput(session, "mp_ref", choices = lv, selected = lv[1])
      shiny::updateSelectInput(session, "mp_test", choices = lv,
                               selected = lv[min(2, length(lv))])
    })
    shiny::observeEvent(input$s_col, {
      shiny::updateSelectInput(session, "s_exclude", choices = levels_of(input$s_col))
    })

    # ---- WGCNA -------------------------------------------------------------
    shiny::observeEvent(input$run_wgcna, {
      shiny::req(state$obj)
      key <- .nz(trimws(input$w_key)) %||% "full"
      out <- run_step(state, paste("WGCNA network", key), {
        run_omni_wgcna(need_model(state),
                       use_masked       = input$w_masked,
                       n_top_genes      = as.integer(input$w_ntop),
                       network_type     = input$w_type,
                       min_module_size  = as.integer(input$w_minsize),
                       merge_cut_height = input$w_merge,
                       key              = key)
      })
      if (!is.null(out)) {
        state$obj <- out
        shiny::updateSelectInput(session, "w_show", choices = names(out$wgcna),
                                 selected = key)
      }
    })

    w_obj <- shiny::reactive({
      obj <- need_obj(state)
      shiny::validate(shiny::need(length(obj$wgcna) > 0 && nzchar(input$w_show),
                                  "Click 'Build network' (requires a fitted model)."))
      obj
    })

    output$w_dendro <- render_smart(safe_plot(
      plot_omni_wgcna_dendrogram(w_obj(), key = input$w_show)))
    # Selected metadata columns -> numeric trait table (factors as level
    # codes). NULL lets the package build its default trait set.
    w_traits <- shiny::reactive({
      cols <- .nz(input$w_traits)
      if (is.null(cols)) return(NULL)
      obj  <- w_obj()
      meta <- omni_metadata(obj)[rownames(obj$wgcna[[input$w_show]]$dat_expr),
                                 cols, drop = FALSE]
      meta[] <- lapply(meta, function(x)
        if (is.numeric(x)) x else as.numeric(as.factor(x)))
      meta
    })

    output$w_trait <- render_smart(safe_plot(
      plot_omni_wgcna_trait_heatmap(w_obj(), key = input$w_show,
                                    traits = w_traits())))
    output$w_hub <- render_smart(safe_plot(
      plot_omni_wgcna_hub_network(w_obj(), key = input$w_show,
                                  n_hubs = input$w_hubs,
                                  plot_type = network_type())))

    w_members <- shiny::reactive({
      entry <- w_obj()$wgcna[[input$w_show]]
      shiny::req(entry)
      data.frame(Feature = names(entry$module_colors) %||%
                   colnames(entry$dat_expr),
                 Module  = unname(entry$module_colors))
    })
    output$w_table <- DT::renderDT(omni_dt(w_members()))
    output$dl_w <- csv_download("wgcna_modules.csv", w_members)

    # ---- Module preservation ----------------------------------------------
    shiny::observeEvent(input$run_multi, {
      shiny::req(state$obj)
      out <- run_step(state, "Per-group WGCNA", {
        run_omni_wgcna_multi(need_model(state),
                             group_col      = input$mp_col,
                             exclude_groups = .nz(input$mp_exclude),
                             n_top_genes    = as.integer(input$mp_ntop))
      })
      if (!is.null(out)) state$obj <- out
    })

    shiny::observeEvent(input$run_pres, {
      shiny::req(state$obj)
      out <- run_step(state, "Module preservation", {
        run_omni_module_preservation(state$obj,
                                     group_col      = input$mp_col,
                                     reference      = input$mp_ref,
                                     test           = input$mp_test,
                                     n_permutations = as.integer(input$mp_perm))
      })
      if (!is.null(out)) state$obj <- out
    })

    output$mp_plot <- render_smart({
      obj <- need_obj(state)
      shiny::validate(shiny::need(
        !is.null(obj$wgcna_multi[[input$mp_col]]$preservation),
        "Run steps 1 and 2 to test module preservation."))
      safe_plot(plot_omni_module_preservation(obj, group_col = input$mp_col))
    })

    # ---- sPLS-DA -----------------------------------------------------------
    shiny::observeEvent(input$run_splsda, {
      shiny::req(state$obj)
      ncomp <- if (is.na(input$s_ncomp)) NULL else as.integer(input$s_ncomp)
      out <- run_step(state, "sPLS-DA", {
        run_omni_splsda(need_model(state),
                        group_col      = input$s_col,
                        exclude_groups = .nz(input$s_exclude),
                        ncomp          = ncomp,
                        use_masked     = input$s_masked)
      })
      if (!is.null(out)) {
        state$obj <- out
        keys <- names(out$splsda)
        shiny::updateSelectInput(session, "s_key", choices = keys,
                                 selected = keys[length(keys)])
      }
    })

    s_obj <- shiny::reactive({
      obj <- need_obj(state)
      shiny::validate(shiny::need(length(obj$splsda) > 0 && nzchar(input$s_key),
                                  "Click 'Fit sPLS-DA' (requires a fitted model)."))
      obj
    })
    s_comp <- shiny::reactive(as.integer(input$s_comp %||% 1))

    output$s_scores <- render_smart(safe_plot(
      plot_omni_splsda_scores(s_obj(), key = input$s_key)))
    output$s_load <- render_smart(safe_plot(
      plot_omni_splsda_loadings(s_obj(), key = input$s_key, comp = s_comp(),
                                n_features = input$s_nfeat)))
    output$s_vip <- render_smart(safe_plot(
      plot_omni_splsda_vip(s_obj(), key = input$s_key, comp = s_comp(),
                           n_features = input$s_nfeat)))
    output$s_roc <- render_smart(safe_plot(
      plot_omni_splsda_roc(s_obj(), key = input$s_key, comp = s_comp())))

    s_sig <- shiny::reactive(
      omni_splsda_signature(s_obj(), key = input$s_key, comp = s_comp()))
    output$s_sig <- DT::renderDT(omni_dt(s_sig()))
    output$dl_sig <- csv_download("splsda_signature.csv", s_sig)
  })
}
