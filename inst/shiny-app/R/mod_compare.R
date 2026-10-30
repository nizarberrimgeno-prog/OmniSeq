# =============================================================================
# Compare contrasts module — two contrasts side by side
# (plot_omni_lfc_scatter_interactive) and formal interaction-term views
# (plot_omni_interaction_scatter / _static, plot_omni_interaction_dumbbell).
# =============================================================================

mod_compare_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 330,
      shiny::selectInput(ns("a"), "Contrast A", choices = NULL),
      shiny::textInput(ns("label_a"), "Label A (optional)"),
      shiny::selectInput(ns("b"), "Contrast B", choices = NULL),
      shiny::textInput(ns("label_b"), "Label B (optional)"),
      shiny::selectInput(ns("pval"), "padj <", c(0.001, 0.01, 0.05, 0.1),
                         selected = 0.05),
      shiny::sliderInput(ns("lfc"), "|log2FC| > (concordance)", 0, 3, 1,
                         step = 0.25),
      shiny::hr(),
      shiny::h6("Interaction term"),
      shiny::selectInput(ns("inter"), "Interaction contrast",
                         choices = NULL),
      shiny::helpText("A cached contrast holding the Wald test for the",
                      "A:B interaction coefficient, with contrasts A and B",
                      "the two simple effects it compares (e.g. stimulus",
                      "response in KO and in WT)."),
      shiny::sliderInput(ns("flip"), "Min |log2FC| to call a reversal",
                         0, 2, 0.3, step = 0.1),
      shiny::textInput(ns("stim"), "Stimulus label (dumbbell axis)"),
      shiny::selectInput(ns("rank_by"), "Dumbbell: rank genes by",
                         c("Interaction padj" = "padj", "Gap |B - A|" = "gap")),
      shiny::sliderInput(ns("n_top"), "Dumbbell: genes shown", 5, 50, 20)
    ),
    bslib::navset_card_tab(
      bslib::nav_panel(
        "log2FC concordance",
        deps_notice("plotly"),
        plotly_output(ns("concord"), height = "620px")
      ),
      bslib::nav_panel(
        "Interaction scatter",
        deps_notice("plotly"),
        plotly_output(ns("inter_int"), height = "620px")
      ),
      bslib::nav_panel(
        "Interaction scatter (static)",
        shiny::plotOutput(ns("inter_static"), height = "620px"),
        shiny::downloadButton(ns("dl_static"), "PNG", class = "btn-sm")
      ),
      bslib::nav_panel(
        "Interaction dumbbell",
        shiny::uiOutput(ns("gap_warning")),
        shiny::plotOutput(ns("dumbbell"), height = "700px"),
        shiny::downloadButton(ns("dl_dumbbell"), "PNG", class = "btn-sm")
      )
    )
  )
}

mod_compare_server <- function(id, state) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(state$obj, {
      cs <- omni_contrasts(state$obj)
      keep <- function(sel, default) {
        cur <- .nz(shiny::isolate(input[[sel]]))
        if (!is.null(cur) && cur %in% cs) cur else default
      }
      shiny::updateSelectInput(session, "a", choices = cs,
                               selected = keep("a", cs[1]))
      shiny::updateSelectInput(session, "b", choices = cs,
                               selected = keep("b", cs[min(2, length(cs))]))
      shiny::updateSelectInput(session, "inter", choices = c("(none)" = "", cs),
                               selected = keep("inter", ""))
    })

    pair <- shiny::reactive({
      obj <- need_contrast(state)
      shiny::validate(
        shiny::need(length(omni_contrasts(obj)) >= 2,
                    "Extract at least two contrasts on the DE tab."),
        shiny::need(nzchar(input$a %||% "") && nzchar(input$b %||% "") &&
                      input$a != input$b, "Pick two different contrasts.")
      )
      list(obj = obj, a = input$a, b = input$b,
           label_a = .nz(trimws(input$label_a)),
           label_b = .nz(trimws(input$label_b)),
           pval = as.numeric(input$pval))
    })

    inter <- shiny::reactive({
      p <- pair()
      shiny::validate(shiny::need(
        nzchar(input$inter %||% "") && !input$inter %in% c(p$a, p$b),
        "Choose an interaction contrast different from A and B."))
      c(p, interaction = input$inter)
    })

    output$concord <- render_plotly({
      p <- pair()
      safe_widget(plot_omni_lfc_scatter_interactive(
        p$obj, contrast_a = p$a, contrast_b = p$b,
        label_a = p$label_a, label_b = p$label_b,
        pval_thresh = p$pval, lfc_thresh = input$lfc,
        width = NULL, height = NULL))
    })

    output$inter_int <- render_plotly({
      p <- inter()
      safe_widget(plot_omni_interaction_scatter(
        p$obj, contrast_a = p$a, contrast_b = p$b,
        interaction_contrast = p$interaction,
        label_a = p$label_a, label_b = p$label_b,
        pval_thresh = p$pval, flip_lfc_thresh = input$flip,
        width = NULL, height = NULL))
    })

    static_fn <- function() {
      p <- inter()
      plot_omni_interaction_scatter_static(
        p$obj, contrast_a = p$a, contrast_b = p$b,
        interaction_contrast = p$interaction,
        label_a = p$label_a, label_b = p$label_b,
        pval_thresh = p$pval, flip_lfc_thresh = input$flip)
    }
    dumbbell_fn <- function() {
      p <- inter()
      plot_omni_interaction_dumbbell(
        p$obj, contrast_a = p$a, contrast_b = p$b,
        interaction_contrast = p$interaction,
        label_a = p$label_a, label_b = p$label_b,
        stim_label = input$stim %||% "", rank_by = input$rank_by,
        n_top = input$n_top)
    }

    output$inter_static <- shiny::renderPlot(safe_plot(static_fn()))
    output$dumbbell     <- shiny::renderPlot(safe_plot(dumbbell_fn()))

    # The package documents that ranking by gap on a null interaction picks
    # noise by construction; surface that next to the plot.
    output$gap_warning <- shiny::renderUI({
      if (!identical(input$rank_by, "gap")) return(NULL)
      shiny::div(class = "alert alert-warning small",
                 "Ranking by gap is only meaningful when the interaction has",
                 "significant genes; otherwise it selects noise. Prefer",
                 "'Interaction padj'.")
    })

    dl <- function(prefix, fn) shiny::downloadHandler(
      filename = function() paste0(prefix, ".png"),
      content  = function(file) {
        grDevices::png(file, width = 10, height = 8, units = "in", res = 300)
        on.exit(grDevices::dev.off())
        show_plot(fn())
      })
    output$dl_static   <- dl("interaction_scatter", static_fn)
    output$dl_dumbbell <- dl("interaction_dumbbell", dumbbell_fn)
  })
}
