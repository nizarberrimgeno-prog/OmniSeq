# =============================================================================
# OmniSeq Shiny App — shared helpers
#
# Files in inst/shiny-app/R/ are sourced automatically by shiny::runApp()
# before app.R. Nothing here performs analysis — these helpers only run
# package functions safely and render their output.
# =============================================================================

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Character vector of choices -> NULL when empty (selectInput returns "").
.nz <- function(x) if (is.null(x) || length(x) == 0 || all(!nzchar(x))) NULL else x

# TRUE if every package in `pkgs` can be loaded.
pkgs_ok <- function(pkgs) {
  all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
}

# A UI notice listing missing packages, or NULL when everything is installed.
# Used at the top of each optional-module tab.
deps_notice <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(NULL)
  bslib::card(
    class = "border-warning",
    bslib::card_body(
      shiny::strong("Missing packages: "),
      paste(missing, collapse = ", "),
      shiny::br(),
      "This module is unavailable until they are installed, e.g. ",
      shiny::code(sprintf("BiocManager::install(c(%s))",
                          paste0('"', missing, '"', collapse = ", ")))
    )
  )
}

# Run one analysis step: progress indicator, message capture into the shared
# log, warnings/errors as notifications. Returns the step's value, or NULL
# on error (callers keep the previous object in that case).
run_step <- function(state, label, expr) {
  log_line <- function(txt) {
    txt <- gsub("^\n+|\n+$", "", txt)
    if (nzchar(txt))
      state$log <- c(state$log, paste0(format(Sys.time(), "%H:%M:%S"),
                                       "  ", txt))
  }
  log_line(paste0(">>> ", label))

  shiny::withProgress(message = label, value = 0.3, {
    out <- tryCatch(
      withCallingHandlers(
        expr,
        message = function(m) {
          log_line(conditionMessage(m))
          invokeRestart("muffleMessage")
        },
        warning = function(w) {
          log_line(paste("Warning:", conditionMessage(w)))
          shiny::showNotification(conditionMessage(w), type = "warning",
                                  duration = 8)
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        log_line(paste("Error:", conditionMessage(e)))
        shiny::showNotification(paste0(label, " failed: ", conditionMessage(e)),
                                type = "error", duration = NULL)
        NULL
      }
    )
    shiny::setProgress(1)
  })

  if (!is.null(out))
    shiny::showNotification(paste(label, "complete."), type = "message")
  out
}

# Draw whatever a plot_omni_*() function produced. ggplot/patchwork objects
# are returned invisibly by the package, so renderPlot() would not print them
# on its own; pheatmap and base-graphics plots have already drawn themselves.
show_plot <- function(p) {
  if (inherits(p, c("gg", "ggplot", "patchwork")))
    print(p)
  invisible(p)
}

# Wrap a plot call so missing results produce a readable message instead of
# a red error in the plot area.
safe_plot <- function(expr) {
  tryCatch(
    show_plot(suppressMessages(expr)),
    error = function(e) shiny::validate(shiny::need(FALSE, conditionMessage(e)))
  )
}

# Same as safe_plot() for htmlwidget-returning functions (plotly, DT).
safe_widget <- function(expr) {
  tryCatch(
    suppressMessages(expr),
    error = function(e) shiny::validate(shiny::need(FALSE, conditionMessage(e)))
  )
}

# plotly is a suggested dependency: only build plotly outputs/renderers when
# it is installed, so the app still starts without it.
plotly_output <- function(id, height = "550px") {
  if (pkgs_ok("plotly")) plotly::plotlyOutput(id, height = height)
}

render_plotly <- function(expr, env = parent.frame()) {
  if (!pkgs_ok("plotly")) return(shiny::renderUI(NULL))
  plotly::renderPlotly(substitute(expr), env = env, quoted = TRUE)
}

# ---- Interactive plot slots ----------------------------------------------------
# Every plot in the app goes through plot_slot() / render_smart():
#   - ggplot results are converted to plotly (hover, zoom, pan); ggrepel
#     label layers become plain text so labels survive the conversion
#   - htmlwidgets (plotly, visNetwork) are shown as they are
#   - anything else (pheatmap, base graphics, patchwork, ggraph) is drawn to
#     a PNG and shown as an image
# Without plotly installed every plot falls back to the image.

plot_slot <- function(id, height = "450px") {
  shiny::div(class = "omni-plot-slot", style = sprintf("height: %s;", height),
             shiny::uiOutput(id, style = "height: 100%;"))
}

.omni_static_classes <- c("patchwork", "ggraph", "ggmatrix")

# Network plots (STRING, TF network, WGCNA hubs) have an interactive
# visNetwork mode in the package; use it when visNetwork is installed.
network_type <- function() if (pkgs_ok("visNetwork")) "interactive" else "static"

.omni_repel_to_text <- function(p) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, c("GeomTextRepel", "GeomLabelRepel"))) {
      keep <- intersect(names(p$layers[[i]]$aes_params),
                        c("size", "colour", "alpha", "family", "fontface"))
      p$layers[[i]]$geom        <- ggplot2::GeomText
      p$layers[[i]]$geom_params <- list(na.rm = TRUE)
      p$layers[[i]]$aes_params  <- p$layers[[i]]$aes_params[keep]
    }
  }
  p
}

# Multi-panel ggplot results (patchwork, enrichplot's gglist) -> plotly
# subplots, one converted panel at a time.
.omni_panels_to_plotly <- function(val) {
  panels <- if (inherits(val, "patchwork")) lapply(seq_along(val), function(i) val[[i]])
            else unclass(val)
  if (length(panels) < 2 || !all(vapply(panels, inherits, logical(1), "ggplot")))
    return(NULL)
  ws <- lapply(panels, function(p) plotly::ggplotly(.omni_repel_to_text(p)))
  # One legend: the panels usually share their colour mapping.
  ws[-length(ws)] <- lapply(ws[-length(ws)], plotly::style, showlegend = FALSE)
  stacked <- inherits(val, "gglist")  # enrichplot running-score layout
  out <- plotly::subplot(ws, nrows = if (stacked) length(ws) else 1,
                         shareX = stacked, titleX = TRUE, titleY = TRUE,
                         margin = 0.05)
  if (stacked) return(out)
  # Side-by-side panels: subplot keeps only one title, so put each panel's
  # own title above it.
  titles <- vapply(panels, function(p) p$labels$title %||% "", character(1))
  n <- length(panels)
  plotly::layout(out, title = list(text = ""), margin = list(t = 60),
                 annotations = lapply(seq_len(n), function(i) list(
                   text = titles[i], x = (i - 0.5) / n, y = 1.04,
                   xref = "paper", yref = "paper", showarrow = FALSE,
                   xanchor = "center", yanchor = "bottom", font = list(size = 14))))
}

# pheatmap only draws; render_smart() traces pheatmap::pheatmap() to keep the
# matrix and settings it was called with, then redraws them as a plotly
# heatmap in pheatmap's clustered row/column order.
.omni_trace_pheatmap <- function(cap) {
  if (!pkgs_ok("pheatmap")) return(FALSE)
  suppressMessages(trace(
    "pheatmap", where = asNamespace("pheatmap"), print = FALSE,
    tracer = bquote(assign("args", list(
      mat = mat, scale = scale, color = color, main = main),
      envir = .(cap)))))
  TRUE
}

.omni_pheatmap_to_plotly <- function(val, args) {
  m <- as.matrix(args$mat)
  if (identical(args$scale, "row")) m <- t(scale(t(m)))
  if (identical(args$scale, "column")) m <- scale(m)
  ord <- function(tree, n) if (is.list(tree) && !is.null(tree$order)) tree$order else seq_len(n)
  m <- m[ord(val$tree_row, nrow(m)), ord(val$tree_col, ncol(m)), drop = FALSE]
  scaled <- args$scale %in% c("row", "column")
  plotly::plot_ly(
    x = colnames(m), y = rownames(m), z = m, type = "heatmap",
    colors = args$color, zmid = if (scaled) 0,
    hovertemplate = paste0("%{y}<br>%{x}<br>", if (scaled) "z-score" else "value",
                           ": %{z:.3f}<extra></extra>")) |>
    plotly::layout(title = list(text = gsub("\n", "<br>", args$main %||% "")),
                   yaxis = list(autorange = "reversed", title = ""),
                   xaxis = list(title = "", tickangle = -45))
}

.omni_png_tag <- function(file) {
  shiny::tags$img(src = base64enc::dataURI(file = file, mime = "image/png"),
                  style = "width: 100%; height: 100%; object-fit: contain;")
}

# Internal column names that ggplotly copies into hover text
.omni_hover_labels <- c(".group" = "Contrast", "neglog10padj" = "-log10(padj)",
                        "log2FoldChange" = "log2FC", "padj" = "padj")

.omni_fill <- function(w) {
  w$width <- "100%"; w$height <- "100%"
  if (inherits(w, "plotly")) {
    w <- plotly::plotly_build(w)
    w$x$data <- lapply(w$x$data, function(tr) {
      if (is.character(tr$text))
        for (k in names(.omni_hover_labels))
          tr$text <- gsub(paste0("(^|<br />)", k, ":"),
                          paste0("\\1", .omni_hover_labels[[k]], ":"),
                          tr$text, perl = TRUE)
      tr
    })
    # Toolbar only on hover, so it does not sit on top of the plot title
    w <- plotly::config(w, displaylogo = FALSE, displayModeBar = "hover")
  }
  w
}

render_smart <- function(expr, env = parent.frame()) {
  force(env)  # capture the caller's frame now, not inside renderUI
  expr <- substitute(expr)
  shiny::renderUI({
    png_file <- tempfile(fileext = ".png")
    grDevices::png(png_file, width = 1100, height = 760, res = 110)
    dev <- grDevices::dev.cur()
    # Package plot helpers sometimes print() widgets; keep that from opening
    # a browser window on the machine running the app.
    old <- options(viewer = function(...) invisible(NULL),
                   browser = function(...) invisible(NULL))
    cap    <- new.env()
    traced <- pkgs_ok("plotly") && .omni_trace_pheatmap(cap)
    val <- tryCatch(
      eval(expr, env),
      finally = {
        options(old)
        if (traced) suppressMessages(untrace("pheatmap", where = asNamespace("pheatmap")))
        if (dev %in% grDevices::dev.list()) grDevices::dev.off(dev)
      })

    if (inherits(val, "htmlwidget"))
      return(.omni_fill(val))

    if (inherits(val, "pheatmap") && !is.null(cap$args)) {
      w <- tryCatch(.omni_pheatmap_to_plotly(val, cap$args), error = function(e) NULL)
      if (!is.null(w)) return(.omni_fill(w))
    }

    if (inherits(val, c("patchwork", "gglist")) && pkgs_ok("plotly")) {
      w <- tryCatch(suppressWarnings(suppressMessages(.omni_panels_to_plotly(val))),
                    error = function(e) NULL)
      if (!is.null(w)) return(.omni_fill(w))
    }

    is_gg <- inherits(val, c("gg", "ggplot")) &&
      !inherits(val, .omni_static_classes)
    if (is_gg && pkgs_ok("plotly")) {
      w <- tryCatch(suppressWarnings(suppressMessages(
        plotly::ggplotly(.omni_repel_to_text(val)))), error = function(e) NULL)
      if (!is.null(w)) return(.omni_fill(w))
    }

    # Image fallback. A ggplot that was only returned (not printed) has not
    # been drawn on the capture device yet.
    if (inherits(val, c("gg", "ggplot", "patchwork"))) {
      grDevices::png(png_file, width = 1100, height = 760, res = 110)
      tryCatch(print(val), finally = grDevices::dev.off())
    }
    shiny::validate(shiny::need(file.exists(png_file) && file.size(png_file) > 2000,
                                "Nothing to plot."))
    .omni_png_tag(png_file)
  })
}

# DE results table for a contrast as a plain data.frame with a Feature column.
de_table <- function(obj, contrast = NULL) {
  res <- omni_results(obj, contrast = contrast)
  if (is.null(res)) return(NULL)
  res <- as.data.frame(res)
  cbind(Feature = rownames(res), res, row.names = NULL)
}

# Enrichment result (clusterProfiler) -> data.frame, or NULL.
as_result_df <- function(x) {
  if (is.null(x)) return(NULL)
  df <- tryCatch(as.data.frame(x), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}

# Standard DT table with export buttons left to download handlers.
omni_dt <- function(df, page_length = 15) {
  num <- vapply(df, is.numeric, logical(1))
  dt  <- DT::datatable(df, rownames = FALSE, filter = "top",
                       options = list(pageLength = page_length,
                                      scrollX = TRUE))
  if (any(num))
    dt <- DT::formatSignif(dt, columns = names(df)[num], digits = 4)
  dt
}

csv_download <- function(filename, data_fn) {
  shiny::downloadHandler(
    filename = function() filename,
    content  = function(file) {
      df <- data_fn()
      shiny::req(df)
      utils::write.csv(df, file, row.names = FALSE)
    }
  )
}

# Read a CSV/TSV/TXT upload with the first column as row names.
read_table_upload <- function(path, name) {
  sep <- if (grepl("\\.(tsv|txt)$", name, ignore.case = TRUE)) "\t" else ","
  df  <- utils::read.delim(path, sep = sep, check.names = FALSE,
                           stringsAsFactors = FALSE)
  if (anyDuplicated(df[[1]]))
    stop("First column of '", name, "' contains duplicated IDs.", call. = FALSE)
  rownames(df) <- as.character(df[[1]])
  df[-1]
}

# Validate counts + metadata before init_omniseq(). Returns a character
# vector of problems (empty when fine) and the coerced count matrix.
check_inputs <- function(counts, meta) {
  problems <- character(0)
  non_num  <- names(counts)[!vapply(counts, is.numeric, logical(1))]
  if (length(non_num) > 0)
    problems <- c(problems, paste0("Non-numeric count columns: ",
                                   paste(utils::head(non_num, 5), collapse = ", ")))
  missing_meta <- setdiff(colnames(counts), rownames(meta))
  if (length(missing_meta) > 0)
    problems <- c(problems, paste0(length(missing_meta),
                                   " count-matrix sample(s) absent from metadata: ",
                                   paste(utils::head(missing_meta, 5), collapse = ", ")))
  if (length(problems) > 0) return(list(problems = problems, counts = NULL))

  mat <- as.matrix(counts)
  if (any(mat < 0, na.rm = TRUE))
    problems <- c(problems, "Count matrix contains negative values.")
  if (anyNA(mat))
    problems <- c(problems, "Count matrix contains missing values.")
  if (length(problems) == 0 && any(mat != round(mat)))
    problems <- c(problems, paste0(
      "Count matrix contains non-integer values. DESeq2 needs raw counts; ",
      "tick 'Round counts' if these are estimated counts (e.g. salmon/RSEM)."))
  list(problems = problems, counts = mat)
}

# Simulated demo dataset: 4 vs 4 samples, two batches, ~10% DE features plus
# an interferon response (the package's ISG signatures, up in Treated) so the
# masking, sensor-score and mediation modules have signal to work with.
# Uses real human gene symbols when org.Hs.eg.db is installed so enrichment
# modules have something to map; otherwise generic gene IDs.
demo_data <- function(n_genes = 3000, seed = 1) {
  set.seed(seed)
  isgs    <- unique(c(omniseq_isgs_28, omniseq_isgs_extended))
  n_other <- n_genes - length(isgs)
  symbols <- NULL
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    all_sym <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "SYMBOL")
    all_sym <- all_sym[grepl("^[A-Z][A-Z0-9]+$", all_sym) &
                         !grepl("^LOC[0-9]", all_sym)]
    all_sym <- setdiff(all_sym, isgs)
    if (length(all_sym) >= n_other) symbols <- sample(all_sym, n_other)
  }
  symbols <- c(isgs, symbols %||% sprintf("GENE%04d", seq_len(n_other)))

  samples <- paste0("S", 1:8)
  meta <- data.frame(
    Condition = factor(rep(c("Control", "Treated"), each = 4),
                       levels = c("Control", "Treated")),
    Batch     = factor(rep(c("B1", "B2"), times = 4)),
    row.names = samples
  )

  base_mean <- stats::rlnorm(n_genes, meanlog = 5, sdlog = 1.5)
  lfc       <- numeric(n_genes)
  de_idx    <- sample(n_genes, round(n_genes * 0.1))
  lfc[de_idx] <- stats::rnorm(length(de_idx), 0, 2)
  # Graded ISG response: varies per treated sample so ISG scores correlate
  # with other per-sample measures rather than being a clean step.
  isg_idx   <- seq_along(isgs)
  isg_dose  <- c(0, 0, 0, 0, 1, 1.5, 2, 2.5)
  batch_eff <- stats::rnorm(n_genes, 0, 0.2)

  mu <- sapply(seq_along(samples), function(j) {
    eff <- lfc * (meta$Condition[j] == "Treated")
    eff[isg_idx] <- isg_dose[j]
    base_mean * 2^(eff + batch_eff * (meta$Batch[j] == "B2"))
  })
  counts <- matrix(stats::rnbinom(length(mu), mu = mu, size = 10),
                   nrow = n_genes, dimnames = list(symbols, samples))
  list(counts = counts, meta = meta)
}

# Common UI fragment: a card holding a plot with a fixed height.
plot_card <- function(title, output_id, height = "450px", ...) {
  bslib::card(
    full_screen = TRUE,
    bslib::card_header(title),
    ...,
    plot_slot(output_id, height = height)
  )
}

# Wrap a plot_omni_*() call in a downloadable PNG.
png_download <- function(filename, plot_fn, width = 10, height = 8) {
  shiny::downloadHandler(
    filename = function() filename,
    content  = function(file) {
      grDevices::png(file, width = width, height = height, units = "in",
                     res = 300)
      on.exit(grDevices::dev.off())
      show_plot(plot_fn())
    }
  )
}

# Needs-an-object guards used inside render functions.
need_obj <- function(state) {
  shiny::validate(shiny::need(state$obj, "Load data on the Data tab first."))
  state$obj
}

need_model <- function(state) {
  obj <- need_obj(state)
  shiny::validate(shiny::need(!is.null(omni_model(obj)),
                              "Run differential expression first (DE tab)."))
  obj
}

need_contrast <- function(state) {
  obj <- need_obj(state)
  shiny::validate(shiny::need(length(omni_contrasts(obj)) > 0,
                              "Run differential expression first (DE tab)."))
  obj
}

# Sample-level metadata columns usable as groupings (few levels).
group_columns <- function(obj, max_levels = 20) {
  meta <- omni_metadata(obj)
  keep <- vapply(meta, function(x) {
    n <- length(unique(stats::na.omit(x)))
    n >= 2 && n <= max_levels
  }, logical(1))
  names(meta)[keep]
}
