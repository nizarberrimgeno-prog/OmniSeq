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
    shiny::plotOutput(output_id, height = height)
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
