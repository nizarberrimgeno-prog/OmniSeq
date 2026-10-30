# =============================================================================
# OmniSeq — 17_Interactive_Reporting.R
#
# Interactive (plotly/DT) counterparts to the static plot_omni_*() functions,
# built for embedding directly in an rmarkdown+flatly report chunk. Developed
# for the TLR7/TLR9 pDC "Results Organized by Research Question" report and
# promoted here so future analyses can call them on any OmniSeq object
# instead of hand-rolling per-project versions.
#
#   plot_omni_volcano_interactive()   Volcano w/ hover (gene/LFC/padj), same
#                                      padj+LFC significance rule and color
#                                      scheme as plot_omni_volcano().
#   plot_omni_lfc_scatter_interactive() Two contrasts' log2FC against each
#                                      other, colored Shared/A-specific/
#                                      B-specific/Not significant.
#   plot_omni_interaction_scatter()   The two contrasts behind a genotype x
#                                      stimulus (or any) interaction term,
#                                      colored by the INTERACTION test's own
#                                      significance (not either contrast's).
#   omni_gene_search_table()          Full DE table as a searchable/sortable
#                                      DT::datatable().
#   plot_omni_enrichment_interactive() Interactive dot plot + linked gene
#                                      table for a clusterProfiler enrichment
#                                      result (enrichResult / compareCluster).
#
# Threshold convention (matches plot_omni_volcano() / the fixed
# plot_omni_heatmap()/plot_omni_top_genes()): significant = padj < pval_thresh
# AND |log2FoldChange| > lfc_thresh. Every plot states its exact contrast name
# and thresholds in the title so it never has to be inferred.
#
# Performance: background "Not significant" points are downsampled to
# max_bg_points (significant points are always kept in full) to keep the
# rendered HTML from becoming enormous/slow when many such plots sit in one
# document (esp. inside {.tabset} panes, which all render into the DOM at
# load regardless of which tab is visible). Plain SVG scatter (not scattergl)
# is used deliberately -- a report with a dozen+ such plots can exceed the
# browser's simultaneous WebGL-context limit with scattergl.
#
# IMPORTANT for plot_omni_enrichment_interactive(): return its result as the
# chunk's sole, un-printed expression (e.g. the last line of the chunk, not
# wrapped in print()). Printing plotly/DT widgets via separate print() calls
# in one chunk silently breaks htmlwidgets' dependency injection in knitted
# HTML -- verified in isolation. Returning one combined htmltools::tagList()
# as the single auto-printed value is the only combination confirmed to work.
# =============================================================================

#' @keywords internal
.omni_downsample_bg <- function(df, max_bg_points) {
  if (!is.null(max_bg_points) && nrow(df) > max_bg_points)
    df <- df[sample(nrow(df), max_bg_points), ]
  df
}

#' Add one scatter trace per category level, each with its own marker colour
#' @keywords internal
.omni_add_category_traces <- function(p, d, col, levels, cols, x, y,
                                      hover_txt, size, opacity) {
  for (lv in levels) {
    dd <- d[as.character(d[[col]]) == lv, , drop = FALSE]
    if (nrow(dd) == 0) next
    p <- plotly::add_trace(p, data = dd, x = x, y = y,
                           type = "scatter", mode = "markers",
                           marker = list(color = unname(cols[lv]), size = size,
                                         opacity = opacity),
                           text = hover_txt(dd), hoverinfo = "text",
                           name = paste0(lv, " (n=", nrow(dd), ")"))
  }
  p
}

#' @keywords internal
.omni_require_interactive_pkgs <- function() {
  .require_pkg("plotly")
  .require_pkg("htmltools")
}

# =============================================================================
# 1. INTERACTIVE VOLCANO
# =============================================================================

#' Interactive Volcano Plot (plotly)
#'
#' Same significance rule and color scheme as [plot_omni_volcano()]
#' (`padj < pval_thresh & |log2FoldChange| > lfc_thresh`), rendered as an
#' interactive plotly scatter with per-gene hover text (gene / log2FC / padj)
#' instead of a static ggplot. Intended for direct use inside an
#' rmarkdown HTML chunk.
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param contrast Contrast name to pull from `results_cache`. Defaults to
#'   `obj$active_contrast`.
#' @param pval_thresh Adjusted p-value threshold. Default `0.05`.
#' @param lfc_thresh Absolute log2FoldChange threshold. Default `1`.
#' @param max_bg_points Downsample the "Not significant" background cloud to
#'   at most this many points (significant points are always kept in full).
#'   `NULL` disables downsampling. Default `1500`.
#' @param title Optional custom title. Defaults to a title stating the
#'   contrast and both thresholds explicitly.
#' @param width,height Plot pixel dimensions (fixed, not responsive — needed
#'   so sizing doesn't depend on a hidden tab pane's computed dimensions).
#'
#' @return A plotly object.
#' @export
plot_omni_volcano_interactive <- function(obj,
                                          contrast      = NULL,
                                          pval_thresh   = 0.05,
                                          lfc_thresh    = 1,
                                          max_bg_points = 1500,
                                          title         = NULL,
                                          width         = 750,
                                          height         = 550) {

  validate_omniseq(obj)
  .omni_require_interactive_pkgs()

  key <- contrast %||% obj$active_contrast
  res <- omni_results(obj, contrast = key)
  if (is.null(res))
    stop("No DE results for contrast '", key, "'. Run run_omni_de()/",
         "run_omni_contrast() first.", call. = FALSE)

  df      <- as.data.frame(res)
  df$Gene <- rownames(df)
  df      <- df[!is.na(df$log2FoldChange), ]
  df$padj <- ifelse(is.na(df$padj), 1, df$padj)

  df$Direction <- dplyr::case_when(
    df$padj < pval_thresh & df$log2FoldChange >  lfc_thresh ~ "Up",
    df$padj < pval_thresh & df$log2FoldChange < -lfc_thresh ~ "Down",
    TRUE ~ "Not significant"
  )
  n_up   <- sum(df$Direction == "Up")
  n_down <- sum(df$Direction == "Down")

  df$neglog10padj <- -log10(pmax(df$padj, 1e-300))
  cols <- c("Up" = "#E74C3C", "Down" = "#2980B9", "Not significant" = "grey75")

  d_bg <- .omni_downsample_bg(df[df$Direction == "Not significant", ], max_bg_points)
  d_fg <- df[df$Direction != "Not significant", ]

  hover_txt <- function(dd) paste0("Gene: ", dd$Gene,
                                    "<br>log2FC: ", round(dd$log2FoldChange, 2),
                                    "<br>padj: ", signif(dd$padj, 3))

  title <- title %||% paste0(key, "  (padj<", pval_thresh, " & |log2FC|>", lfc_thresh,
                              "; ", n_up, " Up / ", n_down, " Down)")

  p <- plotly::plot_ly(width = width, height = height) |>
    plotly::add_trace(data = d_bg, x = ~log2FoldChange, y = ~neglog10padj,
                       type = "scatter", mode = "markers",
                       marker = list(color = cols["Not significant"], size = 4, opacity = 0.35),
                       text = hover_txt(d_bg), hoverinfo = "text", name = "Not significant")
  # One trace per direction with an explicit marker colour: `colors =` is a
  # plot_ly() argument and is ignored by add_trace(), which left plotly's
  # default palette in place of the Up/Down colours.
  p <- .omni_add_category_traces(p, d_fg, "Direction", c("Up", "Down"), cols,
                                 x = ~log2FoldChange, y = ~neglog10padj,
                                 hover_txt = hover_txt, size = 6, opacity = 0.75)
  p |>
    plotly::layout(
      title  = title,
      xaxis  = list(title = "log2 Fold Change"),
      yaxis  = list(title = "-log10 adjusted p-value"),
      shapes = list(
        list(type = "line", x0 = -lfc_thresh, x1 = -lfc_thresh, y0 = 0, y1 = 1,
             yref = "paper", line = list(dash = "dot", color = "black", width = 1)),
        list(type = "line", x0 = lfc_thresh, x1 = lfc_thresh, y0 = 0, y1 = 1,
             yref = "paper", line = list(dash = "dot", color = "black", width = 1)),
        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
             y0 = -log10(pval_thresh), y1 = -log10(pval_thresh),
             line = list(dash = "dot", color = "black", width = 1))
      ),
      legend = list(title = list(text = ""))
    )
}

# =============================================================================
# 2. INTERACTIVE LFC-vs-LFC SCATTER (two contrasts)
# =============================================================================

#' Interactive Fold-Change Concordance Scatter (plotly)
#'
#' Plots one contrast's log2FoldChange against another's, one point per gene,
#' colored by whether each gene is significant in contrast A only, contrast B
#' only, both ("Shared"), or neither. Useful for shared-vs-specific
#' transcriptional program comparisons (e.g. TLR7 vs TLR9 stimulation).
#'
#' @param obj A validated OmniSeq object with DE results for both contrasts.
#' @param contrast_a,contrast_b Contrast names to pull from `results_cache`.
#' @param label_a,label_b Axis/legend labels. Default to the contrast names.
#' @param pval_thresh Adjusted p-value threshold. Default `0.05`.
#' @param lfc_thresh Absolute log2FoldChange threshold. Default `1`.
#' @param max_bg_points Downsample the "Not significant" points. Default `1500`.
#' @param title Optional custom title.
#' @param width,height Plot pixel dimensions.
#'
#' @return A plotly object.
#' @export
plot_omni_lfc_scatter_interactive <- function(obj,
                                              contrast_a,
                                              contrast_b,
                                              label_a       = NULL,
                                              label_b       = NULL,
                                              pval_thresh   = 0.05,
                                              lfc_thresh    = 1,
                                              max_bg_points = 1500,
                                              title         = NULL,
                                              width         = 800,
                                              height        = 550) {

  validate_omniseq(obj)
  .omni_require_interactive_pkgs()

  label_a <- label_a %||% contrast_a
  label_b <- label_b %||% contrast_b

  ra <- omni_results(obj, contrast = contrast_a)
  rb <- omni_results(obj, contrast = contrast_b)
  if (is.null(ra)) stop("No DE results for contrast '", contrast_a, "'.", call. = FALSE)
  if (is.null(rb)) stop("No DE results for contrast '", contrast_b, "'.", call. = FALSE)

  da <- as.data.frame(ra); da$Gene <- rownames(da); da$padj[is.na(da$padj)] <- 1
  db <- as.data.frame(rb); db$Gene <- rownames(db); db$padj[is.na(db$padj)] <- 1
  d  <- merge(da[, c("Gene", "log2FoldChange", "padj")],
              db[, c("Gene", "log2FoldChange", "padj")],
              by = "Gene", suffixes = c("_A", "_B"))

  sig_a <- d$padj_A < pval_thresh & abs(d$log2FoldChange_A) > lfc_thresh
  sig_b <- d$padj_B < pval_thresh & abs(d$log2FoldChange_B) > lfc_thresh
  d$Category <- ifelse(sig_a & sig_b, "Shared",
                 ifelse(sig_a, paste0(label_a, "-specific"),
                 ifelse(sig_b, paste0(label_b, "-specific"), "Not significant")))
  cat_levels <- c("Shared", paste0(label_a, "-specific"), paste0(label_b, "-specific"), "Not significant")
  d$Category <- factor(d$Category, levels = cat_levels)
  cols <- stats::setNames(c("#8172B2", "#4C72B0", "#DD8452", "grey80"), cat_levels)

  d_bg <- .omni_downsample_bg(d[d$Category == "Not significant", ], max_bg_points)
  d_fg <- d[d$Category != "Not significant", ]

  hover_txt <- function(dd) paste0("Gene: ", dd$Gene,
                                    "<br>", label_a, " log2FC: ", round(dd$log2FoldChange_A, 2),
                                    "<br>", label_b, " log2FC: ", round(dd$log2FoldChange_B, 2))

  title <- title %||% paste0(label_a, " vs ", label_b,
                              "  (padj<", pval_thresh, " & |log2FC|>", lfc_thresh, ")")

  p <- plotly::plot_ly(width = width, height = height) |>
    plotly::add_trace(data = d_bg, x = ~log2FoldChange_A, y = ~log2FoldChange_B,
                       type = "scatter", mode = "markers",
                       marker = list(color = cols["Not significant"], size = 4, opacity = 0.25),
                       text = hover_txt(d_bg), hoverinfo = "text", name = "Not significant")
  p <- .omni_add_category_traces(p, d_fg, "Category", setdiff(cat_levels, "Not significant"),
                                 cols, x = ~log2FoldChange_A, y = ~log2FoldChange_B,
                                 hover_txt = hover_txt, size = 6, opacity = 0.75)
  p |>
    plotly::layout(title = title,
                    xaxis = list(title = paste0("log2FC, ", label_a)),
                    yaxis = list(title = paste0("log2FC, ", label_b)),
                    legend = list(title = list(text = "Category")))
}

# =============================================================================
# 3. INTERACTION SCATTER (formal genotype/condition x stimulus test)
# =============================================================================

#' Interactive Interaction-Term Scatter (plotly)
#'
#' For a formal interaction test (e.g. from a `~ ... + FactorA * FactorB`
#' model), plots the two simple-effect contrasts that the interaction
#' coefficient contrasts against each other — one point per gene — colored by
#' whether the **interaction coefficient itself** is significant (padj
#' alone, no LFC filter — the established convention for interaction
#' coefficients), and, among significant genes, further split into
#' amplification/dampening (same sign in both contrasts) vs. a genuine
#' **direction reversal** (opposite sign, gated on a minimum magnitude in
#' both contrasts so near-zero noise isn't miscounted as "reversed" — see
#' `flip_lfc_thresh`). A gene on the dotted y=x line has an identical effect
#' in both arms (interaction = 0); its distance from that line is the
#' interaction log2FC. Axes are forced symmetric around 0 so the y=x
#' diagonal always runs corner-to-corner. Interactive twin of
#' [plot_omni_interaction_scatter_static()]; keep both in sync.
#'
#' @param obj A validated OmniSeq object.
#' @param contrast_a,contrast_b The two simple-effect contrast names (e.g.
#'   `"KO_StimX vs KO_NS"` and `"WT_StimX vs WT_NS"`) that the interaction
#'   term is the difference of.
#' @param interaction_contrast Name of the cached interaction-term contrast
#'   (its DE results should carry the Wald test for the FactorA:FactorB
#'   coefficient from the joint interaction model).
#' @param label_a,label_b Axis labels. Default to the contrast names.
#' @param pval_thresh Adjusted p-value threshold for the interaction test.
#'   Default `0.05`.
#' @param flip_lfc_thresh Minimum absolute log2FC required in BOTH contrasts
#'   for a sign-flipped significant gene to be colored as "reversed" rather
#'   than amplification noise. Default `0.3`.
#' @param title Optional custom title.
#' @param width,height Plot pixel dimensions.
#'
#' @return A plotly object.
#' @export
plot_omni_interaction_scatter <- function(obj,
                                          contrast_a,
                                          contrast_b,
                                          interaction_contrast,
                                          label_a         = NULL,
                                          label_b         = NULL,
                                          pval_thresh     = 0.05,
                                          flip_lfc_thresh = 0.3,
                                          title           = NULL,
                                          width           = 800,
                                          height          = 550) {

  validate_omniseq(obj)
  .omni_require_interactive_pkgs()

  label_a <- label_a %||% contrast_a
  label_b <- label_b %||% contrast_b

  ra <- omni_results(obj, contrast = contrast_a)
  rb <- omni_results(obj, contrast = contrast_b)
  ri <- omni_results(obj, contrast = interaction_contrast)
  if (is.null(ra)) stop("No DE results for contrast '", contrast_a, "'.", call. = FALSE)
  if (is.null(rb)) stop("No DE results for contrast '", contrast_b, "'.", call. = FALSE)
  if (is.null(ri)) stop("No DE results for interaction contrast '",
                        interaction_contrast, "'.", call. = FALSE)

  da <- as.data.frame(ra); da$Gene <- rownames(da)
  db <- as.data.frame(rb); db$Gene <- rownames(db)
  di <- as.data.frame(ri); di$Gene <- rownames(di); di$padj[is.na(di$padj)] <- 1

  d <- merge(da[, c("Gene", "log2FoldChange")], db[, c("Gene", "log2FoldChange")],
             by = "Gene", suffixes = c("_A", "_B"))
  d <- merge(d, di[, c("Gene", "log2FoldChange", "padj")], by = "Gene")
  names(d)[names(d) == "log2FoldChange"] <- "Interaction_LFC"
  names(d)[names(d) == "padj"]           <- "Interaction_padj"

  # Direction-reversal flag mirrors plot_omni_interaction_scatter_static():
  # sign(A) != sign(B) alone massively overcounts "reversed" genes from
  # near-zero noise -- gate on a minimum magnitude in BOTH contrasts too
  # (naive sign-only rule flagged 142-1,690 genes/contrast on real data;
  # the magnitude-gated rule found the true number was 0-1).
  sig <- d$Interaction_padj < pval_thresh
  flipped <- sig & sign(d$log2FoldChange_A) != sign(d$log2FoldChange_B) &
    abs(d$log2FoldChange_A) > flip_lfc_thresh & abs(d$log2FoldChange_B) > flip_lfc_thresh

  d$Interaction_sig <- ifelse(flipped, "Significant (reversed direction)",
                        ifelse(sig, "Significant", "Not significant"))
  cols <- stats::setNames(c("grey80", "#C44E52", "#8172B2"),
                           c("Not significant", "Significant", "Significant (reversed direction)"))

  d_bg <- d[d$Interaction_sig == "Not significant", ]
  d_fg <- d[d$Interaction_sig != "Not significant", ]

  hover_txt <- function(dd) paste0("Gene: ", dd$Gene,
                                    "<br>", label_a, " log2FC: ", round(dd$log2FoldChange_A, 2),
                                    "<br>", label_b, " log2FC: ", round(dd$log2FoldChange_B, 2),
                                    "<br>Interaction log2FC: ", round(dd$Interaction_LFC, 2),
                                    "<br>Interaction padj: ", signif(dd$Interaction_padj, 3))

  # Axes forced symmetric around 0 (not just min/max) so the y=x diagonal
  # always runs corner-to-corner and distance from it reads consistently.
  M   <- max(abs(c(d$log2FoldChange_A, d$log2FoldChange_B)), na.rm = TRUE)
  rng <- c(-M, M)

  n_sig  <- sum(sig)
  n_flip <- sum(flipped)
  title  <- title %||% paste0(label_a, " vs ", label_b,
                               "  (interaction padj<", pval_thresh, ", ", n_sig,
                               " genes", if (n_flip > 0) paste0(", ", n_flip, " reversed") else "", ")")

  p <- plotly::plot_ly(width = width, height = height) |>
    plotly::add_trace(x = rng, y = rng, type = "scatter", mode = "lines",
                       line = list(color = "black", dash = "dot", width = 1),
                       hoverinfo = "skip", showlegend = FALSE) |>
    plotly::add_trace(data = d_bg, x = ~log2FoldChange_A, y = ~log2FoldChange_B,
                       type = "scatter", mode = "markers",
                       marker = list(color = cols["Not significant"], size = 4, opacity = 0.25),
                       text = hover_txt(d_bg), hoverinfo = "text", name = "Not significant")

  for (cat in c("Significant", "Significant (reversed direction)")) {
    dd <- d_fg[d_fg$Interaction_sig == cat, ]
    if (nrow(dd) == 0) next
    p <- p |>
      plotly::add_trace(data = dd, x = ~log2FoldChange_A, y = ~log2FoldChange_B,
                         type = "scatter", mode = "markers",
                         marker = list(color = cols[cat], size = 6, opacity = 0.8),
                         text = hover_txt(dd), hoverinfo = "text",
                         name = paste0(cat, " (n=", nrow(dd), ")"))
  }

  p |>
    plotly::layout(title = title,
                    xaxis = list(title = paste0("log2FC, ", label_a), range = rng, scaleanchor = "y"),
                    yaxis = list(title = paste0("log2FC, ", label_b), range = rng),
                    legend = list(title = list(text = "Formal interaction test")))
}

# =============================================================================
# 4. SEARCHABLE GENE TABLE
# =============================================================================

#' Searchable/Sortable DE Results Table (DT)
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param contrast Contrast name. Defaults to `obj$active_contrast`.
#' @param pval_thresh If given, pre-filter to `padj < pval_thresh` before
#'   display. `NULL` (default) shows all genes.
#' @param digits Rounding for log2FoldChange; p-values/padj use `signif`.
#' @param page_length Rows per page. Default `10`.
#' @param caption Optional table caption.
#'
#' @return A `DT::datatable` htmlwidget.
#' @export
omni_gene_search_table <- function(obj,
                                   contrast    = NULL,
                                   pval_thresh = NULL,
                                   digits      = 3,
                                   page_length = 10,
                                   caption     = NULL) {

  validate_omniseq(obj)
  .require_pkg("DT")

  key <- contrast %||% obj$active_contrast
  res <- omni_results(obj, contrast = key)
  if (is.null(res))
    stop("No DE results for contrast '", key, "'.", call. = FALSE)

  d <- as.data.frame(res)
  d$Gene <- rownames(d)
  d <- d[, c("Gene", setdiff(names(d), "Gene"))]
  d$padj[is.na(d$padj)] <- 1
  if (!is.null(pval_thresh)) d <- d[d$padj < pval_thresh, ]

  if ("log2FoldChange" %in% names(d)) d$log2FoldChange <- round(d$log2FoldChange, digits)
  if ("pvalue" %in% names(d)) d$pvalue <- signif(d$pvalue, 3)
  d$padj <- signif(d$padj, 3)

  caption <- caption %||% key
  padj_col <- which(names(d) == "padj") - 1
  DT::datatable(d, filter = "top", rownames = FALSE, caption = caption,
                options = list(pageLength = page_length, scrollX = TRUE,
                                order = list(list(padj_col, "asc"))))
}

# =============================================================================
# 5. INTERACTIVE ENRICHMENT DOT PLOT + LINKED GENE TABLE
# =============================================================================

#' Interactive Enrichment Dot Plot with Linked Gene Table (plotly + DT)
#'
#' Takes a clusterProfiler enrichment result (`enrichResult`,
#' `compareClusterResult`, or a plain data.frame with equivalent columns) and
#' renders an interactive dot plot (hover shows term/cluster/padj/gene count)
#' plus a searchable table of the underlying genes per term, combined into a
#' single `htmltools::tagList()`.
#'
#' @param enrich_result An `enrichResult`/`compareClusterResult` object, or a
#'   data.frame with columns matching `cluster_col`/`term_col`/`padj_col`/
#'   `gene_col` (default clusterProfiler column names).
#' @param cluster_col Column naming the x-axis groups (e.g. "Cluster",
#'   contrast/direction labels). If the data has only one implicit group,
#'   pass a column that is constant, or set to `NULL` to use a single
#'   "All" group.
#' @param term_col Column with the pathway/term description. Default
#'   `"Description"`.
#' @param padj_col Column with adjusted p-values. Default `"p.adjust"`.
#' @param gene_col Column with `/`-separated gene symbols per term. Default
#'   `"geneID"` (requires `readable = TRUE` in the original
#'   `enrichGO`/`compareCluster` call).
#' @param count_col Column with per-term gene counts, used for point size.
#'   Default `"Count"`.
#' @param n_top Max terms to keep, ranked by padj (kept separately per
#'   cluster). Default `25`.
#' @param color_low,color_high Two-color gradient for -log10(padj), light to
#'   dark (avoid grey — reserved for "not significant" elsewhere in this
#'   package's plots). Defaults to a blue gradient; pass a red gradient
#'   (e.g. `"#F4B8AE"`/`"#B01E2B"`) when plotting an "Up-regulated" subset
#'   alongside this function's default for "Down-regulated", to keep
#'   direction visually distinct across the two calls.
#' @param title Optional plot title.
#' @param width,height Plot pixel dimensions.
#'
#' @return An `htmltools::tagList()` combining the plotly dot plot and a
#'   `DT::datatable` of terms x genes. **Return this as the chunk's sole
#'   expression** (don't wrap in `print()`, don't call `print()` on the plot
#'   and table separately) — see file header for why.
#' @export
plot_omni_enrichment_interactive <- function(enrich_result,
                                             cluster_col = "Cluster",
                                             term_col    = "Description",
                                             padj_col    = "p.adjust",
                                             gene_col    = "geneID",
                                             count_col   = "Count",
                                             n_top       = 25,
                                             color_low   = "#AEC6E8",
                                             color_high  = "#1F4E8C",
                                             title       = NULL,
                                             width       = 800,
                                             height      = 550) {

  .require_pkg("plotly"); .require_pkg("DT"); .require_pkg("htmltools")

  ## clusterProfiler's as.data.frame() methods for enrichResult/
  ## compareClusterResult are S4 methods registered by DOSE/clusterProfiler;
  ## base as.data.frame() only finds them if that package is *attached*
  ## (library(), not just requireNamespace()), which we don't want to
  ## require of the caller. Pull the underlying data.frame slot directly
  ## instead -- works regardless of what's attached.
  d <- if (is.data.frame(enrich_result)) {
    enrich_result
  } else if (methods::.hasSlot(enrich_result, "compareClusterResult")) {
    methods::slot(enrich_result, "compareClusterResult")
  } else if (methods::.hasSlot(enrich_result, "result")) {
    methods::slot(enrich_result, "result")
  } else {
    .require_pkg("clusterProfiler")
    as.data.frame(enrich_result)
  }
  if (nrow(d) == 0) stop("Enrichment result is empty.", call. = FALSE)

  if (is.null(cluster_col) || !cluster_col %in% names(d)) {
    d$.cluster <- "All"; cluster_col <- ".cluster"
  }
  d[[padj_col]][is.na(d[[padj_col]])] <- 1

  ## Keep, per cluster, the n_top terms with the smallest padj (ties broken by
  ## first occurrence) -- base R instead of a dplyr NSE chain to avoid any
  ## .data-pronoun import fragility for a function meant to be called broadly.
  min_padj_per_term <- stats::aggregate(d[[padj_col]],
                                        by = list(cluster = d[[cluster_col]], term = d[[term_col]]),
                                        FUN = min)
  keep_terms <- unlist(lapply(split(min_padj_per_term, min_padj_per_term$cluster), function(sub) {
    sub$term[order(sub$x)][seq_len(min(n_top, nrow(sub)))]
  }), use.names = FALSE)
  d <- d[d[[term_col]] %in% keep_terms, ]

  d$neglog10padj <- -log10(pmax(d[[padj_col]], 1e-300))
  d$row_id       <- seq_len(nrow(d))

  term_min_padj <- stats::aggregate(d[[padj_col]], by = list(term = d[[term_col]]), FUN = min)
  ord <- term_min_padj$term[order(-term_min_padj$x)]

  size_vals <- if (count_col %in% names(d)) d[[count_col]] else rep(8, nrow(d))
  hover_txt <- paste0("Term: ", d[[term_col]],
                       "<br>Cluster: ", d[[cluster_col]],
                       "<br>padj: ", signif(d[[padj_col]], 3),
                       if (count_col %in% names(d)) paste0("<br>Gene count: ", d[[count_col]]) else "")

  title <- title %||% "Enrichment"

  p <- plotly::plot_ly(d, x = d[[cluster_col]], y = d[[term_col]],
                        type = "scatter", mode = "markers", key = ~row_id,
                        width = width, height = height,
                        marker = list(size = ~ (8 + 2 * sqrt(size_vals)),
                                      color = ~neglog10padj,
                                      colorscale = list(list(0, color_low), list(1, color_high)),
                                      showscale = TRUE,
                                      colorbar = list(title = "-log10(padj)")),
                        text = hover_txt, hoverinfo = "text") |>
    plotly::layout(title = title,
                    xaxis = list(title = ""),
                    yaxis = list(title = "", categoryorder = "array", categoryarray = ord),
                    margin = list(l = 250))

  gene_tbl <- if (gene_col %in% names(d)) {
    gt <- d[, c(term_col, cluster_col, padj_col, gene_col)]
    names(gt) <- c("Term", "Cluster", "padj", "Genes")
    gt$padj <- signif(gt$padj, 3)
    gt$Genes <- gsub("/", ", ", gt$Genes)
    DT::datatable(gt, filter = "top", rownames = FALSE,
                  options = list(pageLength = 8, scrollX = TRUE,
                                  order = list(list(2, "asc"))))
  } else NULL

  if (is.null(gene_tbl)) htmltools::tagList(p) else htmltools::tagList(p, gene_tbl)
}
