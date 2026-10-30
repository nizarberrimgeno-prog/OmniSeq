# =============================================================================
# OmniSeq — 18_Static_Comparison_Plots.R
#
# Static (ggplot2, publication/presentation-ready) comparison plots developed
# for the TLR7/TLR9 pDC Q4 formal-interaction-test presentation materials and
# promoted here so future analyses get them as package functions instead of
# per-project scripts.
#
#   plot_omni_interaction_scatter_static()  Two contrasts' log2FC against each
#                                            other (one point per gene),
#                                            colored by the interaction test's
#                                            own significance AND direction.
#   plot_omni_enrichment_comparison()       Term x cluster dot-plot comparing
#                                            enrichment across multiple
#                                            compareCluster results in one
#                                            figure (e.g. the same ontology
#                                            run on several contrasts).
#   plot_omni_interaction_dumbbell()        Paired WT/KO dot-and-line plot
#                                            per gene -- READ THE WARNING in
#                                            its docs before using `rank_by =
#                                            "gap"` on a non-significant
#                                            interaction.
#
# Design choices baked in from hard-won lessons on the TLR7/TLR9 project:
#  - Scatter axes are forced symmetric around 0 (not just min/max of the
#    data) so the y=x diagonal always runs corner-to-corner and distance
#    from it reads consistently regardless of where the data happen to sit.
#  - A third color (not just significant/not) flags genes whose interaction
#    is a genuine DIRECTION REVERSAL (opposite sign in both contrasts) vs.
#    the far more common amplification/dampening pattern (same direction,
#    different magnitude) -- these are mechanistically different findings
#    and were visually indistinguishable before this was added. Reversal
#    requires |value| > a minimum magnitude in BOTH contrasts, not just
#    opposite signs -- without that, near-zero noise on either side of zero
#    registers as "reversed" by sign alone, wildly inflating the count (this
#    was caught and fixed on the real dataset: a naive sign-only rule flagged
#    142-1,690 "reversed" genes per contrast; a magnitude-gated rule found
#    the true number was 0-1).
#  - No correlation coefficient is shown on interaction scatters by default:
#    r is easy to misread as "how much do the genotypes differ" when it
#    actually measures something closer to the opposite ("how proportionally
#    coupled are the two responses", which stays high even under a large,
#    real, consistent amplification) -- removed after this confusion
#    recurred in practice.
# =============================================================================

#' Static Interaction-Term Scatter (ggplot2)
#'
#' Static twin of [plot_omni_interaction_scatter()], for presentations/
#' publications. Plots the two contrasts behind a formal interaction term
#' against each other (one point per gene), colored by whether the
#' **interaction coefficient itself** is significant, and further split by
#' whether a significant gene is an amplification/dampening (same sign in
#' both contrasts) or a genuine direction reversal (opposite sign, gated on
#' a minimum magnitude in both to exclude near-zero noise -- see file
#' header).
#'
#' @inheritParams plot_omni_interaction_scatter
#' @param flip_lfc_thresh Minimum absolute log2FC required in BOTH contrasts
#'   for a sign-flipped gene to be colored as "reversed" rather than treated
#'   as amplification noise. Default `0.3`.
#' @param max_bg_points Downsample the "Not significant" background cloud to
#'   at most this many points. `NULL` disables downsampling. Default `4000`
#'   (static plots tolerate more points than the interactive ones since
#'   there's no per-point tooltip payload).
#' @param point_size,point_alpha Aesthetics for significant points.
#'
#' @return A ggplot object (returned invisibly, also printed).
#' @export
plot_omni_interaction_scatter_static <- function(obj,
                                                 contrast_a,
                                                 contrast_b,
                                                 interaction_contrast,
                                                 label_a          = NULL,
                                                 label_b          = NULL,
                                                 pval_thresh      = 0.05,
                                                 flip_lfc_thresh  = 0.3,
                                                 title            = NULL,
                                                 max_bg_points    = 4000,
                                                 point_size       = 1.2,
                                                 point_alpha      = 0.7) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

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
  d <- merge(d, di[, c("Gene", "padj")], by = "Gene")
  names(d)[names(d) == "padj"] <- "Interaction_padj"

  sig <- d$Interaction_padj < pval_thresh
  flipped <- sig & sign(d$log2FoldChange_A) != sign(d$log2FoldChange_B) &
    abs(d$log2FoldChange_A) > flip_lfc_thresh & abs(d$log2FoldChange_B) > flip_lfc_thresh

  d$Category <- ifelse(flipped, "Significant (reversed direction)",
                 ifelse(sig, "Significant", "Not significant"))
  cols <- c("Not significant" = "grey70", "Significant" = "#C44E52",
            "Significant (reversed direction)" = "#8172B2")

  M   <- max(abs(c(d$log2FoldChange_A, d$log2FoldChange_B)), na.rm = TRUE)
  rng <- c(-M, M)

  bg <- d[d$Category == "Not significant", ]
  set.seed(1)
  if (!is.null(max_bg_points) && nrow(bg) > max_bg_points)
    bg <- bg[sample(nrow(bg), max_bg_points), ]
  fg <- d[d$Category != "Not significant", ]

  n_sig  <- sum(sig)
  n_flip <- sum(flipped)
  title  <- title %||% paste0(label_a, " vs ", label_b,
                               "  (interaction padj<", pval_thresh, ", ", n_sig,
                               " genes", if (n_flip > 0) paste0(", ", n_flip, " reversed") else "", ")")

  p <- ggplot2::ggplot() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    ggplot2::geom_point(data = bg, ggplot2::aes(x = log2FoldChange_A, y = log2FoldChange_B),
                         color = cols["Not significant"], alpha = 0.25, size = point_size * 0.7) +
    ggplot2::geom_point(data = fg, ggplot2::aes(x = log2FoldChange_A, y = log2FoldChange_B, color = Category),
                         alpha = point_alpha, size = point_size) +
    ggplot2::scale_color_manual(values = cols, name = NULL, breaks = c("Significant", "Significant (reversed direction)")) +
    ggplot2::coord_equal(xlim = rng, ylim = rng) +
    ggplot2::labs(title = title, x = paste0("log2FC, ", label_a), y = paste0("log2FC, ", label_b)) +
    ggplot2::theme_bw(base_size = 13)

  print(p)
  invisible(p)
}

#' Interaction Dot-Plot Enrichment Comparison Across Multiple Clusters/Contrasts
#'
#' Term (row) x cluster (column) dot plot comparing enrichment results across
#' several `compareCluster`/`enrichResult` objects in one figure -- e.g. the
#' same ontology run separately on several interaction terms or contrasts,
#' shown side by side so shared vs contrast-specific terms are visible at a
#' glance. Point size = significance, color = direction.
#'
#' @param results A **named list** of `compareClusterResult`/`enrichResult`
#'   objects (or data.frames with equivalent columns) -- one per
#'   cluster/contrast/interaction to compare. Names become the x-axis labels.
#' @param cluster_col Column distinguishing Up/Down (or your own two
#'   directions) WITHIN each result. Default `"Cluster"`.
#' @param term_col,padj_col As in [plot_omni_enrichment_interactive()].
#' @param direction_labels Named character vector of length 2 renaming the
#'   two `cluster_col` values for the legend, e.g.
#'   `c(Up = "KO > WT", Down = "KO < WT")` -- use this instead of the raw
#'   "Up"/"Down" labels whenever the underlying value is a difference-of-
#'   differences (interaction coefficient), where "Up"/"Down" would wrongly
#'   imply raw expression direction. Default keeps "Up"/"Down" as-is.
#' @param n_top Max terms to keep per result (by padj), unioned across all
#'   results so shared terms are directly comparable. Default `12`.
#' @param direction_colors Named character vector matching the (possibly
#'   relabeled) direction values to colors. Default a red/blue pair.
#' @param title,subtitle Plot labels.
#'
#' @return A ggplot object (returned invisibly, also printed).
#' @export
plot_omni_enrichment_comparison <- function(results,
                                            cluster_col      = "Cluster",
                                            term_col         = "Description",
                                            padj_col         = "p.adjust",
                                            direction_labels = NULL,
                                            n_top            = 12,
                                            direction_colors = NULL,
                                            title            = "Enrichment comparison",
                                            subtitle         = NULL) {

  .require_pkg("ggplot2")
  if (is.null(names(results)) || any(names(results) == ""))
    stop("`results` must be a named list (names become the x-axis labels).", call. = FALSE)

  extract_df <- function(x) {
    if (is.data.frame(x)) return(x)
    if (methods::.hasSlot(x, "compareClusterResult")) return(methods::slot(x, "compareClusterResult"))
    if (methods::.hasSlot(x, "result")) return(methods::slot(x, "result"))
    as.data.frame(x)
  }

  d <- do.call(rbind, lapply(names(results), function(nm) {
    x <- extract_df(results[[nm]])
    if (nrow(x) == 0) return(NULL)
    x$.group <- nm
    x
  }))
  if (is.null(d) || nrow(d) == 0) stop("No non-empty results to plot.", call. = FALSE)

  d[[padj_col]][is.na(d[[padj_col]])] <- 1
  d$neglog10padj <- -log10(pmax(d[[padj_col]], 1e-300))

  top_terms <- unique(unlist(lapply(names(results), function(nm) {
    sub <- d[d$.group == nm, ]
    sub[[term_col]][order(sub[[padj_col]])][seq_len(min(n_top, nrow(sub)))]
  })))
  d <- d[d[[term_col]] %in% top_terms, ]

  term_stats <- stats::aggregate(d[[padj_col]], by = list(term = d[[term_col]]), FUN = min)
  term_n     <- table(d[[term_col]])
  term_stats$n <- as.integer(term_n[term_stats$term])
  term_stats   <- term_stats[order(-term_stats$n, term_stats$x), ]
  d[[term_col]] <- factor(d[[term_col]], levels = rev(term_stats$term))
  d$.group      <- factor(d$.group, levels = names(results))

  if (!is.null(direction_labels)) {
    d$Direction <- direction_labels[d[[cluster_col]]]
    dir_levels  <- unname(direction_labels)
  } else {
    d$Direction <- d[[cluster_col]]
    dir_levels  <- unique(d$Direction)
  }
  direction_colors <- direction_colors %||%
    stats::setNames(c("#C44E52", "#4C72B0")[seq_along(dir_levels)], dir_levels)

  n_terms <- length(unique(d[[term_col]]))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .group, y = .data[[term_col]],
                                       color = Direction, size = neglog10padj)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = direction_colors, name = "Direction") +
    ggplot2::scale_size_continuous(name = "-log10(padj)", range = c(2, 10)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold"),
                   panel.grid.major.x = ggplot2::element_line(color = "grey85"))

  print(p)
  invisible(p)
}

#' Paired WT/KO Dot-and-Line Plot for Interaction Genes ("Dumbbell")
#'
#' For a formal interaction term, shows each gene as two dots (contrast A's
#' own log2FC, contrast B's own log2FC) connected by a line -- a wide gap
#' reads as a real difference, overlapping dots as no difference, without
#' requiring the reader to mentally reconstruct direction from a single
#' difference-of-differences value the way an interaction volcano does.
#'
#' @section Gene-selection warning (read before using `rank_by = "gap"`):
#' Ranking genes by the raw gap between the two contrasts (`|B - A|`) is
#' only meaningful when there IS a real signal to rank. For an interaction
#' term with few or zero significant genes, ranking by gap selects for
#' **noise by construction** -- some genes will show a large apparent gap
#' purely by chance even under the null, and a plot built that way visually
#' suggests a pattern the statistics explicitly reject. This was caught on
#' real data: a "top 20 by gap" plot for a 0-significant-gene interaction
#' looked like it showed real divergence; none of it held up. Use
#' `rank_by = "padj"` (the default) whenever the interaction has real
#' significant genes; for a null/near-null interaction, don't use this
#' function at all -- use [plot_omni_interaction_scatter_static()] on the
#' full gene set instead, which can't manufacture a false pattern because
#' it doesn't select a subset.
#'
#' @param obj A validated OmniSeq object.
#' @param contrast_a,contrast_b The two contrasts to compare (e.g. the KO's
#'   own simple-effect contrast and WT's).
#' @param interaction_contrast Name of the cached interaction-term contrast,
#'   used only for ranking/selecting genes (by `padj`), not for the plotted
#'   values themselves.
#' @param label_a,label_b Legend labels (e.g. `"WT"`, `"TLR7-KO"`).
#' @param stim_label Used in the x-axis title ("log2 fold-change induced by
#'   `stim_label`").
#' @param rank_by `"padj"` (default, recommended -- rank by the interaction
#'   term's own significance) or `"gap"` (see warning above).
#' @param n_top Number of genes to show. Default `20`.
#'
#' @return A ggplot object (returned invisibly, also printed).
#' @export
plot_omni_interaction_dumbbell <- function(obj,
                                           contrast_a,
                                           contrast_b,
                                           interaction_contrast,
                                           label_a     = NULL,
                                           label_b     = NULL,
                                           stim_label  = "",
                                           rank_by     = c("padj", "gap"),
                                           n_top       = 20) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  rank_by <- match.arg(rank_by)

  label_a <- label_a %||% contrast_a
  label_b <- label_b %||% contrast_b

  ra <- omni_results(obj, contrast = contrast_a)
  rb <- omni_results(obj, contrast = contrast_b)
  ri <- omni_results(obj, contrast = interaction_contrast)
  if (is.null(ra) || is.null(rb) || is.null(ri))
    stop("Missing DE results for one of the three contrasts.", call. = FALSE)

  da <- as.data.frame(ra); da$Gene <- rownames(da)
  db <- as.data.frame(rb); db$Gene <- rownames(db)
  di <- as.data.frame(ri); di$Gene <- rownames(di); di$padj[is.na(di$padj)] <- 1

  d <- merge(da[, c("Gene", "log2FoldChange")], db[, c("Gene", "log2FoldChange")],
             by = "Gene", suffixes = c("_A", "_B"))
  d <- merge(d, di[, c("Gene", "padj")], by = "Gene")

  if (rank_by == "padj") {
    d <- d[order(d$padj), ]
  } else {
    message("[plot_omni_interaction_dumbbell] rank_by='gap' selected -- see ",
            "?plot_omni_interaction_dumbbell for why this can be misleading ",
            "when the interaction has few/no significant genes.")
    d$gap <- abs(d$log2FoldChange_B - d$log2FoldChange_A)
    d <- d[order(-d$gap), ]
  }
  d <- head(d, n_top)
  d <- d[order(d$log2FoldChange_A), ]
  d$Gene <- factor(d$Gene, levels = d$Gene)

  long <- data.frame(Gene = rep(d$Gene, 2),
                      LFC  = c(d$log2FoldChange_A, d$log2FoldChange_B),
                      Group = factor(rep(c(label_a, label_b), each = nrow(d)), levels = c(label_a, label_b)))
  seg <- data.frame(Gene = d$Gene, x = d$log2FoldChange_A, xend = d$log2FoldChange_B)

  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = x, xend = xend, y = Gene, yend = Gene),
                          color = "grey60", linewidth = 0.6) +
    ggplot2::geom_point(data = long, ggplot2::aes(x = LFC, y = Gene, color = Group), size = 3) +
    ggplot2::scale_color_manual(values = stats::setNames(c("#4C72B0", "#C44E52"), c(label_a, label_b))) +
    ggplot2::labs(x = paste0("log2 fold-change induced by ", stim_label, " (vs. unstimulated)"),
                  y = NULL, color = NULL) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(legend.position = "top")

  print(p)
  invisible(p)
}
