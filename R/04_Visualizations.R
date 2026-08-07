# =============================================================================
# OmniSeq — 04_Visualizations.R
#
# Core plotting engine — 8 functions:
#   plot_omni_pca()          PCA with covariate correction + ellipses + auto Sex shapes
#   plot_omni_volcano()      Volcano plot with labeled top genes
#   plot_omni_ma()           MA plot
#   plot_omni_heatmap()      Annotated heatmap, auto-subsets to active contrast
#   plot_omni_qc_basics()    Library size + detected genes barplots
#   plot_omni_qc_advanced()  Distance heatmap, mean-variance, before/after PCA
#   plot_omni_pvalue_hist()  P-value histogram (model diagnostic)
#   plot_omni_top_genes()    Normalized count dotplots for top N DEGs
#
# Fixes applied vs previous version:
#   - apeglm shrinkage drops zero-count genes, reducing nrow(res) vs nrow(counts)
#     All functions now build df from as.data.frame(res) and replace NA padj
#     with 1 rather than filtering them out, preventing row mismatch errors
#   - plot_omni_pca() auto-detects Sex column and maps it to point shape
#   - plot_omni_top_genes() now infers the correct grouping column from the
#     active contrast (e.g. Genotype for "TLR7 vs UNC93B1") instead of always
#     defaulting to the primary condition column
#   - .omni_contrast_column() helper added to derive the relevant metadata
#     column from any active contrast name
# =============================================================================


# =============================================================================
# 1. PCA
# =============================================================================

#' Plot PCA with Optional Covariate Correction
#'
#' Computes PCA on the top variable features from the cached VST matrix.
#' If covariates are present in the model, they are regressed out using
#' `limma::removeBatchEffect()` before PCA so the plot reflects the primary
#' biological condition. Point shape is automatically mapped to Sex if a
#' Sex column is detected in the metadata.
#'
#' @param obj A validated OmniSeq object with a fitted model.
#' @param ntop Integer. Number of top variable features to use. Default: 500.
#' @param label Logical. Label sample points with their names. Default: `TRUE`.
#' @param ellipse Logical. Draw 95% confidence ellipses per group. Suppressed
#'   automatically if any group has fewer than 3 samples. Default: `TRUE`.
#' @param color_by Character. Metadata column to color points by. Defaults to
#'   the primary condition.
#' @param shape_by Character. Metadata column to map to point shape. Auto-
#'   detected from Sex/sex/gender columns if `NULL`. Set to `NA` to disable.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_pca <- function(obj,
                          ntop     = 500,
                          label    = TRUE,
                          ellipse  = TRUE,
                          color_by = NULL,
                          shape_by = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .require_pkg("ggrepel")

  vsd <- omni_vst(obj)
  if (is.null(vsd))
    stop("VST not cached. Run run_omni_qc() or run_omni_de() first.",
         call. = FALSE)

  mat      <- SummarizedExperiment::assay(vsd)
  meta     <- omni_metadata(obj)
  cond_col <- obj$params$condition
  covars   <- obj$params$covariates
  color_by <- color_by %||% cond_col

  # Auto-detect Sex column for shape if not specified
  if (is.null(shape_by)) {
    sex_col  <- intersect(c("Sex", "sex", "SEX", "gender", "Gender"),
                          colnames(meta))
    shape_by <- if (length(sex_col) > 0) sex_col[1] else NULL
  } else if (identical(shape_by, NA)) {
    shape_by <- NULL
  }

  # Covariate correction via limma (only affects this plot, not stored data)
  if (length(covars) > 0) {
    .require_pkg("limma")
    design_bio <- stats::model.matrix(
      stats::as.formula(paste("~", cond_col)), data = meta
    )
    b1   <- if (length(covars) >= 1) meta[[covars[1]]] else NULL
    b2   <- if (length(covars) >= 2) meta[[covars[2]]] else NULL
    covs <- if (length(covars) >= 3)
      as.matrix(meta[, covars[3:length(covars)], drop = FALSE]) else NULL

    mat <- limma::removeBatchEffect(mat,
                                    batch      = b1,
                                    batch2     = b2,
                                    covariates = covs,
                                    design     = design_bio)
  }

  # Select top variable features
  rv     <- apply(mat, 1, stats::var)
  select <- order(rv, decreasing = TRUE)[1:min(ntop, nrow(mat))]
  pca    <- stats::prcomp(t(mat[select, ]))
  pct    <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))

  df <- data.frame(
    PC1    = pca$x[, 1],
    PC2    = pca$x[, 2],
    Color  = meta[[color_by]],
    Sample = rownames(meta),
    stringsAsFactors = FALSE
  )

  if (!is.null(shape_by) && shape_by %in% colnames(meta))
    df$Shape <- as.factor(meta[[shape_by]])

  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2,
                                        color = Color, fill = Color))

  min_n <- min(table(df$Color))
  if (ellipse && min_n >= 3) {
    p <- p + ggplot2::stat_ellipse(geom = "polygon", alpha = 0.12, color = NA)
  } else if (ellipse && min_n < 3) {
    message("[plot_omni_pca] Some groups have <3 samples — ellipses suppressed.")
  }

  if (!is.null(shape_by) && "Shape" %in% colnames(df)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(shape = Shape),
                                 size = 4, stroke = 0.8)
  } else {
    p <- p + ggplot2::geom_point(size = 4, stroke = 0.8)
  }

  if (label)
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = Sample),
      size         = 3,
      color        = "black",
      show.legend  = FALSE,
      box.padding  = 0.4,
      max.overlaps = 20
    )

  p <- p +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::scale_fill_viridis_d(option  = "turbo") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title    = "PCA — Covariate Corrected",
      subtitle = paste(
        "Condition:", cond_col,
        if (length(covars) > 0)
          paste("| Regressed:", paste(covars, collapse = ", "))
      ),
      x     = paste0("PC1: ", pct[1], "% variance"),
      y     = paste0("PC2: ", pct[2], "% variance"),
      color = color_by,
      fill  = color_by
    )

  if (!is.null(shape_by) && "Shape" %in% colnames(df))
    p <- p + ggplot2::labs(shape = shape_by)

  print(p)
  return(invisible(p))
}


# =============================================================================
# 2. VOLCANO PLOT
# =============================================================================

#' Plot Volcano Plot
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param lfc_thresh Log2 fold change threshold. Default: `1`.
#' @param pval_thresh Adjusted p-value threshold. Default: `0.05`.
#' @param top_n Integer. Number of top significant genes to label. Default: `15`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_volcano <- function(obj,
                              lfc_thresh  = 1,
                              pval_thresh = 0.05,
                              top_n       = 15) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .require_pkg("ggrepel")

  res <- omni_results(obj)
  if (is.null(res))
    stop("No DE results found. Run run_omni_de() first.", call. = FALSE)

  df      <- as.data.frame(res)
  df$Gene <- rownames(df)
  df      <- df[!is.na(df$log2FoldChange), ]
  df$padj <- ifelse(is.na(df$padj), 1, df$padj)

  df$Direction <- dplyr::case_when(
    df$padj < pval_thresh & df$log2FoldChange >  lfc_thresh ~ "Up",
    df$padj < pval_thresh & df$log2FoldChange < -lfc_thresh ~ "Down",
    TRUE ~ "NS"
  )

  top_up   <- df[df$Direction == "Up",   ]
  top_down <- df[df$Direction == "Down", ]
  top_up   <- utils::head(top_up  [order(top_up$padj),   ], ceiling(top_n / 2))
  top_down <- utils::head(top_down[order(top_down$padj), ], floor(top_n / 2))
  label_df <- rbind(top_up, top_down)

  n_up   <- sum(df$Direction == "Up")
  n_down <- sum(df$Direction == "Down")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange,
                                        y = -log10(padj),
                                        color = Direction)) +
    ggplot2::geom_point(alpha = 0.55, size = 1.2) +
    ggrepel::geom_text_repel(
      data         = label_df,
      ggplot2::aes(label = Gene),
      size         = 3,
      color        = "black",
      box.padding  = 0.5,
      max.overlaps = 30,
      show.legend  = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("Up" = "#E74C3C", "Down" = "#2980B9", "NS" = "grey75"),
      labels = c(
        "Up"   = paste0("Up (", n_up, ")"),
        "Down" = paste0("Down (", n_down, ")"),
        "NS"   = "Not significant"
      )
    ) +
    ggplot2::geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
                        linetype = "dashed", color = "black", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(pval_thresh),
                        linetype = "dashed", color = "black", linewidth = 0.4) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title    = "Volcano Plot",
      subtitle = paste("Contrast:", obj$active_contrast),
      x        = "Log\u2082 Fold Change",
      y        = "-Log\u2081\u2080 Adjusted P-value",
      color    = NULL
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 3. MA PLOT
# =============================================================================

#' Plot MA Plot
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param pval_thresh Adjusted p-value threshold for coloring. Default: `0.05`.
#' @param top_n Integer. Number of top significant genes to label. Default: `10`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_ma <- function(obj,
                         pval_thresh = 0.05,
                         top_n       = 10) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .require_pkg("ggrepel")

  res <- omni_results(obj)
  if (is.null(res))
    stop("No DE results found. Run run_omni_de() first.", call. = FALSE)

  df      <- as.data.frame(res)
  df$Gene <- rownames(df)
  df      <- df[!is.na(df$log2FoldChange) & !is.na(df$baseMean), ]
  df$padj <- ifelse(is.na(df$padj), 1, df$padj)

  df$Significant <- df$padj < pval_thresh
  top_genes      <- utils::head(df[order(df$padj), ], top_n)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = baseMean,
                                        y = log2FoldChange,
                                        color = Significant)) +
    ggplot2::geom_point(alpha = 0.5, size = 1) +
    ggrepel::geom_text_repel(
      data         = top_genes,
      ggplot2::aes(label = Gene),
      size         = 3,
      color        = "black",
      box.padding  = 0.4,
      max.overlaps = 20,
      show.legend  = FALSE
    ) +
    ggplot2::scale_x_log10(labels = scales::comma) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "#E74C3C", "FALSE" = "grey70"),
      labels = c("TRUE" = paste0("padj < ", pval_thresh), "FALSE" = "NS")
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "black", linewidth = 0.6) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title    = "MA Plot",
      subtitle = paste("Contrast:", obj$active_contrast),
      x        = "Mean Expression (baseMean, log scale)",
      y        = "Log\u2082 Fold Change",
      color    = NULL
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 4. HEATMAP
# =============================================================================

#' Plot Annotated Heatmap of Top DEGs
#'
#' Generates a pheatmap of the top significant DEGs from the active contrast,
#' using the cached VST matrix. Automatically subsets samples to only those
#' involved in the active contrast when possible.
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param ntop Integer. Number of top significant genes to plot. Default: `50`.
#' @param pval_thresh Adjusted p-value threshold. Default: `0.05`.
#' @param color_up Color for high expression. Default: `"#B2182B"`.
#' @param color_down Color for low expression. Default: `"#2166AC"`.
#' @param cluster_rows Logical. Cluster genes. Default: `TRUE`.
#' @param cluster_cols Logical. Cluster samples. Default: `TRUE`.
#' @param show_colnames Logical. Show sample names. Default: `TRUE`.
#' @param custom_genes Optional character vector of specific genes to plot.
#' @param title Optional custom title string.
#'
#' @return A pheatmap object (returned invisibly).
#' @export
plot_omni_heatmap <- function(obj,
                              ntop          = 50,
                              pval_thresh   = 0.05,
                              color_up      = "#B2182B",
                              color_down    = "#2166AC",
                              cluster_rows  = TRUE,
                              cluster_cols  = TRUE,
                              show_colnames = TRUE,
                              custom_genes  = NULL,
                              title         = NULL) {

  validate_omniseq(obj)
  .require_pkg("pheatmap")

  res <- omni_results(obj)
  vsd <- omni_vst(obj)

  if (is.null(res)) stop("No DE results. Run run_omni_de() first.", call. = FALSE)
  if (is.null(vsd)) stop("VST not cached. Run run_omni_de() first.", call. = FALSE)

  meta     <- omni_metadata(obj)
  cond_col <- obj$params$condition

  if (!is.null(custom_genes)) {
    gene_names <- intersect(custom_genes,
                            rownames(SummarizedExperiment::assay(vsd)))
    if (length(gene_names) == 0)
      stop("None of custom_genes were found in the VST matrix.", call. = FALSE)

  } else {
    # Always select genes by DE significance ranking from the ACTIVE contrast
    # (which may be a masked contrast — ISGs will be absent from masked results)
    res_df      <- as.data.frame(res)
    res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)
    sig         <- res_df[res_df$padj < pval_thresh, ]

    if (nrow(sig) == 0)
      stop("No significant genes at padj < ", pval_thresh, call. = FALSE)

    sig        <- sig[order(sig$padj), ]
    gene_names <- rownames(sig)[1:min(ntop, nrow(sig))]
  }

  vst_mat <- SummarizedExperiment::assay(vsd)[gene_names, , drop = FALSE]

  samples_to_keep <- .omni_subset_samples(obj, meta, cond_col)
  vst_mat         <- vst_mat[, samples_to_keep, drop = FALSE]
  meta_sub        <- meta[samples_to_keep, , drop = FALSE]

  annot_cols <- intersect(c(cond_col, obj$params$covariates), colnames(meta_sub))
  df_annot   <- meta_sub[, annot_cols, drop = FALSE]
  df_annot[] <- lapply(df_annot,
                       function(x) if (is.factor(x)) droplevels(x) else x)

  palette <- grDevices::colorRampPalette(c(color_down, "white", color_up))(100)

  if (is.null(title))
    title <- paste0("Top ", length(gene_names), " DEGs — ", obj$active_contrast)

  fontsize_row <- dplyr::case_when(
    length(gene_names) > 100 ~ 4,
    length(gene_names) > 60  ~ 5,
    TRUE                     ~ 7
  )

  p <- pheatmap::pheatmap(
    vst_mat,
    annotation_col = df_annot,
    color          = palette,
    scale          = "row",
    cluster_rows   = cluster_rows,
    cluster_cols   = cluster_cols,
    show_colnames  = show_colnames,
    show_rownames  = TRUE,
    border_color   = NA,
    fontsize_row   = fontsize_row,
    main           = title
  )

  return(invisible(p))
}


# =============================================================================
# 5. QC BASICS
# =============================================================================

#' Plot Basic QC Metrics
#'
#' Barplots of library size or detected feature count per sample.
#'
#' @param obj A validated OmniSeq object with QC metrics.
#' @param plot_type Character. Either `"library_size"` or `"detected_genes"`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_qc_basics <- function(obj,
                                plot_type = c("library_size", "detected_genes")) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  qc <- omni_qc(obj)
  if (is.null(qc))
    stop("QC metrics not found. Run run_omni_qc() first.", call. = FALSE)

  plot_type <- match.arg(plot_type)
  df        <- qc$stats
  cond_col  <- obj$params$condition

  if (plot_type == "library_size") {
    p <- ggplot2::ggplot(
      df, ggplot2::aes(x = reorder(Sample, Library_Size),
                       y = Library_Size,
                       fill = .data[[cond_col]])) +
      ggplot2::geom_col(color = "white", linewidth = 0.2) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(title = "Library Size per Sample",
                    x = NULL, y = "Total Assigned Reads",
                    fill = cond_col)
  } else {
    p <- ggplot2::ggplot(
      df, ggplot2::aes(x = reorder(Sample, Detected_Genes),
                       y = Detected_Genes,
                       fill = .data[[cond_col]])) +
      ggplot2::geom_col(color = "white", linewidth = 0.2) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(title = "Detected Features per Sample",
                    x = NULL, y = "Features with Count > 0",
                    fill = cond_col)
  }

  p <- p +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::theme_bw(base_size = 12)

  print(p)
  return(invisible(p))
}


# =============================================================================
# 6. QC ADVANCED
# =============================================================================

#' Plot Advanced QC Diagnostics
#'
#' @param obj A validated OmniSeq object with QC metrics.
#' @param plot_type One of `"distance"`, `"mean_var"`, `"pca_compare"`.
#'
#' @return A plot object (returned invisibly).
#' @export
plot_omni_qc_advanced <- function(obj,
                                  plot_type = c("distance",
                                                "mean_var",
                                                "pca_compare")) {

  validate_omniseq(obj)
  plot_type <- match.arg(plot_type)

  qc   <- omni_qc(obj)
  meta <- omni_metadata(obj)

  if (is.null(qc))
    stop("QC metrics not found. Run run_omni_qc() first.", call. = FALSE)

  if (plot_type == "distance") {
    .require_pkg("pheatmap")

    if (is.null(qc$distances))
      stop("Distance matrix not found. Re-run run_omni_qc(calc_distances = TRUE).",
           call. = FALSE)

    dm           <- qc$distances
    rownames(dm) <- colnames(dm)

    annot_cols <- intersect(c(obj$params$condition, obj$params$covariates),
                            colnames(meta))
    df_annot   <- meta[, annot_cols, drop = FALSE]

    pal <- grDevices::colorRampPalette(
      rev(c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF", "#F7FBFF"))
    )(255)

    p <- pheatmap::pheatmap(
      dm,
      clustering_distance_rows = stats::as.dist(dm),
      clustering_distance_cols = stats::as.dist(dm),
      annotation_col = df_annot,
      annotation_row = df_annot,
      col  = pal,
      main = "Sample-to-Sample Euclidean Distance"
    )
    return(invisible(p))
  }

  if (plot_type == "mean_var") {
    .require_pkg("ggplot2")

    mv <- qc$mv_stats
    p  <- ggplot2::ggplot(mv, ggplot2::aes(x = Mean, y = SD)) +
      ggplot2::geom_point(alpha = 0.15, size = 0.4, color = "#2C3E50") +
      ggplot2::geom_smooth(method    = "loess",
                           formula   = y ~ x,
                           color     = "#E74C3C",
                           se        = FALSE,
                           linewidth = 1) +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::labs(
        title    = "Mean-Variance Trend (VST Normalised)",
        subtitle = "A flat loess line confirms successful variance stabilisation",
        x        = "Mean Expression (VST)",
        y        = "Standard Deviation"
      )
    print(p)
    return(invisible(p))
  }

  if (plot_type == "pca_compare") {
    .require_pkg("ggplot2")
    .require_pkg("patchwork")
    .require_pkg("limma")

    vsd <- omni_vst(obj)
    if (is.null(vsd))
      stop("VST not cached. Run run_omni_qc() or run_omni_de() first.",
           call. = FALSE)

    vst_mat  <- SummarizedExperiment::assay(vsd)
    cond_col <- obj$params$condition
    covars   <- obj$params$covariates

    .make_pca_df <- function(mat) {
      rv  <- apply(mat, 1, stats::var)
      sub <- mat[order(rv, decreasing = TRUE)[1:min(500, nrow(mat))], ]
      pca <- stats::prcomp(t(sub))
      pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
      list(
        df  = data.frame(PC1       = pca$x[, 1],
                         PC2       = pca$x[, 2],
                         Condition = meta[[cond_col]]),
        pct = pct
      )
    }

    before   <- .make_pca_df(vst_mat)
    p_before <- ggplot2::ggplot(
      before$df, ggplot2::aes(PC1, PC2, color = Condition)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::scale_color_viridis_d(option = "turbo") +
      ggplot2::labs(title = "Before correction",
                    x = paste0("PC1: ", before$pct[1], "%"),
                    y = paste0("PC2: ", before$pct[2], "%"))

    mat_corr <- vst_mat
    if (length(covars) > 0) {
      design_bio <- stats::model.matrix(
        stats::as.formula(paste("~", cond_col)), data = meta)
      b1   <- if (length(covars) >= 1) meta[[covars[1]]] else NULL
      b2   <- if (length(covars) >= 2) meta[[covars[2]]] else NULL
      covs <- if (length(covars) >= 3)
        as.matrix(meta[, covars[3:length(covars)], drop = FALSE]) else NULL
      mat_corr <- limma::removeBatchEffect(vst_mat,
                                           batch      = b1,
                                           batch2     = b2,
                                           covariates = covs,
                                           design     = design_bio)
    }

    after   <- .make_pca_df(mat_corr)
    p_after <- ggplot2::ggplot(
      after$df, ggplot2::aes(PC1, PC2, color = Condition)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::scale_color_viridis_d(option = "turbo") +
      ggplot2::labs(
        title = paste("After correction:",
                      paste(covars, collapse = ", ")),
        x = paste0("PC1: ", after$pct[1], "%"),
        y = paste0("PC2: ", after$pct[2], "%")
      )

    p <- p_before + p_after +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = "PCA: Before vs After Covariate Correction"
      )

    print(p)
    return(invisible(p))
  }
}


# =============================================================================
# 7. P-VALUE HISTOGRAM
# =============================================================================

#' Plot P-value Histogram (Model Diagnostic)
#'
#' A well-specified model shows a flat distribution with a spike near 0.
#' A spike near 1 suggests model problems.
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param bins Integer. Number of histogram bins. Default: `40`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_pvalue_hist <- function(obj, bins = 40) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  res <- omni_results(obj)
  if (is.null(res))
    stop("No DE results. Run run_omni_de() first.", call. = FALSE)

  df <- as.data.frame(res)
  df <- df[!is.na(df$pvalue) & !is.na(df$log2FoldChange), ]

  n_sig <- sum(!is.na(df$padj) & df$padj < 0.05)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins      = bins,
                            fill      = "#3182BD",
                            color     = "white",
                            linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = nrow(df) / bins,
                        linetype   = "dashed",
                        color      = "#E74C3C",
                        linewidth  = 0.8) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title    = "P-value Distribution (Model Diagnostic)",
      subtitle = paste0(
        "Contrast: ", obj$active_contrast,
        "  |  Significant (padj < 0.05): ", scales::comma(n_sig),
        "\nDashed line = expected uniform baseline"
      ),
      x = "Raw P-value",
      y = "Count"
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 8. TOP GENES COUNT PLOT
# =============================================================================

#' Plot Normalised Counts for Top DEGs
#'
#' Individual sample points for the top N most significant DEGs, one panel
#' per gene. More honest than a heatmap for communicating single-gene effects.
#' X-axis automatically uses the grouping column relevant to the active
#' contrast — e.g. Genotype for "TLR7 vs UNC93B1", not the primary condition
#' column.
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param n Integer. Number of top genes to plot. Default: `12`.
#' @param pval_thresh Adjusted p-value threshold. Default: `0.05`.
#' @param color_by Character. Metadata column for point color. If `NULL`
#'   (default), auto-detected from the active contrast name.
#' @param shape_by Optional metadata column for point shape.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_top_genes <- function(obj,
                                n           = 12,
                                pval_thresh = 0.05,
                                color_by    = NULL,
                                shape_by    = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  res <- omni_results(obj)
  vsd <- omni_vst(obj)

  if (is.null(res)) stop("No DE results. Run run_omni_de() first.", call. = FALSE)
  if (is.null(vsd)) stop("VST not cached. Run run_omni_de() first.", call. = FALSE)

  meta     <- omni_metadata(obj)
  cond_col <- obj$params$condition

  # FIXED: infer the correct grouping column from the active contrast
  # so that "TLR7 vs UNC93B1" shows Genotype on the x-axis,
  # not the primary Condition (Patient/Control)
  color_by <- color_by %||% .omni_contrast_column(obj, meta, cond_col)

  res_df      <- as.data.frame(res)
  res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)
  sig         <- res_df[res_df$padj < pval_thresh, ]

  if (nrow(sig) == 0)
    stop("No significant genes at padj < ", pval_thresh, call. = FALSE)

  top_genes <- rownames(sig[order(sig$padj), ])[1:min(n, nrow(sig))]
  vst_mat   <- SummarizedExperiment::assay(vsd)[top_genes, , drop = FALSE]

  samples_to_keep <- .omni_subset_samples(obj, meta, cond_col)
  vst_mat         <- vst_mat[, samples_to_keep, drop = FALSE]
  meta_sub        <- meta[samples_to_keep, , drop = FALSE]

  # If color_by column is not present in meta_sub, fall back to cond_col
  if (!color_by %in% colnames(meta_sub)) color_by <- cond_col

  plot_df <- data.frame(
    VST    = as.vector(t(vst_mat)),
    Gene   = rep(top_genes, each = ncol(vst_mat)),
    Sample = rep(colnames(vst_mat), times = length(top_genes)),
    stringsAsFactors = FALSE
  )

  plot_df$Color <- as.character(meta_sub[plot_df$Sample, color_by])
  plot_df$Gene  <- factor(plot_df$Gene, levels = top_genes)

  if (!is.null(shape_by) && shape_by %in% colnames(meta_sub))
    plot_df$Shape <- as.factor(meta_sub[plot_df$Sample, shape_by])

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = Color, y = VST, color = Color))

  if (!is.null(shape_by) && "Shape" %in% colnames(plot_df)) {
    p <- p + ggplot2::geom_jitter(ggplot2::aes(shape = Shape),
                                  width = 0.2, size = 2.5, alpha = 0.85)
  } else {
    p <- p + ggplot2::geom_jitter(width = 0.2, size = 2.5, alpha = 0.85)
  }

  p <- p +
    ggplot2::stat_summary(fun       = mean,
                          geom      = "crossbar",
                          width     = 0.45,
                          linewidth = 0.5,
                          color     = "black",
                          fatten    = 1) +
    ggplot2::facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.text      = ggplot2::element_text(face = "italic", size = 9),
      axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title    = paste("Top", length(top_genes),
                       "DEGs — Normalised Counts (VST)"),
      subtitle = paste("Contrast:", obj$active_contrast),
      x        = color_by,
      y        = "VST-normalised expression"
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Infer the metadata column relevant to the active contrast
#'
#' For a contrast like "TLR7 vs UNC93B1", this returns "Genotype".
#' For the primary contrast or masked contrasts, returns `cond_col`.
#' Used by `plot_omni_top_genes()` to color by the correct grouping variable.
#'
#' @keywords internal
.omni_contrast_column <- function(obj, meta, cond_col) {
  contrast_name <- obj$active_contrast
  if (is.null(contrast_name)) return(cond_col)

  # Strip masked suffix e.g. " (excl. Isgs extended)"
  clean <- gsub("\\s*\\(excl\\..*\\)$", "", contrast_name)
  clean <- trimws(clean)

  # Parse groups from contrast name
  if (grepl(" vs ", clean, fixed = TRUE)) {
    groups <- trimws(strsplit(clean, " vs ", fixed = TRUE)[[1]])
  } else if (grepl("_vs_", clean, fixed = TRUE)) {
    clean2 <- sub(paste0("^", cond_col, "_"), "", clean)
    groups <- trimws(strsplit(clean2, "_vs_", fixed = TRUE)[[1]])
  } else {
    return(cond_col)
  }

  # Search ALL metadata columns for a column that contains both groups
  for (col in colnames(meta)) {
    vals <- unique(as.character(meta[[col]]))
    if (all(groups %in% vals)) return(col)
  }

  return(cond_col)
}


#' Subset sample names to those relevant to the active contrast
#' @keywords internal
.omni_subset_samples <- function(obj, meta, cond_col) {
  contrast_name <- obj$active_contrast
  if (is.null(contrast_name)) return(rownames(meta))

  if (grepl(" vs ", contrast_name, fixed = TRUE)) {
    groups <- trimws(strsplit(contrast_name, " vs ", fixed = TRUE)[[1]])
  } else if (grepl("_vs_", contrast_name, fixed = TRUE)) {
    clean  <- sub(paste0("^", cond_col, "_"), "", contrast_name)
    groups <- trimws(strsplit(clean, "_vs_", fixed = TRUE)[[1]])
  } else {
    return(rownames(meta))
  }

  groups <- gsub("\\s*\\(excl\\..*\\)$", "", groups)

  # Search ALL metadata columns for these groups
  for (col in colnames(meta)) {
    col_vals <- unique(as.character(meta[[col]]))
    if (all(groups %in% col_vals)) {
      return(rownames(meta)[as.character(meta[[col]]) %in% groups])
    }
  }

  return(rownames(meta))
}
