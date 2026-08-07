# =============================================================================
# OmniSeq — 10_TE_Sensors.R
#
# Innate immune sensor activation scoring from TE expression data.
#
# SCIENTIFIC BASIS:
#   Three innate immune sensing pathways are activated by specific TE families:
#
#   1. cGAS-STING (cytoplasmic dsDNA sensing)
#      Activated by L1HS and L1PA2 retrotransposition intermediates (cDNA)
#      and ERV/HERV long terminal repeat elements.
#      Reference: Mathavarajah & Dellaire, Biochem Cell Biol 2024.
#
#   2. RIG-I / MDA5 (cytoplasmic dsRNA sensing)
#      Activated by inverted Alu repeats forming dsRNA structures (AluSx,
#      AluY, AluJb, AluSp, AluSq families). ADAR1 normally edits these;
#      when editing fails, MDA5 is triggered.
#      Reference: Athanasiadis 2021 PLOS Genet; Crow lab, Nat Genet 2014.
#
#   3. TLR7 / TLR8 (endosomal ssRNA sensing)
#      Activated by LINE-1 ssRNA and specific SINE-R elements.
#      Directly relevant to TLR7-mutant patients in this cohort.
#      Reference: EMBO Reports 2025 (STING restricts LINE-1 via lysosomal
#      sorting independent of cGAS).
#
# DESIGN:
#   - Requires a TE OmniSeq object with VST computed
#   - Scores computed as sum of VST expression for relevant loci,
#     then z-scored across samples (so scores are comparable between sensors)
#   - IRG-S score computed from gene OmniSeq object using omniseq_isgs_28
#   - All results stored in gene_obj$sensor_scores
#   - Families matched by substring to accommodate versioned names
#     (e.g. "AluSx1" matches "AluSx")
#
# FUNCTIONS:
#   run_omni_sensor_score()       Compute all three sensor scores + IRG-S
#   plot_omni_sensor_scores()     Violin/boxplot per sensor per group
#   plot_omni_sensor_irgs_corr()  Scatter: each sensor vs IRG-S score
#   plot_omni_sensor_heatmap()    Per-patient score heatmap
# =============================================================================


# =============================================================================
# 1. COMPUTE SENSOR SCORES
# =============================================================================

#' Compute Innate Immune Sensor Activation Scores from TE Expression
#'
#' Computes three per-patient sensor activation scores from TE VST expression:
#' cGAS-STING (cytoplasmic DNA), RIG-I/MDA5 (cytoplasmic dsRNA), and
#' TLR7/TLR8 (endosomal ssRNA). Scores are summed VST expression of relevant
#' TE loci, z-scored across patients for comparability.
#'
#' Also computes the IRG-S score from the gene object (mean VST of the
#' 28-gene interferon response panel) for downstream correlation.
#'
#' @param gene_obj OmniSeq object containing standard gene VST and metadata.
#' @param te_obj OmniSeq object containing TE VST expression.
#' @param group_col Character. Metadata column for grouping in plots.
#'   Default: `"Condition"`.
#' @param exclude_groups Character vector. Groups to exclude from scoring.
#'   Default: `NULL`.
#' @param cgas_families Character vector. TE family names for cGAS-STING
#'   score. Defaults to built-in curated list. Set to `NULL` to use default.
#' @param rigi_families Character vector. TE family names for RIG-I/MDA5
#'   score. Defaults to built-in curated list.
#' @param tlr7_families Character vector. TE family names for TLR7/TLR8
#'   score. Defaults to built-in curated list.
#' @param min_loci Integer. Minimum number of loci required for a sensor
#'   score to be computed. Default: `5`.
#'
#' @return The gene_obj with sensor scores stored in `$sensor_scores`.
#' @export
#'
#' @examples
#' \dontrun{
#' gene_obj <- run_omni_sensor_score(gene_obj, te_obj,
#'                                   group_col = "Phenotype",
#'                                   exclude_groups = "SPEN")
#' plot_omni_sensor_scores(gene_obj)
#' plot_omni_sensor_irgs_corr(gene_obj)
#' }
run_omni_sensor_score <- function(gene_obj,
                                  te_obj,
                                  group_col      = "Condition",
                                  exclude_groups = NULL,
                                  cgas_families  = NULL,
                                  rigi_families  = NULL,
                                  tlr7_families  = NULL,
                                  min_loci       = 5L,
                                  sc_deconv      = NULL) {

  validate_omniseq(gene_obj)
  validate_omniseq(te_obj)

  # ------------------------------------------------------------------
  # 1. DEFAULT TE FAMILY LISTS (from data_signatures.R or hardcoded)
  # ------------------------------------------------------------------
  if (is.null(cgas_families))
    cgas_families <- .omni_default_sensor_families("cgas_sting")
  if (is.null(rigi_families))
    rigi_families <- .omni_default_sensor_families("rigi_mda5")
  if (is.null(tlr7_families))
    tlr7_families <- .omni_default_sensor_families("tlr7")

  # ------------------------------------------------------------------
  # 2. EXTRACT TE MATRIX
  # ------------------------------------------------------------------
  vsd_te <- omni_vst(te_obj)
  if (is.null(vsd_te))
    stop("TE VST not cached. Run run_omni_de() on te_obj first.", call. = FALSE)

  mat_te  <- SummarizedExperiment::assay(vsd_te)
  te_info <- .omni_parse_te_ids(rownames(mat_te))

  message("\n[run_omni_sensor_score] TE matrix: ", nrow(mat_te),
          " loci x ", ncol(mat_te), " samples")

  # ------------------------------------------------------------------
  # 3. COMPUTE PER-SENSOR SCORES
  # ------------------------------------------------------------------
  .compute_score <- function(families, mat, te_info, sensor_name) {
    # Match families by substring — handles AluSx, AluSx1, AluSx2 etc.
    hit_idx <- which(sapply(te_info$Family, function(f) {
      any(sapply(families, function(fam) grepl(fam, f, fixed = TRUE)))
    }))

    if (length(hit_idx) < min_loci) {
      # Suggest the correct TE class object based on which sensor failed
      te_class_hint <- switch(sensor_name,
                              "cGAS-STING" = "Use an ERV or LINE te_obj (these classes contain L1HS, HERVK loci).",
                              "RIG-I/MDA5" = "Use a SINE te_obj (these classes contain Alu loci).",
                              "TLR7/TLR8"  = "Use a LINE or SINE te_obj (L1HS, SVA, Alu loci required).",
                              "Insufficient loci matched."
      )
      warning("[run_omni_sensor_score] Sensor '", sensor_name,
              "' matched only ", length(hit_idx), " loci (min = ", min_loci,
              "). Score will be NA. ", te_class_hint, call. = FALSE)
      return(rep(NA_real_, ncol(mat)))
    }

    message("  ", sensor_name, ": ", length(hit_idx), " loci matched.")

    # Sum VST across matched loci per sample
    raw_score <- colSums(mat[hit_idx, , drop = FALSE], na.rm = TRUE)
    # Z-score across samples
    z_score   <- as.numeric(scale(raw_score))
    names(z_score) <- colnames(mat)
    return(z_score)
  }

  cgas_score <- .compute_score(cgas_families, mat_te, te_info, "cGAS-STING")
  rigi_score <- .compute_score(rigi_families, mat_te, te_info, "RIG-I/MDA5")
  tlr7_score <- .compute_score(tlr7_families, mat_te, te_info, "TLR7/TLR8")

  # ------------------------------------------------------------------
  # 4. IRG-S SCORE FROM GENE OBJECT
  # ------------------------------------------------------------------
  vsd_gene <- omni_vst(gene_obj)
  if (is.null(vsd_gene))
    stop("Gene VST not cached. Run run_omni_de() on gene_obj first.",
         call. = FALSE)

  mat_gene  <- SummarizedExperiment::assay(vsd_gene)
  valid_irgs <- intersect(omniseq_isgs_28, rownames(mat_gene))
  message("  IRG-S: ", length(valid_irgs), "/28 genes found.")

  irgs_raw   <- colMeans(mat_gene[valid_irgs, , drop = FALSE], na.rm = TRUE)
  irgs_score <- as.numeric(scale(irgs_raw))
  names(irgs_score) <- colnames(mat_gene)

  # ------------------------------------------------------------------
  # 5. BUILD SCORE DATA FRAME
  # ------------------------------------------------------------------
  # Common samples across gene and TE objects
  # Use TE matrix column names directly — sensor scores may be NA
  # if a sensor has no matching loci (e.g. cGAS-STING on SINE-only data)
  common_samps <- intersect(colnames(mat_gene), colnames(mat_te))
  if (length(common_samps) < 3)
    stop("Fewer than 3 common samples between gene and TE objects.",
         call. = FALSE)

  meta     <- omni_metadata(gene_obj)

  # Apply group exclusion
  if (!is.null(exclude_groups) && group_col %in% colnames(meta)) {
    excl_samps <- rownames(meta)[
      as.character(meta[[group_col]]) %in% exclude_groups]
    common_samps <- setdiff(common_samps, excl_samps)
    message("  Excluded ", length(excl_samps),
            " samples from groups: ",
            paste(exclude_groups, collapse = ", "))
  }

  score_df <- data.frame(
    Sample      = common_samps,
    cGAS_STING  = cgas_score[common_samps],
    RIG_I_MDA5  = rigi_score[common_samps],
    TLR7_TLR8   = tlr7_score[common_samps],
    IRG_S_Score = irgs_score[common_samps],
    stringsAsFactors = FALSE
  )

  # Attach metadata — always store Group (the grouping column used)
  # AND store all key metadata columns under their original names
  # so plot functions can use any column regardless of what group_col was
  if (group_col %in% colnames(meta))
    score_df$Group <- as.character(meta[common_samps, group_col])

  for (extra_col in c("Genotype", "Phenotype", "Condition", "Batch", "Sex")) {
    if (extra_col %in% colnames(meta))
      score_df[[extra_col]] <- as.character(meta[common_samps, extra_col])
  }

  # ------------------------------------------------------------------
  # 5B. IFN MONOCYTE BURDEN SCORE (scRNA hook — SAVI integration)
  # When sc_deconv is provided (MuSiC deconvolution using SAVI as reference),
  # compute the IFN-like monocyte state burden as a weighted sum of
  # IFN-high monocyte fractions. These sub-states are invisible to xCell
  # but are the key disease-relevant cell population in interferonopathies.
  # ------------------------------------------------------------------
  if (!is.null(sc_deconv)) {
    # sc_deconv should be a matrix/df: samples x cell_states
    # with rownames = sample IDs
    ifn_mono_cols <- grep("IFN|ifn|monocyte|mono|inflammatory",
                          colnames(sc_deconv), ignore.case = TRUE, value = TRUE)
    if (length(ifn_mono_cols) == 0) {
      message("[run_omni_sensor_score] sc_deconv provided but no IFN/monocyte ",
              "columns detected. Skipping IFN monocyte burden score.")
    } else {
      sc_samps <- intersect(common_samps, rownames(sc_deconv))
      if (length(sc_samps) >= 3) {
        ifn_burden_raw <- rowSums(
          as.matrix(sc_deconv[sc_samps, ifn_mono_cols, drop = FALSE]),
          na.rm = TRUE
        )
        ifn_burden_z <- rep(NA_real_, length(common_samps))
        names(ifn_burden_z) <- common_samps
        ifn_burden_z[sc_samps] <- as.numeric(scale(ifn_burden_raw))
        score_df$IFN_Monocyte_Burden <- ifn_burden_z
        message("[run_omni_sensor_score] IFN monocyte burden score computed from ",
                length(ifn_mono_cols), " cell state column(s): ",
                paste(ifn_mono_cols, collapse = ", "))
      } else {
        message("[run_omni_sensor_score] Fewer than 3 common samples between ",
                "sc_deconv and TE data. Skipping IFN monocyte burden score.")
      }
    }
  }

  # ------------------------------------------------------------------
  # 6. STORE
  # ------------------------------------------------------------------
  gene_obj$sensor_scores <- list(
    scores         = score_df,
    group_col      = group_col,
    exclude_groups = exclude_groups,
    cgas_families  = cgas_families,
    rigi_families  = rigi_families,
    tlr7_families  = tlr7_families,
    n_loci         = list(
      cgas_sting = sum(sapply(te_info$Family, function(f)
        any(sapply(cgas_families, function(fam) grepl(fam, f, fixed = TRUE))))),
      rigi_mda5  = sum(sapply(te_info$Family, function(f)
        any(sapply(rigi_families, function(fam) grepl(fam, f, fixed = TRUE))))),
      tlr7       = sum(sapply(te_info$Family, function(f)
        any(sapply(tlr7_families, function(fam) grepl(fam, f, fixed = TRUE)))))
    ),
    n_irgs_genes   = length(valid_irgs),
    # DIABLO-ready block: sensor scores as a sample x feature matrix
    # Pass directly as a block to run_omni_diablo()
    diablo_ready   = TRUE,
    diablo_block   = {
      sensor_cols <- intersect(c("cGAS_STING", "RIG_I_MDA5", "TLR7_TLR8",
                                 "IFN_Monocyte_Burden"),
                               colnames(score_df))
      dm <- as.matrix(score_df[, sensor_cols, drop = FALSE])
      rownames(dm) <- score_df$Sample
      dm[is.na(dm)] <- 0  # replace NA sensors with 0 for DIABLO
      dm
    },
    timestamp      = Sys.time()
  )

  message("\n[run_omni_sensor_score] Scores computed for ",
          nrow(score_df), " samples.")
  message("  Access with gene_obj$sensor_scores$scores")
  message("  Next: plot_omni_sensor_scores(gene_obj)")
  message("        plot_omni_sensor_irgs_corr(gene_obj)")

  return(gene_obj)
}


# =============================================================================
# 2. VIOLIN / BOXPLOT PER SENSOR
# =============================================================================

#' Plot Sensor Score Distributions Across Groups
#'
#' Three-panel violin + boxplot comparing cGAS-STING, RIG-I/MDA5, and
#' TLR7/TLR8 sensor scores across patient groups. Points are individual
#' patients, colored by group.
#'
#' @param obj A validated OmniSeq object with sensor scores.
#' @param group_col Character. Grouping variable for x-axis. If `NULL`,
#'   uses the group_col stored in `$sensor_scores`.
#' @param add_points Logical. Overlay individual sample points. Default: `TRUE`.
#' @param title Optional custom title.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_sensor_scores <- function(obj,
                                    group_col  = NULL,
                                    add_points = TRUE,
                                    title      = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  if (is.null(obj$sensor_scores))
    stop("No sensor scores. Run run_omni_sensor_score() first.", call. = FALSE)

  score_df  <- obj$sensor_scores$scores
  grp_col   <- group_col %||% obj$sensor_scores$group_col %||% "Group"

  # If the requested group_col is not a column name (e.g. "Phenotype"), check
  # if the data frame has a "Group" column storing that grouping
  if (!grp_col %in% colnames(score_df)) {
    if ("Group" %in% colnames(score_df)) {
      message("[plot_omni_sensor_scores] Column '", grp_col,
              "' not found — using 'Group' column (set during run_omni_sensor_score).")
      grp_col <- "Group"
    } else {
      stop("Group column '", grp_col, "' not found in score data frame.",
           call. = FALSE)
    }
  }

  # Reshape to long format
  sensors <- c("cGAS_STING", "RIG_I_MDA5", "TLR7_TLR8")
  sensor_labels <- c(
    "cGAS_STING" = "cGAS-STING\n(cytoplasmic DNA)",
    "RIG_I_MDA5" = "RIG-I / MDA5\n(cytoplasmic dsRNA)",
    "TLR7_TLR8"  = "TLR7 / TLR8\n(endosomal ssRNA)"
  )

  long_df <- do.call(rbind, lapply(sensors, function(s) {
    data.frame(
      Sample = score_df$Sample,
      Group  = score_df[[grp_col]],
      Sensor = sensor_labels[s],
      Score  = score_df[[s]],
      stringsAsFactors = FALSE
    )
  }))
  long_df$Sensor <- factor(long_df$Sensor, levels = sensor_labels)

  if (is.null(title))
    title <- "Innate Immune Sensor Activation Scores"

  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = Group, y = Score, fill = Group)) +
    ggplot2::geom_violin(alpha = 0.4, trim = FALSE, linewidth = 0.5) +
    ggplot2::geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA,
                          linewidth = 0.5)

  if (add_points)
    p <- p + ggplot2::geom_jitter(
      ggplot2::aes(color = Group),
      width = 0.08, size = 2.5, alpha = 0.85, show.legend = FALSE
    )

  p <- p +
    ggplot2::facet_wrap(~ Sensor, scales = "free_y", ncol = 3) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey50", linewidth = 0.4) +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1),
      strip.text       = ggplot2::element_text(face = "bold", size = 10),
      legend.position  = "none"
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = paste("Z-scored VST sum of sensor-specific TE loci |",
                       "Group:", grp_col),
      x        = NULL,
      y        = "Sensor score (z-scored)"
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 3. SENSOR vs IRG-S SCATTER MATRIX
# =============================================================================

#' Plot Sensor Score vs IRG-S Correlation Scatter Matrix
#'
#' Three-panel scatter plot: each innate sensor score vs IRG-S score.
#' Points colored by genotype and shaped by phenotype simultaneously.
#' Includes Pearson R and p-value per panel.
#'
#' @param obj A validated OmniSeq object with sensor scores.
#' @param color_col Character. Metadata column for point color.
#'   Default: `"Genotype"`.
#' @param shape_col Character. Metadata column for point shape.
#'   Default: `"Phenotype"`.
#' @param add_fit Logical. Add linear regression line. Default: `TRUE`.
#' @param title Optional custom title.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_sensor_irgs_corr <- function(obj,
                                       color_col = "Genotype",
                                       shape_col = "Phenotype",
                                       add_fit   = TRUE,
                                       title     = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .require_pkg("ggrepel")

  if (is.null(obj$sensor_scores))
    stop("No sensor scores. Run run_omni_sensor_score() first.", call. = FALSE)

  score_df <- obj$sensor_scores$scores

  sensors <- c("cGAS_STING", "RIG_I_MDA5", "TLR7_TLR8")
  sensor_labels <- c(
    "cGAS_STING" = "cGAS-STING Score",
    "RIG_I_MDA5" = "RIG-I/MDA5 Score",
    "TLR7_TLR8"  = "TLR7/TLR8 Score"
  )

  # Build color and shape columns
  color_var <- if (color_col %in% colnames(score_df)) color_col else
    obj$sensor_scores$group_col %||% "Group"
  shape_var <- if (!is.null(shape_col) && shape_col %in% colnames(score_df))
    shape_col else NULL

  # Reshape to long
  long_df <- do.call(rbind, lapply(sensors, function(s) {
    df <- data.frame(
      Sample    = score_df$Sample,
      IRG_S     = score_df$IRG_S_Score,
      SensorScore = score_df[[s]],
      Sensor    = sensor_labels[s],
      stringsAsFactors = FALSE
    )
    df[[color_var]] <- score_df[[color_var]]
    if (!is.null(shape_var))
      df[[shape_var]] <- score_df[[shape_var]]
    df
  }))
  long_df$Sensor <- factor(long_df$Sensor, levels = sensor_labels)

  # Compute per-panel Pearson R and p-value
  stats_df <- do.call(rbind, lapply(sensors, function(s) {
    sub <- score_df[!is.na(score_df[[s]]) & !is.na(score_df$IRG_S_Score), ]
    if (nrow(sub) < 3) return(NULL)
    ct <- stats::cor.test(sub[[s]], sub$IRG_S_Score, method = "pearson")
    data.frame(
      Sensor = sensor_labels[s],
      R      = round(as.numeric(ct$estimate), 3),
      P      = signif(as.numeric(ct$p.value), 2),
      stringsAsFactors = FALSE
    )
  }))
  stats_df$Label   <- paste0("R = ", stats_df$R, "\np = ", stats_df$P)
  stats_df$Sensor  <- factor(stats_df$Sensor, levels = sensor_labels)

  if (is.null(title))
    title <- "TE Innate Sensor Scores vs IRG-S Score"

  # Base plot
  aes_args <- if (!is.null(shape_var)) {
    ggplot2::aes(x = SensorScore, y = IRG_S,
                 color = .data[[color_var]],
                 shape = .data[[shape_var]])
  } else {
    ggplot2::aes(x = SensorScore, y = IRG_S,
                 color = .data[[color_var]])
  }

  p <- ggplot2::ggplot(long_df, aes_args) +
    ggplot2::geom_point(size = 4, alpha = 0.85)

  if (add_fit)
    p <- p + ggplot2::geom_smooth(method = "lm", color = "black",
                                  linetype = "dashed", fill = "grey80",
                                  formula = y ~ x, linewidth = 0.8,
                                  se = TRUE)

  p <- p +
    ggplot2::geom_text(
      data = stats_df,
      ggplot2::aes(x = -Inf, y = Inf, label = Label),
      inherit.aes = FALSE,
      hjust = -0.1, vjust = 1.3,
      size = 3.5, color = "black", fontface = "italic"
    ) +
    ggplot2::facet_wrap(~ Sensor, scales = "free_x", ncol = 3) +
    ggplot2::scale_color_viridis_d(option = "turbo", name = color_var) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 10)
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = paste("Each point = one patient |",
                       "Color =", color_var,
                       if (!is.null(shape_var)) paste("| Shape =", shape_var)),
      x        = "Sensor score (z-scored)",
      y        = "IRG-S Score (z-scored)"
    )

  if (!is.null(shape_var))
    p <- p + ggplot2::labs(shape = shape_var)

  print(p)
  return(invisible(p))
}


# =============================================================================
# 4. PER-PATIENT SENSOR HEATMAP
# =============================================================================

#' Plot Per-Patient Sensor Score Heatmap
#'
#' Heatmap of all four scores (cGAS-STING, RIG-I/MDA5, TLR7/TLR8, IRG-S)
#' across patients. Annotated by Genotype, Phenotype, Batch, and IRG-S
#' quartile. Reveals which patients have which sensing pathways active.
#'
#' @param obj A validated OmniSeq object with sensor scores.
#' @param cluster_cols Logical. Cluster patients. Default: `TRUE`.
#' @param annotation_cols Character vector. Metadata columns to show as
#'   annotation bars. Default: auto-detected from metadata.
#' @param title Optional custom title.
#'
#' @return A pheatmap object (returned invisibly).
#' @export
plot_omni_sensor_heatmap <- function(obj,
                                     cluster_cols     = TRUE,
                                     annotation_cols  = NULL,
                                     title            = NULL) {

  validate_omniseq(obj)
  .require_pkg("pheatmap")

  if (is.null(obj$sensor_scores))
    stop("No sensor scores. Run run_omni_sensor_score() first.", call. = FALSE)

  score_df <- obj$sensor_scores$scores
  meta     <- omni_metadata(obj)

  # Build the score matrix (rows = sensors, cols = samples)
  score_mat <- t(as.matrix(score_df[, c("cGAS_STING", "RIG_I_MDA5",
                                        "TLR7_TLR8", "IRG_S_Score")]))
  colnames(score_mat) <- score_df$Sample
  rownames(score_mat) <- c("cGAS-STING", "RIG-I / MDA5",
                           "TLR7 / TLR8", "IRG-S Score")

  # Handle NA scores — replace with 0 for display
  score_mat[is.na(score_mat)] <- 0

  # Build annotation
  if (is.null(annotation_cols)) {
    annotation_cols <- intersect(
      c("Genotype", "Phenotype", "Condition", "Batch", "Sex"),
      colnames(meta)
    )
  }

  ann_df <- meta[score_df$Sample, annotation_cols, drop = FALSE]
  ann_df[] <- lapply(ann_df, function(x) {
    if (is.factor(x)) droplevels(x) else x
  })

  # Add IRG-S quartile annotation
  irgs_quartile <- cut(score_df$IRG_S_Score,
                       breaks = stats::quantile(score_df$IRG_S_Score,
                                                probs = 0:4 / 4, na.rm = TRUE),
                       labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
                       include.lowest = TRUE)
  ann_df$IRG_S_Quartile <- irgs_quartile

  pal <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

  if (is.null(title))
    title <- "Innate Sensor Scores per Patient"

  p <- pheatmap::pheatmap(
    score_mat,
    annotation_col   = ann_df,
    color            = pal,
    cluster_rows     = FALSE,
    cluster_cols     = cluster_cols,
    border_color     = NA,
    show_colnames    = TRUE,
    fontsize_col     = 8,
    fontsize_row     = 11,
    cellheight       = 30,
    main             = title
  )

  return(invisible(p))
}


# =============================================================================
# INTERNAL: DEFAULT SENSOR FAMILY LISTS
# =============================================================================

#' Return default TE family lists for each innate sensor
#'
#' These lists are curated from the literature and match the family names
#' used in RepeatMasker / Telescope annotations (the format used in
#' locus IDs like `100000:SINE:AluSx`).
#'
#' @keywords internal
.omni_default_sensor_families <- function(sensor) {

  switch(sensor,

         # cGAS-STING: activated by LINE-1 retrotransposition cDNA intermediates
         # and ERV/HERV elements producing cytoplasmic dsDNA.
         # Primary families: L1HS (most active), L1PA2, L1PA3 (recently active),
         # HERV-K (HML-2), HERV-H, HERV-W (cytoplasmic DNA from LTR elements).
         # Source: Mathavarajah & Dellaire 2024; De Cecco et al. Nature 2019.
         "cgas_sting" = c(
           "L1HS", "L1PA2", "L1PA3", "L1PA4",
           "HERVK", "HERV-K", "HML-2",
           "HERVH", "HERV-H",
           "HERVW", "HERV-W", "HERV17",
           "LTR12", "LTR7", "LTR5"
         ),

         # RIG-I / MDA5: activated by inverted Alu repeat pairs forming dsRNA.
         # Alu families with high dsRNA-forming potential: AluSx, AluY, AluJb,
         # AluSp, AluSq, AluJr. MIR elements also contribute dsRNA in some contexts.
         # Source: Ahmad et al. Cell 2018; Athanasiadis PLOS Genet 2021;
         #         Crow lab Nat Genet 2014 (ADAR/MDA5 AGS).
         "rigi_mda5" = c(
           "AluSx", "AluY", "AluJb", "AluSp", "AluSq", "AluJr",
           "AluSc", "AluSg", "AluSz", "AluPa",
           "MIRb", "MIRc", "MIR3",
           "ERVL", "ERVL-MaLR"
         ),

         # TLR7 / TLR8: activated by ssRNA in endosomes.
         # LINE-1 RNA (L1Md_A, L1Md_T in mice; L1HS, L1PA in humans),
         # SVA elements (SINE-R/VNTR/Alu), and specific SINE families.
         # Directly relevant to TLR7-mutant patients: TLR7 recognises G- and
         # U-rich ssRNA including LINE-1 transcripts.
         # Source: EMBO Reports 2025 (LINE-1 RNA sensed by TLR7);
         #         Young et al. 2012.
         "tlr7" = c(
           "L1HS", "L1PA2", "L1PA3",
           "SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F",
           "AluSx", "AluY",
           "MIRb", "MIR3",
           "THE1", "THE1B", "THE1C"
         ),

         stop("Unknown sensor: '", sensor,
              "'. Use 'cgas_sting', 'rigi_mda5', or 'tlr7'.", call. = FALSE)
  )
}
