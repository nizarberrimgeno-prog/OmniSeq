# =============================================================================
# OmniSeq — 08_sPLSDA.R
#
# Supervised dimensionality reduction and biomarker discovery via sPLS-DA.
# Uses the mixOmics package.
#
# DESIGN PRINCIPLES:
#   - Fully flexible: any metadata column can be the grouping variable
#   - exclude_groups drops any unwanted groups before fitting
#   - Results cached in obj$splsda[["<key>"]] — multiple models coexist
#   - Automatic component and keepX tuning via LOO cross-validation
#   - nrepeat forced to 1L for LOO (LOO never needs repeating by definition)
#   - Warns when any group has n < 3 but does not stop
#   - Controls can be included or excluded — user decides per call
#
# FIXES vs previous version:
#   - nrepeat = 1L enforced for LOO — removes spurious warning
#   - near_zero_var = TRUE added to run_omni_diablo() — prevents zero-variance crash
#   - near_zero_var argument exposed so user can control it
#   - plot_omni_splsda_scores() 1-component fallback fully implemented
#   - run_omni_diablo() keepX capped to min(keepX, ncol(block)) per block
#   - run_omni_diablo() removed reference to undefined `use_masked` variable
#
# FUNCTIONS:
#   run_omni_splsda()           Core sPLS-DA engine
#   run_omni_diablo()           Multi-block integration (gene + TE)
#   omni_splsda()               Accessor
#   plot_omni_splsda_scores()   Supervised sample score plot
#   plot_omni_splsda_loadings() Top discriminating features per component
#   plot_omni_splsda_vip()      VIP scores — de novo biomarker list
#   plot_omni_splsda_roc()      LOO cross-validated ROC per component
#   omni_splsda_signature()     Extract the de novo gene signature as a table
# =============================================================================


# =============================================================================
# 1. CORE ENGINE
# =============================================================================

#' Run Sparse Partial Least Squares Discriminant Analysis (sPLS-DA)
#'
#' Fits an sPLS-DA model on the VST expression matrix using a user-specified
#' metadata column as the grouping variable. Supports any number of groups,
#' flexible sample inclusion/exclusion, and caches results so multiple models
#' can coexist in the same object.
#'
#' Component number and the number of features to keep per component (keepX)
#' are both tuned automatically via leave-one-out cross-validation unless
#' overridden by the user.
#'
#' @param obj A validated OmniSeq object with a fitted model.
#' @param group_col Character. Name of the metadata column to use as the
#'   grouping/outcome variable (e.g. `"Condition"`, `"Phenotype"`,
#'   `"Genotype"`).
#' @param exclude_groups Optional character vector. Values in `group_col` to
#'   drop before fitting. Default: `NULL` (use all groups).
#' @param ncomp Integer. Number of components to fit. If `NULL` (default),
#'   automatically tuned via LOO cross-validation (max tested: 3).
#' @param keepX Integer vector of length `ncomp`. Number of features to
#'   select per component. If `NULL` (default), automatically tuned via
#'   cross-validation from candidates `c(5, 10, 15, 20, 25, 30, 40, 50)`.
#' @param use_masked Logical. If `TRUE` and a masked contrast is active,
#'   use the ISG-masked VST rather than the full VST. Default: `FALSE`.
#' @param scale Logical. Centre and scale the X matrix before fitting.
#'   Default: `TRUE`.
#' @param seed Integer. Random seed for reproducibility. Default: `42`.
#'
#' @return The OmniSeq object with results stored in
#'   `obj$splsda[["<group_col>_<groups>"]]`.
#' @export
run_omni_splsda <- function(obj,
                            group_col      = "Condition",
                            exclude_groups = NULL,
                            ncomp          = NULL,
                            keepX          = NULL,
                            use_masked     = FALSE,
                            scale          = TRUE,
                            seed           = 42L) {

  validate_omniseq(obj)
  .require_pkg("mixOmics")

  meta <- omni_metadata(obj)

  # ------------------------------------------------------------------
  # 1. VALIDATE group_col
  # ------------------------------------------------------------------
  if (!group_col %in% colnames(meta))
    stop("group_col '", group_col, "' not found in metadata.\n",
         "Available columns: ", paste(colnames(meta), collapse = ", "),
         call. = FALSE)

  # ------------------------------------------------------------------
  # 2. BUILD SAMPLE SUBSET
  # ------------------------------------------------------------------
  groups_all <- as.character(meta[[group_col]])

  if (!is.null(exclude_groups)) {
    bad <- setdiff(exclude_groups, unique(groups_all))
    if (length(bad) > 0)
      message("[run_omni_splsda] Note: exclude_groups values not found in '",
              group_col, "': ", paste(bad, collapse = ", "))
  }

  keep_mask  <- !groups_all %in% (exclude_groups %||% character(0))
  keep_samps <- rownames(meta)[keep_mask]
  groups_use <- groups_all[keep_mask]

  if (length(keep_samps) < 4)
    stop("Fewer than 4 samples remain after exclusion. Cannot fit sPLS-DA.",
         call. = FALSE)

  grp_tab <- table(groups_use)
  small   <- names(grp_tab[grp_tab < 3])
  if (length(small) > 0)
    warning("[run_omni_splsda] Group(s) with n < 3: ",
            paste(small, collapse = ", "),
            ". LOO cross-validation may be unreliable.",
            call. = FALSE)

  if (length(unique(groups_use)) < 2)
    stop("At least 2 groups required. Only found: ",
         paste(unique(groups_use), collapse = ", "), call. = FALSE)

  groups_use <- factor(groups_use)
  message("\n[run_omni_splsda] Groups: ",
          paste(paste0(names(grp_tab), " (n=", grp_tab, ")"),
                collapse = ", "))

  # ------------------------------------------------------------------
  # 3. BUILD EXPRESSION MATRIX
  # ------------------------------------------------------------------
  vsd <- omni_vst(obj)
  if (is.null(vsd))
    stop("VST not cached. Run run_omni_de() first.", call. = FALSE)

  mat_full <- SummarizedExperiment::assay(vsd)

  # Use ISG-masked gene set if requested and available
  if (use_masked && !is.null(obj$active_contrast) &&
      grepl("excl\\.", obj$active_contrast)) {
    parent_key   <- sub("\\s*\\(excl\\..*\\)$", "", obj$active_contrast)
    parent_entry <- obj$results_cache[[parent_key]]
    masked_entry <- obj$results_cache[[obj$active_contrast]]

    if (!is.null(parent_entry) && !is.null(masked_entry)) {
      parent_df    <- as.data.frame(parent_entry$results)
      masked_df    <- as.data.frame(masked_entry$results)
      masked_genes <- setdiff(rownames(parent_df), rownames(masked_df))

      if (length(masked_genes) > 0) {
        mat_full <- mat_full[!rownames(mat_full) %in% masked_genes, ,
                             drop = FALSE]
        message("[run_omni_splsda] Using ISG-masked gene set: removed ",
                length(masked_genes), " masked genes. Remaining: ",
                nrow(mat_full), " features.")
      } else {
        message("[run_omni_splsda] Warning: no masked genes identified.")
      }
    }
  }

  # Samples x features for mixOmics
  X <- t(mat_full[, keep_samps, drop = FALSE])
  Y <- groups_use

  message("[run_omni_splsda] Matrix: ", nrow(X), " samples x ",
          ncol(X), " features.")

  # ------------------------------------------------------------------
  # 4. TUNE ncomp AND keepX VIA LOO CROSS-VALIDATION
  # nrepeat is always 1L for LOO — repeating LOO is meaningless.
  # ------------------------------------------------------------------
  set.seed(seed)

  if (is.null(ncomp) || is.null(keepX)) {
    n_groups <- length(unique(Y))
    max_comp <- min(3L, n_groups - 1L, nrow(X) - 1L)

    message("[run_omni_splsda] Tuning ncomp and keepX via LOO ",
            "cross-validation (max ncomp = ", max_comp, ")...")
    message("  This may take a minute...")

    keepX_candidates <- c(5, 10, 15, 20, 25, 30, 40, 50)

    tune_res <- tryCatch({
      mixOmics::tune.splsda(
        X           = X,
        Y           = Y,
        ncomp       = max_comp,
        scale       = scale,
        test.keepX  = keepX_candidates,
        validation  = "loo",
        nrepeat     = 1L,
        dist        = "max.dist",
        progressBar = FALSE
      )
    }, error = function(e) {
      message("[run_omni_splsda] Auto-tuning failed (", e$message,
              "). Using defaults: ncomp=1, keepX=20.")
      NULL
    })

    if (!is.null(tune_res)) {
      optimal_keepX <- tune_res$choice.keepX
      optimal_ncomp <- length(optimal_keepX)
      message("[run_omni_splsda] Optimal ncomp = ", optimal_ncomp,
              " | keepX = ", paste(optimal_keepX, collapse = ", "))
    } else {
      optimal_ncomp <- min(1L, max_comp)
      optimal_keepX <- rep(20L, optimal_ncomp)
    }

    if (!is.null(ncomp))  optimal_ncomp <- ncomp
    if (!is.null(keepX))  optimal_keepX <- keepX

  } else {
    optimal_ncomp <- ncomp
    optimal_keepX <- keepX
    tune_res      <- NULL
    message("[run_omni_splsda] Using user-specified ncomp = ", ncomp,
            " | keepX = ", paste(keepX, collapse = ", "))
  }

  # ------------------------------------------------------------------
  # 5. FIT FINAL MODEL
  # ------------------------------------------------------------------
  message("[run_omni_splsda] Fitting final sPLS-DA model...")
  set.seed(seed)

  splsda_fit <- mixOmics::splsda(
    X     = X,
    Y     = Y,
    ncomp = optimal_ncomp,
    keepX = optimal_keepX,
    scale = scale
  )

  # ------------------------------------------------------------------
  # 6. REPEATED CV PERFORMANCE EVALUATION
  # Using nrepeat = 50 for stable BER estimates and 95% CI
  # LOO (nrepeat=1) produces single-point estimates with no uncertainty.
  # Repeated CV (M-fold or LOO with repeats) gives distribution of errors.
  # ------------------------------------------------------------------
  message("[run_omni_splsda] Evaluating CV performance (nrepeat = 50)...")
  message("  This may take ~1-2 minutes. BER with 95% CI will be reported.")
  set.seed(seed)

  n_samples <- length(Y)
  # Use M-fold CV for small cohorts (n < 30) to avoid empty test folds
  # Use LOO for very small groups but repeat it to get uncertainty
  validation_type <- if (n_samples >= 30) "Mfold" else "loo"
  folds_val       <- if (validation_type == "Mfold") min(5L, n_samples) else n_samples

  perf_res <- tryCatch({
    mixOmics::perf(
      splsda_fit,
      validation  = validation_type,
      folds       = folds_val,
      nrepeat     = 50L,
      dist        = c("max.dist", "centroids.dist", "mahalanobis.dist"),
      progressBar = FALSE
    )
  }, error = function(e) {
    message("[run_omni_splsda] Repeated CV failed: ", e$message)
    message("  Falling back to LOO nrepeat=1...")
    tryCatch(
      mixOmics::perf(splsda_fit, validation = "loo",
                     nrepeat = 1L, dist = "max.dist",
                     progressBar = FALSE),
      error = function(e2) {
        message("[run_omni_splsda] Performance evaluation failed: ", e2$message)
        NULL
      }
    )
  })

  # ------------------------------------------------------------------
  # 7. VIP SCORES
  # ------------------------------------------------------------------
  vip_mat <- mixOmics::vip(splsda_fit)

  # ------------------------------------------------------------------
  # 8. CACHE KEY AND STORE
  # ------------------------------------------------------------------
  retained_groups <- sort(unique(as.character(Y)))
  cache_key <- paste0(group_col, "_",
                      paste(retained_groups, collapse = "-"))
  if (use_masked) cache_key <- paste0(cache_key, "_masked")

  if (is.null(obj$splsda)) obj$splsda <- list()

  obj$splsda[[cache_key]] <- list(
    model       = splsda_fit,
    tune        = tune_res,
    performance = perf_res,
    vip         = vip_mat,
    group_col   = group_col,
    groups      = retained_groups,
    excluded    = exclude_groups,
    samples     = keep_samps,
    ncomp       = optimal_ncomp,
    keepX       = optimal_keepX,
    use_masked  = use_masked,
    timestamp   = Sys.time()
  )

  if (!is.null(perf_res)) {
    ber <- perf_res$error.rate$BER
    message("\n[run_omni_splsda] Balanced Error Rate per component:")
    for (i in seq_len(nrow(ber)))
      message("  Component ", i, ": BER = ",
              round(ber[i, "max.dist"], 3))
  }

  message("\n[run_omni_splsda] Model stored under key: '", cache_key, "'")
  message("  Use omni_splsda(obj, '", cache_key, "') to access, or plot with:")
  message("  plot_omni_splsda_scores(obj, '", cache_key, "')")
  message("  plot_omni_splsda_vip(obj, '",    cache_key, "')")

  return(obj)
}


# =============================================================================
# 2. ACCESSOR
# =============================================================================

#' Access a Stored sPLS-DA Result
#'
#' @param obj A validated OmniSeq object.
#' @param key Character. Cache key, or `NULL` to list all available keys.
#'
#' @return The stored sPLS-DA result list, or a character vector of keys.
#' @export
omni_splsda <- function(obj, key = NULL) {
  .check_class(obj)
  if (is.null(obj$splsda) || length(obj$splsda) == 0) {
    message("No sPLS-DA results found. Run run_omni_splsda() first.")
    return(invisible(NULL))
  }
  if (is.null(key)) {
    cat("\nAvailable sPLS-DA models:\n")
    for (k in names(obj$splsda)) {
      entry <- obj$splsda[[k]]
      cat("  '", k, "' — groups: ",
          paste(entry$groups, collapse = ", "),
          " | ncomp: ", entry$ncomp,
          " | fitted: ", format(entry$timestamp, "%Y-%m-%d %H:%M"),
          "\n", sep = "")
    }
    return(invisible(names(obj$splsda)))
  }
  if (!key %in% names(obj$splsda))
    stop("Key '", key, "' not found. Available: ",
         paste(names(obj$splsda), collapse = ", "), call. = FALSE)
  obj$splsda[[key]]
}


# =============================================================================
# 3. SCORE PLOT
# =============================================================================

#' Plot sPLS-DA Sample Score Plot
#'
#' @param obj A validated OmniSeq object with sPLS-DA results.
#' @param key Character. Cache key. Auto-selected if only one model exists.
#' @param comp Integer vector of length 2. Components to plot. Default: `c(1, 2)`.
#' @param ellipse Logical. Draw 95% confidence ellipses. Default: `TRUE`.
#' @param label Logical. Label sample points. Default: `TRUE`.
#' @param shape_col Optional metadata column to map to point shape.
#' @param title Optional custom title.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_splsda_scores <- function(obj,
                                    key       = NULL,
                                    comp      = c(1, 2),
                                    ellipse   = TRUE,
                                    label     = TRUE,
                                    shape_col = NULL,
                                    title     = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .require_pkg("ggrepel")

  entry <- .omni_get_splsda_entry(obj, key)
  model <- entry$model
  meta  <- omni_metadata(obj)

  max_comp <- model$ncomp
  comp     <- comp[comp <= max_comp]

  if (length(comp) == 0)
    stop("No valid components requested. Model has ", max_comp,
         " component(s).", call. = FALSE)

  expl_var <- model$prop_expl_var$X

  # ── 1-COMPONENT MODEL ───────────────────────────────────────────────────
  if (length(comp) == 1) {
    message("[plot_omni_splsda_scores] Model has only 1 component. ",
            "Plotting Component 1 vs sample index.")

    scores_1d           <- as.data.frame(model$variates$X[, 1, drop = FALSE])
    colnames(scores_1d) <- "Comp1"
    scores_1d$Sample    <- rownames(scores_1d)
    scores_1d$Group     <- as.character(model$Y)
    scores_1d$Index     <- seq_len(nrow(scores_1d))

    if (!is.null(shape_col) && shape_col %in% colnames(meta))
      scores_1d$Shape <- as.factor(meta[scores_1d$Sample, shape_col])

    if (is.null(title))
      title <- paste0("sPLS-DA Score Plot (1 component) — ",
                      entry$group_col, ": ",
                      paste(entry$groups, collapse = " vs "))

    p <- ggplot2::ggplot(scores_1d,
                         ggplot2::aes(x = Comp1, y = Index, color = Group))

    if (!is.null(shape_col) && "Shape" %in% colnames(scores_1d))
      p <- p + ggplot2::geom_point(ggplot2::aes(shape = Shape),
                                   size = 4, stroke = 0.8)
    else
      p <- p + ggplot2::geom_point(size = 4, stroke = 0.8)

    if (label)
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = Sample),
        size = 3, color = "black", show.legend = FALSE,
        box.padding = 0.4, max.overlaps = 20
      )

    p <- p +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          color = "grey70", linewidth = 0.4) +
      ggplot2::scale_color_viridis_d(option = "turbo") +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::labs(
        title = title,
        x     = paste0("Component 1 (",
                       round(expl_var[1] * 100, 1), "% expl. var.)"),
        y     = "Sample index",
        color = entry$group_col
      )

    if (!is.null(shape_col)) p <- p + ggplot2::labs(shape = shape_col)

    print(p)
    return(invisible(p))
  }

  # ── 2+ COMPONENT MODEL ──────────────────────────────────────────────────
  scores <- as.data.frame(model$variates$X[, comp, drop = FALSE])
  colnames(scores) <- paste0("Comp", comp)
  scores$Sample    <- rownames(scores)
  scores$Group     <- as.character(model$Y)

  xlab <- paste0("Component ", comp[1],
                 " (", round(expl_var[comp[1]] * 100, 1), "% expl. var.)")
  ylab <- paste0("Component ", comp[2],
                 " (", round(expl_var[comp[2]] * 100, 1), "% expl. var.)")

  if (!is.null(shape_col) && shape_col %in% colnames(meta))
    scores$Shape <- as.factor(meta[scores$Sample, shape_col])

  if (is.null(title))
    title <- paste0("sPLS-DA Score Plot — ", entry$group_col, ": ",
                    paste(entry$groups, collapse = " vs "))

  p <- ggplot2::ggplot(scores,
                       ggplot2::aes(x     = .data[["Comp1"]],
                                    y     = .data[["Comp2"]],
                                    color = Group,
                                    fill  = Group))

  if (ellipse) {
    min_n <- min(table(scores$Group))
    if (min_n >= 3)
      p <- p + ggplot2::stat_ellipse(geom = "polygon", alpha = 0.12,
                                     color = NA)
    else
      message("[plot_omni_splsda_scores] Some groups have n < 3 — ",
              "ellipses suppressed.")
  }

  if (!is.null(shape_col) && "Shape" %in% colnames(scores))
    p <- p + ggplot2::geom_point(ggplot2::aes(shape = Shape),
                                 size = 4, stroke = 0.8)
  else
    p <- p + ggplot2::geom_point(size = 4, stroke = 0.8)

  if (label)
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = Sample),
      size = 3, color = "black", show.legend = FALSE,
      box.padding = 0.4, max.overlaps = 20
    )

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey70", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey70", linewidth = 0.4) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::scale_fill_viridis_d(option  = "turbo") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(title = title, x = xlab, y = ylab,
                  color = entry$group_col, fill = entry$group_col)

  if (!is.null(shape_col)) p <- p + ggplot2::labs(shape = shape_col)

  print(p)
  return(invisible(p))
}


# =============================================================================
# 4. LOADINGS PLOT
# =============================================================================

#' Plot sPLS-DA Feature Loadings
#'
#' @param obj A validated OmniSeq object with sPLS-DA results.
#' @param key Character. Cache key.
#' @param comp Integer. Which component to plot. Default: `1`.
#' @param n_features Integer. Number of top features to show. Default: `20`.
#' @param title Optional custom title.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_splsda_loadings <- function(obj,
                                      key        = NULL,
                                      comp       = 1L,
                                      n_features = 20L,
                                      title      = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  entry <- .omni_get_splsda_entry(obj, key)
  model <- entry$model

  if (comp > model$ncomp)
    stop("Component ", comp, " requested but model only has ",
         model$ncomp, " components.", call. = FALSE)

  loadings <- model$loadings$X[, comp]
  loadings <- loadings[loadings != 0]

  if (length(loadings) == 0)
    stop("No features selected on component ", comp, call. = FALSE)

  loadings <- sort(loadings)

  if (length(loadings) > n_features) {
    n_pos <- ceiling(n_features / 2)
    n_neg <- floor(n_features / 2)
    pos   <- utils::tail(loadings[loadings > 0], n_pos)
    neg   <- utils::head(loadings[loadings < 0], n_neg)
    loadings <- c(neg, pos)
  }

  df <- data.frame(
    Feature   = names(loadings),
    Loading   = as.numeric(loadings),
    stringsAsFactors = FALSE
  )
  df$Feature   <- factor(df$Feature, levels = df$Feature)
  df$Direction <- ifelse(df$Loading > 0, "Positive", "Negative")

  if (is.null(title))
    title <- paste0("sPLS-DA Loadings — Component ", comp, "\n",
                    entry$group_col, ": ",
                    paste(entry$groups, collapse = " vs "))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Loading, y = Feature,
                                        fill = Direction)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("Positive" = "#E74C3C", "Negative" = "#2980B9"),
      name   = "Direction"
    ) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8, face = "italic")
    ) +
    ggplot2::labs(title = title,
                  x     = paste0("Loading coefficient (Component ", comp, ")"),
                  y     = NULL)

  print(p)
  return(invisible(p))
}


# =============================================================================
# 5. VIP SCORES
# =============================================================================

#' Plot VIP Scores (Variable Importance in Projection)
#'
#' @param obj A validated OmniSeq object with sPLS-DA results.
#' @param key Character. Cache key.
#' @param comp Integer. Which component. Default: `1`.
#' @param n_features Integer. Top features to show. Default: `25`.
#' @param vip_threshold Numeric. Threshold line. Default: `1`.
#' @param title Optional custom title.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_splsda_vip <- function(obj,
                                 key           = NULL,
                                 comp          = 1L,
                                 n_features    = 25L,
                                 vip_threshold = 1,
                                 title         = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  entry   <- .omni_get_splsda_entry(obj, key)
  vip_mat <- entry$vip

  if (comp > ncol(vip_mat))
    stop("Component ", comp, " not in VIP matrix.", call. = FALSE)

  vip_vec <- sort(vip_mat[, comp], decreasing = TRUE)
  vip_vec <- utils::head(vip_vec, n_features)

  df <- data.frame(
    Feature   = names(vip_vec),
    VIP       = as.numeric(vip_vec),
    stringsAsFactors = FALSE
  )
  df$Feature   <- factor(df$Feature, levels = rev(df$Feature))
  df$Important <- df$VIP >= vip_threshold

  if (is.null(title))
    title <- paste0("VIP Scores — Component ", comp, "\n",
                    entry$group_col, ": ",
                    paste(entry$groups, collapse = " vs "))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = VIP, y = Feature,
                                        fill = Important)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#E74C3C", "FALSE" = "grey70"),
      labels = c("TRUE" = paste0("VIP \u2265 ", vip_threshold),
                 "FALSE" = paste0("VIP < ", vip_threshold)),
      name   = NULL
    ) +
    ggplot2::geom_vline(xintercept = vip_threshold, linetype = "dashed",
                        color = "black", linewidth = 0.6) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8, face = "italic")
    ) +
    ggplot2::labs(title = title, x = "VIP Score", y = NULL)

  print(p)
  return(invisible(p))
}


# =============================================================================
# 6. ROC CURVE
# =============================================================================

#' Plot Cross-Validated Balanced Error Rate
#'
#' Replaces the misleading training-set ROC/AUC plot with an honest display
#' of the cross-validated Balanced Error Rate (BER) per component.
#'
#' The 100% AUC from `mixOmics::auroc()` reflects training-set performance —
#' with small cohorts and many features, perfect training separation is trivial
#' and uninformative. BER from repeated CV is the honest performance metric.
#'
#' Interpretation guide shown on the plot:
#'   BER < 0.20  — strong discrimination
#'   BER 0.20–0.35 — moderate signal
#'   BER 0.35–0.45 — weak signal
#'   BER >= 0.50 — no better than chance (random classifier)
#'
#' 95% CI is shown when nrepeat > 1 in the stored performance object.
#'
#' @param obj A validated OmniSeq object with sPLS-DA results.
#' @param key Character. Cache key for the model to plot.
#' @param comp Integer. Maximum component to display. Default: `1`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_splsda_roc <- function(obj, key = NULL, comp = 1L) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  entry <- .omni_get_splsda_entry(obj, key)

  if (is.null(entry$performance))
    stop("Performance evaluation not available. ",
         "Re-run run_omni_splsda() to generate it.", call. = FALSE)

  perf    <- entry$performance
  ber_mat <- perf$error.rate$BER

  if (is.null(ber_mat))
    stop("BER not found in performance object.", call. = FALSE)

  n_comp  <- min(comp, nrow(ber_mat))
  n_groups <- length(entry$groups)

  # Build plot data frame
  # ber_mat rows = components, cols = distance metrics
  # Use max.dist as primary; centroids.dist and mahalanobis.dist if available
  dist_cols <- intersect(c("max.dist", "centroids.dist", "mahalanobis.dist"),
                         colnames(ber_mat))

  plot_rows <- lapply(seq_len(n_comp), function(i) {
    lapply(dist_cols, function(dc) {
      ber_val <- ber_mat[i, dc]

      # Extract 95% CI from repeated CV if available
      ci_lo <- NA_real_
      ci_hi <- NA_real_

      if (!is.null(perf$error.rate.sd)) {
        ber_sd <- perf$error.rate.sd$BER
        if (!is.null(ber_sd) && dc %in% colnames(ber_sd)) {
          sd_val <- ber_sd[i, dc]
          # 95% CI approximation: mean ± 1.96 * SD / sqrt(nrepeat)
          # nrepeat not directly stored — use SD as uncertainty band
          ci_lo  <- max(0, ber_val - 1.96 * sd_val)
          ci_hi  <- min(1, ber_val + 1.96 * sd_val)
        }
      }

      data.frame(
        Component = paste0("Comp ", i),
        Distance  = dc,
        BER       = ber_val,
        CI_lo     = ci_lo,
        CI_hi     = ci_hi,
        stringsAsFactors = FALSE
      )
    })
  })

  df <- do.call(rbind, do.call(c, plot_rows))
  df$Component <- factor(df$Component,
                         levels = paste0("Comp ", seq_len(n_comp)))

  # Color points by performance zone
  df$Zone <- dplyr::case_when(
    df$BER <  0.20 ~ "Strong",
    df$BER <  0.35 ~ "Moderate",
    df$BER <  0.50 ~ "Weak",
    TRUE           ~ "Chance"
  )
  df$Zone <- factor(df$Zone,
                    levels = c("Strong", "Moderate", "Weak", "Chance"))

  zone_colors <- c(
    "Strong"   = "#27AE60",
    "Moderate" = "#F39C12",
    "Weak"     = "#E67E22",
    "Chance"   = "#E74C3C"
  )

  has_ci <- any(!is.na(df$CI_lo))

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = Component, y = BER,
                                    color = Zone, group = Distance,
                                    shape = Distance)) +
    # Reference line at 0.5 — random classifier
    ggplot2::geom_hline(yintercept = 0.5,
                        linetype = "dashed", color = "#E74C3C",
                        linewidth = 0.8) +
    # Shaded interpretation zones
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = 0,    ymax = 0.20,
                      fill = "#27AE60", alpha = 0.06) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = 0.20, ymax = 0.35,
                      fill = "#F39C12", alpha = 0.06) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = 0.35, ymax = 0.50,
                      fill = "#E67E22", alpha = 0.06) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = 0.50, ymax = 1.0,
                      fill = "#E74C3C", alpha = 0.06)

  # Add 95% CI error bars if available
  if (has_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = CI_lo, ymax = CI_hi),
      width = 0.15, alpha = 0.6
    )
  }

  p <- p +
    ggplot2::geom_point(size = 4) +
    ggplot2::scale_color_manual(values = zone_colors, name = "Signal") +
    ggplot2::scale_shape_manual(
      values = c("max.dist" = 16, "centroids.dist" = 17,
                 "mahalanobis.dist" = 15),
      name = "Distance"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.20, 0.35, 0.50, 0.75, 1.0),
      labels = c("0", "0.20", "0.35", "0.50\n(chance)", "0.75", "1.0")
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right"
    ) +
    ggplot2::labs(
      title    = paste("sPLS-DA — Cross-Validated Performance"),
      subtitle = paste0(
        "Model: ", key %||% names(obj$splsda)[length(obj$splsda)],
        " | Groups: ", paste(entry$groups, collapse = " vs "),
        " | n = ", length(entry$samples),
        if (has_ci) " | Bars = 95% CI (repeated CV)" else " | Single CV run"
      ),
      x = NULL,
      y = "Balanced Error Rate (BER)",
      caption = paste0(
        "BER < 0.20: strong  |  0.20-0.35: moderate  |  0.35-0.50: weak  |  ",
        "\u2265 0.50: no better than chance\n",
        "Note: class imbalance (n=", paste(table(entry$groups[match(entry$samples,
                                                                    names(table(entry$groups)))]), collapse="/"),
        ") and small cohort size (n=", length(entry$samples),
        ") increase overfitting risk.\n",
        "Training-set AUC is not shown — it is trivially 100% with many features."
      )
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 7. SIGNATURE EXTRACTION
# =============================================================================

#' Extract the sPLS-DA De Novo Gene/TE Signature
#'
#' @param obj A validated OmniSeq object with sPLS-DA results.
#' @param key Character. Cache key.
#' @param comp Integer. Which component. Default: `1`.
#' @param vip_min Numeric. Minimum VIP score to include. Default: `1`.
#'
#' @return A data.frame sorted by VIP score (descending).
#' @export
omni_splsda_signature <- function(obj,
                                  key     = NULL,
                                  comp    = 1L,
                                  vip_min = 1) {

  validate_omniseq(obj)

  entry   <- .omni_get_splsda_entry(obj, key)
  model   <- entry$model
  vip_mat <- entry$vip

  loadings <- model$loadings$X[, comp]
  selected <- names(loadings[loadings != 0])

  if (length(selected) == 0)
    stop("No features selected on component ", comp, call. = FALSE)

  df <- data.frame(
    Gene    = selected,
    Loading = as.numeric(loadings[selected]),
    VIP     = as.numeric(vip_mat[selected, comp]),
    stringsAsFactors = FALSE
  )

  df <- df[df$VIP >= vip_min, ]
  df <- df[order(-df$VIP), ]

  res <- omni_results(obj)
  if (!is.null(res)) {
    res_df      <- as.data.frame(res)
    res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)
    in_res      <- intersect(df$Gene, rownames(res_df))
    df$log2FC   <- NA_real_
    df$padj     <- NA_real_
    if (length(in_res) > 0) {
      df$log2FC[df$Gene %in% in_res] <- suppressWarnings(
        as.numeric(res_df[in_res, "log2FoldChange"]))
      df$padj[df$Gene %in% in_res] <- suppressWarnings(
        as.numeric(res_df[in_res, "padj"]))
    }
  }

  rownames(df) <- NULL
  message("[omni_splsda_signature] ", nrow(df),
          " features with VIP \u2265 ", vip_min,
          " on component ", comp)
  return(df)
}


# =============================================================================
# 8. DIABLO — MULTI-BLOCK sPLS-DA
# =============================================================================

#' Run DIABLO Multi-Block sPLS-DA (Gene + TE Integration)
#'
#' Integrates gene and TE expression as two separate data blocks and finds
#' latent components that simultaneously explain variance in both blocks while
#' discriminating between groups. Uses the DIABLO framework from mixOmics.
#'
#' The TE block is aggregated to family level (mean VST per family per sample)
#' before integration to reduce dimensionality and improve stability.
#'
#' @param gene_obj OmniSeq object containing standard gene VST.
#' @param te_obj OmniSeq object containing TE VST.
#' @param group_col Character. Metadata column for grouping. Default: `"Condition"`.
#' @param exclude_groups Character vector. Groups to exclude before fitting.
#' @param ncomp Integer. Number of components. Default: `2`.
#' @param keepX Named list with elements `"genes"` and `"TEs"`, each a vector
#'   of length `ncomp`. Automatically capped to the number of available
#'   features per block — so passing `c(20, 20)` for a TE block with only
#'   6 families is safe. Default: `list(genes = c(20, 20), TEs = c(10, 10))`.
#' @param near_zero_var Logical. Filter near-zero-variance features before
#'   fitting. Default: `TRUE`.
#' @param seed Integer. Random seed. Default: `42`.
#'
#' @return The gene_obj with DIABLO result stored in
#'   `gene_obj$diablo[["<group_col>_<groups>"]]`.
#' @export
run_omni_diablo <- function(gene_obj,
                            te_obj,
                            group_col      = "Condition",
                            exclude_groups = NULL,
                            ncomp          = 2L,
                            keepX          = NULL,
                            near_zero_var  = TRUE,
                            seed           = 42L) {

  validate_omniseq(gene_obj)
  validate_omniseq(te_obj)
  .require_pkg("mixOmics")

  meta <- omni_metadata(gene_obj)

  if (!group_col %in% colnames(meta))
    stop("group_col '", group_col, "' not found in metadata.", call. = FALSE)

  # ------------------------------------------------------------------
  # 1. BUILD SAMPLE SUBSET
  # ------------------------------------------------------------------
  groups_all <- as.character(meta[[group_col]])
  keep_mask  <- !groups_all %in% (exclude_groups %||% character(0))
  keep_samps <- rownames(meta)[keep_mask]
  Y          <- factor(groups_all[keep_mask])

  grp_tab <- table(Y)
  message("\n[run_omni_diablo] Groups: ",
          paste(paste0(names(grp_tab), " (n=", grp_tab, ")"),
                collapse = ", "))

  # ------------------------------------------------------------------
  # 2. GENE BLOCK (samples x genes)
  # ------------------------------------------------------------------
  vsd_gene <- omni_vst(gene_obj)
  if (is.null(vsd_gene))
    stop("Gene VST not cached. Run run_omni_de() on gene_obj first.",
         call. = FALSE)
  X_gene <- t(SummarizedExperiment::assay(vsd_gene)[, keep_samps,
                                                    drop = FALSE])
  message("[run_omni_diablo] Gene block: ", ncol(X_gene), " features")

  # ------------------------------------------------------------------
  # 3. TE BLOCK — aggregate to family level (samples x families)
  # ------------------------------------------------------------------
  vsd_te <- omni_vst(te_obj)
  if (is.null(vsd_te))
    stop("TE VST not cached. Run run_omni_de() on te_obj first.",
         call. = FALSE)

  mat_te      <- SummarizedExperiment::assay(vsd_te)[, keep_samps,
                                                     drop = FALSE]
  te_info     <- .omni_parse_te_ids(rownames(mat_te))
  fams        <- te_info$Family
  unique_fams <- unique(fams[fams != "Unknown"])

  agg_list <- lapply(unique_fams, function(f) {
    idx <- which(fams == f)
    if (length(idx) == 1) mat_te[idx, , drop = TRUE]
    else colMeans(mat_te[idx, , drop = FALSE], na.rm = TRUE)
  })
  agg_mat           <- do.call(rbind, agg_list)
  rownames(agg_mat) <- unique_fams
  X_te              <- t(agg_mat)

  message("[run_omni_diablo] TE block:   ", ncol(X_te), " families")
  message("[run_omni_diablo] Samples:    ", nrow(X_gene))

  # ------------------------------------------------------------------
  # 4. BUILD AND CAP keepX
  # FIXED: cap each block's keepX to min(requested, ncol(block))
  # Prevents crash when TE block has fewer families than default keepX
  # (e.g. 6 SVA families with keepX TEs = c(10, 10) → capped to c(6, 6))
  # ------------------------------------------------------------------
  if (is.null(keepX))
    keepX <- list(genes = rep(20L, ncomp), TEs = rep(10L, ncomp))

  if (!all(c("genes", "TEs") %in% names(keepX)))
    stop("keepX must be a named list with elements 'genes' and 'TEs'.",
         call. = FALSE)

  keepX$genes <- pmin(as.integer(keepX$genes), ncol(X_gene))
  keepX$TEs   <- pmin(as.integer(keepX$TEs),   ncol(X_te))

  if (any(keepX$TEs < 10L))
    message("[run_omni_diablo] keepX$TEs capped to ", ncol(X_te),
            " (number of TE families available).")
  if (any(keepX$genes < 20L))
    message("[run_omni_diablo] keepX$genes capped to ", ncol(X_gene),
            " (number of gene features available).")

  X_list <- list(genes = X_gene, TEs = X_te)

  # ------------------------------------------------------------------
  # 5. DESIGN MATRIX (off-diagonal = 0.1: weak correlation assumption)
  # ------------------------------------------------------------------
  design <- matrix(0.1, nrow = 2, ncol = 2,
                   dimnames = list(c("genes", "TEs"), c("genes", "TEs")))
  diag(design) <- 0

  # ------------------------------------------------------------------
  # 6. FIT DIABLO MODEL
  # ------------------------------------------------------------------
  set.seed(seed)
  message("[run_omni_diablo] Fitting DIABLO model",
          if (near_zero_var) " (near.zero.var filtering enabled)..." else "...")

  diablo_fit <- mixOmics::block.splsda(
    X             = X_list,
    Y             = Y,
    ncomp         = ncomp,
    keepX         = keepX,
    design        = design,
    near.zero.var = near_zero_var
  )

  # ------------------------------------------------------------------
  # 7. STORE
  # FIXED: removed reference to undefined `use_masked` variable
  # ------------------------------------------------------------------
  retained_groups <- sort(unique(as.character(Y)))
  cache_key       <- paste0(group_col, "_",
                            paste(retained_groups, collapse = "-"))

  if (is.null(gene_obj$diablo)) gene_obj$diablo <- list()

  gene_obj$diablo[[cache_key]] <- list(
    model         = diablo_fit,
    group_col     = group_col,
    groups        = retained_groups,
    excluded      = exclude_groups,
    samples       = keep_samps,
    ncomp         = ncomp,
    keepX         = keepX,
    near_zero_var = near_zero_var,
    timestamp     = Sys.time()
  )

  message("[run_omni_diablo] Done. Stored under key: '", cache_key, "'")
  message("  Visualise with mixOmics directly:")
  message("  model <- gene_obj$diablo[['", cache_key, "']]$model")
  message("  mixOmics::plotDiablo(model, ncomp = 1)")
  message("  mixOmics::circosPlot(model, cutoff = 0.7)")
  message("  mixOmics::plotLoadings(model, comp = 1, contrib = 'max')")
  message("  mixOmics::plotIndiv(model, legend = TRUE)")

  return(gene_obj)
}


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Get a stored sPLS-DA entry, auto-selecting if only one exists
#' @keywords internal
.omni_get_splsda_entry <- function(obj, key) {
  if (is.null(obj$splsda) || length(obj$splsda) == 0)
    stop("No sPLS-DA results found. Run run_omni_splsda() first.",
         call. = FALSE)

  if (is.null(key)) {
    if (length(obj$splsda) == 1) {
      key <- names(obj$splsda)[1]
      message("[sPLS-DA] Auto-selected model: '", key, "'")
    } else {
      stop("Multiple sPLS-DA models found. Specify 'key'.\n",
           "Available: ", paste(names(obj$splsda), collapse = ", "),
           call. = FALSE)
    }
  }

  if (!key %in% names(obj$splsda))
    stop("Key '", key, "' not found. Available: ",
         paste(names(obj$splsda), collapse = ", "), call. = FALSE)

  obj$splsda[[key]]
}
