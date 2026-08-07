# =============================================================================
# OmniSeq — 11_Causal.R
#
# Causal and mediation analysis for interferon-TE-mutation pathway testing.
#
# SCIENTIFIC RATIONALE:
#   The causal hypothesis is: mutation type (TLR7 vs UNC93B1) drives aberrant
#   TE activation, which in turn triggers innate immune sensing (cGAS-STING,
#   RIG-I/MDA5, or TLR7), which produces the interferon storm (IRG-S).
#
#   Mediation analysis formally tests:
#     mutation_type → TE_sensor_score → IRG-S_score
#
#   Decomposing the total effect into:
#     - Average Causal Mediation Effect (ACME): the indirect path through TEs
#     - Average Direct Effect (ADE): mutation effect NOT through TEs
#     - Proportion mediated: how much of the mutation effect is TE-mediated
#
#   If the ACME is significant, TE activation is a statistically supported
#   mediator between mutation type and the interferon response.
#
# DESIGN:
#   - Requires gene_obj (for IRG-S score) and sensor scores from
#     run_omni_sensor_score() or a pre-computed sensor score data frame
#   - Mutation type encoded as binary (e.g. TLR7=1, UNC93B1=0)
#   - SPEN patients excluded from mediation (n=2, insufficient)
#   - Bootstrapped confidence intervals (default 1000 iterations)
#   - Reverse causality model tested as sensitivity analysis
#   - Results stored in obj$mediation[["<mediator_name>"]]
#
# FUNCTIONS:
#   run_omni_mediation()        Single mediator (one sensor)
#   run_omni_mediation_multi()  All three sensors automatically
#   plot_omni_mediation()       Forest plot of ACME, ADE, total effect
#   plot_omni_mediation_path()  Path diagram with effect sizes on arrows
# =============================================================================


# =============================================================================
# 1. SINGLE MEDIATOR MEDIATION ANALYSIS
# =============================================================================

#' Run Mediation Analysis for a Single TE Sensor Score
#'
#' Tests whether a TE innate sensor score mediates the relationship between
#' mutation type and IRG-S score. Uses the `mediation` R package with
#' bootstrapped confidence intervals.
#'
#' The causal model tested is:
#' ```
#' Mutation type → TE sensor score → IRG-S score
#' ```
#'
#' Both a forward model (TE mediates mutation → IFN) and a reverse model
#' (IFN mediates mutation → TE) are fitted as a sensitivity analysis.
#'
#' @param gene_obj A validated OmniSeq gene object. Used to compute IRG-S
#'   score from the VST matrix.
#' @param sensor_scores A data frame from `run_omni_sensor_score()` containing
#'   per-sample sensor scores and metadata. Must contain columns for the
#'   mediator, `IRG_S_Score`, and a mutation type variable.
#' @param mediator Character. Column name in `sensor_scores` to use as the
#'   mediating variable. One of `"cGAS_STING"`, `"RIG_I_MDA5"`, `"TLR7_TLR8"`,
#'   or any other numeric column.
#' @param treatment_col Character. Column name encoding mutation type as a
#'   binary variable (0/1). Default: `"Genotype"`. The function will attempt
#'   to encode it as binary automatically.
#' @param treatment_vals Named character vector of length 2 mapping column
#'   values to 0 and 1. E.g. `c("0" = "UNC93B1", "1" = "TLR7")`.
#'   If `NULL`, the first two unique values (alphabetically) are used.
#' @param exclude_groups Character vector. Values in `treatment_col` to
#'   exclude. Default: `c("SPEN", "Ctrl", "Control")`.
#' @param outcome Character. Column name for the outcome variable.
#'   Default: `"IRG_S_Score"`.
#' @param covariates Character vector. Additional covariates to include in
#'   both mediator and outcome models. Default: `NULL`.
#' @param boot Logical. Use bootstrap for confidence intervals. Default:
#'   `TRUE`.
#' @param boot_sims Integer. Number of bootstrap simulations. Default: `1000`.
#' @param seed Integer. Random seed. Default: `42`.
#'
#' @return The gene_obj with mediation results stored in
#'   `gene_obj$mediation[["<mediator>"]]`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Run sensor scoring first
#' gene_obj <- run_omni_sensor_score(gene_obj, te_obj = sine_obj)
#' scores   <- gene_obj$sensor_scores$scores
#'
#' # Test whether cGAS-STING mediates TLR7 vs UNC93B1 → IRG-S
#' gene_obj <- run_omni_mediation(
#'   gene_obj,
#'   sensor_scores  = scores,
#'   mediator       = "cGAS_STING",
#'   treatment_col  = "Genotype",
#'   treatment_vals = c("0" = "UNC93B1", "1" = "TLR7"),
#'   exclude_groups = c("SPEN", "Ctrl")
#' )
#'
#' plot_omni_mediation(gene_obj)
#' }
run_omni_mediation <- function(gene_obj,
                               sensor_scores,
                               mediator        = "cGAS_STING",
                               treatment_col   = "Genotype",
                               treatment_vals  = NULL,
                               exclude_groups  = c("SPEN", "Ctrl", "Control"),
                               outcome         = "IRG_S_Score",
                               covariates      = NULL,
                               cell_state_scores = NULL,
                               boot            = TRUE,
                               boot_sims       = 1000L,
                               seed            = 42L) {

  validate_omniseq(gene_obj)
  .require_pkg("mediation")

  # ------------------------------------------------------------------
  # 1. VALIDATE INPUTS
  # ------------------------------------------------------------------
  if (!mediator %in% colnames(sensor_scores))
    stop("Mediator '", mediator, "' not found in sensor_scores.\n",
         "Available columns: ",
         paste(colnames(sensor_scores), collapse = ", "), call. = FALSE)

  if (!outcome %in% colnames(sensor_scores))
    stop("Outcome '", outcome, "' not found in sensor_scores.\n",
         "Available columns: ",
         paste(colnames(sensor_scores), collapse = ", "), call. = FALSE)

  if (!treatment_col %in% colnames(sensor_scores))
    stop("treatment_col '", treatment_col,
         "' not found in sensor_scores.", call. = FALSE)

  # ------------------------------------------------------------------
  # 2. PREPARE DATA
  # ------------------------------------------------------------------
  df <- sensor_scores

  # If scRNA-derived cell state scores are provided (e.g. IFN-like monocyte
  # burden from MuSiC/SAVI deconvolution), merge them into df so they can
  # be used as mediators — architecture hook for future SAVI integration
  if (!is.null(cell_state_scores)) {
    if (!is.data.frame(cell_state_scores))
      stop("cell_state_scores must be a data frame with a Sample column.",
           call. = FALSE)
    common_samps <- intersect(df$Sample, cell_state_scores$Sample)
    if (length(common_samps) == 0)
      stop("No common samples between sensor_scores and cell_state_scores.",
           call. = FALSE)
    # Merge new columns — do not overwrite existing sensor columns
    new_cols <- setdiff(colnames(cell_state_scores), c("Sample", colnames(df)))
    if (length(new_cols) > 0) {
      df <- merge(df,
                  cell_state_scores[, c("Sample", new_cols), drop = FALSE],
                  by = "Sample", all.x = TRUE)
      message("[run_omni_mediation] Merged ", length(new_cols),
              " cell state score column(s): ",
              paste(new_cols, collapse = ", "))
    }
  }

  # Exclude unwanted groups
  if (!is.null(exclude_groups)) {
    keep <- !as.character(df[[treatment_col]]) %in% exclude_groups
    n_excluded <- sum(!keep)
    df <- df[keep, ]
    if (n_excluded > 0)
      message("[run_omni_mediation] Excluded ", n_excluded,
              " samples from groups: ",
              paste(exclude_groups, collapse = ", "))
  }

  # Encode treatment as binary
  trt_vals <- as.character(df[[treatment_col]])
  unique_vals <- sort(unique(trt_vals))

  if (length(unique_vals) < 2)
    stop("At least 2 unique values required in '", treatment_col,
         "' after exclusion. Found: ", paste(unique_vals, collapse = ", "),
         call. = FALSE)

  if (length(unique_vals) > 2) {
    warning("[run_omni_mediation] More than 2 unique values in '",
            treatment_col, "': ", paste(unique_vals, collapse = ", "),
            ". Using first two: '", unique_vals[1], "' (0) and '",
            unique_vals[2], "' (1). Specify treatment_vals to override.",
            call. = FALSE)
    unique_vals <- unique_vals[1:2]
  }

  if (!is.null(treatment_vals)) {
    if (!all(c("0", "1") %in% names(treatment_vals)))
      stop("treatment_vals must be a named vector with names '0' and '1'.",
           call. = FALSE)
    val0 <- treatment_vals[["0"]]
    val1 <- treatment_vals[["1"]]
  } else {
    val0 <- unique_vals[1]
    val1 <- unique_vals[2]
  }

  df$Treatment <- ifelse(trt_vals == val1, 1L, 0L)
  df           <- df[trt_vals %in% c(val0, val1), ]

  n0 <- sum(df$Treatment == 0)
  n1 <- sum(df$Treatment == 1)

  message("\n[run_omni_mediation] Mediation analysis")
  message("  Mediator  : ", mediator)
  message("  Treatment : ", treatment_col,
          " (", val0, " = 0, n=", n0, " | ",
          val1, " = 1, n=", n1, ")")
  message("  Outcome   : ", outcome)
  message("  Samples   : ", nrow(df))

  if (nrow(df) < 6)
    stop("Fewer than 6 samples after filtering. Cannot run mediation.",
         call. = FALSE)

  if (any(is.na(df[[mediator]])))
    stop("Mediator '", mediator,
         "' contains NA values. Check run_omni_sensor_score() output.",
         call. = FALSE)

  # ------------------------------------------------------------------
  # 3. BUILD MODEL FORMULAS
  # ------------------------------------------------------------------
  # Mediator model: M ~ Treatment + covariates
  # Outcome model:  Y ~ Treatment + M + covariates

  covar_str <- if (!is.null(covariates) && length(covariates) > 0)
    paste("+", paste(covariates, collapse = " + "))
  else ""

  med_formula <- stats::as.formula(
    paste(mediator, "~ Treatment", covar_str)
  )
  out_formula <- stats::as.formula(
    paste(outcome, "~ Treatment +", mediator, covar_str)
  )

  message("[run_omni_mediation] Mediator model: ",
          deparse(med_formula))
  message("[run_omni_mediation] Outcome model:  ",
          deparse(out_formula))

  # ------------------------------------------------------------------
  # 4. FIT MODELS
  # ------------------------------------------------------------------
  # mediation::mediate() re-evaluates model$call internally in a clean
  # environment, so formula objects must be embedded via the model environment.
  # Solution: set the model environment to the current local environment so
  # all local variables (including the formula) remain accessible.
  local_env <- environment()

  med_model <- stats::lm(med_formula, data = df)
  out_model <- stats::lm(out_formula, data = df)

  # Patch model environments so mediation::mediate() can re-evaluate calls
  environment(med_model$terms) <- local_env
  environment(out_model$terms) <- local_env
  med_model$call$formula <- med_formula
  out_model$call$formula <- out_formula
  med_model$call$data    <- quote(df)
  out_model$call$data    <- quote(df)

  # ------------------------------------------------------------------
  # 5. RUN MEDIATION
  # ------------------------------------------------------------------
  set.seed(seed)
  message("[run_omni_mediation] Running mediation analysis (",
          if (boot) paste0(boot_sims, " bootstrap simulations")
          else "no bootstrap", ")...")

  med_result <- local({
    mediation::mediate(
      model.m      = med_model,
      model.y      = out_model,
      treat        = "Treatment",
      mediator     = mediator,
      boot         = boot,
      sims         = boot_sims,
      boot.ci.type = "perc"
    )
  })

  # ------------------------------------------------------------------
  # 6. REVERSE CAUSALITY (SENSITIVITY)
  # ------------------------------------------------------------------
  message("[run_omni_mediation] Testing reverse causality (sensitivity)...")
  set.seed(seed)

  rev_med_formula <- stats::as.formula(
    paste(outcome, "~ Treatment", covar_str)
  )
  rev_out_formula <- stats::as.formula(
    paste(mediator, "~ Treatment +", outcome, covar_str)
  )

  rev_med_model <- stats::lm(rev_med_formula, data = df)
  rev_out_model <- stats::lm(rev_out_formula, data = df)
  rev_med_model$call$formula <- rev_med_formula
  rev_out_model$call$formula <- rev_out_formula
  rev_med_model$call$data    <- quote(df)
  rev_out_model$call$data    <- quote(df)

  rev_result <- tryCatch({
    local({
      mediation::mediate(
        model.m      = rev_med_model,
        model.y      = rev_out_model,
        treat        = "Treatment",
        mediator     = outcome,
        boot         = boot,
        sims         = min(boot_sims, 500L),
        boot.ci.type = "perc"
      )
    })
  }, error = function(e) {
    message("[run_omni_mediation] Reverse causality test failed: ",
            e$message)
    NULL
  })

  # ------------------------------------------------------------------
  # 7. EXTRACT SUMMARY
  # ------------------------------------------------------------------
  med_sum <- summary(med_result)

  # Guard: proportion mediated (n.avg) can be outside [0,1] when effects
  # are small or go in opposing directions — flag clearly rather than silently
  # storing a nonsensical value
  prop_med_raw <- med_result$n.avg
  if (!is.na(prop_med_raw) && (prop_med_raw < 0 || prop_med_raw > 1))
    message("[run_omni_mediation] Warning: Proportion mediated = ",
            round(prop_med_raw, 3),
            " (outside [0,1]). This occurs when ACME and total effect have ",
            "opposite signs or total effect is near zero. Interpret with caution.")

  result_summary <- data.frame(
    Effect            = c("ACME (indirect)", "ADE (direct)",
                          "Total effect", "Proportion mediated"),
    Estimate          = c(med_result$d.avg,  med_result$z.avg,
                          med_result$tau.coef, prop_med_raw),
    CI_Lower          = c(med_result$d.avg.ci[1], med_result$z.avg.ci[1],
                          med_result$tau.ci[1],    med_result$n.avg.ci[1]),
    CI_Upper          = c(med_result$d.avg.ci[2], med_result$z.avg.ci[2],
                          med_result$tau.ci[2],    med_result$n.avg.ci[2]),
    P_value           = c(med_result$d.avg.p,  med_result$z.avg.p,
                          med_result$tau.p,     med_result$n.avg.p),
    stringsAsFactors  = FALSE
  )

  message("\n[run_omni_mediation] Results:")
  for (i in seq_len(nrow(result_summary))) {
    message(sprintf("  %-22s: %.3f [%.3f, %.3f]  p = %.3f",
                    result_summary$Effect[i],
                    result_summary$Estimate[i],
                    result_summary$CI_Lower[i],
                    result_summary$CI_Upper[i],
                    result_summary$P_value[i]))
  }

  # ------------------------------------------------------------------
  # 8. STORE
  # ------------------------------------------------------------------
  if (is.null(gene_obj$mediation)) gene_obj$mediation <- list()

  gene_obj$mediation[[mediator]] <- list(
    result          = med_result,
    reverse_result  = rev_result,
    summary         = result_summary,
    mediator        = mediator,
    treatment_col   = treatment_col,
    treatment_vals  = c("0" = val0, "1" = val1),
    outcome         = outcome,
    covariates      = covariates,
    data            = df,
    boot            = boot,
    boot_sims       = boot_sims,
    n_samples       = nrow(df),
    timestamp       = Sys.time()
  )

  message("\n[run_omni_mediation] Stored under: gene_obj$mediation[['",
          mediator, "']]")

  # Summarise reverse causality result alongside forward result
  if (!is.null(rev_result)) {
    message("[run_omni_mediation] Reverse causality check (IFN → TE):")
    message(sprintf("  Reverse ACME: %.3f  p = %.3f %s",
                    rev_result$d.avg,
                    rev_result$d.avg.p,
                    if (rev_result$d.avg.p < 0.05)
                      "*** significant — reverse path cannot be excluded"
                    else "(ns — forward direction supported)"))
  }

  message("  Next: plot_omni_mediation(gene_obj)")
  message("        plot_omni_mediation_path(gene_obj)")

  return(gene_obj)
}


# =============================================================================
# 2. MULTI-SENSOR MEDIATION
# =============================================================================

#' Run Mediation Analysis for All Three TE Sensor Scores
#'
#' Convenience wrapper that calls `run_omni_mediation()` once per innate
#' sensor (cGAS-STING, RIG-I/MDA5, TLR7/TLR8). Results for all three
#' mediators are stored in `gene_obj$mediation` under their respective names.
#'
#' @param gene_obj A validated OmniSeq gene object.
#' @param sensor_scores A data frame from `run_omni_sensor_score()`.
#' @param treatment_col Character. Column encoding mutation type. Default:
#'   `"Genotype"`.
#' @param treatment_vals Named character vector. See `run_omni_mediation()`.
#' @param exclude_groups Character vector. Groups to exclude. Default:
#'   `c("SPEN", "Ctrl", "Control")`.
#' @param outcome Character. Outcome variable. Default: `"IRG_S_Score"`.
#' @param covariates Character vector. Additional covariates.
#' @param boot Logical. Bootstrap CIs. Default: `TRUE`.
#' @param boot_sims Integer. Bootstrap simulations per mediator. Default:
#'   `1000`.
#' @param seed Integer. Random seed. Default: `42`.
#'
#' @return The gene_obj with all three mediation results stored.
#' @export
#'
#' @examples
#' \dontrun{
#' gene_obj <- run_omni_mediation_multi(
#'   gene_obj,
#'   sensor_scores  = gene_obj$sensor_scores$scores,
#'   treatment_col  = "Genotype",
#'   treatment_vals = c("0" = "UNC93B1", "1" = "TLR7"),
#'   exclude_groups = c("SPEN", "Ctrl")
#' )
#'
#' # Plot all three mediators
#' plot_omni_mediation(gene_obj)
#' }
run_omni_mediation_multi <- function(gene_obj,
                                     sensor_scores,
                                     treatment_col   = "Genotype",
                                     treatment_vals  = NULL,
                                     exclude_groups  = c("SPEN", "Ctrl",
                                                         "Control"),
                                     outcome         = "IRG_S_Score",
                                     covariates      = NULL,
                                     boot            = TRUE,
                                     boot_sims       = 1000L,
                                     seed            = 42L) {

  validate_omniseq(gene_obj)

  # Determine which sensors are available in the score data frame
  available_sensors <- intersect(
    c("cGAS_STING", "RIG_I_MDA5", "TLR7_TLR8"),
    colnames(sensor_scores)
  )

  if (length(available_sensors) == 0)
    stop("No sensor score columns found in sensor_scores.\n",
         "Expected: cGAS_STING, RIG_I_MDA5, TLR7_TLR8\n",
         "Available: ", paste(colnames(sensor_scores), collapse = ", "),
         call. = FALSE)

  # Check for NA sensors (score was NA because too few loci matched)
  na_sensors <- available_sensors[
    sapply(available_sensors, function(s) all(is.na(sensor_scores[[s]])))
  ]
  if (length(na_sensors) > 0) {
    message("[run_omni_mediation_multi] Skipping sensors with all NA values: ",
            paste(na_sensors, collapse = ", "))
    available_sensors <- setdiff(available_sensors, na_sensors)
  }

  message("\n[run_omni_mediation_multi] Running mediation for ",
          length(available_sensors), " sensor(s): ",
          paste(available_sensors, collapse = ", "))

  for (sens in available_sensors) {
    message("\n", strrep("-", 50))
    gene_obj <- run_omni_mediation(
      gene_obj       = gene_obj,
      sensor_scores  = sensor_scores,
      mediator       = sens,
      treatment_col  = treatment_col,
      treatment_vals = treatment_vals,
      exclude_groups = exclude_groups,
      outcome        = outcome,
      covariates     = covariates,
      boot           = boot,
      boot_sims      = boot_sims,
      seed           = seed
    )
  }

  message("\n", strrep("=", 50))
  message("[run_omni_mediation_multi] All mediation analyses complete.")
  message("  Stored mediators: ",
          paste(names(gene_obj$mediation), collapse = ", "))
  message("  Next: plot_omni_mediation(gene_obj)")

  return(gene_obj)
}


# =============================================================================
# 3. FOREST PLOT
# =============================================================================

#' Plot Mediation Effects as a Forest Plot
#'
#' Displays ACME (indirect effect), ADE (direct effect), total effect, and
#' proportion mediated for one or all mediators. Points show estimates,
#' horizontal lines show 95% confidence intervals. A vertical dashed line
#' at 0 indicates the null hypothesis.
#'
#' @param obj A validated OmniSeq object with mediation results.
#' @param mediators Character vector. Which mediators to plot. Default: all
#'   stored mediators.
#' @param effects Character vector. Which effects to show. Default:
#'   `c("ACME (indirect)", "ADE (direct)", "Total effect")`.
#'   `"Proportion mediated"` can be added but is on a different scale.
#' @param title Optional custom title.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_mediation <- function(obj,
                                mediators = NULL,
                                effects   = c("ACME (indirect)",
                                              "ADE (direct)",
                                              "Total effect"),
                                title     = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  if (is.null(obj$mediation) || length(obj$mediation) == 0)
    stop("No mediation results. Run run_omni_mediation() first.",
         call. = FALSE)

  if (is.null(mediators))
    mediators <- names(obj$mediation)

  bad_med <- setdiff(mediators, names(obj$mediation))
  if (length(bad_med) > 0)
    stop("Mediator(s) not found: ", paste(bad_med, collapse = ", "),
         "\nAvailable: ", paste(names(obj$mediation), collapse = ", "),
         call. = FALSE)

  # Friendly sensor names
  sensor_labels <- c(
    "cGAS_STING" = "cGAS-STING\n(cytoplasmic DNA)",
    "RIG_I_MDA5" = "RIG-I / MDA5\n(cytoplasmic dsRNA)",
    "TLR7_TLR8"  = "TLR7 / TLR8\n(endosomal ssRNA)"
  )

  # Combine summaries across mediators
  plot_df <- do.call(rbind, lapply(mediators, function(med) {
    entry <- obj$mediation[[med]]
    df    <- entry$summary
    df    <- df[df$Effect %in% effects, ]
    df$Mediator <- sensor_labels[med] %||% med
    df
  }))

  plot_df$Mediator <- factor(plot_df$Mediator,
                             levels = unique(plot_df$Mediator))
  plot_df$Effect   <- factor(plot_df$Effect, levels = effects)
  plot_df$Sig      <- plot_df$P_value < 0.05

  if (is.null(title))
    title <- paste0("Mediation Analysis — ",
                    unique(sapply(mediators, function(m)
                      obj$mediation[[m]]$treatment_col)),
                    " \u2192 TE Sensor \u2192 IRG-S")

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = Estimate, y = Effect,
                                    color = Sig,
                                    xmin  = CI_Lower,
                                    xmax  = CI_Upper)) +
    ggplot2::geom_pointrange(size = 0.8, linewidth = 1) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey40", linewidth = 0.5) +
    ggplot2::facet_wrap(~ Mediator, ncol = length(mediators),
                        scales = "free_x") +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "#E74C3C", "FALSE" = "grey50"),
      labels = c("TRUE" = "p < 0.05", "FALSE" = "p \u2265 0.05"),
      name   = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(face = "bold", size = 10),
      axis.text.y      = ggplot2::element_text(size = 10),
      legend.position  = "bottom"
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = paste("Bootstrapped 95% CI |",
                       unique(sapply(mediators, function(m) {
                         e <- obj$mediation[[m]]
                         paste0(e$treatment_vals["0"], " vs ",
                                e$treatment_vals["1"])
                       })), "| n =",
                       unique(sapply(mediators, function(m)
                         obj$mediation[[m]]$n_samples))),
      x        = "Effect estimate (z-score units)",
      y        = NULL
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 4. PATH DIAGRAM
# =============================================================================

#' Plot Mediation Path Diagram
#'
#' Produces a simple path diagram showing the mediation structure:
#' Treatment → Mediator → Outcome, with effect sizes on the arrows.
#' Solid arrows indicate significant effects (p < 0.05), dashed arrows
#' indicate non-significant effects.
#'
#' @param obj A validated OmniSeq object with mediation results.
#' @param mediator Character. Which mediator to plot. If `NULL` and only
#'   one mediator exists, uses it automatically.
#' @param title Optional custom title.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_mediation_path <- function(obj,
                                     mediator = NULL,
                                     title    = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  if (is.null(obj$mediation) || length(obj$mediation) == 0)
    stop("No mediation results. Run run_omni_mediation() first.",
         call. = FALSE)

  # Auto-select mediator
  if (is.null(mediator)) {
    if (length(obj$mediation) == 1) {
      mediator <- names(obj$mediation)[1]
      message("[plot_omni_mediation_path] Auto-selected mediator: '",
              mediator, "'")
    } else {
      stop("Multiple mediators available. Specify 'mediator'.\n",
           "Available: ", paste(names(obj$mediation), collapse = ", "),
           call. = FALSE)
    }
  }

  if (!mediator %in% names(obj$mediation))
    stop("Mediator '", mediator, "' not found.", call. = FALSE)

  entry <- obj$mediation[[mediator]]
  summ  <- entry$summary

  # Friendly labels
  sensor_labels <- c(
    "cGAS_STING" = "cGAS-STING\n(TE sensor)",
    "RIG_I_MDA5" = "RIG-I / MDA5\n(TE sensor)",
    "TLR7_TLR8"  = "TLR7 / TLR8\n(TE sensor)"
  )
  med_label <- sensor_labels[mediator] %||% mediator
  trt_label <- paste0(entry$treatment_vals["1"], " vs ",
                      entry$treatment_vals["0"])
  out_label <- "IRG-S Score\n(IFN response)"

  # Extract effects
  acme_row <- summ[summ$Effect == "ACME (indirect)", ]
  ade_row  <- summ[summ$Effect == "ADE (direct)", ]
  tot_row  <- summ[summ$Effect == "Total effect", ]
  prop_row <- summ[summ$Effect == "Proportion mediated", ]

  # Node positions
  nodes <- data.frame(
    label = c(trt_label, med_label, out_label),
    x     = c(0, 1, 2),
    y     = c(0, 1, 0),
    stringsAsFactors = FALSE
  )

  # Arrows: a→b path (treatment to mediator), b→c path (mediator to outcome),
  # a→c path (direct, treatment to outcome)
  get_effect <- function(row) {
    if (nrow(row) == 0) return(list(est = NA, p = NA, sig = FALSE))
    list(
      est = round(row$Estimate, 3),
      p   = round(row$P_value, 3),
      sig = row$P_value < 0.05
    )
  }

  # We extract path-specific coefficients from the stored lm models
  # rather than the mediation summary which gives total ACME/ADE
  med_model_coef <- tryCatch(
    stats::coef(entry$result$model.m),
    error = function(e) NULL
  )
  out_model_coef <- tryCatch(
    stats::coef(entry$result$model.y),
    error = function(e) NULL
  )

  a_path <- if (!is.null(med_model_coef) && "Treatment" %in% names(med_model_coef))
    round(med_model_coef["Treatment"], 3) else NA
  b_path <- if (!is.null(out_model_coef) && mediator %in% names(out_model_coef))
    round(out_model_coef[mediator], 3) else NA
  c_path <- get_effect(tot_row)
  cp_path <- get_effect(ade_row)

  # ACME = a * b (approximately)
  acme_eff <- get_effect(acme_row)
  prop_med  <- if (nrow(prop_row) > 0) round(prop_row$Estimate * 100, 1) else NA

  if (is.null(title))
    title <- paste0("Mediation Path Diagram — ", mediator)

  # Build path diagram using ggplot2 annotations
  p <- ggplot2::ggplot() +
    ggplot2::xlim(-0.3, 2.3) +
    ggplot2::ylim(-0.5, 1.5) +
    ggplot2::coord_fixed() +

    # Nodes
    ggplot2::annotate("rect",
                      xmin = nodes$x - 0.35, xmax = nodes$x + 0.35,
                      ymin = nodes$y - 0.25, ymax = nodes$y + 0.25,
                      fill = c("#2980B9", "#E74C3C", "#27AE60"),
                      alpha = 0.85, color = "white", linewidth = 1.5) +
    ggplot2::annotate("text",
                      x = nodes$x, y = nodes$y,
                      label = nodes$label,
                      size = 3.5, color = "white", fontface = "bold",
                      lineheight = 0.9) +

    # Arrow: Treatment → Mediator (a path)
    ggplot2::annotate("segment",
                      x = 0.35, xend = 0.65,
                      y = 0.15, yend = 0.85,
                      arrow = grid::arrow(length = grid::unit(3, "mm"),
                                          type = "closed"),
                      color = "black", linewidth = 0.8) +
    ggplot2::annotate("text",
                      x = 0.42, y = 0.55,
                      label = paste0("a = ", a_path),
                      size = 3.5, color = "black",
                      angle = 55, fontface = "italic") +

    # Arrow: Mediator → Outcome (b path)
    ggplot2::annotate("segment",
                      x = 1.35, xend = 1.65,
                      y = 0.85, yend = 0.15,
                      arrow = grid::arrow(length = grid::unit(3, "mm"),
                                          type = "closed"),
                      color = "black", linewidth = 0.8) +
    ggplot2::annotate("text",
                      x = 1.58, y = 0.55,
                      label = paste0("b = ", b_path),
                      size = 3.5, color = "black",
                      angle = -55, fontface = "italic") +

    # Arrow: Treatment → Outcome (c' direct path)
    ggplot2::annotate("segment",
                      x = 0.35, xend = 1.65,
                      y = -0.05, yend = -0.05,
                      arrow = grid::arrow(length = grid::unit(3, "mm"),
                                          type = "closed"),
                      color   = if (isTRUE(cp_path$sig)) "#E74C3C" else "grey50",
                      linetype = if (isTRUE(cp_path$sig)) "solid" else "dashed",
                      linewidth = 0.8) +
    ggplot2::annotate("text",
                      x = 1, y = -0.18,
                      label = paste0("c' = ", cp_path$est,
                                     " (p=", cp_path$p, ")"),
                      size = 3.2, color = "black", fontface = "italic") +

    # ACME annotation
    ggplot2::annotate("text",
                      x = 1, y = 1.35,
                      label = paste0("ACME = ", acme_eff$est,
                                     " (p=", acme_eff$p, ")",
                                     if (!is.na(prop_med))
                                       paste0("\n", prop_med, "% mediated")
                                     else ""),
                      size = 3.5, color = "#8E44AD", fontface = "bold",
                      lineheight = 0.9) +

    ggplot2::theme_void(base_size = 13) +
    ggplot2::labs(
      title    = title,
      subtitle = paste0(
        "ACME: p=", acme_eff$p,
        if (isTRUE(acme_eff$sig)) " \u2605 significant" else " (ns)",
        " | Direct: p=", cp_path$p,
        if (isTRUE(cp_path$sig)) " \u2605 significant" else " (ns)"
      ),
      caption = paste0(
        "Red arrow = significant (p<0.05) | Grey dashed = non-significant\n",
        "a: treatment \u2192 mediator | b: mediator \u2192 outcome | ",
        "c': direct effect"
      )
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11, color = "grey40"),
      plot.caption  = ggplot2::element_text(size = 9, color = "grey50")
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# PLANNED: run_omni_mediation_diablo()
# =============================================================================
# When DIABLO is run with gene + TE blocks, Component 1 latent scores
# represent the multi-omics axis separating patients from controls.
# A planned extension will test:
#   mutation_type -> DIABLO_Component1_score -> IRG-S_score
# using the same mediation framework above.
#
# Architecture:
#   diablo_scores <- as.data.frame(
#     gene_obj$diablo[[key]]$model$variates$genes[, 1, drop = FALSE]
#   )
#   diablo_scores$Sample <- rownames(diablo_scores)
#   colnames(diablo_scores)[1] <- "DIABLO_Comp1"
#
#   gene_obj <- run_omni_mediation(
#     gene_obj,
#     sensor_scores  = merge(sensor_scores, diablo_scores, by = "Sample"),
#     mediator       = "DIABLO_Comp1",
#     treatment_col  = "Genotype",
#     treatment_vals = c("0" = "UNC93B1", "1" = "TLR7")
#   )
#
# This will be implemented as run_omni_mediation_diablo() once the
# combined gene + TE DIABLO model is finalized.
# =============================================================================
