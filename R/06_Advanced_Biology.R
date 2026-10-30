# =============================================================================
# OmniSeq — 06_Advanced_Biology.R  (Part 4B)
#
# Advanced multi-omics and systems biology functions:
#
# TRANSCRIPTION FACTOR ENRICHMENT:
#   run_omni_tf()              TF activity inference via decoupleR
#                              Supports DoRothEA and CollecTRI at runtime
#                              Regulons cached in obj$params$tf_regulons
#   plot_omni_tf_activity()    TF activity barplot (top activated/repressed)
#   plot_omni_tf_network()     TF regulatory network (static + interactive)
#                              Reads cached regulons — no re-download
#
# CELL TYPE DECONVOLUTION:
#   run_omni_deconvolution()   Flexible: immunedeconv (bulk) or MuSiC (scRNA)
#   plot_omni_deconv_heatmap() Cell type fractions heatmap per sample
#   plot_omni_deconv_corr()    Cell fractions vs DEGs or pathway scores
#                              common_samples bug fixed
#
# STRING NETWORK:
#   run_omni_string()          Query STRING for top DEG interactions
#                              Persistent local cache directory supported
#   plot_omni_string_network() Static (ggraph) + interactive (visNetwork)
#
# Fixes vs previous version:
#   - Regulons fetched ONCE and reused from obj$params$tf_regulons for all
#     subsequent contrasts. Re-download only if database/organism changes.
#   - Stat fingerprint (median, sum, range) printed per contrast to prove
#     each contrast produces a genuinely different input to ULM.
#   - OmnipathR deprecation warning suppressed (internal, uncontrollable)
#   - common_samples defined before correlation loop in plot_omni_deconv_corr()
#   - n_total unused variable removed from plot_omni_tf_activity()
#   - STRING input_directory supports persistent local cache path
#   - New object slots (tf_results, deconv, string_network) documented
# =============================================================================


# =============================================================================
# 1. TF ENRICHMENT — run_omni_tf()
# =============================================================================

#' Infer Transcription Factor Activity
#'
#' Uses `decoupleR` to infer TF activity using one of two modes:
#'
#' - `"contrast"` (default, recommended): Uses Wald t-statistics from the
#'   active contrast as input to ULM. Produces one score per TF per contrast —
#'   genuinely different across contrasts. Requires `stat` column in
#'   `results_cache` (populated by `run_omni_de()` / `run_omni_contrast()`).
#'
#' - `"samples"`: Uses the VST expression matrix restricted to a specified
#'   sample group (e.g. patients only). Produces per-sample scores.
#'   Never mixes patients and controls — the previous behaviour of mixing
#'   all 24 samples was incorrect and is no longer the default.
#'
#' Results are stored in:
#' - `obj$tf_contrast_cache[[contrast_key]]` — one score per TF (contrast mode)
#' - `obj$tf_sample_activity` — per-sample scores (samples mode)
#' - `obj$tf_results` — legacy slot (backward compat, always updated)
#'
#' @param obj A validated OmniSeq object with a fitted model.
#' @param mode Character. `"contrast"` (t-stat based, default) or `"samples"`
#'   (VST-based, restricted to `sample_group`).
#' @param sample_group Character. For `mode = "samples"`: which condition group
#'   to restrict to (e.g. `"Patient"`). Defaults to the non-reference group.
#' @param database Character. TF regulon database: `"collectri"`,
#'   `"dorothea"`, or `"both"`. Default: `"collectri"`.
#' @param method Character. decoupleR method: `"ulm"`, `"mlm"`, or `"wsum"`.
#'   Default: `"ulm"`.
#' @param species Character. Inherited from `obj$params$species` if `NULL`.
#' @param dorothea_levels Character vector. DoRothEA confidence levels.
#'   Default: `c("A", "B", "C")`.
#' @param interactive Logical. Prompt user to choose database and method.
#'
#' @return The OmniSeq object with TF results cached.
#' @export
run_omni_tf <- function(obj,
                        mode            = c("contrast", "samples"),
                        sample_group    = NULL,
                        database        = c("collectri", "dorothea", "both"),
                        method          = c("ulm", "mlm", "wsum"),
                        species         = NULL,
                        dorothea_levels = c("A", "B", "C"),
                        interactive     = FALSE) {

  validate_omniseq(obj)
  .require_pkg("decoupleR")

  mode     <- match.arg(mode)
  database <- match.arg(database)
  method   <- match.arg(method)

  # Resolve species
  raw_species <- species %||% obj$params$species
  organism <- if (grepl("sapiens|human", raw_species, ignore.case = TRUE)) {
    "human"
  } else if (grepl("musculus|mouse", raw_species, ignore.case = TRUE)) {
    "mouse"
  } else {
    message("[run_omni_tf] Unknown species '", raw_species,
            "' — defaulting to human.")
    "human"
  }

  if (interactive && base::interactive()) {
    cat("\n[TF Enrichment Setup]\n")
    ans_mode <- utils::menu(
      c("Contrast mode — t-statistics per contrast (recommended)",
        "Samples mode  — VST per sample, restricted to one group"),
      title = "Select TF inference mode:"
    )
    mode <- c("contrast", "samples")[ans_mode]

    ans_db <- utils::menu(
      c("CollecTRI (newer, more complete)",
        "DoRothEA (classic, confidence-filtered)",
        "Both (run and merge)"),
      title = "Select TF regulon database:"
    )
    database <- switch(ans_db, "collectri", "dorothea", "both")

    ans_met <- utils::menu(
      c("ULM — Univariate Linear Model (fast, recommended)",
        "MLM — Multivariate Linear Model",
        "WSUM — Weighted Sum"),
      title = "Select statistical method:"
    )
    method <- switch(ans_met, "ulm", "mlm", "wsum")
  }

  # ================================================================
  # REGULON CACHING — fetch ONCE, reuse for all subsequent contrasts
  # Only re-fetch if database or organism changed since last call.
  # This avoids 3 OmnipathR HTTP requests per contrast (21 total
  # for 7 contrasts → now 3 total for the first call only).
  # ================================================================
  message("\n[run_omni_tf] Mode: ", mode,
          " | Database: ", database,
          " | Method: ", method,
          " | Organism: ", organism)

  cached_db  <- obj$params$tf_database
  cached_org <- obj$params$tf_organism

  if (!is.null(obj$params$tf_regulons) &&
      nrow(obj$params$tf_regulons) > 0 &&
      identical(cached_db, database) &&
      identical(cached_org, organism)) {
    # Reuse cached regulons — no network call
    regulons <- obj$params$tf_regulons
    message("[run_omni_tf] Reusing cached regulons (",
            length(unique(regulons$source)), " TFs, ",
            length(unique(regulons$target)), " targets). No re-download.")
  } else {
    # First call or database/organism changed — fetch and cache
    message("[run_omni_tf] Fetching regulons from ", database, "...")
    regulons <- suppressWarnings(
      .omni_get_regulons(database, organism, dorothea_levels)
    )

    if (is.null(regulons) || nrow(regulons) == 0)
      stop("No regulons loaded. Check database/species arguments.",
           call. = FALSE)

    message("[run_omni_tf] Loaded ", length(unique(regulons$source)),
            " TF regulons covering ", length(unique(regulons$target)),
            " target genes.")

    # Cache regulons and parameters
    obj$params$tf_regulons <- regulons
    obj$params$tf_database <- database
    obj$params$tf_organism <- organism
  }

  obj$params$tf_method   <- method
  obj$params$tf_contrast <- obj$active_contrast

  # ==================================================================
  # MODE 1: CONTRAST — use Wald t-statistics as input
  # Produces one score per TF for the active contrast.
  # Genuinely different across contrasts — correct architecture.
  #
  # ROOT FIX: extract stat directly from the DESeq2 results object

  # stored in results_cache. Do NOT rely on a cached stat slot that
  # may be empty for contrasts built before the stat-preservation fix.
  # The results object is always present and always has the correct
  # genes for that specific contrast.
  #
  # NOTE ON GENE COUNTS: all contrasts from the same count matrix
  # will show a similar (often identical) gene count after NA removal.
  # This is EXPECTED — DESeq2's independent filtering depends on
  # gene-level mean counts (model-wide), not the specific contrast.
  # The stat VALUES are different; the stat fingerprint proves this.
  # ==================================================================
  if (mode == "contrast") {

    stat_vec    <- NULL
    cache_entry <- obj$results_cache[[obj$active_contrast]]

    if (is.null(cache_entry))
      stop("No results found for contrast '", obj$active_contrast,
           "' in results_cache. Run run_omni_de() or run_omni_contrast() first.",
           call. = FALSE)

    res_df <- tryCatch(as.data.frame(cache_entry$results),
                       error = function(e) NULL)

    if (!is.null(res_df)) {
      # Priority 1: stat column (Wald t-statistic — best input for ULM)
      if ("stat" %in% colnames(res_df)) {
        sv        <- as.numeric(res_df$stat)
        names(sv) <- rownames(res_df)
        sv        <- sv[!is.na(sv) & is.finite(sv)]
        if (length(sv) > 100) stat_vec <- sv
      }

      # Priority 2: LFC / lfcSE approximation
      # Used when stat was stripped by ashr/apeglm shrinkage
      if (is.null(stat_vec) &&
          all(c("log2FoldChange", "lfcSE") %in% colnames(res_df))) {
        lfc   <- as.numeric(res_df$log2FoldChange)
        lfcse <- as.numeric(res_df$lfcSE)
        ok    <- !is.na(lfc) & !is.na(lfcse) & lfcse > 0
        if (sum(ok) > 100) {
          sv        <- lfc[ok] / lfcse[ok]
          names(sv) <- rownames(res_df)[ok]
          sv        <- sv[is.finite(sv)]
          stat_vec  <- sv
          message("[run_omni_tf] stat column absent for '",
                  obj$active_contrast,
                  "'. Using LFC/lfcSE approximation (",
                  length(stat_vec), " genes).")
        }
      }
    }

    # Priority 3: raw fold changes — least accurate but never NULL
    if (is.null(stat_vec) && !is.null(cache_entry$fold_changes)) {
      stat_vec <- cache_entry$fold_changes
      message("[run_omni_tf] stat not available for '",
              obj$active_contrast,
              "'. Using fold_changes as proxy.")
    }

    if (is.null(stat_vec))
      stop("No stat, LFC, or fold_changes found for contrast '",
           obj$active_contrast, "'. Run run_omni_de() first.",
           call. = FALSE)

    stat_vec <- stat_vec[!is.na(stat_vec) & is.finite(stat_vec)]

    # ── STAT FINGERPRINT ─────────────────────────────────────────
    # Print a summary of the stat vector to PROVE each contrast
    # produces different input. Even when gene counts are identical
    # (expected — same count matrix), the median/sum/range differ.
    # ─────────────────────────────────────────────────────────────
    message("[run_omni_tf] Stat fingerprint for '", obj$active_contrast, "':",
            "\n  n_genes = ", length(stat_vec),
            " | median = ", round(stats::median(stat_vec), 4),
            " | sum = ", round(sum(stat_vec), 1),
            " | range = [", round(min(stat_vec), 2),
            ", ", round(max(stat_vec), 2), "]")

    # Build a 1-column matrix: genes x 1 contrast
    stat_mat <- matrix(stat_vec, ncol = 1,
                       dimnames = list(names(stat_vec),
                                       obj$active_contrast))

    message("[run_omni_tf] Running contrast-mode ULM on ",
            length(stat_vec), " genes...")

    tf_contrast <- decoupleR::run_ulm(
      mat     = stat_mat,
      network = regulons,
      .source = "source",
      .target = "target",
      .mor    = "mor",
      minsize = 5
    )

    # Store in per-contrast cache — each contrast fully independent
    if (is.null(obj$tf_contrast_cache)) obj$tf_contrast_cache <- list()
    obj$tf_contrast_cache[[obj$active_contrast]] <- tf_contrast

    # Legacy: store as tf_results in a format compatible with existing plots
    tf_res           <- tf_contrast
    tf_res$Condition <- obj$active_contrast

    if (is.null(obj$tf_cache)) obj$tf_cache <- list()
    obj$tf_cache[[obj$active_contrast]] <- tf_res
    obj$tf_results <- tf_res

    n_tfs <- length(unique(tf_res$source))
    message("[run_omni_tf] Contrast mode: scored ", n_tfs,
            " TFs for contrast '", obj$active_contrast, "'.")
    message("  Access with obj$tf_contrast_cache[['",
            obj$active_contrast, "']]")
  }

  # ==================================================================
  # MODE 2: SAMPLES — VST matrix restricted to one group
  # Never mixes patients and controls.
  # ==================================================================
  if (mode == "samples") {

    vsd <- omni_vst(obj)
    if (is.null(vsd))
      stop("VST not cached. Run run_omni_de() first.", call. = FALSE)

    meta     <- omni_metadata(obj)
    cond_col <- obj$params$condition
    groups   <- levels(meta[[cond_col]])

    # Default to non-reference (last) group — patients
    if (is.null(sample_group))
      sample_group <- groups[length(groups)]

    if (!sample_group %in% groups)
      stop("sample_group '", sample_group, "' not found. Available: ",
           paste(groups, collapse = ", "), call. = FALSE)

    # Restrict VST matrix to the requested group only
    group_samples <- rownames(meta)[meta[[cond_col]] == sample_group]
    mat_group     <- SummarizedExperiment::assay(vsd)[, group_samples,
                                                      drop = FALSE]

    message("[run_omni_tf] Sample mode: restricted to group '", sample_group,
            "' (", length(group_samples), " samples)...")

    tf_res <- switch(method,
                     "ulm"  = decoupleR::run_ulm(
                       mat = mat_group, network = regulons,
                       .source = "source", .target = "target",
                       .mor = "mor", minsize = 5),
                     "mlm"  = decoupleR::run_mlm(
                       mat = mat_group, network = regulons,
                       .source = "source", .target = "target",
                       .mor = "mor", minsize = 5),
                     "wsum" = decoupleR::run_wsum(
                       mat = mat_group, network = regulons,
                       .source = "source", .target = "target",
                       .mor = "mor", minsize = 5)
    )

    # Annotate with sample metadata
    matched_idx      <- match(tf_res$condition, rownames(meta))
    tf_res$Condition <- meta[[cond_col]][matched_idx]
    tf_res$Group     <- sample_group

    # Store in dedicated sample activity slot
    obj$tf_sample_activity <- tf_res

    # Also update legacy slots
    if (is.null(obj$tf_cache)) obj$tf_cache <- list()
    obj$tf_cache[[obj$active_contrast]] <- tf_res
    obj$tf_results <- tf_res

    n_tfs <- length(unique(tf_res$source))
    message("[run_omni_tf] Sample mode: scored ", n_tfs,
            " TFs across ", length(group_samples), " samples.")
    message("  Access with obj$tf_sample_activity")
  }

  return(obj)
}


# =============================================================================
# 2. TF ACTIVITY BARPLOT
# =============================================================================

#' Plot Top Activated and Repressed Transcription Factors
#'
#' For contrast mode: draws from `tf_contrast_cache` — genuinely different
#' per contrast. For samples mode: draws from `tf_sample_activity`.
#' Color scale is clamped to the range -5 to 5 to prevent extreme outliers
#' from compressing all other TFs to grey.
#'
#' @param obj A validated OmniSeq object with TF results.
#' @param n Integer. Number of top TFs per direction to show. Default: `15`.
#' @param group Character. For samples mode: condition group to plot.
#'   Defaults to the non-reference group. Ignored in contrast mode.
#' @param min_regulon Integer. Minimum regulon size to include a TF.
#'   Removes noise TFs with very few targets (e.g. 5-gene regulons that score
#'   high by chance). Default: `15`. Set to `0` to disable.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_tf_activity <- function(obj, n = 30, group = NULL,
                                  min_regulon = 15L) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  # Prefer contrast cache (one score per TF, genuinely per-contrast)
  # Fall back to sample-level cache for backward compat
  contrast_res <- if (!is.null(obj$tf_contrast_cache))
    obj$tf_contrast_cache[[obj$active_contrast]]
  else
    NULL

  if (!is.null(contrast_res)) {
    # --- Contrast mode: one score per TF ---
    tf_df <- data.frame(
      TF        = contrast_res$source,
      MeanScore = contrast_res$score,
      stringsAsFactors = FALSE
    )
    tf_df <- tf_df[!is.na(tf_df$MeanScore), ]
    tf_df <- tf_df[order(tf_df$MeanScore), ]

    subtitle_text <- paste0("Contrast: ", obj$active_contrast,
                            " | DB: ", obj$params$tf_database %||% "?",
                            " | Method: ", obj$params$tf_method %||% "ulm")

  } else {
    # --- Samples mode fallback ---
    tf_res <- .omni_get_tf(obj)
    if (is.null(tf_res))
      stop("No TF results for contrast '", obj$active_contrast,
           "'. Run run_omni_tf() first.", call. = FALSE)

    available_conditions <- unique(as.character(tf_res$Condition))
    available_conditions <- available_conditions[!is.na(available_conditions)]

    if (is.null(group)) {
      non_ctrl <- available_conditions[!grepl("Control|Ctrl|ctrl",
                                              available_conditions)]
      group <- if (length(non_ctrl) > 0) non_ctrl[1] else
        available_conditions[1]
    }

    if (!group %in% available_conditions)
      stop("Group '", group, "' not found in TF results. Available: ",
           paste(available_conditions, collapse = ", "), call. = FALSE)

    df_group <- tf_res[tf_res$Condition == group, ]
    tf_mean  <- tapply(df_group$score, df_group$source, mean, na.rm = TRUE)

    tf_df <- data.frame(
      TF        = names(tf_mean),
      MeanScore = as.numeric(tf_mean),
      stringsAsFactors = FALSE
    )
    tf_df <- tf_df[!is.na(tf_df$MeanScore), ]
    tf_df <- tf_df[order(tf_df$MeanScore), ]

    subtitle_text <- paste0("Group: ", group,
                            " | DB: ", obj$params$tf_database %||% "?",
                            " | Method: ", obj$params$tf_method %||% "ulm")
  }

  # Apply minimum regulon size filter
  if (min_regulon > 0 && !is.null(obj$params$tf_regulons)) {
    reg_sizes <- table(obj$params$tf_regulons$source)
    keep_tfs  <- names(reg_sizes)[reg_sizes >= min_regulon]
    tf_df     <- tf_df[tf_df$TF %in% keep_tfs, ]
    if (nrow(tf_df) == 0)
      stop("No TFs remaining after min_regulon filter (", min_regulon,
           "). Lower min_regulon.", call. = FALSE)
  }

  # Top N activated + top N repressed
  top_act <- utils::tail(tf_df, n)
  top_rep <- utils::head(tf_df, n)
  plot_df <- rbind(top_rep, top_act)
  plot_df <- plot_df[!duplicated(plot_df$TF), ]
  plot_df$TF        <- factor(plot_df$TF, levels = plot_df$TF)
  plot_df$Direction <- ifelse(plot_df$MeanScore > 0, "Activated", "Repressed")

  # Dynamic clamping
  .require_pkg("scales")
  all_scores_full <- tf_df$MeanScore
  clamp_abs <- max(abs(stats::quantile(all_scores_full, c(0.02, 0.98),
                                       na.rm = TRUE)))
  plot_df$ScoreClamped <- scales::squish(plot_df$MeanScore,
                                         range = c(-clamp_abs, clamp_abs))

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = MeanScore, y = TF,
                                    fill = Direction)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("Activated" = "#E74C3C", "Repressed" = "#2980B9")
    ) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5) +
    ggplot2::coord_cartesian(xlim = c(-clamp_abs, clamp_abs)) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9, face = "italic")
    ) +
    ggplot2::labs(
      title    = paste("TF Activity —", obj$active_contrast),
      subtitle = subtitle_text,
      x        = "Activity Score",
      y        = NULL,
      fill     = NULL
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 3. TF NETWORK PLOT
# =============================================================================

#' Plot TF Regulatory Network
#'
#' Visualises the top activated/repressed TFs and their significant DEG
#' target genes as a network. Reads regulons from the object cache —
#' no re-download required.
#'
#' @param obj A validated OmniSeq object with TF and DE results.
#' @param n_tfs Integer. Number of top TFs to include (split equally between
#'   activated and repressed). Default: `10`.
#' @param n_targets Integer. Top target genes per TF to show. Default: `5`.
#' @param plot_type Character. `"static"` (ggraph) or `"interactive"`
#'   (visNetwork). Default: `"static"`.
#' @param group Character. Condition group to use for TF ranking.
#' @param min_regulon Integer. Minimum regulon size filter. Default: `15`.
#'
#' @return A plot object (returned invisibly).
#' @export
plot_omni_tf_network <- function(obj,
                                 n_tfs       = 25,
                                 n_targets   = 5,
                                 plot_type   = c("static", "interactive"),
                                 group       = NULL,
                                 min_regulon = 15L) {

  validate_omniseq(obj)
  plot_type <- match.arg(plot_type)

  # Prefer contrast-level scores; fall back to sample-level
  contrast_res <- if (!is.null(obj$tf_contrast_cache))
    obj$tf_contrast_cache[[obj$active_contrast]]
  else
    NULL

  tf_res <- .omni_get_tf(obj)
  res    <- omni_results(obj)

  if (is.null(tf_res) && is.null(contrast_res))
    stop("No TF results for contrast '", obj$active_contrast,
         "'. Run run_omni_tf() first.", call. = FALSE)
  if (is.null(res))
    stop("No DE results. Run run_omni_de() first.", call. = FALSE)

  meta     <- omni_metadata(obj)
  cond_col <- obj$params$condition
  groups   <- levels(meta[[cond_col]])
  if (is.null(group)) group <- groups[length(groups)]

  # Get TF scores — prefer contrast cache
  if (!is.null(contrast_res)) {
    tf_mean <- stats::setNames(contrast_res$score, contrast_res$source)
  } else {
    df_group <- tf_res[tf_res$Condition == group, ]
    tf_mean  <- tapply(df_group$score, df_group$source, mean, na.rm = TRUE)
  }

  # Dynamic clamping
  .require_pkg("scales")

  # Apply minimum regulon size filter before selecting top TFs
  if (min_regulon > 0 && !is.null(obj$params$tf_regulons)) {
    reg_sizes_net <- table(obj$params$tf_regulons$source)
    keep_net      <- names(reg_sizes_net)[reg_sizes_net >= min_regulon]
    tf_mean       <- tf_mean[names(tf_mean) %in% keep_net]
  }

  clamp_abs       <- max(abs(stats::quantile(tf_mean, c(0.02, 0.98),
                                             na.rm = TRUE)))
  tf_mean_clamped <- scales::squish(tf_mean, range = c(-clamp_abs, clamp_abs))

  # Select TFs: top N activated + top N repressed separately
  tf_sorted    <- sort(tf_mean, decreasing = TRUE)
  top_act_tfs  <- names(utils::head(tf_sorted, ceiling(n_tfs / 2)))
  top_rep_tfs  <- names(utils::tail(tf_sorted, floor(n_tfs / 2)))
  top_tfs      <- unique(c(top_act_tfs, top_rep_tfs))

  # Use cached regulons — no re-download
  regulons <- obj$params$tf_regulons
  if (is.null(regulons)) {
    message("[plot_omni_tf_network] Regulon cache empty — re-fetching...")
    regulons <- suppressWarnings(
      .omni_get_regulons(
        obj$params$tf_database  %||% "collectri",
        obj$params$tf_organism  %||% "human",
        c("A", "B", "C")
      )
    )
    obj$params$tf_regulons <- regulons
  }

  # DE results for target gene coloring
  res_df      <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)

  # Build edge list
  edges <- do.call(rbind, lapply(top_tfs, function(tf) {
    targets     <- regulons$target[regulons$source == tf]
    sig_targets <- targets[targets %in%
                             res_df$Gene[res_df$padj < 0.05]]
    if (length(sig_targets) == 0) return(NULL)
    lfc_vals    <- abs(res_df$log2FoldChange[match(sig_targets, res_df$Gene)])
    sig_targets <- sig_targets[order(lfc_vals, decreasing = TRUE)]
    sig_targets <- utils::head(sig_targets, n_targets)
    data.frame(from = tf, to = sig_targets, stringsAsFactors = FALSE)
  }))

  if (is.null(edges) || nrow(edges) == 0)
    stop("No TF-target edges found among significant DEGs.", call. = FALSE)

  # Node attributes
  all_nodes  <- unique(c(edges$from, edges$to))
  node_type  <- ifelse(all_nodes %in% top_tfs, "TF", "Gene")
  node_score <- sapply(all_nodes, function(n) {
    if (n %in% names(tf_mean_clamped)) tf_mean_clamped[n]
    else {
      lfc <- res_df$log2FoldChange[match(n, res_df$Gene)] %||% 0
      scales::squish(lfc, range = c(-clamp_abs, clamp_abs))
    }
  })

  # ── STATIC (ggraph) ──────────────────────────────────────────────────────
  if (plot_type == "static") {
    .require_pkg("igraph")
    .require_pkg("ggraph")

    g <- igraph::graph_from_data_frame(
      edges,
      directed = TRUE,
      vertices = data.frame(name  = all_nodes,
                            type  = node_type,
                            score = node_score)
    )

    p <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(
        ggplot2::aes(),
        arrow     = grid::arrow(length = grid::unit(3, "mm"),
                                type = "closed"),
        end_cap   = ggraph::circle(4, "mm"),
        color     = "grey60",
        linewidth = 0.4,
        alpha     = 0.7
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(color = score, size = type),
        show.legend = TRUE
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label    = name,
                     fontface = ifelse(type == "TF", "bold", "plain")),
        repel = TRUE, size = 3
      ) +
      ggplot2::scale_color_gradient2(
        low      = "#2980B9",
        mid      = "grey90",
        high     = "#E74C3C",
        midpoint = 0,
        limits   = c(-clamp_abs, clamp_abs),
        oob      = scales::squish,
        name     = sprintf("Score / LFC\n(clamped \u00b1%.1f)", clamp_abs)
      ) +
      ggplot2::scale_size_manual(
        values = c("TF" = 6, "Gene" = 3), name = "Node type"
      ) +
      ggraph::theme_graph(base_size = 12, base_family = "sans") +
      ggplot2::labs(
        title   = paste("TF Regulatory Network —", obj$active_contrast),
        subtitle = paste0("Top ", n_tfs, " TFs | ", nrow(edges),
                          " TF-target edges | Color clamped to \u00b1",
                          round(clamp_abs, 1)),
        caption = paste0("TF nodes (large): activity score | ",
                         "Gene nodes (small): LFC | Red = up, Blue = down")
      )

    print(p)
    return(invisible(p))
  }

  # ── INTERACTIVE (visNetwork) ─────────────────────────────────────────────
  if (plot_type == "interactive") {
    .require_pkg("visNetwork")

    .score_to_hex <- function(s) {
      if (is.na(s) || s == 0) return("#CCCCCC")
      if (s > 0) "#E74C3C" else "#2980B9"
    }

    vis_nodes <- data.frame(
      id        = all_nodes,
      label     = all_nodes,
      group     = node_type,
      value     = abs(node_score) + 1,
      color     = sapply(node_score, .score_to_hex),
      font.bold = ifelse(node_type == "TF", TRUE, FALSE),
      title     = paste0("<b>", all_nodes, "</b><br>",
                         "Type: ",  node_type, "<br>",
                         "Score: ", round(node_score, 3)),
      stringsAsFactors = FALSE
    )

    vis_edges <- data.frame(
      from   = edges$from,
      to     = edges$to,
      arrows = "to",
      color  = "grey70",
      width  = 1.2,
      stringsAsFactors = FALSE
    )

    p <- visNetwork::visNetwork(
      vis_nodes, vis_edges,
      main   = paste("TF Network —", group),
      height = "600px",
      width  = "100%"
    ) |>
      visNetwork::visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1),
        nodesIdSelection = TRUE
      ) |>
      visNetwork::visPhysics(
        solver           = "forceAtlas2Based",
        forceAtlas2Based = list(gravitationalConstant = -50)
      ) |>
      visNetwork::visEdges(
        smooth = list(type = "curvedCW", roundness = 0.2)
      ) |>
      visNetwork::visLegend()

    print(p)
    return(invisible(p))
  }
}


# =============================================================================
# 4. TF GENOTYPE ANALYSIS — MULTI-PANEL (redesigned v2)
#    Panels: A (pattern), B (shared repressed), C (shared activated per pair),
#            D (discordant per pair + combined), E (heatmap),
#            F (individual genotype profiles), G (genotype-specific TFs)
# =============================================================================

#' Plot TF Activity by Genotype — Independent Panels
#'
#' Returns a named list of independent ggplot objects. Each is printed,
#' and optionally saved as a separate PNG.
#'
#' @param obj A validated OmniSeq object with `tf_contrast_cache`.
#' @param main_contrast Character. Auto-detected if NULL.
#' @param geno_contrasts Named list. Auto-detected if NULL.
#' @param n_top Integer. Top TFs per panel. Default: 20.
#' @param min_regulon Integer. Minimum regulon size. Default: 15.
#' @param p_thresh Numeric. Default: 0.05.
#' @param use_msigdb Logical. Default: TRUE.
#' @param save_dir Character or NULL. Save each plot as PNG.
#'
#' @return A named list of ggplot objects (invisibly).
#' @export
plot_omni_tf_genotype <- function(obj,
                                  main_contrast  = NULL,
                                  geno_contrasts = NULL,
                                  n_top          = 20L,
                                  min_regulon    = 15L,
                                  p_thresh       = 0.05,
                                  use_msigdb     = TRUE,
                                  save_dir       = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  cache <- obj$tf_contrast_cache
  if (is.null(cache) || length(cache) == 0)
    stop("No TF contrast cache. Run run_omni_tf() first.", call. = FALSE)

  # ── Resolve main contrast ──────────────────────────────────────────
  if (is.null(main_contrast)) {
    keys <- names(cache)
    main_contrast <- keys[grepl("Patient.*Control|Condition", keys,
                                ignore.case = TRUE)][1]
    if (is.na(main_contrast)) main_contrast <- keys[1]
  }
  if (!main_contrast %in% names(cache))
    stop("main_contrast '", main_contrast,
         "' not in tf_contrast_cache.", call. = FALSE)

  message("[plot_omni_tf_genotype] Main contrast: ", main_contrast)

  # ── Resolve genotype vs-control contrasts ──────────────────────────
  if (is.null(geno_contrasts)) {
    other_keys <- setdiff(names(cache), main_contrast)
    ctrl_keys <- other_keys[grepl("vs.Ctrl$|vs.Control$",
                                  other_keys, ignore.case = TRUE)]
    ctrl_keys <- ctrl_keys[!grepl("excl|mask|schoggins|isgs",
                                  ctrl_keys, ignore.case = TRUE)]
    if (length(ctrl_keys) >= 2) {
      geno_labels <- gsub("\\s*vs\\s*(Ctrl|Control).*$", "",
                          ctrl_keys, ignore.case = TRUE)
      geno_labels <- trimws(geno_labels)
      geno_contrasts <- as.list(stats::setNames(ctrl_keys, geno_labels))
      message("[plot_omni_tf_genotype] Auto-detected: ",
              paste(ctrl_keys, collapse = ", "))
    }
  }

  if (is.null(geno_contrasts) || length(geno_contrasts) < 2)
    stop("Need >= 2 genotype vs control contrasts.", call. = FALSE)
  if (is.null(names(geno_contrasts)) || any(!nzchar(names(geno_contrasts))))
    stop("geno_contrasts must be named (genotype label = contrast name), ",
         "e.g. list(KO = \"KO vs Ctrl\", Het = \"Het vs Ctrl\").",
         call. = FALSE)

  gt_labels <- names(geno_contrasts)
  message("[plot_omni_tf_genotype] Genotypes: ",
          paste(gt_labels, collapse = ", "))

  # ── Regulon sizes ──────────────────────────────────────────────────
  reg       <- obj$params$tf_regulons
  reg_sizes <- if (!is.null(reg)) table(reg$source) else NULL
  keep_tfs_global <- if (!is.null(reg_sizes) && min_regulon > 0)
    names(reg_sizes[reg_sizes >= min_regulon]) else NULL

  # ── Main contrast TF scores ───────────────────────────────────────
  main_raw <- as.data.frame(cache[[main_contrast]])
  colnames(main_raw)[colnames(main_raw) == "source"] <- "TF"

  if (!is.null(keep_tfs_global))
    main_raw <- main_raw[main_raw$TF %in% keep_tfs_global, ]

  message(sprintf("[plot_omni_tf_genotype] After min_regulon=%d: %d TFs",
                  min_regulon, nrow(main_raw)))

  sig_df <- main_raw[!is.na(main_raw$p_value) &
                       main_raw$p_value < p_thresh, ]
  sig_df <- sig_df[order(sig_df$score), ]
  message(sprintf("[plot_omni_tf_genotype] Significant (p<%g): %d TFs",
                  p_thresh, nrow(sig_df)))

  # ── MSigDB IFN cross-reference ────────────────────────────────────
  ifn_universe <- character(0)
  ifn_core <- c("STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B",
                "STAT6","IRF1","IRF2","IRF3","IRF4","IRF5","IRF7",
                "IRF8","IRF9")

  if (use_msigdb && requireNamespace("msigdbr", quietly = TRUE)) {
    tryCatch({
      ms <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
      ifn_universe <- unique(c(
        ms$gene_symbol[ms$gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE"],
        ms$gene_symbol[ms$gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE"]
      ))
      message(sprintf("[plot_omni_tf_genotype] MSigDB IFN: %d genes",
                      length(ifn_universe)))
    }, error = function(e)
      message("[plot_omni_tf_genotype] msigdbr failed: ",
              conditionMessage(e)))
  }

  sig_df$db_class <- dplyr::case_when(
    sig_df$TF %in% ifn_core     ~ "Core IFN (IRF/STAT)",
    sig_df$TF %in% ifn_universe ~ "IFN-related (MSigDB)",
    TRUE                         ~ "Novel"
  )
  sig_df$direction <- ifelse(sig_df$score > 0, "Activated", "Repressed")

  # ── Per-genotype scores ────────────────────────────────────────────
  gt_all_scores <- do.call(rbind, lapply(gt_labels, function(gt) {
    key <- geno_contrasts[[gt]]
    if (!key %in% names(cache)) return(NULL)
    df <- as.data.frame(cache[[key]])
    colnames(df)[colnames(df) == "source"] <- "TF"
    if (!is.null(keep_tfs_global))
      df <- df[df$TF %in% keep_tfs_global, ]
    df$Genotype  <- gt
    df$sig       <- !is.na(df$p_value) & df$p_value < p_thresh
    df$direction <- ifelse(df$score > 0, "Activated", "Repressed")
    df[, c("TF", "score", "p_value", "sig", "direction", "Genotype")]
  }))

  if (is.null(gt_all_scores) || nrow(gt_all_scores) == 0)
    stop("Could not extract per-genotype TF scores.", call. = FALSE)

  # Restrict to significant TFs from main contrast for C/D/E panels
  gt_scores <- gt_all_scores[gt_all_scores$TF %in% sig_df$TF, ]

  # ── Wide format for pattern classification ─────────────────────────
  gt_wide <- tidyr::pivot_wider(
    gt_scores[, c("TF", "Genotype", "score", "direction")],
    names_from  = "Genotype",
    values_from = c("score", "direction"),
    names_sep   = "_"
  )
  dir_cols <- grep("^direction_", colnames(gt_wide), value = TRUE)

  gt_wide$pattern <- apply(gt_wide, 1, function(r) {
    dirs <- sapply(dir_cols, function(col) {
      v <- r[col]; if (is.na(v)) "ND" else v
    })
    n_act <- sum(dirs == "Activated")
    n_rep <- sum(dirs == "Repressed")
    n_gt  <- length(dirs)
    if (n_act == n_gt) return("Shared activated")
    if (n_rep == n_gt) return("Shared repressed")
    if (n_act > 0 && n_rep > 0) return("Discordant")
    return("Mixed / ND")
  })

  pattern_df <- merge(gt_wide,
                      sig_df[, c("TF", "db_class", "score")],
                      by = "TF", suffixes = c("_gt", "_main"))

  # ── Colour constants ───────────────────────────────────────────────
  DB_COLORS <- c("Core IFN (IRF/STAT)"  = "#E74C3C",
                 "IFN-related (MSigDB)" = "#E67E22",
                 "Novel"                = "#2166AC")
  GT_COLORS <- c("#E74C3C","#E67E22","#8E44AD","#27AE60",
                 "#2980B9","#16A085")
  gt_col_vals <- stats::setNames(GT_COLORS[seq_along(gt_labels)],
                                 gt_labels)
  PAT_COLORS <- c("Shared repressed" = "#2166AC",
                  "Shared activated" = "#B2182B",
                  "Discordant"       = "#E67E22",
                  "Mixed / ND"       = "#AAAAAA")

  plots <- list()
  pairs <- utils::combn(gt_labels, 2, simplify = FALSE)

  # ════════════════════════════════════════════════════════════════════
  # A — PATTERN DISTRIBUTION
  # ════════════════════════════════════════════════════════════════════
  pat_counts <- as.data.frame(table(Pattern = pattern_df$pattern))
  pat_counts <- pat_counts[pat_counts$Freq > 0, ]
  pat_counts$Pattern <- reorder(pat_counts$Pattern, pat_counts$Freq)

  plots$A_pattern_distribution <- ggplot2::ggplot(
    pat_counts, ggplot2::aes(x = Freq, y = Pattern, fill = Pattern)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = Freq),
                       hjust = -0.2, size = 5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = PAT_COLORS, guide = "none") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    ) +
    ggplot2::labs(
      title = paste0("TF Genotype Pattern Distribution (n=",
                     nrow(sig_df), " significant TFs)"),
      subtitle = paste0("Genotypes: ", paste(gt_labels, collapse = ", ")),
      x = "Number of TFs", y = NULL
    )

  # ════════════════════════════════════════════════════════════════════
  # B — TOP SHARED REPRESSED
  # ════════════════════════════════════════════════════════════════════
  shared_rep <- pattern_df[pattern_df$pattern == "Shared repressed", ]
  shared_rep <- merge(shared_rep, sig_df[, c("TF", "score", "db_class")],
                      by = "TF", suffixes = c("", "_m2"))
  shared_rep <- shared_rep[order(shared_rep$score), ]
  shared_rep <- utils::head(shared_rep, n_top)

  if (nrow(shared_rep) > 0) {
    shared_rep$TF <- reorder(shared_rep$TF, shared_rep$score)
    plots$B_top_shared_repressed <- ggplot2::ggplot(
      shared_rep, ggplot2::aes(x = score, y = TF, fill = db_class)
    ) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
      ggplot2::scale_fill_manual(values = DB_COLORS, name = NULL) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(face = "italic", size = 10),
        plot.title  = ggplot2::element_text(face = "bold", size = 14),
        legend.position = "bottom"
      ) +
      ggplot2::labs(
        title = paste0("Top ", nrow(shared_rep),
                       " Shared Repressed TFs (all genotypes)"),
        x = "ULM activity score", y = NULL
      )
  }

  # ════════════════════════════════════════════════════════════════════
  # C — SHARED ACTIVATED PER PAIR
  # ════════════════════════════════════════════════════════════════════
  for (pair in pairs) {
    g1 <- pair[1]; g2 <- pair[2]
    pair_label <- paste0(g1, "_and_", g2)
    dir_g1 <- paste0("direction_", g1); dir_g2 <- paste0("direction_", g2)
    sc_g1  <- paste0("score_", g1);     sc_g2  <- paste0("score_", g2)

    if (!all(c(dir_g1, dir_g2, sc_g1, sc_g2) %in% colnames(gt_wide))) next

    act_pair <- gt_wide[!is.na(gt_wide[[dir_g1]]) &
                          !is.na(gt_wide[[dir_g2]]) &
                          gt_wide[[dir_g1]] == "Activated" &
                          gt_wide[[dir_g2]] == "Activated", ]

    if (nrow(act_pair) == 0) {
      plots[[paste0("C_shared_activated_", pair_label)]] <-
        ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::annotate("text", x = .5, y = .5, size = 5,
                          label = paste0("No shared activated TFs: ",
                                         g1, " & ", g2))
      next
    }

    act_pair$mean_score <- (as.numeric(act_pair[[sc_g1]]) +
                              as.numeric(act_pair[[sc_g2]])) / 2
    act_pair <- act_pair[order(-act_pair$mean_score), ]
    act_pair <- utils::head(act_pair, n_top)

    act_long <- data.frame(
      TF       = rep(act_pair$TF, 2),
      Score    = c(as.numeric(act_pair[[sc_g1]]),
                   as.numeric(act_pair[[sc_g2]])),
      Genotype = rep(c(g1, g2), each = nrow(act_pair)),
      stringsAsFactors = FALSE
    )
    act_long$TF <- factor(act_long$TF,
                          levels = act_pair$TF[order(act_pair$mean_score)])

    plots[[paste0("C_shared_activated_", pair_label)]] <-
      ggplot2::ggplot(act_long,
                      ggplot2::aes(x = Score, y = TF, color = Genotype)) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50") +
      ggplot2::geom_line(ggplot2::aes(group = TF),
                         color = "grey70", linewidth = 0.5) +
      ggplot2::geom_point(size = 3.5, alpha = 0.9) +
      ggplot2::scale_color_manual(values = gt_col_vals) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(face = "italic", size = 10),
        plot.title  = ggplot2::element_text(face = "bold", size = 14),
        legend.position = "top"
      ) +
      ggplot2::labs(
        title = paste0("Top ", nrow(act_pair),
                       " Shared Activated: ", g1, " & ", g2),
        x = "ULM activity score", y = NULL
      )
  }

  # ════════════════════════════════════════════════════════════════════
  # D — DISCORDANT TFs PER PAIR
  # ════════════════════════════════════════════════════════════════════
  d_plots <- list()

  for (pair in pairs) {
    g1 <- pair[1]; g2 <- pair[2]
    pair_label <- paste0(g1, "_vs_", g2)
    dir_g1 <- paste0("direction_", g1); dir_g2 <- paste0("direction_", g2)
    sc_g1  <- paste0("score_", g1);     sc_g2  <- paste0("score_", g2)

    if (!all(c(dir_g1, dir_g2) %in% colnames(gt_wide))) next

    disc_pair <- gt_wide[
      !is.na(gt_wide[[dir_g1]]) & !is.na(gt_wide[[dir_g2]]) &
        gt_wide[[dir_g1]] != gt_wide[[dir_g2]], ]

    if (nrow(disc_pair) == 0) {
      p_empty <- ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::annotate("text", x = .5, y = .5, size = 5,
                          label = paste0("No discordant TFs: ",
                                         g1, " vs ", g2))
      plots[[paste0("D_discordant_", pair_label)]] <- p_empty
      d_plots[[pair_label]] <- p_empty
      next
    }

    disc_pair$diff <- abs(as.numeric(disc_pair[[sc_g1]]) -
                            as.numeric(disc_pair[[sc_g2]]))
    disc_pair <- disc_pair[order(-disc_pair$diff), ]
    disc_pair <- utils::head(disc_pair, n_top)

    disc_long <- gt_scores[gt_scores$TF %in% disc_pair$TF &
                             gt_scores$Genotype %in% c(g1, g2), ]

    g1_sc <- stats::setNames(as.numeric(disc_pair[[sc_g1]]), disc_pair$TF)
    disc_long$TF <- factor(disc_long$TF, levels = names(sort(g1_sc)))

    p_disc <- ggplot2::ggplot(
      disc_long, ggplot2::aes(x = score, y = TF, color = Genotype)
    ) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50") +
      ggplot2::geom_line(ggplot2::aes(group = TF),
                         color = "grey70", linewidth = 0.5) +
      ggplot2::geom_point(ggplot2::aes(shape = sig),
                          size = 3.5, alpha = 0.9) +
      ggplot2::scale_color_manual(values = gt_col_vals) +
      ggplot2::scale_shape_manual(
        values = c("TRUE" = 16, "FALSE" = 1),
        labels = c("TRUE" = "p<0.05", "FALSE" = "ns"), name = "Sig.") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(face = "italic", size = 10),
        plot.title  = ggplot2::element_text(face = "bold", size = 14),
        legend.position = "top"
      ) +
      ggplot2::labs(
        title = paste0("Discordant TFs: ", g1, " vs ", g2,
                       " (top ", nrow(disc_pair), ")"),
        subtitle = "Opposite direction between genotypes | Filled = p<0.05",
        x = "ULM activity score (vs controls)", y = NULL
      )

    plots[[paste0("D_discordant_", pair_label)]] <- p_disc
    d_plots[[pair_label]] <- p_disc
  }

  # Combined discordant assembly
  if (length(d_plots) >= 2 &&
      requireNamespace("patchwork", quietly = TRUE)) {
    plots$D_discordant_combined <- Reduce(`+`, d_plots) +
      patchwork::plot_layout(ncol = min(length(d_plots), 3)) +
      patchwork::plot_annotation(
        title = "Discordant TFs — All Pairwise Comparisons",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 16)))
  }

  # ════════════════════════════════════════════════════════════════════
  # E — PER-GENOTYPE HEATMAP
  # ════════════════════════════════════════════════════════════════════
  top_act_tfs <- utils::tail(sig_df$TF[order(sig_df$score)], n_top)
  top_rep_tfs <- utils::head(sig_df$TF[order(sig_df$score)], n_top)
  focus_tfs   <- unique(c(top_rep_tfs, top_act_tfs))

  heat_df <- gt_scores[gt_scores$TF %in% focus_tfs, ]
  tf_order_e <- sig_df$TF[sig_df$TF %in% focus_tfs]
  heat_df$TF <- factor(heat_df$TF, levels = rev(tf_order_e))
  heat_df$Genotype <- factor(heat_df$Genotype, levels = gt_labels)
  clamp_val <- stats::quantile(abs(heat_df$score), 0.95, na.rm = TRUE)
  heat_df$score_c <- pmax(pmin(heat_df$score, clamp_val), -clamp_val)

  plots$E_per_genotype_heatmap <- ggplot2::ggplot(
    heat_df, ggplot2::aes(x = Genotype, y = TF, fill = score_c)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(!is.na(score), round(score, 1), ""),
        color = ifelse(sig, "black", "grey60")),
      size = 2.5, fontface = "bold") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "ULM score", na.value = "grey92") +
    ggplot2::scale_color_identity() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 11),
      axis.text.y = ggplot2::element_text(face = "italic", size = 7),
      plot.title  = ggplot2::element_text(face = "bold", size = 14)
    ) +
    ggplot2::labs(
      title = paste0("Per-Genotype TF Scores (top ", n_top,
                     " activated + repressed)"),
      subtitle = "Black = p<0.05 | Grey = ns", x = NULL, y = NULL
    )

  # ════════════════════════════════════════════════════════════════════
  # F — INDIVIDUAL GENOTYPE PROFILES (NEW)
  #   For each genotype: top N activated + top N repressed TFs
  #   from that genotype's own contrast vs controls
  # ════════════════════════════════════════════════════════════════════
  for (gt in gt_labels) {
    key <- geno_contrasts[[gt]]
    if (!key %in% names(cache)) next

    gt_df <- as.data.frame(cache[[key]])
    colnames(gt_df)[colnames(gt_df) == "source"] <- "TF"
    if (!is.null(keep_tfs_global))
      gt_df <- gt_df[gt_df$TF %in% keep_tfs_global, ]

    gt_df <- gt_df[order(gt_df$score), ]
    f_rep <- utils::head(gt_df, n_top)
    f_act <- utils::tail(gt_df, n_top)
    f_df  <- rbind(f_rep, f_act)
    f_df  <- f_df[!duplicated(f_df$TF), ]
    f_df$TF <- factor(f_df$TF, levels = f_df$TF)
    f_df$Direction <- ifelse(f_df$score > 0, "Activated", "Repressed")

    # Add db_class
    f_df <- merge(f_df, sig_df[, c("TF", "db_class")],
                  by = "TF", all.x = TRUE)
    f_df$db_class[is.na(f_df$db_class)] <- "Novel"
    f_df <- f_df[order(f_df$score), ]
    f_df$TF <- factor(f_df$TF, levels = f_df$TF)

    plots[[paste0("F_profile_", gt)]] <- ggplot2::ggplot(
      f_df, ggplot2::aes(x = score, y = TF, fill = db_class)
    ) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
      ggplot2::scale_fill_manual(values = DB_COLORS, name = NULL) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(face = "italic", size = 9),
        plot.title  = ggplot2::element_text(face = "bold", size = 14),
        legend.position = "bottom"
      ) +
      ggplot2::labs(
        title = paste0("TF Profile: ", gt, " vs Controls"),
        subtitle = paste0("Top ", n_top,
                          " activated + repressed | Contrast: ", key),
        x = "ULM activity score", y = NULL
      )
  }

  # ════════════════════════════════════════════════════════════════════
  # G — GENOTYPE-SPECIFIC TFs (NEW)
  #   TFs significant in ONE genotype but NOT significant in any other.
  #   Identifies mutation-specific regulatory programs.
  # ════════════════════════════════════════════════════════════════════
  for (gt in gt_labels) {
    key <- geno_contrasts[[gt]]
    if (!key %in% names(cache)) next

    gt_df <- as.data.frame(cache[[key]])
    colnames(gt_df)[colnames(gt_df) == "source"] <- "TF"
    if (!is.null(keep_tfs_global))
      gt_df <- gt_df[gt_df$TF %in% keep_tfs_global, ]

    # TFs significant in this genotype
    sig_in_gt <- gt_df$TF[!is.na(gt_df$p_value) &
                            gt_df$p_value < p_thresh]

    # Must NOT be significant in any other genotype
    other_gts <- setdiff(gt_labels, gt)
    specific_tfs <- sig_in_gt

    for (ogt in other_gts) {
      okey <- geno_contrasts[[ogt]]
      if (!okey %in% names(cache)) next
      odf <- as.data.frame(cache[[okey]])
      colnames(odf)[colnames(odf) == "source"] <- "TF"
      # TFs NOT significant in this other genotype
      nonsig_other <- odf$TF[is.na(odf$p_value) | odf$p_value >= p_thresh]
      specific_tfs <- intersect(specific_tfs, nonsig_other)
    }

    if (length(specific_tfs) == 0) {
      plots[[paste0("G_specific_", gt)]] <-
        ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::annotate("text", x = .5, y = .5, size = 5,
                          label = paste0("No genotype-specific TFs\n",
                                         "unique to ", gt))
      next
    }

    spec_df <- gt_df[gt_df$TF %in% specific_tfs, ]
    spec_df <- spec_df[order(spec_df$score), ]
    spec_df <- rbind(utils::head(spec_df, n_top),
                     utils::tail(spec_df, n_top))
    spec_df <- spec_df[!duplicated(spec_df$TF), ]

    # Add db_class
    spec_df <- merge(spec_df, sig_df[, c("TF", "db_class")],
                     by = "TF", all.x = TRUE)
    spec_df$db_class[is.na(spec_df$db_class)] <- "Novel"
    spec_df <- spec_df[order(spec_df$score), ]
    spec_df$TF <- factor(spec_df$TF, levels = spec_df$TF)
    spec_df$Direction <- ifelse(spec_df$score > 0, "Activated", "Repressed")

    plots[[paste0("G_specific_", gt)]] <- ggplot2::ggplot(
      spec_df, ggplot2::aes(x = score, y = TF, fill = db_class)
    ) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
      ggplot2::scale_fill_manual(values = DB_COLORS, name = NULL) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(face = "italic", size = 9),
        plot.title  = ggplot2::element_text(face = "bold", size = 14),
        legend.position = "bottom"
      ) +
      ggplot2::labs(
        title = paste0("Genotype-Specific TFs: ", gt,
                       " (", nrow(spec_df), " TFs)"),
        subtitle = paste0("Significant in ", gt,
                          " but NOT in any other genotype (p<", p_thresh, ")"),
        x = "ULM activity score", y = NULL
      )
  }

  # ════════════════════════════════════════════════════════════════════
  # SAVE INDIVIDUALLY
  # ════════════════════════════════════════════════════════════════════
  if (!is.null(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    for (nm in names(plots)) {
      w <- if (grepl("^E_", nm)) 12
      else if (grepl("combined", nm)) min(length(d_plots) * 10, 30)
      else 10
      h <- if (grepl("^E_", nm)) max(8, length(focus_tfs) * 0.3 + 2)
      else if (grepl("^A_", nm)) 6
      else if (grepl("combined", nm)) 10
      else 10
      tryCatch({
        ggplot2::ggsave(file.path(save_dir,
                                  paste0("tf_genotype_", nm, ".png")),
                        plots[[nm]], width = w, height = h, dpi = 180)
      }, error = function(e)
        message("  [skip] ", nm, " — ", e$message))
    }
    message("[plot_omni_tf_genotype] Saved ", length(plots),
            " plots to: ", save_dir)
  }

  # ════════════════════════════════════════════════════════════════════
  # PRINT & SUMMARY
  # ════════════════════════════════════════════════════════════════════
  message("\n[plot_omni_tf_genotype] Generated ", length(plots), " plots:")
  for (nm in names(plots)) {
    message("  - ", nm)
    tryCatch(print(plots[[nm]]), error = function(e) NULL)
  }

  n_sig_act <- sum(sig_df$direction == "Activated")
  n_sig_rep <- sum(sig_df$direction == "Repressed")
  n_novel   <- sum(sig_df$db_class  == "Novel")
  message(sprintf("\n  %d significant TFs (%d activated, %d repressed, %d novel)",
                  nrow(sig_df), n_sig_act, n_sig_rep, n_novel))

  return(invisible(plots))
}

# =============================================================================
# 5. CELL TYPE DECONVOLUTION
# =============================================================================

#' Run Cell Type Deconvolution
#'
#' Estimates cell type composition per sample using either:
#' - `immunedeconv` — bulk RNA-seq methods (xCell, CIBERSORT, TIMER, EPIC,
#'   QUANTISEQ). No single-cell reference required.
#' - `MuSiC` — requires a `SingleCellExperiment` reference. More accurate
#'   when a matched reference is available.
#'
#' @param obj A validated OmniSeq object with a fitted model.
#' @param method Character. `"immunedeconv"` or `"music"`. Default:
#'   `"immunedeconv"`.
#' @param deconv_method Character. For immunedeconv: `"xcell"`, `"epic"`,
#'   `"quantiseq"`, `"timer"`, or `"cibersort"`. Default: `"xcell"`.
#' @param sc_ref A `SingleCellExperiment` object. Required for MuSiC.
#' @param cell_type_col Character. colData column with cell type labels
#'   (MuSiC only).
#' @param sample_col Character. colData column with sample IDs (MuSiC only).
#' @param interactive Logical. Prompt user to choose method and algorithm.
#'
#' @return The OmniSeq object with `$deconv` populated (cell types x samples).
#' @export
run_omni_deconvolution <- function(obj,
                                   method        = c("immunedeconv", "music"),
                                   deconv_method = "xcell",
                                   sc_ref        = NULL,
                                   cell_type_col = NULL,
                                   sample_col    = NULL,
                                   interactive   = FALSE) {

  validate_omniseq(obj)
  method <- match.arg(method)

  if (interactive && base::interactive()) {
    ans <- utils::menu(
      c(paste0("immunedeconv — bulk methods ",
               "(xCell, CIBERSORT, EPIC, TIMER, QUANTISEQ)"),
        "MuSiC — single-cell reference required"),
      title = "Select deconvolution method:"
    )
    method <- if (ans == 1L) "immunedeconv" else "music"

    if (method == "immunedeconv") {
      ans2 <- utils::menu(
        c("xCell (broad immune + stromal, recommended)",
          "EPIC (tumour microenvironment)",
          "QUANTISEQ (immune fraction)",
          "TIMER (TCGA-trained)",
          "CIBERSORT (requires local binary)"),
        title = "Select immunedeconv algorithm:"
      )
      deconv_method <- c("xcell", "epic", "quantiseq",
                         "timer", "cibersort")[ans2]
    }
  }

  vsd <- omni_vst(obj)
  if (is.null(vsd))
    stop("VST not cached. Run run_omni_de() first.", call. = FALSE)

  mat <- SummarizedExperiment::assay(vsd)

  # ── immunedeconv ──────────────────────────────────────────────────────────
  if (method == "immunedeconv") {
    .require_pkg("immunedeconv")

    message("\n[run_omni_deconvolution] Running immunedeconv: ",
            deconv_method, " on ", ncol(mat), " samples...")

    if (deconv_method == "xcell") {
      if (!requireNamespace("xCell", quietly = TRUE))
        stop("Package 'xCell' is required for xcell deconvolution.\n",
             "Install with: devtools::install_github('dviraran/xCell')",
             call. = FALSE)
      suppressPackageStartupMessages(library(xCell))
    }

    result <- suppressMessages(
      immunedeconv::deconvolute(mat, method = deconv_method)
    )

    result_df           <- as.data.frame(result)
    rownames(result_df) <- result_df$cell_type
    result_df$cell_type <- NULL
    deconv_mat          <- as.matrix(result_df)

    message("[run_omni_deconvolution] Estimated ", nrow(deconv_mat),
            " cell types.")
  }

  # ── MuSiC ────────────────────────────────────────────────────────────────
  else if (method == "music") {
    .require_pkg("MuSiC")
    .require_pkg("SingleCellExperiment")

    if (is.null(sc_ref))
      stop("sc_ref (SingleCellExperiment) is required for MuSiC.",
           call. = FALSE)
    if (is.null(cell_type_col) || is.null(sample_col))
      stop("cell_type_col and sample_col must be specified for MuSiC.",
           call. = FALSE)

    message("\n[run_omni_deconvolution] Running MuSiC deconvolution...")

    bulk_mat  <- omni_counts(obj)
    music_res <- MuSiC::music_prop(
      bulk.mtx = bulk_mat,
      sc.sce   = sc_ref,
      clusters = cell_type_col,
      samples  = sample_col,
      verbose  = FALSE
    )

    deconv_mat <- t(music_res$Est.prop.weighted)
    message("[run_omni_deconvolution] Estimated ", nrow(deconv_mat),
            " cell types.")
  }

  obj$deconv               <- deconv_mat
  obj$params$deconv_method <- paste0(
    method,
    if (method == "immunedeconv") paste0(" (", deconv_method, ")")
  )
  obj$params$deconv_contrast <- obj$active_contrast

  message("[run_omni_deconvolution] Done. Access with obj$deconv")
  return(obj)
}


# =============================================================================
# 6. DECONVOLUTION HEATMAP
# =============================================================================

#' Plot Cell Type Fraction Heatmap
#'
#' Displays estimated cell type fractions across all samples as an annotated
#' heatmap. Samples are ordered by condition by default.
#'
#' @param obj A validated OmniSeq object with deconvolution results.
#' @param cluster_rows Logical. Cluster cell types. Default: `TRUE`.
#' @param cluster_cols Logical. Cluster samples. Default: `FALSE`.
#' @param color_high Color for high fraction. Default: `"#B2182B"`.
#' @param color_low Color for low fraction. Default: `"white"`.
#' @param title Optional custom title.
#' @param aggregate_main Logical. Aggregate subtypes to 13 main cell type
#'   categories. Recommended for xCell (39 subtypes). Default: `TRUE`.
#'
#' @return A pheatmap object (returned invisibly).
#' @export
plot_omni_deconv_heatmap <- function(obj,
                                     cluster_rows   = TRUE,
                                     cluster_cols   = FALSE,
                                     color_high     = "#B2182B",
                                     color_low      = "white",
                                     title          = NULL,
                                     aggregate_main = TRUE) {

  validate_omniseq(obj)
  .require_pkg("pheatmap")

  deconv <- obj$deconv
  if (is.null(deconv))
    stop("No deconvolution results. Run run_omni_deconvolution() first.",
         call. = FALSE)

  meta     <- omni_metadata(obj)
  cond_col <- obj$params$condition

  # ── Aggregate subtypes to main cell types ──────────────────────────────────
  if (aggregate_main) {
    main_type_map <- list(
      "B cell"        = c("B cell", "B-cell", "naive B", "memory B",
                          "plasma B", "Class-switched", "pro-B", "pre-B"),
      "CD4+ T cell"   = c("T cell CD4", "CD4\\+", "Th1", "Th2", "Tfh",
                          "T cell regulatory", "Treg"),
      "CD8+ T cell"   = c("T cell CD8", "CD8\\+", "cytotoxic"),
      "NK cell"       = c("NK cell", "natural killer", "NKT"),
      "Neutrophil"    = c("Neutrophil"),
      "Monocyte"      = c("Monocyte", "monocyte"),
      "Dendritic"     = c("dendritic", "pDC", "mDC", "DC"),
      "Macrophage"    = c("Macrophage", "macrophage", "M1", "M2"),
      "Mast cell"     = c("Mast cell", "mast"),
      "Eosinophil"    = c("Eosinophil"),
      "Stem / prog."  = c("Hematopoietic stem", "progenitor", "HSC",
                          "CMP", "GMP", "MEP", "CLP"),
      "Fibroblast"    = c("fibroblast", "stroma"),
      "Endothelial"   = c("Endothelial", "endothelial")
    )

    rn      <- rownames(deconv)
    used    <- character(0)
    agg_mat <- matrix(NA_real_,
                      nrow     = length(main_type_map),
                      ncol     = ncol(deconv),
                      dimnames = list(names(main_type_map),
                                      colnames(deconv)))

    for (main_nm in names(main_type_map)) {
      patterns <- main_type_map[[main_nm]]
      pat      <- paste(patterns, collapse = "|")
      matched  <- rn[grepl(pat, rn, ignore.case = TRUE)]
      if (length(matched) > 0) {
        agg_mat[main_nm, ] <- colMeans(
          deconv[matched, , drop = FALSE], na.rm = TRUE)
        used <- union(used, matched)
      }
    }
    keep_rows  <- rowSums(!is.na(agg_mat)) > 0
    # Cell types outside the immune map (e.g. a tissue reference used with
    # MuSiC) are kept as they are rather than dropped.
    unmatched  <- setdiff(rn, used)
    deconv_use <- rbind(agg_mat[keep_rows, , drop = FALSE],
                        as.matrix(deconv[unmatched, , drop = FALSE]))
    subtitle   <- paste0("Aggregated to ", sum(keep_rows),
                         " main cell types",
                         if (length(unmatched) > 0)
                           paste0(" + ", length(unmatched), " unmapped"))
  } else {
    deconv_use <- deconv
    subtitle   <- paste0(nrow(deconv), " cell types")
  }

  # Align sample order to metadata
  sample_order <- rownames(meta)[order(meta[[cond_col]])]
  sample_order <- intersect(sample_order, colnames(deconv_use))
  deconv_plot  <- deconv_use[, sample_order, drop = FALSE]
  meta_ordered <- meta[sample_order, , drop = FALSE]

  annot_cols <- intersect(c(cond_col, obj$params$covariates),
                          colnames(meta_ordered))
  df_annot   <- meta_ordered[, annot_cols, drop = FALSE]
  df_annot[] <- lapply(df_annot,
                       function(x) if (is.factor(x)) droplevels(x) else x)

  pal <- grDevices::colorRampPalette(c(color_low, color_high))(100)

  if (is.null(title))
    title <- paste0("Cell Type Fractions — ",
                    obj$params$deconv_method %||% "Deconvolution")

  # Drop rows that are all zero or all NA
  keep_rows   <- apply(deconv_plot, 1,
                       function(x) any(!is.na(x) & x > 0))
  deconv_plot <- deconv_plot[keep_rows, , drop = FALSE]

  if (nrow(deconv_plot) == 0)
    stop("No non-zero cell type fractions to plot.", call. = FALSE)

  p <- pheatmap::pheatmap(
    deconv_plot,
    annotation_col = df_annot,
    color          = pal,
    cluster_rows   = cluster_rows && nrow(deconv_plot) >= 2,
    cluster_cols   = cluster_cols && ncol(deconv_plot) >= 2,
    border_color   = NA,
    show_colnames  = TRUE,
    fontsize_row   = 10,
    main           = paste0(title, "\n", subtitle)
  )

  return(invisible(p))
}


# =============================================================================
# 7. DECONVOLUTION CORRELATION
# =============================================================================

#' Correlate Cell Type Fractions with a Gene or Pathway
#'
#' For each cell type, computes Pearson correlation against either a gene's
#' VST expression or an ssGSEA pathway score. Displays results as a dot plot
#' (x = R, y = cell type, size = -log10 FDR).
#'
#' @param obj A validated OmniSeq object with deconvolution results.
#' @param feature Character. Gene name or pathway name to correlate against.
#' @param feature_type Character. `"gene"` or `"pathway"`. Auto-detected if
#'   `NULL`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_deconv_corr <- function(obj,
                                  feature      = NULL,
                                  feature_type = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  deconv <- obj$deconv
  if (is.null(deconv))
    stop("No deconvolution results. Run run_omni_deconvolution() first.",
         call. = FALSE)

  if (is.null(feature))
    stop("Specify a gene or pathway name via the 'feature' argument.",
         call. = FALSE)

  if (is.null(feature_type)) {
    vsd    <- omni_vst(obj)
    in_vst <- !is.null(vsd) &&
      feature %in% rownames(SummarizedExperiment::assay(vsd))
    in_sg  <- !is.null(obj$advanced_ssgsea) &&
      any(grepl(feature, rownames(obj$advanced_ssgsea),
                ignore.case = TRUE))

    feature_type <- if (in_vst) "gene" else if (in_sg) "pathway" else
      stop("Feature '", feature,
           "' not found in VST matrix or ssGSEA scores.", call. = FALSE)

    message("[plot_omni_deconv_corr] Auto-detected feature type: ",
            feature_type)
  }

  if (feature_type == "gene") {
    vsd <- omni_vst(obj)
    if (is.null(vsd)) stop("VST not cached.", call. = FALSE)
    mat <- SummarizedExperiment::assay(vsd)
    if (!feature %in% rownames(mat))
      stop("Gene '", feature, "' not found in VST matrix.", call. = FALSE)
    feature_vec <- mat[feature, ]
    x_label     <- paste("VST Expression:", feature)

  } else {
    scores <- obj$advanced_ssgsea
    if (is.null(scores))
      stop("No ssGSEA scores. Run run_omni_ssgsea() first.", call. = FALSE)
    idx <- grep(feature, rownames(scores), ignore.case = TRUE)[1]
    if (is.na(idx))
      stop("Pathway '", feature, "' not found in ssGSEA scores.",
           call. = FALSE)
    feature_vec <- scores[idx, ]
    x_label     <- paste("ssGSEA Score:", rownames(scores)[idx])
  }

  common_samples <- intersect(colnames(deconv), names(feature_vec))
  if (length(common_samples) < 3)
    stop("Fewer than 3 common samples between deconv and feature.",
         call. = FALSE)

  message("[plot_omni_deconv_corr] Correlating ", nrow(deconv),
          " cell types across ", length(common_samples), " samples...")

  cor_res <- lapply(rownames(deconv), function(ct) {
    x <- as.numeric(deconv[ct, common_samples])
    y <- as.numeric(feature_vec[common_samples])

    if (stats::var(x, na.rm = TRUE) == 0) return(NULL)

    ct_res <- stats::cor.test(x, y, method = "pearson")
    data.frame(
      CellType = ct,
      R        = as.numeric(ct_res$estimate),
      P        = as.numeric(ct_res$p.value),
      stringsAsFactors = FALSE
    )
  })

  cor_df          <- do.call(rbind, Filter(Negate(is.null), cor_res))
  cor_df$P_adj    <- stats::p.adjust(cor_df$P, method = "BH")
  cor_df$LogP     <- -log10(cor_df$P_adj + 1e-10)
  cor_df          <- cor_df[order(cor_df$R), ]
  cor_df$CellType <- factor(cor_df$CellType, levels = cor_df$CellType)

  p <- ggplot2::ggplot(
    cor_df,
    ggplot2::aes(x = R, y = CellType, size = LogP, color = R)
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_gradient2(
      low = "#2980B9", mid = "grey80", high = "#E74C3C",
      midpoint = 0, name = "Pearson R"
    ) +
    ggplot2::scale_size_continuous(name = "-log10(FDR)", range = c(2, 8)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        linewidth = 0.4) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(
      title    = paste("Cell Type Correlation with", feature),
      subtitle = paste("Method:",
                       obj$params$deconv_method %||% "deconvolution"),
      x        = paste("Pearson R vs", x_label),
      y        = NULL
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 8. STRING NETWORK QUERY
# =============================================================================

#' Query STRING Database for Top DEG Interactions
#'
#' Takes the top N significant DEGs and queries the STRING protein-protein
#' interaction database via `STRINGdb`. To avoid re-downloading the large
#' STRING files on every session, set `cache_dir` to a persistent local path.
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param n Integer. Number of top DEGs to query. Default: `50`.
#' @param pval_thresh Adjusted p-value threshold. Default: `0.05`.
#' @param string_score Integer. Minimum STRING interaction score (0–1000).
#'   Default: `400` (medium confidence).
#' @param species_id Integer. NCBI taxonomy ID. Default: `9606` (human).
#'   Use `10090` for mouse.
#' @param cache_dir Character. Local directory for STRING cache files.
#'   Default: `"~/.omniseq_cache/string"`. Set this to avoid re-downloading
#'   the 70MB STRING files on every new R session.
#'
#' @return The OmniSeq object with `$string_network` populated.
#' @export
run_omni_string <- function(obj,
                            n            = 50,
                            pval_thresh  = 0.05,
                            string_score = 400,
                            species_id   = 9606,
                            cache_dir    = "~/.omniseq_cache/string") {

  validate_omniseq(obj)
  .require_pkg("STRINGdb")

  res <- omni_results(obj)
  if (is.null(res))
    stop("No DE results. Run run_omni_de() first.", call. = FALSE)

  res_df      <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)

  sig_df <- res_df[res_df$padj < pval_thresh &
                     !is.na(res_df$log2FoldChange), ]

  if (nrow(sig_df) == 0)
    stop("No significant genes found.", call. = FALSE)

  sig_df <- sig_df[order(sig_df$padj), ]
  top_df <- utils::head(sig_df, n)

  cache_dir <- path.expand(cache_dir)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    message("[run_omni_string] Created STRING cache directory: ", cache_dir)
  }

  message("\n[run_omni_string] Querying STRING for ", nrow(top_df),
          " genes (score threshold: ", string_score, ")...")
  message("  Cache: ", cache_dir)

  string_db <- STRINGdb::STRINGdb$new(
    version         = "11.5",
    species         = species_id,
    score_threshold = string_score,
    network_type    = "full",
    input_directory = cache_dir
  )

  mapped <- string_db$map(top_df, "Gene", removeUnmappedRows = TRUE)

  if (nrow(mapped) == 0)
    stop("No genes could be mapped to STRING IDs.", call. = FALSE)

  message("[run_omni_string] Mapped ", nrow(mapped), " / ",
          nrow(top_df), " genes to STRING.")

  interactions <- string_db$get_interactions(mapped$STRING_id)

  if (is.null(interactions) || nrow(interactions) == 0)
    stop("No interactions found at score threshold ", string_score,
         ". Try lowering string_score.", call. = FALSE)

  id_to_gene             <- stats::setNames(mapped$Gene, mapped$STRING_id)
  interactions$from_gene <- id_to_gene[interactions$from]
  interactions$to_gene   <- id_to_gene[interactions$to]
  interactions           <- interactions[
    !is.na(interactions$from_gene) & !is.na(interactions$to_gene), ]

  lfc_map               <- stats::setNames(sig_df$log2FoldChange,
                                           sig_df$Gene)
  mapped$log2FoldChange <- lfc_map[mapped$Gene]
  mapped$padj           <- sig_df$padj[match(mapped$Gene, sig_df$Gene)]

  # STRINGdb can map one symbol to several STRING proteins (e.g. CT45A
  # paralogs). Every protein ID was used above to collect edges; keep one
  # node per gene and one edge per gene pair (highest score) so the
  # gene-level graph has unique vertex names.
  mapped       <- mapped[!duplicated(mapped$Gene), ]
  interactions <- .omni_collapse_gene_edges(interactions)

  obj$string_network <- list(
    nodes         = mapped,
    interactions  = interactions,
    score_thresh  = string_score,
    n_genes       = nrow(mapped),
    cache_dir     = cache_dir
  )

  message("[run_omni_string] Found ", nrow(interactions),
          " interactions among ", nrow(mapped), " hub genes.")
  message("  Access with obj$string_network")

  return(obj)
}


# =============================================================================
# 9. STRING NETWORK PLOT
# =============================================================================

#' Plot STRING Protein Interaction Network
#'
#' Visualises the STRING PPI network with nodes colored by log2 fold change
#' and sized by interaction degree. Supports both static (ggraph) and
#' interactive (visNetwork) output.
#'
#' @param obj A validated OmniSeq object with STRING results.
#' @param plot_type Character. `"static"` (ggraph) or `"interactive"`
#'   (visNetwork). Default: `"static"`.
#' @param min_degree Integer. Minimum interactions for a node to be shown.
#'   Default: `1`.
#' @param label_hubs Logical. Label only hub genes (degree >= median).
#'   Default: `TRUE`.
#'
#' @return A plot object (returned invisibly).
#' @export
plot_omni_string_network <- function(obj,
                                     plot_type  = c("static", "interactive"),
                                     min_degree = 1,
                                     label_hubs = TRUE) {

  validate_omniseq(obj)
  plot_type <- match.arg(plot_type)

  net <- obj$string_network
  if (is.null(net))
    stop("No STRING network. Run run_omni_string() first.", call. = FALSE)

  .require_pkg("igraph")

  interactions <- net$interactions
  nodes        <- net$nodes

  # Normalise column names defensively
  if (!"combined_score" %in% colnames(interactions)) {
    score_col <- grep("score|weight", colnames(interactions),
                      ignore.case = TRUE, value = TRUE)
    if (length(score_col) > 0)
      interactions$combined_score <- interactions[[score_col[1]]]
    else
      interactions$combined_score <- 1
  }

  if (!"from_gene" %in% colnames(interactions) &&
      "from" %in% colnames(interactions))
    interactions$from_gene <- interactions$from
  if (!"to_gene" %in% colnames(interactions) &&
      "to" %in% colnames(interactions))
    interactions$to_gene <- interactions$to

  # Remove self-loops and NA edges
  interactions <- interactions[
    !is.na(interactions$from_gene) & !is.na(interactions$to_gene) &
      interactions$from_gene != interactions$to_gene, ]

  if (nrow(interactions) == 0)
    stop("No valid interactions to plot after filtering.", call. = FALSE)

  # Objects created before run_omni_string() collapsed multi-mapped genes
  # can still hold duplicate gene nodes/edges.
  nodes        <- nodes[!duplicated(nodes$Gene), ]
  interactions <- .omni_collapse_gene_edges(interactions)

  g <- igraph::graph_from_data_frame(
    d        = interactions[, c("from_gene", "to_gene", "combined_score")],
    directed = FALSE,
    vertices = nodes[, c("Gene", "log2FoldChange", "padj")]
  )

  deg  <- igraph::degree(g)
  keep <- names(deg)[deg >= min_degree]
  g    <- igraph::induced_subgraph(g, keep)

  if (igraph::vcount(g) == 0)
    stop("No nodes remaining after degree filter.", call. = FALSE)

  deg2       <- igraph::degree(g)
  hub_thresh <- stats::median(deg2)
  lfc_vals   <- igraph::V(g)$log2FoldChange
  lfc_vals[is.na(lfc_vals)] <- 0

  # Dynamic symmetric color scale
  lfc_abs_max <- max(abs(lfc_vals), na.rm = TRUE)
  if (lfc_abs_max < 0.1) lfc_abs_max <- 1
  lfc_limits  <- c(-lfc_abs_max, lfc_abs_max)

  # ── STATIC (ggraph) ──────────────────────────────────────────────────────
  if (plot_type == "static") {
    .require_pkg("ggraph")

    label_nodes <- if (label_hubs)
      ifelse(deg2 >= hub_thresh, igraph::V(g)$name, "")
    else
      igraph::V(g)$name

    p <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(
        ggplot2::aes(width = combined_score),
        color = "grey75", alpha = 0.6
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(color = lfc_vals, size = deg2)
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label = label_nodes),
        repel    = TRUE,
        size     = 3,
        fontface = "italic"
      ) +
      ggplot2::scale_color_gradient2(
        low      = "#2166AC",
        mid      = "grey90",
        high     = "#B2182B",
        midpoint = 0,
        limits   = lfc_limits,
        oob      = scales::squish,
        name     = sprintf("Log\u2082 FC\n(\u00b1%.1f)", lfc_abs_max)
      ) +
      ggplot2::scale_size_continuous(name = "Degree", range = c(2, 8)) +
      ggraph::scale_edge_width_continuous(
        name = "STRING score", range = c(0.3, 1.5)
      ) +
      ggraph::theme_graph(base_size = 12, base_family = "sans") +
      ggplot2::labs(
        title    = "STRING Protein Interaction Network",
        subtitle = paste0(
          "Contrast: ", obj$active_contrast,
          " | Score \u2265 ", net$score_thresh,
          " | ", igraph::vcount(g), " nodes, ",
          igraph::ecount(g), " edges"
        )
      )

    print(p)
    return(invisible(p))
  }

  # ── INTERACTIVE (visNetwork) ─────────────────────────────────────────────
  if (plot_type == "interactive") {
    .require_pkg("visNetwork")

    .lfc_to_color <- function(lfc) {
      if (is.na(lfc) || lfc == 0) return("#AAAAAA")
      intensity <- min(255L, as.integer(abs(lfc) / 3 * 200 + 55))
      if (lfc > 0)
        sprintf("#FF%02X%02X", 255L - intensity, 255L - intensity)
      else
        sprintf("#%02X%02XFF", 255L - intensity, 255L - intensity)
    }

    vis_nodes <- data.frame(
      id    = igraph::V(g)$name,
      label = igraph::V(g)$name,
      value = deg2,
      color = sapply(lfc_vals, .lfc_to_color),
      title = paste0(
        "<b>", igraph::V(g)$name, "</b><br>",
        "LFC: ",    round(lfc_vals, 2), "<br>",
        "Degree: ", deg2
      ),
      stringsAsFactors = FALSE
    )

    el           <- igraph::as_edgelist(g)
    edge_weights <- igraph::E(g)$combined_score %||%
      rep(1, igraph::ecount(g))

    vis_edges <- data.frame(
      from  = el[, 1],
      to    = el[, 2],
      value = edge_weights / 1000,
      color = "rgba(150,150,150,0.5)",
      stringsAsFactors = FALSE
    )

    p <- visNetwork::visNetwork(
      vis_nodes, vis_edges,
      main   = "STRING PPI Network",
      height = "650px",
      width  = "100%"
    ) |>
      visNetwork::visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1),
        nodesIdSelection = TRUE
      ) |>
      visNetwork::visPhysics(
        solver           = "forceAtlas2Based",
        forceAtlas2Based = list(gravitationalConstant = -80,
                                centralGravity        = 0.01)
      ) |>
      visNetwork::visEdges(smooth = FALSE) |>
      visNetwork::visLayout(randomSeed = 42)

    print(p)
    return(invisible(p))
  }
}


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Retrieve TF results for the active contrast
#' Falls back to obj$tf_results for backward compatibility
#' @keywords internal
.omni_get_tf <- function(obj, contrast = NULL) {
  key <- contrast %||% obj$active_contrast
  if (!is.null(obj$tf_cache) && !is.null(obj$tf_cache[[key]]))
    return(obj$tf_cache[[key]])
  obj$tf_results
}

#' Fetch TF regulons from DoRothEA or CollecTRI via decoupleR
#' @keywords internal
.omni_get_regulons <- function(database, organism, dorothea_levels) {

  .require_pkg("decoupleR")

  regs <- switch(database,
                 "collectri" = decoupleR::get_collectri(
                   organism         = organism,
                   split_complexes  = FALSE),
                 "dorothea"  = decoupleR::get_dorothea(
                   organism = organism,
                   levels   = dorothea_levels),
                 "both"      = {
                   r1 <- decoupleR::get_collectri(
                     organism        = organism,
                     split_complexes = FALSE)
                   r2 <- decoupleR::get_dorothea(
                     organism = organism,
                     levels   = dorothea_levels)
                   r1 <- r1[, intersect(c("source", "target", "mor"),
                                        colnames(r1)), drop = FALSE]
                   r2 <- r2[, intersect(c("source", "target", "mor"),
                                        colnames(r2)), drop = FALSE]
                   combined <- rbind(r1, r2)
                   dup_key  <- paste(combined$source, combined$target,
                                     sep = "___")
                   combined[!duplicated(dup_key), ]
                 },
                 stop("Unknown database: ", database, call. = FALSE)
  )

  if (!"mor" %in% colnames(regs)) regs$mor <- 1
  regs <- regs[, c("source", "target", "mor")]
  regs <- regs[!is.na(regs$source) & !is.na(regs$target), ]

  return(regs)
}


#' Collapse STRING edges to one per unordered gene pair (max combined_score)
#' @keywords internal
.omni_collapse_gene_edges <- function(interactions) {
  if (nrow(interactions) == 0) return(interactions)
  a   <- pmin(interactions$from_gene, interactions$to_gene)
  b   <- pmax(interactions$from_gene, interactions$to_gene)
  ord <- order(-interactions$combined_score)
  interactions <- interactions[ord, ]
  key <- paste(a[ord], b[ord], sep = "___")
  interactions[!duplicated(key), ]
}
