# =============================================================================
# OmniSeq — 03_DE_engine.R
#
# Fixes vs previous version:
#   - lfcShrink() applied in the MAIN engine (not just custom contrasts)
#   - Results stored in results_cache[["contrast_name"]], never overwritten
#   - active_contrast pointer updated so accessors always return the right result
#   - VST cached in model_info$vst after fitting — never recomputed downstream
#   - fold_changes stored per-contrast inside results_cache entry, not in params
#   - interactive() guard on coefficient selection
#   - "Control" reference logic removed from here — handled at init
#   - Shrinkage method is configurable and recorded in model_info$lfc_shrinkage
#
# FIXES applied in this version:
#   - stat column preserved after ashr/apeglm shrinkage (both engines)
#     decoupleR ULM uses Wald t-statistics — losing stat forces fallback to LFC
#   - stat stored in results_cache entry so run_omni_mask() propagates it
#     to masked contrasts correctly
#   - Rebuilt DDS cached per column in model_info$cached_column_dds so that
#     multiple genotype contrasts (same column) reuse the same model fit
#     instead of re-running DESeq() for each contrast (major speedup)
#   - Enrichment slots NOT wiped in run_omni_contrast()
# =============================================================================


#' Run Differential Expression Analysis
#'
#' Fits a DESeq2 model using the design stored in the OmniSeq object, extracts
#' results for a chosen coefficient, applies log-fold-change shrinkage, and
#' caches both the results and the VST matrix for downstream use.
#'
#' Results are stored in `obj$results_cache[["<contrast_name>"]]` so that
#' multiple contrasts can coexist. Use [omni_results()] to retrieve them and
#' [omni_contrasts()] to list all available contrasts.
#'
#' @param obj A validated OmniSeq object produced by [init_omniseq()].
#' @param shrinkage Character. LFC shrinkage method. One of `"apeglm"`
#'   (recommended, default), `"ashr"`, or `"none"`. `"apeglm"` requires the
#'   `apeglm` package; if unavailable it falls back to `"ashr"` automatically.
#' @param coef Integer or character. Which DESeq2 result coefficient to extract.
#'   If `NULL` (default), the last coefficient in the model (the main condition
#'   effect) is used automatically. In interactive mode, the user is prompted.
#' @param interactive Logical. If `TRUE` and the session is interactive,
#'   prompts the user to select a coefficient from the fitted model.
#'   Default: `FALSE`.
#'
#' @return The OmniSeq object with the following updated:
#'   - `$model_info$dds` — the fitted DESeqDataSet
#'   - `$model_info$vst` — cached VST (blind = FALSE)
#'   - `$model_info$extracted_coef` — name of the active coefficient
#'   - `$model_info$lfc_shrinkage` — shrinkage method used
#'   - `$results_cache[["<coef_name>"]]` — a list with `$results` and `$fold_changes`
#'   - `$active_contrast` — set to the extracted coefficient name
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- run_omni_de(obj)
#' obj <- run_omni_de(obj, shrinkage = "ashr")
#' obj <- run_omni_de(obj, interactive = TRUE)  # prompts for coefficient
#'
#' # Retrieve results
#' res <- omni_results(obj)
#' head(as.data.frame(res))
#' }
run_omni_de <- function(obj,
                        shrinkage   = c("apeglm", "ashr", "none"),
                        coef        = NULL,
                        interactive = FALSE) {

  validate_omniseq(obj)
  .require_pkg("DESeq2")
  .require_pkg("SummarizedExperiment")

  shrinkage <- match.arg(shrinkage)

  message("\n", strrep("=", 52))
  message("  OmniSeq: DESeq2 Differential Expression Engine")
  message(strrep("=", 52))

  counts_mat <- omni_counts(obj)
  meta       <- omni_metadata(obj)
  formula    <- obj$params$formula

  message("[run_omni_de] Design formula: ", deparse(formula))

  # ------------------------------------------------------------------
  # 1. BUILD DESEQ2 DATASET
  # ------------------------------------------------------------------
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = meta,
    design    = formula
  )

  # ------------------------------------------------------------------
  # 2. FIT THE MODEL
  # ------------------------------------------------------------------
  message("[run_omni_de] Fitting DESeq2 model (this may take a moment)...")
  dds <- DESeq2::DESeq(dds)

  # ------------------------------------------------------------------
  # 3. CACHE VST IMMEDIATELY (blind = FALSE uses the fitted dispersions)
  # ------------------------------------------------------------------
  message("[run_omni_de] Caching VST normalised matrix...")
  vsd <- DESeq2::vst(dds, blind = FALSE)
  obj$model_info$dds <- dds
  obj$model_info$vst <- vsd

  # ------------------------------------------------------------------
  # 4. COEFFICIENT SELECTION
  # ------------------------------------------------------------------
  res_names   <- DESeq2::resultsNames(dds)
  target_coef <- res_names[length(res_names)]

  if (interactive && base::interactive()) {
    cat("\nAvailable model coefficients:\n")
    for (i in seq_along(res_names))
      cat("  ", i, ": ", res_names[i], "\n", sep = "")

    ans <- readline(prompt = paste0(
      "Select coefficient [default = ", length(res_names), "]: "
    ))

    if (nzchar(trimws(ans))) {
      idx <- suppressWarnings(as.integer(ans))
      if (is.na(idx) || idx < 1L || idx > length(res_names)) {
        warning("Invalid selection — using default coefficient: ",
                target_coef, call. = FALSE)
      } else {
        target_coef <- res_names[idx]
      }
    }

  } else if (!is.null(coef)) {
    if (is.numeric(coef)) {
      coef <- as.integer(coef)
      if (coef < 1L || coef > length(res_names))
        stop("coef index ", coef, " is out of range (1 to ",
             length(res_names), ").", call. = FALSE)
      target_coef <- res_names[coef]
    } else if (is.character(coef)) {
      if (!coef %in% res_names)
        stop("coef '", coef, "' not found in model. Available: ",
             paste(res_names, collapse = ", "), call. = FALSE)
      target_coef <- coef
    }
  }

  message("[run_omni_de] Extracting results for: ", target_coef)

  # ------------------------------------------------------------------
  # 5. EXTRACT RESULTS & PRESERVE STAT COLUMN
  # ------------------------------------------------------------------
  res       <- DESeq2::results(dds, name = target_coef)
  orig_stat <- as.numeric(res$stat)
  names(orig_stat) <- rownames(res)

  # ------------------------------------------------------------------
  # 6. LFC SHRINKAGE
  # ------------------------------------------------------------------
  shrinkage_used <- "none"

  if (shrinkage != "none") {
    apeglm_ok <- requireNamespace("apeglm", quietly = TRUE)
    ashr_ok   <- requireNamespace("ashr",   quietly = TRUE)

    if (shrinkage == "apeglm" && !apeglm_ok) {
      if (ashr_ok) {
        warning("Package 'apeglm' not installed — falling back to 'ashr'.",
                call. = FALSE)
        shrinkage <- "ashr"
      } else {
        warning("Neither 'apeglm' nor 'ashr' installed — skipping.",
                call. = FALSE)
        shrinkage <- "none"
      }
    }

    if (shrinkage == "apeglm" && apeglm_ok) {
      message("[run_omni_de] Applying apeglm LFC shrinkage...")
      res <- suppressMessages(
        DESeq2::lfcShrink(dds, coef = target_coef, type = "apeglm",
                          quiet = TRUE)
      )
      shrinkage_used <- "apeglm"

    } else if (shrinkage == "ashr" && ashr_ok) {
      message("[run_omni_de] Applying ashr LFC shrinkage...")
      res <- suppressMessages(
        DESeq2::lfcShrink(dds, coef = target_coef, type = "ashr",
                          quiet = TRUE)
      )
      shrinkage_used <- "ashr"
    }
  }

  # Restore stat column if shrinkage dropped it
  if (is.null(res$stat) || all(is.na(res$stat))) {
    res$stat <- orig_stat[rownames(res)]
  }

  obj$model_info$lfc_shrinkage  <- shrinkage_used
  obj$model_info$extracted_coef <- target_coef

  # ------------------------------------------------------------------
  # 7. STORE IN results_cache
  # ------------------------------------------------------------------
  fc        <- res$log2FoldChange
  names(fc) <- rownames(res)
  fc        <- sort(fc[!is.na(fc)], decreasing = TRUE)

  stat_vec        <- as.numeric(res$stat)
  names(stat_vec) <- rownames(res)
  stat_vec        <- stat_vec[!is.na(stat_vec)]

  obj$results_cache[[target_coef]] <- list(
    results      = res,
    fold_changes = fc,
    stat         = stat_vec,
    shrinkage    = shrinkage_used,
    timestamp    = Sys.time()
  )

  obj$active_contrast <- target_coef

  # ------------------------------------------------------------------
  # 8. SUMMARY
  # ------------------------------------------------------------------
  res_df <- as.data.frame(res)
  n_sig  <- sum(!is.na(res_df$padj) & res_df$padj < 0.05, na.rm = TRUE)
  n_up   <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 &
                  res_df$log2FoldChange > 1,  na.rm = TRUE)
  n_down <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 &
                  res_df$log2FoldChange < -1, na.rm = TRUE)

  message("\n[run_omni_de] DE complete!")
  message("  Contrast      : ", target_coef)
  message("  Shrinkage     : ", shrinkage_used)
  message("  Significant   : ", n_sig,  " (padj < 0.05)")
  message("  Up (|LFC|>1)  : ", n_up)
  message("  Down (|LFC|>1): ", n_down)
  message("\n  Tip: use omni_results(obj) to access results,")
  message("       omni_contrasts(obj) to list all contrasts.")

  return(obj)
}


# =============================================================================
# CUSTOM CONTRAST ENGINE
# =============================================================================

#' Extract a Pairwise Contrast from a Fitted OmniSeq Model
#'
#' Extracts a specific pairwise comparison (e.g. GroupA vs GroupB) from the
#' already-fitted DESeq2 model without refitting. If the requested condition
#' column differs from the one used to build the current model, the model is
#' automatically rebuilt with the new column as the grouping variable.
#'
#' The rebuilt model is cached in `model_info$cached_column_dds[[key]]` so
#' that multiple contrasts using the same column (e.g. UNC93B1 vs Ctrl,
#' SPEN vs Ctrl, TLR7 vs Ctrl — all using the Genotype column) reuse the
#' same fitted model. This avoids redundant DESeq() calls.
#'
#' Results are added to `results_cache` under the key `"GroupA vs GroupB"`,
#' leaving all previous contrasts intact.
#'
#' @param obj A validated OmniSeq object with a fitted model.
#' @param column Character. The metadata column containing the groups to
#'   compare. Defaults to the primary condition set at init.
#' @param group1 Character. The numerator group (e.g. `"Treated"`).
#' @param group2 Character. The denominator / baseline group (e.g. `"Control"`).
#' @param shrinkage Character. LFC shrinkage method: `"apeglm"`, `"ashr"`,
#'   or `"none"`. Defaults to `"ashr"` (always works for contrast-based
#'   extraction; apeglm requires a named coefficient).
#' @param new_covariates Optional character vector. If the model needs to be
#'   rebuilt (column differs from current condition), use these covariates.
#'   Defaults to the covariates already stored in the object.
#'
#' @return The OmniSeq object with a new entry in `$results_cache` and
#'   `$active_contrast` updated.
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- run_omni_contrast(obj, group1 = "Treated", group2 = "Control")
#' obj <- run_omni_contrast(obj, column = "Genotype",
#'                          group1 = "KO", group2 = "WT",
#'                          shrinkage = "ashr")
#' }
run_omni_contrast <- function(obj,
                              column         = NULL,
                              group1         = NULL,
                              group2         = NULL,
                              shrinkage      = c("ashr", "apeglm", "none"),
                              new_covariates = NULL) {

  validate_omniseq(obj)
  .require_pkg("DESeq2")

  shrinkage <- match.arg(shrinkage)
  column    <- column %||% obj$params$condition
  meta      <- omni_metadata(obj)

  if (!column %in% colnames(meta))
    stop("Column '", column, "' not found in metadata.", call. = FALSE)

  available_groups <- unique(as.character(meta[[column]]))
  available_groups <- available_groups[!is.na(available_groups)]

  if (is.null(group1) || is.null(group2))
    stop("Both group1 and group2 must be specified.\n",
         "Available groups in '", column, "': ",
         paste(available_groups, collapse = ", "), call. = FALSE)

  for (g in c(group1, group2)) {
    if (!g %in% available_groups)
      stop("Group '", g, "' not found in column '", column, "'.\n",
           "Available: ", paste(available_groups, collapse = ", "),
           call. = FALSE)
  }

  contrast_key <- paste(group1, "vs", group2)
  message("\n[run_omni_contrast] Extracting contrast: ", contrast_key)

  # ------------------------------------------------------------------
  # DETERMINE WHICH DDS TO USE
  # If column == primary condition, use the main model (already fitted).
  # If column differs, we need a model fitted with that column.
  # Cache rebuilt models in model_info$cached_column_dds so that
  # multiple contrasts using the same column reuse the same fit.
  # ------------------------------------------------------------------
  active_dds <- NULL

  if (column == obj$params$condition && !is.null(obj$model_info$dds)) {
    # Same column as the main model — reuse it directly
    active_dds <- obj$model_info$dds
    message("[run_omni_contrast] Reusing main model (column: ", column, ")")

  } else {
    # Different column — check cache first, rebuild if needed
    covariates     <- new_covariates %||% obj$params$covariates
    covariates     <- setdiff(covariates, column)
    design_formula <- .build_formula(covariates, column)
    cache_key      <- paste0(column, "|||", deparse(design_formula))

    if (!is.null(obj$model_info$cached_column_dds) &&
        !is.null(obj$model_info$cached_column_dds[[cache_key]])) {
      # Reuse cached model — no re-fit needed
      active_dds <- obj$model_info$cached_column_dds[[cache_key]]
      message("[run_omni_contrast] Reusing cached model for column: ",
              column)

    } else {
      # Rebuild from scratch
      message("[run_omni_contrast] Rebuilding DESeq2 model for column: ",
              column)
      message("[run_omni_contrast] Formula: ", deparse(design_formula))

      counts_mat      <- omni_counts(obj)
      meta[[column]]  <- as.factor(meta[[column]])

      suppressMessages({
        dds_new <- DESeq2::DESeqDataSetFromMatrix(
          countData = counts_mat,
          colData   = meta,
          design    = design_formula
        )
        dds_new <- DESeq2::DESeq(dds_new, quiet = TRUE)
      })

      # Cache the rebuilt model for reuse by subsequent contrasts
      if (is.null(obj$model_info$cached_column_dds))
        obj$model_info$cached_column_dds <- list()
      obj$model_info$cached_column_dds[[cache_key]] <- dds_new

      active_dds <- dds_new
      message("[run_omni_contrast] Model rebuilt and cached.")
    }
  }

  # ------------------------------------------------------------------
  # EXTRACT RESULTS & PRESERVE STAT COLUMN
  # ------------------------------------------------------------------
  res       <- DESeq2::results(active_dds,
                               contrast = c(column, group1, group2))
  orig_stat <- as.numeric(res$stat)
  names(orig_stat) <- rownames(res)

  # ------------------------------------------------------------------
  # LFC SHRINKAGE
  # ------------------------------------------------------------------
  shrinkage_used <- "none"

  if (shrinkage != "none") {
    ashr_ok   <- requireNamespace("ashr",   quietly = TRUE)
    apeglm_ok <- requireNamespace("apeglm", quietly = TRUE)

    # Try apeglm via named coefficient
    coef_name <- paste0(column, "_", group1, "_vs_", group2)

    if (shrinkage == "apeglm" &&
        coef_name %in% DESeq2::resultsNames(active_dds) &&
        apeglm_ok) {
      message("[run_omni_contrast] Applying apeglm shrinkage...")
      res <- suppressMessages(
        DESeq2::lfcShrink(active_dds, coef = coef_name, type = "apeglm",
                          res = res, quiet = TRUE)
      )
      shrinkage_used <- "apeglm"

    } else if (ashr_ok) {
      if (shrinkage == "apeglm")
        message("[run_omni_contrast] apeglm unavailable for this ",
                "contrast — using ashr.")
      else
        message("[run_omni_contrast] Applying ashr shrinkage...")

      res <- suppressMessages(
        DESeq2::lfcShrink(active_dds,
                          contrast = c(column, group1, group2),
                          type     = "ashr",
                          res      = res,
                          quiet    = TRUE)
      )
      shrinkage_used <- "ashr"

    } else {
      warning("No shrinkage package available. Returning unshrunk LFCs.",
              call. = FALSE)
    }
  }

  # Restore stat column if shrinkage dropped it
  if (is.null(res$stat) || all(is.na(res$stat))) {
    res$stat <- orig_stat[rownames(res)]
  }

  obj$model_info$lfc_shrinkage  <- shrinkage_used
  obj$model_info$extracted_coef <- contrast_key

  # ------------------------------------------------------------------
  # STORE IN results_cache — each contrast fully independent
  # stat stored alongside so masked contrasts and run_omni_tf() can
  # read the Wald t-statistics without accessing the DDS
  # ------------------------------------------------------------------
  fc        <- res$log2FoldChange
  names(fc) <- rownames(res)
  fc        <- sort(fc[!is.na(fc)], decreasing = TRUE)

  stat_vec        <- as.numeric(res$stat)
  names(stat_vec) <- rownames(res)
  stat_vec        <- stat_vec[!is.na(stat_vec)]

  obj$results_cache[[contrast_key]] <- list(
    results      = res,
    fold_changes = fc,
    stat         = stat_vec,
    shrinkage    = shrinkage_used,
    timestamp    = Sys.time()
  )

  obj$active_contrast <- contrast_key

  # NOTE: enrichment slots are NOT wiped here.

  res_df <- as.data.frame(res)
  n_sig  <- sum(!is.na(res_df$padj) & res_df$padj < 0.05, na.rm = TRUE)

  message("[run_omni_contrast] Done. Significant (padj<0.05): ", n_sig)
  message("  Stored under key: '", contrast_key, "'")
  message("  Shrinkage: ", shrinkage_used)
  message("  Stat vector: n=", length(stat_vec),
          " | median=", round(stats::median(stat_vec), 4),
          " | sum=", round(sum(stat_vec), 1))

  return(obj)
}


# =============================================================================
# GENE MASKING (post-DE)
# =============================================================================

#' Mask Specific Features from the Active Contrast
#'
#' Removes a list of features (genes, TEs) from the active contrast's results
#' and fold-change vector to reveal secondary biology in downstream plots.
#' The original results are preserved in `results_cache` under the original
#' key; the masked version is stored under a new key and becomes the new
#' `active_contrast`.
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param features Optional character vector of gene/feature names to mask.
#' @param signature Optional character string. Shorthand name of a built-in
#'   OmniSeq signature to mask. Options: `"isgs_28"`, `"isgs_extended"`,
#'   `"ribosomal"`, `"mitochondrial"`, `"immunoglobulin"`. Run
#'   `omni_list_signatures()` for details.
#' @param mask_label A short label describing what was masked. Used to name
#'   the new cache entry. Auto-set from `signature` if not provided.
#'
#' @return The OmniSeq object with a new masked contrast in `results_cache`
#'   and `active_contrast` updated to the masked key.
#'
#' @details
#' Features can be provided in three ways:
#' 1. A character vector of gene names via `features`
#' 2. A built-in signature shorthand via `signature` (e.g. `"isgs_28"`)
#' 3. Both combined — `signature` genes are merged with `features`
#'
#' Run `omni_list_signatures()` to see all available built-in signatures.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Option 1: use a built-in signature by shorthand
#' obj <- run_omni_mask(obj, signature = "isgs_28")
#'
#' # Option 2: use the built-in vector directly
#' obj <- run_omni_mask(obj, features = omniseq_isgs_28, mask_label = "ISGs")
#'
#' # Option 3: combine built-in + custom genes
#' obj <- run_omni_mask(obj,
#'                      features  = c("EXTRA_GENE1", "EXTRA_GENE2"),
#'                      signature = "isgs_28",
#'                      mask_label = "ISGs_Extended")
#'
#' # The original contrast is still accessible
#' omni_contrasts(obj)
#' omni_results(obj, contrast = "Condition_Patient_vs_Control")
#' omni_results(obj, contrast = "Condition_Patient_vs_Control (excl. ISGs_28)")
#' }
run_omni_mask <- function(obj,
                          features   = NULL,
                          signature  = NULL,
                          mask_label = NULL) {

  validate_omniseq(obj)

  # ------------------------------------------------------------------
  # Resolve features from signature shorthand and/or direct vector
  # ------------------------------------------------------------------
  sig_genes <- character(0)

  if (!is.null(signature)) {
    sig_genes <- .omni_get_signature(signature)
    if (is.null(mask_label))
      mask_label <- paste0(
        toupper(substring(signature, 1, 1)),
        substring(gsub("_", " ", signature), 2)
      )
    message("[run_omni_mask] Built-in signature '", signature,
            "': ", length(sig_genes), " genes loaded.")
  }

  all_features <- unique(c(sig_genes, features))

  if (length(all_features) == 0)
    stop("No features to mask. Provide 'features', 'signature', or both.",
         call. = FALSE)

  if (is.null(mask_label)) mask_label <- "Masked"

  # ------------------------------------------------------------------
  # Validate active contrast
  # ------------------------------------------------------------------
  if (is.null(obj$active_contrast))
    stop("No active contrast found. Run run_omni_de() first.",
         call. = FALSE)

  cache_entry <- obj$results_cache[[obj$active_contrast]]
  if (is.null(cache_entry))
    stop("Active contrast '", obj$active_contrast,
         "' not found in results_cache.", call. = FALSE)

  res    <- cache_entry$results
  res_df <- as.data.frame(res)
  hits   <- intersect(all_features, rownames(res_df))

  if (length(hits) == 0L) {
    message("[run_omni_mask] None of the specified features were found. ",
            "Returning object unchanged.")
    return(obj)
  }

  message("[run_omni_mask] Masking ", length(hits), " / ",
          length(all_features), " features (label: '", mask_label, "')")

  if (length(all_features) - length(hits) > 0)
    message("[run_omni_mask] Note: ", length(all_features) - length(hits),
            " feature(s) not found and were skipped.")

  # ------------------------------------------------------------------
  # Build masked results — originals untouched
  # ------------------------------------------------------------------
  keep       <- !rownames(res_df) %in% hits
  masked_res <- res[keep, ]

  masked_fc <- cache_entry$fold_changes
  masked_fc <- masked_fc[!names(masked_fc) %in% hits]

  # Propagate stat vector
  masked_stat <- cache_entry$stat
  if (!is.null(masked_stat))
    masked_stat <- masked_stat[!names(masked_stat) %in% hits]

  new_key <- paste0(obj$active_contrast, " (excl. ", mask_label, ")")

  obj$results_cache[[new_key]] <- list(
    results      = masked_res,
    fold_changes = masked_fc,
    stat         = masked_stat,
    shrinkage    = cache_entry$shrinkage,
    timestamp    = Sys.time(),
    masked       = hits,
    mask_label   = mask_label
  )

  obj$active_contrast <- new_key

  # Wipe enrichment — it belongs to the unmasked contrast
  obj$enrichment      <- NULL
  obj$advanced_gsea   <- NULL
  obj$pathway_results <- NULL

  message("[run_omni_mask] Active contrast is now: '", new_key, "'")
  message("[run_omni_mask] Original preserved under: '",
          gsub(" \\(excl.*", "", new_key), "'")
  message("  Tip: use omni_contrasts(obj) to see all available contrasts.")

  return(obj)
}


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Retrieve fold changes for the active contrast
#' @keywords internal
.omni_fold_changes <- function(obj, contrast = NULL) {
  key   <- contrast %||% obj$active_contrast
  entry <- obj$results_cache[[key]]
  if (is.null(entry)) return(NULL)
  entry$fold_changes
}

#' Retrieve Wald stat vector for the active contrast
#' Falls back to LFC / lfcSE if stat was not cached
#' @keywords internal
.omni_stat_vector <- function(obj, contrast = NULL) {
  key   <- contrast %||% obj$active_contrast
  entry <- obj$results_cache[[key]]
  if (is.null(entry)) return(NULL)

  # Prefer cached stat (saved before shrinkage)
  if (!is.null(entry$stat) && length(entry$stat) > 0)
    return(entry$stat)

  # Fallback: compute from results
  res <- as.data.frame(entry$results)
  if ("stat" %in% colnames(res)) {
    sv        <- res$stat
    names(sv) <- rownames(res)
    sv        <- sv[!is.na(sv)]
    if (length(sv) > 0) return(sv)
  }

  # Last resort: LFC / lfcSE approximation
  if (all(c("log2FoldChange","lfcSE") %in% colnames(res))) {
    message("[.omni_stat_vector] stat not cached for '",
            key, "'. Using LFC/lfcSE approximation.")
    sv        <- res$log2FoldChange / res$lfcSE
    names(sv) <- rownames(res)
    return(sv[!is.na(sv)])
  }

  NULL
}
