# =============================================================================
# OmniSeq — 02_QC.R
#
# Fixes vs FlexiSeq:
#   - VST computed once and cached in model_info$vst (never recomputed)
#   - interactive() guard on all menu() calls
#   - QC stats stored in a clean, documented list structure
#   - distance matrix calculation is optional and clearly separated
# =============================================================================


#' Run Quality Control
#'
#' Calculates per-sample library sizes, detected feature counts, and
#' optionally sample-to-sample Euclidean distances on VST-normalised data.
#' Results are stored in `obj$qc` and the VST matrix is cached in
#' `obj$model_info$vst` for reuse by all downstream functions.
#'
#' @param obj A validated OmniSeq object.
#' @param calc_distances Logical. Whether to compute the sample-to-sample
#'   Euclidean distance matrix (used by the distance heatmap). Slightly
#'   slow on very large datasets. Default: `TRUE`.
#' @param interactive Logical. If `TRUE` and the session is interactive,
#'   prompts the user to confirm distance calculation. Default: `FALSE`.
#'
#' @return The OmniSeq object with `$qc` and `$model_info$vst` populated.
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- run_omni_qc(obj)
#' obj <- run_omni_qc(obj, calc_distances = FALSE)  # faster, skips heatmap data
#' }
run_omni_qc <- function(obj,
                        calc_distances = TRUE,
                        interactive    = FALSE) {

  validate_omniseq(obj)
  .require_pkg("DESeq2")
  .require_pkg("SummarizedExperiment")

  # ------------------------------------------------------------------
  # INTERACTIVE MODE
  # ------------------------------------------------------------------
  if (interactive && base::interactive()) {
    ans <- utils::menu(
      choices = c("Yes — calculate distances (needed for distance heatmap)",
                  "No  — skip distances (faster)"),
      title   = "Calculate sample-to-sample distance matrix?"
    )
    calc_distances <- (ans == 1L)
  }

  message("\n[run_omni_qc] Starting QC calculations...")

  counts_mat <- omni_counts(obj)
  meta       <- omni_metadata(obj)

  # ------------------------------------------------------------------
  # 1. BASIC STATS
  # ------------------------------------------------------------------
  message("[run_omni_qc] Computing library sizes and detected features...")

  lib_sizes      <- colSums(counts_mat)
  detected_feats <- colSums(counts_mat > 0)

  qc_stats <- data.frame(
    Sample         = colnames(counts_mat),
    Library_Size   = lib_sizes,
    Detected_Genes = detected_feats,
    Condition      = meta[[obj$params$condition]],
    row.names      = colnames(counts_mat),
    stringsAsFactors = FALSE
  )

  # ------------------------------------------------------------------
  # 2. VST NORMALISATION
  #    Use cached VST if already available (e.g. model was run first).
  #    Otherwise build a minimal intercept-only DESeqDataSet just for VST.
  #    Result is stored in model_info$vst — ONE cache for the whole session.
  # ------------------------------------------------------------------
  if (!is.null(obj$model_info$vst)) {
    message("[run_omni_qc] Using cached VST from model_info$vst.")
    vsd <- obj$model_info$vst

  } else if (!is.null(obj$model_info$dds)) {
    message("[run_omni_qc] Extracting VST from fitted DESeq2 model...")
    vsd <- .omni_vst(obj$model_info$dds, blind = FALSE)
    obj$model_info$vst <- vsd

  } else {
    message("[run_omni_qc] No model fitted yet — computing blind VST for QC only...")
    dds_qc <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts_mat,
      colData   = meta,
      design    = ~ 1
    )
    vsd <- .omni_vst(dds_qc, blind = TRUE)
    obj$model_info$vst <- vsd
  }

  # ------------------------------------------------------------------
  # 3. MEAN-VARIANCE STATS (for mean-variance trend plot)
  # ------------------------------------------------------------------
  vst_mat  <- SummarizedExperiment::assay(vsd)
  mv_stats <- data.frame(
    Mean = rowMeans(vst_mat),
    SD   = apply(vst_mat, 1, stats::sd)
  )

  # ------------------------------------------------------------------
  # 4. OPTIONAL: SAMPLE DISTANCES
  # ------------------------------------------------------------------
  dist_matrix <- NULL
  if (calc_distances) {
    message("[run_omni_qc] Computing Euclidean sample-to-sample distances...")
    dist_matrix <- as.matrix(stats::dist(t(vst_mat)))
  }

  # ------------------------------------------------------------------
  # 5. STORE IN CLEAN QC SLOT
  # ------------------------------------------------------------------
  obj$qc <- list(
    stats      = qc_stats,    # data.frame: Sample, Library_Size, Detected_Genes, Condition
    mv_stats   = mv_stats,    # data.frame: Mean, SD (VST-normalised, per feature)
    distances  = dist_matrix  # matrix or NULL
  )

  message("[run_omni_qc] QC complete. Access results with omni_qc(obj).")
  return(obj)
}
