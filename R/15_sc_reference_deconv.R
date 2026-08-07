# =============================================================================
# OmniSeq — 15_sc_reference_deconv.R
#
# Single-cell / pseudo-bulk reference-based deconvolution.
# Implements Emilia's two-layer architecture:
#   Layer 1: Cell identity from healthy controls → composition estimates
#   Layer 2: Disease activation states from untreated disease → activation scores
#
# AUTO-BRANCHING LOGIC:
#   Raw scRNA-seq (>1000 cells) → MuSiC deconvolution
#   Pseudo-bulk (<500 profiles)  → Signature scoring + optional BisqueRNA
#
# SECTION INDEX:
#   SECTION 1:  run_omni_sc_deconv()         — auto-branching entry point
#   SECTION 2:  .omni_sc_signature_score()   — pseudo-bulk signature scoring
#   SECTION 3:  .omni_sc_music()             — raw SCE → MuSiC
#   SECTION 4:  .omni_sc_bisque()            — pseudo-bulk → BisqueRNA (NNLS)
#   SECTION 5:  run_omni_disease_state()     — disease activation scoring
#   SECTION 6:  plot_omni_sc_heatmap()       — cell type score heatmap
#   SECTION 7:  plot_omni_sc_genotype()      — per-genotype boxplots + signif
#   SECTION 8:  plot_omni_sc_correlation()   — cell type vs feature bubble
#   SECTION 9:  plot_omni_sc_profile()       — genotype immune profile lines
#   SECTION 10: plot_omni_sc_two_layer()     — composition vs activation figure
#   SECTION 11: plot_omni_sc_benchmark()     — compare xCell vs signature
#   SECTION 12: Internal helpers
# =============================================================================


# =============================================================================
# SECTION 1: AUTO-BRANCHING ENTRY POINT
# =============================================================================

#' Reference-Based Cell Type Deconvolution (Auto-Branching)
#'
#' Automatically detects whether the reference is raw single-cell data
#' or processed pseudo-bulk, and routes to the appropriate method.
#'
#' @param obj A validated OmniSeq object with VST cached.
#' @param sc_ref A SingleCellExperiment object.
#' @param cell_type_col Character. Cell type column. Default: `"cell_type"`.
#' @param sample_col Character. Donor/sample column. Default: `"donor"`.
#' @param group_col Character. Condition column. Default: `"group"`.
#' @param ref_group Character. Reference group. Default: `"Control"`.
#' @param method Character. `"auto"`, `"signature"`, `"music"`, or
#'   `"bisque"`. Default: `"auto"`.
#' @param n_markers Integer. Markers per cell type. Default: `50`.
#' @param min_cpm Numeric. Min CPM for markers. Default: `10`.
#' @param ref_label Character. Label for plots. Default: `"scRNA-seq reference"`.
#' @param run_bisque Logical. Also run BisqueRNA. Default: `FALSE`.
#'
#' @return The OmniSeq object with `$sc_deconv` populated.
#' @export
run_omni_sc_deconv <- function(obj,
                               sc_ref,
                               cell_type_col = "cell_type",
                               sample_col    = "donor",
                               group_col     = "group",
                               ref_group     = "Control",
                               method        = c("auto", "signature",
                                                 "music", "bisque"),
                               n_markers     = 50L,
                               min_cpm       = 10,
                               ref_label     = "scRNA-seq reference",
                               run_bisque    = FALSE) {

  validate_omniseq(obj)
  .require_pkg("SingleCellExperiment")
  .require_pkg("SummarizedExperiment")

  method <- match.arg(method)

  if (!inherits(sc_ref, "SingleCellExperiment"))
    stop("sc_ref must be a SingleCellExperiment object.", call. = FALSE)

  if (!"counts" %in% SummarizedExperiment::assayNames(sc_ref))
    stop("sc_ref must have a 'counts' assay.", call. = FALSE)

  # Validate columns
  cd <- as.data.frame(SummarizedExperiment::colData(sc_ref))
  for (col in c(cell_type_col, sample_col)) {
    if (!col %in% colnames(cd))
      stop("Column '", col, "' not in colData(sc_ref). Available: ",
           paste(colnames(cd), collapse = ", "), call. = FALSE)
  }

  if (!is.null(ref_group) && !is.null(group_col)) {
    if (!group_col %in% colnames(cd))
      stop("group_col '", group_col, "' not in colData.", call. = FALSE)
    if (!ref_group %in% cd[[group_col]])
      stop("ref_group '", ref_group, "' not found. Available: ",
           paste(unique(cd[[group_col]]), collapse = ", "), call. = FALSE)
  }

  n_cells    <- ncol(sc_ref)
  counts_mat <- SummarizedExperiment::assay(sc_ref, "counts")

  # ── Auto-detect: raw single-cell vs pseudo-bulk ────────────────────
  if (method == "auto") {
    max_count  <- max(counts_mat)
    median_lib <- stats::median(colSums(counts_mat))

    is_pseudobulk <- n_cells < 500 && max_count > 100000

    if (is_pseudobulk) {
      method <- "signature"
      message("[run_omni_sc_deconv] Auto-detected: PSEUDO-BULK (",
              n_cells, " profiles, max count = ",
              format(max_count, big.mark = ","), ")")
      message("  -> Routing to signature scoring method")
    } else {
      method <- "music"
      message("[run_omni_sc_deconv] Auto-detected: RAW SINGLE-CELL (",
              n_cells, " cells, median library = ",
              format(round(median_lib), big.mark = ","), ")")
      message("  -> Routing to MuSiC deconvolution")
    }
  }

  # ── Route to method ────────────────────────────────────────────────
  if (method == "signature") {
    obj <- .omni_sc_signature_score(
      obj = obj, sc_ref = sc_ref,
      cell_type_col = cell_type_col, sample_col = sample_col,
      group_col = group_col, ref_group = ref_group,
      n_markers = n_markers, min_cpm = min_cpm,
      ref_label = ref_label
    )

    if (run_bisque) {
      message("\n[run_omni_sc_deconv] Also running BisqueRNA (NNLS)...")
      obj <- tryCatch(
        .omni_sc_bisque(obj, sc_ref = sc_ref,
                        cell_type_col = cell_type_col,
                        sample_col = sample_col,
                        group_col = group_col,
                        ref_group = ref_group),
        error = function(e) {
          message("  BisqueRNA failed: ", e$message); obj })
    }

  } else if (method == "music") {
    obj <- .omni_sc_music(obj, sc_ref = sc_ref,
                          cell_type_col = cell_type_col,
                          sample_col = sample_col)

  } else if (method == "bisque") {
    obj <- .omni_sc_bisque(obj, sc_ref = sc_ref,
                           cell_type_col = cell_type_col,
                           sample_col = sample_col,
                           group_col = group_col,
                           ref_group = ref_group)
  }

  return(obj)
}


# =============================================================================
# SECTION 2: PSEUDO-BULK SIGNATURE SCORING
# =============================================================================

#' @keywords internal
.omni_sc_signature_score <- function(obj, sc_ref,
                                     cell_type_col, sample_col,
                                     group_col, ref_group,
                                     n_markers, min_cpm,
                                     ref_label,
                                     min_profiles = 2L) {

  counts_mat <- SummarizedExperiment::assay(sc_ref, "counts")
  cd <- as.data.frame(SummarizedExperiment::colData(sc_ref))

  message("\n[Signature Scoring] Reference: ", ref_label)
  message("  Profiles: ", ncol(sc_ref), " | Genes: ", nrow(sc_ref))
  message("  Cell types: ", length(unique(cd[[cell_type_col]])))
  message("  Donors: ", length(unique(cd[[sample_col]])))

  # Select reference profiles
  if (!is.null(ref_group) && !is.null(group_col) &&
      group_col %in% colnames(cd)) {
    ref_idx    <- cd[[group_col]] == ref_group
    ref_counts <- counts_mat[, ref_idx, drop = FALSE]
    ref_meta   <- cd[ref_idx, ]
    message("  Using group: '", ref_group, "' (",
            sum(ref_idx), " profiles)")
  } else {
    ref_counts <- counts_mat
    ref_meta   <- cd
    message("  Using all profiles")
  }

  # Build markers per cell type
  cell_types   <- unique(ref_meta[[cell_type_col]])
  markers_list <- list()
  markers_lfc  <- list()

  for (ct in cell_types) {
    is_ct  <- ref_meta[[cell_type_col]] == ct
    not_ct <- !is_ct

    if (sum(is_ct) < min_profiles || sum(not_ct) < min_profiles) {
      message("  [skip] ", ct, " -- ", sum(is_ct), " profiles")
      next
    }

    ct_cpm  <- sweep(ref_counts[, is_ct, drop = FALSE], 2,
                     colSums(ref_counts[, is_ct, drop = FALSE]), "/") * 1e6
    oth_cpm <- sweep(ref_counts[, not_ct, drop = FALSE], 2,
                     colSums(ref_counts[, not_ct, drop = FALSE]), "/") * 1e6

    ct_mean  <- rowMeans(ct_cpm)
    oth_mean <- rowMeans(oth_cpm)
    lfc      <- log2((ct_mean + 1) / (oth_mean + 1))

    candidates  <- lfc[ct_mean > min_cpm]
    top_markers <- names(sort(candidates, decreasing = TRUE))[
      seq_len(min(n_markers, length(candidates)))]
    top_markers <- top_markers[!is.na(top_markers)]

    if (length(top_markers) < 5) {
      message("  [skip] ", ct, " -- only ", length(top_markers), " markers")
      next
    }

    markers_list[[ct]] <- top_markers
    markers_lfc[[ct]]  <- lfc[top_markers]
    message("  ", ct, ": ", length(top_markers), " markers (top 3: ",
            paste(utils::head(top_markers, 3), collapse = ", "), ")")
  }

  if (length(markers_list) == 0)
    stop("No cell types with enough markers.", call. = FALSE)

  # Score bulk samples
  vsd     <- omni_vst(obj)
  if (is.null(vsd))
    stop("VST not cached. Run run_omni_de() first.", call. = FALSE)
  vst_mat <- SummarizedExperiment::assay(vsd)

  bulk_genes <- rownames(vst_mat)
  sc_genes   <- rownames(counts_mat)
  overlap    <- length(intersect(bulk_genes, sc_genes))
  message("  Gene overlap: ", overlap, " / ", length(sc_genes),
          " (", round(overlap / length(sc_genes) * 100, 1), "%)")

  ct_scores  <- matrix(NA_real_, nrow = length(markers_list),
                       ncol = ncol(vst_mat),
                       dimnames = list(names(markers_list),
                                       colnames(vst_mat)))
  genes_used <- list()

  for (ct in names(markers_list)) {
    genes <- intersect(markers_list[[ct]], bulk_genes)
    if (length(genes) < 5) next
    sub_mat <- vst_mat[genes, , drop = FALSE]
    z_mat   <- t(scale(t(sub_mat)))
    ct_scores[ct, ] <- colMeans(z_mat, na.rm = TRUE)
    genes_used[[ct]] <- genes
    message("  Scored ", ct, ": ", length(genes), " genes")
  }

  valid_rows <- !apply(ct_scores, 1, function(x) all(is.na(x)))
  ct_scores  <- ct_scores[valid_rows, , drop = FALSE]

  # IRG-S correlation
  irgs_vec      <- NULL
  cor_with_irgs <- NULL

  isgs <- tryCatch(OmniSeq::omniseq_isgs_28, error = function(e) NULL)
  if (!is.null(isgs)) {
    isgs_in_vst <- intersect(isgs, bulk_genes)
    if (length(isgs_in_vst) >= 5) {
      isg_z    <- scale(t(vst_mat[isgs_in_vst, ]))
      irgs_vec <- rowMeans(isg_z, na.rm = TRUE)
      names(irgs_vec) <- colnames(vst_mat)

      meta     <- omni_metadata(obj)
      cond_col <- obj$params$condition
      ref_lev  <- levels(meta[[cond_col]])[1]

      cor_rows <- lapply(rownames(ct_scores), function(ct) {
        ct_test <- stats::cor.test(ct_scores[ct, ], irgs_vec)
        ctrl_m  <- mean(ct_scores[ct, meta[[cond_col]] == ref_lev],
                        na.rm = TRUE)
        pat_m   <- mean(ct_scores[ct, meta[[cond_col]] != ref_lev],
                        na.rm = TRUE)
        data.frame(CellType = ct,
                   R = as.numeric(ct_test$estimate),
                   P = ct_test$p.value,
                   Ctrl_mean = ctrl_m, Pat_mean = pat_m,
                   Diff = pat_m - ctrl_m,
                   stringsAsFactors = FALSE)
      })
      cor_with_irgs <- do.call(rbind, cor_rows)
      cor_with_irgs <- cor_with_irgs[order(-abs(cor_with_irgs$R)), ]

      message("\n  Cell type correlations with IRG-S:")
      for (i in seq_len(nrow(cor_with_irgs))) {
        r <- cor_with_irgs[i, ]
        sig <- if (r$P < 0.05) " ***" else if (r$P < 0.1) " *" else ""
        message(sprintf("    %-22s R=%+.3f  p=%.4f  diff=%+.3f%s",
                        r$CellType, r$R, r$P, r$Diff, sig))
      }
    }
  }

  # Store
  obj$sc_deconv <- list(
    method         = "signature",
    scores         = ct_scores,
    markers        = markers_list,
    markers_lfc    = markers_lfc,
    genes_used     = genes_used,
    reference      = ref_label,
    ref_group      = ref_group %||% "all",
    n_ref_donors   = length(unique(ref_meta[[sample_col]])),
    n_ref_profiles = nrow(ref_meta),
    cell_types     = rownames(ct_scores),
    gene_overlap   = overlap,
    cor_with_irgs  = cor_with_irgs,
    irgs_scores    = irgs_vec,
    params         = list(cell_type_col = cell_type_col,
                          sample_col = sample_col,
                          group_col = group_col,
                          ref_group = ref_group,
                          n_markers = n_markers,
                          min_cpm = min_cpm)
  )

  message("\n[Signature Scoring] Done. ",
          nrow(ct_scores), " cell types scored.")
  return(obj)
}


# =============================================================================
# SECTION 3: RAW SCE -> MuSiC DECONVOLUTION
# =============================================================================

#' @keywords internal
.omni_sc_music <- function(obj, sc_ref, cell_type_col, sample_col) {

  .require_pkg("MuSiC")

  message("\n[MuSiC] Running deconvolution...")
  message("  Reference: ", ncol(sc_ref), " cells")

  bulk_mat  <- omni_counts(obj)
  music_res <- MuSiC::music_prop(
    bulk.mtx = bulk_mat,
    sc.sce   = sc_ref,
    clusters = cell_type_col,
    samples  = sample_col,
    verbose  = FALSE
  )

  deconv_mat <- t(music_res$Est.prop.weighted)
  message("[MuSiC] Estimated ", nrow(deconv_mat), " cell types.")

  obj$sc_deconv <- list(
    method     = "music",
    scores     = deconv_mat,
    reference  = "MuSiC (raw scRNA-seq)",
    cell_types = rownames(deconv_mat),
    params     = list(cell_type_col = cell_type_col,
                      sample_col = sample_col)
  )

  obj$deconv <- deconv_mat
  obj$params$deconv_method <- "MuSiC (scRNA-seq reference)"

  return(obj)
}


# =============================================================================
# SECTION 4: PSEUDO-BULK -> BisqueRNA (NNLS)
# =============================================================================

#' @keywords internal
.omni_sc_bisque <- function(obj, sc_ref,
                            cell_type_col, sample_col,
                            group_col = NULL, ref_group = NULL) {

  if (!requireNamespace("BisqueRNA", quietly = TRUE)) {
    message("[BisqueRNA] Not installed. Install: install.packages('BisqueRNA')")
    return(obj)
  }
  .require_pkg("Biobase")

  message("\n[BisqueRNA] Running reference-based decomposition...")

  cd <- as.data.frame(SummarizedExperiment::colData(sc_ref))
  counts_mat <- SummarizedExperiment::assay(sc_ref, "counts")

  if (!is.null(ref_group) && !is.null(group_col) &&
      group_col %in% colnames(cd)) {
    ref_idx    <- cd[[group_col]] == ref_group
    counts_mat <- counts_mat[, ref_idx, drop = FALSE]
    cd         <- cd[ref_idx, ]
  }

  sc_pheno <- data.frame(
    cellType    = cd[[cell_type_col]],
    SubjectName = cd[[sample_col]],
    row.names   = colnames(counts_mat)
  )
  sc_eset <- Biobase::ExpressionSet(
    assayData = as.matrix(counts_mat),
    phenoData = Biobase::AnnotatedDataFrame(sc_pheno)
  )

  bulk_mat     <- omni_counts(obj)
  common_genes <- intersect(rownames(bulk_mat), rownames(counts_mat))
  message("  Gene overlap: ", length(common_genes))

  bulk_eset <- Biobase::ExpressionSet(
    assayData = as.matrix(bulk_mat[common_genes, ])
  )

  bisque_res <- tryCatch(
    BisqueRNA::ReferenceBasedDecomposition(
      bulk.eset   = bulk_eset,
      sc.eset     = sc_eset,
      markers     = NULL,
      use.overlap = TRUE
    ),
    error = function(e) {
      message("  BisqueRNA failed: ", e$message); NULL })

  if (!is.null(bisque_res)) {
    bisque_props <- bisque_res$bulk.props
    message("[BisqueRNA] Estimated ", nrow(bisque_props), " cell types.")
    if (is.null(obj$sc_deconv)) obj$sc_deconv <- list()
    obj$sc_deconv$bisque_proportions <- bisque_props
    obj$sc_deconv$bisque_method      <- "BisqueRNA (NNLS)"
  }

  return(obj)
}


# =============================================================================
# SECTION 5: DISEASE ACTIVATION STATE SCORING (Layer 2)
# =============================================================================

#' Score Disease Activation States from Untreated Reference
#'
#' Second layer: cell identities come from healthy controls (Layer 1),
#' activation states come from comparing untreated disease vs healthy
#' WITHIN each cell type. This captures IFN-driven activation that is
#' cell-type-agnostic.
#'
#' @param obj A validated OmniSeq object (ideally with sc_deconv from Layer 1).
#' @param sc_ref The same SingleCellExperiment object.
#' @param disease_group Character. Group label for disease profiles.
#' @param healthy_group Character. Group label for healthy profiles.
#' @param cell_type_col Character. Cell type column.
#' @param group_col Character. Group column.
#' @param n_markers Integer. Markers per state. Default: `30`.
#' @param min_cpm Numeric. Min CPM for marker genes. Default: `5`.
#'
#' @return The OmniSeq object with disease state scores in `$sc_deconv`.
#' @export
run_omni_disease_state <- function(obj,
                                   sc_ref,
                                   disease_group = "SAVI_untreated",
                                   healthy_group = "Control",
                                   cell_type_col = "cell_type",
                                   group_col     = "group",
                                   n_markers     = 30L,
                                   min_cpm       = 5) {

  validate_omniseq(obj)
  .require_pkg("SingleCellExperiment")
  .require_pkg("SummarizedExperiment")

  cd <- as.data.frame(SummarizedExperiment::colData(sc_ref))
  counts_mat <- SummarizedExperiment::assay(sc_ref, "counts")

  if (!disease_group %in% cd[[group_col]])
    stop("disease_group '", disease_group, "' not found. Available: ",
         paste(unique(cd[[group_col]]), collapse = ", "), call. = FALSE)
  if (!healthy_group %in% cd[[group_col]])
    stop("healthy_group '", healthy_group, "' not found.", call. = FALSE)

  message("\n[Disease State Scoring] Layer 2: Activation states")
  message("  Comparing: ", disease_group, " vs ", healthy_group,
          " (within each cell type)")

  vsd     <- omni_vst(obj)
  if (is.null(vsd))
    stop("VST not cached. Run run_omni_de() first.", call. = FALSE)
  vst_mat <- SummarizedExperiment::assay(vsd)

  all_cell_types <- unique(cd[[cell_type_col]])
  state_markers  <- list()
  state_scores   <- matrix(NA_real_,
                           nrow     = 0,
                           ncol     = ncol(vst_mat),
                           dimnames = list(NULL, colnames(vst_mat)))

  for (ct in all_cell_types) {
    ct_dis  <- cd[[group_col]] == disease_group &
      cd[[cell_type_col]] == ct
    ct_ctrl <- cd[[group_col]] == healthy_group &
      cd[[cell_type_col]] == ct

    if (sum(ct_dis) < 1 || sum(ct_ctrl) < 1) {
      message("  [skip] ", ct, " -- no disease/healthy pair")
      next
    }

    dis_cpm  <- sweep(counts_mat[, ct_dis, drop = FALSE], 2,
                      colSums(counts_mat[, ct_dis, drop = FALSE]),
                      "/") * 1e6
    ctrl_cpm <- sweep(counts_mat[, ct_ctrl, drop = FALSE], 2,
                      colSums(counts_mat[, ct_ctrl, drop = FALSE]),
                      "/") * 1e6

    dis_mean  <- rowMeans(dis_cpm)
    ctrl_mean <- rowMeans(ctrl_cpm)
    lfc       <- log2((dis_mean + 1) / (ctrl_mean + 1))

    # Genes upregulated in disease state
    up_genes <- names(sort(lfc[dis_mean > min_cpm],
                           decreasing = TRUE))[
                             seq_len(min(n_markers, sum(dis_mean > min_cpm)))]
    up_genes <- up_genes[!is.na(up_genes)]

    if (length(up_genes) < 3) {
      message("  [skip] ", ct, " -- only ", length(up_genes), " markers")
      next
    }

    # Score bulk samples
    genes_in_bulk <- intersect(up_genes, rownames(vst_mat))
    if (length(genes_in_bulk) < 3) {
      message("  [skip] ", ct, " -- only ", length(genes_in_bulk),
              " in bulk")
      next
    }

    z_mat <- t(scale(t(vst_mat[genes_in_bulk, , drop = FALSE])))
    score <- colMeans(z_mat, na.rm = TRUE)

    state_name <- paste0(ct, "_activated")
    state_markers[[state_name]] <- up_genes

    # Add row to scores matrix
    new_row <- matrix(score, nrow = 1,
                      dimnames = list(state_name, names(score)))
    state_scores <- rbind(state_scores, new_row)

    message("  ", state_name, ": ", length(genes_in_bulk),
            " markers (top 5: ",
            paste(utils::head(up_genes, 5), collapse = ", "), ")")
  }

  if (nrow(state_scores) == 0) {
    message("[Disease State] No states could be scored.")
    return(obj)
  }

  # Compute correlation with IRG-S
  state_cor_irgs <- NULL
  if (!is.null(obj$sc_deconv) && !is.null(obj$sc_deconv$irgs_scores)) {
    irgs <- obj$sc_deconv$irgs_scores
    meta <- omni_metadata(obj)
    cond_col <- obj$params$condition
    ref_lev  <- levels(meta[[cond_col]])[1]

    cor_rows <- lapply(rownames(state_scores), function(st) {
      ct_test <- stats::cor.test(state_scores[st, ], irgs)
      ctrl_m  <- mean(state_scores[st, meta[[cond_col]] == ref_lev],
                      na.rm = TRUE)
      pat_m   <- mean(state_scores[st, meta[[cond_col]] != ref_lev],
                      na.rm = TRUE)
      data.frame(State = st, R = as.numeric(ct_test$estimate),
                 P = ct_test$p.value,
                 Ctrl_mean = ctrl_m, Pat_mean = pat_m,
                 Diff = pat_m - ctrl_m,
                 stringsAsFactors = FALSE)
    })
    state_cor_irgs <- do.call(rbind, cor_rows)
    state_cor_irgs <- state_cor_irgs[order(-abs(state_cor_irgs$R)), ]
  }

  # Store in sc_deconv
  if (is.null(obj$sc_deconv)) obj$sc_deconv <- list()
  obj$sc_deconv$disease_states <- list(
    scores        = state_scores,
    markers       = state_markers,
    disease_group = disease_group,
    healthy_group = healthy_group,
    cor_with_irgs = state_cor_irgs
  )

  message("\n[Disease State] Done. ", nrow(state_scores),
          " activation states scored.")
  return(obj)
}


# =============================================================================
# SECTION 6: CELL TYPE SCORE HEATMAP
# =============================================================================

#' Plot Cell Type Signature Score Heatmap
#'
#' @param obj OmniSeq object with sc_deconv.
#' @param layer Character. `"composition"`, `"activation"`, or `"both"`.
#' @param cluster_cols Logical. Default: `FALSE`.
#' @param annot_cols Character vector. Default: auto.
#'
#' @return A pheatmap object (invisibly).
#' @export
plot_omni_sc_heatmap <- function(obj,
                                 layer        = "composition",
                                 cluster_cols = FALSE,
                                 annot_cols   = NULL) {

  validate_omniseq(obj)
  .require_pkg("pheatmap")

  sc <- obj$sc_deconv
  if (is.null(sc))
    stop("No sc_deconv. Run run_omni_sc_deconv() first.", call. = FALSE)

  meta     <- omni_metadata(obj)
  cond_col <- obj$params$condition

  if (is.null(annot_cols))
    annot_cols <- intersect(c(cond_col, "Genotype", "Phenotype"),
                            colnames(meta))

  annot_df <- meta[, annot_cols, drop = FALSE]
  annot_df[] <- lapply(annot_df,
                       function(x) if (is.factor(x)) droplevels(x) else x)

  pal <- grDevices::colorRampPalette(
    c("#2166AC", "white", "#B2182B"))(100)

  # Select scores matrix
  if (layer == "activation" && !is.null(sc$disease_states)) {
    scores <- sc$disease_states$scores
    title_layer <- "Disease Activation States"
  } else {
    scores <- sc$scores
    title_layer <- "Cell Type Composition"
  }

  sample_order <- rownames(meta)[order(meta[[cond_col]])]
  sample_order <- intersect(sample_order, colnames(scores))

  p <- pheatmap::pheatmap(
    scores[, sample_order],
    annotation_col = annot_df[sample_order, , drop = FALSE],
    scale          = "row",
    cluster_cols   = cluster_cols,
    cluster_rows   = TRUE,
    color          = pal,
    border_color   = NA,
    fontsize_row   = 11,
    main           = paste0(sc$reference %||% "scRNA-seq", " -- ",
                            title_layer, "\n(z-scored, ",
                            ncol(scores), " samples)")
  )

  return(invisible(p))
}


# =============================================================================
# SECTION 7: PER-GENOTYPE BOXPLOTS WITH SIGNIFICANCE
# =============================================================================

#' Plot Cell Type Scores by Genotype with Significance Marks
#'
#' @param obj OmniSeq object with sc_deconv.
#' @param layer Character. `"composition"` or `"activation"`.
#' @param genotype_col Character. Default: `"Genotype"`.
#' @param cell_types Character vector. Default: all.
#' @param ref_genotype Character. Control genotype for Wilcoxon tests.
#'   Default: auto-detected.
#'
#' @return A ggplot object (invisibly).
#' @export
plot_omni_sc_genotype <- function(obj,
                                  layer        = "composition",
                                  genotype_col = "Genotype",
                                  cell_types   = NULL,
                                  ref_genotype = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .require_pkg("tidyr")

  sc <- obj$sc_deconv
  if (is.null(sc))
    stop("No sc_deconv. Run run_omni_sc_deconv() first.", call. = FALSE)

  meta <- omni_metadata(obj)

  if (!genotype_col %in% colnames(meta))
    stop("Column '", genotype_col, "' not in metadata.", call. = FALSE)

  # Select scores
  if (layer == "activation" && !is.null(sc$disease_states)) {
    scores    <- sc$disease_states$scores
    y_label   <- "Activation score (z)"
    title_pfx <- "Disease Activation States"
  } else {
    scores    <- sc$scores
    y_label   <- "Composition score (z)"
    title_pfx <- "Cell Type Composition"
  }

  if (is.null(cell_types)) cell_types <- rownames(scores)
  cell_types <- intersect(cell_types, rownames(scores))

  # Detect control genotype
  cond_col <- obj$params$condition
  if (is.null(ref_genotype)) {
    ref_lev <- levels(meta[[cond_col]])[1]
    ref_genotype <- unique(
      meta[[genotype_col]][meta[[cond_col]] == ref_lev])[1]
  }

  # Build data frame
  plot_data <- data.frame(
    Sample   = colnames(scores),
    Genotype = meta[colnames(scores), genotype_col],
    stringsAsFactors = FALSE
  )
  for (ct in cell_types) plot_data[[ct]] <- scores[ct, ]

  ct_long <- tidyr::pivot_longer(plot_data, cols = all_of(cell_types),
                                 names_to = "CellType",
                                 values_to = "Score")

  # Clean cell type labels for display
  ct_labels <- gsub("_activated$", "", cell_types)
  ct_labels <- gsub("_", " ", ct_labels)
  ct_long$CellType <- factor(ct_long$CellType,
                             levels = cell_types, labels = ct_labels)

  # Order genotype: ref first
  gt_levels <- unique(plot_data$Genotype)
  gt_levels <- c(ref_genotype,
                 sort(setdiff(gt_levels, ref_genotype)))
  ct_long$Genotype <- factor(ct_long$Genotype, levels = gt_levels)

  # Color palette
  geno_cols <- c("#2980B9", "#E74C3C", "#E67E22", "#8E44AD",
                 "#27AE60", "#16A085")
  fill_vals <- stats::setNames(geno_cols[seq_along(gt_levels)],
                               gt_levels)

  # ── Compute significance: each genotype vs ref ─────────────────────
  signif_df <- do.call(rbind, lapply(seq_along(cell_types), function(i) {
    ct <- cell_types[i]
    ct_label <- ct_labels[i]
    ctrl_vals <- scores[ct, meta[[genotype_col]] == ref_genotype]

    do.call(rbind, lapply(setdiff(gt_levels, ref_genotype), function(gt) {
      gt_vals <- scores[ct, meta[[genotype_col]] == gt]
      if (length(gt_vals) < 2) return(NULL)

      wt <- tryCatch(stats::wilcox.test(gt_vals, ctrl_vals)$p.value,
                     error = function(e) NA)
      if (is.na(wt)) return(NULL)

      label <- if (wt < 0.001) "***"
      else if (wt < 0.01) "**"
      else if (wt < 0.05) "*"
      else if (wt < 0.1) "."
      else ""

      if (nzchar(label)) {
        y_pos <- max(c(gt_vals, ctrl_vals), na.rm = TRUE) + 0.15
        data.frame(CellType = ct_label, Genotype = gt,
                   Label = label, P = wt, Y = y_pos,
                   stringsAsFactors = FALSE)
      } else NULL
    }))
  }))

  # ── Plot ───────────────────────────────────────────────────────────
  ref_label <- sc$reference %||% "scRNA-seq"

  p <- ggplot2::ggplot(
    ct_long, ggplot2::aes(x = Genotype, y = Score, fill = Genotype)
  ) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7) +
    ggplot2::geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    ggplot2::facet_wrap(~ CellType, scales = "free_y",
                        ncol = min(4, length(cell_types))) +
    ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      strip.text  = ggplot2::element_text(face = "bold", size = 10)
    ) +
    ggplot2::labs(
      title    = paste0(title_pfx, " by Genotype"),
      subtitle = paste0(ref_label,
                        " | Wilcoxon vs ", ref_genotype,
                        ": *** p<0.001  ** p<0.01  * p<0.05  . p<0.1"),
      y = y_label, x = NULL
    )

  # Add significance marks
  if (!is.null(signif_df) && nrow(signif_df) > 0) {
    signif_df$CellType <- factor(signif_df$CellType, levels = ct_labels)
    signif_df$Genotype <- factor(signif_df$Genotype, levels = gt_levels)
    p <- p + ggplot2::geom_text(
      data = signif_df,
      ggplot2::aes(x = Genotype, y = Y, label = Label),
      inherit.aes = FALSE, size = 5, fontface = "bold",
      color = "black", vjust = 0)
  }

  print(p)
  return(invisible(p))
}


# =============================================================================
# SECTION 8: CELL TYPE CORRELATION WITH FEATURES
# =============================================================================

#' Correlate Cell Type Scores with a Feature
#'
#' @param obj OmniSeq object with sc_deconv.
#' @param feature Character. `"IRGS"`, gene name, or metadata column.
#' @param layer Character. `"composition"` or `"activation"`.
#' @param method Character. Default: `"pearson"`.
#'
#' @return A ggplot object (invisibly).
#' @export
plot_omni_sc_correlation <- function(obj,
                                     feature = "IRGS",
                                     layer   = "composition",
                                     method  = "pearson") {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  sc <- obj$sc_deconv
  if (is.null(sc))
    stop("No sc_deconv.", call. = FALSE)

  # Select scores
  if (layer == "activation" && !is.null(sc$disease_states))
    ct_scores <- sc$disease_states$scores
  else
    ct_scores <- sc$scores

  meta <- omni_metadata(obj)

  # Resolve feature
  if (toupper(feature) == "IRGS" && !is.null(sc$irgs_scores)) {
    feat_vec   <- sc$irgs_scores
    feat_label <- "IRG-S Score"
  } else if (feature %in% colnames(meta) && is.numeric(meta[[feature]])) {
    feat_vec   <- stats::setNames(meta[[feature]], rownames(meta))
    feat_label <- feature
  } else {
    vsd <- omni_vst(obj)
    if (is.null(vsd)) stop("VST not cached.", call. = FALSE)
    vst_mat <- SummarizedExperiment::assay(vsd)
    if (!feature %in% rownames(vst_mat))
      stop("Feature '", feature, "' not found.", call. = FALSE)
    feat_vec   <- vst_mat[feature, ]
    feat_label <- paste("VST:", feature)
  }

  common <- intersect(colnames(ct_scores), names(feat_vec))

  cor_df <- do.call(rbind, lapply(rownames(ct_scores), function(ct) {
    ct_test <- stats::cor.test(ct_scores[ct, common],
                               feat_vec[common], method = method)
    data.frame(CellType = ct,
               R = as.numeric(ct_test$estimate),
               P = ct_test$p.value,
               stringsAsFactors = FALSE)
  }))

  cor_df$LogP     <- -log10(cor_df$P + 1e-10)
  cor_df          <- cor_df[order(-abs(cor_df$R)), ]
  cor_df$Label    <- gsub("_activated$", "", cor_df$CellType)
  cor_df$Label    <- gsub("_", " ", cor_df$Label)
  cor_df$Label    <- factor(cor_df$Label, levels = cor_df$Label)

  layer_label <- if (layer == "activation") "Activation" else "Composition"

  p <- ggplot2::ggplot(
    cor_df, ggplot2::aes(x = R, y = Label, size = LogP, color = R)
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        linewidth = 0.4) +
    ggplot2::scale_color_gradient2(
      low = "#2166AC", mid = "grey80", high = "#E74C3C",
      midpoint = 0, name = "R") +
    ggplot2::scale_size_continuous(name = "-log10(p)", range = c(3, 10)) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title    = paste(layer_label, "Scores vs", feat_label),
      subtitle = paste0(sc$reference %||% "scRNA-seq", " | ",
                        method, " correlation"),
      x = paste(method, "R"), y = NULL)

  print(p)
  return(invisible(p))
}


# =============================================================================
# SECTION 9: GENOTYPE IMMUNE PROFILE LINES
# =============================================================================

#' Plot Mean Immune Profile per Genotype
#'
#' @param obj OmniSeq object with sc_deconv.
#' @param layer Character. `"composition"` or `"activation"`.
#' @param genotype_col Character. Default: `"Genotype"`.
#'
#' @return A ggplot object (invisibly).
#' @export
plot_omni_sc_profile <- function(obj,
                                 layer        = "composition",
                                 genotype_col = "Genotype") {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .require_pkg("tidyr")

  sc <- obj$sc_deconv
  if (is.null(sc)) stop("No sc_deconv.", call. = FALSE)

  meta <- omni_metadata(obj)
  if (!genotype_col %in% colnames(meta))
    stop("Column '", genotype_col, "' not in metadata.", call. = FALSE)

  if (layer == "activation" && !is.null(sc$disease_states))
    scores <- sc$disease_states$scores
  else
    scores <- sc$scores

  geno_data <- data.frame(
    Genotype = meta[colnames(scores), genotype_col],
    stringsAsFactors = FALSE
  )
  for (ct in rownames(scores)) geno_data[[ct]] <- scores[ct, ]

  profile_df <- stats::aggregate(. ~ Genotype, data = geno_data, FUN = mean)

  profile_long <- tidyr::pivot_longer(
    profile_df, cols = -Genotype,
    names_to = "CellType", values_to = "MeanScore")

  profile_long$Label <- gsub("_activated$", "", profile_long$CellType)
  profile_long$Label <- gsub("_", " ", profile_long$Label)

  gt_levels <- unique(profile_long$Genotype)
  geno_cols <- c("#2980B9", "#E74C3C", "#E67E22", "#8E44AD",
                 "#27AE60", "#16A085")
  col_vals  <- stats::setNames(geno_cols[seq_along(gt_levels)], gt_levels)

  layer_label <- if (layer == "activation") "Activation" else "Composition"

  p <- ggplot2::ggplot(
    profile_long,
    ggplot2::aes(x = Label, y = MeanScore,
                 color = Genotype, group = Genotype)
  ) +
    ggplot2::geom_line(linewidth = 1.2, alpha = 0.8) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey50") +
    ggplot2::scale_color_manual(values = col_vals) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 11),
      legend.position = "top") +
    ggplot2::labs(
      title = paste("Immune", layer_label, "Profile by Genotype"),
      subtitle = sc$reference %||% "scRNA-seq",
      y = "Mean score (z)", x = NULL, color = NULL)

  print(p)
  message("\n  Mean scores by genotype:")
  print(round(profile_df, 3))
  return(invisible(p))
}


# =============================================================================
# SECTION 10: TWO-LAYER COMBINED FIGURE
# =============================================================================

#' Plot Two-Layer Deconvolution: Composition vs Activation
#'
#' Side-by-side comparison of Layer 1 (cell composition from healthy ref)
#' and Layer 2 (IFN activation from disease ref) with significance marks.
#'
#' @param obj OmniSeq object with both layers in sc_deconv.
#' @param genotype_col Character. Default: `"Genotype"`.
#' @param cell_types Character vector. Cell types to show (matched across
#'   both layers). Default: top 4 most informative.
#'
#' @return A patchwork ggplot object (invisibly).
#' @export
plot_omni_sc_two_layer <- function(obj,
                                   genotype_col = "Genotype",
                                   cell_types   = NULL) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .require_pkg("patchwork")

  sc <- obj$sc_deconv
  if (is.null(sc) || is.null(sc$scores))
    stop("No Layer 1 scores.", call. = FALSE)
  if (is.null(sc$disease_states) || is.null(sc$disease_states$scores))
    stop("No Layer 2 scores. Run run_omni_disease_state() first.",
         call. = FALSE)

  comp_scores <- sc$scores
  act_scores  <- sc$disease_states$scores

  # Auto-select cell types present in both layers
  comp_cts <- rownames(comp_scores)
  act_cts  <- gsub("_activated$", "", rownames(act_scores))

  if (is.null(cell_types)) {
    shared <- intersect(comp_cts, act_cts)
    if (length(shared) > 4)
      shared <- shared[seq_len(4)]  # top 4
    cell_types <- shared
  }

  p1 <- plot_omni_sc_genotype(obj, layer = "composition",
                              genotype_col = genotype_col,
                              cell_types = cell_types)

  act_cell_types <- paste0(cell_types, "_activated")
  act_cell_types <- intersect(act_cell_types, rownames(act_scores))

  p2 <- plot_omni_sc_genotype(obj, layer = "activation",
                              genotype_col = genotype_col,
                              cell_types = act_cell_types)

  combined <- p1 / p2 +
    patchwork::plot_annotation(
      title    = "Two-Layer Deconvolution: Composition vs Activation",
      subtitle = paste0(sc$reference %||% "scRNA-seq",
                        " | Layer 1: healthy ref | Layer 2: disease ref"),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 10, color = "grey40")
      )
    )

  print(combined)
  return(invisible(combined))
}


# =============================================================================
# SECTION 11: BENCHMARK — COMPARE METHODS
# =============================================================================

#' Compare Multiple Deconvolution Methods
#'
#' @param obj OmniSeq object with deconv and/or sc_deconv.
#'
#' @return A ggplot object (invisibly).
#' @export
plot_omni_sc_benchmark <- function(obj) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  if (is.null(obj$deconv) && is.null(obj$sc_deconv))
    stop("No deconvolution results.", call. = FALSE)

  meta     <- omni_metadata(obj)
  cond_col <- obj$params$condition
  dfs      <- list()

  # xCell results
  if (!is.null(obj$deconv)) {
    xcell <- obj$deconv
    for (pat in c("[Mm]onocyte|Mono", "NK", "B.cell|B_cell", "T.cell|T_cell")) {
      rows <- grep(pat, rownames(xcell), value = TRUE, ignore.case = TRUE)
      if (length(rows) == 0) next
      label <- gsub("\\[.*", "", pat)
      label <- gsub("\\.", " ", label)
      label <- trimws(gsub("\\|.*", "", label))
      dfs[[paste0("xCell_", label)]] <- data.frame(
        Sample = colnames(xcell),
        Score  = colMeans(xcell[rows, , drop = FALSE], na.rm = TRUE),
        Method = "xCell", CellType = label,
        stringsAsFactors = FALSE)
    }
  }

  # Signature scoring results
  if (!is.null(obj$sc_deconv) && !is.null(obj$sc_deconv$scores)) {
    sc_scores <- obj$sc_deconv$scores
    for (ct in rownames(sc_scores)) {
      dfs[[paste0("Sig_", ct)]] <- data.frame(
        Sample = colnames(sc_scores),
        Score  = sc_scores[ct, ],
        Method = "Signature", CellType = gsub("_", " ", ct),
        stringsAsFactors = FALSE)
    }
  }

  if (length(dfs) == 0)
    stop("No results to compare.", call. = FALSE)

  bench_df <- do.call(rbind, dfs)
  bench_df$Condition <- meta[bench_df$Sample, cond_col]

  p <- ggplot2::ggplot(
    bench_df,
    ggplot2::aes(x = Method, y = Score, fill = Condition)
  ) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::facet_wrap(~ CellType, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("#2980B9", "#E74C3C")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title    = "Deconvolution Method Comparison",
      subtitle = "xCell vs scRNA-seq Signature Scoring",
      y = "Score", x = NULL)

  print(p)
  return(invisible(p))
}


# =============================================================================
# SECTION 12: INTERNAL HELPERS
# =============================================================================

#' @keywords internal
.omni_get_sc_deconv <- function(obj) {
  if (is.null(obj$sc_deconv)) return(NULL)
  obj$sc_deconv$scores
}

#' @keywords internal
.omni_get_sc_markers <- function(obj) {
  if (is.null(obj$sc_deconv)) return(NULL)
  obj$sc_deconv$markers
}

#' @keywords internal
.omni_get_disease_states <- function(obj) {
  if (is.null(obj$sc_deconv) ||
      is.null(obj$sc_deconv$disease_states)) return(NULL)
  obj$sc_deconv$disease_states$scores
}
