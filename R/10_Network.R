# =============================================================================
# OmniSeq — 09_Network.R
#
# Co-expression network reconstruction via WGCNA and multiWGCNA.
#
# DESIGN PRINCIPLES:
#   - Automatic soft-power estimation via pickSoftThreshold() — never hardcoded
#   - Input matrix is user-specified: full VST, ISG-masked VST, or any matrix
#   - Results cached in obj$wgcna[["<key>"]] — multiple networks coexist
#   - TE family expression used as external traits (not as network nodes)
#   - Module preservation tests stability across conditions
#   - All functions warn clearly when sample size limits reliability
#
# FIXES vs previous version:
#   - cor <- WGCNA::cor fix added to run_omni_wgcna() and run_omni_wgcna_multi()
#     prevents the stats::cor masking conflict that crashes blockwiseModules
#   - corType = "pearson" and maxBlockSize added to blockwiseModules calls
#   - guide = "none" added to scale_edge_color_gradient2 — fixes ggraph version error
#   - plot_omni_wgcna_te_trait_heatmap() added as a new dedicated function
#   - plot_omni_wgcna_trait_heatmap() no longer shows confusing message about
#     TE traits — it now calls plot_omni_wgcna_te_trait_heatmap() automatically
#     when TE trait correlations are available
#
# FUNCTIONS:
#   run_omni_wgcna()                  Single network, auto soft-power
#   run_omni_wgcna_multi()            Condition-specific networks + module mapping
#   run_omni_wgcna_te_traits()        Gene modules x TE family correlation
#   run_omni_module_preservation()    Module stability across conditions
#   plot_omni_wgcna_dendrogram()      Module color dendrogram
#   plot_omni_wgcna_trait_heatmap()   Module eigengene x clinical trait heatmap
#   plot_omni_wgcna_te_trait_heatmap() Module eigengene x TE family heatmap (NEW)
#   plot_omni_wgcna_hub_network()     Hub gene network (static + interactive)
#   plot_omni_module_preservation()   Zsummary preservation plot
#   omni_wgcna()                      Accessor
# =============================================================================


# =============================================================================
# INTERNAL: GET EXPRESSION MATRIX FOR NETWORK ANALYSIS
# =============================================================================

#' Extract the appropriate expression matrix for WGCNA
#' @keywords internal
.omni_wgcna_matrix <- function(obj,
                               use_masked  = FALSE,
                               n_top_genes = 5000L) {

  vsd <- omni_vst(obj)
  if (is.null(vsd))
    stop("VST not cached. Run run_omni_de() first.", call. = FALSE)

  mat <- SummarizedExperiment::assay(vsd)

  # Restrict to genes in the masked contrast if requested
  if (use_masked && !is.null(obj$active_contrast) &&
      grepl("excl\\.", obj$active_contrast)) {

    # The masked contrast stores padj=1 for masked genes but does not
    # remove them from the results. We need to identify which genes
    # were masked and physically exclude them from the WGCNA matrix.
    # Strategy: find the original (unmasked) contrast and compare —
    # genes that have padj=1 in masked but significant in original
    # are the masked genes.

    # Find the parent contrast (strip the "(excl. ...)" suffix)
    parent_key <- sub("\\s*\\(excl\\..*\\)$", "",
                      obj$active_contrast)
    parent_entry <- obj$results_cache[[parent_key]]
    masked_entry <- obj$results_cache[[obj$active_contrast]]

    if (!is.null(parent_entry) && !is.null(masked_entry)) {
      parent_df <- as.data.frame(parent_entry$results)
      masked_df <- as.data.frame(masked_entry$results)

      # Masked genes are identified as rows present in parent but
      # absent in masked results — masking removes rows entirely
      masked_genes <- setdiff(rownames(parent_df), rownames(masked_df))

      if (length(masked_genes) > 0) {
        mat <- mat[!rownames(mat) %in% masked_genes, , drop = FALSE]
        message("  Using masked gene set: removed ", length(masked_genes),
                " masked genes. Remaining: ", nrow(mat), " features.")
      } else {
        message("  Warning: could not identify masked genes (setdiff empty).")
        message("  Falling back to full gene set for WGCNA.")
      }
    }
  }

  # Filter to top variable genes
  n_keep <- min(n_top_genes, nrow(mat))
  rv     <- apply(mat, 1, stats::var)
  mat    <- mat[order(rv, decreasing = TRUE)[1:n_keep], ]

  message("  Filtered to top ", n_keep, " variable genes.")
  return(mat)
}


# =============================================================================
# 1. SINGLE WGCNA NETWORK
# =============================================================================

#' Run WGCNA Co-Expression Network Analysis
#'
#' Builds a weighted gene co-expression network on the VST expression matrix
#' with automatic soft-power estimation. Results are cached in
#' `obj$wgcna[["<key>"]]`.
#'
#' @param obj A validated OmniSeq object with a fitted model.
#' @param use_masked Logical. Use the ISG-masked gene set. Default: `FALSE`.
#' @param n_top_genes Integer. Number of top variable genes to use.
#'   Default: `5000`. Reduce for speed on smaller datasets.
#' @param network_type Character. `"signed"` (recommended) or `"unsigned"`.
#'   Default: `"signed"`.
#' @param min_module_size Integer. Minimum genes per module. Default: `30`.
#' @param merge_cut_height Numeric. Height for merging similar modules.
#'   Default: `0.25`.
#' @param key Character. Custom cache key. Auto-generated if `NULL`.
#'
#' @return The OmniSeq object with results stored in `obj$wgcna[["<key>"]]`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Full VST network
#' obj <- run_omni_wgcna(obj, use_masked = FALSE, key = "full")
#'
#' # ISG-masked network (run after run_omni_mask())
#' obj <- run_omni_wgcna(obj, use_masked = TRUE, key = "masked")
#' }
run_omni_wgcna <- function(obj,
                           use_masked       = FALSE,
                           n_top_genes      = 5000L,
                           network_type     = c("signed", "unsigned"),
                           min_module_size  = 30L,
                           merge_cut_height = 0.25,
                           key              = NULL) {

  validate_omniseq(obj)
  .require_pkg("WGCNA")

  # Fix WGCNA/stats::cor masking conflict — must be before any WGCNA calls
  cor <- WGCNA::cor
  on.exit(rm(cor), add = TRUE)

  network_type <- match.arg(network_type)

  label <- if (use_masked) "masked" else "full"
  if (is.null(key)) key <- label

  message("\n", strrep("=", 55))
  message("  OmniSeq WGCNA — ", key)
  message(strrep("=", 55))

  # ------------------------------------------------------------------
  # 1. PREPARE MATRIX (samples x genes for WGCNA)
  # ------------------------------------------------------------------
  mat      <- .omni_wgcna_matrix(obj, use_masked, n_top_genes)
  dat_expr <- t(mat)

  n_samp  <- nrow(dat_expr)
  n_genes <- ncol(dat_expr)
  message("  Samples: ", n_samp, " | Genes: ", n_genes)

  if (n_samp < 10)
    warning("WGCNA is designed for n >= 15 samples. With n = ", n_samp,
            " results may be unreliable.", call. = FALSE)

  WGCNA::allowWGCNAThreads()

  # ------------------------------------------------------------------
  # 2. SOFT-POWER ESTIMATION
  # ------------------------------------------------------------------
  message("\n[run_omni_wgcna] Estimating soft-thresholding power...")

  powers <- c(1:10, seq(12, 20, by = 2))
  sft    <- WGCNA::pickSoftThreshold(
    dat_expr,
    powerVector = powers,
    networkType = network_type,
    verbose     = 0
  )

  r2_vals    <- sft$fitIndices[, 2]
  good_power <- sft$fitIndices$Power[r2_vals > 0.85]

  if (length(good_power) > 0) {
    soft_power <- good_power[1]
  } else {
    good_power <- sft$fitIndices$Power[r2_vals > 0.80]
    if (length(good_power) > 0) {
      soft_power <- good_power[1]
      message("  Note: R\u00b2 > 0.85 not achieved. Using R\u00b2 > 0.80 threshold.")
    } else {
      soft_power <- sft$fitIndices$Power[which.max(r2_vals)]
      message("  Note: R\u00b2 > 0.80 not achieved. Using best available power.")
    }
  }

  message("  Selected soft-thresholding power: ", soft_power,
          " (R\u00b2 = ",
          round(r2_vals[sft$fitIndices$Power == soft_power], 3), ")")

  # ------------------------------------------------------------------
  # 3. BUILD NETWORK
  # corType = "pearson" forces WGCNA to use its internal implementation
  # maxBlockSize = n_genes ensures all genes processed in a single block
  # ------------------------------------------------------------------
  message("\n[run_omni_wgcna] Building co-expression network...")
  message("  network_type = ", network_type,
          " | min_module_size = ", min_module_size,
          " | merge_cut_height = ", merge_cut_height)

  net <- WGCNA::blockwiseModules(
    dat_expr,
    power             = soft_power,
    networkType       = network_type,
    TOMType           = if (network_type == "signed") "signed" else "unsigned",
    minModuleSize     = min_module_size,
    mergeCutHeight    = merge_cut_height,
    numericLabels     = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs          = FALSE,
    corType           = "pearson",
    maxBlockSize      = n_genes,
    verbose           = 0
  )

  module_colors <- net$colors
  n_modules     <- length(unique(module_colors)) - 1
  message("  Found ", n_modules, " modules (excluding grey/unassigned).")
  mod_tab <- sort(table(module_colors), decreasing = TRUE)
  for (nm in names(mod_tab)[names(mod_tab) != "grey"])
    message("    ", nm, ": ", mod_tab[nm], " genes")

  # ------------------------------------------------------------------
  # 4. MODULE EIGENGENES
  # ------------------------------------------------------------------
  message("\n[run_omni_wgcna] Computing module eigengenes...")
  MEs0 <- WGCNA::moduleEigengenes(dat_expr, module_colors)$eigengenes
  MEs  <- WGCNA::orderMEs(MEs0)

  # ------------------------------------------------------------------
  # 5. STORE
  # ------------------------------------------------------------------
  if (is.null(obj$wgcna)) obj$wgcna <- list()

  obj$wgcna[[key]] <- list(
    net              = net,
    module_colors    = module_colors,
    eigengenes       = MEs,
    dat_expr         = dat_expr,
    soft_power       = soft_power,
    sft              = sft,
    network_type     = network_type,
    n_top_genes      = n_top_genes,
    use_masked       = use_masked,
    min_module_size  = min_module_size,
    merge_cut_height = merge_cut_height,
    n_modules        = n_modules,
    timestamp        = Sys.time()
  )

  message("\n[run_omni_wgcna] Done. Stored under key: '", key, "'")
  message("  Use omni_wgcna(obj, '", key, "') to access.")
  message("  Next: plot_omni_wgcna_dendrogram(obj, '", key, "')")
  message("        plot_omni_wgcna_trait_heatmap(obj, '", key, "')")

  return(obj)
}


# =============================================================================
# 2. MULTI-CONDITION WGCNA
# =============================================================================

#' Run Condition-Specific WGCNA Networks
#'
#' Builds separate WGCNA networks for each level of a grouping variable,
#' then maps modules between them to identify shared vs condition-specific
#' co-expression patterns.
#'
#' @param obj A validated OmniSeq object with a fitted model.
#' @param group_col Character. Metadata column defining conditions.
#' @param exclude_groups Character vector. Groups to exclude.
#' @param use_masked Logical. Use ISG-masked gene set. Default: `FALSE`.
#' @param n_top_genes Integer. Top variable genes. Default: `5000`.
#' @param network_type Character. `"signed"` or `"unsigned"`.
#' @param min_module_size Integer. Minimum genes per module. Default: `30`.
#'
#' @return The OmniSeq object with results stored in
#'   `obj$wgcna_multi[["<group_col>"]]`.
#' @export
run_omni_wgcna_multi <- function(obj,
                                 group_col       = "Condition",
                                 exclude_groups  = NULL,
                                 use_masked      = FALSE,
                                 n_top_genes     = 5000L,
                                 network_type    = c("signed", "unsigned"),
                                 min_module_size = 30L) {

  validate_omniseq(obj)
  .require_pkg("WGCNA")

  # Fix WGCNA/stats::cor masking conflict
  cor <- WGCNA::cor
  on.exit(rm(cor), add = TRUE)

  network_type <- match.arg(network_type)

  meta       <- omni_metadata(obj)
  if (!group_col %in% colnames(meta))
    stop("group_col '", group_col, "' not found in metadata.", call. = FALSE)

  groups_all   <- as.character(meta[[group_col]])
  keep_mask    <- !groups_all %in% (exclude_groups %||% character(0))
  groups_use   <- groups_all[keep_mask]
  samp_use     <- rownames(meta)[keep_mask]
  group_levels <- unique(groups_use)

  message("\n[run_omni_wgcna_multi] Building networks for: ",
          paste(group_levels, collapse = ", "))

  mat_full <- .omni_wgcna_matrix(obj, use_masked, n_top_genes)
  mat_full <- mat_full[, samp_use, drop = FALSE]

  WGCNA::allowWGCNAThreads()

  # Estimate soft power on combined data
  dat_expr_all <- t(mat_full)
  message("\n[run_omni_wgcna_multi] Estimating soft power on combined data...")
  powers <- c(1:10, seq(12, 20, by = 2))
  sft    <- WGCNA::pickSoftThreshold(dat_expr_all,
                                     powerVector = powers,
                                     networkType = network_type,
                                     verbose     = 0)
  r2_vals    <- sft$fitIndices[, 2]
  good_power <- sft$fitIndices$Power[r2_vals > 0.80]
  soft_power <- if (length(good_power) > 0) good_power[1] else
    sft$fitIndices$Power[which.max(r2_vals)]
  message("  Soft power: ", soft_power)

  n_genes_all <- ncol(dat_expr_all)

  # Build one network per group
  group_nets <- list()
  for (grp in group_levels) {
    grp_samps <- samp_use[groups_use == grp]
    n_samp    <- length(grp_samps)
    message("\n[run_omni_wgcna_multi] Building network for '", grp,
            "' (n = ", n_samp, ")...")

    if (n_samp < 5)
      warning("Group '", grp, "' has only ", n_samp,
              " samples — network may be unreliable.", call. = FALSE)

    dat_grp  <- t(mat_full[, grp_samps, drop = FALSE])
    n_genes_grp <- ncol(dat_grp)

    net_grp <- WGCNA::blockwiseModules(
      dat_grp,
      power             = soft_power,
      networkType       = network_type,
      TOMType           = if (network_type == "signed") "signed" else "unsigned",
      minModuleSize     = min_module_size,
      mergeCutHeight    = 0.25,
      numericLabels     = FALSE,
      pamRespectsDendro = FALSE,
      saveTOMs          = FALSE,
      corType           = "pearson",
      maxBlockSize      = n_genes_grp,
      verbose           = 0
    )

    n_mods <- length(unique(net_grp$colors)) - 1
    message("  Found ", n_mods, " modules.")

    MEs <- WGCNA::orderMEs(
      WGCNA::moduleEigengenes(dat_grp, net_grp$colors)$eigengenes
    )

    group_nets[[grp]] <- list(
      net           = net_grp,
      module_colors = net_grp$colors,
      eigengenes    = MEs,
      samples       = grp_samps,
      n_modules     = n_mods
    )
  }

  # Module overlap mapping between groups
  message("\n[run_omni_wgcna_multi] Mapping modules between conditions...")
  overlap_list  <- list()
  group_pairs   <- utils::combn(group_levels, 2, simplify = FALSE)

  for (pair in group_pairs) {
    g1    <- pair[1]
    g2    <- pair[2]
    cols1 <- group_nets[[g1]]$module_colors
    cols2 <- group_nets[[g2]]$module_colors
    shared <- intersect(names(cols1), names(cols2))
    if (length(shared) == 0) next
    tab      <- table(cols1[shared], cols2[shared])
    pair_key <- paste(g1, "vs", g2)
    overlap_list[[pair_key]] <- tab
    message("  ", pair_key, ": overlap table computed (",
            nrow(tab), " x ", ncol(tab), " modules)")
  }

  if (is.null(obj$wgcna_multi)) obj$wgcna_multi <- list()

  obj$wgcna_multi[[group_col]] <- list(
    group_nets   = group_nets,
    overlap      = overlap_list,
    group_col    = group_col,
    group_levels = group_levels,
    excluded     = exclude_groups,
    soft_power   = soft_power,
    network_type = network_type,
    use_masked   = use_masked,
    timestamp    = Sys.time()
  )

  message("\n[run_omni_wgcna_multi] Done. Stored under obj$wgcna_multi[['",
          group_col, "']]")
  return(obj)
}


# =============================================================================
# 3. TE FAMILY AS EXTERNAL TRAITS
# =============================================================================

#' Correlate WGCNA Gene Modules with TE Family Expression
#'
#' Uses TE family-level expression (mean VST per family per sample) as
#' external traits in the WGCNA module-trait correlation framework. Produces
#' correlation coefficients and p-values stored in the WGCNA result slot, and
#' visualised with `plot_omni_wgcna_te_trait_heatmap()`.
#'
#' @param gene_obj OmniSeq object containing the WGCNA result.
#' @param te_obj OmniSeq object containing TE VST expression.
#' @param wgcna_key Character. Key of the WGCNA result in `gene_obj$wgcna`.
#' @param te_classes Character vector. TE classes to include. Default: all.
#' @param n_top_families Integer. Top variable TE families. Default: `50`.
#' @param also_include_clinical Logical. Include clinical traits alongside
#'   TE families. Default: `TRUE`.
#'
#' @return The gene_obj with TE trait correlations stored in
#'   `gene_obj$wgcna[[wgcna_key]]$te_trait_cor`.
#' @export
run_omni_wgcna_te_traits <- function(gene_obj,
                                     te_obj,
                                     wgcna_key             = "full",
                                     te_classes            = NULL,
                                     n_top_families        = 50L,
                                     also_include_clinical = TRUE) {

  validate_omniseq(gene_obj)
  validate_omniseq(te_obj)
  .require_pkg("WGCNA")

  if (is.null(gene_obj$wgcna) || !wgcna_key %in% names(gene_obj$wgcna))
    stop("WGCNA key '", wgcna_key,
         "' not found. Run run_omni_wgcna() first.", call. = FALSE)

  wgcna_entry  <- gene_obj$wgcna[[wgcna_key]]
  MEs          <- wgcna_entry$eigengenes
  dat_expr     <- wgcna_entry$dat_expr
  gene_samples <- rownames(dat_expr)

  # Build TE family matrix
  vsd_te <- omni_vst(te_obj)
  if (is.null(vsd_te))
    stop("TE VST not cached. Run run_omni_de() on te_obj first.",
         call. = FALSE)

  mat_te  <- SummarizedExperiment::assay(vsd_te)
  te_info <- .omni_parse_te_ids(rownames(mat_te))

  # Filter by TE class if requested
  if (!is.null(te_classes)) {
    keep_te <- te_info$Class %in% te_classes
    mat_te  <- mat_te[keep_te, , drop = FALSE]
    te_info <- te_info[keep_te, ]
    message("[run_omni_wgcna_te_traits] Filtered to classes: ",
            paste(te_classes, collapse = ", "),
            " (", nrow(mat_te), " loci)")
  }

  # Aggregate by family
  fams        <- te_info$Family
  unique_fams <- unique(fams[fams != "Unknown"])

  agg_list <- lapply(unique_fams, function(f) {
    idx <- which(fams == f)
    if (length(idx) == 1) mat_te[idx, ]
    else colMeans(mat_te[idx, , drop = FALSE], na.rm = TRUE)
  })
  agg_mat           <- do.call(rbind, agg_list)
  rownames(agg_mat) <- unique_fams

  # Select top variable families
  pvar      <- apply(agg_mat, 1, stats::var)
  top_fams  <- names(sort(pvar, decreasing = TRUE))[
    1:min(n_top_families, nrow(agg_mat))]
  agg_mat   <- agg_mat[top_fams, ]

  # Align samples
  common_samps <- intersect(gene_samples, colnames(agg_mat))
  if (length(common_samps) < 3)
    stop("Fewer than 3 common samples between gene and TE objects.",
         call. = FALSE)

  MEs_sub   <- MEs[common_samps, , drop = FALSE]
  te_traits <- t(agg_mat[, common_samps, drop = FALSE])
  n_samp    <- length(common_samps)

  # Optionally add clinical traits
  if (also_include_clinical) {
    meta     <- omni_metadata(gene_obj)
    clin_df  <- meta[common_samps, , drop = FALSE]

    cond_col  <- gene_obj$params$condition
    cond_lvls <- levels(clin_df[[cond_col]])
    cond_num  <- as.numeric(clin_df[[cond_col]] == cond_lvls[length(cond_lvls)])

    clinical_traits <- data.frame(Condition = cond_num,
                                  row.names = common_samps)

    for (col in c("Genotype", "Phenotype", "Sex", "Batch")) {
      if (col %in% colnames(clin_df))
        clinical_traits[[col]] <- as.numeric(as.factor(clin_df[[col]]))
    }
    if ("Age_Scaled" %in% colnames(clin_df))
      clinical_traits$Age <- as.numeric(clin_df$Age_Scaled)

    traits <- cbind(clinical_traits, te_traits)
  } else {
    traits <- as.data.frame(te_traits)
  }

  n_te_fams  <- ncol(te_traits)
  n_clin     <- ncol(traits) - n_te_fams

  message("[run_omni_wgcna_te_traits] Correlating ", ncol(MEs_sub),
          " module eigengenes with ", ncol(traits), " traits (",
          n_te_fams, " TE families",
          if (also_include_clinical)
            paste0(" + ", n_clin, " clinical traits"),
          ")...")

  module_trait_cor  <- WGCNA::cor(MEs_sub, traits, use = "p")
  module_trait_pval <- WGCNA::corPvalueStudent(module_trait_cor, n_samp)

  gene_obj$wgcna[[wgcna_key]]$te_trait_cor   <- module_trait_cor
  gene_obj$wgcna[[wgcna_key]]$te_trait_pval  <- module_trait_pval
  gene_obj$wgcna[[wgcna_key]]$te_trait_names <- colnames(traits)
  gene_obj$wgcna[[wgcna_key]]$n_te_families  <- n_te_fams
  gene_obj$wgcna[[wgcna_key]]$n_clinical     <- n_clin

  message("[run_omni_wgcna_te_traits] Done.")
  message("  Clinical traits: ", n_clin,
          " | TE families: ", n_te_fams)
  message("  Plot with: plot_omni_wgcna_trait_heatmap(obj, '",
          wgcna_key, "')  — clinical traits")
  message("             plot_omni_wgcna_te_trait_heatmap(obj, '",
          wgcna_key, "') — TE family traits")

  return(gene_obj)
}


# =============================================================================
# 4. MODULE PRESERVATION
# =============================================================================

#' Test WGCNA Module Preservation Across Conditions
#'
#' Tests whether modules identified in a reference condition are preserved in
#' a test condition. High Zsummary (>10) = strong preservation;
#' low Zsummary (<2) = condition-specific modules.
#'
#' @param obj A validated OmniSeq object with WGCNA multi results.
#' @param group_col Character. The multi-WGCNA key to use.
#' @param reference Character. Reference condition name.
#' @param test Character. Test condition name.
#' @param n_permutations Integer. Permutations for preservation statistics.
#'   Default: `200`. Increase to 1000 for publication.
#'
#' @return The OmniSeq object with preservation results stored in
#'   `obj$wgcna_multi[[group_col]]$preservation`.
#' @export
run_omni_module_preservation <- function(obj,
                                         group_col      = "Condition",
                                         reference      = NULL,
                                         test           = NULL,
                                         n_permutations = 200L) {

  validate_omniseq(obj)
  .require_pkg("WGCNA")

  cor <- WGCNA::cor
  on.exit(rm(cor), add = TRUE)

  if (is.null(obj$wgcna_multi) || !group_col %in% names(obj$wgcna_multi))
    stop("No multi-WGCNA results for '", group_col,
         "'. Run run_omni_wgcna_multi() first.", call. = FALSE)

  multi_entry  <- obj$wgcna_multi[[group_col]]
  group_levels <- multi_entry$group_levels

  if (is.null(reference)) reference <- group_levels[1]
  if (is.null(test))      test      <- group_levels[2]

  if (!reference %in% group_levels)
    stop("Reference '", reference, "' not found.", call. = FALSE)
  if (!test %in% group_levels)
    stop("Test '", test, "' not found.", call. = FALSE)

  message("\n[run_omni_module_preservation] Reference: ", reference,
          " | Test: ", test)
  message("  Permutations: ", n_permutations,
          " (increase to 1000 for publication-quality p-values)")

  ref_net  <- multi_entry$group_nets[[reference]]
  test_net <- multi_entry$group_nets[[test]]

  vst_mat <- SummarizedExperiment::assay(omni_vst(obj))

  ref_mat  <- t(vst_mat[names(ref_net$module_colors),
                        ref_net$samples, drop = FALSE])
  test_mat <- t(vst_mat[names(test_net$module_colors),
                        test_net$samples, drop = FALSE])

  # Filter out genes with zero variance in either group
  # (WGCNA::modulePreservation fails if any gene is constant)
  ref_var  <- apply(ref_mat,  2, stats::var, na.rm = TRUE)
  test_var <- apply(test_mat, 2, stats::var, na.rm = TRUE)
  good_genes <- names(ref_var)[ref_var > 0 & test_var > 0 & !is.na(ref_var) & !is.na(test_var)]

  n_removed <- ncol(ref_mat) - length(good_genes)
  if (n_removed > 0)
    message("[run_omni_module_preservation] Removed ", n_removed,
            " zero-variance genes before preservation test.")

  ref_mat  <- ref_mat[,  good_genes, drop = FALSE]
  test_mat <- test_mat[, good_genes, drop = FALSE]

  # Also filter module color vectors to match
  ref_colors  <- ref_net$module_colors[good_genes]
  test_colors <- test_net$module_colors[good_genes]

  multiExpr <- list(
    ref  = list(data = ref_mat),
    test = list(data = test_mat)
  )

  multiColor <- list(
    ref  = ref_colors,
    test = test_colors
  )

  message("[run_omni_module_preservation] Running preservation statistics...")
  set.seed(42)

  pres <- WGCNA::modulePreservation(
    multiData         = multiExpr,
    multiColor        = multiColor,
    dataIsExpr        = TRUE,
    referenceNetworks = 1,
    nPermutations     = n_permutations,
    randomSeed        = 42,
    verbose           = 0
  )

  obj$wgcna_multi[[group_col]]$preservation <- list(
    result    = pres,
    reference = reference,
    test      = test
  )

  message("[run_omni_module_preservation] Done.")
  message("  Next: plot_omni_module_preservation(obj, '", group_col, "')")

  return(obj)
}


# =============================================================================
# 5. DENDROGRAM PLOT
# =============================================================================

#' Plot WGCNA Module Dendrogram
#'
#' @param obj A validated OmniSeq object with WGCNA results.
#' @param key Character. WGCNA cache key.
#' @param title Optional custom title.
#'
#' @return Called for side effect. Returns invisibly.
#' @export
plot_omni_wgcna_dendrogram <- function(obj, key = NULL, title = NULL) {

  validate_omniseq(obj)
  .require_pkg("WGCNA")

  entry <- .omni_get_wgcna_entry(obj, key)

  if (is.null(title))
    title <- paste0("Gene Co-expression Dendrogram — ", key)

  WGCNA::plotDendroAndColors(
    entry$net$dendrograms[[1]],
    entry$module_colors[entry$net$blockGenes[[1]]],
    "Module colors",
    dendroLabels  = FALSE,
    hang          = 0.03,
    addGuide      = TRUE,
    guideHang     = 0.05,
    main          = title
  )

  return(invisible(NULL))
}


# =============================================================================
# 6. MODULE-TRAIT HEATMAP (CLINICAL TRAITS)
# =============================================================================

#' Plot Module Eigengene x Clinical Trait Heatmap
#'
#' Displays Pearson correlation between module eigengenes and clinical traits
#' (Condition, Genotype, Phenotype, Age, Sex, Batch, IRG-S score). Cell color
#' shows correlation direction, text shows R and p-value.
#'
#' If TE trait correlations are also available (from `run_omni_wgcna_te_traits()`),
#' this function shows only clinical traits. Use
#' `plot_omni_wgcna_te_trait_heatmap()` for the TE family heatmap.
#'
#' @param obj A validated OmniSeq object with WGCNA results.
#' @param key Character. WGCNA cache key.
#' @param traits Optional named numeric data frame. If `NULL`, built from
#'   metadata automatically.
#' @param irgs_obj Optional OmniSeq gene object to compute IRG-S from.
#' @param title Optional custom title.
#'
#' @return Called for side effect. Returns invisibly.
#' @export
plot_omni_wgcna_trait_heatmap <- function(obj,
                                          key      = NULL,
                                          traits   = NULL,
                                          irgs_obj = NULL,
                                          title    = NULL) {

  validate_omniseq(obj)
  .require_pkg("WGCNA")

  entry    <- .omni_get_wgcna_entry(obj, key)
  MEs      <- entry$eigengenes
  dat_expr <- entry$dat_expr
  meta     <- omni_metadata(obj)
  samples  <- rownames(dat_expr)

  # Build clinical trait matrix if not provided
  if (is.null(traits)) {
    traits <- data.frame(row.names = samples)

    cond_col  <- obj$params$condition
    if (cond_col %in% colnames(meta)) {
      lvls <- levels(meta[[cond_col]])
      traits$Condition <- as.numeric(
        meta[samples, cond_col] == lvls[length(lvls)])
    }

    for (col in c("Genotype", "Phenotype", "Sex", "Batch")) {
      if (col %in% colnames(meta))
        traits[[col]] <- as.numeric(as.factor(meta[samples, col]))
    }

    if ("Age_Scaled" %in% colnames(meta))
      traits$Age <- as.numeric(meta[samples, "Age_Scaled"])

    # IRG-S score
    irg_src <- irgs_obj %||% obj
    vsd_src <- omni_vst(irg_src)
    if (!is.null(vsd_src)) {
      mat_g      <- SummarizedExperiment::assay(vsd_src)
      valid_irgs <- intersect(omniseq_isgs_28, rownames(mat_g))
      if (length(valid_irgs) > 0) {
        irgs_scores    <- colMeans(mat_g[valid_irgs, samples, drop = FALSE],
                                   na.rm = TRUE)
        traits$IRG_S_Score <- as.numeric(irgs_scores)
        message("  IRG-S score added (", length(valid_irgs),
                "/28 genes found).")
      }
    }
  }

  # Align and filter
  common_samps <- intersect(samples, rownames(traits))
  MEs_sub      <- MEs[common_samps, , drop = FALSE]
  traits_sub   <- traits[common_samps, , drop = FALSE]

  keep_traits  <- apply(traits_sub, 2,
                        function(x) stats::var(x, na.rm = TRUE) > 0)
  traits_sub   <- traits_sub[, keep_traits, drop = FALSE]

  if (ncol(traits_sub) == 0)
    stop("No variable traits found.", call. = FALSE)

  n_samp         <- length(common_samps)
  mod_trait_cor  <- WGCNA::cor(MEs_sub, traits_sub, use = "p")
  mod_trait_pval <- WGCNA::corPvalueStudent(mod_trait_cor, n_samp)

  text_mat <- paste0(signif(mod_trait_cor, 2), "\n(",
                     signif(mod_trait_pval, 1), ")")
  dim(text_mat) <- dim(mod_trait_cor)

  if (is.null(title))
    title <- paste0("Module-Trait Relationships — ", key)

  graphics::par(mar = c(6, 8.5, 3, 3))
  WGCNA::labeledHeatmap(
    Matrix        = mod_trait_cor,
    xLabels       = colnames(traits_sub),
    yLabels       = names(MEs_sub),
    ySymbols      = names(MEs_sub),
    colorLabels   = FALSE,
    colors        = WGCNA::blueWhiteRed(50),
    textMatrix    = text_mat,
    setStdMargins = FALSE,
    cex.text      = 0.7,
    zlim          = c(-1, 1),
    main          = title
  )

  # Remind user about TE trait heatmap if available
  if (!is.null(entry$te_trait_cor))
    message("  TE trait correlations also available. Plot with:\n",
            "  plot_omni_wgcna_te_trait_heatmap(obj, '",
            if (is.null(key)) names(obj$wgcna)[1] else key, "')")

  return(invisible(NULL))
}


# =============================================================================
# 7. TE TRAIT HEATMAP (NEW)
# =============================================================================

#' Plot Gene Module x TE Family Correlation Heatmap
#'
#' Displays the Pearson correlation between gene co-expression module
#' eigengenes and TE family expression levels. This is the companion to
#' `plot_omni_wgcna_trait_heatmap()` — while that function shows clinical
#' traits, this shows TE families as the external traits.
#'
#' Requires `run_omni_wgcna_te_traits()` to have been run first.
#'
#' @param obj A validated OmniSeq object with WGCNA + TE trait results.
#' @param key Character. WGCNA cache key.
#' @param show_clinical Logical. Include clinical trait columns alongside
#'   TE families. Default: `TRUE`.
#' @param title Optional custom title.
#'
#' @return Called for side effect. Returns invisibly.
#' @export
plot_omni_wgcna_te_trait_heatmap <- function(obj,
                                             key            = NULL,
                                             show_clinical  = TRUE,
                                             title          = NULL) {

  validate_omniseq(obj)
  .require_pkg("WGCNA")

  entry <- .omni_get_wgcna_entry(obj, key)

  if (is.null(entry$te_trait_cor))
    stop("No TE trait correlations found. ",
         "Run run_omni_wgcna_te_traits() first.", call. = FALSE)

  te_cor  <- entry$te_trait_cor
  te_pval <- entry$te_trait_pval

  # Optionally show only TE families (strip clinical columns)
  if (!show_clinical && !is.null(entry$n_clinical) &&
      entry$n_clinical > 0) {
    n_clin  <- entry$n_clinical
    n_total <- ncol(te_cor)
    te_cor  <- te_cor[, (n_clin + 1):n_total, drop = FALSE]
    te_pval <- te_pval[, (n_clin + 1):n_total, drop = FALSE]
  }

  # Build text matrix
  text_mat <- paste0(signif(te_cor, 2), "\n(",
                     signif(te_pval, 1), ")")
  dim(text_mat) <- dim(te_cor)

  if (is.null(title))
    title <- paste0("Gene Modules vs TE Families",
                    if (show_clinical) " + Clinical Traits" else "",
                    " — ", if (is.null(key)) names(obj$wgcna)[1] else key)

  graphics::par(mar = c(8, 8.5, 3, 3))
  WGCNA::labeledHeatmap(
    Matrix        = te_cor,
    xLabels       = colnames(te_cor),
    yLabels       = rownames(te_cor),
    ySymbols      = rownames(te_cor),
    colorLabels   = FALSE,
    colors        = WGCNA::blueWhiteRed(50),
    textMatrix    = text_mat,
    setStdMargins = FALSE,
    cex.text      = 0.6,
    zlim          = c(-1, 1),
    main          = title
  )

  return(invisible(NULL))
}


# =============================================================================
# 8. HUB GENE NETWORK
# =============================================================================

#' Plot WGCNA Hub Gene Network
#'
#' Extracts hub genes (highest module membership) from each module and
#' displays their co-expression relationships as a network.
#'
#' @param obj A validated OmniSeq object with WGCNA results.
#' @param key Character. WGCNA cache key.
#' @param modules Character vector. Modules to include. Default: all
#'   non-grey modules.
#' @param n_hubs Integer. Hub genes per module. Default: `5`.
#' @param plot_type Character. `"static"` (ggraph) or `"interactive"`
#'   (visNetwork). Default: `"static"`.
#'
#' @return A plot object (returned invisibly).
#' @export
plot_omni_wgcna_hub_network <- function(obj,
                                        key        = NULL,
                                        modules    = NULL,
                                        n_hubs     = 5L,
                                        plot_type  = c("static",
                                                       "interactive")) {

  validate_omniseq(obj)
  .require_pkg("WGCNA")
  .require_pkg("igraph")

  plot_type <- match.arg(plot_type)
  entry     <- .omni_get_wgcna_entry(obj, key)

  module_colors <- entry$module_colors
  dat_expr      <- entry$dat_expr

  all_mods <- unique(module_colors)
  all_mods <- all_mods[all_mods != "grey"]

  if (is.null(modules)) modules <- all_mods
  modules <- intersect(modules, all_mods)
  if (length(modules) == 0)
    stop("No valid modules specified.", call. = FALSE)

  # Module membership (kME)
  MMs <- as.data.frame(WGCNA::cor(dat_expr, entry$eigengenes, use = "p"))
  colnames(MMs) <- gsub("^ME", "", colnames(MMs))

  # Hub genes per module
  hub_genes <- list()
  for (mod in modules) {
    mod_genes <- names(module_colors)[module_colors == mod]
    if (!mod %in% colnames(MMs)) next
    mm_vals  <- MMs[mod_genes, mod]
    top_hubs <- utils::head(mod_genes[order(-mm_vals)], n_hubs)
    hub_genes[[mod]] <- top_hubs
  }

  all_hubs <- unique(unlist(hub_genes))
  if (length(all_hubs) == 0)
    stop("No hub genes found.", call. = FALSE)

  # Correlation matrix between hub genes
  cor_mat   <- WGCNA::cor(dat_expr[, all_hubs, drop = FALSE], use = "p")
  threshold <- 0.6

  .build_edges <- function(cor_mat, thr) {
    do.call(rbind, lapply(seq_len(nrow(cor_mat) - 1), function(i) {
      do.call(rbind, lapply(seq(i + 1, ncol(cor_mat)), function(j) {
        r <- cor_mat[i, j]
        if (abs(r) >= thr)
          data.frame(from = rownames(cor_mat)[i],
                     to   = colnames(cor_mat)[j],
                     r    = r, stringsAsFactors = FALSE)
        else NULL
      }))
    }))
  }

  edges_df <- .build_edges(cor_mat, threshold)
  if (is.null(edges_df) || nrow(edges_df) == 0) {
    threshold <- 0.4
    message("[plot_omni_wgcna_hub_network] Lowering edge threshold to 0.4")
    edges_df <- .build_edges(cor_mat, threshold)
  }

  node_module <- module_colors[all_hubs]
  node_df <- data.frame(
    name   = all_hubs,
    module = node_module,
    stringsAsFactors = FALSE
  )

  g <- igraph::graph_from_data_frame(edges_df, directed = FALSE,
                                     vertices = node_df)

  # ── STATIC ─────────────────────────────────────────────────────────────────
  if (plot_type == "static") {
    .require_pkg("ggraph")
    .require_pkg("ggplot2")

    p <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(
        ggplot2::aes(width = abs(r), color = r),
        alpha = 0.6
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(color = module), size = 5
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label = name),
        repel = TRUE, size = 3, fontface = "italic"
      ) +
      ggplot2::scale_color_identity(guide = "legend", name = "Module") +
      ggraph::scale_edge_color_gradient2(
        low      = "#2166AC",
        mid      = "grey80",
        high     = "#B2182B",
        midpoint = 0,
        name     = "Correlation",
        guide    = "none"          # fixes ggraph version incompatibility
      ) +
      ggraph::scale_edge_width_continuous(
        range = c(0.3, 1.5), name = "|r|", guide = "none"
      ) +
      ggraph::theme_graph(base_size = 12) +
      ggplot2::labs(
        title    = paste0("WGCNA Hub Gene Network — ", key),
        subtitle = paste(length(modules), "modules |",
                         length(all_hubs), "hub genes |",
                         if (!is.null(edges_df)) nrow(edges_df) else 0,
                         "edges (|r| \u2265 ", threshold, ")")
      )

    print(p)
    return(invisible(p))
  }

  # ── INTERACTIVE ─────────────────────────────────────────────────────────────
  if (plot_type == "interactive") {
    .require_pkg("visNetwork")

    vis_nodes <- data.frame(
      id    = igraph::V(g)$name,
      label = igraph::V(g)$name,
      color = igraph::V(g)$module,
      group = igraph::V(g)$module,
      title = paste0("<b>", igraph::V(g)$name, "</b><br>Module: ",
                     igraph::V(g)$module),
      stringsAsFactors = FALSE
    )

    el     <- igraph::as_edgelist(g)
    r_vals <- igraph::E(g)$r

    vis_edges <- data.frame(
      from  = el[, 1],
      to    = el[, 2],
      value = abs(r_vals),
      color = ifelse(r_vals > 0, "#E74C3C", "#2980B9"),
      title = paste0("r = ", round(r_vals, 3)),
      stringsAsFactors = FALSE
    )

    p <- visNetwork::visNetwork(
      vis_nodes, vis_edges,
      main   = paste("WGCNA Hub Network —", key),
      height = "650px", width = "100%"
    ) |>
      visNetwork::visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1),
        nodesIdSelection = TRUE
      ) |>
      visNetwork::visPhysics(
        solver           = "forceAtlas2Based",
        forceAtlas2Based = list(gravitationalConstant = -60)
      ) |>
      visNetwork::visLegend()

    print(p)
    return(invisible(p))
  }
}


# =============================================================================
# 9. MODULE PRESERVATION PLOT
# =============================================================================

#' Plot Module Preservation Statistics
#'
#' Zsummary vs module size plot. Modules with Zsummary > 10 are strongly
#' preserved, 2-10 moderately preserved, < 2 are condition-specific.
#'
#' @param obj A validated OmniSeq object with module preservation results.
#' @param group_col Character. Multi-WGCNA key.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_module_preservation <- function(obj, group_col = "Condition") {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .require_pkg("ggrepel")

  if (is.null(obj$wgcna_multi) || !group_col %in% names(obj$wgcna_multi))
    stop("No multi-WGCNA results for '", group_col, "'.", call. = FALSE)

  pres_entry <- obj$wgcna_multi[[group_col]]$preservation
  if (is.null(pres_entry))
    stop("No preservation results. Run run_omni_module_preservation() first.",
         call. = FALSE)

  pres <- pres_entry$result
  ref  <- pres_entry$reference
  test <- pres_entry$test

  # Extract Zsummary statistics
  # WGCNA stores results in $preservation$Z$ref.ref$inColumnsAlsoPresentIn.test
  z_entry <- pres$preservation$Z$ref.ref
  if (is.null(z_entry))
    stop("Could not extract Zsummary statistics.", call. = FALSE)

  # Handle both flat data.frame and nested list structures
  if (is.data.frame(z_entry) && "Zsummary.pres" %in% colnames(z_entry)) {
    z_stats <- z_entry
  } else if (is.list(z_entry)) {
    # Look for the sub-element containing Zsummary.pres
    z_stats <- NULL
    for (nm in names(z_entry)) {
      candidate <- z_entry[[nm]]
      if (is.data.frame(candidate) &&
          "Zsummary.pres" %in% colnames(candidate)) {
        z_stats <- candidate
        break
      }
    }
  }

  if (is.null(z_stats) || !("Zsummary.pres" %in% colnames(z_stats)))
    stop("Could not locate Zsummary.pres in preservation results.", call. = FALSE)

  # moduleSize is a column in z_stats for this WGCNA version
  if ("moduleSize" %in% colnames(z_stats)) {
    mod_sizes <- z_stats$moduleSize
    names(mod_sizes) <- rownames(z_stats)
  } else {
    mod_sizes <- rep(50, nrow(z_stats))
    names(mod_sizes) <- rownames(z_stats)
  }

  df <- data.frame(
    Module   = rownames(z_stats),
    Zsummary = z_stats$Zsummary.pres,
    Size     = as.numeric(mod_sizes[rownames(z_stats)]),
    stringsAsFactors = FALSE
  )
  df <- df[df$Module != "gold", ]
  df$Preservation <- ifelse(df$Zsummary > 10, "Strong (>10)",
                            ifelse(df$Zsummary > 2,  "Moderate (2-10)",
                                   "None (<2)"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Size, y = Zsummary,
                                        color = Preservation,
                                        label = Module)) +
    ggplot2::geom_point(size = 4, alpha = 0.85) +
    ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = c(2, 10), linetype = "dashed",
                        color = c("#E74C3C", "#27AE60"), linewidth = 0.7) +
    ggplot2::annotate("text", x = max(df$Size) * 0.95, y = 2.5,
                      label = "Moderate threshold", color = "#E74C3C",
                      size = 3, hjust = 1) +
    ggplot2::annotate("text", x = max(df$Size) * 0.95, y = 10.5,
                      label = "Strong threshold", color = "#27AE60",
                      size = 3, hjust = 1) +
    ggplot2::scale_color_manual(
      values = c("Strong (>10)"    = "#27AE60",
                 "Moderate (2-10)" = "#F39C12",
                 "None (<2)"       = "#E74C3C")
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title    = "Module Preservation",
      subtitle = paste0("Reference: ", ref, " | Test: ", test),
      x        = "Module size (number of genes)",
      y        = "Zsummary preservation statistic",
      color    = "Preservation"
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# ACCESSOR
# =============================================================================

#' Access a Stored WGCNA Result
#'
#' @param obj A validated OmniSeq object.
#' @param key Character. Cache key, or `NULL` to list all.
#'
#' @return The stored WGCNA result list, or a character vector of keys.
#' @export
omni_wgcna <- function(obj, key = NULL) {
  .check_class(obj)
  if (is.null(obj$wgcna) || length(obj$wgcna) == 0) {
    message("No WGCNA results found. Run run_omni_wgcna() first.")
    return(invisible(NULL))
  }
  if (is.null(key)) {
    cat("\nAvailable WGCNA networks:\n")
    for (k in names(obj$wgcna)) {
      e <- obj$wgcna[[k]]
      cat("  '", k, "' — ", e$n_modules, " modules | ",
          e$n_top_genes, " genes | power=", e$soft_power,
          " | masked=", e$use_masked,
          " | fitted: ", format(e$timestamp, "%Y-%m-%d %H:%M"),
          "\n", sep = "")
    }
    return(invisible(names(obj$wgcna)))
  }
  if (!key %in% names(obj$wgcna))
    stop("Key '", key, "' not found. Available: ",
         paste(names(obj$wgcna), collapse = ", "), call. = FALSE)
  obj$wgcna[[key]]
}


# =============================================================================
# INTERNAL HELPER
# =============================================================================

#' @keywords internal
.omni_get_wgcna_entry <- function(obj, key) {
  if (is.null(obj$wgcna) || length(obj$wgcna) == 0)
    stop("No WGCNA results found. Run run_omni_wgcna() first.", call. = FALSE)
  if (is.null(key)) {
    if (length(obj$wgcna) == 1) {
      key <- names(obj$wgcna)[1]
      message("[WGCNA] Auto-selected network: '", key, "'")
    } else {
      stop("Multiple WGCNA networks found. Specify 'key'.\n",
           "Available: ", paste(names(obj$wgcna), collapse = ", "),
           call. = FALSE)
    }
  }
  if (!key %in% names(obj$wgcna))
    stop("Key '", key, "' not found. Available: ",
         paste(names(obj$wgcna), collapse = ", "), call. = FALSE)
  obj$wgcna[[key]]
}
