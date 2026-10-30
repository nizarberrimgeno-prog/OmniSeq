# =============================================================================
# OmniSeq — 05_Enrichment.R  (Part 4A)
#
# Core enrichment engine — 12 functions:
#
# ANALYSIS:
#   run_omni_enrichment()          ORA with GO/KEGG + auto redundancy reduction
#   run_omni_gsea()                GSEA with msigdbr gene sets
#   run_omni_ssgsea()              Per-sample pathway scoring (GSVA)
#   run_omni_simplify()            Semantic similarity redundancy filter
#   run_omni_pathway_de()          limma DE on ssGSEA scores between groups
#
# PLOTTING:
#   plot_omni_ora()                Dotplot, cnetplot, emapplot
#   plot_omni_gsea()               Waterfall, ridge, NES barplot
#   plot_omni_ssgsea_heatmap()     Ranked per-sample pathway activity heatmap
#   plot_omni_compare_clusters()   Up vs Down side-by-side dotplot
#   plot_omni_leading_edge()       Leading edge gene heatmap per pathway
#   plot_omni_pathway_lfc()        Per-gene LFC barplot inside a pathway
#   plot_omni_multicontrast()      Bubble heatmap across multiple contrasts
# =============================================================================


# =============================================================================
# 1. ORA
# =============================================================================

#' Run Gene Ontology Over-Representation Analysis
#'
#' Extracts significant DEGs from the active contrast and runs GO enrichment
#' (Biological Process by default). Results are stored in `obj$enrichment`
#' and optionally simplified to remove redundant terms.
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param lfc_thresh Log2 fold change threshold. Default: `1`.
#' @param pval_thresh Adjusted p-value threshold. Default: `0.05`.
#' @param ont Character. GO ontology: `"BP"`, `"MF"`, `"CC"`, or `"ALL"`.
#' @param org_db Character. Annotation database. Inherited from `obj$params$org_db`.
#' @param simplify Logical. Remove redundant GO terms. Default: `TRUE`.
#' @param simplify_cutoff Numeric. Similarity cutoff. Default: `0.7`.
#' @param direction Character. `"up"`, `"down"`, or `"both"`. Default: `"both"`.
#'
#' @return The OmniSeq object with `$enrichment` populated.
#' @export
run_omni_enrichment <- function(obj,
                                lfc_thresh      = 0,
                                pval_thresh     = 0.05,
                                ont             = "BP",
                                org_db          = NULL,
                                simplify        = TRUE,
                                simplify_cutoff = 0.7,
                                direction       = c("both", "up", "down")) {

  validate_omniseq(obj)
  .require_pkg("clusterProfiler")

  direction <- match.arg(direction)
  org_db    <- org_db %||% obj$params$org_db
  .require_pkg(org_db)

  res <- omni_results(obj)
  if (is.null(res))
    stop("No DE results. Run run_omni_de() first.", call. = FALSE)

  message("\n[run_omni_enrichment] Running GO-ORA (", ont, ") for: ",
          obj$active_contrast)

  res_df      <- as.data.frame(res)
  res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)
  res_df$Gene <- rownames(res_df)

  sig <- res_df[res_df$padj < pval_thresh &
                  abs(res_df$log2FoldChange) > lfc_thresh, ]

  if (direction == "up")
    sig <- sig[sig$log2FoldChange > 0, ]
  else if (direction == "down")
    sig <- sig[sig$log2FoldChange < 0, ]

  if (nrow(sig) == 0)
    stop("No significant genes at these thresholds.", call. = FALSE)

  message("[run_omni_enrichment] Mapping ", nrow(sig), " genes to Entrez IDs...")

  gene_ids <- unique(sig$Gene)
  gene_ids <- gene_ids[nchar(gene_ids) <= 25 &
                         grepl("^[A-Za-z][A-Za-z0-9\\-\\./_]*$", gene_ids)]

  if (length(gene_ids) == 0)
    stop("No mappable gene symbols found in significant results.", call. = FALSE)

  message("[run_omni_enrichment] Running enrichGO with ", length(gene_ids),
          " gene symbols...")

  if (!requireNamespace(org_db, quietly = TRUE))
    stop("Package '", org_db, "' is required.\n",
         "Install with: BiocManager::install('", org_db, "')",
         call. = FALSE)
  suppressPackageStartupMessages(
    do.call(library, list(org_db, character.only = TRUE))
  )

  ego <- clusterProfiler::enrichGO(
    gene          = gene_ids,
    keyType       = "SYMBOL",
    OrgDb         = get(org_db),
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = FALSE
  )

  n_terms <- nrow(as.data.frame(ego))
  message("[run_omni_enrichment] Found ", n_terms, " enriched terms.")
  if (n_terms == 0)
    message("[run_omni_enrichment] No significant GO terms found at current thresholds.\n",
            "  Try run_omni_gsea() which is more powerful for smaller gene lists,\n",
            "  or relax thresholds with lfc_thresh = 0 and pval_thresh = 1.")

  if (simplify && n_terms > 0 && ont != "ALL") {
    message("[run_omni_enrichment] Simplifying redundant terms ",
            "(cutoff = ", simplify_cutoff, ")...")
    ego <- clusterProfiler::simplify(ego, cutoff = simplify_cutoff,
                                     by = "p.adjust", select_fun = min)
    message("[run_omni_enrichment] Retained ",
            nrow(as.data.frame(ego)), " terms after simplification.")
  }

  # Store in PER-CONTRAST CACHE keyed by active contrast
  if (is.null(obj$enrichment_cache)) obj$enrichment_cache <- list()
  obj$enrichment_cache[[obj$active_contrast]] <- ego

  # Keep legacy single slot for backward compat
  obj$enrichment                 <- ego
  obj$params$enrichment_contrast <- obj$active_contrast

  message("[run_omni_enrichment] Done. Access with obj$enrichment.")
  return(obj)
}


# =============================================================================
# 1B. ORA, UP/DOWN SPLIT, MULTIPLE DATABASES AT ONCE
# =============================================================================

#' Run Over-Representation Analysis, Split by Direction, Across Multiple Databases
#'
#' The workflow `run_omni_enrichment()` runs for one ontology/direction at a
#' time. This wraps that same significance logic but (1) always computes Up
#' and Down as two parallel clusters in one `compareCluster` call per
#' database (so they can be plotted/compared side by side), and (2) runs
#' GO-BP and MSigDB Hallmark ORA in the same call, since in practice these
#' are almost always wanted together. Each returned element is a
#' `compareClusterResult` ready to pass straight into
#' [plot_omni_enrichment_interactive()] for an interactive dot plot with a
#' linked, searchable gene table directly beneath it.
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param contrast Contrast name to pull from `results_cache`. Defaults to
#'   `obj$active_contrast`.
#' @param pval_thresh Adjusted p-value threshold for the input DEG list.
#'   Default `0.05`.
#' @param lfc_thresh Absolute log2FoldChange threshold for the input DEG
#'   list. Default `1`. Set to `0` to use DESeq2's padj cutoff alone (e.g.
#'   for a formal interaction coefficient, where no LFC filter is typically
#'   applied — see the package's interaction-test functions).
#' @param databases Character vector, any of `"GO-BP"`, `"Hallmark"`.
#'   Default both.
#' @param org_db Character. Annotation database for GO-BP. Inherited from
#'   `obj$params$org_db` if `NULL`.
#' @param species Character. Species name for `msigdbr` (Hallmark). Inherited
#'   from `obj$params$species` if `NULL`.
#' @param pvalue_cutoff Enrichment significance threshold (applied inside
#'   `enrichGO`/`enricher`). Default `0.05`.
#'
#' @return A named list, one element per requested database (`` `GO-BP` ``
#'   and/or `Hallmark`), each a `compareClusterResult` with `Cluster` levels
#'   `"Up"` and `"Down"`. An element is `NULL` (with a message, not an error)
#'   if that database found zero significant terms for this gene list — pass
#'   through this check before plotting.
#'
#' @examples
#' \dontrun{
#' res <- run_omni_ora_updown(obj, contrast = "WT_R848 vs WT_NS")
#' if (!is.null(res[["GO-BP"]]))
#'   plot_omni_enrichment_interactive(res[["GO-BP"]], title = "GO-BP -- WT_R848 vs WT_NS")
#' if (!is.null(res$Hallmark))
#'   plot_omni_enrichment_interactive(res$Hallmark, title = "Hallmark -- WT_R848 vs WT_NS")
#' }
#'
#' @export
run_omni_ora_updown <- function(obj,
                                contrast      = NULL,
                                pval_thresh   = 0.05,
                                lfc_thresh    = 1,
                                databases     = c("GO-BP", "Hallmark"),
                                org_db        = NULL,
                                species       = NULL,
                                pvalue_cutoff = 0.05) {

  validate_omniseq(obj)
  .require_pkg("clusterProfiler")
  databases <- match.arg(databases, c("GO-BP", "Hallmark"), several.ok = TRUE)

  key <- contrast %||% obj$active_contrast
  res <- omni_results(obj, contrast = key)
  if (is.null(res))
    stop("No DE results for contrast '", key, "'. Run run_omni_de()/",
         "run_omni_contrast() first.", call. = FALSE)

  d <- as.data.frame(res)
  d$Gene <- rownames(d)
  d$padj <- ifelse(is.na(d$padj), 1, d$padj)
  sig <- d[d$padj < pval_thresh & (lfc_thresh <= 0 | abs(d$log2FoldChange) > lfc_thresh), ]
  if (nrow(sig) == 0)
    stop("No significant genes at padj<", pval_thresh,
         if (lfc_thresh > 0) paste0(" & |log2FC|>", lfc_thresh) else "",
         " for contrast '", key, "'.", call. = FALSE)

  up_genes   <- sig$Gene[sig$log2FoldChange > 0]
  down_genes <- sig$Gene[sig$log2FoldChange < 0]
  message("[run_omni_ora_updown] ", key, " -- Up: ", length(up_genes),
          "  Down: ", length(down_genes))
  gene_clusters <- list(Up = up_genes, Down = down_genes)
  gene_clusters <- gene_clusters[vapply(gene_clusters, length, integer(1)) > 0]
  if (length(gene_clusters) == 0)
    stop("Neither direction has any genes to test.", call. = FALSE)

  out <- list()

  if ("GO-BP" %in% databases) {
    org_db_r <- org_db %||% obj$params$org_db
    .require_pkg(org_db_r)
    suppressPackageStartupMessages(do.call(library, list(org_db_r, character.only = TRUE)))
    cc <- tryCatch(
      clusterProfiler::compareCluster(
        geneCluster   = gene_clusters,
        ## Pass the function OBJECT, not the string "enrichGO" -- compareCluster
        ## resolves a string fun via a name lookup that requires the package to
        ## be attached (library()'d), which OmniSeq deliberately avoids doing
        ## to dependencies (requireNamespace only). A function object sidesteps
        ## that lookup entirely and works regardless of what's attached (this
        ## was caught failing for fun="enricher" below -- same risk applies
        ## here even where it happened to work).
        fun           = clusterProfiler::enrichGO,
        keyType       = "SYMBOL",
        OrgDb         = get(org_db_r),
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = pvalue_cutoff,
        readable      = FALSE
      ),
      error = function(e) { message("[run_omni_ora_updown] GO-BP failed: ", e$message); NULL }
    )
    if (is.null(cc) || nrow(as.data.frame(cc)) == 0) {
      message("[run_omni_ora_updown] GO-BP: no significant terms.")
      out[["GO-BP"]] <- NULL
    } else {
      out[["GO-BP"]] <- cc
    }
  }

  if ("Hallmark" %in% databases) {
    .require_pkg("msigdbr")
    species_r <- species %||% obj$params$species %||% "Homo sapiens"
    ## enricher() against Hallmark needs gene SYMBOLS, not Entrez -- using
    ## Entrez here silently returns zero hits against msigdbr's
    ## gene_symbol-keyed TERM2GENE (a real bug hit and fixed on a real
    ## project; gene_clusters are already symbols from `sig$Gene` above, so
    ## no conversion is done here -- keep it that way).
    hallmark_t2g <- msigdbr::msigdbr(species = species_r, collection = "H")[, c("gs_name", "gene_symbol")]
    cc <- tryCatch(
      clusterProfiler::compareCluster(
        geneCluster = gene_clusters,
        fun         = clusterProfiler::enricher,
        TERM2GENE   = hallmark_t2g,
        pvalueCutoff = pvalue_cutoff
      ),
      error = function(e) { message("[run_omni_ora_updown] Hallmark failed: ", e$message); NULL }
    )
    if (is.null(cc) || nrow(as.data.frame(cc)) == 0) {
      message("[run_omni_ora_updown] Hallmark: no significant terms.")
      out[["Hallmark"]] <- NULL
    } else {
      cc@compareClusterResult$Description <- gsub("^HALLMARK_", "", cc@compareClusterResult$Description)
      cc@compareClusterResult$Description <- gsub("_", " ", cc@compareClusterResult$Description)
      out[["Hallmark"]] <- cc
    }
  }

  out
}


# =============================================================================
# 2. GSEA
# =============================================================================

#' Run Gene Set Enrichment Analysis
#'
#' Uses the full ranked fold-change vector from the active contrast to run
#' GSEA against MSigDB gene sets.
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param species Character. Inherited from `obj$params$species`.
#' @param category Character. MSigDB collection. Default: `"H"`.
#' @param subcategory Character. Optional subcategory. Default: `NULL`.
#' @param pval_thresh Numeric. p-value cutoff. Default: `0.05`.
#' @param interactive Logical. Prompt user to select category.
#'
#' @return The OmniSeq object with `$advanced_gsea` populated.
#' @export
run_omni_gsea <- function(obj,
                          species     = NULL,
                          category    = "H",
                          subcategory = NULL,
                          pval_thresh = 0.05,
                          interactive = FALSE) {

  validate_omniseq(obj)
  .require_pkg("clusterProfiler")
  .require_pkg("msigdbr")

  species <- species %||% obj$params$species

  if (interactive && base::interactive()) {
    cat("\n[GSEA Setup]\n")
    cat("  H  = Hallmark (recommended first pass)\n")
    cat("  C2 = Curated (KEGG, Reactome, WikiPathways)\n")
    cat("  C5 = Ontology (GO BP, CC, MF)\n")
    cat("  C7 = Immunologic signatures\n\n")
    ans_cat <- readline("Category [default: H]: ")
    if (nzchar(trimws(ans_cat))) category <- toupper(trimws(ans_cat))
    if (category %in% c("C2", "C5", "C7")) {
      ans_sub <- readline("Subcategory (Enter to skip): ")
      if (nzchar(trimws(ans_sub))) subcategory <- toupper(trimws(ans_sub))
    }
  }

  message("\n[run_omni_gsea] Fetching MSigDB: species=", species,
          " | category=", category,
          " | subcategory=", subcategory %||% "ALL")

  m_df <- msigdbr::msigdbr(species     = species,
                           collection  = category,
                           subcategory = subcategory)

  if (nrow(m_df) == 0)
    stop("No gene sets returned. Check species/category/subcategory.",
         call. = FALSE)

  term2gene <- m_df[, c("gs_name", "gene_symbol")]

  fc <- .omni_fold_changes(obj)
  if (is.null(fc) || length(fc) == 0)
    stop("No fold changes found. Run run_omni_de() first.", call. = FALSE)

  message("[run_omni_gsea] Running GSEA on ", length(fc),
          " genes (", nrow(m_df), " gene sets)...")

  set.seed(42)
  gsea_res <- clusterProfiler::GSEA(
    geneList     = fc,
    TERM2GENE    = term2gene,
    pvalueCutoff = pval_thresh,
    verbose      = FALSE
  )

  n_paths <- nrow(as.data.frame(gsea_res))
  message("[run_omni_gsea] Found ", n_paths, " enriched gene sets.")

  # Store in PER-CONTRAST CACHE
  if (is.null(obj$gsea_cache)) obj$gsea_cache <- list()
  cache_key <- paste0(obj$active_contrast, "::GSEA::", category)
  obj$gsea_cache[[cache_key]]           <- gsea_res
  obj$gsea_cache[[obj$active_contrast]] <- gsea_res  # plain key for plots

  # Legacy backward compat
  obj$advanced_gsea           <- gsea_res
  obj$params$gsea_type        <- "GSEA"
  obj$params$gsea_category    <- category
  obj$params$gsea_subcategory <- subcategory
  obj$params$gsea_contrast    <- obj$active_contrast

  return(obj)
}


# =============================================================================
# 3. ssGSEA
# =============================================================================

#' Run Single-Sample GSEA (ssGSEA) via GSVA
#'
#' Computes per-sample pathway enrichment scores using the GSVA package.
#'
#' @param obj A validated OmniSeq object with a fitted model.
#' @param species Character. Inherited from `obj$params$species`.
#' @param category Character. MSigDB collection. Default: `"H"`.
#' @param subcategory Character. Optional subcategory.
#' @param interactive Logical. Prompt for parameters.
#'
#' @return The OmniSeq object with `$advanced_ssgsea` populated.
#' @export
run_omni_ssgsea <- function(obj,
                            species     = NULL,
                            category    = "H",
                            subcategory = NULL,
                            interactive = FALSE) {

  validate_omniseq(obj)
  .require_pkg("GSVA")
  .require_pkg("msigdbr")

  species <- species %||% obj$params$species

  if (interactive && base::interactive()) {
    cat("\n[ssGSEA Setup]\n")
    ans_cat <- readline("MSigDB Category [default: H]: ")
    if (nzchar(trimws(ans_cat))) category <- toupper(trimws(ans_cat))
    if (category %in% c("C2", "C5", "C7")) {
      ans_sub <- readline("Subcategory (Enter to skip): ")
      if (nzchar(trimws(ans_sub))) subcategory <- toupper(trimws(ans_sub))
    }
  }

  message("\n[run_omni_ssgsea] Fetching MSigDB: ", category,
          if (!is.null(subcategory)) paste0(" / ", subcategory))

  m_df <- msigdbr::msigdbr(species     = species,
                           collection  = category,
                           subcategory = subcategory)

  if (nrow(m_df) == 0)
    stop("No gene sets returned.", call. = FALSE)

  gs_list <- split(m_df$gene_symbol, m_df$gs_name)

  vsd <- omni_vst(obj)
  if (is.null(vsd))
    stop("VST not cached. Run run_omni_de() first.", call. = FALSE)

  mat <- SummarizedExperiment::assay(vsd)

  message("[run_omni_ssgsea] Running ssGSEA on ", ncol(mat),
          " samples x ", length(gs_list), " gene sets...")

  if (utils::packageVersion("GSVA") >= "1.50.0") {
    param         <- GSVA::ssgseaParam(mat, gs_list)
    ssgsea_scores <- GSVA::gsva(param, verbose = FALSE)
  } else {
    ssgsea_scores <- suppressWarnings(
      GSVA::gsva(mat, gs_list, method = "ssgsea",
                 kcdf = "Gaussian", verbose = FALSE)
    )
  }

  obj$advanced_ssgsea         <- ssgsea_scores
  obj$params$gsea_type        <- "ssGSEA"
  obj$params$gsea_category    <- category
  obj$params$gsea_subcategory <- subcategory

  message("[run_omni_ssgsea] Scored ", nrow(ssgsea_scores),
          " pathways across ", ncol(ssgsea_scores), " samples.")
  return(obj)
}


# =============================================================================
# 4. SIMPLIFY
# =============================================================================

#' Remove Redundant GO Terms by Semantic Similarity
#'
#' @param obj A validated OmniSeq object with enrichment results.
#' @param cutoff Numeric. Similarity cutoff. Default: `0.7`.
#'
#' @return The OmniSeq object with simplified `$enrichment`.
#' @export
run_omni_simplify <- function(obj, cutoff = 0.7) {

  validate_omniseq(obj)
  .require_pkg("clusterProfiler")

  ora <- .omni_get_enrichment(obj)
  if (is.null(ora))
    stop("No enrichment results. Run run_omni_enrichment() first.",
         call. = FALSE)

  n_before <- nrow(as.data.frame(ora))
  message("[run_omni_simplify] Simplifying ", n_before,
          " terms (cutoff = ", cutoff, ")...")

  ora_simplified <- clusterProfiler::simplify(
    ora,
    cutoff     = cutoff,
    by         = "p.adjust",
    select_fun = min
  )

  # Update both cache and legacy slot
  if (!is.null(obj$enrichment_cache) &&
      !is.null(obj$enrichment_cache[[obj$active_contrast]]))
    obj$enrichment_cache[[obj$active_contrast]] <- ora_simplified
  obj$enrichment <- ora_simplified

  n_after <- nrow(as.data.frame(ora_simplified))
  message("[run_omni_simplify] Retained ", n_after, " / ", n_before,
          " terms after simplification.")
  return(obj)
}


# =============================================================================
# 5. PATHWAY DE
# =============================================================================

#' Differential Pathway Activity Analysis
#'
#' Runs a limma linear model on ssGSEA scores to identify significantly
#' differentially active pathways between two groups.
#'
#' @param obj A validated OmniSeq object with ssGSEA scores.
#' @param group1 Character. Numerator group.
#' @param group2 Character. Denominator group.
#'
#' @return The OmniSeq object with `$pathway_results` populated.
#' @export
run_omni_pathway_de <- function(obj, group1 = NULL, group2 = NULL) {

  validate_omniseq(obj)
  .require_pkg("limma")

  if (is.null(obj$advanced_ssgsea))
    stop("No ssGSEA scores. Run run_omni_ssgsea() first.", call. = FALSE)

  meta     <- omni_metadata(obj)
  cond_col <- obj$params$condition
  groups   <- unique(as.character(meta[[cond_col]]))
  groups   <- groups[!is.na(groups)]

  if (is.null(group1) || is.null(group2)) {
    if (!base::interactive())
      stop("group1 and group2 must be specified in non-interactive mode.",
           call. = FALSE)
    cat("\nAvailable groups:", paste(groups, collapse = ", "), "\n")
    group1 <- readline("Group 1 (numerator) : ")
    group2 <- readline("Group 2 (baseline)  : ")
  }

  for (g in c(group1, group2))
    if (!g %in% groups)
      stop("Group '", g, "' not found. Available: ",
           paste(groups, collapse = ", "), call. = FALSE)

  message("[run_omni_pathway_de] Comparing: ", group1, " vs ", group2)

  g1_clean <- make.names(group1)
  g2_clean <- make.names(group2)

  fmla             <- stats::as.formula(paste("~ 0 +", cond_col))
  design           <- stats::model.matrix(fmla, data = meta)
  colnames(design) <- make.names(gsub(cond_col, "", colnames(design)))

  contrast_list <- list(Contrast = paste0(g1_clean, "-", g2_clean),
                        levels   = design)
  cm   <- do.call(limma::makeContrasts, contrast_list)
  fit  <- limma::lmFit(obj$advanced_ssgsea, design)
  fit2 <- limma::contrasts.fit(fit, cm)
  fit2 <- limma::eBayes(fit2)
  res  <- limma::topTable(fit2, coef = "Contrast", number = Inf)

  obj$pathway_results         <- res
  obj$params$pathway_contrast <- paste(group1, "vs", group2)

  message("[run_omni_pathway_de] Found ",
          sum(res$adj.P.Val < 0.05, na.rm = TRUE),
          " significant pathways (adj.P.Val < 0.05).")
  return(obj)
}


# =============================================================================
# 6. PLOT ORA
# =============================================================================

#' Plot ORA Enrichment Results
#'
#' @param obj A validated OmniSeq object with enrichment results.
#' @param plot_type Character. One of `"dot"`, `"cnet"`, `"emap"`.
#' @param show_n Integer. Number of top terms to display. Default: `15`.
#'
#' @return A ggplot/enrichplot object (returned invisibly).
#' @export
plot_omni_ora <- function(obj,
                          plot_type = c("dot", "cnet", "emap"),
                          show_n    = 15) {

  validate_omniseq(obj)
  .require_pkg("enrichplot")
  .require_pkg("ggplot2")

  plot_type <- match.arg(plot_type)

  ora <- .omni_get_enrichment(obj)
  if (is.null(ora))
    stop("No ORA results. Run run_omni_enrichment() first.", call. = FALSE)

  if (nrow(as.data.frame(ora)) == 0)
    stop("ORA result is empty — no significant terms found.", call. = FALSE)

  # Use full LFC vector (all genes, not just significant ones) for cnet coloring
  # This gives correct red=up / blue=down directionality across the full range
  fc <- .omni_fold_changes(obj)

  p <- switch(plot_type,

              "dot" = {
                enrichplot::dotplot(ora, showCategory = show_n) +
                  ggplot2::ggtitle("GO Over-Representation",
                                   subtitle = obj$active_contrast)
              },

              "cnet" = {
                .require_pkg("scales")
                fc_clamped <- scales::squish(fc, range = c(-3, 3))
                # colorEdge supported in enrichplot >= 1.14 — fall back gracefully
                p_cnet <- tryCatch(
                  enrichplot::cnetplot(ora,
                                       showCategory = show_n,
                                       foldChange   = fc_clamped,
                                       colorEdge    = TRUE),
                  error = function(e)
                    enrichplot::cnetplot(ora,
                                         showCategory = show_n,
                                         foldChange   = fc_clamped)
                )
                p_cnet +
                  ggplot2::scale_color_gradient2(
                    low      = "#2980B9",
                    mid      = "white",
                    high     = "#E74C3C",
                    midpoint = 0,
                    limits   = c(-3, 3),
                    oob      = scales::squish,
                    name     = "Log\u2082 FC"
                  ) +
                  ggplot2::ggtitle("Gene-Concept Network",
                                   subtitle = obj$active_contrast)
              },

              "emap" = {
                ora_sim <- enrichplot::pairwise_termsim(ora)
                enrichplot::emapplot(ora_sim, showCategory = show_n) +
                  ggplot2::ggtitle("Enrichment Map (Clustered Pathways)",
                                   subtitle = obj$active_contrast)
              }
  )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 7. PLOT GSEA
# =============================================================================

#' Plot GSEA Results
#'
#' @param obj A validated OmniSeq object with GSEA results.
#' @param plot_type Character. One of `"waterfall"`, `"ridge"`, `"bar"`.
#' @param show_n Integer. Number of pathways to show. Default: `15`.
#' @param pathway_idx Integer. For `"waterfall"`: which pathway rank to plot.
#'
#' @return A ggplot/enrichplot object (returned invisibly).
#' @export
plot_omni_gsea <- function(obj,
                           plot_type   = c("waterfall", "ridge", "bar"),
                           show_n      = 15,
                           pathway_idx = 1) {

  validate_omniseq(obj)
  .require_pkg("enrichplot")
  .require_pkg("ggplot2")

  plot_type <- match.arg(plot_type)

  gsea <- .omni_get_gsea(obj)
  if (is.null(gsea))
    stop("No GSEA results. Run run_omni_gsea() first.", call. = FALSE)

  gsea_df <- as.data.frame(gsea)
  if (nrow(gsea_df) == 0)
    stop("GSEA result is empty — no significant pathways.", call. = FALSE)

  gsea@result$Description <- gsub(
    "^HALLMARK_|^KEGG_|^REACTOME_|^WP_|^GOBP_|^GOMF_|^GOCC_", "",
    gsea@result$Description
  )
  gsea@result$Description <- gsub("_", " ", gsea@result$Description)

  p <- switch(plot_type,

              "waterfall" = {
                enrichplot::gseaplot2(gsea,
                                      geneSetID = pathway_idx,
                                      title     = gsea@result$Description[pathway_idx])
              },

              "ridge" = {
                enrichplot::ridgeplot(gsea, showCategory = show_n) +
                  ggplot2::labs(title    = "GSEA Ridge Plot",
                                subtitle = obj$active_contrast)
              },

              "bar" = {
                df             <- gsea_df[1:min(show_n, nrow(gsea_df)), ]
                df             <- df[order(df$NES), ]
                df$Description <- factor(df$Description, levels = df$Description)

                ggplot2::ggplot(df,
                                ggplot2::aes(x = NES, y = Description,
                                             fill = p.adjust)) +
                  ggplot2::geom_col() +
                  ggplot2::scale_fill_gradient(low = "#E74C3C", high = "#2980B9",
                                               name = "FDR") +
                  ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
                  ggplot2::theme_minimal(base_size = 12) +
                  ggplot2::labs(title    = "GSEA — Normalised Enrichment Scores",
                                subtitle = obj$active_contrast,
                                x = "NES", y = NULL)
              }
  )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 8. ssGSEA HEATMAP
# =============================================================================

#' Plot ssGSEA Pathway Activity Heatmap
#'
#' Rows are the most variable pathways, columns are samples, fill is the
#' ssGSEA enrichment score (row-scaled).
#'
#' @param obj A validated OmniSeq object with ssGSEA scores.
#' @param show_n Integer. Number of top variable pathways. Default: `30`.
#' @param isolate_contrast Logical. Subset samples to active contrast groups.
#'   Default: `TRUE`.
#' @param cluster_rows Logical. Cluster pathways. Default: `TRUE`.
#' @param cluster_cols Logical. Cluster samples. Default: `TRUE`.
#' @param color_up Color for high activity. Default: `"#B2182B"`.
#' @param color_down Color for low activity. Default: `"#2166AC"`.
#' @param title Optional custom title.
#'
#' @return A pheatmap object (returned invisibly).
#' @export
plot_omni_ssgsea_heatmap <- function(obj,
                                     show_n           = 30,
                                     isolate_contrast = TRUE,
                                     cluster_rows     = TRUE,
                                     cluster_cols     = TRUE,
                                     color_up         = "#B2182B",
                                     color_down       = "#2166AC",
                                     title            = NULL) {

  validate_omniseq(obj)
  .require_pkg("pheatmap")

  scores <- obj$advanced_ssgsea
  if (is.null(scores))
    stop("No ssGSEA scores. Run run_omni_ssgsea() first.", call. = FALSE)

  meta     <- omni_metadata(obj)
  cond_col <- obj$params$condition

  samples_to_keep <- if (isolate_contrast)
    .omni_subset_samples(obj, meta, cond_col)
  else
    rownames(meta)

  mat      <- scores[, samples_to_keep, drop = FALSE]
  meta_sub <- meta[samples_to_keep, , drop = FALSE]

  pvar      <- apply(mat, 1, stats::var)
  top_paths <- names(sort(pvar, decreasing = TRUE))[1:min(show_n, nrow(mat))]
  mat_plot  <- mat[top_paths, ]

  rownames(mat_plot) <- gsub(
    "^HALLMARK_|^KEGG_|^REACTOME_|^WP_|^GOBP_", "",
    rownames(mat_plot)
  )
  rownames(mat_plot) <- gsub("_", " ", rownames(mat_plot))

  annot_cols <- intersect(c(cond_col, obj$params$covariates),
                          colnames(meta_sub))
  df_annot   <- meta_sub[, annot_cols, drop = FALSE]
  df_annot[] <- lapply(df_annot,
                       function(x) if (is.factor(x)) droplevels(x) else x)

  pal <- grDevices::colorRampPalette(c(color_down, "white", color_up))(100)

  if (is.null(title))
    title <- paste0("ssGSEA Pathway Activity — ", obj$active_contrast)

  p <- pheatmap::pheatmap(
    mat_plot,
    annotation_col = df_annot,
    color          = pal,
    scale          = "row",
    border_color   = NA,
    cluster_rows   = cluster_rows,
    cluster_cols   = cluster_cols,
    show_rownames  = TRUE,
    fontsize_row   = if (show_n > 40) 5 else 8,
    main           = title
  )

  return(invisible(p))
}


# =============================================================================
# 9. COMPARE CLUSTERS
# =============================================================================

#' Compare Pathway Enrichment: Upregulated vs Downregulated Genes
#'
#' @param obj A validated OmniSeq object with DE results.
#' @param lfc_thresh Log2 fold change threshold. Default: `1`.
#' @param pval_thresh Adjusted p-value threshold. Default: `0.05`.
#' @param show_n Integer. Top pathways per cluster. Default: `10`.
#' @param org_db Character. Annotation database.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_compare_clusters <- function(obj,
                                       lfc_thresh  = 1,
                                       pval_thresh = 0.05,
                                       show_n      = 10,
                                       org_db      = NULL) {

  validate_omniseq(obj)
  .require_pkg("clusterProfiler")
  .require_pkg("enrichplot")
  .require_pkg("ggplot2")
  .require_pkg("stringr")

  org_db <- org_db %||% obj$params$org_db
  .require_pkg(org_db)

  res <- omni_results(obj)
  if (is.null(res))
    stop("No DE results. Run run_omni_de() first.", call. = FALSE)

  res_df      <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)

  up_genes   <- res_df$Gene[res_df$padj < pval_thresh &
                              res_df$log2FoldChange >  lfc_thresh]
  down_genes <- res_df$Gene[res_df$padj < pval_thresh &
                              res_df$log2FoldChange < -lfc_thresh]

  if (length(up_genes) == 0 && length(down_genes) == 0)
    stop("No significant genes found.", call. = FALSE)

  message("[plot_omni_compare_clusters] Up: ", length(up_genes),
          " | Down: ", length(down_genes))

  .to_entrez <- function(genes) {
    if (length(genes) == 0) return(character(0))
    map <- clusterProfiler::bitr(genes, fromType = "SYMBOL",
                                 toType = "ENTREZID", OrgDb = org_db)
    map$ENTREZID
  }

  gene_clusters <- list(
    Upregulated   = .to_entrez(up_genes),
    Downregulated = .to_entrez(down_genes)
  )
  gene_clusters <- gene_clusters[sapply(gene_clusters, length) > 0]

  message("[plot_omni_compare_clusters] Running compareCluster...")

  comp_res <- clusterProfiler::compareCluster(
    geneCluster  = gene_clusters,
    fun          = "enrichGO",
    OrgDb        = org_db,
    ont          = "BP",
    pvalueCutoff = 0.05,
    readable     = TRUE
  )

  if (is.null(comp_res) || nrow(as.data.frame(comp_res)) == 0)
    stop("No significant pathways found in either cluster.", call. = FALSE)

  p <- enrichplot::dotplot(comp_res, showCategory = show_n) +
    ggplot2::scale_y_discrete(
      labels = function(x) stringr::str_wrap(x, width = 50)
    ) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = if (show_n > 12) 7 else 9),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::ggtitle("Up vs Down — GO Biological Process",
                     subtitle = obj$active_contrast)

  print(p)
  return(invisible(p))
}


# =============================================================================
# 10. LEADING EDGE HEATMAP
# =============================================================================

#' Plot Leading Edge Gene Heatmap
#'
#' For each of the top N enriched GSEA pathways, extracts the leading edge
#' genes and displays their VST expression as an annotated heatmap.
#'
#' @param obj A validated OmniSeq object with GSEA results.
#' @param n_pathways Integer. Number of top pathways to extract. Default: `5`.
#' @param isolate_contrast Logical. Subset samples to contrast groups.
#' @param cluster_rows Logical. Cluster genes. Default: `FALSE`.
#' @param cluster_cols Logical. Cluster samples. Default: `TRUE`.
#' @param show_rownames Logical. Show gene names. Auto-detected if `NULL`:
#'   shown when <= 50 genes, hidden otherwise.
#' @param color_up Color for high expression. Default: `"#B2182B"`.
#' @param color_down Color for low expression. Default: `"#2166AC"`.
#' @param title Optional custom title.
#'
#' @return A pheatmap object (returned invisibly).
#' @export
plot_omni_leading_edge <- function(obj,
                                   n_pathways       = 5,
                                   isolate_contrast = TRUE,
                                   cluster_rows     = FALSE,
                                   cluster_cols     = TRUE,
                                   show_rownames    = NULL,
                                   color_up         = "#B2182B",
                                   color_down       = "#2166AC",
                                   title            = NULL) {

  validate_omniseq(obj)
  .require_pkg("pheatmap")

  gsea <- .omni_get_gsea(obj)
  vsd  <- omni_vst(obj)

  if (is.null(gsea))
    stop("No GSEA results. Run run_omni_gsea() first.", call. = FALSE)
  if (is.null(vsd))
    stop("VST not cached. Run run_omni_de() first.", call. = FALSE)

  gsea_df <- as.data.frame(gsea)
  if (nrow(gsea_df) == 0)
    stop("GSEA result is empty.", call. = FALSE)

  top_paths <- utils::head(gsea_df[order(gsea_df$p.adjust), ], n_pathways)

  leading_genes <- unique(unlist(
    lapply(top_paths$core_enrichment, function(x) strsplit(x, "/")[[1]])
  ))

  message("[plot_omni_leading_edge] Extracted ", length(leading_genes),
          " leading edge genes from top ", n_pathways, " pathways.")

  mat     <- SummarizedExperiment::assay(vsd)
  present <- intersect(leading_genes, rownames(mat))

  if (length(present) == 0)
    stop("No leading edge genes found in the VST matrix.", call. = FALSE)

  meta     <- omni_metadata(obj)
  cond_col <- obj$params$condition

  samples_to_keep <- if (isolate_contrast)
    .omni_subset_samples(obj, meta, cond_col)
  else
    rownames(meta)

  mat_plot <- mat[present, samples_to_keep, drop = FALSE]
  meta_sub <- meta[samples_to_keep, , drop = FALSE]

  annot_cols <- intersect(c(cond_col, obj$params$covariates),
                          colnames(meta_sub))
  df_annot   <- meta_sub[, annot_cols, drop = FALSE]

  gene_path_map <- lapply(seq_len(nrow(top_paths)), function(i) {
    genes <- strsplit(top_paths$core_enrichment[i], "/")[[1]]
    path  <- gsub("_", " ",
                  gsub("^HALLMARK_|^KEGG_|^REACTOME_|^GOBP_", "",
                       top_paths$Description[i]))
    data.frame(Gene = genes, Pathway = path, stringsAsFactors = FALSE)
  })
  gene_path_df           <- do.call(rbind, gene_path_map)
  gene_path_df           <- gene_path_df[!duplicated(gene_path_df$Gene), ]
  rownames(gene_path_df) <- gene_path_df$Gene

  df_row_annot <- data.frame(
    Pathway   = gene_path_df[present, "Pathway"],
    row.names = present
  )

  if (is.null(show_rownames))
    show_rownames <- length(present) <= 50

  pal <- grDevices::colorRampPalette(c(color_down, "white", color_up))(100)

  if (is.null(title))
    title <- paste0("Leading Edge Genes — Top ", n_pathways,
                    " GSEA Pathways\n", obj$active_contrast)

  p <- pheatmap::pheatmap(
    mat_plot,
    annotation_col = df_annot,
    annotation_row = df_row_annot,
    color          = pal,
    scale          = "row",
    border_color   = NA,
    cluster_rows   = cluster_rows,
    cluster_cols   = cluster_cols,
    show_rownames  = show_rownames,
    fontsize_row   = if (show_rownames) 7 else 0,
    main           = title
  )

  return(invisible(p))
}


# =============================================================================
# 11. PATHWAY LFC BARPLOT
# =============================================================================

#' Plot Per-Gene LFC Inside a Pathway
#'
#' For a given pathway, plots each gene's log2 fold change as a ranked
#' horizontal barplot colored by significance.
#'
#' @param obj A validated OmniSeq object with DE and GSEA results.
#' @param pathway_name Character. Pathway name (partial match allowed).
#'   Defaults to the top enriched pathway.
#' @param pval_thresh Adjusted p-value threshold for coloring. Default: `0.05`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_pathway_lfc <- function(obj,
                                  pathway_name = NULL,
                                  pval_thresh  = 0.05) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  gsea <- .omni_get_gsea(obj)
  res  <- omni_results(obj)

  if (is.null(gsea))
    stop("No GSEA results. Run run_omni_gsea() first.", call. = FALSE)
  if (is.null(res))
    stop("No DE results. Run run_omni_de() first.", call. = FALSE)

  gsea_df <- as.data.frame(gsea)

  if (is.null(pathway_name)) {
    pathway_name <- gsea_df$Description[1]
    message("[plot_omni_pathway_lfc] No pathway specified — using top pathway: ",
            pathway_name)
  }

  idx <- grep(pathway_name, gsea_df$Description, ignore.case = TRUE)[1]
  if (is.na(idx))
    stop("Pathway '", pathway_name, "' not found in GSEA results.",
         call. = FALSE)

  pathway_genes <- strsplit(gsea_df$core_enrichment[idx], "/")[[1]]

  res_df      <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)

  plot_df <- res_df[res_df$Gene %in% pathway_genes &
                      !is.na(res_df$log2FoldChange), ]

  if (nrow(plot_df) == 0)
    stop("No genes from this pathway found in DE results.", call. = FALSE)

  plot_df      <- plot_df[order(plot_df$log2FoldChange), ]
  plot_df$Gene <- factor(plot_df$Gene, levels = plot_df$Gene)
  plot_df$Sig  <- plot_df$padj < pval_thresh

  clean_name <- gsub("_", " ",
                     gsub("^HALLMARK_|^KEGG_|^REACTOME_|^GOBP_", "",
                          gsea_df$Description[idx]))

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = log2FoldChange,
                                    y = Gene,
                                    fill = Sig)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#E74C3C", "FALSE" = "grey70"),
      labels = c("TRUE" = paste0("padj < ", pval_thresh), "FALSE" = "NS"),
      name   = NULL
    ) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        size = if (nrow(plot_df) > 40) 6 else 9,
        face = "italic"
      )
    ) +
    ggplot2::labs(
      title    = clean_name,
      subtitle = paste("Contrast:", obj$active_contrast),
      x        = "Log\u2082 Fold Change",
      y        = NULL
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 12. MULTI-CONTRAST PATHWAY COMPARISON
# =============================================================================

#' Compare Pathway Enrichment Across Multiple Contrasts
#'
#' Runs ORA on each contrast in `results_cache` and displays a bubble heatmap
#' where rows are pathways, columns are contrasts, size is gene ratio and
#' color is -log10(FDR).
#'
#' @param obj A validated OmniSeq object with multiple contrasts.
#' @param contrasts Character vector of contrast names. Defaults to all.
#' @param lfc_thresh Log2 fold change threshold. Default: `1`.
#' @param pval_thresh Adjusted p-value threshold. Default: `0.05`.
#' @param show_n Integer. Top N pathways to display. Default: `20`.
#' @param org_db Character. Annotation database.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_multicontrast <- function(obj,
                                    contrasts   = NULL,
                                    lfc_thresh  = 1,
                                    pval_thresh = 0.05,
                                    show_n      = 20,
                                    org_db      = NULL) {

  validate_omniseq(obj)
  .require_pkg("clusterProfiler")
  .require_pkg("ggplot2")
  .require_pkg("stringr")

  org_db    <- org_db %||% obj$params$org_db
  .require_pkg(org_db)

  contrasts <- contrasts %||% omni_contrasts(obj)
  if (length(contrasts) < 2)
    stop("At least 2 contrasts required. Run run_omni_contrast() to add more.",
         call. = FALSE)

  message("[plot_omni_multicontrast] Running ORA for ",
          length(contrasts), " contrasts...")

  all_results <- lapply(contrasts, function(contrast_key) {

    res_entry <- obj$results_cache[[contrast_key]]
    if (is.null(res_entry)) return(NULL)

    res_df      <- as.data.frame(res_entry$results)
    res_df$Gene <- rownames(res_df)
    res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)

    sig_genes <- res_df$Gene[res_df$padj < pval_thresh &
                               abs(res_df$log2FoldChange) > lfc_thresh]
    if (length(sig_genes) == 0) return(NULL)

    gene_map <- tryCatch(
      clusterProfiler::bitr(sig_genes, fromType = "SYMBOL",
                            toType = "ENTREZID", OrgDb = org_db),
      error = function(e) NULL
    )
    if (is.null(gene_map) || nrow(gene_map) == 0) return(NULL)

    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene          = gene_map$ENTREZID,
        OrgDb         = org_db,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        readable      = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(ego)) return(NULL)

    df          <- as.data.frame(ego)
    if (nrow(df) == 0) return(NULL)
    df$Contrast <- contrast_key
    df
  })

  all_results <- Filter(Negate(is.null), all_results)
  if (length(all_results) == 0)
    stop("No enrichment results found for any contrast.", call. = FALSE)

  combined  <- do.call(rbind, all_results)
  min_padj  <- tapply(combined$p.adjust, combined$Description, min)
  top_terms <- names(sort(min_padj)[1:min(show_n, length(min_padj))])
  plot_df   <- combined[combined$Description %in% top_terms, ]

  plot_df$GeneRatioNum <- sapply(plot_df$GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    parts[1] / parts[2]
  })
  plot_df$LogFDR      <- -log10(plot_df$p.adjust + 1e-10)
  plot_df$Description <- stringr::str_wrap(plot_df$Description, width = 40)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = Contrast, y = Description,
                 size = GeneRatioNum, color = LogFDR)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "#FDBB84", high = "#B30000",
                                  name = "-log10(FDR)") +
    ggplot2::scale_size_continuous(name = "Gene ratio", range = c(2, 8)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      axis.text.y = ggplot2::element_text(size  = 8)
    ) +
    ggplot2::labs(
      title    = "Pathway Enrichment Across Contrasts",
      subtitle = paste("Top", length(top_terms), "GO-BP terms"),
      x        = NULL,
      y        = NULL
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Retrieve ORA results for the active contrast
#' Falls back to obj$enrichment for objects created before cache architecture
#' @keywords internal
.omni_get_enrichment <- function(obj, contrast = NULL) {
  key <- contrast %||% obj$active_contrast
  if (!is.null(obj$enrichment_cache) && !is.null(obj$enrichment_cache[[key]]))
    return(obj$enrichment_cache[[key]])
  obj$enrichment
}

#' Retrieve GSEA results for the active contrast
#' @keywords internal
.omni_get_gsea <- function(obj, contrast = NULL) {
  key <- contrast %||% obj$active_contrast
  if (!is.null(obj$gsea_cache) && !is.null(obj$gsea_cache[[key]]))
    return(obj$gsea_cache[[key]])
  obj$advanced_gsea
}

#' @keywords internal
.omni_fold_changes <- function(obj, contrast = NULL) {
  key   <- contrast %||% obj$active_contrast
  entry <- obj$results_cache[[key]]
  if (is.null(entry)) return(NULL)
  entry$fold_changes
}
