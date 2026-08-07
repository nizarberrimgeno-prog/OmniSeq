# =============================================================================
# OmniSeq — 14_Export.R
#
# Export engine: saves all results, plots, and tables to disk in a clean
# structured directory. Also generates a self-contained HTML report.
#
# FUNCTIONS:
#   export_omniseq()        Save RDS, DE tables, VST matrix, all plots as PNG
#   export_omniseq_report() Generate self-contained HTML report (rmarkdown)
#
# FIXES vs previous version:
#   - Every contrast (main / additional / masked) gets unique PNG filenames.
#     Masked plots are named e.g. volcano_masked_isgs_extended.png so they
#     never overwrite the unmasked versions.  (Issues 1, 4, 7)
#   - HTML report now contains full per-contrast sections: volcano, heatmap,
#     top genes, MA plot, TF activity reference, and sPLS-DA.  (Issues 3, 5, 6)
#   - sPLS-DA export loops over ALL stored models, not just the first key.
#     (Issue 6)
#   - Masked sections in HTML carry the mask label in both the section title
#     and every plot title so unmasked and masked are unambiguous.  (Issue 7)
#   - DIABLO tables and plots exported; DIABLO section added to HTML.  (Issue 8)
#   - gene_obj passed through to export and report for TE-gene correlations.
#
# FIXES applied in this version:
#   - has_tf (and all other feature detection flags) defined at the top of
#     export_omniseq(). Previous version used has_tf without defining it,
#     causing the entire TF plot export section to silently fail.
#   - Per-contrast TF activity plots exported for every key in
#     tf_contrast_cache (not just the main contrast).
#   - HTML report per-contrast TF sections wrapped in tryCatch so missing
#     TF results for a specific contrast don't break rendering.
#   - $$ Rmd bug fixed: all generated R code uses [["col"]] not $col
#     to prevent knitr doubling $ into $$ which crashes R parsing.
#   - TF genotype analysis iterates returned plot list — each panel saved
#     as a separate PNG (A through G) instead of one combined image.
#   - F panels (individual genotype profiles) and G panels (genotype-specific
#     TFs) added to both export and HTML report.
# =============================================================================


# =============================================================================
# SECTION 1: MAIN EXPORT FUNCTION
# =============================================================================

#' Export OmniSeq Results to Disk
#'
#' Saves all results from a completed OmniSeq analysis to a structured
#' directory. Creates subdirectories for tables (CSV) and plots (PNG).
#' Skips any section that has not been run without erroring.
#'
#' @param obj A validated OmniSeq object.
#' @param outdir Character. Path to output directory. Created if it
#'   does not exist. Default: `"./omniseq_output"`.
#' @param prefix Character. Optional prefix for all output files.
#'   Default: `NULL` (no prefix).
#' @param gene_obj Optional gene OmniSeq object for TE-gene correlation plots.
#' @param width Numeric. Plot width in inches. Default: `10`.
#' @param height Numeric. Plot height in inches. Default: `8`.
#' @param dpi Integer. Plot resolution. Default: `300`.
#' @param save_rds Logical. Save the full OmniSeq object as RDS.
#'   Default: `TRUE`.
#' @param overwrite Logical. Overwrite existing files. Default: `TRUE`.
#'
#' @return Invisibly returns a named list of paths to all saved files.
#' @export
export_omniseq <- function(obj,
                           outdir    = "./omniseq_output",
                           prefix    = NULL,
                           gene_obj  = NULL,
                           width     = 10,
                           height    = 8,
                           dpi       = 300,
                           save_rds  = TRUE,
                           overwrite = TRUE) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  outdir    <- path.expand(outdir)
  table_dir <- file.path(outdir, "tables")
  plot_dir  <- file.path(outdir, "plots")

  for (d in c(outdir, table_dir, plot_dir))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  pfx   <- if (!is.null(prefix)) paste0(prefix, "_") else ""
  saved <- list()

  message("\n[export_omniseq] Exporting to: ", outdir)
  message("  Prefix: ", if (is.null(prefix)) "(none)" else prefix)

  # ------------------------------------------------------------------
  # SECTION 1A: FEATURE DETECTION FLAGS
  # Defined once at the top so all sections can use them consistently.
  # ------------------------------------------------------------------
  has_de      <- !is.null(obj$results_cache) && length(obj$results_cache) > 0
  has_ora     <- !is.null(obj$enrichment) &&
    tryCatch(nrow(as.data.frame(obj$enrichment)) > 0,
             error = function(e) FALSE)
  has_gsea    <- !is.null(obj$advanced_gsea) &&
    tryCatch(nrow(as.data.frame(obj$advanced_gsea)) > 0,
             error = function(e) FALSE)
  has_ssgsea  <- !is.null(obj$advanced_ssgsea)
  has_tf      <- (!is.null(obj$tf_contrast_cache) &&
                    length(obj$tf_contrast_cache) > 0) ||
    !is.null(obj$tf_results)
  has_string  <- !is.null(obj$string_network)
  has_deconv  <- !is.null(obj$deconv)
  has_splsda  <- !is.null(obj$splsda) && length(obj$splsda) > 0
  has_wgcna   <- !is.null(obj$wgcna) && length(obj$wgcna) > 0
  has_te      <- !is.null(obj$te_family)
  has_te_corr <- has_te && !is.null(gene_obj)
  has_sensors <- !is.null(obj$sensor_scores) ||
    (!is.null(gene_obj) && !is.null(gene_obj$sensor_scores))
  has_diablo  <- !is.null(obj$diablo) && length(obj$diablo) > 0
  has_med     <- !is.null(obj$mediation) && length(obj$mediation) > 0
  has_pathway <- !is.null(obj$pathway_results)

  message("  Feature flags: DE=", has_de,
          " TF=", has_tf,
          " ORA=", has_ora,
          " GSEA=", has_gsea,
          " sPLSDA=", has_splsda,
          " WGCNA=", has_wgcna,
          " DIABLO=", has_diablo)

  # ------------------------------------------------------------------
  # SECTION 1B: HELPERS
  # ------------------------------------------------------------------
  .save_plot <- function(name, expr, w = width, h = height) {
    path <- file.path(plot_dir, paste0(pfx, name, ".png"))
    if (file.exists(path) && !overwrite) {
      message("  [skip] ", name, ".png (exists)")
      return(path)
    }
    tryCatch({
      grDevices::png(path, width = w, height = h, units = "in", res = dpi)
      force(expr)
      grDevices::dev.off()
      message("  [plot] ", name, ".png")
      saved[[name]] <<- path
      path
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      message("  [skip] ", name, ".png — ", e$message)
      NULL
    })
  }

  .save_csv <- function(name, df) {
    if (is.null(df) || (is.data.frame(df) && nrow(df) == 0)) {
      message("  [skip] ", name, ".csv — empty or NULL")
      return(NULL)
    }
    path <- file.path(table_dir, paste0(pfx, name, ".csv"))
    if (file.exists(path) && !overwrite) {
      message("  [skip] ", name, ".csv (exists)")
      return(path)
    }
    tryCatch({
      utils::write.csv(df, path, row.names = TRUE)
      message("  [csv]  ", name, ".csv (", nrow(df), " rows)")
      saved[[name]] <<- path
      path
    }, error = function(e) {
      message("  [skip] ", name, ".csv — ", e$message)
      NULL
    })
  }

  # ------------------------------------------------------------------
  # SECTION 1C: RDS
  # ------------------------------------------------------------------
  if (save_rds) {
    rds_path <- file.path(outdir, paste0(pfx, "omniseq_object.rds"))
    saveRDS(obj, rds_path)
    message("  [rds]  omniseq_object.rds")
    saved$rds <- rds_path
  }

  # ------------------------------------------------------------------
  # SECTION 1D: TABLES
  # ------------------------------------------------------------------
  message("\n[export_omniseq] Saving tables...")

  main_contrast <- obj$active_contrast

  # Main DE results
  res <- tryCatch(omni_results(obj), error = function(e) NULL)
  if (!is.null(res)) {
    res_df      <- as.data.frame(res)
    res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)
    .save_csv("de_results_full", res_df)
    .save_csv("de_results_sig", res_df[res_df$padj < 0.05, ])
  }

  # DE results for every additional contrast (including masked)
  all_contrast_keys <- names(obj$results_cache)
  for (ct in setdiff(all_contrast_keys, main_contrast)) {
    obj$active_contrast <- ct
    ct_res <- tryCatch(omni_results(obj), error = function(e) NULL)
    if (!is.null(ct_res)) {
      safe_ct  <- gsub("[^a-zA-Z0-9]", "_", ct)
      ct_df    <- as.data.frame(ct_res)
      ct_df$padj <- ifelse(is.na(ct_df$padj), 1, ct_df$padj)
      .save_csv(paste0("de_results_sig_", safe_ct),
                ct_df[ct_df$padj < 0.05, ])
    }
  }
  obj$active_contrast <- main_contrast

  # VST matrix
  vst <- tryCatch(omni_vst(obj), error = function(e) NULL)
  if (!is.null(vst)) {
    vst_mat <- as.data.frame(SummarizedExperiment::assay(vst))
    .save_csv("vst_matrix", vst_mat)
  }

  # ORA
  if (has_ora) {
    ora_df <- as.data.frame(obj$enrichment)
    if (nrow(ora_df) > 0) .save_csv("ora_terms", ora_df)
  }

  # GSEA
  if (has_gsea) {
    gsea_df <- as.data.frame(obj$advanced_gsea)
    if (nrow(gsea_df) > 0) .save_csv("gsea_results", gsea_df)
  }

  # Pathway DE
  if (has_pathway)
    .save_csv("pathway_de", obj$pathway_results)

  # TF activity — export all contrast scores as separate CSVs
  if (!is.null(obj$tf_contrast_cache) && length(obj$tf_contrast_cache) > 0) {
    for (ct_key in names(obj$tf_contrast_cache)) {
      tf_df <- tryCatch(as.data.frame(obj$tf_contrast_cache[[ct_key]]),
                        error = function(e) NULL)
      if (!is.null(tf_df) && nrow(tf_df) > 0) {
        safe_key <- gsub("[^a-zA-Z0-9]", "_", ct_key)
        .save_csv(paste0("tf_activity_", safe_key), tf_df)
      }
    }
  } else if (!is.null(obj$tf_results)) {
    # legacy fallback
    tf_wide <- tryCatch({
      tf_df  <- obj$tf_results
      tf_avg <- tapply(tf_df$score, list(tf_df$source, tf_df$condition),
                       mean, na.rm = TRUE)
      as.data.frame(tf_avg)
    }, error = function(e) NULL)
    if (!is.null(tf_wide)) .save_csv("tf_activity", tf_wide)
  }

  # Deconvolution
  if (has_deconv)
    .save_csv("deconv_fractions", as.data.frame(obj$deconv))

  # sPLS-DA signatures — all models
  if (has_splsda) {
    for (splsda_key in names(obj$splsda)) {
      sig <- tryCatch(
        omni_splsda_signature(obj, key = splsda_key, comp = 1),
        error = function(e) NULL
      )
      if (!is.null(sig)) {
        safe_sk <- gsub("[^a-zA-Z0-9]", "_", splsda_key)
        .save_csv(paste0("splsda_signature_", safe_sk), sig)
      }
    }
  }

  # WGCNA modules
  if (has_wgcna) {
    first_wgcna <- obj$wgcna[[1]]
    if (!is.null(first_wgcna$module_colors)) {
      mod_df <- data.frame(
        Gene   = names(first_wgcna$module_colors),
        Module = as.character(first_wgcna$module_colors),
        stringsAsFactors = FALSE
      )
      .save_csv("wgcna_modules", mod_df)
    }
    if (!is.null(first_wgcna$trait_cor))
      .save_csv("wgcna_trait_cor", as.data.frame(first_wgcna$trait_cor))
  }

  # Sensor scores
  if (!is.null(obj$sensor_scores))
    .save_csv("sensor_scores", obj$sensor_scores$scores)

  # Mediation
  if (has_med) {
    med_rows <- lapply(names(obj$mediation), function(nm) {
      s <- obj$mediation[[nm]]$summary
      if (is.null(s)) return(NULL)
      s$Mediator <- nm; s
    })
    med_all <- do.call(rbind, Filter(Negate(is.null), med_rows))
    if (!is.null(med_all)) .save_csv("mediation_summary", med_all)
  }

  # TE family / GREAT
  if (has_te)
    .save_csv("te_family_results", obj$te_family$by_family)
  if (!is.null(obj$te_great))
    .save_csv("te_great_results", obj$te_great$results)

  # DIABLO loadings
  if (has_diablo) {
    first_diablo_model <- obj$diablo[[1]]$model
    tryCatch({
      gene_load <- as.data.frame(first_diablo_model$loadings$genes)
      gene_load <- gene_load[rowSums(gene_load != 0) > 0, , drop = FALSE]
      .save_csv("diablo_gene_loadings", gene_load)
      te_load_d <- as.data.frame(first_diablo_model$loadings$TEs)
      te_load_d <- te_load_d[rowSums(te_load_d != 0) > 0, , drop = FALSE]
      .save_csv("diablo_te_loadings", te_load_d)
    }, error = function(e)
      message("  [skip] diablo loading tables — ", e$message))
  }

  # ------------------------------------------------------------------
  # SECTION 1E: PLOTS — QC
  # ------------------------------------------------------------------
  message("\n[export_omniseq] Saving plots...")

  .save_plot("qc_library_size",
             plot_omni_qc_basics(obj, plot_type = "library_size"))
  .save_plot("qc_detected_genes",
             plot_omni_qc_basics(obj, plot_type = "detected"))
  .save_plot("qc_distance",
             plot_omni_qc_advanced(obj, plot_type = "distance"))
  .save_plot("qc_pca_compare",
             plot_omni_qc_advanced(obj, plot_type = "pca_compare"))
  .save_plot("qc_mean_var",
             plot_omni_qc_advanced(obj, plot_type = "mean_var"))

  # ------------------------------------------------------------------
  # SECTION 1F: PLOTS — per contrast (main + additional + masked)
  # ------------------------------------------------------------------
  all_contrasts     <- names(obj$results_cache)
  non_masked        <- all_contrasts[!grepl("excl\\.", all_contrasts)]
  masked_contrasts  <- all_contrasts[ grepl("excl\\.", all_contrasts)]

  # Main contrast
  obj$active_contrast <- main_contrast
  .save_plot("pca",           plot_omni_pca(obj))
  .save_plot("volcano",
             plot_omni_volcano(obj, lfc_thresh = 1, top_n = 15))
  .save_plot("ma_plot",       plot_omni_ma(obj))
  .save_plot("heatmap_top50", plot_omni_heatmap(obj, ntop = 50))
  .save_plot("top_genes",     plot_omni_top_genes(obj, n = 12))
  .save_plot("pvalue_hist",   plot_omni_pvalue_hist(obj))

  # Additional non-masked contrasts
  for (ct in setdiff(non_masked, main_contrast)) {
    obj$active_contrast <- ct
    safe_ct <- gsub("[^a-zA-Z0-9]", "_", ct)
    .save_plot(paste0("volcano_",   safe_ct),
               plot_omni_volcano(obj, lfc_thresh = 1, top_n = 15))
    .save_plot(paste0("heatmap_",   safe_ct),
               plot_omni_heatmap(obj, ntop = 50))
    .save_plot(paste0("top_genes_", safe_ct),
               plot_omni_top_genes(obj, n = 12))
    .save_plot(paste0("ma_plot_",   safe_ct),
               plot_omni_ma(obj))
  }

  # Masked contrasts
  for (mk in masked_contrasts) {
    obj$active_contrast <- mk
    safe_mk <- gsub("[^a-zA-Z0-9]", "_", mk)
    .save_plot(paste0("volcano_",   safe_mk),
               plot_omni_volcano(obj, lfc_thresh = 1, top_n = 15))
    .save_plot(paste0("heatmap_",   safe_mk),
               plot_omni_heatmap(obj, ntop = 50))
    .save_plot(paste0("top_genes_", safe_mk),
               plot_omni_top_genes(obj, n = 12))
    .save_plot(paste0("ma_plot_",   safe_mk),
               plot_omni_ma(obj))
    if (has_ora) {
      .save_plot(paste0("ora_dotplot_", safe_mk),
                 plot_omni_ora(obj, plot_type = "dot"))
      .save_plot(paste0("ora_emapplot_", safe_mk),
                 plot_omni_ora(obj, plot_type = "emap"))
    }
    if (has_gsea)
      .save_plot(paste0("gsea_waterfall_", safe_mk),
                 plot_omni_gsea(obj, show_n = 10))
  }

  obj$active_contrast <- main_contrast  # restore

  # ------------------------------------------------------------------
  # SECTION 1G: PLOTS — enrichment (main contrast)
  # ------------------------------------------------------------------
  if (has_ora) {
    .save_plot("ora_dotplot",
               plot_omni_ora(obj, plot_type = "dot"))
    .save_plot("ora_cnetplot",
               plot_omni_ora(obj, plot_type = "cnet"))
    .save_plot("ora_emapplot",
               plot_omni_ora(obj, plot_type = "emap"))
    .save_plot("ora_compare_clusters",
               plot_omni_compare_clusters(obj))
  }
  if (has_gsea)
    .save_plot("gsea_waterfall", plot_omni_gsea(obj, show_n = 10))
  if (has_ssgsea)
    .save_plot("ssgsea_heatmap", plot_omni_ssgsea_heatmap(obj))
  if (has_pathway)
    .save_plot("pathway_lfc", plot_omni_pathway_lfc(obj))

  # ------------------------------------------------------------------
  # SECTION 1H: PLOTS — TF (all contrasts + genotype panels)
  # ── FIX: iterates returned list from plot_omni_tf_genotype()
  # ── FIX: each genotype panel (A-G) saved as separate PNG
  # ── FIX: per-contrast TF activity + network for all contrasts
  # ------------------------------------------------------------------
  if (has_tf) {
    # Main contrast TF plots
    obj$active_contrast <- main_contrast
    .save_plot("tf_activity",
               plot_omni_tf_activity(obj, n = 15),
               w = 10, h = 10)
    .save_plot("tf_network",
               plot_omni_tf_network(obj, plot_type = "static"),
               w = 12, h = 10)

    # TF genotype analysis — returns named list of independent plots
    # Each plot (A through G) saved as a separate PNG
    tf_geno_plots <- tryCatch(
      plot_omni_tf_genotype(obj),
      error = function(e) {
        message("  [skip] tf_genotype — ", e$message)
        NULL
      }
    )

    if (!is.null(tf_geno_plots) && length(tf_geno_plots) > 0) {
      for (tg_name in names(tf_geno_plots)) {
        tg_plot <- tf_geno_plots[[tg_name]]
        if (is.null(tg_plot)) next

        # Determine dimensions based on plot type
        tg_w <- if (grepl("combined", tg_name)) 24
        else if (grepl("^E_", tg_name)) 12
        else 10
        tg_h <- if (grepl("combined", tg_name)) 10
        else if (grepl("^E_", tg_name)) 14
        else if (grepl("^A_", tg_name)) 6
        else 10

        .save_plot(paste0("tf_genotype_", tg_name),
                   print(tg_plot),
                   w = tg_w, h = tg_h)
      }
      message("  [tf_genotype] Saved ", length(tf_geno_plots),
              " independent genotype plots")
    }

    # Per-contrast TF activity and network plots
    if (!is.null(obj$tf_contrast_cache)) {
      tf_ct_keys <- setdiff(names(obj$tf_contrast_cache), main_contrast)
      for (tf_ct in tf_ct_keys) {
        obj$active_contrast <- tf_ct
        safe_tf_ct <- gsub("[^a-zA-Z0-9]", "_", tf_ct)
        .save_plot(paste0("tf_activity_", safe_tf_ct),
                   tryCatch(
                     plot_omni_tf_activity(obj, n = 15),
                     error = function(e) {
                       message("  [skip] tf_activity_", safe_tf_ct,
                               ".png — ", e$message)
                       NULL
                     }
                   ),
                   w = 10, h = 10)
        .save_plot(paste0("tf_network_", safe_tf_ct),
                   tryCatch(
                     plot_omni_tf_network(obj, plot_type = "static"),
                     error = function(e) {
                       message("  [skip] tf_network_", safe_tf_ct,
                               ".png — ", e$message)
                       NULL
                     }
                   ),
                   w = 12, h = 10)
      }
      obj$active_contrast <- main_contrast  # restore
    }
  }

  # ------------------------------------------------------------------
  # SECTION 1I: PLOTS — Deconvolution
  # ------------------------------------------------------------------
  if (has_deconv)
    .save_plot("deconv_heatmap", plot_omni_deconv_heatmap(obj))

  # ------------------------------------------------------------------
  # SECTION 1J: PLOTS — sPLS-DA (all models)
  # ------------------------------------------------------------------
  if (has_splsda) {
    for (splsda_key in names(obj$splsda)) {
      safe_sk <- gsub("[^a-zA-Z0-9]", "_", splsda_key)
      .save_plot(paste0("splsda_scores_",  safe_sk),
                 plot_omni_splsda_scores(obj, key = splsda_key))
      .save_plot(paste0("splsda_vip_",     safe_sk),
                 plot_omni_splsda_vip(obj, key = splsda_key))
      .save_plot(paste0("splsda_roc_",     safe_sk),
                 plot_omni_splsda_roc(obj, key = splsda_key, comp = 1))
    }
  }

  # ------------------------------------------------------------------
  # SECTION 1K: PLOTS — WGCNA (all keys)
  # ------------------------------------------------------------------
  if (has_wgcna) {
    for (wgcna_key in names(obj$wgcna)) {
      safe_wk <- gsub("[^a-zA-Z0-9]", "_", wgcna_key)
      .save_plot(paste0("wgcna_dendrogram_",  safe_wk),
                 plot_omni_wgcna_dendrogram(obj, key = wgcna_key),
                 w = 12, h = 7)
      .save_plot(paste0("wgcna_traits_",      safe_wk),
                 plot_omni_wgcna_trait_heatmap(obj, key = wgcna_key),
                 w = 12, h = 9)
      .save_plot(paste0("wgcna_hub_network_", safe_wk),
                 plot_omni_wgcna_hub_network(obj, key = wgcna_key,
                                             plot_type = "static"))
    }
  }

  # ------------------------------------------------------------------
  # SECTION 1I2: PLOTS — SC Reference Deconvolution
  # ------------------------------------------------------------------
  has_sc_deconv <- !is.null(obj$sc_deconv) && !is.null(obj$sc_deconv$scores)
  has_disease_states <- has_sc_deconv &&
    !is.null(obj$sc_deconv$disease_states) &&
    !is.null(obj$sc_deconv$disease_states$scores)

  if (has_sc_deconv) {
    # Layer 1: Composition
    .save_plot("sc_deconv_heatmap_composition",
               plot_omni_sc_heatmap(obj, layer = "composition"),
               w = 12, h = 8)
    .save_plot("sc_deconv_genotype_composition",
               plot_omni_sc_genotype(obj, layer = "composition"),
               w = 14, h = 10)
    .save_plot("sc_deconv_correlation_irgs",
               plot_omni_sc_correlation(obj, feature = "IRGS",
                                        layer = "composition"),
               w = 10, h = 8)
    .save_plot("sc_deconv_profile_composition",
               plot_omni_sc_profile(obj, layer = "composition"),
               w = 10, h = 7)

    # Layer 2: Activation (if available)
    if (has_disease_states) {
      .save_plot("sc_deconv_heatmap_activation",
                 plot_omni_sc_heatmap(obj, layer = "activation"),
                 w = 12, h = 8)
      .save_plot("sc_deconv_genotype_activation",
                 plot_omni_sc_genotype(obj, layer = "activation"),
                 w = 14, h = 10)
      .save_plot("sc_deconv_profile_activation",
                 plot_omni_sc_profile(obj, layer = "activation"),
                 w = 10, h = 7)
      .save_plot("sc_deconv_two_layer",
                 plot_omni_sc_two_layer(obj),
                 w = 14, h = 10)
    }

    # Benchmark vs xCell
    if (!is.null(obj$deconv))
      .save_plot("sc_deconv_benchmark",
                 plot_omni_sc_benchmark(obj),
                 w = 12, h = 8)

    # Tables
    .save_csv("sc_deconv_composition_scores",
              as.data.frame(obj$sc_deconv$scores))
    if (!is.null(obj$sc_deconv$cor_with_irgs))
      .save_csv("sc_deconv_irgs_correlation",
                obj$sc_deconv$cor_with_irgs)
    if (has_disease_states) {
      .save_csv("sc_deconv_activation_scores",
                as.data.frame(obj$sc_deconv$disease_states$scores))
      if (!is.null(obj$sc_deconv$disease_states$cor_with_irgs))
        .save_csv("sc_deconv_activation_irgs_correlation",
                  obj$sc_deconv$disease_states$cor_with_irgs)
    }
  }

  # ------------------------------------------------------------------
  # SECTION 1L: PLOTS — TE
  # ------------------------------------------------------------------
  if (has_te) {
    .save_plot("te_family_bar", plot_omni_te_family(obj))
    .save_plot("te_volcano",    plot_omni_te_volcano(obj))
    .save_plot("te_heatmap",    plot_omni_te_heatmap(obj))
  }
  if (!is.null(obj$te_great))
    .save_plot("te_great", plot_omni_te_great(obj))

  # TE-gene correlation
  if (has_te_corr) {
    tryCatch({
      .save_plot("te_gene_corr_irgs",
                 plot_omni_te_gene_corr(
                   gene_obj = gene_obj, te_obj = obj,
                   marker_genes = OmniSeq::omniseq_isgs_28,
                   te_class = "ALL",
                   score_name = "IRG-S Score"))
      g_res      <- as.data.frame(omni_results(gene_obj))
      g_res$padj <- ifelse(is.na(g_res$padj), 1, g_res$padj)
      non_isg    <- setdiff(
        rownames(g_res)[g_res$padj < 0.05],
        OmniSeq::omniseq_isgs_28)
      if (length(non_isg) >= 3)
        .save_plot("te_gene_corr_nonisg",
                   plot_omni_te_gene_corr(
                     gene_obj = gene_obj, te_obj = obj,
                     marker_genes = non_isg,
                     te_class = "ALL",
                     score_name = "ISG-independent score"))
    }, error = function(e)
      message("  [skip] te_gene_corr — ", e$message))
  }

  # ------------------------------------------------------------------
  # SECTION 1M: PLOTS — Sensors
  # Sensor scores live on gene_obj (not the TE obj) — check both
  # ------------------------------------------------------------------
  sensor_source <- if (!is.null(obj$sensor_scores)) obj
  else if (!is.null(gene_obj) &&
           !is.null(gene_obj$sensor_scores)) gene_obj
  else NULL

  if (!is.null(sensor_source)) {
    .save_plot("sensor_scores",
               plot_omni_sensor_scores(sensor_source))
    .save_plot("sensor_irgs_corr",
               plot_omni_sensor_irgs_corr(sensor_source))
    .save_plot("sensor_heatmap",
               plot_omni_sensor_heatmap(sensor_source))
  }

  # ------------------------------------------------------------------
  # SECTION 1N: PLOTS — DIABLO
  # ------------------------------------------------------------------
  if (has_diablo) {
    first_diablo_key   <- names(obj$diablo)[1]
    first_diablo_model <- obj$diablo[[first_diablo_key]]$model
    .save_plot("diablo_scores",
               mixOmics::plotIndiv(first_diablo_model, legend = TRUE,
                                   title = "DIABLO — Sample Plot"))
    .save_plot("diablo_overview",
               mixOmics::plotDiablo(first_diablo_model, ncomp = 1))
    tryCatch(
      .save_plot("diablo_circos",
                 mixOmics::circosPlot(first_diablo_model, cutoff = 0.7,
                                      line = TRUE, size.labels = 1.5),
                 w = 10, h = 10),
      error = function(e)
        message("  [skip] diablo_circos.png — ", e$message))
    .save_plot("diablo_loadings",
               mixOmics::plotLoadings(first_diablo_model, comp = 1,
                                      contrib = "max", method = "median"))
  }

  # ------------------------------------------------------------------
  # SECTION 1O: PLOTS — Mediation
  # ------------------------------------------------------------------
  if (has_med) {
    .save_plot("mediation_forest", plot_omni_mediation(obj))
    first_med <- names(obj$mediation)[1]
    .save_plot(paste0("mediation_path_", first_med),
               plot_omni_mediation_path(obj, mediator = first_med))
  }

  # ------------------------------------------------------------------
  # SECTION 1P: SUMMARY
  # ------------------------------------------------------------------
  n_csv <- length(list.files(table_dir, pattern = "\\.csv$"))
  n_png <- length(list.files(plot_dir,  pattern = "\\.png$"))

  message("\n[export_omniseq] Export complete.")
  message("  Tables : ", n_csv, " CSV files -> ", table_dir)
  message("  Plots  : ", n_png, " PNG files -> ", plot_dir)
  if (save_rds)
    message("  Object : omniseq_object.rds -> ", outdir)

  return(invisible(saved))
}


# =============================================================================
# SECTION 2: HTML REPORT
# =============================================================================

#' Generate an HTML Report from an OmniSeq Object
#'
#' Renders a self-contained HTML report summarising all completed analyses.
#' Requires the `rmarkdown` and `knitr` packages.
#'
#' @param obj A validated OmniSeq object.
#' @param outdir Character. Output directory for the report.
#'   Default: `"./omniseq_output"`.
#' @param report_file Character. Output filename. Default:
#'   `"omniseq_report.html"`.
#' @param title Character. Report title.
#' @param author Character. Author name for the report header.
#' @param gene_obj Optional gene OmniSeq object for TE-gene correlation.
#' @param open Logical. Open the report in the browser after rendering.
#'   Default: `TRUE`.
#'
#' @return Invisibly returns the path to the HTML file.
#' @export
export_omniseq_report <- function(obj,
                                  outdir      = "./omniseq_output",
                                  report_file = "omniseq_report.html",
                                  title       = "OmniSeq Analysis Report",
                                  author      = "",
                                  gene_obj    = NULL,
                                  open        = TRUE) {

  validate_omniseq(obj)
  .require_pkg("rmarkdown")
  .require_pkg("knitr")

  outdir <- path.expand(outdir)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  html_path <- file.path(outdir, report_file)

  message("\n[export_omniseq_report] Generating HTML report...")
  message("  Output: ", html_path)

  # ------------------------------------------------------------------
  # SECTION 2A: Detect what has been run
  # ------------------------------------------------------------------
  has_de      <- !is.null(obj$results_cache) && length(obj$results_cache) > 0
  has_ora     <- !is.null(obj$enrichment) &&
    tryCatch(nrow(as.data.frame(obj$enrichment)) > 0,
             error = function(e) FALSE)
  has_gsea    <- !is.null(obj$advanced_gsea) &&
    tryCatch(nrow(as.data.frame(obj$advanced_gsea)) > 0,
             error = function(e) FALSE)
  has_ssgsea  <- !is.null(obj$advanced_ssgsea)
  has_tf      <- (!is.null(obj$tf_contrast_cache) &&
                    length(obj$tf_contrast_cache) > 0) ||
    !is.null(obj$tf_results)
  has_string  <- !is.null(obj$string_network)
  has_deconv  <- !is.null(obj$deconv)
  has_splsda  <- !is.null(obj$splsda) && length(obj$splsda) > 0
  has_wgcna   <- !is.null(obj$wgcna) && length(obj$wgcna) > 0
  has_te      <- !is.null(obj$te_family)
  has_te_corr <- has_te && !is.null(gene_obj)
  has_sensors <- !is.null(obj$sensor_scores) ||
    (!is.null(gene_obj) && !is.null(gene_obj$sensor_scores))
  sensor_obj_name <- if (!is.null(obj$sensor_scores)) "obj" else "gene_obj"
  has_diablo  <- !is.null(obj$diablo) && length(obj$diablo) > 0
  has_med     <- !is.null(obj$mediation) && length(obj$mediation) > 0

  main_contrast    <- obj$active_contrast
  all_ct_keys      <- names(obj$results_cache)
  non_masked_cts   <- all_ct_keys[!grepl("excl\\.", all_ct_keys)]
  masked_cts       <- all_ct_keys[ grepl("excl\\.", all_ct_keys)]
  extra_cts        <- setdiff(non_masked_cts, main_contrast)

  n_sig <- tryCatch({
    res <- as.data.frame(omni_results(obj))
    res[["padj"]] <- ifelse(is.na(res[["padj"]]), 1, res[["padj"]])
    sum(res[["padj"]] < 0.05)
  }, error = function(e) NA)

  # Helper to find a sPLS-DA key that matches a contrast
  .splsda_key_for <- function(ct) {
    if (!has_splsda) return(NULL)
    clean <- gsub("\\s*\\(excl\\..*\\)$", "", ct)
    clean <- trimws(clean)
    if (ct %in% names(obj$splsda)) return(ct)
    if (grepl(" vs ", clean, fixed = TRUE)) {
      parts <- trimws(strsplit(clean, " vs ", fixed = TRUE)[[1]])
      for (k in names(obj$splsda)) {
        if (all(sapply(parts, function(p)
          grepl(p, k, ignore.case = TRUE))))
          return(k)
      }
    }
    NULL
  }

  # Genotype labels for TF section (computed once for reuse)
  gt_keys_for_rmd <- tryCatch({
    ok <- setdiff(names(obj$tf_contrast_cache), main_contrast)
    ck <- ok[grepl("vs.Ctrl$|vs.Control$", ok, ignore.case = TRUE)]
    ck <- ck[!grepl("excl|mask", ck, ignore.case = TRUE)]
    trimws(gsub("\\s*vs\\s*(Ctrl|Control).*$", "", ck, ignore.case = TRUE))
  }, error = function(e) character(0))

  # ------------------------------------------------------------------
  # SECTION 2B: Build Rmd — Header + Overview + QC
  # ── FIX: All generated R code uses [["col"]] not $col
  # ------------------------------------------------------------------
  rmd <- c(
    "---",
    paste0('title: "', title, '"'),
    if (nzchar(author)) paste0('author: "', author, '"'),
    paste0('date: "', format(Sys.Date(), "%B %d, %Y"), '"'),
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    theme: flatly",
    "    self_contained: true",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE,",
    "                      fig.width=10, fig.height=7, dpi=150)",
    "library(OmniSeq)",
    "```",
    "",
    "# Dataset Overview",
    "",
    "```{r overview}",
    "print(obj)",
    "```",
    "",
    paste0("**Features:** ", nrow(obj$data),
           " | **Samples:** ", ncol(obj$data),
           " | **Data type:** ", obj$params$data_type),
    "",
    '```{r metadata-table}',
    'knitr::kable(omni_metadata(obj)[,',
    '  intersect(c("Condition","Genotype","Phenotype","Sex","Batch"),',
    '            colnames(omni_metadata(obj))), drop=FALSE],',
    '  caption="Sample metadata")',
    '```',
    "",
    "# Methods Snapshot",
    "",
    paste0("**Design formula:** `", deparse(obj$params$formula), "`  "),
    paste0("**Condition column:** `", obj$params$condition, "`  "),
    if (length(obj$params$covariates) > 0)
      paste0("**Covariates:** `",
             paste(obj$params$covariates, collapse = ", "), "`  ")
    else
      "**Covariates:** none  ",
    paste0("**Features (after filter):** ", nrow(obj$data),
           " | **Protein-coding filter:** ",
           if (!is.null(obj$params$filter_coding) &&
               isTRUE(obj$params$filter_coding))
             "applied" else "not applied"),
    paste0("**Normalisation:** DESeq2 VST (blind = FALSE, model-aware)  "),
    paste0("**DE engine:** DESeq2 (Wald test)  "),
    paste0("**LFC shrinkage:** ",
           obj$model_info$lfc_shrinkage %||% "none", "  "),
    "**Significance threshold:** padj < 0.05 (Benjamini-Hochberg)  ",
    "**sPLS-DA performance:** cross-validated BER reported; training-set AUC is not shown  ",
    "**TF activity:** ULM on contrast Wald t-statistics (contrast mode)  ",
    "",
    "# Known Limitations",
    "",
    paste0("> **Sample size:** n = ", ncol(obj$data),
           " total samples. Small cohort increases risk of overfitting ",
           "in supervised models (sPLS-DA) and reduces statistical ",
           "power for subgroup contrasts."),
    ">  ",
    "> **Group imbalance:** Some subgroups have very few samples ",
    "> (e.g. SPEN n=2, TLR7 n=4). BER for these contrasts should be ",
    "> interpreted with caution — majority-class prediction inflates ",
    "> overall accuracy. BER = 0.5 indicates no discrimination above chance.  ",
    ">  ",
    "> **TE locus ambiguity:** TE loci are multi-mapping by nature. ",
    "> Telescope-style unique assignment is used, but quantification ",
    "> uncertainty remains, especially for young/active families (L1HS, Alu).  ",
    ">  ",
    "> **Compositional effects:** Bulk RNA-seq measures a mixture of cell types. ",
    "> Changes in cell composition can confound both gene and TE DE results. ",
    "> xCell deconvolution partially addresses this but cannot fully correct for it.  ",
    ">  ",
    "> **Causal inference:** Mediation analysis assumes no unmeasured confounders ",
    "> between treatment, mediator, and outcome. With n < 15 per group this ",
    "> assumption cannot be fully tested.  ",
    ">  ",
    "> **TE-specific:** TE sensor scores are based on sum VST of matched loci — ",
    "> a proxy for pathway activation, not a direct functional measure.",
    "",
    "# Quality Control",
    "",
    '```{r qc-library-size, fig.height=6}',
    'tryCatch(plot_omni_qc_basics(obj, plot_type = "library_size"),',
    '  error = function(e) cat("QC metrics not available:", e$message,',
    '    "\\nRun run_omni_qc() and re-export to include QC plots."))',
    '```',
    "",
    '```{r qc-detected, fig.height=6}',
    'tryCatch(plot_omni_qc_basics(obj, plot_type = "detected"),',
    '  error = function(e) invisible(NULL))',
    '```',
    "",
    '```{r qc-distance, fig.height=7}',
    'tryCatch(plot_omni_qc_advanced(obj, plot_type = "distance"),',
    '  error = function(e) invisible(NULL))',
    '```',
    "",
    '```{r qc-pca-compare, fig.height=5}',
    'tryCatch(plot_omni_qc_advanced(obj, plot_type = "pca_compare"),',
    '  error = function(e) invisible(NULL))',
    '```',
    "",
    '```{r qc-mean-var, fig.height=5}',
    'tryCatch(plot_omni_qc_advanced(obj, plot_type = "mean_var"),',
    '  error = function(e) invisible(NULL))',
    '```'
  )

  # ------------------------------------------------------------------
  # SECTION 2C: Helper — build a full contrast section
  # ── FIX: All $col replaced with [["col"]] to avoid $$ Rmd bug
  # ── FIX: All plot/table chunks wrapped in tryCatch
  # ------------------------------------------------------------------
  .contrast_section <- function(ct, level = "##", id_suffix = NULL) {
    safe      <- id_suffix %||% gsub("[^a-zA-Z0-9]", "_", ct)
    is_masked <- grepl("excl\\.", ct)

    ct_res <- obj$results_cache[[ct]]
    n_sig_ct <- tryCatch({
      d <- as.data.frame(ct_res$results)
      d[["padj"]] <- ifelse(is.na(d[["padj"]]), 1, d[["padj"]])
      sum(d[["padj"]] < 0.05, na.rm = TRUE)
    }, error = function(e) "?")

    mask_note <- if (is_masked) {
      mask_label <- gsub(".*\\(excl\\.\\s*(.*)\\)$", "\\1", ct)
      paste0(" *(masked: ", mask_label, ")*")
    } else ""

    lines <- c(
      "",
      paste0(level, " ", ct, mask_note),
      "",
      paste0("**Significant features (padj < 0.05):** ", n_sig_ct),
      "",
      paste0('```{r ct-set-', safe, ', include=FALSE}'),
      paste0('obj$active_contrast <- "', ct, '"'),
      '```',
      "",
      paste0('```{r ct-volcano-', safe, '}'),
      'tryCatch(plot_omni_volcano(obj, lfc_thresh = 1, top_n = 15),',
      '  error = function(e) cat("Volcano not available:", e$message))',
      '```',
      "",
      paste0('```{r ct-heatmap-', safe, '}'),
      'tryCatch(plot_omni_heatmap(obj, ntop = 50),',
      '  error = function(e) cat("Heatmap not available:", e$message))',
      '```',
      "",
      paste0('```{r ct-topgenes-', safe, '}'),
      'tryCatch(plot_omni_top_genes(obj, n = 12),',
      '  error = function(e) cat("Top genes not available:", e$message))',
      '```',
      "",
      paste0('```{r ct-ma-', safe, '}'),
      'tryCatch(plot_omni_ma(obj),',
      '  error = function(e) cat("MA plot not available:", e$message))',
      '```'
    )

    # ── FIX: Table chunk uses [["padj"]] instead of $padj ──
    lines <- c(lines,
               "",
               paste0('```{r ct-table-', safe, '}'),
               'tryCatch({',
               '  res_ct <- as.data.frame(omni_results(obj))',
               '  res_ct[["padj"]] <- ifelse(is.na(res_ct[["padj"]]), 1, res_ct[["padj"]])',
               '  sig_ct <- res_ct[res_ct[["padj"]] < 0.05, ]',
               '  sig_ct <- sig_ct[order(sig_ct[["padj"]]), ]',
               '  if (nrow(sig_ct) > 0) {',
               '    show_cols <- intersect(c("log2FoldChange","pvalue","padj"), colnames(sig_ct))',
               '    knitr::kable(utils::head(sig_ct[, show_cols], 20),',
               paste0('      digits = 4, caption = "Top 20 features: ', ct, '")'),
               '  } else { cat("No significant features.") }',
               '}, error = function(e) cat("Table error:", e$message))',
               '```'
    )

    if (!is_masked) {
      lines <- c(lines,
                 "",
                 paste0('```{r ct-ora-', safe, '}'),
                 'tryCatch({',
                 '  .ora <- .omni_get_enrichment(obj)',
                 '  if (!is.null(.ora) && nrow(as.data.frame(.ora)) > 0) {',
                 '    plot_omni_ora(obj, plot_type = "dot")',
                 '    plot_omni_ora(obj, plot_type = "emap")',
                 '  } else { cat("No enriched GO terms for this contrast.") }',
                 '}, error = function(e) cat("ORA not available:", e$message))',
                 '```'
      )
    } else {
      lines <- c(lines,
                 "",
                 paste0('```{r ct-ora-masked-', safe, '}'),
                 'tryCatch({',
                 '  .ora <- .omni_get_enrichment(obj)',
                 '  if (!is.null(.ora) && nrow(as.data.frame(.ora)) > 0) {',
                 '    plot_omni_ora(obj, plot_type = "dot")',
                 '    plot_omni_ora(obj, plot_type = "emap")',
                 '    plot_omni_compare_clusters(obj)',
                 '  } else { cat("No enriched GO terms after masking.") }',
                 '}, error = function(e) cat("ORA not available:", e$message))',
                 '```',
                 "",
                 paste0('```{r ct-gsea-masked-', safe, '}'),
                 'tryCatch({',
                 '  .gsea <- .omni_get_gsea(obj)',
                 '  if (!is.null(.gsea) && nrow(as.data.frame(.gsea)) > 0) {',
                 '    plot_omni_gsea(obj, show_n = 10)',
                 '  } else { cat("No enriched GSEA terms after masking.") }',
                 '}, error = function(e) cat("GSEA not available:", e$message))',
                 '```'
      )
    }

    splsda_k <- .splsda_key_for(ct)
    if (!is.null(splsda_k)) {
      lines <- c(lines,
                 "",
                 paste0('```{r ct-splsda-', safe, '}'),
                 'tryCatch({',
                 paste0('  plot_omni_splsda_scores(obj, key = "', splsda_k, '")'),
                 paste0('  plot_omni_splsda_vip(obj, key = "', splsda_k, '")'),
                 paste0('  tryCatch(plot_omni_splsda_roc(obj, key = "',
                        splsda_k, '", comp = 1),',
                        ' error = function(e) invisible(NULL))'),
                 '}, error = function(e) cat("sPLS-DA not available:", e$message))',
                 '```'
      )
    }

    # Per-contrast TF activity — wrapped in tryCatch
    if (!is_masked && has_tf) {
      lines <- c(lines,
                 "",
                 paste0('```{r ct-tf-', safe, '}'),
                 'tryCatch({',
                 paste0('  obj$active_contrast <- "', ct, '"'),
                 '  plot_omni_tf_activity(obj, n = 15)',
                 paste0('  obj$active_contrast <- "', main_contrast, '"'),
                 '}, error = function(e) {',
                 paste0('  obj$active_contrast <- "', main_contrast, '"'),
                 '  cat("TF activity not available for this contrast:", e$message)',
                 '})',
                 '```'
      )
    }

    lines <- c(lines,
               "",
               paste0('```{r ct-restore-', safe, ', include=FALSE}'),
               paste0('obj$active_contrast <- "', main_contrast, '"'),
               '```'
    )

    lines
  }

  # ------------------------------------------------------------------
  # SECTION 2D: Main DE section
  # ------------------------------------------------------------------
  if (has_de) {
    rmd <- c(rmd,
             "",
             "# Differential Expression",
             "",
             paste0("**Main contrast:** `", main_contrast, "`  "),
             paste0("**Significant features (padj < 0.05):** ", n_sig),
             "",
             '```{r de-pca}',
             'plot_omni_pca(obj)',
             '```',
             "",
             '```{r de-pval}',
             'plot_omni_pvalue_hist(obj)',
             '```',
             .contrast_section(main_contrast, level = "##",
                               id_suffix = "main")
    )
  }

  # ------------------------------------------------------------------
  # SECTION 2E: Additional contrasts (non-masked)
  # ------------------------------------------------------------------
  if (length(extra_cts) > 0) {
    rmd <- c(rmd, "",
             "# Additional Contrasts",
             "",
             paste0("The following ", length(extra_cts),
                    " additional contrast(s) were run."))
    for (ct in extra_cts)
      rmd <- c(rmd, .contrast_section(ct, level = "##"))
  }

  # ------------------------------------------------------------------
  # SECTION 2F: Enrichment (main contrast, global)
  # ------------------------------------------------------------------
  if (has_ora) {
    rmd <- c(rmd, "",
             "# GO Enrichment (ORA — Main Contrast)",
             "",
             paste0("**Enriched GO terms:** ",
                    nrow(as.data.frame(obj$enrichment))),
             "",
             '```{r ora-dot, fig.height=7}',
             'plot_omni_ora(obj, plot_type = "dot")',
             '```',
             "",
             '```{r ora-cnet, fig.height=8}',
             'plot_omni_ora(obj, plot_type = "cnet")',
             '```',
             "",
             '```{r ora-emap, fig.height=7}',
             'plot_omni_ora(obj, plot_type = "emap")',
             '```',
             "",
             '```{r ora-compare, fig.height=7}',
             'plot_omni_compare_clusters(obj)',
             '```',
             "",
             '```{r ora-table}',
             'ora_df <- as.data.frame(obj[["enrichment"]])',
             'knitr::kable(ora_df[, intersect(c("Description","p.adjust","Count","geneID"), colnames(ora_df))],',
             '  digits = 4, caption = "Enriched GO terms")',
             '```'
    )
  }

  if (has_gsea) {
    rmd <- c(rmd, "",
             "# GSEA",
             "",
             '```{r gsea}',
             'plot_omni_gsea(obj, show_n = 10)',
             '```',
             "",
             '```{r leading-edge}',
             'tryCatch(plot_omni_leading_edge(obj, n_pathways = 5),',
             '  error = function(e) cat("Leading edge not available:", e$message))',
             '```'
    )
  }

  if (has_ssgsea) {
    rmd <- c(rmd, "",
             "# ssGSEA Pathway Scores",
             "",
             '```{r ssgsea}',
             'plot_omni_ssgsea_heatmap(obj)',
             '```'
    )
  }

  # ------------------------------------------------------------------
  # SECTION 2G: TF activity — redesigned with independent sections
  # ── FIX: Each panel (A-G) on its own Rmd subsection
  # ── FIX: F panels (individual genotype profiles) added
  # ── FIX: G panels (genotype-specific TFs) added
  # ------------------------------------------------------------------
  if (has_tf) {
    rmd <- c(rmd, "",
             "# Transcription Factor Activity {.tabset}",
             "",
             paste0("> **Method:** ULM (Univariate Linear Model) on ",
                    "DESeq2 Wald t-statistics"),
             paste0("> using the CollecTRI database. ",
                    "A positive score means the TF's target genes are ",
                    "systematically upregulated. A negative score means repression."),
             "> TFs are cross-referenced against MSigDB Hallmark IFN gene sets.",
             "",
             "## Top activated and repressed TFs (Main Contrast)",
             "",
             '```{r tf-activity, fig.width=10, fig.height=10}',
             'plot_omni_tf_activity(obj, n = 15)',
             '```',
             "",
             "## TF regulatory network (Main Contrast)",
             "",
             '```{r tf-network, fig.width=12, fig.height=10}',
             'plot_omni_tf_network(obj, plot_type = "static")',
             '```'
    )

    # Per-contrast TF activity subsections
    if (!is.null(obj$tf_contrast_cache)) {
      tf_ct_other <- setdiff(names(obj$tf_contrast_cache), main_contrast)
      for (tf_ct in tf_ct_other) {
        safe_tfc <- gsub("[^a-zA-Z0-9]", "_", tf_ct)
        rmd <- c(rmd, "",
                 paste0("## TF Activity: ", tf_ct), "",
                 paste0('```{r tf-act-', safe_tfc,
                        ', fig.width=10, fig.height=10}'),
                 'tryCatch({',
                 paste0('  obj$active_contrast <- "', tf_ct, '"'),
                 '  plot_omni_tf_activity(obj, n = 15)',
                 paste0('  obj$active_contrast <- "', main_contrast, '"'),
                 '}, error = function(e) {',
                 paste0('  obj$active_contrast <- "', main_contrast, '"'),
                 '  cat("TF activity not available:", e$message)',
                 '})',
                 '```'
        )
      }
    }

    # Genotype analysis — compute once, render each panel separately
    rmd <- c(rmd, "",
             "## Genotype: Pattern Distribution", "",
             '```{r tf-geno-compute, fig.width=10, fig.height=6}',
             'tf_plots <- tryCatch(plot_omni_tf_genotype(obj),',
             '  error = function(e) { cat("Genotype TF analysis:", e$message); NULL })',
             'if (!is.null(tf_plots) && "A_pattern_distribution" %in% names(tf_plots))',
             '  print(tf_plots[["A_pattern_distribution"]])',
             '```', "",
             "## Genotype: Top Shared Repressed", "",
             '```{r tf-geno-B, fig.width=10, fig.height=10}',
             'if (exists("tf_plots") && !is.null(tf_plots) &&',
             '    "B_top_shared_repressed" %in% names(tf_plots))',
             '  print(tf_plots[["B_top_shared_repressed"]])',
             '```'
    )

    # Dynamic C/D/F/G panels based on genotype labels
    if (length(gt_keys_for_rmd) >= 2) {
      rpairs <- utils::combn(gt_keys_for_rmd, 2, simplify = FALSE)

      for (rp in rpairs) {
        g1 <- rp[1]; g2 <- rp[2]
        sp <- paste0(gsub("[^a-zA-Z0-9]", "_", g1), "_",
                     gsub("[^a-zA-Z0-9]", "_", g2))
        ck <- paste0("C_shared_activated_", g1, "_and_", g2)
        dk <- paste0("D_discordant_", g1, "_vs_", g2)

        rmd <- c(rmd, "",
                 paste0("## Shared Activated: ", g1, " & ", g2), "",
                 paste0('```{r tf-C-', sp, ', fig.width=10, fig.height=10}'),
                 paste0('if (exists("tf_plots") && !is.null(tf_plots) && "',
                        ck, '" %in% names(tf_plots))'),
                 paste0('  print(tf_plots[["', ck, '"]])'),
                 '```', "",
                 paste0("## Discordant TFs: ", g1, " vs ", g2), "",
                 paste0('```{r tf-D-', sp, ', fig.width=10, fig.height=10}'),
                 paste0('if (exists("tf_plots") && !is.null(tf_plots) && "',
                        dk, '" %in% names(tf_plots))'),
                 paste0('  print(tf_plots[["', dk, '"]])'),
                 '```'
        )
      }

      # D combined
      rmd <- c(rmd, "",
               "## Discordant TFs: All Pairs Combined", "",
               '```{r tf-D-combined, fig.width=24, fig.height=10}',
               'if (exists("tf_plots") && !is.null(tf_plots) &&',
               '    "D_discordant_combined" %in% names(tf_plots))',
               '  print(tf_plots[["D_discordant_combined"]])',
               '```'
      )

      # F panels — individual genotype profiles
      for (gt in gt_keys_for_rmd) {
        sgt <- gsub("[^a-zA-Z0-9]", "_", gt)
        fk  <- paste0("F_profile_", gt)
        rmd <- c(rmd, "",
                 paste0("## TF Profile: ", gt, " vs Controls"), "",
                 paste0('```{r tf-F-', sgt, ', fig.width=10, fig.height=10}'),
                 paste0('if (exists("tf_plots") && !is.null(tf_plots) && "',
                        fk, '" %in% names(tf_plots))'),
                 paste0('  print(tf_plots[["', fk, '"]])'),
                 '```'
        )
      }

      # G panels — genotype-specific TFs
      for (gt in gt_keys_for_rmd) {
        sgt <- gsub("[^a-zA-Z0-9]", "_", gt)
        gk  <- paste0("G_specific_", gt)
        rmd <- c(rmd, "",
                 paste0("## Genotype-Specific TFs: ", gt), "",
                 paste0('```{r tf-G-', sgt, ', fig.width=10, fig.height=10}'),
                 paste0('if (exists("tf_plots") && !is.null(tf_plots) && "',
                        gk, '" %in% names(tf_plots))'),
                 paste0('  print(tf_plots[["', gk, '"]])'),
                 '```'
        )
      }
    }

    # E — heatmap (always last in TF section)
    rmd <- c(rmd, "",
             "## Per-Genotype Score Heatmap", "",
             '```{r tf-geno-E, fig.width=12, fig.height=14}',
             'if (exists("tf_plots") && !is.null(tf_plots) &&',
             '    "E_per_genotype_heatmap" %in% names(tf_plots))',
             '  print(tf_plots[["E_per_genotype_heatmap"]])',
             '```'
    )
  }

  # ------------------------------------------------------------------
  # SECTION 2H: Deconvolution + STRING
  # ------------------------------------------------------------------
  if (has_deconv) {
    rmd <- c(rmd, "",
             "# Cell Type Deconvolution",
             "",
             '```{r deconv}',
             'plot_omni_deconv_heatmap(obj)',
             '```'
    )
  }

  if (has_string) {
    n_nodes <- tryCatch(
      length(unique(c(obj$string_network$interactions$from_gene,
                      obj$string_network$interactions$to_gene))),
      error = function(e) "?"
    )
    n_edges <- nrow(obj$string_network$interactions)
    rmd <- c(rmd, "",
             "# STRING Protein Interaction Network",
             "",
             paste0("Network built from top ",
                    obj$string_network$n_genes %||% "?",
                    " DEGs | **", n_nodes, " proteins** | **",
                    n_edges, " interactions** | Score threshold: ",
                    obj$string_network$score_thresh %||% "?"),
             "",
             '```{r string-network, fig.width=10, fig.height=8}',
             'plot_omni_string_network(obj, plot_type = "static")',
             '```'
    )
  }

  # SC reference deconvolution
  has_sc_rmd <- !is.null(obj$sc_deconv) && !is.null(obj$sc_deconv$scores)
  has_ds_rmd <- has_sc_rmd && !is.null(obj$sc_deconv$disease_states)

  if (has_sc_rmd) {
    rmd <- c(rmd, "",
             "# scRNA-seq Reference Deconvolution {.tabset}", "",
             paste0("> **Reference:** ", obj$sc_deconv$reference %||% "scRNA-seq",
                    " | **Method:** ", obj$sc_deconv$method %||% "signature",
                    " | **Cell types:** ", length(obj$sc_deconv$cell_types),
                    " | **Gene overlap:** ", obj$sc_deconv$gene_overlap %||% "?"),
             "",
             "## Layer 1: Cell Composition (Healthy Reference)", "",
             '```{r sc-heat-comp, fig.width=12, fig.height=8}',
             'plot_omni_sc_heatmap(obj, layer = "composition")',
             '```', "",
             "## Composition by Genotype", "",
             '```{r sc-geno-comp, fig.width=14, fig.height=10}',
             'tryCatch(plot_omni_sc_genotype(obj, layer = "composition"),',
             '  error = function(e) cat("Plot failed:", e$message))',
             '```', "",
             "## Composition vs IRG-S", "",
             '```{r sc-corr-comp, fig.width=10, fig.height=8}',
             'tryCatch(plot_omni_sc_correlation(obj, feature = "IRGS", layer = "composition"),',
             '  error = function(e) cat("Plot failed:", e$message))',
             '```', "",
             "## Composition Profile", "",
             '```{r sc-prof-comp, fig.width=10, fig.height=7}',
             'tryCatch(plot_omni_sc_profile(obj, layer = "composition"),',
             '  error = function(e) cat("Plot failed:", e$message))',
             '```'
    )

    if (has_ds_rmd) {
      rmd <- c(rmd, "",
               "## Layer 2: Disease Activation (Disease Reference)", "",
               paste0("> Disease states derived from **",
                      obj$sc_deconv$disease_states$disease_group %||% "untreated",
                      "** vs **",
                      obj$sc_deconv$disease_states$healthy_group %||% "control",
                      "** within each cell type."),
               "",
               '```{r sc-heat-act, fig.width=12, fig.height=8}',
               'plot_omni_sc_heatmap(obj, layer = "activation")',
               '```', "",
               "## Activation by Genotype", "",
               '```{r sc-geno-act, fig.width=14, fig.height=10}',
               'tryCatch(plot_omni_sc_genotype(obj, layer = "activation"),',
               '  error = function(e) cat("Plot failed:", e$message))',
               '```', "",
               "## Activation Profile", "",
               '```{r sc-prof-act, fig.width=10, fig.height=7}',
               'tryCatch(plot_omni_sc_profile(obj, layer = "activation"),',
               '  error = function(e) cat("Plot failed:", e$message))',
               '```', "",
               "## Two-Layer Combined", "",
               '```{r sc-two-layer, fig.width=14, fig.height=10}',
               'tryCatch(plot_omni_sc_two_layer(obj),',
               '  error = function(e) cat("Plot failed:", e$message))',
               '```'
      )
    }

    if (!is.null(obj$deconv)) {
      rmd <- c(rmd, "",
               "## Benchmark: xCell vs Signature Scoring", "",
               '```{r sc-benchmark, fig.width=12, fig.height=8}',
               'tryCatch(plot_omni_sc_benchmark(obj),',
               '  error = function(e) cat("Benchmark failed:", e$message))',
               '```')
    }
  }

  # ------------------------------------------------------------------
  # SECTION 2I: sPLS-DA — all models
  # ------------------------------------------------------------------
  if (has_splsda) {
    rmd <- c(rmd, "",
             "# sPLS-DA",
             "",
             paste0("**Models stored:** ",
                    paste(names(obj$splsda), collapse = ", "))
    )
    for (sk in names(obj$splsda)) {
      safe_sk <- gsub("[^a-zA-Z0-9]", "_", sk)
      rmd <- c(rmd,
               "",
               paste0("## Model: `", sk, "`"),
               "",
               paste0('```{r splsda-scores-', safe_sk, '}'),
               paste0('plot_omni_splsda_scores(obj, key = "', sk, '")'),
               '```',
               "",
               paste0('```{r splsda-vip-', safe_sk, '}'),
               paste0('plot_omni_splsda_vip(obj, key = "', sk, '")'),
               '```',
               "",
               paste0('```{r splsda-roc-', safe_sk, '}'),
               paste0('tryCatch(plot_omni_splsda_roc(obj, key = "',
                      sk, '", comp = 1),',
                      ' error = function(e) cat("ROC not available:",',
                      ' e$message))'),
               '```',
               "",
               paste0('```{r splsda-sig-', safe_sk, '}'),
               paste0('sig_tab <- tryCatch(',
                      'omni_splsda_signature(obj, key = "',
                      sk, '", comp = 1),',
                      ' error = function(e) NULL)'),
               'if (!is.null(sig_tab)) {',
               '  knitr::kable(utils::head(sig_tab, 20), digits = 4,',
               paste0('    caption = "Top sPLS-DA features: ', sk, '")'),
               '} else { cat("Signature not available.") }',
               '```'
      )
    }
  }

  # ------------------------------------------------------------------
  # SECTION 2J: WGCNA — all keys
  # ------------------------------------------------------------------
  if (has_wgcna) {
    rmd <- c(rmd, "",
             "# WGCNA Co-expression Network",
             "",
             paste0("**Networks stored:** ",
                    paste(names(obj$wgcna), collapse = ", "))
    )
    for (wk in names(obj$wgcna)) {
      safe_wk <- gsub("[^a-zA-Z0-9]", "_", wk)
      rmd <- c(rmd,
               "",
               paste0("## Network: `", wk, "`"),
               "",
               paste0('```{r wgcna-dendro-', safe_wk, ', fig.height=6}'),
               paste0('plot_omni_wgcna_dendrogram(obj, key = "', wk, '")'),
               '```',
               "",
               paste0('```{r wgcna-traits-', safe_wk, ', fig.height=8}'),
               paste0('plot_omni_wgcna_trait_heatmap(obj, key = "', wk, '")'),
               '```',
               "",
               paste0('```{r wgcna-hub-', safe_wk, ', fig.height=8}'),
               paste0('plot_omni_wgcna_hub_network(obj, key = "',
                      wk, '", plot_type = "static")'),
               '```'
      )
    }
  }

  # ------------------------------------------------------------------
  # SECTION 2K: TE Analysis
  # ------------------------------------------------------------------
  if (has_te) {
    rmd <- c(rmd, "",
             "# Transposable Element Analysis",
             "",
             '```{r te-class-pie}',
             'plot_omni_te_class_pie(obj)',
             '```',
             "",
             '```{r te-family}',
             'plot_omni_te_family(obj)',
             '```',
             "",
             '```{r te-lfc-density, fig.height=6}',
             'plot_omni_te_lfc_density(obj, min_loci = 10)',
             '```',
             "",
             '```{r te-volcano}',
             'plot_omni_te_volcano(obj)',
             '```',
             "",
             '```{r te-heatmap}',
             'plot_omni_te_heatmap(obj)',
             '```'
    )
  }

  if (!is.null(obj$te_great)) {
    rmd <- c(rmd, "",
             "## GREAT Genomic Enrichment",
             "",
             '```{r te-great}',
             'plot_omni_te_great(obj, show_n = 15)',
             '```'
    )
  }

  # ── FIX: TE-gene correlation uses [["padj"]] not $padj ──
  if (has_te_corr) {
    rmd <- c(rmd, "",
             "## TE-Gene Correlation",
             "",
             "### IRG-S Score",
             "",
             '```{r te-gene-corr-irgs}',
             'plot_omni_te_gene_corr(gene_obj = gene_obj, te_obj = obj,',
             '  marker_genes = OmniSeq::omniseq_isgs_28,',
             '  te_class = "ALL", score_name = "IRG-S Score")',
             '```',
             "",
             "### ISG-independent signature",
             "",
             '```{r te-gene-corr-nonisg}',
             'g_res <- as.data.frame(omni_results(gene_obj))',
             'g_res[["padj"]] <- ifelse(is.na(g_res[["padj"]]), 1, g_res[["padj"]])',
             'non_isg <- setdiff(rownames(g_res)[g_res[["padj"]] < 0.05], OmniSeq::omniseq_isgs_28)',
             'if (length(non_isg) >= 3) {',
             '  plot_omni_te_gene_corr(gene_obj = gene_obj, te_obj = obj,',
             '    marker_genes = non_isg, te_class = "ALL",',
             '    score_name = "ISG-independent score")',
             '} else {',
             '  cat("Fewer than 3 ISG-independent DEGs available.")',
             '}',
             '```'
    )
  }

  # ------------------------------------------------------------------
  # SECTION 2L: Sensors + DIABLO + Mediation
  # ------------------------------------------------------------------
  if (has_sensors) {
    so <- sensor_obj_name
    rmd <- c(rmd, "",
             "# Innate Immune Sensor Scores",
             "",
             "> Sensor activation scores computed from TE VST expression.",
             "> cGAS-STING (cytoplasmic DNA), RIG-I/MDA5 (dsRNA), TLR7/TLR8 (ssRNA).",
             "> Scores are z-scored across samples for comparability.",
             "",
             paste0('```{r sensors, fig.height=6}'),
             paste0('tryCatch(plot_omni_sensor_scores(', so, '),'),
             '  error = function(e) cat("Sensor score plot failed:", e$message))',
             '```',
             "",
             paste0('```{r sensors-corr, fig.height=6}'),
             paste0('tryCatch(plot_omni_sensor_irgs_corr(', so, '),'),
             '  error = function(e) cat("Sensor-IRGs correlation plot failed:", e$message))',
             '```',
             "",
             paste0('```{r sensors-heatmap, fig.height=7}'),
             paste0('tryCatch(plot_omni_sensor_heatmap(', so, '),'),
             '  error = function(e) cat("Sensor heatmap failed:", e$message))',
             '```'
    )
  }

  # ── FIX: DIABLO uses [["..."]] not $... to avoid $$ bug ──
  if (has_diablo) {
    first_dk   <- names(obj$diablo)[1]
    first_dobj <- obj$diablo[[first_dk]]
    rmd <- c(rmd, "",
             "# DIABLO Multi-block Integration",
             "",
             paste0("**Model:** `", first_dk, "`  "),
             paste0("**Components:** ", first_dobj$ncomp),
             "",
             "## Sample plot",
             "",
             '```{r diablo-indiv}',
             paste0('model <- obj[["diablo"]][["', first_dk, '"]][["model"]]'),
             'mixOmics::plotIndiv(model, legend = TRUE, title = "DIABLO - Sample Plot")',
             '```',
             "",
             "## Block overview",
             "",
             '```{r diablo-overview, fig.width=8, fig.height=8}',
             paste0('model <- obj[["diablo"]][["', first_dk, '"]][["model"]]'),
             'mixOmics::plotDiablo(model, ncomp = 1)',
             'title("DIABLO - Block Correlation Overview (Component 1)",',
             '      outer = TRUE, line = -1.5, cex.main = 1.2)',
             '```',
             "",
             "## Circos correlation plot",
             "",
             '```{r diablo-circos, fig.height=9, fig.width=9}',
             paste0('model <- obj[["diablo"]][["', first_dk, '"]][["model"]]'),
             'tryCatch(mixOmics::circosPlot(model, cutoff = 0.7,',
             '  line = TRUE, size.labels = 1.5),',
             '  error = function(e) cat("Circos plot failed:", e$message))',
             '```',
             "",
             "## Loadings (Component 1)",
             "",
             '```{r diablo-loadings}',
             paste0('model <- obj[["diablo"]][["', first_dk, '"]][["model"]]'),
             'mixOmics::plotLoadings(model, comp = 1, contrib = "max", method = "median")',
             '```',
             "",
             "## Gene loading table",
             "",
             '```{r diablo-gene-table}',
             paste0('model <- obj[["diablo"]][["', first_dk, '"]][["model"]]'),
             'gene_load <- as.data.frame(model[["loadings"]][["genes"]])',
             'gene_load <- gene_load[rowSums(gene_load != 0) > 0, , drop = FALSE]',
             'gene_load <- gene_load[order(-abs(gene_load[, 1])), ]',
             'knitr::kable(utils::head(gene_load, 20), digits = 4,',
             '  caption = "Top gene loadings (Comp 1)")',
             '```',
             "",
             "## TE family loading table",
             "",
             '```{r diablo-te-table}',
             paste0('model <- obj[["diablo"]][["', first_dk, '"]][["model"]]'),
             'te_load <- as.data.frame(model[["loadings"]][["TEs"]])',
             'te_load <- te_load[rowSums(te_load != 0) > 0, , drop = FALSE]',
             'knitr::kable(te_load, digits = 4, caption = "TE family loadings")',
             '```'
    )
  }

  # ── FIX: Mediation uses [["..."]] not $... to avoid $$ bug ──
  if (has_med) {
    rmd <- c(rmd, "",
             "# Mediation Analysis",
             "",
             paste0("**Mediators tested:** ",
                    paste(names(obj$mediation), collapse = ", ")),
             "",
             '```{r mediation-forest}',
             'plot_omni_mediation(obj)',
             '```',
             "",
             '```{r mediation-table}',
             'med_rows <- lapply(names(obj[["mediation"]]), function(nm) {',
             '  s <- obj[["mediation"]][[nm]][["summary"]]',
             '  if (is.null(s)) return(NULL)',
             '  s[["Mediator"]] <- nm; s',
             '})',
             'med_all <- do.call(rbind, Filter(Negate(is.null), med_rows))',
             'knitr::kable(med_all, digits = 3,',
             '  caption = "Mediation analysis summary")',
             '```'
    )
  }

  # ------------------------------------------------------------------
  # SECTION 2M: Masked analysis sections
  # ------------------------------------------------------------------
  if (length(masked_cts) > 0) {
    rmd <- c(rmd, "",
             "# Masked Analyses",
             "",
             "Each section below re-runs DE visualizations and enrichment after",
             "removing a specific gene signature. Plots and titles include the",
             "mask label to distinguish from the unmasked results."
    )
    for (mk in masked_cts) {
      rmd <- c(rmd, .contrast_section(
        mk, level = "##",
        id_suffix = gsub("[^a-zA-Z0-9]", "_", mk)))
    }
  }

  # ------------------------------------------------------------------
  # SECTION 2N: Session info + Render
  # ------------------------------------------------------------------
  rmd <- c(rmd, "",
           "# Session Info",
           "",
           '```{r session}',
           'sessionInfo()',
           '```'
  )

  # Render
  rmd_path <- tempfile(fileext = ".Rmd")
  writeLines(rmd, rmd_path)

  # Save a permanent copy of the Rmd in outdir for debugging
  rmd_perm <- file.path(outdir, sub("\\.html$", ".Rmd", report_file))
  file.copy(rmd_path, rmd_perm, overwrite = TRUE)

  render_env          <- new.env(parent = globalenv())
  render_env$obj      <- obj
  render_env$gene_obj <- gene_obj

  tryCatch({
    rmarkdown::render(
      input       = rmd_path,
      output_file = html_path,
      envir       = render_env,
      quiet       = FALSE
    )
    message("[export_omniseq_report] Report saved: ", html_path)
    if (open && interactive()) utils::browseURL(html_path)
  }, error = function(e) {
    message("[export_omniseq_report] Rendering FAILED.")
    message("  Error   : ", e$message)
    message("  Rmd kept: ", rmd_perm)
    message("  To debug: rmarkdown::render('", rmd_perm, "')")
    message("  Or open the Rmd in RStudio and knit manually.")
  })

  return(invisible(html_path))
}
