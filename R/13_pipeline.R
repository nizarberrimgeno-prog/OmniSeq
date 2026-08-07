# =============================================================================
# OmniSeq — 13_Pipeline.R
#
# High-level pipeline wrappers that run multiple analysis steps in sequence.
#
# SECTION INDEX (for targeted replacement):
#   SECTION 1:  run_omniseq_auto()
#   SECTION 2:  run_omniseq_guided()
#     2A: Init + data type detection
#     2B: QC
#     2C: Differential Expression + Custom Contrasts
#     2D: GO Enrichment (ORA) — with preservation fix
#     2E: GSEA + ssGSEA
#     2F: TF Enrichment
#     2G: Deconvolution
#     2H: STRING
#     2I: sPLS-DA
#     2J: WGCNA
#     2K: TE Sensor Scoring
#     2L: Export
#     2M: Masked Analysis Prompt
#     2N: HTML Report + Final Summary
#   SECTION 3:  run_omniseq_masked()
#   SECTION 4:  .run_omniseq_guided_te()
#
# FIX applied: ORA enrichment preserved across per-contrast loop (Section 2D)
# =============================================================================


# =============================================================================
# SECTION 1: AUTOMATIC PIPELINE
# =============================================================================

#' Run the Full OmniSeq Pipeline Automatically
#'
#' Runs the complete OmniSeq analysis pipeline from raw counts to exported
#' results using sensible defaults. All steps are protected with
#' `tryCatch` — a failure in one step is reported but does not abort
#' subsequent steps.
#'
#' @param count_matrix Matrix or data frame of raw counts (features x samples).
#' @param metadata Data frame with sample metadata (rows = samples).
#' @param condition_col Character. Metadata column for main condition.
#' @param covariates Character vector. Covariate columns for the design.
#' @param ref_level Character. Reference level for the condition factor.
#' @param species Character. Species name. Default: `"Homo sapiens"`.
#' @param org_db Character. Annotation database. Default: `"org.Hs.eg.db"`.
#' @param filter_coding Logical. Filter to protein-coding genes. Default: `TRUE`.
#' @param te_obj Optional OmniSeq TE object for sensor scoring and DIABLO.
#' @param outdir Character. Export directory. Default: `"./omniseq_output"`.
#' @param export Logical. Export results on completion. Default: `TRUE`.
#' @param n_top_genes Integer. Genes for WGCNA. Default: `3000L`.
#' @param gsea_category Character. MSigDB category for GSEA. Default: `"H"`.
#'
#' @return A fully populated OmniSeq object.
#' @export
run_omniseq_auto <- function(count_matrix,
                             metadata,
                             condition_col,
                             covariates    = NULL,
                             ref_level     = NULL,
                             species       = "Homo sapiens",
                             org_db        = "org.Hs.eg.db",
                             filter_coding = TRUE,
                             te_obj        = NULL,
                             outdir        = "./omniseq_output",
                             export        = TRUE,
                             n_top_genes   = 3000L,
                             gsea_category = "H") {

  status <- list()
  .step <- function(name, expr) {
    message("\n", strrep("=", 55))
    message("  STEP: ", name)
    message(strrep("=", 55))
    tryCatch({
      result <- force(expr)
      status[[name]] <<- "PASS"
      result
    }, error = function(e) {
      status[[name]] <<- paste("FAIL —", e$message)
      message("[run_omniseq_auto] SKIP: ", name, " failed — ", e$message)
      NULL
    })
  }

  # ------------------------------------------------------------------
  # 1. INITIALISE
  # ------------------------------------------------------------------
  obj <- .step("init_omniseq", {
    init_omniseq(
      count_matrix  = count_matrix,
      metadata      = metadata,
      condition_col = condition_col,
      covariates    = covariates,
      filter_coding = filter_coding,
      species       = species,
      org_db        = org_db,
      ref_level     = ref_level
    )
  })
  if (is.null(obj)) {
    message("[run_omniseq_auto] Cannot continue — initialisation failed.")
    return(invisible(NULL))
  }

  # ------------------------------------------------------------------
  # 2. QC
  # ------------------------------------------------------------------
  obj <- .step("run_omni_qc", {
    o <- run_omni_qc(obj)
    plot_omni_qc_basics(o, plot_type = "library_size")
    plot_omni_qc_basics(o, plot_type = "detected")
    plot_omni_qc_advanced(o, plot_type = "distance")
    plot_omni_qc_advanced(o, plot_type = "pca_compare")
    plot_omni_qc_advanced(o, plot_type = "mean_var")
    o
  }) %||% obj

  # ------------------------------------------------------------------
  # 3. DIFFERENTIAL EXPRESSION
  # ------------------------------------------------------------------
  obj <- .step("run_omni_de", {
    run_omni_de(obj)
  }) %||% obj

  # ------------------------------------------------------------------
  # 4. VISUALIZATIONS
  # ------------------------------------------------------------------
  .step("visualizations", {
    plot_omni_pca(obj)
    plot_omni_volcano(obj, lfc_thresh = 1, top_n = 15)
    plot_omni_ma(obj)
    plot_omni_heatmap(obj, ntop = 50)
    plot_omni_top_genes(obj, n = 12)
    plot_omni_pvalue_hist(obj)
    TRUE
  })

  # ------------------------------------------------------------------
  # 5. ENRICHMENT
  # ------------------------------------------------------------------
  obj <- .step("run_omni_enrichment", {
    run_omni_enrichment(obj)
  }) %||% obj

  obj <- .step("run_omni_gsea", {
    run_omni_gsea(obj, category = gsea_category)
  }) %||% obj

  obj <- .step("run_omni_ssgsea", {
    run_omni_ssgsea(obj, category = gsea_category)
  }) %||% obj

  obj <- .step("run_omni_pathway_de", {
    meta     <- omni_metadata(obj)
    grps     <- levels(meta[[condition_col]])
    grp1     <- grps[grps != ref_level][1]
    grp2     <- ref_level %||% grps[1]
    run_omni_pathway_de(obj, group1 = grp1, group2 = grp2)
  }) %||% obj

  # ------------------------------------------------------------------
  # 6. TF ENRICHMENT (main contrast + all per-genotype contrasts)
  # ------------------------------------------------------------------
  obj <- .step("run_omni_tf", {
    o <- run_omni_tf(obj, database = "collectri", method = "ulm")
    main_ct   <- o$active_contrast
    other_cts <- setdiff(
      omni_contrasts(o),
      grep("excl|mask|isgs|schoggins", omni_contrasts(o),
           ignore.case = TRUE, value = TRUE)
    )
    other_cts <- setdiff(other_cts, main_ct)
    for (ct in other_cts) {
      o$active_contrast <- ct
      o <- tryCatch(
        run_omni_tf(o, database = "collectri", method = "ulm"),
        error = function(e) {
          message("[run_omni_tf] Skipping contrast '", ct, "': ", e$message)
          o
        }
      )
    }
    o$active_contrast <- main_ct
    o
  }) %||% obj

  # ------------------------------------------------------------------
  # 7. DECONVOLUTION
  # ------------------------------------------------------------------
  obj <- .step("run_omni_deconvolution", {
    run_omni_deconvolution(obj, method = "immunedeconv",
                           deconv_method = "xcell")
  }) %||% obj

  # ------------------------------------------------------------------
  # 8. STRING NETWORK
  # ------------------------------------------------------------------
  obj <- .step("run_omni_string", {
    run_omni_string(obj, n = 50, string_score = 400)
  }) %||% obj

  # ------------------------------------------------------------------
  # 9. sPLS-DA
  # ------------------------------------------------------------------
  obj <- .step("run_omni_splsda", {
    run_omni_splsda(obj, group_col = condition_col)
  }) %||% obj

  # ------------------------------------------------------------------
  # 10. WGCNA
  # ------------------------------------------------------------------
  obj <- .step("run_omni_wgcna", {
    run_omni_wgcna(obj,
                   use_masked  = FALSE,
                   n_top_genes = n_top_genes,
                   key         = "full")
  }) %||% obj

  # ------------------------------------------------------------------
  # 11. SENSOR SCORING (if TE object provided)
  # ------------------------------------------------------------------
  if (!is.null(te_obj)) {
    obj <- .step("run_omni_sensor_score", {
      run_omni_sensor_score(obj,
                            te_obj    = te_obj,
                            group_col = condition_col)
    }) %||% obj

    obj <- .step("run_omni_diablo", {
      run_omni_diablo(obj,
                      te_obj    = te_obj,
                      group_col = condition_col,
                      ncomp     = 2L)
    }) %||% obj
  }

  # ------------------------------------------------------------------
  # 12. EXPORT
  # ------------------------------------------------------------------
  if (export) {
    .step("export_omniseq", {
      export_omniseq(obj, outdir = outdir)
    })
  }

  # ------------------------------------------------------------------
  # SUMMARY
  # ------------------------------------------------------------------
  message("\n", strrep("=", 55))
  message("  run_omniseq_auto — PIPELINE SUMMARY")
  message(strrep("=", 55))
  for (nm in names(status)) {
    icon <- if (startsWith(status[[nm]], "PASS")) "\u2713" else "\u2717"
    message("  ", icon, "  ", nm, ": ", status[[nm]])
  }
  message(strrep("=", 55), "\n")

  return(invisible(obj))
}


# =============================================================================
# SECTION 2: GUIDED PIPELINE
# =============================================================================

#' Run the OmniSeq Pipeline Interactively
#'
#' Guides the user through the complete OmniSeq pipeline step by step,
#' prompting for parameter choices at each stage. Designed for exploratory
#' use in an interactive R session.
#'
#' @param count_matrix Matrix or data frame of raw counts.
#' @param metadata Data frame with sample metadata.
#' @param condition_col Character. Metadata column for main condition.
#' @param covariates Character vector. Covariate columns.
#' @param ref_level Character. Reference level.
#' @param species Character. Species name. Default: `"Homo sapiens"`.
#' @param org_db Character. Annotation database. Default: `"org.Hs.eg.db"`.
#' @param te_obj Optional OmniSeq TE object.
#' @param outdir Character. Export directory. Default: `"./omniseq_output"`.
#'
#' @return A populated OmniSeq object.
#' @export
run_omniseq_guided <- function(count_matrix,
                               metadata,
                               condition_col,
                               covariates = NULL,
                               ref_level  = NULL,
                               species    = "Homo sapiens",
                               org_db     = "org.Hs.eg.db",
                               te_obj     = NULL,
                               outdir     = "./omniseq_output") {

  if (!base::interactive())
    stop("run_omniseq_guided() must be run in an interactive R session.",
         call. = FALSE)

  .confirm <- function(step_name) {
    cat("\n")
    cat(strrep("-", 55), "\n")
    cat("  Next step: ", step_name, "\n")
    cat(strrep("-", 55), "\n")
    ans <- utils::menu(
      c("Run this step",
        "Skip this step",
        "Repeat last step",
        "Stop pipeline"),
      title = "What would you like to do?")
    if (ans == 4L) {
      message("Pipeline stopped by user.")
      return("stop")
    }
    if (ans == 2L) return("skip")
    if (ans == 3L) return("repeat")
    return("run")
  }

  status <- list()
  obj    <- NULL

  # ------------------------------------------------------------------
  # SECTION 2A: INIT + DATA TYPE DETECTION
  # ------------------------------------------------------------------
  cat("\n")
  cat(strrep("=", 55), "\n")
  cat("  OmniSeq Guided Pipeline\n")
  cat(strrep("=", 55), "\n")

  is_te_data <- any(grepl(":", rownames(count_matrix)[1:min(20, nrow(count_matrix))],
                          fixed = TRUE))

  if (is_te_data) {
    cat("\n  [Auto-detected: TE data (locus IDs detected)]\n")
    cat("  Routing to TE analysis pipeline.\n")
    cat("  Protein-coding filter: skipped (not applicable to TE loci)\n")
    filter_coding <- FALSE
  } else {
    filter_coding <- utils::menu(
      c("Yes — filter to protein-coding genes (recommended for GO enrichment)",
        "No  — keep all features"),
      title = "Filter to protein-coding genes?"
    ) == 1L
  }

  obj <- tryCatch({
    init_omniseq(
      count_matrix  = count_matrix,
      metadata      = metadata,
      condition_col = condition_col,
      covariates    = covariates,
      filter_coding = filter_coding,
      species       = species,
      org_db        = org_db,
      ref_level     = ref_level
    )
  }, error = function(e) {
    message("Initialisation failed: ", e$message)
    return(NULL)
  })
  if (is.null(obj)) return(invisible(NULL))
  print(obj)

  # BRANCH: TE DATA
  if (obj$params$data_type == "TE") {
    return(invisible(
      .run_omniseq_guided_te(obj       = obj,
                             condition_col = condition_col,
                             status    = status,
                             .confirm  = .confirm,
                             outdir    = outdir,
                             gene_obj  = te_obj)
    ))
  }

  # ------------------------------------------------------------------
  # SECTION 2B: QC
  # ------------------------------------------------------------------
  decision <- .confirm("Quality Control")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    tryCatch(grDevices::dev.off(), error = function(e) invisible(NULL))
    obj <- tryCatch({
      o <- run_omni_qc(obj)
      plot_omni_qc_basics(o, plot_type = "library_size")
      plot_omni_qc_basics(o, plot_type = "detected")
      plot_omni_qc_advanced(o, plot_type = "distance")
      plot_omni_qc_advanced(o, plot_type = "pca_compare")
      plot_omni_qc_advanced(o, plot_type = "mean_var")
      status[["QC"]] <- "PASS"
      o
    }, error = function(e) { status[["QC"]] <<- paste("FAIL:", e$message); obj })
  }

  # ------------------------------------------------------------------
  # SECTION 2C: DIFFERENTIAL EXPRESSION + CUSTOM CONTRASTS
  # ------------------------------------------------------------------
  decision <- .confirm("Differential Expression (DESeq2)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    shrink_choice <- utils::menu(
      c("apeglm (recommended for named coefficients)",
        "ashr   (more flexible, works with any contrast)",
        "none   (no shrinkage)"),
      title = "LFC shrinkage method:"
    )
    shrinkage <- c("apeglm", "ashr", "none")[shrink_choice]

    obj <- tryCatch({
      o <- run_omni_de(obj, shrinkage = shrinkage)
      plot_omni_pca(o)
      plot_omni_volcano(o, lfc_thresh = 1, top_n = 15)
      plot_omni_heatmap(o, ntop = 50)
      status[["DE"]] <- "PASS"
      o
    }, error = function(e) { status[["DE"]] <<- paste("FAIL:", e$message); obj })
  }

  meta       <- omni_metadata(obj)
  meta_cols  <- colnames(meta)

  add_contrasts <- utils::menu(
    c("Yes", "No"),
    title = "Add custom pairwise contrasts?"
  ) == 1L

  while (add_contrasts && !is.null(obj$model_info$dds)) {
    col_idx <- utils::menu(meta_cols, title = "Select metadata column:")
    col     <- meta_cols[col_idx]
    grps    <- unique(as.character(meta[[col]]))
    grps    <- grps[!is.na(grps)]
    g1_idx  <- utils::menu(grps, title = "Group 1 (numerator):")
    g2_idx  <- utils::menu(grps, title = "Group 2 (baseline):")
    g1      <- grps[g1_idx]
    g2      <- grps[g2_idx]

    obj <- tryCatch({
      run_omni_contrast(obj, column = col, group1 = g1, group2 = g2)
    }, error = function(e) {
      message("Contrast failed: ", e$message); obj
    })
    cts <- omni_contrasts(obj)
    if (length(cts) > 0) obj$active_contrast <- cts[1]

    add_contrasts <- utils::menu(c("Yes", "No"),
                                 title = "Add another contrast?") == 1L
  }

  cts <- omni_contrasts(obj)
  if (length(cts) > 0) obj$active_contrast <- cts[1]

  # ------------------------------------------------------------------
  # SECTION 2D: GO ENRICHMENT (ORA)
  # ── FIX: Main enrichment preserved across per-contrast loop
  # ------------------------------------------------------------------
  decision <- .confirm("GO Enrichment (ORA)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    lfc_thresh <- as.numeric(readline(
      prompt = "LFC threshold for ORA gene list [default: 0]: "))
    if (is.na(lfc_thresh)) lfc_thresh <- 0

    suppressPackageStartupMessages(
      do.call(library, list(org_db, character.only = TRUE))
    )

    # Run ORA on main contrast
    main_ct <- omni_contrasts(obj)[1]
    obj$active_contrast <- main_ct
    obj <- tryCatch({
      o <- run_omni_enrichment(obj, lfc_thresh = lfc_thresh)
      if (!is.null(o$enrichment) &&
          nrow(as.data.frame(o$enrichment)) > 0) {
        plot_omni_ora(o, plot_type = "dot")
        plot_omni_ora(o, plot_type = "cnet")
        plot_omni_ora(o, plot_type = "emap")
        plot_omni_compare_clusters(o)
      }
      status[["ORA"]] <- "PASS"
      o
    }, error = function(e) { status[["ORA"]] <<- paste("FAIL:", e$message); obj })

    # ── FIX: Save main enrichment + GSEA BEFORE per-contrast loop ──
    main_enrichment <- obj$enrichment
    main_gsea       <- obj$advanced_gsea

    # Run ORA on each additional contrast
    all_cts <- omni_contrasts(obj)
    extra_cts <- all_cts[!grepl("excl\\.", all_cts)]
    extra_cts <- setdiff(extra_cts, main_ct)
    for (ct in extra_cts) {
      obj$active_contrast <- ct
      obj <- tryCatch({
        o2 <- run_omni_enrichment(obj, lfc_thresh = lfc_thresh)
        n_t <- nrow(as.data.frame(o2$enrichment))
        message("[ORA] Contrast '", ct, "': ", n_t, " terms found.")
        if (n_t > 0) {
          plot_omni_ora(o2, plot_type = "dot")
          plot_omni_ora(o2, plot_type = "emap")
        }
        o2
      }, error = function(e) {
        message("[ORA] Could not run for '", ct, "': ", e$message)
        obj
      })
    }

    # ── FIX: Restore main enrichment so export detects has_ora=TRUE ──
    obj$enrichment    <- main_enrichment
    obj$advanced_gsea <- main_gsea
    obj$active_contrast <- main_ct

  }

  # ------------------------------------------------------------------
  # SECTION 2E: GSEA + ssGSEA
  # ------------------------------------------------------------------
  decision <- .confirm("GSEA (MSigDB)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    cat("\nMSigDB categories:\n")
    cat("  H  = Hallmark (50 curated gene sets — recommended first pass)\n")
    cat("  C2 = Curated (KEGG, Reactome, WikiPathways)\n")
    cat("  C5 = Ontology (GO BP, CC, MF)\n")
    cat("  C7 = Immunologic signatures (highly relevant for immune diseases)\n")
    cat("  Enter multiple separated by comma e.g.: H,C7\n\n")
    gsea_input <- readline(prompt = "Category/ies [default: H]: ")
    if (!nzchar(trimws(gsea_input))) gsea_input <- "H"
    gsea_cats <- toupper(trimws(strsplit(gsea_input, ",")[[1]]))

    for (cat_i in gsea_cats) {
      obj <- tryCatch({
        o <- run_omni_gsea(obj, category = cat_i)
        n_g <- nrow(as.data.frame(o$advanced_gsea))
        message("[GSEA] Category ", cat_i, ": ", n_g, " terms found.")
        if (n_g > 0) {
          plot_omni_gsea(o, show_n = 10)
          if (n_g >= 3)
            plot_omni_leading_edge(o, n_pathways = min(5, n_g))
        }
        status[[paste0("GSEA_", cat_i)]] <- "PASS"
        o
      }, error = function(e) {
        status[[paste0("GSEA_", cat_i)]] <<- paste("FAIL:", e$message)
        obj
      })
    }
  }

  decision <- .confirm("ssGSEA (per-sample pathway scores)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    obj <- tryCatch({
      o <- run_omni_ssgsea(obj)
      plot_omni_ssgsea_heatmap(o)
      status[["ssGSEA"]] <- "PASS"
      o
    }, error = function(e) { status[["ssGSEA"]] <<- paste("FAIL:", e$message); obj })
  }

  # ------------------------------------------------------------------
  # SECTION 2F: TF ENRICHMENT
  # ------------------------------------------------------------------
  decision <- .confirm("Transcription Factor Activity (decoupleR)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    db_choice <- utils::menu(
      c("CollecTRI (newer, more complete)",
        "DoRothEA  (classic, confidence-filtered)",
        "Both"),
      title = "TF database:"
    )
    db <- c("collectri", "dorothea", "both")[db_choice]

    obj <- tryCatch({
      o <- run_omni_tf(obj, database = db, method = "ulm")

      main_ct   <- o$active_contrast
      other_cts <- setdiff(
        omni_contrasts(o),
        grep("excl|mask|isgs|schoggins", omni_contrasts(o),
             ignore.case = TRUE, value = TRUE)
      )
      other_cts <- setdiff(other_cts, main_ct)

      for (ct in other_cts) {
        o$active_contrast <- ct
        o <- tryCatch(
          run_omni_tf(o, database = db, method = "ulm"),
          error = function(e) {
            message("[TF] Skipping contrast '", ct, "': ", e$message)
            o
          }
        )
      }
      o$active_contrast <- main_ct

      plot_omni_tf_activity(o, n = 15)
      plot_omni_tf_network(o, plot_type = "static")

      tryCatch(
        plot_omni_tf_genotype(o),
        error = function(e)
          message("[TF] Genotype plot skipped: ", e$message)
      )

      status[["TF"]] <- "PASS"
      o
    }, error = function(e) { status[["TF"]] <<- paste("FAIL:", e$message); obj })
  }

  # ------------------------------------------------------------------
  # SECTION 2G: DECONVOLUTION
  # ------------------------------------------------------------------
  decision <- .confirm("Cell Type Deconvolution (immunedeconv)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    method_choice <- utils::menu(
      c("xCell    (broad immune + stromal, 39 cell types)",
        "EPIC     (tumour microenvironment)",
        "QUANTISEQ (immune fraction, 10 cell types)"),
      title = "Deconvolution algorithm:"
    )
    deconv_method <- c("xcell", "epic", "quantiseq")[method_choice]

    obj <- tryCatch({
      o <- run_omni_deconvolution(obj, method = "immunedeconv",
                                  deconv_method = deconv_method)
      plot_omni_deconv_heatmap(o)
      status[["Deconv"]] <- "PASS"
      o
    }, error = function(e) { status[["Deconv"]] <<- paste("FAIL:", e$message); obj })
  }

  # ------------------------------------------------------------------
  # SECTION 2G2: SC REFERENCE DECONVOLUTION (optional)
  # Runs if user provides a SingleCellExperiment reference.
  # Auto-detects raw vs pseudo-bulk and routes accordingly.
  # Layer 1: cell composition from healthy reference
  # Layer 2: disease activation states from untreated reference
  # ------------------------------------------------------------------
  decision <- .confirm("scRNA-seq Reference Deconvolution (requires SCE object)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {

    # Ask user how to provide the SCE reference
    cat("\n")
    cat(strrep("-", 55), "\n")
    cat("  A SingleCellExperiment reference object is required.\n")
    cat("  This can be raw scRNA-seq or pseudo-bulk (donor x cell type).\n")
    cat(strrep("-", 55), "\n")

    sce_source <- utils::menu(
      c("SCE object is already loaded in my R session",
        "Load SCE from an RDS file on disk",
        "Skip — run sc reference deconv later"),
      title = "How would you like to provide the SCE reference?"
    )

    active_sce <- NULL

    if (sce_source == 1L) {
      sce_name <- trimws(readline(
        prompt = "Name of the SCE object in your session [e.g. savi_sce]: "))
      active_sce <- tryCatch(
        get(sce_name, envir = globalenv()),
        error = function(e) NULL)
      if (is.null(active_sce)) {
        message("[SC Deconv] Object '", sce_name, "' not found.")
      } else if (!inherits(active_sce, "SingleCellExperiment")) {
        message("[SC Deconv] '", sce_name,
                "' is not a SingleCellExperiment (class: ",
                class(active_sce)[1], ")")
        active_sce <- NULL
      } else {
        message("[SC Deconv] Found: ", ncol(active_sce), " profiles, ",
                nrow(active_sce), " genes")
      }

    } else if (sce_source == 2L) {
      rds_path <- trimws(readline(prompt = "Path to SCE RDS file: "))
      active_sce <- tryCatch({
        o <- readRDS(rds_path)
        if (!inherits(o, "SingleCellExperiment")) {
          message("[SC Deconv] Loaded object is not SCE (class: ",
                  class(o)[1], ")")
          NULL
        } else {
          message("[SC Deconv] Loaded: ", ncol(o), " profiles, ",
                  nrow(o), " genes")
          o
        }
      }, error = function(e) {
        message("[SC Deconv] Could not load: ", e$message)
        NULL
      })
    }

    if (!is.null(active_sce)) {
      # Show available metadata columns
      sce_cd <- as.data.frame(
        SummarizedExperiment::colData(active_sce))
      cat("\nAvailable metadata columns:\n")
      for (i in seq_along(colnames(sce_cd)))
        cat("  ", i, ": ", colnames(sce_cd)[i], " (",
            length(unique(sce_cd[[i]])), " unique)\n", sep = "")

      ct_col_idx <- utils::menu(colnames(sce_cd),
                                title = "Which column has cell type labels?")
      ct_col <- colnames(sce_cd)[ct_col_idx]

      smp_col_idx <- utils::menu(colnames(sce_cd),
                                 title = "Which column has donor/sample IDs?")
      smp_col <- colnames(sce_cd)[smp_col_idx]

      grp_col_idx <- utils::menu(colnames(sce_cd),
                                 title = "Which column has condition/group labels?")
      grp_col <- colnames(sce_cd)[grp_col_idx]

      # Show available groups
      cat("\nAvailable groups in '", grp_col, "':\n")
      print(table(sce_cd[[grp_col]]))

      ref_grp <- trimws(readline(
        prompt = "Reference group for marker discovery [e.g. Control]: "))

      ref_lbl <- trimws(readline(
        prompt = "Label for this reference [e.g. SAVI scRNA-seq]: "))
      if (!nzchar(ref_lbl)) ref_lbl <- "scRNA-seq reference"

      # Layer 1: Cell composition
      obj <- tryCatch({
        o <- run_omni_sc_deconv(obj,
                                sc_ref        = active_sce,
                                cell_type_col = ct_col,
                                sample_col    = smp_col,
                                group_col     = grp_col,
                                ref_group     = ref_grp,
                                ref_label     = ref_lbl)
        plot_omni_sc_heatmap(o, layer = "composition")
        plot_omni_sc_genotype(o, layer = "composition")
        plot_omni_sc_correlation(o, feature = "IRGS",
                                 layer = "composition")
        plot_omni_sc_profile(o, layer = "composition")
        status[["SC_Deconv_L1"]] <- "PASS"
        o
      }, error = function(e) {
        status[["SC_Deconv_L1"]] <<- paste("FAIL:", e$message)
        message("[SC Deconv Layer 1] Failed: ", e$message)
        obj
      })

      # Layer 2: Disease activation states
      run_l2 <- utils::menu(
        c("Yes — score disease activation states (untreated vs healthy)",
          "No — skip Layer 2"),
        title = "Run disease activation state scoring (Layer 2)?"
      ) == 1L

      if (run_l2) {
        cat("\nAvailable groups for disease reference:\n")
        print(table(sce_cd[[grp_col]]))
        dis_grp <- trimws(readline(
          prompt = "Disease group [e.g. SAVI_untreated]: "))

        obj <- tryCatch({
          o <- run_omni_disease_state(obj,
                                      sc_ref        = active_sce,
                                      disease_group = dis_grp,
                                      healthy_group = ref_grp,
                                      cell_type_col = ct_col,
                                      group_col     = grp_col)
          plot_omni_sc_heatmap(o, layer = "activation")
          plot_omni_sc_genotype(o, layer = "activation")
          plot_omni_sc_profile(o, layer = "activation")
          tryCatch(plot_omni_sc_two_layer(o),
                   error = function(e)
                     message("[SC Deconv] Two-layer plot: ", e$message))
          status[["SC_Deconv_L2"]] <- "PASS"
          o
        }, error = function(e) {
          status[["SC_Deconv_L2"]] <<- paste("FAIL:", e$message)
          message("[SC Deconv Layer 2] Failed: ", e$message)
          obj
        })
      }

      # Benchmark vs xCell
      if (!is.null(obj$deconv) && !is.null(obj$sc_deconv)) {
        tryCatch(plot_omni_sc_benchmark(obj),
                 error = function(e)
                   message("[SC Deconv] Benchmark: ", e$message))
      }
    } else {
      status[["SC_Deconv_L1"]] <- "SKIPPED — no SCE provided"
    }
  }

  # ------------------------------------------------------------------
  # SECTION 2H: STRING
  # ------------------------------------------------------------------
  decision <- .confirm("STRING Protein Interaction Network")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    n_genes <- as.integer(readline(
      prompt = "Number of top DEGs to query [default: 50]: "))
    if (is.na(n_genes)) n_genes <- 50L

    score_thresh <- as.integer(readline(
      prompt = "Minimum STRING score 0-1000 [default: 400]: "))
    if (is.na(score_thresh)) score_thresh <- 400L

    obj <- tryCatch({
      o <- run_omni_string(obj, n = n_genes, string_score = score_thresh)
      plot_omni_string_network(o, plot_type = "static")
      status[["STRING"]] <- "PASS"
      o
    }, error = function(e) { status[["STRING"]] <<- paste("FAIL:", e$message); obj })
  }

  # ------------------------------------------------------------------
  # SECTION 2I: sPLS-DA
  # ------------------------------------------------------------------
  decision <- .confirm("sPLS-DA Dimensionality Reduction")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    obj <- tryCatch({
      o <- run_omni_splsda(obj, group_col = condition_col)
      key <- names(o$splsda)[1]
      plot_omni_splsda_scores(o, key = key)
      plot_omni_splsda_vip(o,    key = key)
      plot_omni_splsda_roc(o,    key = key, comp = 1)
      status[["sPLS-DA"]] <- "PASS"
      o
    }, error = function(e) { status[["sPLS-DA"]] <<- paste("FAIL:", e$message); obj })

    # Visualize additional contrast DE results
    extra_cts_vis <- setdiff(omni_contrasts(obj), omni_contrasts(obj)[1])
    extra_cts_vis <- extra_cts_vis[!grepl("excl\\.", extra_cts_vis)]
    for (ct_vis in extra_cts_vis) {
      obj$active_contrast <- ct_vis
      message("\n[Pipeline] Visualizing contrast: ", ct_vis)
      tryCatch({
        plot_omni_volcano(obj, lfc_thresh = 1, top_n = 15)
        plot_omni_heatmap(obj, ntop = 50)
        plot_omni_top_genes(obj, n = 12)
        plot_omni_ma(obj)
      }, error = function(e) {
        message("[Pipeline] Visualization failed for '", ct_vis, "': ", e$message)
      })
    }
    obj$active_contrast <- omni_contrasts(obj)[1]

    # Run sPLS-DA for each additional contrast
    extra_cts <- setdiff(omni_contrasts(obj),
                         omni_contrasts(obj)[1])
    extra_cts <- extra_cts[!grepl("excl\\.", extra_cts)]
    for (ct in extra_cts) {
      parts <- strsplit(ct, " vs ")[[1]]
      if (length(parts) == 2) {
        meta_tmp <- omni_metadata(obj)
        found_col <- NULL
        for (col in colnames(meta_tmp)) {
          vals <- unique(as.character(meta_tmp[[col]]))
          if (all(parts %in% vals)) { found_col <- col; break }
        }
        if (!is.null(found_col)) {
          obj <- tryCatch({
            o2 <- run_omni_splsda(obj,
                                  group_col      = found_col,
                                  exclude_groups = setdiff(
                                    unique(as.character(meta_tmp[[found_col]])),
                                    parts))
            key2 <- paste0(found_col, "_",
                           paste(sort(parts), collapse = "-"))
            if (key2 %in% names(o2$splsda)) {
              plot_omni_splsda_scores(o2, key = key2)
              plot_omni_splsda_vip(o2,    key = key2)
            }
            message("[sPLS-DA] Completed for contrast: ", ct)
            o2$active_contrast <- omni_contrasts(o2)[1]
            o2
          }, error = function(e) {
            message("[sPLS-DA] Could not run for ", ct, ": ", e$message)
            obj
          })
        }
      }
    }
  }

  # ------------------------------------------------------------------
  # SECTION 2J: WGCNA
  # ------------------------------------------------------------------
  decision <- .confirm("WGCNA Co-expression Network")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    n_genes_wgcna <- as.integer(readline(
      prompt = "Number of top variable genes [default: 3000]: "))
    if (is.na(n_genes_wgcna)) n_genes_wgcna <- 3000L

    obj <- tryCatch({
      o <- run_omni_wgcna(obj, n_top_genes = n_genes_wgcna, key = "full")
      plot_omni_wgcna_dendrogram(o,    key = "full")
      plot_omni_wgcna_trait_heatmap(o, key = "full")
      status[["WGCNA"]] <- "PASS"
      o
    }, error = function(e) { status[["WGCNA"]] <<- paste("FAIL:", e$message); obj })
  }

  # ------------------------------------------------------------------
  # SECTION 2K: TE SENSOR SCORING
  # ------------------------------------------------------------------
  if (!is.null(te_obj)) {
    decision <- .confirm("TE Sensor Scoring (requires TE object)")
    if (decision == "stop") return(invisible(obj))
    if (decision == "run") {
      obj <- tryCatch({
        o <- run_omni_sensor_score(obj, te_obj = te_obj,
                                   group_col = condition_col)
        plot_omni_sensor_scores(o)
        plot_omni_sensor_irgs_corr(o)
        plot_omni_sensor_heatmap(o)
        status[["Sensors"]] <- "PASS"
        o
      }, error = function(e) { status[["Sensors"]] <<- paste("FAIL:", e$message); obj })
    }
  }

  # ------------------------------------------------------------------
  # SECTION 2L: EXPORT
  # ------------------------------------------------------------------
  decision <- .confirm("Export results (CSV tables + PNG plots)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    tryCatch({
      export_omniseq(obj, outdir = outdir)
      status[["Export"]] <- "PASS"
    }, error = function(e) {
      status[["Export"]] <- paste("FAIL:", e$message)
    })
  }

  # ------------------------------------------------------------------
  # SECTION 2M: MASKED ANALYSIS PROMPT
  # ------------------------------------------------------------------
  cat("\n")
  cat(strrep("-", 55), "\n")
  cat("  Would you like to run a masked analysis?\n")
  cat("  Results will be included in the HTML report.\n")
  cat("  (reruns DE plots, ORA, GSEA, sPLS-DA, WGCNA\n")
  cat("   after removing a gene signature)\n")
  cat(strrep("-", 55), "\n")

  mask_choice <- utils::menu(
    c("isgs_28 — mask the 28-gene NIH interferon panel",
      "isgs_extended — mask the 51-gene extended IFN signature",
      "isgs_schoggins — mask the 288-gene Schoggins ISG list",
      "ribosomal — mask ribosomal protein genes",
      "immunoglobulin — mask immunoglobulin variable region genes",
      "Custom — enter gene names manually",
      "No masking — skip"),
    title = "Run masked analysis?"
  )

  if (mask_choice < 7L) {
    mask_sig <- c("isgs_28", "isgs_extended", "isgs_schoggins",
                  "ribosomal", "immunoglobulin", NULL)[mask_choice]
    custom_features <- NULL
    if (mask_choice == 6L) {
      mask_sig    <- NULL
      raw_input   <- readline(prompt = "Enter gene names separated by commas: ")
      custom_features <- trimws(strsplit(raw_input, ",")[[1]])
      custom_features <- custom_features[nzchar(custom_features)]
      if (length(custom_features) == 0) {
        message("[Masking] No genes entered — skipping masked analysis.")
        mask_choice <- 7L
      } else {
        message("[Masking] Custom gene list: ",
                length(custom_features), " genes.")
      }
    }

    obj <- tryCatch({
      run_omniseq_masked(
        obj             = obj,
        signature       = mask_sig,
        features        = custom_features,
        gsea_categories = "H",
        n_wgcna         = 3000L,
        org_db          = org_db,
        outdir          = NULL
      )
    }, error = function(e) {
      message("[Masking] Failed: ", e$message)
      obj
    })
    status[["Masking"]] <- "PASS"
  }

  # ------------------------------------------------------------------
  # SECTION 2N: HTML REPORT + FINAL SUMMARY
  # ------------------------------------------------------------------
  decision <- .confirm("Generate HTML report")
  if (decision != "stop" && decision == "run") {
    tryCatch({
      export_omniseq_report(obj, outdir = outdir)
      status[["Report"]] <- "PASS"
    }, error = function(e) {
      status[["Report"]] <- paste("FAIL:", e$message)
    })
  }

  cat("\n")
  cat(strrep("=", 55), "\n")
  cat("  OmniSeq Guided Pipeline — COMPLETE\n")
  cat(strrep("=", 55), "\n")
  for (nm in names(status)) {
    icon <- if (startsWith(status[[nm]], "PASS")) "\u2713" else "\u2717"
    cat("  ", icon, " ", nm, ": ", status[[nm]], "\n", sep = "")
  }
  cat(strrep("=", 55), "\n\n")
  cat("  Tip: To run masked analysis later:\n")
  cat("  gene_obj <- run_omniseq_masked(gene_obj, signature = 'isgs_28')\n\n")

  print(obj)

  return(invisible(obj))
}


# =============================================================================
# SECTION 3: POST-PIPELINE MASKING ANALYSIS
# =============================================================================

#' Run Focused Masked Analysis After the Full Pipeline
#'
#' After running run_omniseq_guided() or run_omniseq_auto() on the full
#' gene set, this function applies a gene signature mask and reruns only
#' the steps that are biologically meaningful to compare:
#' - DE results (volcano, heatmap, top genes)
#' - GO ORA enrichment
#' - GSEA
#' - sPLS-DA
#' - WGCNA
#'
#' @param obj A completed OmniSeq object (after run_omniseq_guided).
#' @param signature Character. Built-in signature to mask. Options:
#'   "isgs_28", "isgs_extended", "ribosomal", "mitochondrial",
#'   "immunoglobulin". Default: "isgs_28".
#' @param features Character vector. Additional genes to mask beyond the
#'   signature. Default: NULL.
#' @param mask_label Character. Label for the masked contrast.
#'   Default: auto-generated.
#' @param lfc_thresh Numeric. LFC threshold for ORA. Default: 0.
#' @param gsea_categories Character vector. MSigDB categories for GSEA.
#'   Default: c("H") (Hallmark).
#' @param n_wgcna Integer. Number of top variable genes for masked WGCNA.
#'   Default: 3000.
#' @param org_db Character. Annotation database. Default: "org.Hs.eg.db".
#' @param outdir Character. Directory for masked results export.
#'   Default: NULL (uses same directory as original export).
#'
#' @return The OmniSeq object updated with masked results.
#'
#' @details
#' Steps NOT rerun (masking does not change them):
#' - QC (cell-level metrics unchanged)
#' - Deconvolution (cell type fractions unchanged)
#' - STRING network (topology unchanged)
#' - TF activity (regulon-based scores unchanged)
#' - ssGSEA (per-sample pathway scores use all genes)
#'
#' Results are stored inside the original object under masked contrast keys
#' and are automatically included in the HTML report.
#' @export
run_omniseq_masked <- function(obj,
                               signature       = "isgs_28",
                               features        = NULL,
                               mask_label      = NULL,
                               lfc_thresh      = 0,
                               gsea_categories = "H",
                               n_wgcna         = 3000L,
                               org_db          = "org.Hs.eg.db",
                               outdir          = NULL) {

  validate_omniseq(obj)

  if (!base::interactive() &&
      is.null(obj$model_info$dds))
    stop("No fitted model found. Run run_omni_de() first.", call. = FALSE)

  # ------------------------------------------------------------------
  # 1. APPLY MASK
  # ------------------------------------------------------------------
  cat("\n")
  cat(strrep("=", 55), "\n")
  cat("  OmniSeq — Masked Analysis\n")
  cat(strrep("=", 55), "\n")

  obj <- tryCatch({
    if (!is.null(features) && !is.null(signature)) {
      run_omni_mask(obj, signature = signature,
                    features = features, mask_label = mask_label)
    } else if (!is.null(signature)) {
      run_omni_mask(obj, signature = signature,
                    mask_label = mask_label)
    } else if (!is.null(features)) {
      run_omni_mask(obj, features = features,
                    mask_label = mask_label)
    } else {
      stop("Provide at least one of: signature or features.", call. = FALSE)
    }
  }, error = function(e) {
    stop("Masking failed: ", e$message, call. = FALSE)
  })

  masked_key <- obj$active_contrast
  message("\n[run_omniseq_masked] Masked contrast: '", masked_key, "'")

  res_m      <- as.data.frame(omni_results(obj))
  res_m$padj <- ifelse(is.na(res_m$padj), 1, res_m$padj)
  n_sig_m    <- sum(res_m$padj < 0.05)
  message("[run_omniseq_masked] Significant genes after masking: ", n_sig_m)

  status <- list()

  # ------------------------------------------------------------------
  # 2. DE VISUALIZATIONS ON MASKED RESULTS
  # ------------------------------------------------------------------
  message("\n[run_omniseq_masked] Step 1/5 — DE visualizations...")
  tryCatch({
    plot_omni_pca(obj)
    plot_omni_volcano(obj, lfc_thresh = 1, top_n = 15)
    plot_omni_heatmap(obj, ntop = 50)
    plot_omni_top_genes(obj, n = 12)
    plot_omni_ma(obj)
    plot_omni_pvalue_hist(obj)
    status[["DE_plots"]] <- "PASS"
  }, error = function(e) {
    status[["DE_plots"]] <<- paste("FAIL:", e$message)
    message("[run_omniseq_masked] DE plots failed: ", e$message)
  })

  # ------------------------------------------------------------------
  # 3. GO ORA ON MASKED RESULTS
  # ------------------------------------------------------------------
  message("\n[run_omniseq_masked] Step 2/5 — GO ORA...")
  tryCatch({
    suppressPackageStartupMessages(
      do.call(library, list(org_db, character.only = TRUE))
    )
    obj <- run_omni_enrichment(obj, lfc_thresh = lfc_thresh)
    n_ora <- nrow(as.data.frame(obj$enrichment))
    message("[run_omniseq_masked] ORA found ", n_ora, " terms.")
    if (n_ora > 0) {
      plot_omni_ora(obj, plot_type = "dot")
      plot_omni_ora(obj, plot_type = "cnet")
      plot_omni_ora(obj, plot_type = "emap")
      plot_omni_compare_clusters(obj)
    }
    status[["ORA"]] <- paste("PASS —", n_ora, "terms")
  }, error = function(e) {
    status[["ORA"]] <<- paste("FAIL:", e$message)
    message("[run_omniseq_masked] ORA failed: ", e$message)
  })

  # ------------------------------------------------------------------
  # 4. GSEA ON MASKED RESULTS
  # ------------------------------------------------------------------
  message("\n[run_omniseq_masked] Step 3/5 — GSEA...")
  for (cat_i in gsea_categories) {
    tryCatch({
      obj <- run_omni_gsea(obj, category = cat_i)
      n_g <- nrow(as.data.frame(obj$advanced_gsea))
      message("[run_omniseq_masked] GSEA ", cat_i, ": ", n_g, " terms.")
      if (n_g > 0) {
        plot_omni_gsea(obj, show_n = 10)
        if (n_g >= 3)
          plot_omni_leading_edge(obj, n_pathways = min(5, n_g))
      }
      status[[paste0("GSEA_", cat_i)]] <- paste("PASS —", n_g, "terms")
    }, error = function(e) {
      status[[paste0("GSEA_", cat_i)]] <<- paste("FAIL:", e$message)
      message("[run_omniseq_masked] GSEA ", cat_i, " failed: ", e$message)
    })
  }

  # ------------------------------------------------------------------
  # 5. sPLS-DA ON MASKED RESULTS
  # ------------------------------------------------------------------
  message("\n[run_omniseq_masked] Step 4/5 — sPLS-DA...")
  tryCatch({
    cond_col <- obj$params$condition
    obj <- run_omni_splsda(obj, group_col = cond_col, use_masked = TRUE)
    splsda_keys <- names(obj$splsda)
    new_key     <- splsda_keys[length(splsda_keys)]
    plot_omni_splsda_scores(obj,   key = new_key)
    plot_omni_splsda_vip(obj,      key = new_key)
    plot_omni_splsda_loadings(obj, key = new_key, comp = 1)
    sig_masked <- omni_splsda_signature(obj, key = new_key, comp = 1)
    message("[run_omniseq_masked] sPLS-DA signature: ",
            nrow(sig_masked), " features (VIP >= 1)")
    ber_val <- tryCatch({
      ber_mat <- obj$splsda[[new_key]]$performance$error.rate$BER
      round(ber_mat[nrow(ber_mat), "max.dist"], 3)
    }, error = function(e) "unknown")
    status[["sPLS-DA"]] <- paste("PASS — BER =", ber_val)
  }, error = function(e) {
    status[["sPLS-DA"]] <<- paste("FAIL:", e$message)
    message("[run_omniseq_masked] sPLS-DA failed: ", e$message)
  })

  # ------------------------------------------------------------------
  # 6. WGCNA ON MASKED RESULTS
  # ------------------------------------------------------------------
  message("\n[run_omniseq_masked] Step 5/5 — WGCNA...")
  tryCatch({
    obj <- run_omni_wgcna(obj,
                          use_masked  = TRUE,
                          n_top_genes = n_wgcna,
                          key         = "masked")
    plot_omni_wgcna_dendrogram(obj,    key = "masked")
    plot_omni_wgcna_trait_heatmap(obj, key = "masked")
    plot_omni_wgcna_hub_network(obj,   key = "masked",
                                plot_type = "static")
    n_mod <- obj$wgcna[["masked"]]$n_modules
    status[["WGCNA"]] <- paste("PASS —", n_mod, "modules")
  }, error = function(e) {
    status[["WGCNA"]] <<- paste("FAIL:", e$message)
    message("[run_omniseq_masked] WGCNA failed: ", e$message)
  })

  # ------------------------------------------------------------------
  # 7. EXPORT MASKED RESULTS
  # ------------------------------------------------------------------
  if (!is.null(outdir)) {
    masked_dir <- file.path(outdir, "masked")
    tryCatch({
      export_omniseq(obj,
                     outdir = masked_dir,
                     prefix = paste0("masked_", gsub("[^A-Za-z0-9]", "_",
                                                     signature %||% "custom")))
      message("[run_omniseq_masked] Masked results exported to: ", masked_dir)
      status[["Export"]] <- "PASS"
    }, error = function(e) {
      status[["Export"]] <<- paste("FAIL:", e$message)
    })
  }

  # ------------------------------------------------------------------
  # 8. REGENERATE HTML REPORT WITH MASKED SECTIONS
  # ------------------------------------------------------------------
  if (!is.null(outdir)) {
    tryCatch({
      export_omniseq_report(obj,
                            outdir      = outdir,
                            report_file = "omniseq_report_with_masked.html",
                            title       = paste("OmniSeq Analysis — Full + Masked (",
                                                signature %||% "custom", ")"),
                            open        = TRUE)
      status[["Report"]] <- "PASS"
    }, error = function(e) {
      status[["Report"]] <<- paste("FAIL:", e$message)
    })
  }

  # ------------------------------------------------------------------
  # SUMMARY
  # ------------------------------------------------------------------
  cat("\n")
  cat(strrep("=", 55), "\n")
  cat("  run_omniseq_masked — COMPLETE\n")
  cat(strrep("=", 55), "\n")
  cat("  Masked contrast: '", masked_key, "'\n", sep = "")
  cat("  Significant genes remaining: ", n_sig_m, "\n\n", sep = "")
  for (nm in names(status)) {
    icon <- if (startsWith(status[[nm]], "PASS")) "\u2713" else "\u2717"
    cat("  ", icon, " ", nm, ": ", status[[nm]], "\n", sep = "")
  }
  cat(strrep("=", 55), "\n\n")

  obj$active_contrast <- omni_contrasts(obj)[1]

  return(invisible(obj))
}


# =============================================================================
# SECTION 4: TE-SPECIFIC GUIDED PIPELINE (internal)
# Called automatically when data_type == "TE" is detected
# =============================================================================

#' @keywords internal
.run_omniseq_guided_te <- function(obj, condition_col, status, .confirm,
                                   outdir, gene_obj = NULL) {

  cat("\n")
  cat(strrep("=", 55), "\n")
  cat("  OmniSeq TE Analysis Pipeline\n")
  cat(strrep("=", 55), "\n")
  cat("  Data type : TE (", obj$params$data_type, ")\n", sep = "")
  cat("  Features  : ", nrow(obj$data), " loci\n", sep = "")
  cat("  Samples   : ", ncol(obj$data), "\n", sep = "")
  cat("\n")
  cat("  Steps in this pipeline:\n")
  cat("    1. Quality Control\n")
  cat("    2. Differential Expression (DESeq2 on raw TE counts)\n")
  cat("    3. TE Family & Class Analysis\n")
  cat("    4. TE Visualization (volcano, heatmap, pie, density)\n")
  cat("    5. GREAT Genomic Enrichment (proximity-based GO)\n")
  cat("    6. sPLS-DA (does TE expression classify patients?)\n")
  cat("    7. WGCNA TE Co-expression Network\n")
  cat("    8. Multi-omics Integration (Sensor Scoring + DIABLO)\n")
  if (!is.null(gene_obj))
    cat("    9. TE-Gene Correlation (vs IRG-S / gene signatures)\n")
  cat("\n")
  cat("  [SKIPPED — not applicable to TE loci]:\n")
  cat("    ORA, GSEA, ssGSEA, TF activity, Deconvolution, STRING\n")
  cat(strrep("=", 55), "\n")

  # ------------------------------------------------------------------
  # 1. QC
  # ------------------------------------------------------------------
  decision <- .confirm("Quality Control")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    tryCatch(grDevices::dev.off(), error = function(e) invisible(NULL))
    obj <- tryCatch({
      o <- run_omni_qc(obj)
      plot_omni_qc_basics(o, plot_type = "library_size")
      plot_omni_qc_basics(o, plot_type = "detected")
      plot_omni_qc_advanced(o, plot_type = "distance")
      plot_omni_qc_advanced(o, plot_type = "mean_var")
      status[["QC"]] <- "PASS"
      o
    }, error = function(e) { status[["QC"]] <<- paste("FAIL:", e$message); obj })
  }

  # ------------------------------------------------------------------
  # 2. DIFFERENTIAL EXPRESSION
  # ------------------------------------------------------------------
  decision <- .confirm("Differential Expression (DESeq2 on TE loci)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    shrink_choice <- utils::menu(
      c("apeglm (recommended)",
        "ashr   (more flexible)",
        "none   (no shrinkage)"),
      title = "LFC shrinkage method:"
    )
    shrinkage <- c("apeglm", "ashr", "none")[shrink_choice]

    obj <- tryCatch({
      o <- run_omni_de(obj, shrinkage = shrinkage)
      plot_omni_pca(o)
      plot_omni_volcano(o, lfc_thresh = 1, top_n = 15)

      res_df <- as.data.frame(omni_results(o))
      res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)
      n_sig  <- sum(res_df$padj < 0.05)
      n_up   <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1)
      n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1)
      message("\n[TE DE] Significant loci (padj<0.05): ", n_sig)
      message("[TE DE] Derepressed  (LFC>1):          ", n_up)
      message("[TE DE] Silenced     (LFC<-1):          ", n_down)
      message("[TE DE] Fraction derepressed: ",
              round(n_up / max(nrow(res_df), 1) * 100, 2), "%")

      status[["DE"]] <- paste0("PASS — ", n_sig, " sig. loci (", n_up, " up)")
      o
    }, error = function(e) { status[["DE"]] <<- paste("FAIL:", e$message); obj })
  }

  # Optional custom contrasts
  meta      <- omni_metadata(obj)
  meta_cols <- colnames(meta)
  add_contrasts <- utils::menu(c("Yes", "No"),
                               title = "Add custom pairwise contrasts? (e.g. TLR7 vs UNC93B1)") == 1L

  while (add_contrasts && !is.null(obj$model_info$dds)) {
    col_idx <- utils::menu(meta_cols, title = "Select metadata column:")
    col     <- meta_cols[col_idx]
    grps    <- unique(as.character(meta[[col]]))
    grps    <- grps[!is.na(grps)]
    g1_idx  <- utils::menu(grps, title = "Group 1 (numerator):")
    g2_idx  <- utils::menu(grps, title = "Group 2 (baseline):")
    g1 <- grps[g1_idx]; g2 <- grps[g2_idx]
    obj <- tryCatch(
      run_omni_contrast(obj, column = col, group1 = g1, group2 = g2),
      error = function(e) { message("Contrast failed: ", e$message); obj })
    cts <- omni_contrasts(obj)
    if (length(cts) > 0) obj$active_contrast <- cts[1]
    add_contrasts <- utils::menu(c("Yes", "No"),
                                 title = "Add another contrast?") == 1L
  }
  cts <- omni_contrasts(obj)
  if (length(cts) > 0) obj$active_contrast <- cts[1]

  # ------------------------------------------------------------------
  # 3. TE FAMILY & CLASS ANALYSIS
  # ------------------------------------------------------------------
  decision <- .confirm("TE Family & Class Analysis")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    lfc_te <- as.numeric(readline(
      prompt = "LFC threshold for derepressed loci [default: 1]: "))
    if (is.na(lfc_te)) lfc_te <- 1

    obj <- tryCatch({
      o <- run_omni_te_family(obj, pval_thresh = 0.05, lfc_thresh = lfc_te)

      if (!is.null(o$te_family) && !is.null(o$te_family$by_family)) {
        fam_df <- o$te_family$by_family
        cat("\n  Top 10 derepressed TE families (by N_Sig):\n")
        top_up <- utils::head(fam_df[order(-fam_df$Mean_LFC), ], 10)
        for (i in seq_len(nrow(top_up))) {
          cat("    ", top_up$Family[i],
              " class=", top_up$Class[i],
              " mean_LFC=", round(top_up$Mean_LFC[i], 2),
              " n_sig=",    top_up$N_Sig[i],
              " n_loci=",   top_up$N_Loci[i], "\n")
        }
        cat("\n  Class-level summary:\n")
        cls_df <- o$te_family$by_class
        for (i in seq_len(nrow(cls_df))) {
          cat("    ", cls_df$Class[i],
              " n_loci=",   cls_df$N_Loci[i],
              " n_sig=",    cls_df$N_Sig[i],
              " frac_sig=", cls_df$Frac_Sig[i],
              " mean_LFC=", round(cls_df$Mean_LFC[i], 2), "\n")
        }
      }

      status[["TE_Family"]] <- "PASS"
      o
    }, error = function(e) {
      status[["TE_Family"]] <<- paste("FAIL:", e$message); obj })
  }

  # ------------------------------------------------------------------
  # 4. TE VISUALIZATION
  # ------------------------------------------------------------------
  decision <- .confirm("TE Visualization (volcano, heatmap, class plots, density)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    obj <- tryCatch({
      plot_omni_te_volcano(obj, lfc_thresh = 1, top_n = 15)
      plot_omni_te_class_pie(obj)
      plot_omni_te_family(obj, n_families = 20)
      plot_omni_te_heatmap(obj, n_families = 20)
      plot_omni_te_lfc_density(obj, min_loci = 10)
      status[["TE_Plots"]] <- "PASS"
      obj
    }, error = function(e) {
      status[["TE_Plots"]] <<- paste("FAIL:", e$message); obj })
  }

  # ------------------------------------------------------------------
  # 5. GREAT GENOMIC ENRICHMENT
  # ------------------------------------------------------------------
  decision <- .confirm("GREAT Genomic Enrichment (proximity-based GO)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    lfc_great <- as.numeric(readline(
      prompt = "Min LFC for GREAT input loci [default: 1]: "))
    if (is.na(lfc_great)) lfc_great <- 1

    cat("\n  GREAT requires a coordinate lookup table built from your GTF files.\n")
    cat("  If you have not built one yet, run this first:\n")
    cat("    te_coords <- omni_build_te_coords(gtf_files = c('path/to/te.gtf'))\n")
    cat("    saveRDS(te_coords, '~/te_coords.rds')\n\n")

    has_coords <- utils::menu(
      c("Yes — I have te_coords.rds saved on disk",
        "No  — skip GREAT for now"),
      title = "Do you have a te_coords file?"
    ) == 1L

    if (has_coords) {
      coords_path <- readline(prompt = "Path to te_coords.rds: ")
      coords_path <- trimws(coords_path)
      obj <- tryCatch({
        te_coords <- readRDS(coords_path)
        message("[GREAT] Loaded te_coords: ", nrow(te_coords), " loci")
        obj <- run_omni_te_proximity(obj,
                                     te_coords  = te_coords,
                                     lfc_thresh = lfc_great)
        obj <- run_omni_te_great(obj)
        plot_omni_te_great(obj, show_n = 15)
        status[["GREAT"]] <- "PASS"
        obj
      }, error = function(e) {
        status[["GREAT"]] <<- paste("FAIL:", e$message)
        message("[GREAT] Failed: ", e$message)
        obj
      })
    } else {
      status[["GREAT"]] <- "SKIPPED — no te_coords"
      message("[GREAT] Skipped. Build te_coords to enable this step.")
    }
  }

  # ------------------------------------------------------------------
  # 6. TE-GENE CORRELATION (if gene_obj provided)
  # ------------------------------------------------------------------
  if (!is.null(gene_obj)) {
    decision <- .confirm("TE-Gene Correlation (TE expression vs IRG-S / gene scores)")
    if (decision == "stop") return(invisible(obj))
    if (decision == "run") {
      obj <- tryCatch({
        cat("\n  Running TE-gene correlation — two passes:\n")
        cat("  Pass 1: TE load vs IRG-S (28-gene IFN panel)\n")
        cat("  Pass 2: TE load vs ISG-independent disease genes\n")
        cat("  Comparing both reveals whether TE-IFN coupling is\n")
        cat("  purely an IFN artifact or reflects broader biology.\n\n")

        message("[TE-Gene] Pass 1 — IRG-S Score (28-gene IFN panel)")
        plot_omni_te_gene_corr(
          gene_obj     = gene_obj,
          te_obj       = obj,
          marker_genes = OmniSeq::omniseq_isgs_28,
          te_class     = "ALL",
          score_name   = "IRG-S Score (full IFN panel)"
        )

        non_isg_genes <- tryCatch({
          g_res         <- as.data.frame(omni_results(gene_obj))
          g_res$padj    <- ifelse(is.na(g_res$padj), 1, g_res$padj)
          sig_all       <- rownames(g_res)[g_res$padj < 0.05]
          non_isg       <- setdiff(sig_all, OmniSeq::omniseq_isgs_28)
          message("[TE-Gene] Pass 2 — ", length(non_isg),
                  " ISG-independent significant genes as signature")
          non_isg
        }, error = function(e) {
          message("[TE-Gene] Could not extract non-ISG genes: ", e$message,
                  " — using extended ISG panel instead")
          setdiff(OmniSeq::omniseq_isgs_extended %||%
                    OmniSeq::omniseq_isgs_28, OmniSeq::omniseq_isgs_28)
        })

        if (length(non_isg_genes) >= 3) {
          plot_omni_te_gene_corr(
            gene_obj     = gene_obj,
            te_obj       = obj,
            marker_genes = non_isg_genes,
            te_class     = "ALL",
            score_name   = "ISG-independent disease score"
          )
        } else {
          message("[TE-Gene] Pass 2 skipped — fewer than 3 non-ISG DEGs available.")
        }

        te_classes <- unique(
          OmniSeq:::.omni_parse_te_ids(rownames(obj$data))$Class)
        te_classes <- te_classes[te_classes != "Unknown" & !is.na(te_classes)]
        if (length(te_classes) > 1) {
          message("[TE-Gene] Running per-class correlations (",
                  paste(te_classes, collapse=", "), ")")
          for (cls in te_classes) {
            tryCatch({
              plot_omni_te_gene_corr(
                gene_obj     = gene_obj,
                te_obj       = obj,
                marker_genes = OmniSeq::omniseq_isgs_28,
                te_class     = cls,
                score_name   = paste("IRG-S —", cls)
              )
            }, error = function(e)
              message("[TE-Gene] Skipping class ", cls, ": ", e$message))
          }
        }

        status[["TE_Gene_Corr"]] <- paste0(
          "PASS — IRG-S + ISG-independent (",
          length(non_isg_genes), " non-ISG genes)")
        obj
      }, error = function(e) {
        status[["TE_Gene_Corr"]] <<- paste("FAIL:", e$message); obj })
    }
  }

  # ------------------------------------------------------------------
  # 7. sPLS-DA
  # ------------------------------------------------------------------
  decision <- .confirm("sPLS-DA (does TE expression classify patients vs controls?)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    obj <- tryCatch({
      o <- run_omni_splsda(obj, group_col = condition_col)
      key <- names(o$splsda)[1]
      plot_omni_splsda_scores(o, key = key)
      plot_omni_splsda_vip(o, key = key)
      plot_omni_splsda_roc(o, key = key, comp = 1)

      ber <- tryCatch(
        round(o$splsda[[key]]$performance$error.rate$BER[1, "max.dist"], 3),
        error = function(e) "unknown")
      message("[TE sPLS-DA] BER = ", ber,
              "  (compare to gene BER=0.144 to assess TE added value)")
      status[["sPLS-DA"]] <- paste("PASS — BER =", ber)
      o
    }, error = function(e) { status[["sPLS-DA"]] <<- paste("FAIL:", e$message); obj })
  }

  # ------------------------------------------------------------------
  # 8. WGCNA TE CO-EXPRESSION NETWORK
  # ------------------------------------------------------------------
  decision <- .confirm("WGCNA TE Co-expression Network")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    n_loci_wgcna <- as.integer(readline(
      prompt = "Number of top variable TE loci [default: 3000]: "))
    if (is.na(n_loci_wgcna)) n_loci_wgcna <- 3000L

    obj <- tryCatch({
      o <- run_omni_wgcna(obj, n_top_genes = n_loci_wgcna, key = "full")
      plot_omni_wgcna_dendrogram(o, key = "full")
      plot_omni_wgcna_trait_heatmap(o, key = "full")
      status[["WGCNA"]] <- paste("PASS —", o$wgcna[["full"]]$n_modules, "modules")
      o
    }, error = function(e) { status[["WGCNA"]] <<- paste("FAIL:", e$message); obj })
  }
  # ------------------------------------------------------------------
  # 8b. SC REFERENCE DECONVOLUTION (optional, for gene_obj)
  # In the TE pipeline, SC deconv runs on gene_obj (not TE obj).
  # Only offered if gene_obj is available.
  # ------------------------------------------------------------------
  if (!is.null(gene_obj)) {
    decision <- .confirm("scRNA-seq Reference Deconvolution (on gene object)")
    if (decision == "stop") return(invisible(obj))
    if (decision == "run") {

      cat("\n  SC reference deconvolution runs on your gene_obj.\n")
      cat("  You need a SingleCellExperiment reference object.\n\n")

      sce_source <- utils::menu(
        c("SCE object is loaded in my session",
          "Load from RDS file",
          "Skip"),
        title = "Provide SCE reference?"
      )

      active_sce <- NULL
      if (sce_source == 1L) {
        sce_name <- trimws(readline(prompt = "SCE object name: "))
        active_sce <- tryCatch(
          get(sce_name, envir = globalenv()),
          error = function(e) NULL)
        if (!is.null(active_sce) &&
            !inherits(active_sce, "SingleCellExperiment"))
          active_sce <- NULL
      } else if (sce_source == 2L) {
        rds_path <- trimws(readline(prompt = "Path to SCE RDS: "))
        active_sce <- tryCatch({
          o <- readRDS(rds_path)
          if (inherits(o, "SingleCellExperiment")) o else NULL
        }, error = function(e) NULL)
      }

      if (!is.null(active_sce) && !is.null(gene_obj)) {
        sce_cd <- as.data.frame(
          SummarizedExperiment::colData(active_sce))
        cat("\nSCE metadata columns:\n")
        for (i in seq_along(colnames(sce_cd)))
          cat("  ", i, ": ", colnames(sce_cd)[i], "\n", sep = "")

        ct_col  <- colnames(sce_cd)[utils::menu(colnames(sce_cd),
                                                title = "Cell type column?")]
        smp_col <- colnames(sce_cd)[utils::menu(colnames(sce_cd),
                                                title = "Donor/sample column?")]
        grp_col <- colnames(sce_cd)[utils::menu(colnames(sce_cd),
                                                title = "Condition column?")]
        print(table(sce_cd[[grp_col]]))
        ref_grp <- trimws(readline(prompt = "Reference group: "))
        ref_lbl <- trimws(readline(prompt = "Label [default: scRNA-seq]: "))
        if (!nzchar(ref_lbl)) ref_lbl <- "scRNA-seq reference"

        gene_obj <- tryCatch({
          g <- run_omni_sc_deconv(gene_obj, sc_ref = active_sce,
                                  cell_type_col = ct_col,
                                  sample_col = smp_col,
                                  group_col = grp_col,
                                  ref_group = ref_grp,
                                  ref_label = ref_lbl)
          plot_omni_sc_heatmap(g, layer = "composition")
          plot_omni_sc_genotype(g, layer = "composition")
          status[["SC_Deconv"]] <- "PASS"
          g
        }, error = function(e) {
          status[["SC_Deconv"]] <<- paste("FAIL:", e$message)
          gene_obj
        })

        # Layer 2
        run_l2 <- utils::menu(
          c("Yes", "No"),
          title = "Run disease activation scoring (Layer 2)?") == 1L
        if (run_l2) {
          print(table(sce_cd[[grp_col]]))
          dis_grp <- trimws(readline(prompt = "Disease group: "))
          gene_obj <- tryCatch({
            g <- run_omni_disease_state(gene_obj, sc_ref = active_sce,
                                        disease_group = dis_grp,
                                        healthy_group = ref_grp,
                                        cell_type_col = ct_col,
                                        group_col = grp_col)
            plot_omni_sc_genotype(g, layer = "activation")
            tryCatch(plot_omni_sc_two_layer(g),
                     error = function(e) NULL)
            status[["SC_Deconv_L2"]] <- "PASS"
            g
          }, error = function(e) {
            status[["SC_Deconv_L2"]] <<- paste("FAIL:", e$message)
            gene_obj
          })
        }
      }
    }
  }

  # ------------------------------------------------------------------
  # 9. MULTI-OMICS INTEGRATION
  # ------------------------------------------------------------------
  decision <- .confirm("Multi-omics Integration (Sensor Scoring + DIABLO)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {

    active_gene_obj <- gene_obj

    if (!is.null(active_gene_obj)) {
      message("[Multi-omics] Using gene_obj provided at pipeline start.")
      message("  Features: ", nrow(active_gene_obj$data),
              " | Samples: ", ncol(active_gene_obj$data))
    }

    if (is.null(active_gene_obj)) {
      cat("\n")
      cat(strrep("-", 55), "\n")
      cat("  A completed gene OmniSeq object is required.\n")
      cat("  This should be gene_obj from the gene expression pipeline.\n")
      cat(strrep("-", 55), "\n")

      gene_source <- utils::menu(
        c("gene_obj is already loaded in my R session",
          "Load gene_obj from an RDS file on disk",
          "Skip — run multi-omics integration later"),
        title = "How would you like to provide the gene object?"
      )

      if (gene_source == 1L) {
        active_gene_obj <- tryCatch(
          get("gene_obj", envir = globalenv()),
          error = function(e) NULL
        )
        if (is.null(active_gene_obj)) {
          message("[Multi-omics] 'gene_obj' not found in global environment.")
          message("  Tip: pass it at pipeline start via te_obj = gene_obj")
          status[["Sensors"]] <- "SKIPPED — gene_obj not found in session"
          status[["DIABLO"]]  <- "SKIPPED — gene_obj not found in session"
          gene_source <- 3L
        } else {
          message("[Multi-omics] Found gene_obj in global environment.")
          message("  Features: ", nrow(active_gene_obj$data),
                  " | Samples: ", ncol(active_gene_obj$data))
        }

      } else if (gene_source == 2L) {
        rds_path <- trimws(readline(prompt = "Path to gene_obj RDS file: "))
        active_gene_obj <- tryCatch({
          o <- readRDS(rds_path)
          message("[Multi-omics] Loaded gene_obj from: ", rds_path)
          message("  Features: ", nrow(o$data), " | Samples: ", ncol(o$data))
          o
        }, error = function(e) {
          message("[Multi-omics] Could not load RDS: ", e$message)
          NULL
        })
        if (is.null(active_gene_obj)) {
          status[["Sensors"]] <- "SKIPPED — could not load gene_obj"
          status[["DIABLO"]]  <- "SKIPPED — could not load gene_obj"
          gene_source <- 3L
        }
      }

      if (gene_source == 3L) {
        message("[Multi-omics] Skipped.")
        status[["Sensors"]] <- status[["Sensors"]] %||% "SKIPPED"
        status[["DIABLO"]]  <- status[["DIABLO"]]  %||% "SKIPPED"
        active_gene_obj <- NULL
      }
    }

    # --- Sensor Scoring ---
    if (!is.null(active_gene_obj)) {
      cat("\n[Multi-omics] Running Sensor Scoring...\n")
      updated_gene_obj <- tryCatch({
        o <- run_omni_sensor_score(active_gene_obj, te_obj = obj,
                                   group_col = condition_col)
        plot_omni_sensor_scores(o)
        plot_omni_sensor_irgs_corr(o)
        plot_omni_sensor_heatmap(o)
        status[["Sensors"]] <- "PASS"
        message("[Multi-omics] Sensor scoring complete.")
        o
      }, error = function(e) {
        status[["Sensors"]] <<- paste("FAIL:", e$message)
        message("[Multi-omics] Sensor scoring failed: ", e$message)
        active_gene_obj
      })

      # --- DIABLO ---
      cat("\n[Multi-omics] Running DIABLO multi-block sPLS-DA...\n")
      updated_gene_obj <- tryCatch({
        o <- run_omni_diablo(updated_gene_obj, te_obj = obj,
                             group_col = condition_col, ncomp = 2L)
        diablo_key   <- names(o$diablo)[length(o$diablo)]
        diablo_model <- o$diablo[[diablo_key]]$model
        message("[Multi-omics] DIABLO complete.")
        message("  Key: '", diablo_key, "'")
        message("  Visualise:")
        message("    model <- gene_obj$diablo[['", diablo_key, "']]$model")
        message("    mixOmics::plotDiablo(model, ncomp = 1)")
        message("    mixOmics::circosPlot(model, cutoff = 0.7)")
        message("    mixOmics::plotIndiv(model, legend = TRUE)")
        tryCatch({ mixOmics::plotDiablo(diablo_model, ncomp = 1) },
                 error = function(e) NULL)
        tryCatch({ mixOmics::plotIndiv(diablo_model, legend = TRUE) },
                 error = function(e) NULL)
        status[["DIABLO"]] <- "PASS"
        o
      }, error = function(e) {
        status[["DIABLO"]] <<- paste("FAIL:", e$message)
        message("[Multi-omics] DIABLO failed: ", e$message)
        updated_gene_obj
      })

      # Offer to save updated gene_obj
      cat("\n")
      save_gene <- utils::menu(
        c("Yes — save updated gene_obj to disk now",
          "No  — I will save it manually later"),
        title = "gene_obj has been updated with sensor scores and DIABLO results. Save now?"
      )
      if (save_gene == 1L) {
        save_path <- trimws(readline(
          prompt = "Save path [e.g. ~/results/gene/gene_obj_complete.rds]: "))
        tryCatch({
          saveRDS(updated_gene_obj, save_path)
          message("[Multi-omics] Saved updated gene_obj to: ", save_path)
        }, error = function(e) message("[Multi-omics] Save failed: ", e$message))
      }

      obj$.gene_obj_ref <- updated_gene_obj
    }
  }

  # ------------------------------------------------------------------
  # 10. EXPORT
  # ------------------------------------------------------------------
  decision <- .confirm("Export results (CSV tables + PNG plots)")
  if (decision == "stop") return(invisible(obj))
  if (decision == "run") {
    export_gene_obj <- obj$.gene_obj_ref %||% gene_obj
    tryCatch({
      export_omniseq(obj, outdir = outdir, gene_obj = export_gene_obj)
      status[["Export"]] <- "PASS"
    }, error = function(e) {
      status[["Export"]] <<- paste("FAIL:", e$message) })
  }

  # HTML report
  decision <- .confirm("Generate HTML report")
  if (decision != "stop" && decision == "run") {
    export_gene_obj <- obj$.gene_obj_ref %||% gene_obj
    tryCatch({
      export_omniseq_report(obj, outdir = outdir, gene_obj = export_gene_obj)
      status[["Report"]] <- "PASS"
    }, error = function(e) {
      status[["Report"]] <<- paste("FAIL:", e$message) })
  }

  # ------------------------------------------------------------------
  # FINAL SUMMARY
  # ------------------------------------------------------------------
  cat("\n")
  cat(strrep("=", 55), "\n")
  cat("  OmniSeq TE Pipeline — COMPLETE\n")
  cat(strrep("=", 55), "\n")
  res_df  <- as.data.frame(omni_results(obj))
  res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)
  n_sig   <- sum(res_df$padj < 0.05)
  n_up    <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1)
  te_info <- .omni_parse_te_ids(rownames(res_df))
  n_class <- length(unique(te_info$Class))
  cat("  Total loci tested    :", nrow(res_df), "\n")
  cat("  Significant (padj<0.05):", n_sig, "\n")
  cat("  Derepressed (LFC>1)  :", n_up, "\n")
  cat("  TE classes in dataset:", n_class, "\n\n")
  for (nm in names(status)) {
    icon <- if (startsWith(status[[nm]], "PASS")) "\u2713" else "\u2717"
    cat("  ", icon, " ", nm, ": ", status[[nm]], "\n", sep = "")
  }
  cat(strrep("=", 55), "\n\n")

  print(obj)
  return(invisible(obj))
}
