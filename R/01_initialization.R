# =============================================================================
# OmniSeq — 01_initialization.R
#
# Fixes vs FlexiSeq:
#   - data_type ("gene" vs "TE") detected ONCE here, stored in params
#   - params / model_info / results_cache are separate slots (not one grab-bag)
#   - interactive() guard on ALL menu()/readline() calls
#   - covariates input accepts both programmatic vector and interactive menu
#   - species + org_db set at init and propagated — not re-asked downstream
#   - exclude_genes stored in params, not silently dropped after filtering
#   - returns a validated OmniSeq object
# =============================================================================


#' Initialize an OmniSeq Object
#'
#' Takes a raw count matrix and sample metadata, aligns them, detects whether
#' the data are standard genes or transposable elements (TEs), sets up the
#' statistical design, and returns a validated OmniSeq object ready for QC
#' and differential expression.
#'
#' @param count_matrix A matrix of **raw integer counts** (genes/TEs in rows,
#'   samples in columns). Row names must be feature identifiers; column names
#'   must match `rownames(metadata)`.
#' @param metadata A `data.frame` of sample metadata. Row names must match
#'   `colnames(count_matrix)`.
#' @param condition_col Character. Name of the primary biological condition
#'   column in `metadata`. If `NULL` and the session is interactive, the user
#'   is prompted to choose.
#' @param covariates Character vector of covariate column names to include in
#'   the model formula. Defaults to none. In interactive mode, the user is
#'   prompted if this is `NULL`.
#' @param exclude_genes Optional character vector of gene/feature identifiers
#'   to permanently remove from the count matrix **before** any analysis.
#'   Useful for masking immunoglobulin genes, spike-ins, etc.
#' @param filter_coding Logical. If `TRUE`, keep only features whose names
#'   look like standard gene symbols (drops lncRNA/novel IDs that cannot be
#'   mapped to GO terms). TE locus IDs are always kept. Default: `FALSE`.
#' @param species Character. Species name passed to `msigdbr` and used to
#'   select the correct annotation database. Default: `"Homo sapiens"`.
#' @param org_db Character. Bioconductor annotation package for GO mapping.
#'   Default: `"org.Hs.eg.db"`. For mouse, use `"org.Mm.eg.db"`.
#' @param ref_level Optional character. The reference (baseline) level for the
#'   condition factor. If `NULL`, R's default alphabetical ordering is used,
#'   unless a level named `"Control"` or `"control"` is detected automatically.
#' @param interactive Logical. If `TRUE` (and `R` session is interactive),
#'   prompts the user with menus to select the design. Ignored in
#'   non-interactive sessions. Default: `FALSE`.
#'
#' @return A validated `OmniSeq` object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Programmatic (fast-track) mode
#' obj <- init_omniseq(
#'   count_matrix  = my_counts,
#'   metadata      = my_meta,
#'   condition_col = "Genotype",
#'   covariates    = "Batch",
#'   species       = "Homo sapiens",
#'   org_db        = "org.Hs.eg.db"
#' )
#'
#' # Interactive mode (runs menu prompts in the console)
#' obj <- init_omniseq(my_counts, my_meta, interactive = TRUE)
#' }
init_omniseq <- function(count_matrix,
                         metadata,
                         condition_col  = NULL,
                         covariates     = NULL,
                         exclude_genes  = NULL,
                         filter_coding  = FALSE,
                         species        = "Homo sapiens",
                         org_db         = "org.Hs.eg.db",
                         ref_level      = NULL,
                         interactive    = FALSE) {

  # -------------------------------------------------------------------
  # 0. PREREQUISITE CHECK
  # -------------------------------------------------------------------
  .require_pkg("SummarizedExperiment")

  # -------------------------------------------------------------------
  # 1. INPUT VALIDATION
  # -------------------------------------------------------------------
  if (!is.matrix(count_matrix))
    count_matrix <- as.matrix(count_matrix)

  if (is.null(rownames(count_matrix)))
    stop("count_matrix must have row names (feature identifiers).", call. = FALSE)

  if (is.null(colnames(count_matrix)))
    stop("count_matrix must have column names (sample identifiers).", call. = FALSE)

  if (!is.data.frame(metadata))
    metadata <- as.data.frame(metadata)

  if (is.null(rownames(metadata)))
    stop("metadata must have row names matching colnames(count_matrix).", call. = FALSE)

  # Check for sample alignment
  missing_in_counts <- setdiff(rownames(metadata), colnames(count_matrix))
  if (length(missing_in_counts) > 0)
    stop(
      "The following metadata row names are not in count_matrix columns: ",
      paste(head(missing_in_counts, 5), collapse = ", "),
      if (length(missing_in_counts) > 5) paste0(" ... (", length(missing_in_counts) - 5, " more)"),
      call. = FALSE
    )

  # Align column order to metadata row order
  count_matrix <- count_matrix[, rownames(metadata), drop = FALSE]

  # -------------------------------------------------------------------
  # 2. OPTIONAL GENE EXCLUSION (pre-DE)
  # -------------------------------------------------------------------
  excluded_names <- character(0)
  if (!is.null(exclude_genes) && length(exclude_genes) > 0) {
    hits <- intersect(exclude_genes, rownames(count_matrix))
    if (length(hits) > 0) {
      count_matrix  <- count_matrix[!rownames(count_matrix) %in% hits, , drop = FALSE]
      excluded_names <- hits
      message("[init_omniseq] Permanently removed ", length(hits),
              " feature(s) from count matrix before analysis.")
    } else {
      message("[init_omniseq] Note: none of the genes listed for exclusion ",
              "were found in the count matrix.")
    }
  }

  # -------------------------------------------------------------------
  # 2B. OPTIONAL PROTEIN-CODING FILTER
  #     When filter_coding = TRUE, removes non-coding features whose names
  #     are not standard HGNC gene symbols. This is recommended when the
  #     count matrix was generated from a full GENCODE annotation that
  #     includes lncRNAs, miRNAs, novel transcripts and pseudogenes, as
  #     these features cannot be mapped to GO terms and dilute enrichment.
  #     All TE locus IDs are automatically preserved regardless.
  # -------------------------------------------------------------------
  if (filter_coding) {
    data_type_pre <- .detect_data_type(count_matrix)
    if (data_type_pre == "gene") {
      n_before  <- nrow(count_matrix)
      keep_mask <- .filter_coding_genes(rownames(count_matrix))
      count_matrix <- count_matrix[keep_mask, , drop = FALSE]
      n_removed <- n_before - nrow(count_matrix)
      message("[init_omniseq] filter_coding: removed ", n_removed,
              " non-coding/unannotated features, kept ",
              nrow(count_matrix), " gene symbols.")
    } else {
      message("[init_omniseq] filter_coding ignored for TE data.")
    }
  }

  # -------------------------------------------------------------------
  # 3. AUTO-DETECT DATA TYPE (gene vs TE)
  #    Stored in params once here — never re-detected downstream.
  # -------------------------------------------------------------------
  data_type <- .detect_data_type(count_matrix)
  message("[init_omniseq] Data type detected: ", data_type)

  # -------------------------------------------------------------------
  # 4. DESIGN SETUP
  #    Interactive path only fires when the user explicitly opts in AND
  #    the R session is actually interactive.
  # -------------------------------------------------------------------
  cols <- colnames(metadata)

  if (interactive && base::interactive()) {
    result      <- .interactive_design_setup(cols)
    condition_col <- result$condition
    covariates    <- result$covariates
  } else {
    # Programmatic / fast-track path
    if (is.null(condition_col))
      stop(
        "condition_col must be provided in non-interactive mode. ",
        "Set interactive = TRUE to use the menu-driven setup.",
        call. = FALSE
      )
    if (!(condition_col %in% cols))
      stop("condition_col '", condition_col, "' not found in metadata columns. ",
           "Available columns: ", paste(cols, collapse = ", "), call. = FALSE)

    # Validate any covariates provided
    covariates <- covariates %||% character(0)
    bad_covs   <- setdiff(covariates, cols)
    if (length(bad_covs) > 0)
      stop("The following covariates were not found in metadata: ",
           paste(bad_covs, collapse = ", "), call. = FALSE)

    # Warn if a covariate is the same as the condition
    if (condition_col %in% covariates)
      warning("condition_col '", condition_col,
              "' is also listed as a covariate. It will be removed from covariates.",
              call. = FALSE)
    covariates <- setdiff(covariates, condition_col)
  }

  # -------------------------------------------------------------------
  # 5. FACTOR SETUP & REFERENCE LEVEL
  # -------------------------------------------------------------------
  metadata[[condition_col]] <- as.factor(metadata[[condition_col]])

  if (!is.null(ref_level)) {
    if (!ref_level %in% levels(metadata[[condition_col]]))
      stop("ref_level '", ref_level, "' is not a level of condition '",
           condition_col, "'. Available levels: ",
           paste(levels(metadata[[condition_col]]), collapse = ", "),
           call. = FALSE)
    metadata[[condition_col]] <- stats::relevel(metadata[[condition_col]], ref = ref_level)
    message("[init_omniseq] Reference level set to: '", ref_level, "'")

  } else {
    # Auto-detect common control labels and set them as reference
    lvls        <- levels(metadata[[condition_col]])
    control_hit <- lvls[tolower(lvls) %in% c("control", "ctrl", "wt",
                                             "wildtype", "untreated", "vehicle")]
    if (length(control_hit) == 1) {
      metadata[[condition_col]] <- stats::relevel(metadata[[condition_col]],
                                                  ref = control_hit)
      message("[init_omniseq] Auto-set reference level to: '", control_hit, "'")
    }
  }

  # Coerce covariates that are character to factor for lm compatibility
  for (cv in covariates) {
    if (is.character(metadata[[cv]])) {
      metadata[[cv]] <- as.factor(metadata[[cv]])
    }
  }

  # -------------------------------------------------------------------
  # 6. BUILD THE SUMMARIZEDEXPERIMENT
  # -------------------------------------------------------------------
  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = count_matrix),
    colData = metadata
  )

  # -------------------------------------------------------------------
  # 7. CONSTRUCT, VALIDATE & RETURN
  # -------------------------------------------------------------------
  obj <- new_omniseq(
    data           = se,
    data_type      = data_type,
    condition      = condition_col,
    covariates     = covariates,
    species        = species,
    org_db         = org_db,
    excluded_genes = excluded_names
  )

  validate_omniseq(obj)

  # Print a clean summary to the console
  message("\n[OmniSeq] Object created successfully.")
  message("  Features : ", nrow(se))
  message("  Samples  : ", ncol(se))
  message("  Data type: ", data_type)
  message("  Formula  : ", deparse(obj$params$formula))
  if (length(excluded_names) > 0)
    message("  Excluded : ", length(excluded_names), " feature(s) removed pre-analysis")

  return(obj)
}


# =============================================================================
# INTERNAL: Interactive design setup (only called when interactive = TRUE)
# =============================================================================

#' @keywords internal
.interactive_design_setup <- function(cols) {

  cat("\n", strrep("=", 50), "\n", sep = "")
  cat("  OmniSeq — Interactive Design Setup\n")
  cat(strrep("=", 50), "\n\n", sep = "")

  # --- Step 1: Primary condition ---
  cat("[STEP 1] Select the PRIMARY biological condition column:\n\n")
  cond_idx <- utils::menu(cols, title = "Primary condition:")
  if (cond_idx == 0L) stop("Setup aborted by user.", call. = FALSE)
  condition <- cols[[cond_idx]]
  cat("\n  -> Selected condition: '", condition, "'\n\n", sep = "")

  # --- Step 2: Covariates ---
  remaining <- setdiff(cols, condition)
  covariates <- character(0)

  if (length(remaining) > 0) {
    cat("[STEP 2] Select COVARIATE columns to include in the model.\n")
    cat("         Press 0 when you are done (or straight away for no covariates).\n\n")

    repeat {
      available <- setdiff(remaining, covariates)
      if (length(available) == 0) break

      label <- if (length(covariates) == 0) "No covariates selected yet"
      else paste("Selected:", paste(covariates, collapse = ", "))

      cov_idx <- utils::menu(available, title = paste0("Add covariate (", label, "):"))
      if (cov_idx == 0L) break
      covariates <- c(covariates, available[[cov_idx]])
    }
  }

  if (length(covariates) > 0)
    cat("\n  -> Covariates: ", paste(covariates, collapse = ", "), "\n\n", sep = "")
  else
    cat("\n  -> No covariates selected.\n\n")

  list(condition = condition, covariates = covariates)
}
