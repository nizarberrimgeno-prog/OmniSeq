# =============================================================================
# OmniSeq — 00_object.R
# Core S3 class: definition, constructor, validator, accessors, print/summary
# =============================================================================


# -----------------------------------------------------------------------------
# 1. INTERNAL CONSTRUCTOR
#    Called only by init_omniseq(). Users never call this directly.
# -----------------------------------------------------------------------------

#' @keywords internal
new_omniseq <- function(data,
                        data_type  = c("gene", "TE"),
                        condition  = character(1),
                        covariates = character(0),
                        species    = "Homo sapiens",
                        org_db     = "org.Hs.eg.db",
                        excluded_genes = character(0)) {

  data_type <- match.arg(data_type)

  # Build the three cleanly-separated slots
  params <- list(
    condition      = condition,
    covariates     = covariates,
    formula        = .build_formula(covariates, condition),
    data_type      = data_type,          # "gene" or "TE" — set ONCE here
    species        = species,
    org_db         = org_db,
    excluded_genes = excluded_genes
  )

  model_info <- list(
    dds            = NULL,   # The fitted DESeqDataSet
    vst            = NULL,   # Cached VST — computed ONCE, reused everywhere
    extracted_coef = NULL,   # Name of the active contrast
    lfc_shrinkage  = "none"  # "apeglm" | "ashr" | "none"
  )

  results_cache <- list(
    # Keyed by contrast name so multiple contrasts can coexist
    # e.g. results_cache[["Treatment_vs_Control"]] <- DESeqResults object
  )

  structure(
    list(
      data          = data,          # SummarizedExperiment
      params        = params,
      model_info    = model_info,
      results_cache = results_cache,
      active_contrast = NA_character_, # Key into results_cache — never NULL (NULL removes the slot)
      qc            = NULL,
      enrichment    = NULL,
      advanced_gsea = NULL,
      pathway_results = NULL
    ),
    class = "OmniSeq"
  )
}


# -----------------------------------------------------------------------------
# 2. VALIDATOR
#    Called by new_omniseq() and can be called by the user to sanity-check.
# -----------------------------------------------------------------------------

#' Validate an OmniSeq object
#'
#' Checks structural integrity of an OmniSeq object and returns it invisibly
#' if valid, or stops with a descriptive error.
#'
#' @param obj An OmniSeq object.
#' @return The object, invisibly.
#' @export
validate_omniseq <- function(obj) {

  if (!inherits(obj, "OmniSeq"))
    stop("Object is not of class 'OmniSeq'.", call. = FALSE)

  required_slots <- c("data", "params", "model_info", "results_cache",
                      "active_contrast", "qc")
  missing_slots <- setdiff(required_slots, names(obj))
  if (length(missing_slots) > 0)
    stop("OmniSeq object is malformed. Missing slots: ",
         paste(missing_slots, collapse = ", "), call. = FALSE)

  required_params <- c("condition", "covariates", "formula",
                       "data_type", "species", "org_db")
  missing_params <- setdiff(required_params, names(obj$params))
  if (length(missing_params) > 0)
    stop("OmniSeq $params slot is malformed. Missing fields: ",
         paste(missing_params, collapse = ", "), call. = FALSE)

  if (!inherits(obj$data, "SummarizedExperiment"))
    stop("$data slot must be a SummarizedExperiment object.", call. = FALSE)

  if (!obj$params$data_type %in% c("gene", "TE"))
    stop("$params$data_type must be 'gene' or 'TE'.", call. = FALSE)

  invisible(obj)
}


# -----------------------------------------------------------------------------
# 3. ACCESSOR FUNCTIONS
#    Users and internal functions should use these, not `obj$...` directly.
#    This insulates downstream code from slot renames.
# -----------------------------------------------------------------------------

#' Get the active DE results from an OmniSeq object
#'
#' Returns the DESeqResults for the currently active contrast.
#'
#' @param obj An OmniSeq object.
#' @param contrast Optional contrast name. Defaults to the active contrast.
#' @return A DESeqResults object, or NULL if not yet run.
#' @export
omni_results <- function(obj, contrast = NULL) {
  .check_class(obj)
  key <- contrast %||% obj$active_contrast
  if (is.null(key)) return(NULL)
  obj$results_cache[[key]]$results
}

#' Get the fitted DESeq2 model from an OmniSeq object
#'
#' @param obj An OmniSeq object.
#' @return A DESeqDataSet, or NULL if the model has not been fitted.
#' @export
omni_model <- function(obj) {
  .check_class(obj)
  obj$model_info$dds
}

#' Get the cached VST matrix from an OmniSeq object
#'
#' @param obj An OmniSeq object.
#' @return A DESeqTransform (VST) object, or NULL.
#' @export
omni_vst <- function(obj) {
  .check_class(obj)
  obj$model_info$vst
}

#' Get all available contrast names
#'
#' @param obj An OmniSeq object.
#' @return A character vector of contrast names stored in the results cache.
#' @export
omni_contrasts <- function(obj) {
  .check_class(obj)
  names(obj$results_cache)
}

#' Set the active contrast
#'
#' Switches which cached contrast downstream functions (plots, enrichment,
#' TF activity, networks) operate on. No results are recomputed.
#'
#' @param obj An OmniSeq object.
#' @param contrast Character. A contrast name from [omni_contrasts()].
#' @return The OmniSeq object with `active_contrast` updated.
#' @export
omni_set_contrast <- function(obj, contrast) {
  .check_class(obj)
  if (!is.character(contrast) || length(contrast) != 1 ||
      !contrast %in% names(obj$results_cache))
    stop("Contrast '", contrast, "' not found. Available: ",
         paste(names(obj$results_cache), collapse = ", "), call. = FALSE)
  obj$active_contrast <- contrast
  obj
}

#' Get the QC metrics slot
#'
#' @param obj An OmniSeq object.
#' @return A list of QC metrics, or NULL.
#' @export
omni_qc <- function(obj) {
  .check_class(obj)
  obj$qc
}

#' Get sample metadata as a data frame
#'
#' @param obj An OmniSeq object.
#' @return A data.frame of sample metadata.
#' @export
omni_metadata <- function(obj) {
  .check_class(obj)
  as.data.frame(SummarizedExperiment::colData(obj$data))
}

#' Get the raw count matrix
#'
#' @param obj An OmniSeq object.
#' @return A matrix of raw integer counts.
#' @export
omni_counts <- function(obj) {
  .check_class(obj)
  SummarizedExperiment::assay(obj$data)
}


# -----------------------------------------------------------------------------
# 4. S3 METHODS: print, summary, [
# -----------------------------------------------------------------------------

#' Print method for OmniSeq objects
#'
#' @param x An OmniSeq object.
#' @param ... Ignored.
#' @export
print.OmniSeq <- function(x, ...) {
  validate_omniseq(x)

  n_samples   <- ncol(x$data)
  n_features  <- nrow(x$data)
  has_model   <- !is.null(x$model_info$dds)
  has_vst     <- !is.null(x$model_info$vst)
  n_contrasts <- length(x$results_cache)

  cat("-- OmniSeq object ", packageVersion_safe(), " --\n", sep = "")
  cat("  Features : ", n_features, "\n", sep = "")
  cat("  Samples  : ", n_samples,  "\n", sep = "")
  cat("  Data type: ", x$params$data_type, "\n", sep = "")
  cat("  Species  : ", x$params$species, "\n", sep = "")
  cat("\n")
  cat("  Design formula : ", deparse(x$params$formula), "\n", sep = "")
  cat("  Condition      : ", x$params$condition, "\n", sep = "")

  if (length(x$params$covariates) > 0)
    cat("  Covariates     : ",
        paste(x$params$covariates, collapse = ", "), "\n", sep = "")

  cat("\n")
  cat("  Model fitted : ", if (has_model) "yes" else "no", "\n", sep = "")
  cat("  VST cached   : ", if (has_vst)   "yes" else "no", "\n", sep = "")
  cat("  Contrasts run: ", n_contrasts, "\n", sep = "")

  if (n_contrasts > 0) {
    active <- x$active_contrast
    for (nm in names(x$results_cache)) {
      entry   <- x$results_cache[[nm]]
      is_act  <- !is.null(active) && nm == active
      is_mask <- grepl("\\(excl\\.", nm)
      tag     <- if (is_act) " <- active" else ""
      prefix  <- if (is_mask) "    ~ " else "    - "
      cat(prefix, nm, tag, "\n", sep = "")
    }
  }

  # Advanced results summary
  has_ora    <- !is.null(x$enrichment)
  has_gsea   <- !is.null(x$advanced_gsea)
  has_ssgsea <- !is.null(x$advanced_ssgsea)
  has_tf     <- !is.null(x$tf_results)
  has_deconv <- !is.null(x$deconv)
  has_string <- !is.null(x$string_network)
  has_te     <- !is.null(x$te_family)

  advanced <- c(
    if (has_ora)    "ORA",
    if (has_gsea)   "GSEA",
    if (has_ssgsea) "ssGSEA",
    if (has_tf)     "TF activity",
    if (has_deconv) paste0("Deconv (", x$params$deconv_method %||% "?", ")"),
    if (has_string) "STRING network",
    if (has_te)     "TE family"
  )

  if (length(advanced) > 0) {
    cat("\n  Results available: ",
        paste(advanced, collapse = ", "), "\n", sep = "")
  }

  invisible(x)
}

#' Summary method for OmniSeq objects
#'
#' @param object An OmniSeq object.
#' @param ... Ignored.
#' @export
summary.OmniSeq <- function(object, ...) {
  print(object)

  res <- omni_results(object)
  if (!is.null(res)) {
    res_df <- as.data.frame(res)
    sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
    up   <- sum(sig$log2FoldChange > 1,  na.rm = TRUE)
    down <- sum(sig$log2FoldChange < -1, na.rm = TRUE)
    cat("\n  Active contrast results (padj < 0.05, |LFC| > 1):\n")
    cat("    Up-regulated  : ", up,   "\n", sep = "")
    cat("    Down-regulated: ", down, "\n", sep = "")
  }

  if (!is.null(object$qc)) {
    stats <- object$qc$stats
    cat("\n  QC summary:\n")
    cat("    Median library size  : ",
        format(median(stats$Library_Size), big.mark = ","), "\n", sep = "")
    cat("    Median detected genes: ",
        format(median(stats$Detected_Genes), big.mark = ","), "\n", sep = "")
  }

  invisible(object)
}


# -----------------------------------------------------------------------------
# 5. INTERNAL HELPERS
#    Prefixed with . so they are clearly private.
# -----------------------------------------------------------------------------

#' @keywords internal
.check_class <- function(obj) {
  if (!inherits(obj, "OmniSeq"))
    stop("Expected an OmniSeq object.", call. = FALSE)
}

#' @keywords internal
.build_formula <- function(covariates, condition) {
  terms <- c(covariates, condition)
  as.formula(paste("~", paste(terms, collapse = " + ")))
}

#' @keywords internal
.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required but not installed. ",
         "Install it with: ", .omni_install_hint(pkg),
         call. = FALSE)
}

#' @keywords internal
.require_slot <- function(obj, slot_name, run_fn) {
  if (is.null(obj[[slot_name]]))
    stop("Slot '", slot_name, "' is empty. Run ", run_fn, "() first.",
         call. = FALSE)
}

#' @keywords internal
`%||%` <- function(x, y) if (!is.null(x)) x else y

#' Compute a VST matrix, falling back below DESeq2::vst()'s row minimum
#'
#' DESeq2::vst() requires at least `nsub` (default 1000) rows for its fast
#' dispersion-subsampling path and errors below that — a common situation
#' with filtered gene panels or TE locus counts. Falls back to
#' DESeq2::varianceStabilizingTransformation(), which has no row minimum.
#'
#' @param dds A DESeqDataSet.
#' @param blind Logical, passed through to the DESeq2 transform.
#' @return A DESeqTransform (VST) object.
#' @keywords internal
.omni_vst <- function(dds, blind) {
  if (nrow(dds) < 1000) {
    DESeq2::varianceStabilizingTransformation(dds, blind = blind)
  } else {
    DESeq2::vst(dds, blind = blind)
  }
}

#' @keywords internal
packageVersion_safe <- function() {
  tryCatch(
    paste0("v", utils::packageVersion("OmniSeq")),
    error = function(e) "(development)"
  )
}

# TE class vocabulary shared by .detect_data_type() and .filter_coding_genes()
# — keep this as the single source of truth so the two can't drift apart.
.OMNISEQ_TE_CLASSES <- c("SINE", "LINE", "ERV", "Retroposon", "DNA", "RC",
                         "Satellite", "rRNA", "scRNA", "snRNA", "tRNA", "Unknown")

#' @keywords internal
.detect_data_type <- function(count_matrix) {
  # Strict TE locus ID format: "<number>:<TEClass>:<Family>"
  # e.g. "1000:SINE:MIRb", "1000014:ERV:MSTB"
  # Simple colon check is too broad — gene annotation strings also contain colons.
  sample_rows <- head(rownames(count_matrix), 100)
  te_pattern  <- paste0("^[0-9]+:(", paste(.OMNISEQ_TE_CLASSES, collapse = "|"), "):")
  n_match     <- sum(grepl(te_pattern, sample_rows))
  if (n_match / length(sample_rows) >= 0.8) "TE" else "gene"
}

#' Clean gene symbol strings for enrichment mapping
#'
#' Filters out gene names that are not standard HGNC symbols — removes
#' Ensembl IDs, long annotation strings, pseudogene descriptions, and
#' other names that will fail bitr() Entrez ID mapping.
#'
#' @param genes Character vector of gene names.
#' @return Character vector with non-symbol names replaced by NA.
#' @keywords internal
.clean_gene_symbols <- function(genes) {
  # Keep only strings that look like HGNC gene symbols:
  # - 1-25 characters
  # - Only letters, digits, hyphens, dots, underscores
  # - Does NOT start with ENSG (Ensembl)
  # - Does NOT contain brackets, parentheses, or Source: strings
  is_clean <- nchar(genes) <= 25 &
    grepl("^[A-Za-z][A-Za-z0-9\\-\\./_]*$", genes) &
    !grepl("^ENSG|^ENSM|\\[|\\(|Source:|novel|pseudogene",
           genes, ignore.case = TRUE)
  genes[!is_clean] <- NA_character_
  genes
}

#' Filter count matrix row names to protein-coding gene symbols
#'
#' Returns a logical mask keeping only features whose names look like
#' standard HGNC gene symbols. Removes lncRNAs, miRNAs, novel transcripts,
#' pseudogenes, Ensembl IDs, and accession-style names that cannot be
#' mapped to GO terms. TE locus IDs (number:Class:Family format) are
#' always preserved.
#'
#' @param gene_names Character vector of feature names.
#' @return Logical vector of length equal to gene_names.
#' @keywords internal
.filter_coding_genes <- function(gene_names) {
  # Always keep TE locus IDs
  te_pattern <- paste0("^[0-9]+(:(", paste(.OMNISEQ_TE_CLASSES, collapse = "|"), ")):")
  is_te      <- grepl(te_pattern, gene_names)

  # Keep features that:
  # - are 1-25 characters
  # - start with a letter
  # - contain only letters, digits, hyphens, dots, underscores
  # - do NOT match known non-coding patterns
  is_clean <- nchar(gene_names) <= 25 &
    grepl("^[A-Za-z][A-Za-z0-9\\-\\._]*$", gene_names) &
    !grepl(paste0(
      "^ENSG|^ENSM|",                      # Ensembl IDs
      "novel|pseudogene|\\[Source|",        # annotation strings
      "Y_RNA|snoRNA|snRNA|miRNA|",          # small ncRNAs
      "MIR[0-9]|RNU[0-9]|SNOR|",           # miRNA/snoRNA symbols
      "LINC[0-9]|",                         # lncRNA accessions
      "^AL[0-9]{5}|^AC[0-9]{5}|",          # lncRNA accessions
      "^AP[0-9]{5}|^FO[0-9]{5}|",          # lncRNA accessions
      "^BX[0-9]{5}|^CR[0-9]{5}|",          # lncRNA accessions
      "^Z[0-9]{5}|^AJ[0-9]{5}"             # lncRNA accessions
    ), gene_names, ignore.case = TRUE)

  is_te | is_clean
}
