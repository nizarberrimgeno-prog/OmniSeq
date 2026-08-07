# =============================================================================
# OmniSeq — 07_TE_Analysis.R
#
# Transposable Element (TE) specific analysis functions.
# These functions activate automatically when obj$params$data_type == "TE"
# but can also be called directly on any OmniSeq object with TE locus IDs.
#
# TE locus ID format: "1000014:ERV:MSTB"
#   - Field 1: Numeric locus ID assigned by HERVd / featureCounts
#   - Field 2: TE Class  (e.g. ERV, LINE, SINE, Retroposon)
#   - Field 3: TE Family (e.g. MSTB, L1PA6, AluSx, SVA_A)
#
# Genomic coordinates are NOT embedded in the locus ID.
# Use omni_build_te_coords(hervd_dir) to build a coordinate lookup table
# from the HERVd GTF files, then pass it to run_omni_te_proximity().
#
# ANALYSIS:
#   run_omni_te_family()        Aggregate DE results by TE family/class
#   run_omni_te_proximity()     Prepare BED coordinates for GREAT analysis
#   run_omni_te_great()         Run rGREAT proximity enrichment
#
# PLOTTING:
#   plot_omni_te_family()       TE family enrichment barplot
#   plot_omni_te_heatmap()      TE family activity heatmap (VST aggregated)
#   plot_omni_te_class_pie()    TE class composition pie/donut chart
#   plot_omni_te_volcano()      Volcano colored by TE class
#   plot_omni_te_lfc_density()  LFC density curves per TE class
#   plot_omni_te_great()        GREAT proximity enrichment dotplot
#
# MULTI-OMICS:
#   plot_omni_te_gene_corr()    Correlate TE load vs gene signature score
#
# Design rules:
#   - All parsers use the colon-delimited locus ID format
#   - TE results stored in obj$te_family (aggregated stats)
#   - GREAT results stored in obj$te_great
#   - All functions work on any OmniSeq object — data_type check is a warning
#     not a hard stop, so users can run these on gene data if needed
# =============================================================================


# =============================================================================
# INTERNAL: TE ID PARSER
# =============================================================================

#' Parse TE locus IDs into a structured data frame
#' @keywords internal
.omni_parse_te_ids <- function(ids) {
  # Expected format: "chr1-10000-20000:ERV:MSTB"
  # Returns data.frame with columns: Locus_ID, Coords, Class, Family
  parts <- strsplit(ids, ":", fixed = TRUE)

  data.frame(
    Locus_ID = ids,
    Coords   = sapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_),
    Class    = sapply(parts, function(x) if (length(x) >= 2) x[2] else "Unknown"),
    Family   = sapply(parts, function(x) if (length(x) >= 3) x[3] else "Unknown"),
    stringsAsFactors = FALSE
  )
}

#' Warn if object does not appear to be TE data
#' @keywords internal
.omni_check_te <- function(obj) {
  if (!is.null(obj$params$data_type) && obj$params$data_type != "TE") {
    message("[TE] Note: obj$params$data_type is '", obj$params$data_type,
            "'. TE functions expect locus IDs in 'coords:Class:Family' format.")
  }
}


# =============================================================================
# 1. AGGREGATE DE RESULTS BY TE FAMILY / CLASS
# =============================================================================

#' Aggregate DE Results by TE Family and Class
#'
#' Parses the locus IDs in the active contrast results, aggregates by TE
#' Family and Class, and returns summary statistics (count, mean LFC,
#' fraction significant). Results are stored in `obj$te_family`.
#'
#' @param obj A validated OmniSeq object with TE DE results.
#' @param lfc_thresh Log2 fold change threshold. Default: `1`.
#' @param pval_thresh Adjusted p-value threshold. Default: `0.05`.
#' @param direction Character. `"up"`, `"down"`, or `"both"`. Default: `"both"`.
#'
#' @return The OmniSeq object with `$te_family` populated.
#' @export
run_omni_te_family <- function(obj,
                               lfc_thresh  = 1,
                               pval_thresh = 0.05,
                               direction   = c("both", "up", "down")) {

  validate_omniseq(obj)
  .omni_check_te(obj)
  direction <- match.arg(direction)

  res <- omni_results(obj)
  if (is.null(res))
    stop("No DE results. Run run_omni_de() first.", call. = FALSE)

  res_df      <- as.data.frame(res)
  res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)

  # Parse TE IDs
  te_df <- .omni_parse_te_ids(rownames(res_df))
  te_df <- cbind(te_df, res_df[, c("log2FoldChange", "padj",
                                   "baseMean"), drop = FALSE])

  # Filter by significance and direction
  sig_mask <- te_df$padj < pval_thresh

  if (direction == "up")
    sig_mask <- sig_mask & te_df$log2FoldChange > lfc_thresh
  else if (direction == "down")
    sig_mask <- sig_mask & te_df$log2FoldChange < -lfc_thresh
  else
    sig_mask <- sig_mask & abs(te_df$log2FoldChange) > lfc_thresh

  te_df$Significant <- sig_mask

  # Aggregate by Family
  fam_stats <- do.call(rbind, lapply(
    split(te_df, te_df$Family), function(df) {
      data.frame(
        Family       = df$Family[1],
        Class        = df$Class[1],
        N_Loci       = nrow(df),
        N_Sig        = sum(df$Significant, na.rm = TRUE),
        Frac_Sig     = round(mean(df$Significant, na.rm = TRUE), 3),
        Mean_LFC     = round(mean(df$log2FoldChange, na.rm = TRUE), 3),
        Mean_baseMean = round(mean(df$baseMean, na.rm = TRUE), 1),
        stringsAsFactors = FALSE
      )
    }
  ))

  fam_stats <- fam_stats[order(-fam_stats$N_Sig), ]
  rownames(fam_stats) <- NULL

  # Also aggregate by Class
  class_stats <- do.call(rbind, lapply(
    split(te_df, te_df$Class), function(df) {
      data.frame(
        Class    = df$Class[1],
        N_Loci   = nrow(df),
        N_Sig    = sum(df$Significant, na.rm = TRUE),
        Frac_Sig = round(mean(df$Significant, na.rm = TRUE), 3),
        Mean_LFC = round(mean(df$log2FoldChange, na.rm = TRUE), 3),
        stringsAsFactors = FALSE
      )
    }
  ))
  class_stats <- class_stats[order(-class_stats$N_Sig), ]
  rownames(class_stats) <- NULL

  obj$te_family <- list(
    by_family    = fam_stats,
    by_class     = class_stats,
    all_loci     = te_df,
    lfc_thresh   = lfc_thresh,
    pval_thresh  = pval_thresh,
    direction    = direction,
    contrast     = obj$active_contrast
  )

  n_sig_fam <- sum(fam_stats$N_Sig > 0)
  message("[run_omni_te_family] Parsed ", nrow(te_df), " TE loci across ",
          nrow(fam_stats), " families and ", nrow(class_stats), " classes.")
  message("  Families with >= 1 significant locus: ", n_sig_fam)
  message("  Access with obj$te_family")

  return(obj)
}


# =============================================================================
# 2. BUILD TE COORDINATE TABLE FROM USER-PROVIDED GTF FILES
# =============================================================================

#' Build a TE Coordinate Lookup Table from GTF Files
#'
#' Loads one or more TE annotation GTF files and builds a coordinate table
#' that maps every locus ID in your count matrices to its exact genomic
#' location. Each locus is matched by the `transcript_id` field in the GTF.
#'
#' This function is fully generic — it works with any TE annotation GTF
#' that uses `transcript_id` as the locus identifier (e.g. HERVd, Telescope,
#' TEtranscripts, RepeatMasker-derived GTFs). The only requirement is that
#' the `transcript_id` values match the row names of your count matrix.
#'
#' **Run once and save.** This function can take several minutes for large
#' GTF files. Save the result with `save_path` and reload with `readRDS()`
#' in future sessions.
#'
#' @param gtf_files Character vector. One or more paths to TE annotation
#'   GTF files. Each file is loaded, collapsed to locus-level coordinates,
#'   and combined into a single table.
#' @param locus_id_col Character. Name of the GTF attribute column that
#'   contains the locus identifier matching your count matrix row names.
#'   Default: `"transcript_id"`. Change to `"gene_id"` if needed.
#' @param save_path Character. Optional path to save the result as an RDS
#'   file for fast reloading in future sessions. Default: `NULL`.
#' @param genome Character. Genome build for filtering to standard
#'   chromosomes. `"human"` (default) uses chr1-22, X, Y.
#'   `"mouse"` uses chr1-19, X, Y. `"all"` keeps all chromosomes.
#'
#' @return A data.frame with columns: Locus_ID, Chr, Start, End, Strand,
#'   Class, Family. Class and Family are parsed from the locus ID if it
#'   follows the `number:Class:Family` format, otherwise set to NA.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Single GTF file
#' te_coords <- omni_build_te_coords(
#'   gtf_files = "/path/to/TE_annotation.gtf",
#'   save_path = "~/te_coords.rds"
#' )
#'
#' # Multiple GTF files (one per TE class)
#' te_coords <- omni_build_te_coords(
#'   gtf_files = c(
#'     "/path/to/sine.gtf",
#'     "/path/to/line.gtf",
#'     "/path/to/erv.gtf",
#'     "/path/to/retroposon.gtf"
#'   ),
#'   save_path = "~/te_coords.rds"
#' )
#'
#' # Reload fast in future sessions
#' te_coords <- readRDS("~/te_coords.rds")
#'
#' # Use in proximity analysis
#' obj <- run_omni_te_proximity(obj, te_coords = te_coords)
#' }
omni_build_te_coords <- function(gtf_files,
                                 locus_id_col = "transcript_id",
                                 save_path    = NULL,
                                 genome       = c("human", "mouse", "all")) {

  .require_pkg("rtracklayer")
  .require_pkg("dplyr")

  genome <- match.arg(genome)

  # Standard chromosome sets
  valid_chrs <- switch(genome,
                       "human" = paste0("chr", c(1:22, "X", "Y")),
                       "mouse" = paste0("chr", c(1:19, "X", "Y")),
                       "all"   = NULL
  )

  # Validate input files
  missing_files <- gtf_files[!file.exists(gtf_files)]
  if (length(missing_files) > 0)
    stop("GTF file(s) not found:\n",
         paste(" ", missing_files, collapse = "\n"), call. = FALSE)

  message("[omni_build_te_coords] Processing ", length(gtf_files),
          " GTF file(s)...")

  coord_list <- lapply(gtf_files, function(gtf_path) {

    message("[omni_build_te_coords] Loading: ", basename(gtf_path), "...")
    gtf <- as.data.frame(rtracklayer::readGFF(gtf_path))

    # Validate locus_id_col exists
    if (!locus_id_col %in% colnames(gtf))
      stop("Column '", locus_id_col, "' not found in GTF: ",
           basename(gtf_path), ".\n",
           "Available columns: ", paste(colnames(gtf), collapse = ", "),
           call. = FALSE)

    message("[omni_build_te_coords] ", nrow(gtf),
            " rows loaded. Collapsing to locus level...")

    # Filter to standard chromosomes if requested
    if (!is.null(valid_chrs))
      gtf <- gtf[gtf$seqid %in% valid_chrs, ]

    # Collapse multiple exon rows per locus to one coordinate
    # using min(start) and max(end) — the outer boundaries
    gtf$locus_id_field <- gtf[[locus_id_col]]

    coords <- gtf %>%
      dplyr::group_by(locus_id_field, seqid) %>%
      dplyr::summarise(
        Start  = min(start),
        End    = max(end),
        Strand = as.character(dplyr::first(strand)),
        .groups = "drop"
      ) %>%
      dplyr::rename(Locus_ID = locus_id_field,
                    Chr      = seqid) %>%
      as.data.frame()

    message("[omni_build_te_coords]   ",
            nrow(coords), " unique loci.")
    return(coords)
  })

  # Combine all files
  te_coords <- do.call(rbind, coord_list)
  rownames(te_coords) <- NULL

  # Remove duplicates (same locus in multiple files)
  n_before <- nrow(te_coords)
  te_coords <- te_coords[!duplicated(te_coords$Locus_ID), ]
  n_removed <- n_before - nrow(te_coords)
  if (n_removed > 0)
    message("[omni_build_te_coords] Removed ", n_removed,
            " duplicate locus IDs across files.")

  # Parse Class and Family from locus ID if format is "number:Class:Family"
  # Works for HERVd-style IDs like "1000:SINE:MIRb"
  # Sets NA if format is different
  parts          <- strsplit(te_coords$Locus_ID, ":", fixed = TRUE)
  te_coords$Class  <- sapply(parts, function(x)
    if (length(x) >= 2) x[2] else NA_character_)
  te_coords$Family <- sapply(parts, function(x)
    if (length(x) >= 3) x[3] else NA_character_)

  message("\n[omni_build_te_coords] Complete.")
  message("  Total loci: ", nrow(te_coords))
  if (!all(is.na(te_coords$Class)))
    message("  Classes:    ",
            paste(names(table(te_coords$Class)),
                  table(te_coords$Class),
                  sep = " (", collapse = "), "), ")")

  if (!is.null(save_path)) {
    saveRDS(te_coords, save_path)
    message("  Saved to:   ", save_path)
    message("  Reload with: te_coords <- readRDS('", save_path, "')")
  }

  return(te_coords)
}


# =============================================================================
# 3. PREPARE BED COORDINATES FOR GREAT
# =============================================================================

#' Prepare TE Coordinates for GREAT Proximity Analysis
#'
#' Filters significant TEs from the active contrast, looks up their exact
#' genomic coordinates from the HERVd coordinate table, and returns a BED
#' file ready for GREAT submission.
#'
#' Requires a coordinate table built by `omni_build_te_coords()`. This
#' separates the slow GTF-loading step (done once) from the fast
#' filtering step (done per analysis).
#'
#' @param obj A validated OmniSeq object with TE DE results.
#' @param te_coords A data.frame from `omni_build_te_coords()`. Must have
#'   columns: Locus_ID, Chr, Start, End, Strand.
#' @param lfc_thresh Log2 fold change threshold. Default: `1`.
#' @param pval_thresh Adjusted p-value threshold. Default: `0.05`.
#' @param direction Character. `"up"`, `"down"`, or `"both"`. Default: `"up"`.
#' @param save_bed Logical. Save the BED file to disk. Default: `FALSE`.
#' @param bed_path Character. Path to save BED file.
#'
#' @return The OmniSeq object with `$te_bed` populated.
#' @export
#'
#' @examples
#' \dontrun{
#' # Build coordinate table once
#' te_coords <- omni_build_te_coords(
#'   gtf_files = c("/path/to/sine.gtf", "/path/to/line.gtf")
#' )
#'
#' # Then use for each analysis
#' obj <- run_omni_te_proximity(obj, te_coords = te_coords,
#'                               direction = "up")
#' obj <- run_omni_te_great(obj, genome = "hg38")
#' }
run_omni_te_proximity <- function(obj,
                                  te_coords   = NULL,
                                  lfc_thresh  = 1,
                                  pval_thresh = 0.05,
                                  direction   = c("up", "both", "down"),
                                  save_bed    = FALSE,
                                  bed_path    = NULL) {

  validate_omniseq(obj)
  .omni_check_te(obj)
  direction <- match.arg(direction)

  # Validate te_coords
  if (is.null(te_coords))
    stop(paste0(
      "te_coords is required. Build it once with:\n",
      "  te_coords <- omni_build_te_coords(\n",
      "    gtf_files = c('/path/to/te.gtf')\n",
      "  )\n",
      "Then pass it to this function via te_coords = te_coords."
    ), call. = FALSE)

  required_cols <- c("Locus_ID", "Chr", "Start", "End")
  missing_cols  <- setdiff(required_cols, colnames(te_coords))
  if (length(missing_cols) > 0)
    stop("te_coords is missing columns: ",
         paste(missing_cols, collapse = ", "),
         ". Build with omni_build_te_coords().", call. = FALSE)

  # Get DE results
  res <- omni_results(obj)
  if (is.null(res))
    stop("No DE results. Run run_omni_de() first.", call. = FALSE)

  res_df      <- as.data.frame(res)
  res_df$padj <- ifelse(is.na(res_df$padj), 1, res_df$padj)

  # Filter by significance and direction
  sig_mask <- res_df$padj < pval_thresh & !is.na(res_df$log2FoldChange)

  if (direction == "up")
    sig_mask <- sig_mask & res_df$log2FoldChange > lfc_thresh
  else if (direction == "down")
    sig_mask <- sig_mask & res_df$log2FoldChange < -lfc_thresh
  else
    sig_mask <- sig_mask & abs(res_df$log2FoldChange) > lfc_thresh

  sig_df <- res_df[sig_mask, ]

  if (nrow(sig_df) == 0)
    stop("No significant TEs found at lfc_thresh=", lfc_thresh,
         ", pval_thresh=", pval_thresh,
         ", direction='", direction, "'.", call. = FALSE)

  message("[run_omni_te_proximity] ", nrow(sig_df),
          " significant TEs found.")

  # Look up exact coordinates by locus ID
  locus_ids   <- rownames(sig_df)
  coord_match <- te_coords[te_coords$Locus_ID %in% locus_ids, ]

  n_matched <- length(unique(coord_match$Locus_ID))
  n_missing <- length(locus_ids) - n_matched

  message("[run_omni_te_proximity] Coordinate lookup: ",
          n_matched, " / ", length(locus_ids), " loci matched.")

  if (n_missing > 0)
    message("[run_omni_te_proximity] ", n_missing,
            " loci not found in te_coords — these will be excluded.",
            " Ensure te_coords covers all TE classes in this object.")

  if (nrow(coord_match) == 0)
    stop("No loci could be matched to coordinates. ",
         "Check that te_coords was built from the correct GTF files.",
         call. = FALSE)

  # Build BED data frame — one row per unique locus
  bed_df <- coord_match[!duplicated(coord_match$Locus_ID), ]
  bed_df$Score <- sig_df[bed_df$Locus_ID, "log2FoldChange"]

  # Keep standard columns
  bed_df <- data.frame(
    Chr    = bed_df$Chr,
    Start  = bed_df$Start,
    End    = bed_df$End,
    Name   = bed_df$Locus_ID,
    Score  = round(bed_df$Score, 3),
    Strand = bed_df$Strand,
    stringsAsFactors = FALSE
  )

  # Keep only standard chromosomes
  valid_chrs <- paste0("chr", c(1:22, "X", "Y"))
  bed_df     <- bed_df[bed_df$Chr %in% valid_chrs, ]

  message("[run_omni_te_proximity] BED file ready: ",
          nrow(bed_df), " loci on standard chromosomes.")

  if (save_bed) {
    bed_path <- bed_path %||% paste0(
      "TE_significant_", direction, "_",
      gsub("[^A-Za-z0-9]", "_", obj$active_contrast), ".bed"
    )
    utils::write.table(bed_df, file = bed_path, sep = "\t",
                       quote = FALSE, row.names = FALSE,
                       col.names = FALSE)
    message("[run_omni_te_proximity] BED file saved: ", bed_path)
  }

  obj$te_bed <- bed_df
  message("  Access with obj$te_bed")
  message("  Next: run_omni_te_great(obj, genome = 'hg38')")

  return(obj)
}


# =============================================================================
# 4. rGREAT PROXIMITY ENRICHMENT  [FIXED: uses local great() — no web server]
# =============================================================================

#' Run rGREAT Proximity-Based Functional Enrichment
#'
#' Runs GREAT proximity enrichment locally via the `rGREAT::great()` function
#' (rGREAT >= 1.99.0). This avoids the deprecated web-server API which no
#' longer supports named ontology categories in v4.
#'
#' Requires `run_omni_te_proximity()` to be run first.
#'
#' @param obj A validated OmniSeq object with TE BED coordinates.
#' @param genome Character. Genome build: `"hg38"` or `"hg19"`.
#'   Default: `"hg38"`.
#' @param show_n Integer. Number of top GO terms to return. Default: `20`.
#'
#' @return The OmniSeq object with `$te_great` populated.
#' @export
run_omni_te_great <- function(obj,
                              genome = "hg38",
                              show_n = 20) {

  validate_omniseq(obj)
  .require_pkg("rGREAT")
  .require_pkg("GenomicRanges")
  .require_pkg("IRanges")

  bed_df <- obj$te_bed
  if (is.null(bed_df))
    stop("No BED coordinates. Run run_omni_te_proximity() first.",
         call. = FALSE)

  message("\n[run_omni_te_great] Submitting ", nrow(bed_df),
          " TE loci to GREAT (genome: ", genome, ")...")

  # Build GRanges
  gr <- GenomicRanges::GRanges(
    seqnames = bed_df$Chr,
    ranges   = IRanges::IRanges(start = bed_df$Start,
                                end   = bed_df$End)
  )

  # Use local great() — fully offline, no web server dependency.
  # submitGreatJob() is deprecated in rGREAT v4 (availableCategories returns "").
  res <- rGREAT::great(gr,
                       gene_sets  = "GO:BP",
                       tss_source = genome,
                       background = NULL,
                       cores      = 1)

  message("[run_omni_te_great] Retrieving enrichment table...")
  tb <- rGREAT::getEnrichmentTable(res)

  if (is.null(tb) || nrow(tb) == 0)
    stop("GREAT returned no GO Biological Process terms.", call. = FALSE)

  tb <- tb[order(tb$p_adjust), ]
  tb <- utils::head(tb, show_n)

  obj$te_great <- list(
    results  = tb,
    genome   = genome,
    n_loci   = nrow(bed_df),
    contrast = obj$active_contrast
  )

  message("[run_omni_te_great] ", nrow(tb),
          " GO terms retained (top ", show_n, ").")
  message("  Access with obj$te_great")

  return(obj)
}


# =============================================================================
# 5. TE FAMILY BARPLOT
# =============================================================================

#' Plot TE Family Enrichment Barplot
#'
#' Horizontal barplot of the top N most enriched TE families by count of
#' significant loci, colored by TE class.
#'
#' @param obj A validated OmniSeq object with TE family results.
#' @param n_families Integer. Number of top families to display. Default: `20`.
#' @param min_loci Integer. Minimum significant loci to include a family.
#'   Default: `1`.
#' @param color_by Character. `"class"` (color by TE class) or `"lfc"`
#'   (color by mean LFC). Default: `"class"`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_te_family <- function(obj,
                                n_families = 20,
                                min_loci   = 1,
                                color_by   = c("class", "lfc")) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .omni_check_te(obj)

  color_by <- match.arg(color_by)

  if (is.null(obj$te_family))
    stop("No TE family results. Run run_omni_te_family() first.", call. = FALSE)

  fam_df <- obj$te_family$by_family
  fam_df <- fam_df[fam_df$N_Sig >= min_loci, ]

  if (nrow(fam_df) == 0)
    stop("No families with >= ", min_loci, " significant loci.", call. = FALSE)

  top_fams <- utils::head(fam_df[order(-fam_df$N_Sig), ], n_families)
  top_fams$Family <- factor(top_fams$Family,
                            levels = rev(top_fams$Family))

  if (color_by == "class") {
    p <- ggplot2::ggplot(top_fams,
                         ggplot2::aes(x = N_Sig, y = Family,
                                      fill = Class)) +
      ggplot2::geom_col(show.legend = TRUE) +
      ggplot2::scale_fill_viridis_d(option = "turbo", name = "TE Class")

  } else {
    p <- ggplot2::ggplot(top_fams,
                         ggplot2::aes(x = N_Sig, y = Family,
                                      fill = Mean_LFC)) +
      ggplot2::geom_col(show.legend = TRUE) +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                    high = "#B2182B", midpoint = 0,
                                    name = "Mean LFC")
  }

  p <- p +
    ggplot2::geom_text(ggplot2::aes(label = N_Sig),
                       hjust = -0.2, size = 3, color = "grey30") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(
      title    = "TE Family Enrichment",
      subtitle = paste("Contrast:", obj$active_contrast,
                       "| Direction:", obj$te_family$direction),
      x        = "Number of Significant Loci",
      y        = NULL
    ) +
    ggplot2::expand_limits(x = max(top_fams$N_Sig) * 1.12)

  print(p)
  return(invisible(p))
}


# =============================================================================
# 6. TE FAMILY ACTIVITY HEATMAP
# =============================================================================

#' Plot TE Family Activity Heatmap
#'
#' Aggregates VST-normalised TE expression by Family (mean per family per
#' sample) and displays the most variable families as a heatmap. Reveals
#' which TE families are most active across patient groups.
#'
#' @param obj A validated OmniSeq object with TE data and a fitted model.
#' @param n_families Integer. Number of top variable families. Default: `25`.
#' @param cluster_rows Logical. Default: `TRUE`.
#' @param cluster_cols Logical. Default: `TRUE`.
#' @param color_up Color for high expression. Default: `"#B2182B"`.
#' @param color_down Color for low expression. Default: `"#2166AC"`.
#' @param title Optional custom title.
#'
#' @return A pheatmap object (returned invisibly).
#' @export
plot_omni_te_heatmap <- function(obj,
                                 n_families   = 25,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 color_up     = "#B2182B",
                                 color_down   = "#2166AC",
                                 title        = NULL) {

  validate_omniseq(obj)
  .require_pkg("pheatmap")
  .omni_check_te(obj)

  vsd <- omni_vst(obj)
  if (is.null(vsd))
    stop("VST not cached. Run run_omni_de() first.", call. = FALSE)

  mat  <- SummarizedExperiment::assay(vsd)
  meta <- omni_metadata(obj)

  # Parse families from row names
  te_info <- .omni_parse_te_ids(rownames(mat))
  fams    <- te_info$Family

  # Aggregate: mean VST per family per sample
  unique_fams <- unique(fams[fams != "Unknown"])
  agg_list    <- lapply(unique_fams, function(f) {
    idx <- which(fams == f)
    if (length(idx) == 1) mat[idx, ]
    else colMeans(mat[idx, ], na.rm = TRUE)
  })

  agg_mat           <- do.call(rbind, agg_list)
  rownames(agg_mat) <- unique_fams

  # Select top variable families
  pvar      <- apply(agg_mat, 1, stats::var)
  top_fams  <- names(sort(pvar, decreasing = TRUE))[
    1:min(n_families, nrow(agg_mat))]
  mat_plot  <- agg_mat[top_fams, ]

  # Column annotation
  cond_col   <- obj$params$condition
  annot_cols <- intersect(c(cond_col, obj$params$covariates), colnames(meta))
  df_annot   <- meta[, annot_cols, drop = FALSE]
  df_annot[] <- lapply(df_annot,
                       function(x) if (is.factor(x)) droplevels(x) else x)

  pal <- grDevices::colorRampPalette(c(color_down, "white", color_up))(100)

  if (is.null(title))
    title <- paste0("TE Family Activity — ", obj$active_contrast)

  p <- pheatmap::pheatmap(
    mat_plot,
    annotation_col = df_annot,
    color          = pal,
    scale          = "row",
    cluster_rows   = cluster_rows,
    cluster_cols   = cluster_cols,
    border_color   = NA,
    show_rownames  = TRUE,
    fontsize_row   = if (n_families > 40) 6 else 9,
    main           = title
  )

  return(invisible(p))
}


# =============================================================================
# 7. TE CLASS COMPOSITION DONUT
# =============================================================================

#' Plot TE Class Composition Donut Chart
#'
#' Shows the proportion of significant TE loci belonging to each TE class
#' (ERV, LINE, SINE, DNA, etc.) as a donut chart.
#'
#' @param obj A validated OmniSeq object with TE family results.
#' @param show_counts Logical. Show locus counts inside segments. Default: `TRUE`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_te_class_pie <- function(obj, show_counts = TRUE) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .omni_check_te(obj)

  if (is.null(obj$te_family))
    stop("No TE family results. Run run_omni_te_family() first.", call. = FALSE)

  class_df <- obj$te_family$by_class
  class_df <- class_df[class_df$N_Sig > 0, ]

  if (nrow(class_df) == 0)
    stop("No significant TE loci found in any class.", call. = FALSE)

  class_df <- class_df[order(-class_df$N_Sig), ]
  class_df$Fraction <- class_df$N_Sig / sum(class_df$N_Sig)
  class_df$ymax     <- cumsum(class_df$Fraction)
  class_df$ymin     <- c(0, utils::head(class_df$ymax, -1))
  class_df$LabelPos <- (class_df$ymax + class_df$ymin) / 2
  class_df$Label    <- if (show_counts)
    paste0(class_df$Class, "\n(", class_df$N_Sig, ")")
  else
    class_df$Class

  p <- ggplot2::ggplot(class_df,
                       ggplot2::aes(ymax = ymax, ymin = ymin,
                                    xmax = 4, xmin = 2.5,
                                    fill = Class)) +
    ggplot2::geom_rect() +
    ggplot2::geom_text(ggplot2::aes(x = 3.5, y = LabelPos, label = Label),
                       size = 3.5, color = "white", fontface = "bold") +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(c(0, 4)) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(
      title    = "Significant TE Loci by Class",
      subtitle = paste("Contrast:", obj$active_contrast,
                       "| Total:", sum(class_df$N_Sig), "loci"),
      fill     = "TE Class"
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 8. TE VOLCANO COLORED BY CLASS
# =============================================================================

#' Plot TE Volcano Colored by Class
#'
#' Standard volcano plot but with points colored by TE class rather than
#' direction, making it easy to see which classes drive the signal.
#'
#' @param obj A validated OmniSeq object with TE DE results.
#' @param lfc_thresh Log2 fold change threshold. Default: `1`.
#' @param pval_thresh Adjusted p-value threshold. Default: `0.05`.
#' @param top_n Integer. Number of top loci to label. Default: `10`.
#' @param label_col Character. `"family"` or `"locus"`. What to use as
#'   point labels. Default: `"family"`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_te_volcano <- function(obj,
                                 lfc_thresh  = 1,
                                 pval_thresh = 0.05,
                                 top_n       = 10,
                                 label_col   = c("family", "locus")) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .require_pkg("ggrepel")
  .omni_check_te(obj)

  label_col <- match.arg(label_col)

  res <- omni_results(obj)
  if (is.null(res))
    stop("No DE results. Run run_omni_de() first.", call. = FALSE)

  df        <- as.data.frame(res)
  df$padj   <- ifelse(is.na(df$padj), 1, df$padj)
  df        <- df[!is.na(df$log2FoldChange), ]

  # Parse TE info
  te_info  <- .omni_parse_te_ids(rownames(df))
  df$Class  <- te_info$Class
  df$Family <- te_info$Family
  df$Locus  <- rownames(df)

  df$Significant <- df$padj < pval_thresh &
    abs(df$log2FoldChange) > lfc_thresh

  # Label top significant loci
  top_df    <- df[df$Significant, ]
  top_df    <- utils::head(top_df[order(top_df$padj), ], top_n)
  top_df$Label <- if (label_col == "family") top_df$Family else top_df$Locus

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange,
                                        y = -log10(padj),
                                        color = Class,
                                        alpha = Significant)) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.2),
                                guide = "none") +
    ggrepel::geom_text_repel(
      data         = top_df,
      ggplot2::aes(label = Label),
      size         = 3,
      color        = "black",
      box.padding  = 0.4,
      max.overlaps = 20,
      show.legend  = FALSE
    ) +
    ggplot2::scale_color_viridis_d(option = "turbo", name = "TE Class") +
    ggplot2::geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
                        linetype = "dashed", color = "black",
                        linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(pval_thresh),
                        linetype = "dashed", color = "black",
                        linewidth = 0.4) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title    = "TE Volcano Plot",
      subtitle = paste("Contrast:", obj$active_contrast),
      x        = "Log\u2082 Fold Change",
      y        = "-Log\u2081\u2080 Adjusted P-value"
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 9. LFC DENSITY PER TE CLASS
# =============================================================================

#' Plot LFC Density Curves per TE Class
#'
#' Overlapping density curves of log2 fold change, one per TE class.
#' Useful for comparing the overall shift in expression across classes.
#'
#' @param obj A validated OmniSeq object with TE DE results.
#' @param min_loci Integer. Minimum loci for a class to be plotted.
#'   Default: `50`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_te_lfc_density <- function(obj, min_loci = 50) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")
  .omni_check_te(obj)

  res <- omni_results(obj)
  if (is.null(res))
    stop("No DE results. Run run_omni_de() first.", call. = FALSE)

  df      <- as.data.frame(res)
  df      <- df[!is.na(df$log2FoldChange), ]
  te_info <- .omni_parse_te_ids(rownames(df))
  df$Class <- te_info$Class

  # Keep only classes with enough loci
  class_counts <- table(df$Class)
  keep_classes <- names(class_counts[class_counts >= min_loci])
  df           <- df[df$Class %in% keep_classes, ]

  if (nrow(df) == 0)
    stop("No TE classes with >= ", min_loci, " loci.", call. = FALSE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange,
                                        fill = Class,
                                        color = Class)) +
    ggplot2::geom_density(alpha = 0.25, linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "black", linewidth = 0.5) +
    ggplot2::scale_fill_viridis_d(option  = "turbo") +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title    = "LFC Distribution by TE Class",
      subtitle = paste("Contrast:", obj$active_contrast,
                       "| Classes with \u2265", min_loci, "loci"),
      x        = "Log\u2082 Fold Change",
      y        = "Density",
      fill     = "TE Class",
      color    = "TE Class"
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 10. GREAT PROXIMITY DOTPLOT  [FIXED: column names from local great()]
# =============================================================================

#' Plot rGREAT Proximity Enrichment Results
#'
#' Dotplot of top GO Biological Process terms from GREAT proximity analysis.
#' Uses column names produced by `rGREAT::great()` (local mode):
#'   `description`, `p_adjust`, `fold_enrichment`.
#'
#' @param obj A validated OmniSeq object with GREAT results.
#' @param show_n Integer. Number of top terms to display. Default: `15`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_te_great <- function(obj, show_n = 15) {

  validate_omniseq(obj)
  .require_pkg("ggplot2")

  great <- obj$te_great
  if (is.null(great))
    stop("No GREAT results. Run run_omni_te_great() first.", call. = FALSE)

  df      <- great$results
  df      <- utils::head(df[order(df$p_adjust), ], show_n)
  df$term <- factor(df$description, levels = rev(df$description))

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x     = -log10(p_adjust),
                                    y     = term,
                                    size  = fold_enrichment,
                                    color = p_adjust)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low  = "#B2182B",
                                  high = "#4393C3",
                                  name = "FDR") +
    ggplot2::scale_size_continuous(name  = "Fold enrichment",
                                   range = c(2, 8)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(
      title    = "TE Host-Gene Proximity (GREAT)",
      subtitle = paste0("Genome: ", great$genome,
                        " | ", great$n_loci, " TE loci submitted",
                        " | Contrast: ", great$contrast),
      x        = "-Log\u2081\u2080(FDR)",
      y        = "GO Biological Process"
    )

  print(p)
  return(invisible(p))
}


# =============================================================================
# 11. TE LOAD vs GENE SIGNATURE CORRELATION
# =============================================================================

#' Correlate TE Load with a Gene Signature Score
#'
#' Computes a per-sample TE load (sum of VST expression for a given TE class
#' or all TEs) and correlates it against a gene signature score (mean VST of
#' a gene set). Produces a scatter plot with Pearson R and linear fit.
#'
#' @param gene_obj A FlexiSeq/OmniSeq object containing standard gene VST.
#' @param te_obj A FlexiSeq/OmniSeq object containing TE VST.
#' @param marker_genes Character vector of gene symbols for the signature.
#' @param te_class Character. TE class to sum for load: `"ERV"`, `"LINE"`,
#'   `"SINE"`, etc. Use `"ALL"` to sum all TEs. Default: `"ALL"`.
#' @param score_name Character. Y-axis label. Default: `"Marker Score"`.
#'
#' @return A ggplot object (returned invisibly).
#' @export
plot_omni_te_gene_corr <- function(gene_obj,
                                   te_obj,
                                   marker_genes,
                                   te_class   = "ALL",
                                   score_name = "Marker Score") {

  validate_omniseq(gene_obj)
  validate_omniseq(te_obj)
  .require_pkg("ggplot2")
  .require_pkg("ggpubr")

  # 1. Gene signature score (mean VST of marker genes)
  vsd_genes <- omni_vst(gene_obj)
  if (is.null(vsd_genes))
    stop("VST not cached in gene_obj. Run run_omni_de() first.", call. = FALSE)

  mat_genes     <- SummarizedExperiment::assay(vsd_genes)
  valid_markers <- intersect(marker_genes, rownames(mat_genes))

  if (length(valid_markers) == 0)
    stop("None of the marker genes found in gene_obj VST matrix.", call. = FALSE)

  message("[plot_omni_te_gene_corr] Scoring ", length(valid_markers),
          " / ", length(marker_genes), " marker genes...")

  marker_scores <- colMeans(mat_genes[valid_markers, , drop = FALSE],
                            na.rm = TRUE)

  # 2. TE load (sum of VST for target class)
  vsd_te <- omni_vst(te_obj)
  if (is.null(vsd_te))
    stop("VST not cached in te_obj. Run run_omni_de() first.", call. = FALSE)

  mat_te <- SummarizedExperiment::assay(vsd_te)

  if (te_class != "ALL") {
    message("[plot_omni_te_gene_corr] Filtering to TE class: ", te_class)
    te_info  <- .omni_parse_te_ids(rownames(mat_te))
    keep_idx <- te_info$Class == te_class
    mat_te   <- mat_te[keep_idx, , drop = FALSE]
    if (nrow(mat_te) == 0)
      stop("No TEs found for class '", te_class, "'.", call. = FALSE)
  }

  te_load <- colSums(mat_te, na.rm = TRUE)
  message("[plot_omni_te_gene_corr] TE load computed from ",
          nrow(mat_te), " loci.")

  # 3. Merge
  meta     <- omni_metadata(gene_obj)
  cond_col <- gene_obj$params$condition
  common   <- intersect(names(marker_scores), names(te_load))

  if (length(common) < 3)
    stop("Fewer than 3 common samples between gene and TE objects.",
         call. = FALSE)

  df <- data.frame(
    Sample       = common,
    Marker_Score = marker_scores[common],
    TE_Load      = te_load[common],
    Condition    = meta[common, cond_col],
    stringsAsFactors = FALSE
  )

  # 4. Plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = TE_Load,
                                        y = Marker_Score,
                                        color = Condition)) +
    ggplot2::geom_point(size = 4, alpha = 0.85) +
    ggplot2::geom_smooth(method = "lm", color = "black",
                         linetype = "dashed", fill = "grey85",
                         formula  = y ~ x) +
    ggpubr::stat_cor(method         = "pearson",
                     label.x.npc    = "left",
                     label.y.npc    = "top",
                     size           = 4.5,
                     color          = "black") +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      title    = paste(score_name, "vs TE Load"),
      subtitle = paste("TE class:", te_class,
                       "| Markers:", length(valid_markers), "genes"),
      x        = paste0("TE Load (VST sum) — ", te_class),
      y        = paste0(score_name, " (mean VST)")
    )

  print(p)
  return(invisible(p))
}
