# =============================================================================
# Exercise every exported OmniSeq function on real public data.
# Run from the repo root after prepare_test_data.R:
#   Rscript dev/real-data-tests/test_all_functions.R
# Writes workspace/test_results/{function_results.csv, plots/, logs}.
# =============================================================================

suppressPackageStartupMessages({
  library(OmniSeq)
  library(SummarizedExperiment)
})

set.seed(1)
res_dir  <- "workspace/test_results"
plot_dir <- file.path(res_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
proc     <- "workspace/data/processed"

results <- data.frame(section = character(), fn = character(), call = character(),
                      status = character(), seconds = numeric(),
                      n_warnings = integer(), detail = character())

# Run one check. `plot = TRUE` draws the value into a PNG (ggplot/patchwork
# printed; pheatmap/base plots draw during evaluation). Returns the value, or
# NULL on error.
chk <- function(section, fn, call_label, expr, plot = FALSE, widget = FALSE) {
  warns <- character()
  t0    <- Sys.time()
  png_file <- NULL
  if (plot) {
    png_file <- file.path(plot_dir, paste0(sprintf("%03d", nrow(results) + 1), "_",
                                           gsub("[^A-Za-z0-9]+", "_", call_label), ".png"))
    grDevices::png(png_file, width = 10, height = 8, units = "in", res = 90)
  }
  val <- tryCatch(
    withCallingHandlers(
      {
        v <- force(expr)
        if (plot && inherits(v, c("gg", "ggplot", "patchwork"))) {
          print(v)
          # Also render to PDF: some font settings only fail on pdf devices.
          cur <- grDevices::dev.cur()
          grDevices::pdf(tempfile(fileext = ".pdf"))
          tryCatch(print(v), finally = grDevices::dev.off())
          grDevices::dev.set(cur)
        }
        v
      },
      warning = function(w) { warns <<- c(warns, conditionMessage(w))
                              invokeRestart("muffleWarning") },
      message = function(m) invokeRestart("muffleMessage")),
    error = function(e) structure(conditionMessage(e), class = "chk_error"))
  if (plot) grDevices::dev.off()
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  if (inherits(val, "chk_error")) {
    status <- "FAIL"; detail <- gsub("\\s+", " ", unclass(val))
    if (!is.null(png_file)) unlink(png_file)
    out <- NULL
  } else {
    status <- "PASS"
    detail <- if (widget) paste("htmlwidget:", paste(class(val)[1:2], collapse = "/")) else
      if (!is.null(png_file)) basename(png_file) else
        paste(class(val)[1], collapse = "")
    out <- val
  }
  if (length(warns))
    detail <- paste0(detail, " | warnings: ",
                     paste(unique(substr(warns, 1, 120)), collapse = " || "))
  results[nrow(results) + 1, ] <<- list(section, fn, call_label, status, secs,
                                        length(warns), substr(detail, 1, 600))
  cat(sprintf("[%s] %-38s %-5s %6.1fs %s\n", section, call_label, status, secs,
              if (status == "FAIL") substr(detail, 1, 160) else ""))
  invisible(out)
}
skip <- function(section, fn, why) {
  results[nrow(results) + 1, ] <<- list(section, fn, fn, "SKIP", 0, 0L, why)
  cat(sprintf("[%s] %-38s SKIP  %s\n", section, fn, why))
}
keep <- function(new, old) if (is.null(new)) old else new

# =============================================================================
# A. OVARIAN GENE DATA (GSE182119)
# =============================================================================
ov <- readRDS(file.path(proc, "ovarian_gene_te.rds"))
S  <- "A-gene"

chk(S, "omni_check_deps", "omni_check_deps()", utils::capture.output(omni_check_deps()))
chk(S, "omni_list_signatures", "omni_list_signatures()", omni_list_signatures())
for (sig in c("omniseq_isgs_28", "omniseq_isgs_extended", "omniseq_isgs_schoggins",
              "omniseq_ribosomal", "omniseq_mitochondrial", "omniseq_immunoglobulin"))
  chk(S, sig, sig, { x <- get(sig); stopifnot(is.character(x), length(x) > 0); length(x) })

g <- chk(S, "init_omniseq", "init_omniseq(gene)",
         init_omniseq(ov$gene_counts, ov$meta, condition_col = "Treatment",
                      covariates = "CellLine", ref_level = "Mock",
                      filter_coding = TRUE))
chk(S, "validate_omniseq", "validate_omniseq()", validate_omniseq(g))
chk(S, "print.OmniSeq", "print(obj)", utils::capture.output(print(g)))
chk(S, "summary.OmniSeq", "summary(obj)", utils::capture.output(summary(g)))
chk(S, "omni_counts", "omni_counts()", dim(omni_counts(g)))
chk(S, "omni_metadata", "omni_metadata()", dim(omni_metadata(g)))

g <- keep(chk(S, "run_omni_qc", "run_omni_qc()", run_omni_qc(g)), g)
chk(S, "omni_qc", "omni_qc()", names(omni_qc(g)))
chk(S, "plot_omni_qc_basics", "qc_basics(library_size)",
    plot_omni_qc_basics(g, "library_size"), plot = TRUE)
chk(S, "plot_omni_qc_basics", "qc_basics(detected_genes)",
    plot_omni_qc_basics(g, "detected_genes"), plot = TRUE)
for (pt in c("distance", "mean_var", "pca_compare"))
  chk(S, "plot_omni_qc_advanced", paste0("qc_advanced(", pt, ")"),
      plot_omni_qc_advanced(g, pt), plot = TRUE)

g <- keep(chk(S, "run_omni_de", "run_omni_de(apeglm)",
              run_omni_de(g, shrinkage = "apeglm")), g)
chk(S, "omni_model", "omni_model()", class(omni_model(g)))
chk(S, "omni_vst", "omni_vst()", class(omni_vst(g)))
chk(S, "omni_results", "omni_results()", nrow(omni_results(g)))
main_key <- g$active_contrast

for (pr in list(c("AZA", "Mock"), c("ITF", "Mock"), c("AZA_ITF", "ITF"),
                c("AZA_ITF", "AZA")))
  g <- keep(chk(S, "run_omni_contrast", paste0("run_omni_contrast(", pr[1], " vs ", pr[2], ")"),
                run_omni_contrast(g, group1 = pr[1], group2 = pr[2], shrinkage = "ashr")), g)
g <- keep(chk(S, "run_omni_contrast", "run_omni_contrast(new column CellLine)",
              run_omni_contrast(g, column = "CellLine", group1 = "TykNu",
                                group2 = "A2780", shrinkage = "none")), g)
chk(S, "omni_contrasts", "omni_contrasts()", omni_contrasts(g))
g <- keep(chk(S, "omni_set_contrast", "omni_set_contrast(AZA vs Mock)",
              omni_set_contrast(g, "AZA vs Mock")), g)

chk(S, "plot_omni_volcano", "plot_omni_volcano()", plot_omni_volcano(g), plot = TRUE)
chk(S, "plot_omni_ma", "plot_omni_ma()", plot_omni_ma(g), plot = TRUE)
chk(S, "plot_omni_heatmap", "plot_omni_heatmap(lfc 1)",
    plot_omni_heatmap(g, ntop = 40, lfc_thresh = 1), plot = TRUE)
chk(S, "plot_omni_heatmap", "plot_omni_heatmap(custom ISGs)",
    plot_omni_heatmap(g, custom_genes = omniseq_isgs_28), plot = TRUE)
chk(S, "plot_omni_pca", "plot_omni_pca(color/shape)",
    plot_omni_pca(g, color_by = "Treatment", shape_by = "CellLine"), plot = TRUE)
chk(S, "plot_omni_top_genes", "plot_omni_top_genes(lfc 1)",
    plot_omni_top_genes(g, lfc_thresh = 1, color_by = "Treatment"), plot = TRUE)
chk(S, "plot_omni_pvalue_hist", "plot_omni_pvalue_hist()", plot_omni_pvalue_hist(g), plot = TRUE)

# ---- Interaction term: ~ CellLine + AZA * ITF (2x2 drugs, blocked by line) --
# The package has no constructor for interaction contrasts; build one with
# DESeq2 the way a user would and cache it in results_cache.
inter_key <- "AZA x ITF interaction"
chk(S, "(setup)", "build AZA:ITF interaction contrast", {
  dds <- DESeq2::DESeqDataSetFromMatrix(omni_counts(g), omni_metadata(g),
                                        ~ CellLine + AZA * ITF)
  dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
  r   <- DESeq2::results(dds, name = "AZAyes.ITFyes")
  fc  <- sort(stats::setNames(r$log2FoldChange, rownames(r))[!is.na(r$log2FoldChange)],
              decreasing = TRUE)
  st  <- stats::setNames(r$stat, rownames(r))[!is.na(r$stat)]
  g$results_cache[[inter_key]] <- list(results = r, fold_changes = fc, stat = st,
                                        shrinkage = "none", timestamp = Sys.time())
  sum(r$padj < 0.05, na.rm = TRUE)
})

chk(S, "plot_omni_volcano_interactive", "plot_omni_volcano_interactive()",
    plot_omni_volcano_interactive(g), widget = TRUE)
chk(S, "plot_omni_lfc_scatter_interactive", "lfc_scatter_interactive(AZA vs ITF)",
    plot_omni_lfc_scatter_interactive(g, "AZA vs Mock", "ITF vs Mock"), widget = TRUE)
chk(S, "plot_omni_interaction_scatter", "interaction_scatter()",
    plot_omni_interaction_scatter(g, "AZA_ITF vs ITF", "AZA vs Mock", inter_key,
                                  label_a = "AZA | ITF", label_b = "AZA | no ITF"),
    widget = TRUE)
chk(S, "plot_omni_interaction_scatter_static", "interaction_scatter_static()",
    plot_omni_interaction_scatter_static(g, "AZA_ITF vs ITF", "AZA vs Mock", inter_key),
    plot = TRUE)
chk(S, "plot_omni_interaction_dumbbell", "interaction_dumbbell(padj)",
    plot_omni_interaction_dumbbell(g, "AZA vs Mock", "AZA_ITF vs ITF", inter_key,
                                   label_a = "no ITF", label_b = "with ITF",
                                   stim_label = "AZA"), plot = TRUE)
chk(S, "omni_gene_search_table", "omni_gene_search_table()",
    omni_gene_search_table(g, pval_thresh = 0.05), widget = TRUE)

g <- keep(chk(S, "run_omni_mask", "run_omni_mask(isgs_28)",
              run_omni_mask(g, signature = "isgs_28")), g)
masked_key <- g$active_contrast
g <- omni_set_contrast(g, "AZA vs Mock")

# ---- Enrichment ----------------------------------------------------------------
g <- keep(chk(S, "run_omni_enrichment", "run_omni_enrichment(GO BP)",
              run_omni_enrichment(g, lfc_thresh = 1, simplify = TRUE)), g)
g <- keep(chk(S, "run_omni_simplify", "run_omni_simplify()", run_omni_simplify(g)), g)
for (pt in c("dot", "cnet", "emap"))
  chk(S, "plot_omni_ora", paste0("plot_omni_ora(", pt, ")"),
      plot_omni_ora(g, plot_type = pt), plot = TRUE)
# AZA alone is weak in these heterogeneous lines (139 DEGs; ashr shrinks most
# LFCs to ~0, leaving GSEA a flat ranking), so GSEA runs on AZA+ITF vs Mock.
g <- omni_set_contrast(g, main_key)
g <- keep(chk(S, "run_omni_gsea", "run_omni_gsea(H, AZA_ITF vs Mock)",
              run_omni_gsea(g, category = "H")), g)
chk(S, "run_omni_gsea", "run_omni_gsea(C2 CP:REACTOME)",
    run_omni_gsea(g, category = "C2", subcategory = "CP:REACTOME"))
for (pt in c("waterfall", "ridge", "bar"))
  chk(S, "plot_omni_gsea", paste0("plot_omni_gsea(", pt, ")"),
      plot_omni_gsea(g, plot_type = pt), plot = TRUE)
chk(S, "plot_omni_leading_edge", "plot_omni_leading_edge()",
    plot_omni_leading_edge(g, n_pathways = 4), plot = TRUE)
chk(S, "plot_omni_pathway_lfc", "plot_omni_pathway_lfc(top pathway)",
    plot_omni_pathway_lfc(g, pathway_name = as.data.frame(g$gsea_cache[[main_key]])$ID[1]),
    plot = TRUE)
g <- omni_set_contrast(g, "AZA vs Mock")
g <- keep(chk(S, "run_omni_ssgsea", "run_omni_ssgsea(H)", run_omni_ssgsea(g)), g)
chk(S, "plot_omni_ssgsea_heatmap", "plot_omni_ssgsea_heatmap()",
    plot_omni_ssgsea_heatmap(g, show_n = 25), plot = TRUE)
g <- keep(chk(S, "run_omni_pathway_de", "run_omni_pathway_de(AZA vs Mock)",
              run_omni_pathway_de(g, group1 = "AZA", group2 = "Mock")), g)
chk(S, "plot_omni_compare_clusters", "plot_omni_compare_clusters()",
    plot_omni_compare_clusters(g), plot = TRUE)
chk(S, "plot_omni_multicontrast", "plot_omni_multicontrast(3 contrasts)",
    plot_omni_multicontrast(g, contrasts = c("AZA vs Mock", "ITF vs Mock", main_key)),
    plot = TRUE)
ud_main <- chk(S, "run_omni_ora_updown", "run_omni_ora_updown(AZA_ITF vs Mock)",
               run_omni_ora_updown(g, contrast = main_key))
ud_itf  <- chk(S, "run_omni_ora_updown", "run_omni_ora_updown(AZA_ITF vs ITF)",
               run_omni_ora_updown(g, contrast = "AZA_ITF vs ITF"))
chk(S, "run_omni_ora_updown", "run_omni_ora_updown(ITF vs Mock: 0 DEGs -> clear error)",
    tryCatch(run_omni_ora_updown(g, contrast = "ITF vs Mock"),
             error = function(e) if (grepl("No significant genes", conditionMessage(e)))
               "clear error (expected)" else stop(e)))
for (db in c("GO-BP", "Hallmark"))
  chk(S, "plot_omni_enrichment_interactive", paste0("enrichment_interactive(", db, ")"),
      plot_omni_enrichment_interactive(ud_main[[db]], title = db), widget = TRUE)
chk(S, "plot_omni_enrichment_comparison", "enrichment_comparison(GO-BP x2)",
    plot_omni_enrichment_comparison(list(`AZA+ITF vs Mock` = ud_main[["GO-BP"]],
                                         `AZA+ITF vs ITF`  = ud_itf[["GO-BP"]])),
    plot = TRUE)

# ---- TF activity (CollecTRI fetched from the OmniPath web API) ---------------
# OmnipathR cannot reach Ensembl from this container, so seed run_omni_tf()'s
# own regulon cache with CollecTRI downloaded directly from omnipathdb.org.
ct <- utils::read.delim("workspace/data/raw/collectri_human.tsv")
g$params$tf_regulons <- data.frame(
  source = ct$source_genesymbol, target = ct$target_genesymbol,
  mor = ifelse(ct$is_inhibition == "True" & ct$is_stimulation != "True", -1, 1))
g$params$tf_regulons <- g$params$tf_regulons[!duplicated(
  paste(g$params$tf_regulons$source, g$params$tf_regulons$target)), ]
g$params$tf_database <- "collectri"
g$params$tf_organism <- "human"

for (k in c("AZA vs Mock", "ITF vs Mock", "AZA_ITF vs ITF")) {
  g <- omni_set_contrast(g, k)
  g <- keep(chk(S, "run_omni_tf", paste0("run_omni_tf(contrast: ", k, ")"),
                run_omni_tf(g, mode = "contrast")), g)
}
chk(S, "plot_omni_tf_activity", "plot_omni_tf_activity(contrast)",
    plot_omni_tf_activity(g, n = 20), plot = TRUE)
chk(S, "plot_omni_tf_network", "plot_omni_tf_network(static)",
    plot_omni_tf_network(g, plot_type = "static"), plot = TRUE)
chk(S, "plot_omni_tf_network", "plot_omni_tf_network(interactive)",
    plot_omni_tf_network(g, plot_type = "interactive"), widget = TRUE)
chk(S, "plot_omni_tf_genotype", "plot_omni_tf_genotype()",
    plot_omni_tf_genotype(g, main_contrast = "AZA vs Mock",
                          geno_contrasts = list(ITF = "ITF vs Mock",
                                                `AZA+ITF` = "AZA_ITF vs ITF"),
                          save_dir = file.path(plot_dir, "tf_genotype")), plot = TRUE)
g <- keep(chk(S, "run_omni_tf", "run_omni_tf(samples: AZA)",
              run_omni_tf(g, mode = "samples", sample_group = "AZA")), g)
chk(S, "plot_omni_tf_activity", "plot_omni_tf_activity(samples)",
    plot_omni_tf_activity(g, group = "AZA"), plot = TRUE)
g <- omni_set_contrast(g, "AZA vs Mock")

# ---- STRING ---------------------------------------------------------------------
g <- keep(chk(S, "run_omni_string", "run_omni_string(top 60)",
              run_omni_string(g, n = 60, cache_dir = "workspace/data/raw/string_cache")), g)
chk(S, "plot_omni_string_network", "plot_omni_string_network(static)",
    plot_omni_string_network(g, plot_type = "static"), plot = TRUE)
chk(S, "plot_omni_string_network", "plot_omni_string_network(interactive)",
    plot_omni_string_network(g, plot_type = "interactive"), widget = TRUE)

# ---- Deconvolution (bulk) --------------------------------------------------------
g2 <- chk(S, "run_omni_deconvolution", "run_omni_deconvolution(immunedeconv quantiseq)",
          run_omni_deconvolution(g, method = "immunedeconv", deconv_method = "quantiseq"))
if (!is.null(g2)) {
  g <- g2
  chk(S, "plot_omni_deconv_heatmap", "plot_omni_deconv_heatmap()",
      plot_omni_deconv_heatmap(g), plot = TRUE)
  chk(S, "plot_omni_deconv_corr", "plot_omni_deconv_corr(IFI27)",
      plot_omni_deconv_corr(g, feature = "IFI27"), plot = TRUE)
}

# ---- WGCNA -------------------------------------------------------------------------
g <- keep(chk(S, "run_omni_wgcna", "run_omni_wgcna(3000 genes)",
              run_omni_wgcna(g, n_top_genes = 3000, min_module_size = 20, key = "full")), g)
chk(S, "omni_wgcna", "omni_wgcna('full')", names(omni_wgcna(g, "full")))
chk(S, "plot_omni_wgcna_dendrogram", "plot_omni_wgcna_dendrogram()",
    plot_omni_wgcna_dendrogram(g, key = "full"), plot = TRUE)
chk(S, "plot_omni_wgcna_trait_heatmap", "plot_omni_wgcna_trait_heatmap()",
    plot_omni_wgcna_trait_heatmap(g, key = "full"), plot = TRUE)
chk(S, "plot_omni_wgcna_hub_network", "plot_omni_wgcna_hub_network(static)",
    plot_omni_wgcna_hub_network(g, key = "full", plot_type = "static"), plot = TRUE)
chk(S, "plot_omni_wgcna_hub_network", "plot_omni_wgcna_hub_network(interactive)",
    plot_omni_wgcna_hub_network(g, key = "full", plot_type = "interactive"), widget = TRUE)
g <- keep(chk(S, "run_omni_wgcna_multi", "run_omni_wgcna_multi(AZA)",
              run_omni_wgcna_multi(g, group_col = "AZA", n_top_genes = 3000,
                                   min_module_size = 20)), g)
g <- keep(chk(S, "run_omni_module_preservation", "run_omni_module_preservation(no->yes)",
              run_omni_module_preservation(g, group_col = "AZA", reference = "no",
                                           test = "yes", n_permutations = 50)), g)
chk(S, "plot_omni_module_preservation", "plot_omni_module_preservation()",
    plot_omni_module_preservation(g, group_col = "AZA"), plot = TRUE)

# ---- sPLS-DA -------------------------------------------------------------------------
g <- keep(chk(S, "run_omni_splsda", "run_omni_splsda(AZA)",
              run_omni_splsda(g, group_col = "AZA")), g)
skey <- names(g$splsda)[1]
chk(S, "omni_splsda", "omni_splsda(key)", names(omni_splsda(g, skey)))
chk(S, "plot_omni_splsda_scores", "plot_omni_splsda_scores()",
    plot_omni_splsda_scores(g, key = skey), plot = TRUE)
chk(S, "plot_omni_splsda_loadings", "plot_omni_splsda_loadings()",
    plot_omni_splsda_loadings(g, key = skey), plot = TRUE)
chk(S, "plot_omni_splsda_vip", "plot_omni_splsda_vip()",
    plot_omni_splsda_vip(g, key = skey), plot = TRUE)
chk(S, "plot_omni_splsda_roc", "plot_omni_splsda_roc()",
    plot_omni_splsda_roc(g, key = skey), plot = TRUE)
chk(S, "omni_splsda_signature", "omni_splsda_signature()",
    nrow(omni_splsda_signature(g, key = skey)))

saveRDS(g, file.path(res_dir, "ovarian_gene_obj.rds"))

# =============================================================================
# B. OVARIAN TE DATA + GENE x TE MULTI-OMICS
# =============================================================================
S  <- "B-TE"
te <- chk(S, "init_omniseq", "init_omniseq(TE loci)",
          init_omniseq(ov$te_counts, ov$meta, condition_col = "Treatment",
                       covariates = "CellLine", ref_level = "Mock"))
chk(S, "(check)", "data_type detected as TE",
    stopifnot(identical(te$params$data_type, "TE")))
te <- keep(chk(S, "run_omni_de", "run_omni_de(TE)", run_omni_de(te, shrinkage = "none")), te)
te <- keep(chk(S, "run_omni_contrast", "run_omni_contrast(TE AZA vs Mock)",
               run_omni_contrast(te, group1 = "AZA", group2 = "Mock", shrinkage = "ashr")), te)
te <- keep(chk(S, "run_omni_te_family", "run_omni_te_family()",
               run_omni_te_family(te, lfc_thresh = 1)), te)
chk(S, "plot_omni_te_family", "plot_omni_te_family()", plot_omni_te_family(te), plot = TRUE)
chk(S, "plot_omni_te_heatmap", "plot_omni_te_heatmap()", plot_omni_te_heatmap(te), plot = TRUE)
chk(S, "plot_omni_te_class_pie", "plot_omni_te_class_pie()",
    plot_omni_te_class_pie(te), plot = TRUE)
chk(S, "plot_omni_te_volcano", "plot_omni_te_volcano()", plot_omni_te_volcano(te), plot = TRUE)
chk(S, "plot_omni_te_lfc_density", "plot_omni_te_lfc_density()",
    plot_omni_te_lfc_density(te, min_loci = 20), plot = TRUE)
coords <- chk(S, "omni_build_te_coords", "omni_build_te_coords(Telescope GTF)",
              omni_build_te_coords(file.path(proc, "ovarian_te_loci.gtf")))
if (!is.null(coords)) {
  te <- keep(chk(S, "run_omni_te_proximity", "run_omni_te_proximity(up)",
                 run_omni_te_proximity(te, te_coords = coords, direction = "up")), te)
  te <- keep(chk(S, "run_omni_te_great", "run_omni_te_great(hg38)",
                 run_omni_te_great(te, genome = "hg38")), te)
  chk(S, "plot_omni_te_great", "plot_omni_te_great()", plot_omni_te_great(te), plot = TRUE)
}

gm <- omni_set_contrast(g, "AZA vs Mock")
chk(S, "plot_omni_te_gene_corr", "plot_omni_te_gene_corr(ISGs, ALL)",
    plot_omni_te_gene_corr(gm, te, marker_genes = omniseq_isgs_28), plot = TRUE)

# Sensor scoring: package defaults first (they target RepeatMasker family
# names), then families named as in the Telescope annotation.
gs <- chk(S, "run_omni_sensor_score", "run_omni_sensor_score(defaults)",
          run_omni_sensor_score(gm, te, group_col = "Treatment"))
gs <- keep(chk(S, "run_omni_sensor_score", "run_omni_sensor_score(Telescope families)",
               run_omni_sensor_score(gm, te, group_col = "Treatment",
                                     cgas_families = c("L1FLnI", "L1FLI", "HML2", "HERVH", "HERVW"),
                                     rigi_families = c("ERV316A3", "ERVLE", "ERVLB4", "MER4", "HERVL"),
                                     tlr7_families = c("L1FLI", "L1ORF2", "HERVK", "HML3"))), gs)
if (!is.null(gs)) {
  chk(S, "plot_omni_sensor_scores", "plot_omni_sensor_scores()",
      plot_omni_sensor_scores(gs), plot = TRUE)
  chk(S, "plot_omni_sensor_irgs_corr", "plot_omni_sensor_irgs_corr()",
      plot_omni_sensor_irgs_corr(gs, color_col = "Treatment", shape_col = "CellLine"), plot = TRUE)
  chk(S, "plot_omni_sensor_heatmap", "plot_omni_sensor_heatmap()",
      plot_omni_sensor_heatmap(gs), plot = TRUE)
  sc_df <- gs$sensor_scores$scores   # now carries every metadata column
  gs <- keep(chk(S, "run_omni_mediation", "run_omni_mediation(AZA -> cGAS -> IRG-S)",
                 run_omni_mediation(gs, sensor_scores = sc_df, mediator = "cGAS_STING",
                                    treatment_col = "AZA",
                                    treatment_vals = c("0" = "no", "1" = "yes"),
                                    exclude_groups = NULL, boot_sims = 200)), gs)
  gs <- keep(chk(S, "run_omni_mediation_multi", "run_omni_mediation_multi(AZA)",
                 run_omni_mediation_multi(gs, sensor_scores = sc_df, treatment_col = "AZA",
                                          treatment_vals = c("0" = "no", "1" = "yes"),
                                          exclude_groups = NULL, boot_sims = 200)), gs)
  chk(S, "plot_omni_mediation", "plot_omni_mediation()", plot_omni_mediation(gs), plot = TRUE)
  chk(S, "plot_omni_mediation_path", "plot_omni_mediation_path()",
      plot_omni_mediation_path(gs, mediator = "cGAS_STING"), plot = TRUE)
}
gd <- chk(S, "run_omni_diablo", "run_omni_diablo(gene + TE, AZA)",
          run_omni_diablo(gm, te, group_col = "AZA", ncomp = 2))
gw <- chk(S, "run_omni_wgcna_te_traits", "run_omni_wgcna_te_traits(full)",
          run_omni_wgcna_te_traits(gm, te, wgcna_key = "full", n_top_families = 30))
if (!is.null(gw))
  chk(S, "plot_omni_wgcna_te_trait_heatmap", "plot_omni_wgcna_te_trait_heatmap()",
      plot_omni_wgcna_te_trait_heatmap(gw, key = "full"), plot = TRUE)

# ---- Export ------------------------------------------------------------------
out_dir <- file.path(res_dir, "export")
chk(S, "export_omniseq", "export_omniseq(gene)", {
  export_omniseq(gm, outdir = file.path(out_dir, "gene"), dpi = 72)
  n <- length(list.files(file.path(out_dir, "gene"), recursive = TRUE))
  stopifnot(n > 10)
  paste(n, "files written")
})
chk(S, "export_omniseq", "export_omniseq(TE, gene_obj)",
    export_omniseq(te, outdir = file.path(out_dir, "te"), gene_obj = gm, dpi = 72))
chk(S, "export_omniseq_report", "export_omniseq_report(gene, relative outdir)", {
  html <- export_omniseq_report(gm, outdir = file.path(out_dir, "report"), open = FALSE)
  stopifnot(!is.null(html), file.exists(html), file.size(html) > 1e5)
  html
})
saveRDS(te, file.path(res_dir, "ovarian_te_obj.rds"))

# ---- Pipelines -----------------------------------------------------------------
S <- "C-pipelines"
chk(S, "run_omniseq_auto", "run_omniseq_auto(gene + TE)",
    run_omniseq_auto(ov$gene_counts, ov$meta, condition_col = "Treatment",
                     covariates = "CellLine", ref_level = "Mock", te_obj = te,
                     outdir = file.path(res_dir, "auto"), n_top_genes = 2000))
chk(S, "run_omniseq_masked", "run_omniseq_masked(isgs_28)",
    run_omniseq_masked(gm, signature = "isgs_28", n_wgcna = 2000,
                       outdir = file.path(res_dir, "masked")))
chk(S, "run_omniseq_guided", "run_omniseq_guided() non-interactive guard",
    tryCatch(run_omniseq_guided(ov$gene_counts, ov$meta, condition_col = "Treatment"),
             error = function(e) if (grepl("interactive", conditionMessage(e)))
               "refuses non-interactive session (expected)" else stop(e)))
skip(S, "omni_install_deps", "installs packages system-wide; not run")
skip(S, "omni_launch_app", "tested separately by the app smoke tests")

# =============================================================================
# D. PANCREATIC ISLETS: BULK + scRNA REFERENCE (GSE50244 + Xin et al.)
# =============================================================================
S   <- "D-singlecell"
pan <- readRDS(file.path(proc, "pancreas_bulk.rds"))
sce <- readRDS(file.path(proc, "pancreas_sc_xin.rds"))
p <- chk(S, "init_omniseq", "init_omniseq(islets T2D)",
         init_omniseq(pan$counts, pan$meta, condition_col = "T2D",
                      covariates = "gender", ref_level = "Normal"))
p <- keep(chk(S, "run_omni_de", "run_omni_de(islets)", run_omni_de(p, shrinkage = "none")), p)
for (m in c("signature", "bisque", "music", "auto"))   # music/auto need MuSiC
  p <- keep(chk(S, "run_omni_sc_deconv", paste0("run_omni_sc_deconv(", m, ")"),
                run_omni_sc_deconv(p, sce, cell_type_col = "cell_type",
                                   sample_col = "donor", group_col = "group",
                                   ref_group = "Control", method = m,
                                   ref_label = "Xin 2016 islets")), p)
p <- keep(chk(S, "run_omni_disease_state", "run_omni_disease_state(T2D vs Control)",
              run_omni_disease_state(p, sce, disease_group = "T2D",
                                     healthy_group = "Control",
                                     cell_type_col = "cell_type", group_col = "group")), p)
chk(S, "plot_omni_sc_heatmap", "plot_omni_sc_heatmap()", plot_omni_sc_heatmap(p), plot = TRUE)
chk(S, "plot_omni_sc_genotype", "plot_omni_sc_genotype(T2D)",
    plot_omni_sc_genotype(p, genotype_col = "T2D"), plot = TRUE)
chk(S, "plot_omni_sc_correlation", "plot_omni_sc_correlation(INS)",
    plot_omni_sc_correlation(p, feature = "INS"), plot = TRUE)
chk(S, "(check)", "bisque proportions stored",
    round(rowMeans(p$sc_deconv$bisque_proportions), 3))
chk(S, "plot_omni_sc_profile", "plot_omni_sc_profile(T2D)",
    plot_omni_sc_profile(p, genotype_col = "T2D"), plot = TRUE)
chk(S, "plot_omni_sc_two_layer", "plot_omni_sc_two_layer(T2D)",
    plot_omni_sc_two_layer(p, genotype_col = "T2D"), plot = TRUE)
chk(S, "plot_omni_sc_benchmark", "plot_omni_sc_benchmark()",
    plot_omni_sc_benchmark(p), plot = TRUE)
p2 <- chk(S, "run_omni_deconvolution", "run_omni_deconvolution(music, sc_ref)",
          run_omni_deconvolution(p, method = "music", sc_ref = sce,
                                 cell_type_col = "cell_type", sample_col = "donor"))
if (!is.null(p2)) {
  chk(S, "plot_omni_deconv_heatmap", "plot_omni_deconv_heatmap(music)",
      plot_omni_deconv_heatmap(p2), plot = TRUE)
  chk(S, "plot_omni_deconv_corr", "plot_omni_deconv_corr(INS)",
      plot_omni_deconv_corr(p2, feature = "INS"), plot = TRUE)
}
saveRDS(p, file.path(res_dir, "pancreas_obj.rds"))

# =============================================================================
write.csv(results, file.path(res_dir, "function_results.csv"), row.names = FALSE)
exported <- getNamespaceExports("OmniSeq")
covered  <- unique(results$fn)
cat("\n==== SUMMARY ====\n")
print(table(results$status))
cat("Exported functions not exercised:",
    paste(sort(setdiff(exported, covered)), collapse = ", "), "\n")
