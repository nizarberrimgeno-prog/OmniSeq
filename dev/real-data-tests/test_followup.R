# Re-run the checks that failed or were not reached in test_all_functions.R,
# against the patched package. Reuses chk() from that script.
suppressPackageStartupMessages({ library(OmniSeq); library(SummarizedExperiment) })
src <- readLines("dev/real-data-tests/test_all_functions.R")
eval(parse(text = src[1:grep("^keep <- function", src)]))   # setup + chk()/skip()/keep()
plot_dir <- file.path(res_dir, "plots_followup")
dir.create(plot_dir, showWarnings = FALSE)

g  <- readRDS(file.path(res_dir, "ovarian_gene_obj.rds"))
te <- readRDS(file.path(res_dir, "ovarian_te_obj.rds"))
ov <- readRDS(file.path(proc, "ovarian_gene_te.rds"))
S  <- "followup"

# GSEA / Up-Down ORA on contrasts that carry signal (AZA alone is weak here).
main_key <- "Treatment_AZA_ITF_vs_Mock"
g <- omni_set_contrast(g, main_key)
g <- keep(chk(S, "run_omni_gsea", "run_omni_gsea(H, AZA_ITF vs Mock)",
              run_omni_gsea(g, category = "H")), g)
chk(S, "run_omni_gsea", "run_omni_gsea(C2 CP:REACTOME subcollection)",
    nrow(as.data.frame(run_omni_gsea(g, category = "C2",
                                     subcategory = "CP:REACTOME")$advanced_gsea)))
for (pt in c("waterfall", "ridge", "bar"))
  chk(S, "plot_omni_gsea", paste0("plot_omni_gsea(", pt, ")"),
      plot_omni_gsea(g, plot_type = pt), plot = TRUE)
chk(S, "plot_omni_leading_edge", "plot_omni_leading_edge()",
    plot_omni_leading_edge(g, n_pathways = 4), plot = TRUE)
top_path <- as.data.frame(g$gsea_cache[[main_key]])$ID[1]
chk(S, "plot_omni_pathway_lfc", paste0("plot_omni_pathway_lfc(", top_path, ")"),
    plot_omni_pathway_lfc(g, pathway_name = top_path), plot = TRUE)

ud_main <- chk(S, "run_omni_ora_updown", "run_omni_ora_updown(AZA_ITF vs Mock)",
               run_omni_ora_updown(g, contrast = main_key))
ud_itf  <- chk(S, "run_omni_ora_updown", "run_omni_ora_updown(AZA_ITF vs ITF)",
               run_omni_ora_updown(g, contrast = "AZA_ITF vs ITF"))
for (db in c("GO-BP", "Hallmark"))
  if (!is.null(ud_main[[db]]))
    chk(S, "plot_omni_enrichment_interactive", paste0("enrichment_interactive(", db, ")"),
        plot_omni_enrichment_interactive(ud_main[[db]], title = db), widget = TRUE)
chk(S, "plot_omni_enrichment_comparison", "enrichment_comparison(GO-BP, 2 contrasts)",
    plot_omni_enrichment_comparison(list(`AZA+ITF vs Mock` = ud_main[["GO-BP"]],
                                         `AZA+ITF vs ITF`  = ud_itf[["GO-BP"]])),
    plot = TRUE)
chk(S, "run_omni_ora_updown", "run_omni_ora_updown(ITF vs Mock: 0 DEGs, expect error)",
    tryCatch(run_omni_ora_updown(g, contrast = "ITF vs Mock"),
             error = function(e) if (grepl("No significant genes", conditionMessage(e)))
               "clear error for a contrast without DEGs (expected)" else stop(e)))

# TF genotype: geno_contrasts must be a *named* list
chk(S, "plot_omni_tf_genotype", "plot_omni_tf_genotype(unnamed -> clear error)",
    tryCatch(plot_omni_tf_genotype(g, main_contrast = "AZA vs Mock",
                                   geno_contrasts = list("ITF vs Mock", "AZA_ITF vs ITF")),
             error = function(e) if (grepl("must be named", conditionMessage(e)))
               "clear error (expected)" else stop(e)))
chk(S, "plot_omni_tf_genotype", "plot_omni_tf_genotype(named list)",
    plot_omni_tf_genotype(g, main_contrast = "AZA vs Mock",
                          geno_contrasts = list(ITF = "ITF vs Mock",
                                                `AZA+ITF` = "AZA_ITF vs ITF"),
                          save_dir = file.path(plot_dir, "tf_genotype")), plot = TRUE)

# STRING with the duplicate-mapping fix (re-query + plot an old object too)
g_old <- g
chk(S, "plot_omni_string_network", "string_network(object saved before fix)",
    plot_omni_string_network(g_old, plot_type = "static"), plot = TRUE)
g <- omni_set_contrast(g, "AZA vs Mock")
g <- keep(chk(S, "run_omni_string", "run_omni_string(top 60, fixed)",
              run_omni_string(g, n = 60, cache_dir = "workspace/data/raw/string_cache")), g)
chk(S, "(check)", "STRING nodes unique",
    stopifnot(!anyDuplicated(g$string_network$nodes$Gene)))
chk(S, "plot_omni_string_network", "plot_omni_string_network(static)",
    plot_omni_string_network(g, plot_type = "static"), plot = TRUE)
chk(S, "plot_omni_string_network", "plot_omni_string_network(interactive)",
    plot_omni_string_network(g, plot_type = "interactive"), widget = TRUE)

# TE coordinates -> proximity -> GREAT
coords <- chk(S, "omni_build_te_coords", "omni_build_te_coords(Telescope GTF)",
              omni_build_te_coords(file.path(proc, "ovarian_te_loci.gtf")))
chk(S, "(check)", "TE coords cover object loci",
    sprintf("%d / %d", sum(rownames(te$data) %in% coords$Locus_ID), nrow(te$data)))
te <- keep(chk(S, "run_omni_te_proximity", "run_omni_te_proximity(up)",
               run_omni_te_proximity(te, te_coords = coords, direction = "up",
                                     lfc_thresh = 1)), te)
te <- keep(chk(S, "run_omni_te_great", "run_omni_te_great(hg38)",
               run_omni_te_great(te, genome = "hg38")), te)
chk(S, "plot_omni_te_great", "plot_omni_te_great()", plot_omni_te_great(te), plot = TRUE)

# Sensors with the metadata fix: plot by any metadata column
gm <- omni_set_contrast(g, "AZA vs Mock")
gs <- chk(S, "run_omni_sensor_score", "run_omni_sensor_score(Telescope families)",
          run_omni_sensor_score(gm, te, group_col = "Treatment",
                                cgas_families = c("L1FLnI", "L1FLI", "HML2", "HERVH", "HERVW"),
                                rigi_families = c("ERV316A3", "ERVLE", "ERVLB4", "MER4", "HERVL"),
                                tlr7_families = c("L1FLI", "L1ORF2", "HERVK", "HML3")))
chk(S, "(check)", "sensor scores keep all metadata columns",
    stopifnot(all(c("Treatment", "CellLine", "AZA", "ITF") %in%
                    names(gs$sensor_scores$scores))))
chk(S, "plot_omni_sensor_irgs_corr", "plot_omni_sensor_irgs_corr(Treatment/CellLine)",
    plot_omni_sensor_irgs_corr(gs, color_col = "Treatment", shape_col = "CellLine"),
    plot = TRUE)
gs <- keep(chk(S, "run_omni_mediation", "run_omni_mediation(AZA, metadata col direct)",
               run_omni_mediation(gs, sensor_scores = gs$sensor_scores$scores,
                                  mediator = "cGAS_STING", treatment_col = "AZA",
                                  treatment_vals = c("0" = "no", "1" = "yes"),
                                  exclude_groups = NULL, boot_sims = 200)), gs)
chk(S, "plot_omni_mediation_path", "plot_omni_mediation_path(cGAS_STING)",
    plot_omni_mediation_path(gs, mediator = "cGAS_STING"), plot = TRUE)

# Single-cell: Bisque now installed; sc_profile fix
pan <- readRDS(file.path(proc, "pancreas_bulk.rds"))
sce <- readRDS(file.path(proc, "pancreas_sc_xin.rds"))
p   <- readRDS(file.path(res_dir, "pancreas_obj.rds"))
chk(S, "plot_omni_sc_profile", "plot_omni_sc_profile(T2D)",
    plot_omni_sc_profile(p, genotype_col = "T2D"), plot = TRUE)
pb <- chk(S, "run_omni_sc_deconv", "run_omni_sc_deconv(bisque)",
          run_omni_sc_deconv(p, sce, cell_type_col = "cell_type", sample_col = "donor",
                             group_col = "group", ref_group = "Control",
                             method = "bisque"))
chk(S, "(check)", "bisque stored results",
    stopifnot(!is.null(pb$sc_deconv$bisque_proportions)))
chk(S, "plot_omni_sc_benchmark", "plot_omni_sc_benchmark(with bisque)",
    plot_omni_sc_benchmark(pb), plot = TRUE)

write.csv(results, file.path(res_dir, "function_results_followup.csv"), row.names = FALSE)
cat("\n==== FOLLOW-UP SUMMARY ====\n"); print(table(results$status))

# plot_omni_deconv_heatmap()/_corr() only need a cell-type x sample matrix in
# obj$deconv; immunedeconv/MuSiC are unavailable here, so feed them the real
# BisqueRNA islet proportions computed above.
S <- "followup-deconv-plots"
pd <- pb
pd$deconv <- pb$sc_deconv$bisque_proportions
pd$params$deconv_method <- "BisqueRNA (via run_omni_sc_deconv)"
chk(S, "plot_omni_deconv_heatmap", "plot_omni_deconv_heatmap(no aggregation)",
    plot_omni_deconv_heatmap(pd, aggregate_main = FALSE), plot = TRUE)
chk(S, "plot_omni_deconv_heatmap", "plot_omni_deconv_heatmap(aggregate_main)",
    plot_omni_deconv_heatmap(pd), plot = TRUE)
chk(S, "plot_omni_deconv_corr", "plot_omni_deconv_corr(INS)",
    plot_omni_deconv_corr(pd, feature = "INS"), plot = TRUE)
chk(S, "plot_omni_deconv_corr", "plot_omni_deconv_corr(GCG)",
    plot_omni_deconv_corr(pd, feature = "GCG"), plot = TRUE)
write.csv(results, file.path(res_dir, "function_results_followup.csv"), row.names = FALSE)
cat("\n==== FOLLOW-UP SUMMARY (final) ====\n"); print(table(results$status))
