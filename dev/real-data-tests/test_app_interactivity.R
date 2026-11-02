# Render every plot output of the Shiny app on the real test objects and
# report whether it comes out as an interactive plotly chart, an interactive
# network (visNetwork) or a static image. Run from the repo root after
# test_all_functions.R.
suppressPackageStartupMessages({ library(shiny); library(OmniSeq) })
for (f in list.files("inst/shiny-app/R", full.names = TRUE)) source(f)

res_dir <- "workspace/test_results"
gene_rds <- file.path(res_dir, "ovarian_gene_obj.rds")
te_rds   <- file.path(res_dir, "ovarian_te_obj.rds")
out <- data.frame(tab = character(), output = character(), kind = character())

kind_of <- function(x) {
  html <- if (is.list(x)) paste(unlist(x$html), collapse = "") else paste(x, collapse = "")
  if (grepl("visNetwork", html)) "interactive network"
  else if (grepl("plotly", html)) "interactive plotly"
  else if (grepl("<img", html)) "static image"
  else if (grepl("<table|datatables", html)) "table"
  else substr(gsub("\\s+", " ", html), 1, 80)
}
check <- function(session, tab, ids, label = NULL) for (id in ids) {
  k <- tryCatch(kind_of(session$getOutput(id)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  lab <- if (is.null(label)) id else paste0(id, " (", label, ")")
  out[nrow(out) + 1, ] <<- list(tab, lab, k)
  cat(sprintf("%-14s %-32s %s\n", tab, lab, k))
}

holder <- NULL
testServer(function(input, output, session) {
  state <- reactiveValues(obj = NULL, partner = NULL, log = character(0),
                          ora_updown = list())
  holder <<- state
  mod_data_server("data", state); mod_qc_server("qc", state)
  mod_de_server("de", state); mod_enrichment_server("enrich", state)
  mod_compare_server("compare", state); mod_biology_server("bio", state)
  mod_multivariate_server("multi", state); mod_te_server("te", state)
  mod_multiomics_server("omics", state)
}, {
  st <- holder
  session$setInputs(`data-source` = "rds", `data-round_counts` = FALSE)
  st$obj <- readRDS(gene_rds)
  st$obj <- omni_set_contrast(st$obj, "Treatment_AZA_ITF_vs_Mock")
  session$flushReact()

  check(session, "QC", paste0("qc-", c("libsize", "detected", "distance", "mean_var", "pca")))

  session$setInputs(`de-lfc_thresh` = 1, `de-pval_thresh` = "0.05", `de-top_n` = 10,
                    `de-heat_n` = 40, `de-heat_genes` = "", `de-color_by` = "Treatment",
                    `de-shape_by` = "CellLine")
  check(session, "DE", paste0("de-", c("volcano", "ma", "heatmap", "pca", "top_genes", "phist")))

  session$setInputs(`enrich-gsea_plot` = "waterfall", `enrich-gsea_n` = 15,
                    `enrich-ss_n` = 25, `enrich-ss_isolate` = TRUE)
  check(session, "Enrichment", "enrich-gsea_plot_out")
  for (pt in c("ridge", "bar")) { session$setInputs(`enrich-gsea_plot` = pt)
    check(session, "Enrichment", "enrich-gsea_plot_out", pt) }
  check(session, "Enrichment", "enrich-ss_plot")

  st$obj <- omni_set_contrast(st$obj, "AZA vs Mock"); session$flushReact()
  session$setInputs(`enrich-ora_plot` = "dot", `enrich-ora_n` = 15, `enrich-ora_ont` = "BP")
  check(session, "Enrichment", "enrich-ora_plot_out")
  for (pt in c("cnet", "emap")) { session$setInputs(`enrich-ora_plot` = pt)
    check(session, "Enrichment", "enrich-ora_plot_out", pt) }

  st$obj <- omni_set_contrast(st$obj, "Treatment_AZA_ITF_vs_Mock")
  session$setInputs(`enrich-ud_db` = c("GO-BP", "Hallmark"), `enrich-ud_lfc` = 1,
                    `enrich-ud_pval` = "0.05", `enrich-ud_show_db` = "GO-BP",
                    `enrich-ud_n` = 15, `enrich-run_ud` = 1)
  st$obj <- omni_set_contrast(st$obj, "AZA_ITF vs ITF")
  session$setInputs(`enrich-run_ud` = 2)
  session$setInputs(`enrich-ud_cmp` = names(st$ora_updown))
  check(session, "Enrichment", c("enrich-ud_plot", "enrich-ud_compare"))
  session$setInputs(`enrich-cmp_contrasts` = c("Treatment_AZA_ITF_vs_Mock", "AZA_ITF vs ITF"),
                    `enrich-cmp_lfc` = 1, `enrich-cmp_pval` = "0.05", `enrich-cmp_n` = 8,
                    `enrich-run_cmp` = 1)
  check(session, "Enrichment", "enrich-cmp_plot")

  session$setInputs(`compare-a` = "AZA_ITF vs ITF", `compare-b` = "AZA vs Mock",
                    `compare-label_a` = "", `compare-label_b` = "", `compare-pval` = "0.05",
                    `compare-lfc` = 1, `compare-inter` = "AZA x ITF interaction",
                    `compare-flip` = 0.3, `compare-stim` = "AZA", `compare-rank_by` = "padj",
                    `compare-n_top` = 20)
  check(session, "Compare", paste0("compare-", c("concord", "inter_int", "inter_static", "dumbbell")))

  st$obj <- omni_set_contrast(st$obj, "AZA vs Mock"); session$flushReact()
  session$setInputs(`bio-tf_mode` = "contrast", `bio-tf_n` = 25, `bio-tf_min_reg` = 15,
                    `bio-tf_targets` = 5, `bio-string_degree` = 1, `bio-string_hubs` = TRUE)
  check(session, "Biology", paste0("bio-", c("tf_activity", "tf_network", "string_plot")))

  session$setInputs(`multi-w_show` = "full", `multi-w_traits` = NULL, `multi-w_hubs` = 5,
                    `multi-s_key` = names(st$obj$splsda)[1], `multi-s_comp` = 1,
                    `multi-s_nfeat` = 15, `multi-mp_col` = "AZA")
  check(session, "Multivariate", paste0("multi-", c("w_dendro", "w_trait", "w_hub", "mp_plot",
                                                    "s_scores", "s_load", "s_vip", "s_roc")))

  st$partner <- readRDS(te_rds); session$flushReact()
  session$setInputs(`omics-sens_col` = "Treatment", `omics-sens_excl` = NULL, `omics-sens_min` = 5,
                    `omics-fam_cgas` = "L1FLnI, L1FLI, HML2, HERVH, HERVW",
                    `omics-fam_rigi` = "ERV316A3, ERVLE, ERVLB4, MER4, HERVL",
                    `omics-fam_tlr7` = "L1FLI, L1ORF2, HERVK, HML3", `omics-run_sens` = 1,
                    `omics-irg_color` = "Treatment", `omics-irg_shape` = "CellLine",
                    `omics-corr_genes` = "IFI27\nIFI44L\nISG15\nMX1\nRSAD2\nOAS1",
                    `omics-corr_class` = "ALL", `omics-corr_name` = "ISG score")
  check(session, "Multi-omics", paste0("omics-", c("sens_plot", "sens_heat", "sens_irg", "corr_plot")))

  st$obj <- omni_set_contrast(readRDS(te_rds), "AZA vs Mock"); session$flushReact()
  session$setInputs(`te-lfc` = 1, `te-pval` = "0.05", `te-n_fam` = 20, `te-top_n` = 10,
                    `te-label_col` = "family", `te-color_by` = "class", `te-min_loci` = 20)
  check(session, "TE", paste0("te-", c("volcano", "family", "pie", "density", "heatmap", "great")))
})

write.csv(out, file.path(res_dir, "app_interactivity.csv"), row.names = FALSE)
cat("\n"); print(table(out$kind))
