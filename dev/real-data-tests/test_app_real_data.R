# Drive the Shiny app's modules headlessly on the real ovarian object
# produced by test_all_functions.R (run from repo root).
suppressPackageStartupMessages({library(shiny); library(OmniSeq)})
app_dir <- "inst/shiny-app"
for (f in list.files(file.path(app_dir, "R"), full.names = TRUE)) source(f)

ok <- function(label, expr) {
  r <- tryCatch({force(expr); "OK"}, error = function(e) paste("ERR:", conditionMessage(e)))
  cat(sprintf("%-34s %s\n", label, r))
}
rds <- "workspace/test_results/ovarian_gene_obj.rds"
te_csv <- "workspace/data/processed/ovarian_te_counts.csv"
holder <- NULL

testServer(function(input, output, session) {
  state <- reactiveValues(obj = NULL, partner = NULL, log = character(0),
                          ora_updown = list())
  holder <<- state
  mod_data_server("data", state); mod_de_server("de", state)
  mod_enrichment_server("enrich", state); mod_compare_server("compare", state)
  mod_te_server("te", state); mod_multiomics_server("omics", state)
}, {
  st <- holder
  session$setInputs(`data-source` = "rds",
                    `data-rds_file` = list(datapath = rds, name = basename(rds)))
  session$setInputs(`data-load_rds` = 1)
  cat("loaded contrasts:", paste(omni_contrasts(st$obj), collapse = " | "), "\n")
  st$obj <- omni_set_contrast(st$obj, "AZA vs Mock")

  session$setInputs(`de-lfc_thresh` = 1, `de-pval_thresh` = "0.05", `de-top_n` = 10,
                    `de-heat_n` = 40, `de-heat_genes` = "", `de-color_by` = "Treatment",
                    `de-shape_by` = "CellLine", `de-sig_only` = TRUE)
  for (o in c("volcano", "heatmap", "top_genes", "table"))
    ok(paste("de", o), session$getOutput(paste0("de-", o)))

  st$obj <- omni_set_contrast(st$obj, "Treatment_AZA_ITF_vs_Mock")
  session$setInputs(`enrich-ud_db` = c("GO-BP", "Hallmark"), `enrich-ud_lfc` = 1,
                    `enrich-ud_pval` = "0.05", `enrich-ud_show_db` = "Hallmark",
                    `enrich-ud_n` = 15, `enrich-run_ud` = 1)
  ok("ora up/down plot", session$getOutput("enrich-ud_plot")); ok("ora up/down table", session$getOutput("enrich-ud_table"))
  st$obj <- omni_set_contrast(st$obj, "AZA_ITF vs ITF")
  session$setInputs(`enrich-run_ud` = 2)
  session$setInputs(`enrich-ud_show_db` = "GO-BP")
  session$setInputs(`enrich-ud_cmp` = names(st$ora_updown))
  cat("ora_updown contrasts:", paste(names(st$ora_updown), collapse = " | "), "\n")
  ok("ora up/down comparison", session$getOutput("enrich-ud_compare"))

  session$setInputs(`compare-a` = "AZA_ITF vs ITF", `compare-b` = "AZA vs Mock",
                    `compare-label_a` = "AZA | ITF", `compare-label_b` = "AZA alone",
                    `compare-pval` = "0.05", `compare-lfc` = 1,
                    `compare-inter` = "AZA x ITF interaction", `compare-flip` = 0.3,
                    `compare-stim` = "AZA", `compare-rank_by` = "padj", `compare-n_top` = 20)
  for (o in c("concord", "inter_int", "inter_static", "dumbbell"))
    ok(paste("compare", o), session$getOutput(paste0("compare-", o)))
  session$setInputs(`compare-rank_by` = "gap")
  ok("compare gap warning", stopifnot(!is.null(session$getOutput("compare-gap_warning"))))

  # Multi-omics: pair with TE counts and score sensors with custom families
  session$setInputs(`omics-pair_src` = "counts",
                    `omics-pair_counts` = list(datapath = te_csv, name = basename(te_csv)))
  session$setInputs(`omics-load_pair` = 1)
  cat("partner type:", st$partner$params$data_type, "\n")
  ok("omics family hint", session$getOutput("omics-fam_hint"))
  session$setInputs(`omics-sens_col` = "Treatment", `omics-sens_excl` = NULL,
                    `omics-sens_min` = 5, `omics-fam_cgas` = "L1FLnI, HML2, HERVH",
                    `omics-fam_rigi` = "ERV316A3, ERVLE", `omics-fam_tlr7` = "L1FLI, L1ORF2",
                    `omics-run_sens` = 1)
  ok("omics sensor plot", session$getOutput("omics-sens_plot"))
  cat("\nlog errors:\n"); cat(grep("Error", st$log, value = TRUE), sep = "\n")
})
