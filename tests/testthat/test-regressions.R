# Regression tests for bugs found by running every exported function on
# public data (workspace/scripts/test_all_functions.R).

make_obj <- function(meta) {
  counts <- matrix(10L, nrow = 3, ncol = nrow(meta),
                   dimnames = list(c("G1", "G2", "G3"), rownames(meta)))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), colData = meta)
  new_omniseq(data = se, data_type = "gene", condition = names(meta)[1])
}

test_that(".omni_collapse_gene_edges keeps one edge per gene pair", {
  edges <- data.frame(from_gene = c("A", "B", "A", "C"),
                      to_gene   = c("B", "A", "B", "A"),
                      combined_score = c(500, 900, 700, 400))
  out <- .omni_collapse_gene_edges(edges)
  expect_equal(nrow(out), 2)
  expect_equal(max(out$combined_score[out$from_gene %in% c("A", "B") &
                                        out$to_gene %in% c("A", "B")]), 900)
})

test_that("omni_build_te_coords works without dplyr/magrittr attached", {
  skip_if_not_installed("rtracklayer")
  gtf <- withr::local_tempfile(fileext = ".gtf")
  writeLines(c(
    'chr1\tx\texon\t100\t200\t.\t+\t.\ttranscript_id "1:LINE:L1HS";',
    'chr1\tx\texon\t300\t400\t.\t+\t.\ttranscript_id "1:LINE:L1HS";',
    'chr2\tx\texon\t50\t80\t.\t-\t.\ttranscript_id "2:ERV:HERVK";'), gtf)
  expect_false("package:dplyr" %in% search())
  coords <- suppressMessages(omni_build_te_coords(gtf))
  expect_equal(nrow(coords), 2)
  l1 <- coords[coords$Locus_ID == "1:LINE:L1HS", ]
  expect_equal(c(l1$Start, l1$End), c(100, 400))
  expect_equal(coords$Family[coords$Locus_ID == "2:ERV:HERVK"], "HERVK")
})

test_that("plot_omni_sc_genotype labels a factor reference group", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  meta <- data.frame(Status = factor(rep(c("Normal", "T2D"), each = 3),
                                     levels = c("Normal", "T2D")),
                     row.names = paste0("S", 1:6))
  obj <- make_obj(meta)
  obj$sc_deconv <- list(scores = rbind(alpha = c(0, .1, .2, 1, 1.1, 1.2),
                                       beta  = c(1, 1.1, .9, 0, .1, -.1)))
  colnames(obj$sc_deconv$scores) <- rev(rownames(meta))  # order differs from meta
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  p <- suppressMessages(plot_omni_sc_genotype(obj, genotype_col = "Status"))
  expect_equal(levels(p$data$Genotype), c("Normal", "T2D"))
  expect_false(anyNA(p$data$Genotype))
})

test_that("plot_omni_sc_profile does not fail after drawing", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  meta <- data.frame(Status = factor(rep(c("Normal", "T2D"), each = 3)),
                     row.names = paste0("S", 1:6))
  obj <- make_obj(meta)
  obj$sc_deconv <- list(scores = rbind(alpha = 1:6 / 6, beta = 6:1 / 6))
  colnames(obj$sc_deconv$scores) <- rownames(meta)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  p <- NULL
  utils::capture.output(
    p <- suppressMessages(plot_omni_sc_profile(obj, genotype_col = "Status")))
  expect_s3_class(p, "ggplot")
})

test_that("plot_omni_tf_genotype explains an unnamed geno_contrasts", {
  meta <- data.frame(G = factor(c("a", "a", "b", "b")),
                     row.names = paste0("S", 1:4))
  obj <- make_obj(meta)
  obj$tf_contrast_cache <- list(m = data.frame(), x = data.frame(), y = data.frame())
  expect_error(plot_omni_tf_genotype(obj, main_contrast = "m",
                                     geno_contrasts = list("x", "y")),
               "must be named")
})

test_that("interactive volcano uses the documented Up/Down colours", {
  skip_if_not_installed("plotly")
  skip_if_not_installed("htmltools")
  meta <- data.frame(G = factor(rep(c("a", "b"), each = 2)),
                     row.names = paste0("S", 1:4))
  obj <- make_obj(meta)
  res <- data.frame(log2FoldChange = c(3, -3, 0.1), padj = c(1e-5, 1e-5, 0.9),
                    row.names = c("G1", "G2", "G3"))
  obj$results_cache[["b vs a"]] <- list(results = res)
  obj$active_contrast <- "b vs a"
  built <- plotly::plotly_build(plot_omni_volcano_interactive(obj))
  cols  <- vapply(built$x$data, function(t) t$marker$color, character(1))
  names(cols) <- vapply(built$x$data, function(t) sub(" \\(n=.*", "", t$name),
                        character(1))
  expect_equal(unname(cols[c("Up", "Down")]), c("#E74C3C", "#2980B9"))
})

test_that("plot_omni_deconv_heatmap keeps non-immune cell types", {
  skip_if_not_installed("pheatmap")
  meta <- data.frame(G = factor(rep(c("a", "b"), each = 3)),
                     row.names = paste0("S", 1:6))
  obj <- make_obj(meta)
  obj$deconv <- rbind(alpha = c(.5, .4, .6, .3, .2, .3),
                      beta  = c(.4, .5, .3, .6, .7, .6),
                      Monocyte = rep(.1, 6))
  colnames(obj$deconv) <- rownames(meta)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  p <- plot_omni_deconv_heatmap(obj)
  expect_setequal(p$tree_row$labels, c("Monocyte", "alpha", "beta"))
})
