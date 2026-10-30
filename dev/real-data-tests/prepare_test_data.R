# =============================================================================
# Build OmniSeq-ready test datasets from public data (run from repo root).
#
# 1. GSE182119 — 4 ovarian cancer cell lines (A2780, Hey, Kuramochi, TykNu)
#    x {Mock, ITF-2357 (HDACi), 5-azacytidine (DNMTi), AZA+ITF}.
#    - Gene counts: NCBI-computed raw counts (Entrez IDs -> symbols).
#    - TE counts:   per-sample Telescope locus reports (L1 + HERV loci),
#                   renamed to OmniSeq's "<id>:<Class>:<Family>" format.
#    DNMT inhibition de-represses TEs and triggers an interferon "viral
#    mimicry" response, which suits the TE, sensor and mediation modules.
#    The 2x2 drug design (blocked by cell line) gives a real AZA x ITF
#    interaction term for the interaction plots.
#
# 2. Pancreatic islets (MuSiC tutorial data):
#    - Bulk: GSE50244 (89 donors, HbA1c-defined T2D vs normal).
#    - scRNA-seq: Xin et al. (healthy + T2D donors) for MuSiC / disease state.
# =============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(Biobase)
})

raw <- "workspace/data/raw"
out <- "workspace/data/processed"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- 1a. Ovarian: sample metadata -------------------------------------------
te_files <- list.files(file.path(raw, "GSE182119"),
                       pattern = "telescope_report\\.tsv\\.gz$", full.names = TRUE)
gsm   <- sub("_.*", "", basename(te_files))
line  <- sub(".*RNAseq\\.([A-Za-z0-9]+)\\.library.*", "\\1", basename(te_files))
treat <- sub(".*treat[0-9]\\.(.*)-telescope.*", "\\1", basename(te_files))
treat <- c(mock = "Mock", ITF = "ITF", AZA = "AZA", AZA.ITF = "AZA_ITF")[treat]

meta <- data.frame(
  CellLine  = factor(line),
  Treatment = factor(treat, levels = c("Mock", "ITF", "AZA", "AZA_ITF")),
  AZA       = factor(ifelse(grepl("AZA", treat), "yes", "no"), levels = c("no", "yes")),
  ITF       = factor(ifelse(grepl("ITF", treat), "yes", "no"), levels = c("no", "yes")),
  row.names = gsm
)

# ---- 1b. Ovarian: gene counts -----------------------------------------------
g <- utils::read.delim(file.path(raw, "GSE182119", "ncbi_counts.tsv.gz"),
                       check.names = FALSE)
sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, as.character(g$GeneID),
                             "SYMBOL", "ENTREZID")
keep <- !is.na(sym)
gene_counts <- rowsum(as.matrix(g[keep, rownames(meta)]), sym[keep])
gene_counts <- gene_counts[rowSums(gene_counts >= 10) >= 3, ]

# ---- 1c. Ovarian: TE locus counts -------------------------------------------
reports <- lapply(te_files, function(f) {
  d <- utils::read.delim(f, skip = 1, check.names = FALSE)
  stats::setNames(d$final_count, d$transcript)
})
loci <- setdiff(Reduce(union, lapply(reports, names)), "__no_feature")
te_counts <- sapply(reports, function(x) { v <- x[loci]; v[is.na(v)] <- 0; v })
dimnames(te_counts) <- list(loci, gsm)
te_counts <- te_counts[rowSums(te_counts >= 5) >= 3, ]

family  <- sub("_.*", "", rownames(te_counts))
te_cls  <- ifelse(grepl("^L1", family), "LINE", "ERV")
te_ids  <- paste0(1000000 + seq_along(family), ":", te_cls, ":", family)
te_map  <- data.frame(Locus_ID = te_ids, Telescope_ID = rownames(te_counts),
                      Class = te_cls, Family = family)
rownames(te_counts) <- te_ids

# ---- 2. Pancreas bulk + scRNA reference --------------------------------------
eset  <- readRDS(file.path(raw, "music", "GSE50244bulkeset.rds"))
pmeta <- pData(eset)
pmeta$T2D <- factor(ifelse(!is.na(pmeta$hba1c) & pmeta$hba1c >= 6.5, "T2D",
                           "Normal"), levels = c("Normal", "T2D"))
keep_s <- !is.na(pmeta$hba1c)
pan_counts <- exprs(eset)[, keep_s]
pan_meta   <- droplevels(pmeta[keep_s, , drop = FALSE])
pan_counts <- pan_counts[rowSums(pan_counts >= 10) >= 5, ]

sc <- readRDS(file.path(raw, "music", "XinT2Dsce.rds"))
sc$cell_type <- as.character(sc$cellType)
sc$donor     <- as.character(sc$SubjectName)
sc$group     <- ifelse(sc$Disease == "Healthy", "Control", "T2D")

# ---- Save --------------------------------------------------------------------
saveRDS(list(gene_counts = gene_counts, te_counts = te_counts, meta = meta,
             te_map = te_map), file.path(out, "ovarian_gene_te.rds"))
saveRDS(list(counts = pan_counts, meta = pan_meta), file.path(out, "pancreas_bulk.rds"))
saveRDS(sc, file.path(out, "pancreas_sc_xin.rds"))

write.csv(cbind(Gene = rownames(gene_counts), gene_counts),
          file.path(out, "ovarian_gene_counts.csv"), row.names = FALSE)
write.csv(cbind(Locus = rownames(te_counts), te_counts),
          file.path(out, "ovarian_te_counts.csv"), row.names = FALSE)
write.csv(cbind(Sample = rownames(meta), meta),
          file.path(out, "ovarian_metadata.csv"), row.names = FALSE)

cat("Ovarian genes:", dim(gene_counts), "| TE loci:", dim(te_counts), "\n")
print(table(meta$CellLine, meta$Treatment))
cat("Pancreas bulk:", dim(pan_counts), "\n"); print(table(pan_meta$T2D))
cat("Xin sc:", dim(sc), "\n")

# ---- TE coordinates: Telescope retro.hg38.v1 GTF with OmniSeq locus IDs -----
# (source: github.com/mlbendall/telescope_annotation_db, builds/retro.hg38.v1)
gtf <- readLines(file.path(raw, "retro.hg38.gtf"))
tid <- sub('.*transcript_id "([^"]+)".*', "\\1", gtf)
hit <- tid %in% te_map$Telescope_ID
new <- te_map$Locus_ID[match(tid[hit], te_map$Telescope_ID)]
gtf <- gtf[hit]
gtf <- vapply(seq_along(gtf), function(i)
  sub('transcript_id "[^"]+"', paste0('transcript_id "', new[i], '"'), gtf[i],
      fixed = FALSE), character(1))
writeLines(gtf, file.path(out, "ovarian_te_loci.gtf"))
cat("TE GTF rows written:", length(gtf), "for", length(unique(new)), "loci\n")
