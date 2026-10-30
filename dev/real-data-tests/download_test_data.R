# =============================================================================
# Download the public datasets used by the real-data function tests.
# Run from the repo root:  Rscript dev/real-data-tests/download_test_data.R
# Files land in workspace/data/raw/ (gitignored). ~110 MB in total.
# =============================================================================

raw <- "workspace/data/raw"
dir.create(file.path(raw, "GSE182119"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(raw, "music"), recursive = TRUE, showWarnings = FALSE)
options(timeout = 1200)

get <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 0) return(invisible(dest))
  message("Downloading ", basename(dest))
  utils::download.file(url, dest, mode = "wb", quiet = TRUE)
}

# 1. GSE182119: ovarian cancer lines x {Mock, ITF-2357, 5-AZA, AZA+ITF}
#    Per-sample Telescope TE locus reports (series RAW tar) ...
tar <- file.path(raw, "GSE182119", "GSE182119_RAW.tar")
get("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182119/suppl/GSE182119_RAW.tar", tar)
utils::untar(tar, exdir = file.path(raw, "GSE182119"))
#    ... and NCBI-computed raw gene counts for the same 16 samples.
get(paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE182119",
           "&format=file&file=GSE182119_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
    file.path(raw, "GSE182119", "ncbi_counts.tsv.gz"))

# 2. Telescope retro.hg38.v1 annotation (TE locus coordinates; Git LFS file)
get(paste0("https://media.githubusercontent.com/media/mlbendall/telescope_annotation_db/",
           "master/builds/retro.hg38.v1/transcripts.gtf"),
    file.path(raw, "retro.hg38.gtf"))

# 3. Pancreatic islets (MuSiC tutorial data): GSE50244 bulk + Xin 2016 scRNA
for (f in c("GSE50244bulkeset.rds", "XinT2Dsce.rds", "EMTABsce_healthy.rds"))
  get(paste0("https://xuranw.github.io/MuSiC/data/", f), file.path(raw, "music", f))

# 4. CollecTRI TF regulons straight from the OmniPath web service (used to
#    seed run_omni_tf()'s regulon cache where OmnipathR cannot run)
get(paste0("https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes",
           "&organisms=9606&fields=is_stimulation,is_inhibition"),
    file.path(raw, "collectri_human.tsv"))

message("Done. Next: Rscript dev/real-data-tests/prepare_test_data.R")
