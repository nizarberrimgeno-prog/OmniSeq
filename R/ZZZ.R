# =============================================================================
# OmniSeq — zzz.R
#
# Package startup message and dependency installation helper.
# Named zzz.R so it is sourced last, after all other R files.
# =============================================================================


# =============================================================================
# STARTUP MESSAGE
# =============================================================================

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n",
    "  ╔══════════════════════════════════════════════════╗\n",
    "  ║              OmniSeq v", utils::packageVersion("OmniSeq"),
    "                         ║\n",
    "  ║  Multi-omics RNA-seq + TE integration pipeline   ║\n",
    "  ╚══════════════════════════════════════════════════╝\n",
    "\n",
    "  First time? Run:  omni_install_deps()\n",
    "  List signatures:  omni_list_signatures()\n",
    "  Full pipeline:    ?run_omniseq_auto\n",
    "  Guided pipeline:  ?run_omniseq_guided\n"
  )
}


# =============================================================================
# DEPENDENCY INSTALLATION HELPER
# =============================================================================

#' Install All OmniSeq Optional Dependencies
#'
#' Installs all packages required by OmniSeq modules that are not part of
#' the hard `Imports` (those are installed automatically with OmniSeq).
#' Packages are installed from CRAN, Bioconductor, or GitHub as appropriate.
#'
#' You only need to run this once after installing OmniSeq. Safe to re-run —
#' already-installed packages are skipped.
#'
#' @param modules Character vector. Which module groups to install.
#'   Options: `"all"`, `"core"`, `"enrichment"`, `"network"`,
#'   `"multiomics"`, `"te"`, `"export"`. Default: `"all"`.
#' @param ask Logical. Ask before installing each group. Default: `FALSE`.
#' @param bioc_version Character. Bioconductor version to use. Default: `NULL`
#'   (uses currently installed version).
#'
#' @return Invisibly returns a data frame summarising installation status.
#' @export
omni_install_deps <- function(modules    = "all",
                              ask        = FALSE,
                              bioc_version = NULL) {

  # ------------------------------------------------------------------
  # Package registry — grouped by module and source
  # ------------------------------------------------------------------
  deps <- list(

    core = list(
      cran  = c("pheatmap", "ggpubr", "patchwork", "stringr",
                "scales", "dplyr", "ggrepel"),
      bioc  = c("DESeq2", "SummarizedExperiment", "limma",
                "apeglm", "ashr"),
      github = character(0)
    ),

    enrichment = list(
      cran  = c("msigdbr"),
      bioc  = c("clusterProfiler", "enrichplot",
                "org.Hs.eg.db", "org.Mm.eg.db", "GSVA"),
      github = character(0)
    ),

    network = list(
      cran  = c("igraph", "ggraph", "visNetwork", "WGCNA"),
      bioc  = c("STRINGdb"),
      github = character(0)
    ),

    multiomics = list(
      cran  = c("mediation"),
      bioc  = c("mixOmics", "decoupleR", "immunedeconv"),
      github = c("dviraran/xCell")
    ),

    te = list(
      cran  = c(),
      bioc  = c("rGREAT", "rtracklayer", "GenomicRanges", "IRanges"),
      github = character(0)
    ),

    export = list(
      cran  = c("rmarkdown", "knitr"),
      bioc  = character(0),
      github = character(0)
    )
  )

  # Resolve module selection
  if ("all" %in% modules) {
    selected <- names(deps)
  } else {
    selected <- intersect(modules, names(deps))
    if (length(selected) == 0)
      stop("No valid module names. Choose from: all, ",
           paste(names(deps), collapse = ", "), call. = FALSE)
  }

  # ------------------------------------------------------------------
  # Install helper
  # ------------------------------------------------------------------
  status_rows <- list()

  .install_pkg <- function(pkg, source) {
    is_installed <- requireNamespace(pkg, quietly = TRUE)
    if (is_installed) {
      message("  [skip]    ", pkg, " (already installed)")
      return(data.frame(Package = pkg, Source = source,
                        Status = "already installed",
                        stringsAsFactors = FALSE))
    }

    message("  [install] ", pkg, " from ", source, "...")
    result <- tryCatch({
      if (source == "CRAN") {
        utils::install.packages(pkg, quiet = TRUE)
      } else if (source == "Bioconductor") {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          utils::install.packages("BiocManager", quiet = TRUE)
        BiocManager::install(pkg, ask = FALSE, quiet = TRUE,
                             version = bioc_version,
                             update = FALSE)
      } else if (source == "GitHub") {
        if (!requireNamespace("remotes", quietly = TRUE))
          utils::install.packages("remotes", quiet = TRUE)
        remotes::install_github(pkg, quiet = TRUE, upgrade = "never")
      }
      "installed"
    }, error = function(e) {
      paste("FAILED:", e$message)
    })

    data.frame(Package = pkg, Source = source,
               Status = result, stringsAsFactors = FALSE)
  }

  # ------------------------------------------------------------------
  # Run installation per module group
  # ------------------------------------------------------------------
  for (mod in selected) {
    group <- deps[[mod]]

    if (ask) {
      cat("\n")
      ans <- utils::menu(c("Install", "Skip"),
                         title = paste0("Install '", mod, "' dependencies?"))
      if (ans != 1L) {
        message("[omni_install_deps] Skipping module: ", mod)
        next
      }
    }

    message("\n[omni_install_deps] Module: ", toupper(mod))

    for (pkg in group$cran)
      status_rows[[length(status_rows) + 1]] <- .install_pkg(pkg, "CRAN")
    for (pkg in group$bioc)
      status_rows[[length(status_rows) + 1]] <- .install_pkg(pkg, "Bioconductor")
    for (pkg in group$github)
      status_rows[[length(status_rows) + 1]] <- .install_pkg(pkg, "GitHub")
  }

  # ------------------------------------------------------------------
  # Summary table
  # ------------------------------------------------------------------
  status_df <- do.call(rbind, status_rows)

  n_installed <- sum(status_df$Status == "installed")
  n_skipped   <- sum(status_df$Status == "already installed")
  n_failed    <- sum(startsWith(status_df$Status, "FAILED"))

  message("\n[omni_install_deps] Done.")
  message("  Installed : ", n_installed)
  message("  Already   : ", n_skipped)
  if (n_failed > 0) {
    message("  Failed    : ", n_failed)
    failed_pkgs <- status_df$Package[startsWith(status_df$Status, "FAILED")]
    message("  Failed packages: ", paste(failed_pkgs, collapse = ", "))
  }

  message("\n  Restart R and run library(OmniSeq) to begin.\n")

  return(invisible(status_df))
}


#' Check Which OmniSeq Dependencies Are Installed
#'
#' Prints a table showing which optional packages are installed and which
#' are missing, grouped by OmniSeq module.
#'
#' @return Invisibly returns a data frame with installation status per package.
#' @export
omni_check_deps <- function() {

  all_deps <- list(
    "Core (always needed)"     = c("DESeq2", "SummarizedExperiment",
                                   "limma", "ggplot2", "ggrepel",
                                   "pheatmap", "patchwork", "dplyr"),
    "Enrichment"               = c("clusterProfiler", "enrichplot",
                                   "org.Hs.eg.db", "msigdbr", "GSVA"),
    "Network & TF"             = c("igraph", "ggraph", "visNetwork",
                                   "WGCNA", "STRINGdb", "decoupleR"),
    "Multi-omics"              = c("mixOmics", "immunedeconv", "xCell",
                                   "mediation", "apeglm", "ashr"),
    "TE Analysis"              = c("rtracklayer", "GenomicRanges",
                                   "IRanges", "rGREAT"),
    "Export & Report"          = c("rmarkdown", "knitr")
  )

  cat("\n OmniSeq Dependency Status\n")
  cat(strrep("-", 55), "\n")
  cat(sprintf(" %-30s %-12s\n", "Package", "Status"))
  cat(strrep("-", 55), "\n")

  rows <- list()
  for (group in names(all_deps)) {
    cat("\n ", group, "\n")
    for (pkg in all_deps[[group]]) {
      installed <- requireNamespace(pkg, quietly = TRUE)
      status    <- if (installed) "installed" else "MISSING"
      symbol    <- if (installed) "\u2713" else "\u2717"
      cat(sprintf("   %s  %-28s %s\n", symbol, pkg, status))
      rows[[length(rows) + 1]] <- data.frame(
        Group     = group,
        Package   = pkg,
        Installed = installed,
        stringsAsFactors = FALSE
      )
    }
  }

  cat("\n", strrep("-", 55), "\n")
  result <- do.call(rbind, rows)
  n_miss <- sum(!result$Installed)
  if (n_miss > 0) {
    cat(" ", n_miss, "package(s) missing. Run omni_install_deps() to install.\n")
  } else {
    cat("  All dependencies installed.\n")
  }
  cat("\n")

  return(invisible(result))
}
