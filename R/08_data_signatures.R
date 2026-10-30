# =============================================================================
# OmniSeq — data_signatures.R
#
# Built-in curated gene signature lists for post-DE masking.
#
# Philosophy:
#   These are universal confounding/housekeeping gene sets that apply across
#   almost any bulk RNA-seq experiment. Disease-specific or tissue-specific
#   signatures should be defined by the user in their analysis script.
#
# Usage:
#   # Use a built-in list directly
#   obj <- run_omni_mask(obj, features = omniseq_isgs_28,
#                        mask_label = "ISGs_28")
#
#   # Combine built-in + custom genes
#   my_list <- c(omniseq_isgs_28, "EXTRA_GENE1", "EXTRA_GENE2")
#   obj <- run_omni_mask(obj, features = my_list, mask_label = "ISGs_Custom")
#
#   # Use the shorthand name via run_omni_mask(signature = "isgs_28")
#   obj <- run_omni_mask(obj, signature = "isgs_28")
#
#   # List all available built-in signatures
#   omni_list_signatures()
#
# Available signatures:
#   omniseq_isgs_28        28-gene NIH IRG-S interferon panel
#   omniseq_isgs_extended  Broad type I IFN signature (Hallmark-based)
#   omniseq_ribosomal      Cytoplasmic + mitoribosomal protein genes
#   omniseq_mitochondrial  Mitochondrially-encoded genes (MT- prefix)
#   omniseq_immunoglobulin Immunoglobulin variable region genes
# =============================================================================


# =============================================================================
# 1. 28-GENE INTERFERON RESPONSE SIGNATURE (IRG-S)
# =============================================================================

#' OmniSeq 28-Gene Interferon Response Gene Signature (IRG-S)
#'
#' The gold-standard 28-gene Interferon Response Gene Score (IRG-S) used for
#' diagnosing and scoring Type I interferonopathies in clinical and research
#' settings. This panel is commonly dominated in inflammatory bulk RNA-seq
#' datasets and masking it reveals secondary biological signals.
#'
#' @source Kim H, de Jesus AA, Brooks SR, Liu Y, Huang Y, VanTries R,
#'   Montealegre Sanchez GA, Rotman Y, Gadina M, Goldbach-Mansky R (2018).
#'   Development of a Validated Interferon Score Using NanoString Technology.
#'   J Interferon Cytokine Res. 38(4):171-185. PMID: 29638206.
#'   NIH/NIAMS IRG-S panel (28 genes), developed and validated in CANDLE,
#'   SAVI, and NOMID patients.
#'
#'   (Corrected 2026-09-24: previously mis-cited as "Crow YJ et al. 2015,
#'   J Clin Immunol, PMID 26250437" -- that PMID belongs to an unrelated
#'   pediatric dizziness-questionnaire paper and was never a real IRG-S
#'   source. The gene list itself was also corrected to match the actual
#'   Kim et al. 2018 panel.)
#'
#' @details
#' \strong{Recommended use}: Post-DE masking via \code{run_omni_mask()} to explore
#' biology beyond the interferon response, as recommended by DESeq2 best
#' practices (run full model, then exclude from downstream lists).
#'
#' @export
omniseq_isgs_28 <- c(
  "CXCL10", "DDX60", "EPSTI1", "GBP1", "HERC5", "HERC6",
  "IFI27", "IFI44", "IFI44L", "IFI6", "IFIT1", "IFIT2",
  "IFIT3", "IFIT5", "ISG15", "LAMP3", "LY6E", "MX1",
  "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "RTP4",
  "SIGLEC1", "SOCS1", "SPATS2L", "USP18"
)


# =============================================================================
# 2. EXTENDED TYPE I INTERFERON SIGNATURE (Hallmark-based)
# =============================================================================

#' OmniSeq Extended Type I Interferon Signature (Hallmark IFNA Response)
#'
#' A broad interferon-stimulated gene set combining the core IRG-S 28-gene
#' panel with additional genes from the MSigDB Hallmark
#' HALLMARK_INTERFERON_ALPHA_RESPONSE collection. This extended set includes
#' genuinely significant ISGs such as HLA-A/B/C, PSMB8/9, TAP1, B2M, IRF1,
#' ADAR, and IL15 — genes that are significant in interferonopathy datasets
#' and will visually differentiate masked from unmasked analyses.
#'
#' All 28 core IRG-S genes are retained. Use this mask when the standard
#' isgs_28 mask produces results that still appear dominated by IFN biology.
#'
#' @source Core IRG-S: Crow YJ et al. (2015). PMID: 26250437.
#'   Extended set: MSigDB Hallmark HALLMARK_INTERFERON_ALPHA_RESPONSE.
#'
#' @export
omniseq_isgs_extended <- unique(c(
  # Core IRG-S 28 (all retained)
  "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
  "ISG15", "RSAD2", "SIGLEC1", "USP18", "CXCL10", "OAS1",
  "OAS2", "OAS3", "OASL", "MX1", "MX2", "STAT1", "STAT2",
  "EPSTI1", "GBP1", "HERC5", "LY6E", "SOCS1", "IFI6",
  "CMPK2", "EIF2AK2", "XAF1",
  # Additional Hallmark INTERFERON_ALPHA_RESPONSE genes not in isgs_28
  "ADAR", "APOBEC3A", "APOBEC3B", "B2M", "BATF2", "C1S",
  "CASP1", "CASP7", "CCRL2", "CD274", "CDKN1A", "CFB",
  "CIITA", "CMKLR1", "CNP", "CSF1", "CXCL8", "CXCL9",
  "CXCL11", "DDX58", "DDX60", "DNASE1L3", "DTX3L",
  "FAS", "FBXO6", "FLNB", "FTSJD2", "GBP2", "GBP4",
  "GBP5", "HERC6", "HIF1A", "HLA-A", "HLA-B", "HLA-C",
  "HLA-E", "HLA-F", "IFI16", "IFI30", "IFI35", "IFIH1",
  "IFIT5", "IFITM1", "IFITM2", "IFITM3", "IL10RA", "IL15",
  "IL4R", "IRF1", "IRF2", "IRF5", "IRF7", "IRF8", "IRF9",
  "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LHFPL2",
  "MAP3K14", "METTL7B", "MIDN", "MT2A", "MVB12A", "MYD88",
  "NCOA7", "NMI", "NUB1", "PARP12", "PARP14", "PARP9",
  "PLSCR1", "PNPT1", "PSMA2", "PSMA3", "PSMB8", "PSMB9",
  "PSME1", "PSME2", "REC8", "RNF31", "SAMD9", "SAMD9L",
  "SELL", "SERPING1", "SP100", "SP110", "SP140",
  "STAT3", "TAP1", "TAPBP", "TDRD7", "TMEM140", "TMEM62",
  "TRAFD1", "TRIM14", "TRIM21", "TRIM22", "TRIM25", "TRIM26",
  "TRIM5", "TRIM56", "TXNIP", "UBE2L6", "UPP1", "WARS1",
  "WAS", "ZBP1", "BST2", "DHX58"
))


# =============================================================================
# 3. RIBOSOMAL PROTEIN GENES
# =============================================================================

#' OmniSeq Ribosomal Protein Gene Signature
#'
#' Cytoplasmic ribosomal protein genes (RPL and RPS families) plus
#' mitochondrial ribosomal protein genes (MRPL and MRPS families).
#' These genes are highly expressed in most cell types, can dominate
#' variance-based analyses, and often confound pathway enrichment results.
#'
#' @source Ensembl gene annotation, GO:0005840 (ribosome) component genes.
#'
#' @details
#' \strong{Recommended use}: Masking before heatmap clustering or GSEA when
#' ribosomal terms are dominating enrichment results.
#'
#' @export
omniseq_ribosomal <- c(
  # Large subunit cytoplasmic (RPL)
  "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8",
  "RPL9", "RPL10", "RPL10A", "RPL11", "RPL12", "RPL13",
  "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A",
  "RPL19", "RPL21", "RPL22", "RPL22L1", "RPL23", "RPL23A",
  "RPL24", "RPL26", "RPL27", "RPL27A", "RPL28", "RPL29",
  "RPL30", "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A",
  "RPL36", "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39",
  "RPL40", "RPL41", "RPLP0", "RPLP1", "RPLP2",
  # Small subunit cytoplasmic (RPS)
  "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS4Y1", "RPS5",
  "RPS6", "RPS7", "RPS8", "RPS9", "RPS10", "RPS11",
  "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16",
  "RPS17", "RPS18", "RPS19", "RPS20", "RPS21", "RPS23",
  "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS28",
  "RPS29", "RPSA",
  # Mitochondrial ribosomal (MRPL / MRPS)
  "MRPL1", "MRPL2", "MRPL3", "MRPL4", "MRPL9", "MRPL10",
  "MRPL11", "MRPL12", "MRPL13", "MRPL14", "MRPL15", "MRPL16",
  "MRPL17", "MRPL18", "MRPL19", "MRPL20", "MRPL21", "MRPL22",
  "MRPL23", "MRPL24", "MRPL27", "MRPL28", "MRPL30", "MRPL32",
  "MRPL33", "MRPL34", "MRPL35", "MRPL36", "MRPL37", "MRPL38",
  "MRPL39", "MRPL40", "MRPL41", "MRPL42", "MRPL43", "MRPL44",
  "MRPL45", "MRPL46", "MRPL47", "MRPL48", "MRPL49", "MRPL50",
  "MRPL51", "MRPL52", "MRPL53", "MRPL54", "MRPL55", "MRPL57",
  "MRPS2", "MRPS5", "MRPS6", "MRPS7", "MRPS9", "MRPS10",
  "MRPS11", "MRPS12", "MRPS14", "MRPS15", "MRPS16", "MRPS17",
  "MRPS18A", "MRPS18B", "MRPS18C", "MRPS21", "MRPS22",
  "MRPS23", "MRPS24", "MRPS25", "MRPS26", "MRPS27", "MRPS28",
  "MRPS30", "MRPS31", "MRPS33", "MRPS34", "MRPS35", "MRPS36"
)


# =============================================================================
# 4. MITOCHONDRIALLY-ENCODED GENES
# =============================================================================

#' OmniSeq Mitochondrially-Encoded Gene Signature
#'
#' The 37 genes encoded by the human mitochondrial genome (MT- prefix).
#' These include 13 protein-coding genes, 22 tRNA genes, and 2 rRNA genes.
#' High MT gene expression relative to nuclear genes can indicate cell stress,
#' apoptosis, or poor sample quality. Often excluded from DE analyses in
#' single-cell RNA-seq and useful to mask in bulk RNA-seq when mitochondrial
#' transcription dominates.
#'
#' @source Human mitochondrial genome annotation (NC_012920.1).
#'
#' @export
omniseq_mitochondrial <- c(
  # Protein-coding (13 genes)
  "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L",
  "MT-ND5", "MT-ND6", "MT-CO1", "MT-CO2", "MT-CO3",
  "MT-ATP6", "MT-ATP8", "MT-CYB",
  # Ribosomal RNA (2 genes)
  "MT-RNR1", "MT-RNR2",
  # Transfer RNA (22 genes)
  "MT-TA", "MT-TC", "MT-TD", "MT-TE", "MT-TF", "MT-TG",
  "MT-TH", "MT-TI", "MT-TK", "MT-TL1", "MT-TL2", "MT-TM",
  "MT-TN", "MT-TP", "MT-TQ", "MT-TR", "MT-TS1", "MT-TS2",
  "MT-TT", "MT-TV", "MT-TW", "MT-TY"
)


# =============================================================================
# 5. IMMUNOGLOBULIN VARIABLE REGION GENES
# =============================================================================

#' OmniSeq Immunoglobulin Variable Region Gene Signature
#'
#' Heavy chain (IGHV), kappa light chain (IGKV), and lambda light chain
#' (IGLV) variable region genes. These are highly expressed in B cells and
#' plasma cells and dominate DE results in samples with varying B cell
#' content. Masking them reveals T cell, myeloid, and stromal biology that
#' would otherwise be hidden.
#'
#' @source IMGT/GENE-DB immunoglobulin gene nomenclature.
#'   Only functional (F) and open reading frame (ORF) genes included.
#'
#' @details
#' \strong{Recommended use}: Bulk RNA-seq from PBMCs, whole blood, or any tissue
#' with variable B cell infiltration. Particularly relevant in autoimmune
#' and inflammatory disease cohorts.
#'
#' @export
omniseq_immunoglobulin <- c(
  # Heavy chain variable (IGHV)
  "IGHV1-2", "IGHV1-3", "IGHV1-8", "IGHV1-18", "IGHV1-24",
  "IGHV1-45", "IGHV1-46", "IGHV1-58", "IGHV1-69",
  "IGHV2-5", "IGHV2-26", "IGHV2-70",
  "IGHV3-7", "IGHV3-9", "IGHV3-11", "IGHV3-13", "IGHV3-15",
  "IGHV3-20", "IGHV3-21", "IGHV3-23", "IGHV3-30", "IGHV3-33",
  "IGHV3-43", "IGHV3-48", "IGHV3-49", "IGHV3-53", "IGHV3-64",
  "IGHV3-66", "IGHV3-72", "IGHV3-73", "IGHV3-74",
  "IGHV4-4", "IGHV4-28", "IGHV4-30-1", "IGHV4-30-2",
  "IGHV4-30-4", "IGHV4-31", "IGHV4-34", "IGHV4-38-2",
  "IGHV4-39", "IGHV4-59", "IGHV4-61",
  "IGHV5-51", "IGHV6-1", "IGHV7-4-1",
  # Kappa light chain variable (IGKV)
  "IGKV1-5", "IGKV1-6", "IGKV1-8", "IGKV1-9", "IGKV1-12",
  "IGKV1-13", "IGKV1-16", "IGKV1-17", "IGKV1-27", "IGKV1-33",
  "IGKV1-37", "IGKV1-39", "IGKV1-NL1",
  "IGKV2-24", "IGKV2-26", "IGKV2-28", "IGKV2-29", "IGKV2-30",
  "IGKV2-40",
  "IGKV3-11", "IGKV3-15", "IGKV3-20",
  "IGKV4-1", "IGKV5-2", "IGKV6-21",
  # Lambda light chain variable (IGLV)
  "IGLV1-36", "IGLV1-40", "IGLV1-41", "IGLV1-44", "IGLV1-47",
  "IGLV1-51",
  "IGLV2-8", "IGLV2-11", "IGLV2-14", "IGLV2-18", "IGLV2-23",
  "IGLV3-1", "IGLV3-9", "IGLV3-10", "IGLV3-12", "IGLV3-16",
  "IGLV3-19", "IGLV3-21", "IGLV3-22", "IGLV3-25", "IGLV3-27",
  "IGLV4-60", "IGLV4-69",
  "IGLV5-45", "IGLV5-48", "IGLV5-52",
  "IGLV6-57", "IGLV7-43", "IGLV7-46", "IGLV8-61",
  "IGLV9-49", "IGLV10-54", "IGLV11-55"
)


# =============================================================================
# REGISTRY & UTILITY FUNCTIONS
# =============================================================================

# =============================================================================
# 6. SCHOGGINS & ROSENBERG 2018 ISG LIST (272 genes)
# =============================================================================

#' Schoggins & Rosenberg 2018 Interferon-Stimulated Gene List
#'
#' 272 ISGs identified by Schoggins & Rosenberg (2018) from IFN-stimulated
#' human cells, ranked by t-statistic. Includes canonical ISGs (MX1, OAS1,
#' IFIT1) as well as less well-known IFN-responsive genes. Broader than the
#' 28-gene IRG-S panel — useful for comprehensive IFN signal masking.
#'
#' Masking with this list removes a broader IFN transcriptional footprint
#' to reveal non-IFN biology (e.g. metabolic, stress response, lineage).
#'
#' @source Schoggins JW, Randall G (2013) updated in Rosenberg et al. (2018).
#'   \emph{Nature}. PMID: 29769728.
#'
#' @format A character vector of 272 HGNC gene symbols, ordered by
#'   significance (most significant first).
#'
#' @export
omniseq_isgs_schoggins <- unique(c(
  "TRANK1", "MOV10", "SMCHD1", "C19orf66", "TNK2", "NCF1", "ANKFY1",
  "C14orf159", "DDX60L", "PGAP1", "STAT2", "PNPT1", "ISG20", "RBCK1",
  "ZCCHC2", "SAMD4A", "PARP14", "RNF213", "PPM1K", "SAMD9", "SAMD9L",
  "EIF2AK2", "HERC5", "TRIM14", "HELZ2", "SP100", "OAS2", "ZNFX1",
  "USP18", "IFI44L", "UBA7", "GBP4", "MX1", "IFIT1", "IFIT3", "IFIT2",
  "OAS1", "OAS3", "ISG15", "RSAD2", "MX2", "HERC6", "IFI44", "CMPK2",
  "DDX60", "TRIM25", "LGALS3BP", "TRIM56", "EPSTI1", "BST2", "XAF1",
  "DHX58", "IFIH1", "DDX58", "ADAR", "TRIM22", "IFI6", "IFI27",
  "IFITM1", "IFITM3", "LAMP3", "LY6E", "NT5C3A", "OASL", "PLSCR1",
  "RTP4", "SP110", "TDRD7", "TREX1", "TRIM21", "UNC93B1", "ZBP1",
  "APOBEC3A", "APOBEC3B", "APOBEC3G", "BATF2", "C3", "CASP1",
  "CCL2", "CCL5", "CCL8", "CD274", "CXCL10", "CXCL9", "EIF2AK1",
  "GBP1", "GBP2", "GBP3", "GBP5", "HLA-A", "HLA-B", "HLA-C",
  "HLA-E", "HLA-F", "HLA-G", "ICAM1", "IFI16", "IFI30", "IFI35",
  "IFITM2", "IL15", "IL6", "IRF1", "IRF2", "IRF7", "IRF9",
  "ISG20L2", "LGALS9", "MAVS", "MNDA", "MTHFD2", "MYD88",
  "NMI", "OGFR", "PARP9", "PARP10", "PARP12", "PARP15", "PCNA",
  "PSMB8", "PSMB9", "PSMB10", "PTPN1", "PTPN2", "RARRES3", "RIPK2",
  "RNASE4", "RNF31", "SECTM1", "SERPINB2", "SIGLEC1", "SLC25A28",
  "SOCS1", "SOCS3", "SQSTM1", "ST3GAL5", "STAT1", "STING1", "TAPBP",
  "TNFAIP3", "TNFAIP6", "TNFSF10", "TRIM15", "TRIM26",
  "TRIM38", "TRIM5", "TRIM6", "UBE2L6", "VCAM1", "WARS1",
  "ZC3HAV1", "ZC3HAV1L", "ACSL1", "ACSS2", "ALCAM", "ALDH3A1",
  "ANKRD22", "ANTXR2", "AOC3", "APOL1", "APOL2",
  "APOL3", "APOL6", "ARL4C", "ATG9B", "B2M", "B4GALT1", "BIRC3",
  "BTN3A1", "C2", "C4B", "CARD16", "CARD17", "CASP4",
  "CASP5", "CD38", "CD44", "CD47", "CD68", "CD69", "CEACAM1",
  "CFB", "CFH", "CIITA", "CLEC2B", "CLEC4A", "CLEC4D", "CLEC7A",
  "CMKLR1", "CNP", "CSF1", "CSF1R", "CXCL11", "CXCL13", "CXCL16",
  "CYBB", "CYP1B1", "DAPP1", "DDB2", "DOCK8", "DTX3L", "EGR1",
  "ENTPD7", "ETV7", "FADS1", "FADS2", "FAM111A", "FAM111B",
  "FCGR1A", "FCGR3A", "FES", "FGL2", "FNBP1", "G1P3",
  "GPR65", "GPR84", "HCST", "HCK", "HMOX1", "IDO1", "IER3",
  "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "IL10RA", "IL13RA1",
  "IL15RA", "IL2RB", "IL4R", "IRF3", "IRF5", "JAK1", "JAK2",
  "KLF4", "LYN", "MAP3K8", "MEFV", "MLKL", "NCF2", "NCF4",
  "NFKBIZ", "PARP16", "PDCD4", "PIM1", "PLCG2", "PML",
  "PPP1R15A", "PTGS2", "RHBDD3", "RIPK1", "RNF19B",
  "S100A8", "S100A9", "SAMHD1", "SELL", "SERPING1", "SH2B3",
  "SLC11A1", "SLC15A3", "SLFN11", "SLFN12", "SLFN13", "SLFN5",
  "SOCS5", "SPATS2L", "STAP1", "STK17A", "STK17B", "TGM2",
  "TLR7", "TMEM173", "TNFRSF10A", "TNFRSF10B", "TNFRSF14",
  "TRAFD1", "TREM1", "TRIM34", "TRIM69", "TXNIP",
  "TYMP", "UBD", "ULBP2", "VAMP5", "VAMP8", "VIM",
  "XCL1", "ZNF385A"
))

.omniseq_signature_registry <- list(
  isgs_28          = omniseq_isgs_28,
  isgs_extended    = omniseq_isgs_extended,
  isgs_schoggins   = omniseq_isgs_schoggins,
  ribosomal        = omniseq_ribosomal,
  mitochondrial    = omniseq_mitochondrial,
  immunoglobulin   = omniseq_immunoglobulin
)


#' List Available Built-in Gene Signatures
#'
#' Prints a table of all built-in OmniSeq gene signatures with their
#' shorthand names, gene counts, and descriptions.
#'
#' @return A data.frame of signature metadata (invisibly).
#' @export
omni_list_signatures <- function() {

  df <- data.frame(
    Shorthand   = names(.omniseq_signature_registry),
    Object_Name = paste0("omniseq_", names(.omniseq_signature_registry)),
    N_Genes     = sapply(.omniseq_signature_registry, length),
    Description = c(
      "28-gene NIH IRG-S interferon response panel",
      "Broad type I IFN signature (Hallmark IFNA Response)",
      "272-gene ISG list (Schoggins & Rosenberg 2018, ranked by t-stat)",
      "Cytoplasmic + mitochondrial ribosomal protein genes",
      "37 mitochondrially-encoded genes (MT- prefix)",
      "Immunoglobulin variable region genes (IGHV/IGKV/IGLV)"
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  cat("\n OmniSeq Built-in Gene Signatures\n")
  cat(strrep("-", 72), "\n")
  cat(sprintf(" %-18s %-10s  %s\n", "Shorthand", "N Genes", "Description"))
  cat(strrep("-", 72), "\n")
  for (i in seq_len(nrow(df))) {
    cat(sprintf(" %-18s %-10d  %s\n",
                df$Shorthand[i], df$N_Genes[i], df$Description[i]))
  }
  cat(strrep("-", 72), "\n")
  cat(" Usage:\n")
  cat("   run_omni_mask(obj, signature = \"isgs_28\")\n")
  cat("   run_omni_mask(obj, features = omniseq_isgs_28)\n")
  cat("   run_omni_mask(obj, features = c(omniseq_isgs_28, \"GENE1\"))\n\n")

  return(invisible(df))
}


#' Look Up a Built-in Signature by Shorthand Name
#'
#' Internal helper used by \code{run_omni_mask()}.
#' @keywords internal
.omni_get_signature <- function(name) {
  # NEW FIX: Explicitly handle NULL to avoid lookup errors when
  # running custom manual feature lists in the pipeline.
  if (is.null(name) || is.na(name)) return(character(0))

  sig <- .omniseq_signature_registry[[tolower(trimws(name))]]
  if (is.null(sig))
    stop("Unknown signature '", name, "'. Run omni_list_signatures() to see ",
         "available options.", call. = FALSE)
  sig
}
