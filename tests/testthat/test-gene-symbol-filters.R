test_that(".clean_gene_symbols keeps valid HGNC symbols", {
  genes <- c(
    "TP53", "ENSG00000141510", "IFI27",
    "gene [Source:HGNC]", "novel_transcript_1"
  )

  result <- .clean_gene_symbols(genes)

  expect_equal(result, c("TP53", NA, "IFI27", NA, NA))
})

test_that(".filter_coding_genes keeps gene symbols and TE locus IDs, drops the rest", {
  genes <- c(
    "TP53",               # keep: clean HGNC symbol
    "1000014:ERV:MSTB",   # keep: TE locus ID
    "ENSG00000141510",    # drop: Ensembl ID
    "LINC00152",          # drop: lncRNA accession
    "RNU6-1",             # drop: small ncRNA
    "AL123456.1"          # drop: lncRNA accession
  )

  mask <- .filter_coding_genes(genes)

  expect_equal(mask, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
})
