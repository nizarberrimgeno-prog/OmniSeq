test_that(".detect_data_type identifies TE locus IDs", {
  te_ids <- c(
    "1000:SINE:MIRb", "1000014:ERV:MSTB", "2000:LINE:L1",
    "3000:DNA:hAT", "4000:Satellite:ALR", "5000:RC:Helitron",
    "6000:rRNA:5S", "7000:tRNA:Ala", "8000:scRNA:7SL", "9000:Unknown:X"
  )
  mat <- matrix(1, nrow = length(te_ids), ncol = 1,
                dimnames = list(te_ids, "S1"))

  expect_equal(.detect_data_type(mat), "TE")
})

test_that(".detect_data_type identifies gene symbols", {
  genes <- c("TP53", "GAPDH", "ACTB", "IFI27", "MX1", "STAT1")
  mat <- matrix(1, nrow = length(genes), ncol = 1,
                dimnames = list(genes, "S1"))

  expect_equal(.detect_data_type(mat), "gene")
})

test_that(".detect_data_type falls back to gene below the 80% TE threshold", {
  mixed <- c("1000:SINE:MIRb", "TP53", "GAPDH", "ACTB", "MX1")
  mat <- matrix(1, nrow = length(mixed), ncol = 1,
                dimnames = list(mixed, "S1"))

  expect_equal(.detect_data_type(mat), "gene")
})
