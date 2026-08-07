test_that("new_omniseq + validate_omniseq build a well-formed object", {
  skip_if_not_installed("SummarizedExperiment")

  counts <- matrix(1:12, nrow = 3, ncol = 4,
                    dimnames = list(c("TP53", "GAPDH", "ACTB"),
                                    c("S1", "S2", "S3", "S4")))
  meta <- data.frame(Genotype = c("WT", "WT", "KO", "KO"),
                      row.names = c("S1", "S2", "S3", "S4"))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), colData = meta
  )

  obj <- new_omniseq(data = se, data_type = "gene", condition = "Genotype")

  expect_s3_class(obj, "OmniSeq")
  expect_silent(validate_omniseq(obj))
  expect_length(omni_contrasts(obj), 0)
  expect_null(omni_model(obj))
  expect_null(omni_vst(obj))
  expect_null(omni_results(obj))
  expect_equal(dim(omni_counts(obj)), c(3, 4))
  expect_equal(omni_metadata(obj)$Genotype, meta$Genotype)
})

test_that("validate_omniseq rejects non-OmniSeq and malformed objects", {
  expect_error(validate_omniseq(list()), "not of class 'OmniSeq'")

  broken <- structure(list(data = "not a SummarizedExperiment"), class = "OmniSeq")
  expect_error(validate_omniseq(broken), "Missing slots")
})
