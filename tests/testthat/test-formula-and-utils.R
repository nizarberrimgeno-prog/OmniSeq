test_that(".build_formula includes covariates and condition", {
  f <- .build_formula(covariates = "Batch", condition = "Genotype")

  expect_s3_class(f, "formula")
  expect_equal(all.vars(f), c("Batch", "Genotype"))
})

test_that(".build_formula works with no covariates", {
  f <- .build_formula(covariates = character(0), condition = "Genotype")

  expect_equal(all.vars(f), "Genotype")
})

test_that("%||% falls back only when the left side is NULL", {
  expect_equal(5 %||% 10, 5)
  expect_equal(NULL %||% 10, 10)
  expect_null(NULL %||% NULL)
})

test_that(".check_class accepts OmniSeq objects and rejects everything else", {
  fake <- structure(list(), class = "OmniSeq")

  expect_silent(.check_class(fake))
  expect_error(.check_class(list()), "Expected an OmniSeq object")
})
