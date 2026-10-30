test_that("Shiny app sources and its demo data builds a valid object", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("SummarizedExperiment")

  app_dir <- system.file("shiny-app", package = "OmniSeq")
  skip_if(app_dir == "", "App directory not installed")

  env <- new.env(parent = asNamespace("OmniSeq"))
  for (f in list.files(file.path(app_dir, "R"), full.names = TRUE))
    sys.source(f, envir = env)

  d <- env$demo_data(n_genes = 500)
  expect_equal(dim(d$counts), c(500L, 8L))
  expect_true(all(omniseq_isgs_28 %in% rownames(d$counts)))
  expect_length(env$check_inputs(as.data.frame(d$counts), d$meta)$problems, 0)

  bad <- as.data.frame(d$counts)
  bad[1, 1] <- 1.5
  expect_match(env$check_inputs(bad, d$meta)$problems, "non-integer")

  obj <- suppressMessages(init_omniseq(d$counts, d$meta,
                                       condition_col = "Condition"))
  expect_s3_class(obj, "OmniSeq")
})
