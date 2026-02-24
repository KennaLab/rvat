# generate sv results for testing
data(GTsmall, envir = environment())
sv_results <- assocTest(
  GTsmall,
  covar = paste0("PC", 1:4),
  test = "scoreSPA",
  singlevar = TRUE,
  pheno = "pheno",
  verbose = FALSE
)

test_that("rvbResult read/write roundtrip preserves data", {
  data(rvbresults, envir = environment())

  resultfile <- withr::local_tempfile()
  writeResult(rvbresults, file = resultfile)
  rvbresults_read <- rvbResult(resultfile)

  expect_equal(rvbresults, rvbresults_read)
})

test_that("rvbResult appending works correctly", {
  data(rvbresults, envir = environment())
  resultfile <- withr::local_tempfile()

  # write results in chunks using append
  writeResult(rvbresults[1:100, ], file = resultfile)
  writeResult(
    rvbresults[101:nrow(rvbresults), ],
    file = resultfile,
    append = TRUE
  )

  # should match writing in one go
  resultfile_onechunk <- withr::local_tempfile()
  writeResult(rvbresults, file = resultfile_onechunk)

  expect_equal(
    rvbResult(resultfile),
    rvbResult(resultfile_onechunk)
  )
})

test_that("readResults infers rvbResult type correctly", {
  data(rvbresults, envir = environment())

  resultfile <- withr::local_tempfile()
  writeResult(rvbresults, file = resultfile)

  # read results without specifying type, should be inferred
  rvbresults_read <- readResults(resultfile)
  expect_s4_class(rvbresults_read, "rvbResult")
  expect_equal(rvbresults, rvbresults_read)

  # write result w/o header info
  resultfile_no_header <- withr::local_tempfile()
  readr::write_tsv(
    as.data.frame(rvbresults),
    file = resultfile_no_header,
    progress = FALSE
  )

  ## expect warning when reading
  expect_warning(
    rvbresults_read <- rvbResult(resultfile_no_header),
    regexp = "The following metadata fields are missing"
  )

  ## results should be identical
  expect_equal(
    as.data.frame(rvbresults),
    as.data.frame(rvbresults_read)
  )
})

test_that("singlevarResult read/write roundtrip preserves data", {
  resultfile <- withr::local_tempfile()
  writeResult(sv_results, file = resultfile)
  sv_results_read <- singlevarResult(resultfile)

  expect_equal(sv_results, sv_results_read)
})

test_that("readResults infers singlevarResult type correctly", {
  resultfile <- withr::local_tempfile()
  writeResult(sv_results, file = resultfile)

  sv_results_read <- readResults(resultfile)
  expect_s4_class(sv_results_read, "singlevarResult")
  expect_equal(sv_results, sv_results_read)

  # read result w/o header info
  resultfile_no_header <- withr::local_tempfile()
  readr::write_tsv(
    as.data.frame(sv_results),
    file = resultfile_no_header,
    progress = FALSE
  )

  expect_warning(
    sv_read <- singlevarResult(resultfile_no_header),
    regexp = "The following metadata"
  )

  expect_equal(
    as.data.frame(sv_results),
    as.data.frame(sv_read)
  )
})


test_that("reading and writing rvatResults input validation works", {
  # reading rvbResult as singlevarResult fails
  data(rvbresults, envir = environment())

  resultfile <- withr::local_tempfile()
  writeResult(rvbresults, file = resultfile)

  expect_error(singlevarResult(resultfile), regexp = "Unexpected filetype")

  # reading singlevarResult as rvbResult fails
  resultfile <- withr::local_tempfile()
  writeResult(sv_results, file = resultfile)
  expect_error(rvbResult(resultfile), regexp = "Unexpected filetype")

  # expect error when required columns are missing
  resultfile <- withr::local_tempfile()
  writeResult(
    rvbresults[, colnames(rvbresults) != "P"],
    file = resultfile
  )
  expect_error(
    {
      rvbresult_read <- rvbResult(resultfile)
    },
    regexp = "The results file is missing the following required"
  )

  ##  expect error when indistinguishable
  resultfile_mixed <- withr::local_tempfile()
  rvbresults_mixed <- rvbresults
  rvbresults_mixed$VAR_id <- "1"
  rvbresults_mixed$caseCarriers <- NULL
  writeResult(rvbresults_mixed, resultfile_mixed)
  expect_error(
    {
      checkClassrvatResult(resultfile_mixed)
    },
    regexp = "Could not infer result type"
  )

  expect_error(
    {
      checkClassrvatResult(resultfile_mixed, type = "rvbResult")
    },
    regexp = "The file is missing required"
  )
})
