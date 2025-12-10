data(GTsmall, envir = environment())
sv_results <- assocTest(
  GTsmall,
  covar = paste0("PC", 1:4),
  test = "scoreSPA",
  singlevar = TRUE,
  pheno = "pheno",
  verbose = FALSE
)

test_that("rvatResult summary methods work", {
  data(rvbresults, envir = environment())

  rvb_summary <- capture_output(summary(rvbresults))
  expect_true(stringr::str_detect(rvb_summary, "n units = "))

  sv_summary <- capture_output(summary(sv_results))
  expect_true(stringr::str_detect(sv_summary, "n vars = "))
})

test_that("empty result constructors work", {
  singlevar_empty <- singlevarResult()
  expect_equal(nrow(singlevar_empty), 0)
  expect_s4_class(singlevar_empty, "singlevarResult")

  rvbresults_empty <- rvbResult()
  expect_equal(nrow(rvbresults_empty), 0)
  expect_s4_class(rvbresults_empty, "rvbResult")

  gsaresults_empty <- gsaResult()
  expect_equal(nrow(gsaresults_empty), 0)
  expect_s4_class(gsaresults_empty, "gsaResult")
})

test_that("rvatResult merging works", {
  data(rvbresults, envir = environment())

  # merge results
  rvbresults_merge <- merge(
    rvbresults[, colnames(rvbresults)[1:24]],
    dplyr::distinct(as.data.frame(rvbresults)[, colnames(rvbresults)[c(
      1,
      25:ncol(rvbresults)
    )]]),
    by = "unit"
  )
  expect_equal(rvbresults, rvbresults_merge)

  # merge results (DataFrame)
  rvbresults_merge_DFrame <- merge(
    rvbresults[, colnames(rvbresults)[1:24]],
    DataFrame(dplyr::distinct(as.data.frame(rvbresults)[, colnames(
      rvbresults
    )[c(
      1,
      25:ncol(rvbresults)
    )]])),
    by = "unit"
  )
  expect_equal(rvbresults_merge, rvbresults_merge_DFrame)
})


test_that("misc rvatResult methods work", {
  data(rvbresults, envir = environment())
  id_col <- getIdCol(rvbresults)
  expect_equal(id_col, "unit")
})


test_that("rvatResult input validation works", {

  data(rvbresults, envir = environment())
  # expect error when object is of invalid type
  expect_error(
    {
      rvbResult(5)
    },
    "should be either a"
  )

  # checkClassrvatResult
  resultfile <- withr::local_tempfile()
  writeResult(rvbresults, resultfile)

  # expect error when type is invalid
  expect_error(
    {
      checkClassrvatResult(resultfile, type = "invalid")
    },
    regexp = "Type should be one of the following"
  )

  # expect error when object is of invalid type
  expect_error(
    {
      checkClassrvatResult(TRUE, type = "rvbResult")
    },
    regexp = "Input must be an rvatResult object or a single file path."
  )

  # validity
  rvbresults_invalid <- rvbresults
  rvbresults_invalid$P <- as.character(rvbresults_invalid$P)
  expect_error(
    {
      validObject(rvbresults_invalid)
    },
    regexp = "The following columns should be of type: P: numeric"
  )
  rvbresults_invalid$P <- NULL
  expect_error(
    {
      validObject(rvbresults_invalid)
    },
    regexp = "The following columns are missing: P"
  )
})
