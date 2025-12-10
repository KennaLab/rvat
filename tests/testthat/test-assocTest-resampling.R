data(GT)
data(GTsmall)

test_that("resampled assocTests snapshots are identical", {
  # snapshot resampled tests
  set.seed(10)
  warnings <- capture_warnings({
    test <- assocTest(
      GT,
      pheno = "pheno",
      covar = paste0("PC", 1:4),
      test = c(
        "skat",
        "skat_robust",
        "skat_burden",
        "skat_burden_robust",
        "skato_robust",
        "acatv"
      ),
      nResampling = 10,
      methodResampling = "permutation",
      verbose = FALSE
    )
  })
  metadata(test)$creationDate <- NA_character_
  metadata(test)$gdbPath <- NA_character_
  metadata(test)$rvatVersion <- NA_character_
  expect_snapshot_value(test, style = "serialize")
})

test_that("resamplingFile and resamplingMatrix are identical", {
  set.seed(10)
  resamplingfile_path <- withr::local_tempfile(fileext = ".gz")
  buildResamplingFile(
    nSamples = ncol(GTsmall),
    nResampling = 10,
    methodResampling = "permutation",
    output = resamplingfile_path
  )
  resamplingfile <- resamplingFile(resamplingfile_path)
  expect_true(stringr::str_detect(
    capture_output({
      show(resamplingfile)
    }),
    "resamplingFile object"
  ))

  test_resamplingfile <- assocTest(
    GTsmall,
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    test = "skat",
    resamplingFile = resamplingfile,
    verbose = FALSE,
    outputResampling = TRUE
  )
  ## expect identical results when using resamplingMatrix
  set.seed(10)
  resamplingmatrix <- buildResamplingFile(
    nSamples = ncol(GTsmall),
    nResampling = 10,
    methodResampling = "permutation"
  )
  test_resamplingmatrix <- assocTest(
    GTsmall,
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    test = "skat",
    resamplingMatrix = resamplingmatrix,
    verbose = FALSE,
    outputResampling = TRUE
  )
  metadata(test_resamplingfile)$creationDate <- NA_character_
  metadata(test_resamplingmatrix)$creationDate <- NA_character_
  expect_equal(test_resamplingfile, test_resamplingmatrix)
  }
)

test_that("outputting resamplings works", {
  
  # output resamplings to file
  outputresampling <- withr::local_tempfile()
  set.seed(10)
  test <- assocTest(
    GTsmall,
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    test = c("skat"),
    methodResampling = "permutation",
    nResampling = 10,
    outputResampling = outputresampling,
    verbose = FALSE
  )
  outputresampling <- readr::read_tsv(
    outputresampling,
    show_col_types = FALSE,
    col_types = list(MAFweight = "character")
  )

  set.seed(10)
  test <- assocTest(
    GTsmall,
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    test = c("skat"),
    methodResampling = "permutation",
    nResampling = 10,
    outputResampling = TRUE,
    verbose = FALSE
  )
  expect_equal(as.data.frame(outputresampling), as.data.frame(test))
})
