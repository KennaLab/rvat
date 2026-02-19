data(GTsmall)

test_that("assocTest input validation work", {
  GT_test <- GTsmall

  # expect a warning when a covariate with zero covariance is included
  GT_test$covar_novar <- 1
  expect_warning(
    {
      test <- assocTest(
        GT_test,
        test = "glm",
        covar = c(paste0("PC", 1:4), "covar_novar"),
        pheno = "pheno",
        verbose = FALSE
      )
    },
    regexp = "The following covariate\\(s\\) have zero covariance: covar_novar"
  )
  ## also: the covariate should not be present in the output
  expect_identical(as.character(test$covar), "PC1,PC2,PC3,PC4")

  # expect error when covariate is not a a vector/list
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "glm",
        covar = TRUE,
        pheno = "pheno",
        verbose = FALSE
      )
    },
    regexp = "`covar` must be either a character vector or a list"
  )

  # expect error when covariates are not present in GT
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "glm",
        covar = c("hello", "world"),
        pheno = "pheno",
        verbose = FALSE
      )
    },
    regexp = "the following `covar` fields were not found in colData"
  )

  # expect error when binary phenotype is specified for non-binary phenotype
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "glm",
        covar = paste0("PC", 1:4),
        pheno = "PC1",
        verbose = FALSE
      )
    },
    regexp = "should be coded"
  )

  # expect error when phenotype is non-numeric
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "glm",
        covar = paste0("PC", 1:4),
        pheno = "superPop",
        verbose = FALSE
      )
    },
    regexp = "Phenotype 'superPop' should be numeric"
  )

  # keep list
  keep <- colnames(GT_test)[1:100]
  expect_true(any(stringr::str_detect(
    capture_messages(
      {
        assocTest(
          GT_test,
          keep = keep,
          test = "glm",
          covar = paste0("PC", 1:4),
          pheno = "pheno",
          verbose = TRUE
        )
      }
    ),
    sprintf("100/%s", ncol(GT_test))
  )))

  # expect warning when missing weights are included
  tmp <- GT_test
  rowData(tmp)$w[1:5] <- NA_real_
  expect_warning(
    {
      test <- assocTest(
        tmp,
        test = c("glm"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        verbose = FALSE
      )
    },
    regexp = "are missing"
  )

  # expect warning when negative weights are included
  tmp <- GT_test
  rowData(tmp)$w[1:5] <- -1
  expect_warning(
    {
      test <- assocTest(
        tmp,
        test = c("glm"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        verbose = FALSE
      )
    },
    regexp = "are < 0"
  )
  rm(tmp)

  # expect empty results when no variants pass thresholds
  test <- assocTest(
    GT_test,
    test = c("glm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    minMAF = 0.011,
    maxMAF = 0.011,
    verbose = FALSE
  )
  expect_true(nrow(test) == 0 && is(test, "rvatResult"))

  # expect empty results when no samples pass thresholds
  test <- assocTest(
    GT_test,
    test = c("glm"),
    covar = c(paste0("PC", 1:4)),
    minCallrateSM = 0.9,
    maxCallrateSM = 0.9,
    pheno = "pheno",
    verbose = FALSE
  )
  expect_true(nrow(test) == 0 && is(test, "rvatResult"))

  # expect an error when a non-implemented test is specified
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = c("glm", "test"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        verbose = FALSE
      )
    },
    regexp = "The following tests are not valid: test"
  )

  # expect a warning when tests are specified that are only valid for binary tests or vice versa
  expect_warning(
    {
      test <- assocTest(
        GT_test,
        test = c("glm", "lm"),
        covar = c(paste0("PC", 2:4)),
        pheno = "PC1",
        verbose = FALSE,
        continuous = TRUE
      )
    },
    regexp = "The following tests were excluded since they are not available for continuous rvb tests: glm"
  )

  # expect an error when a non-existing phenotype is specified
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = c("glm"),
        covar = c(paste0("PC", 1:4)),
        pheno = "test",
        verbose = FALSE
      )
    },
    regexp = "'test' is not present in `colData\\(GT\\)`"
  )

  # expect error the genoMatrix is set to dominant/recessive
  # and another geneticModel is specified
  expect_error(
    {
      test <- assocTest(
        recode(GT_test, geneticModel = "recessive"),
        test = c("glm"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        geneticModel = "dominant",
        verbose = FALSE
      )
    },
    regexp = "in order to apply"
  )

  # expect warning the genoMatrix is set to dominant/recessive
  # and 'allelic' geneticModel is specified
  expect_error(
    {
      test <- assocTest(
        recode(GT_test, geneticModel = "recessive"),
        test = c("glm"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        geneticModel = "allelic",
        verbose = FALSE
      )
    },
    regexp = "Current geneticModel should be"
  )

  # expect error when MAC/MAF filters are set while geneticModel != allelic
  expect_error(
    {
      test <- assocTest(
        recode(GT_test, geneticModel = "recessive"),
        test = c("glm"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        geneticModel = "dominant",
        maxMAC = 10,
        verbose = FALSE
      )
    },
    regexp = "do not apply"
  )

  # expect error when geneticModel is not implemented
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = c("glm"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        geneticModel = "non-existent",
        verbose = FALSE
      )
    },
    regexp = "`geneticModel` should be"
  )

  # export error when non-implemented imputeMethod is specified
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = c("glm"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        imputeMethod = "hello!",
        verbose = FALSE
      )
    },
    regexp = "should be either"
  )

  # warning on non-implemented tests
  expect_warning(
    {
      test <- assocTest(
        GT_test,
        test = c("glm", "lm"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        verbose = FALSE
      )
    },
    regexp = "The following tests were excluded"
  )

  expect_warning(
    {
      test <- assocTest(
        GT_test,
        test = c("glm", "lm"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        singlevar = TRUE,
        verbose = FALSE
      )
    },
    regexp = "The following tests were excluded"
  )

  expect_warning(
    {
      test <- assocTest(
        GT_test,
        test = c("glm", "lm"),
        covar = c(paste0("PC", 2:4)),
        pheno = "PC1",
        continuous = TRUE,
        verbose = FALSE
      )
    },
    regexp = "The following tests were excluded"
  )

  expect_warning(
    {
      test <- assocTest(
        GT_test,
        test = c("glm", "lm"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        singlevar = TRUE,
        continuous = TRUE,
        verbose = FALSE
      )
    },
    regexp = "The following tests were excluded"
  )

  # errror if no applicable tests are left
  expect_error(
    {
      suppressWarnings(
        test <- assocTest(
          GT_test,
          test = c("lm"),
          covar = c(paste0("PC", 1:4)),
          pheno = "pheno",
          verbose = FALSE
        )
      )
    },
    regexp = "No applicable tests"
  )

  # only one offset can be specified
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "firth",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        offset = c("PC", 5:6),
        verbose = FALSE
      )
    },
    regexp = "Currently at most"
  )

  # offset should be available
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "firth",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        offset = "hello!",
        verbose = FALSE
      )
    },
    regexp = "not available"
  )

  # warning
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "firth",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        offset = "hello!",
        verbose = FALSE
      )
    },
    regexp = "not available"
  )

  # resampling

  ## not implemented for singlevar tetss
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "firth",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        methodResampling = "permutation",
        verbose = FALSE,
        singlevar = TRUE
      )
    },
    regexp = "not implemented"
  )

  ## throw warning if some tests are not valid
  expect_warning(
    {
      test <- assocTest(
        GT_test,
        test = c("skat", "firth"),
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        methodResampling = "permutation",
        nResampling = 2,
        verbose = FALSE
      )
    },
    regexp = "not implemented"
  )

  ## throw error if no valid tests are implemented
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "firth",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        methodResampling = "permutation",
        verbose = FALSE
      )
    },
    regexp = "not implemented"
  )

  ## throw error if no valid permutation methods are specified
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "skat",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        methodResampling = "hello!",
        verbose = FALSE
      )
    },
    regexp = "accepted option"
  )

  ## expect errror when resamplingMatrix is of wrong type
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "skat",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        resamplingMatrix = TRUE,
        verbose = FALSE
      )
    },
    regexp = "`resamplingMatrix` should either be `NULL`"
  )

  ## expect errror when resamplingFile is of wrong type
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "skat",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        resamplingFile = TRUE,
        verbose = FALSE
      )
    },
    regexp = "should be an object of class"
  )

  ## expect error when outputResampling is of wrong type
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "skat",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        methodResampling = "permutation",
        outputResampling = data.frame(),
        verbose = FALSE
      )
    },
    regexp = "`outputResampling` should be a filepath or a boolean"
  )

  ## expect error when number of rows in resamplingMatrix doesn't match
  resamplingmatrix <- buildResamplingFile(
    nSamples = 1000,
    nResampling = 1000,
    methodResampling = "permutation"
  )
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "skat",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        resamplingMatrix = resamplingmatrix,
        verbose = FALSE
      )
    },
    regexp = "Number of rows in resamplingMatrix should match"
  )

  ## expect error when number of rows in resamplingFile doesn't match
  resamplingfile <- withr::local_tempfile()
  buildResamplingFile(
    nSamples = 1000,
    nResampling = 1000,
    methodResampling = "permutation",
    output = resamplingfile
  )
  expect_error(
    {
      test <- assocTest(
        GTsmall,
        test = "skat",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        resamplingFile = resamplingFile(resamplingfile),
        verbose = FALSE
      )
    },
    regexp = "Number of samples in resamplingFile does not match"
  )

  ## expect error when both a resamplingMatrix and methodResampling are specified
  resamplingmatrix <- buildResamplingFile(
    nSamples = ncol(GT_test),
    nResampling = 100,
    methodResampling = "permutation"
  )
  expect_error(
    {
      test <- assocTest(
        GT_test,
        test = "skat",
        covar = c(paste0("PC", 1:4)),
        pheno = "pheno",
        resamplingMatrix = resamplingmatrix,
        methodResampling = "permutation",
        verbose = FALSE
      )
    },
    regexp = "Only one of"
  )

  # expect error when supplying varSetFile from different gdb
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
  test_gdb <- withr::local_tempfile(fileext = ".gdb")
  buildGdb(
    test_path("data/vcf_multiallelic.vcf"),
    output = test_gdb,
    genomeBuild = "GRCh38",
    verbose = FALSE
  )
  test_gdb <- gdb(test_gdb)
  expect_error(
    {
      assocTest(
        test_gdb,
        varSet = varsetfile,
        cohort = "SM",
        pheno = "pheno",
        test = "firth",
        verbose = FALSE
      )
    },
    regexp = "generated from a different gdb"
  )
})

