data(GT)
data(GTsmall)
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))

# run_tests_GT etc -> defined in helper-functions.R
test_that("assocTest-GT snapshots are identical", {
  params <- generate_param()
  warnings <- capture_warnings(
    run_test <- run_tests_GT(GT, params, var_sv = rownames(GT)[1:25])
  )
  expect_true(all(
    c(
      "There are p-values that are exactly 0!",
      "no non-missing arguments to min; returning Inf"
    ) %in%
      warnings
  ))
  expect_snapshot_value(run_test, style = "serialize")
})

test_that("assocTest-gdb works", {
  gdb <- create_example_gdb()
  moderate_varsetfile <- withr::local_tempfile()
  merged_varsetfile <- withr::local_tempfile()

  moderate_varsets <- getVarSet(
    varsetfile,
    varSetName = "Moderate",
    unit = c("CYP19A1", "FUS", "OPTN")
  )
  merged_varsets <- getVarSet(
    varsetfile,
    varSetName = c("Moderate", "High", "CADD"),
    unit = c("CYP19A1", "FUS", "OPTN")
  )
  write(moderate_varsets, moderate_varsetfile)
  write(merged_varsets, merged_varsetfile)

  ## looping through either a varSetList or a varSetFile should yield identical result
  varset <- varSetFile(moderate_varsetfile)
  test1 <- assocTest(
    gdb,
    cohort = "pheno",
    varSet = varset,
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    test = "scoreSPA",
    verbose = FALSE
  )
  metadata(test1)$creationDate <- NA_character_

  test2 <- assocTest(
    gdb,
    cohort = "pheno",
    varSet = getVarSet(varset, unit = listUnits(varset)),
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    test = "scoreSPA",
    verbose = FALSE
  )
  metadata(test2)$creationDate <- NA_character_
  expect_identical(test1, test2)

  ## looping through varSet should give same results as separate varSets (to ensure re-using pheno cohort etc. works as it should)
  varset <- varSetFile(merged_varsetfile)
  warnings <- capture_warnings({
    test1 <- assocTest(
      gdb,
      cohort = "pheno",
      varSet = varset,
      pheno = "pheno",
      geneticModel = c("allelic", "dominant", "recessive"),
      covar = paste0("PC", 1:4),
      maxMAF = 0.0001,
      test = "scoreSPA",
      verbose = FALSE
    )
  })
  expect_true(all(
    warnings %in%
      c(
        "NAs introduced by coercion",
        "There are p-values that are exactly 0!",
        "no non-missing arguments to min; returning Inf",
        "Less than two samples have a non-zero burden score, skipping tests."
      ) |
      grepl("weights are missing", warnings)
  ))

  test2 <- list()
  for (i in 1:length(varset)) {
    warnings <- capture_warnings({
      test <- assocTest(
        gdb,
        cohort = "pheno",
        varSet = getVarSet(
          varset,
          unit = listUnits(varset)[i],
          varSetName = listVarSets(varset)[i]
        ),
        pheno = "pheno",
        geneticModel = c("allelic", "dominant", "recessive"),
        covar = paste0("PC", 1:4),
        maxMAF = 0.0001,
        test = "scoreSPA",
        verbose = FALSE
      )
    })
    expect_true(all(
      warnings %in%
        c(
          "NAs introduced by coercion",
          "There are p-values that are exactly 0!",
          "no non-missing arguments to min; returning Inf",
          "Less than two samples have a non-zero burden score, skipping tests."
        ) |
        grepl("weights are missing", warnings)
    ))
    test2[[i]] <- test
  }
  test2 <- do.call(rbind, test2)
  test1$ID <- paste(test1$varSetName, test1$unit, test1$geneticModel)
  test2$ID <- paste(test2$varSetName, test2$unit, test2$geneticModel)
  test2 <- test2[match(test1$ID, test2$ID), ]
  metadata(test1)$creationDate <- NA_character_
  metadata(test2)$creationDate <- NA_character_
  expect_identical(test1, test2)

  # also compare with running directly on genoMatrix
  test3 <- list()
  j <- 1
  for (i in 1:length(varset)) {
    varSet <- getVarSet(
      varset,
      unit = listUnits(varset)[i],
      varSetName = listVarSets(varset)[i]
    )
    gt <- suppressWarnings(getGT(
      gdb,
      varSet = varSet,
      cohort = "pheno",
      verbose = FALSE
    ))

    for (geneticmodel in c("allelic", "dominant", "recessive")) {
      test <- suppressWarnings(assocTest(
        gt,
        pheno = "pheno",
        geneticModel = geneticmodel,
        covar = paste0("PC", 1:4),
        maxMAF = 0.0001,
        test = "scoreSPA",
        verbose = FALSE
      ))
      test3[[j]] <- test
      j <- j + 1
    }
  }
  test3 <- do.call(rbind, test3)
  test3$ID <- paste(test3$varSetName, test3$unit, test3$geneticModel)
  test3 <- test3[match(test2$ID, test3$ID), ]
  metadata(test3)$creationDate <- NA_character_
  expect_identical(test2, test3)

  # test writing to output
  varset <- varSetFile(moderate_varsetfile)
  output_assoctest <- withr::local_tempfile()
  test_output <- assocTest(
    gdb,
    cohort = "pheno",
    varSet = varset,
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    test = "scoreSPA",
    verbose = FALSE,
    output = output_assoctest
  )
  test_from_file <- rvbResult(output_assoctest)
  expect_equal(as.data.frame(test_output), as.data.frame(test_from_file))
  
  ## also check for singlevar
  varset <- varSetFile(moderate_varsetfile)
  output_assoctest <- withr::local_tempfile()
  test_output <- assocTest(
    gdb,
    cohort = "pheno",
    varSet = varset,
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    test = "scoreSPA",
    singlevar = TRUE,
    verbose = FALSE,
    output = output_assoctest
  )
  test_from_file <- singlevarResult(output_assoctest)
  expect_equal(as.data.frame(test_output), as.data.frame(test_from_file))
})

test_that("assocTest-gdb and assocTest-GT result in identical output", {
  params <- generate_param()
  warnings <- capture_warnings(
    run_test_gdb <- run_tests_gdb(
      gdb(rvat_example("rvatData.gdb")),
      params,
      VAR_id = rownames(GT),
      var_sv = rownames(GT)[1:25]
    )
  )
  expect_true(all(
    c(
      "There are p-values that are exactly 0!",
      "no non-missing arguments to min; returning Inf"
    ) %in%
      warnings
  ))
  warnings <- capture_warnings(
    run_test_GT <- run_tests_GT(GT, params, var_sv = rownames(GT)[1:25])
  )
  expect_true(all(
    c(
      "There are p-values that are exactly 0!",
      "no non-missing arguments to min; returning Inf"
    ) %in%
      warnings
  ))
  run_test_gdb$rvb$binary$varSetName = "unnamed"
  run_test_gdb$rvb$binary$unit = "unnamed"
  run_test_gdb$rvb$continuous$varSetName = "unnamed"
  run_test_gdb$rvb$continuous$unit = "unnamed"
  run_test_gdb$sv$binary$varSetName = "unnamed"
  run_test_gdb$sv$continuous$varSetName = "unnamed"
  run_test_gdb$sv$binary$caseMAF = NA_real_
  run_test_gdb$sv$continuous$caseMAF = NA_real_
  run_test_gdb$sv$binary$ctrlMAF = NA_real_
  run_test_gdb$sv$continuous$ctrlMAF = NA_real_
  run_test_GT$sv$binary$caseMAF = NA_real_
  run_test_GT$sv$continuous$caseMAF = NA_real_
  run_test_GT$sv$binary$ctrlMAF = NA_real_
  run_test_GT$sv$continuous$ctrlMAF = NA_real_
  rownames(run_test_GT$sv$continuous) <- NULL
  rownames(run_test_GT$sv$binary) <- NULL
  rownames(run_test_gdb$sv$continuous) <- NULL
  rownames(run_test_gdb$sv$binary) <- NULL
  expect_equal(run_test_gdb, run_test_GT, tolerance = 1.490116e-08)
})
