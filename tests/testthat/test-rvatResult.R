data(rvbresults)
data(GTsmall)
rvbresults_test <- rvbresults

test_that("reading and writing rvbResults works" ,{
  resultfile <- withr::local_tempfile()
  # write results
  writeResult(rvbresults_test, file = resultfile)
  rvbresults_read <- rvbResult(resultfile)
  expect_identical(rvbresults_test, rvbresults_read, tolerance = 1e-10)
  
  # readresults without specifying type
  rvbresults_read <- readResults(resultfile)
  expect_identical(rvbresults_test, rvbresults_read, tolerance = 1e-10)
  
  # load results without a header
  rvbresults_test$CHROM <- rvbresults_test$POS <- rvbresults_test$start <- rvbresults_test$end <- NULL
  resultfile_noheader <- withr::local_tempfile()
  writeResult(rvbresults_test, file = resultfile_noheader, col.names = FALSE)
  expect_error(rvbResult(resultfile_noheader))
  resultfile_noheader_check <- readResults(resultfile_noheader, header = FALSE, type = "rvbResult")
  expect_equal(resultfile_noheader_check, rvbresults_test)
  
  # add some fields
  rvbresults_test$gene <- sample(LETTERS, nrow(rvbresults_test), replace = TRUE)
  rvbresults_test$transcript <- sample(LETTERS, nrow(rvbresults_test), replace = TRUE)
  rvbresults_test$score <- rnorm(n = nrow(rvbresults_test))
  writeResult(rvbresults_test, file = resultfile)
  rvbresults_read <- rvbResult(resultfile)
  expect_identical(rvbresults_test, rvbresults_read, tolerance = 1e-10)
  
  resultfile_noheader <- withr::local_tempfile()
  writeResult(rvbresults_test, file = resultfile_noheader, col.names = FALSE)
  expect_error(rvbResult(resultfile_noheader))
  expect_message({resultfile_noheader_check <- readResults(resultfile_noheader, header = FALSE, type = "rvbResult")},
                 regexp = "Number of columns is larger than the default")
  colnames(resultfile_noheader_check)[25:27] <- c("gene", "transcript", "score")
  expect_equal(resultfile_noheader_check, rvbresults_test)

  # expect error when reading as singlevarResult
  expect_error(singlevarResult(resultfile))
  
  # svresults
  sv_results <- assocTest(
    GTsmall,
    covar = paste0("PC", 1:4),
    test = "scoreSPA",
    singlevar = TRUE,
    pheno = "pheno",
    verbose = FALSE
  )
  writeResult(sv_results, file = resultfile)
  sv_results_read <- singlevarResult(resultfile)
  expect_identical(sv_results, sv_results_read, tolerance = 1e-10)
  
  # readresults without specifying type
  sv_results_read <- readResults(resultfile)
  expect_identical(sv_results, sv_results_read, tolerance = 1e-10)
  
  # results w/o header
  writeResult(sv_results, file = resultfile_noheader, col.names = FALSE)
  expect_error(singlevarResult(resultfile_noheader))
  resultfile_noheader_check <- readResults(resultfile_noheader, header = FALSE, type = "singlevarResult")
  expect_equal(resultfile_noheader_check, sv_results)
  
  
  }
)

test_that("misc rvatResult methods work" ,{
  # check summary
  expect_true(stringr::str_detect(capture_output({summary(rvbresults)}), "n units = "))
  
  # check summary
  sv <- assocTest(
    GTsmall,
    test = "scoreSPA",
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    singlevar = TRUE,
    verbose = FALSE
  )
  expect_true(stringr::str_detect(capture_output({summary(sv)}), "n vars = "))
  
  # empty results
  singlevar_empty <- singlevarResult()
  expect_equal(nrow(singlevar_empty), 0)
  rvbresults_empty <- rvbResult()
  expect_equal(nrow(rvbresults_empty), 0)
  gsaresults_empty <- gsaResult()
  expect_equal(nrow(gsaresults_empty), 0)
  
  # merge results
  rvbresults_check <- merge(rvbresults[,colnames(rvbresults)[1:24]], 
                            distinct(as.data.frame(rvbresults[,colnames(rvbresults)[c(1,25:ncol(rvbresults))]])),
                            by = "unit"
                            )
  expect_equal(rvbresults, rvbresults_check)
  
}
)

test_that("plots work", {
  expect_no_error(qqplot(rvbresults[rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",]))
  expect_no_error(manhattan(rvbresults[rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",], label = "unit", contigs = "GRCh38"))
  genesetlist <- buildGeneSet(list("genesetA" = sample(rvbresults$unit, 50)))
  expect_no_error({
    suppressMessages(densityPlot(rvbresults[rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",], 
                "genesetA",
                genesetlist, 
                showMeans = TRUE, 
                title = "test"))
    })
  }
)


test_that("ACAT works", {
    rvbresults_ACAT <- suppressMessages(ACAT(
      rvbresults,
      aggregate = list("test", "varSetName"),
      fixpval = TRUE,
      fixpval_method = "Liu",
      fixpval_maxP = 0.99,
      fixpval_minP = 1e-16
    ))
    metadata(rvbresults_ACAT)$gdbPath <- NA_character_
    metadata(rvbresults_ACAT)$creationDate <- NA_character_
    expect_snapshot_value(rvbresults_ACAT, style = "serialize")
}
)

