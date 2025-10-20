gdb <- gdb(rvat_example("rvatData.gdb"))
moderate <- withr::local_tempfile()
varsetfile <- buildVarSet(
            object = gdb,
            output = moderate,
            varSetName = "Moderate",
            unitTable = "varInfo",
            unitName = "gene_name",
            where = "(ModerateImpact = 1 or HighImpact = 1)",
            verbose = FALSE
            )
varsetfile <- varSetFile(moderate)
aggdb <-  withr::local_tempfile()

# compare assocTest-aggregatefile with assocTest-GT
test_that("assocTest-aggregateFile works", {

  # generate aggregates
  expect_no_error(
    aggregate(x = gdb, 
              varSet = varsetfile,
              maxMAF = 0.001,
              output = aggdb,
              verbose = FALSE,
              signif = 12
  )
  )

  genesetlist <- buildGeneSet(
    list("geneset1" = listUnits(varsetfile)[c(1,3,4,9)],
         "geneset2" = listUnits(varsetfile)[c(2,4,5,6,10)],
         "geneset3" = listUnits(varsetfile)[c(1,2)],
         "geneset4" = listUnits(varsetfile)[c(7,9,10)],
         "geneset5" = listUnits(varsetfile)
    ))
  genesetfile <- withr::local_tempfile()
  write(genesetlist,
        file = genesetfile)
  aggAssoc <- assocTest(
    aggdb(aggdb),
    gdb = gdb,
    test = c("glm", "firth"),
    cohort = "pheno",
    pheno = "pheno",
    geneSet = genesetlist,
    covar = paste0("PC", 1:4),
    verbose = FALSE
  )

  # run assocTest on GT directly
  assoc <- list()
  for(geneset in listGeneSets(genesetlist)) {
    test <- assocTest(
      object = gdb,
      test = c("glm", "firth"),
      varSet = collapseVarSetList(getVarSet(varsetfile, unit = listUnits(getGeneSet(genesetlist, geneSet = geneset))),drop=FALSE),
      cohort = "pheno",
      pheno="pheno",
      maxMAF = 0.001,
      covar = paste0("PC", 1:4),
      verbose = FALSE
    )
    test$unit <- geneset
    assoc[[geneset]] <- test
  }
  assoc <- do.call(rbind, assoc)
  rownames(aggAssoc) <- NULL
  expect_equal(
    as.data.frame(assoc)[,c("meanCaseScore", "meanCtrlScore", "effect", "effectSE", "effectCIlower", "effectCIupper", "OR", "P")],
    as.data.frame(aggAssoc)[,c("meanCaseScore", "meanCtrlScore", "effect", "effectSE", "effectCIlower", "effectCIupper", "OR", "P")]
  )

  # check if writing to disk results in identical results
  output <- withr::local_tempfile()
  aggAssoc_output <- assocTest(
    aggdb(aggdb),
    gdb = gdb,
    test = c("glm", "firth"),
    cohort = "pheno",
    pheno = "pheno",
    geneSet = genesetlist,
    covar = paste0("PC", 1:4),
    verbose = FALSE,
    output = output
  )
  aggAssoc_output <- readr::read_tsv(output, show_col_types = FALSE)
  expect_equal(as.data.frame(aggAssoc_output), as.data.frame(aggAssoc), tolerance = 1e-6)

  # check keeping subset of samples
  expect_warning({aggAssoc <- suppressMessages(assocTest(
    aggdb(aggdb),
    gdb = gdb,
    test = c("glm"),
    cohort = "pheno",
    pheno="pheno",
    geneSet = genesetlist[1],
    covar = paste0("PC", 1:4),
    keep = listSamples(aggdb(aggdb))[1:15000],
    verbose = FALSE
  ))})
  ## check if caseN/ctrlN match
  expect_identical(unique(aggAssoc$caseN), 5000L)
  expect_identical(unique(aggAssoc$ctrlN), 10000L)

  # check continuous results
  aggAssoc_cont <- assocTest(
    aggdb(aggdb),
    gdb = gdb,
    test = "lm",
    cohort = "pheno",
    pheno = "PC1",
    geneSet = genesetlist,
    covar = paste0("PC", 2:4),
    verbose = FALSE,
    continuous = TRUE
  )
  assoc_cont <- list()
  for(geneset in listGeneSets(genesetlist)) {
    test <- assocTest(
      object = gdb,
      test = c("lm"),
      varSet = collapseVarSetList(getVarSet(varsetfile, unit = listUnits(getGeneSet(genesetlist, geneSet = geneset))),drop=FALSE),
      cohort = "pheno",
      pheno = "PC1",
      maxMAF = 0.001,
      covar = paste0("PC", 2:4),
      verbose = FALSE,
      continuous = TRUE
    )
    test$unit <- geneset
    assoc_cont[[geneset]] <- test
  }
    assoc_cont <- do.call(rbind, assoc_cont)
    expect_equal(
      as.data.frame(assoc_cont)[,c("meanCaseScore", "meanCtrlScore", "effect", "effectSE", "effectCIlower", "effectCIupper", "P")],
      as.data.frame(aggAssoc_cont)[,c("meanCaseScore", "meanCtrlScore", "effect", "effectSE", "effectCIlower", "effectCIupper", "P")]
    )
    
    # check expected messages when verbose=TRUE
    messages <- capture_messages({aggAssoc <- assocTest(
      aggdb(aggdb),
      gdb = gdb,
      test = c("glm"),
      cohort = "pheno",
      pheno="pheno",
      geneSet = genesetlist[1],
      covar = paste0("PC", 1:4),
      verbose = TRUE
    )})
    expect_equal(messages[1], "Analysing geneset1\n")
    expect_equal(messages[2], "   4/4 units in the geneSet are present in the aggdb\n")
    expect_equal(messages[3], "   Analysing pheno: pheno / covar: PC1,PC2,PC3,PC4\n")
    
    # units in the geneSet are present in the aggregateFile
    # check adding categorical variables
    aggAssoc_dummy <- assocTest(
      aggdb(aggdb),
      gdb = gdb,
      test = c("lm"),
      cohort = "pheno",
      pheno = "PC1",
      geneSet = genesetlist[1:2],
      covar = c(paste0("PC", 2:4), "superPop"),
      verbose = FALSE,
      continuous = TRUE
    )
    assoc_dummy <- list()
    for(geneset in listGeneSets(genesetlist)[1:2]) {
      test <- assocTest(
        object = gdb,
        test = c("lm"),
        varSet = collapseVarSetList(getVarSet(varsetfile, unit = listUnits(getGeneSet(genesetlist, geneSet = geneset))),drop=FALSE),
        cohort = "pheno",
        pheno = "PC1",
        maxMAF = 0.001,
        covar = c(paste0("PC", 2:4), "superPop"),
        verbose = FALSE,
        continuous = TRUE
      )
      test$unit <- geneset
      assoc_dummy[[geneset]] <- test
    }
    assoc_dummy <- do.call(rbind, assoc_dummy)
    expect_equal(
      as.data.frame(assoc_dummy)[,c("meanCaseScore", "meanCtrlScore", "effect", "effectSE", "effectCIlower", "effectCIupper", "P")],
      as.data.frame(aggAssoc_dummy)[,c("meanCaseScore", "meanCtrlScore", "effect", "effectSE", "effectCIlower", "effectCIupper", "P")]
    )
    
    # expect error when:
    
    ## invalid test is specified
    expect_error(
      {
        aggAssoc <- assocTest(
          aggdb(aggdb),
          gdb = gdb,
          test = "test",
          cohort = "pheno",
          pheno="pheno",
          geneSet = genesetlist,
          covar = paste0("PC", 1:4),
          verbose = FALSE
        )
      }
    )
    
    ## non-existing phenotypes are specified
    expect_error(
      {
        aggAssoc <- assocTest(
          aggdb(aggdb),
          gdb = gdb,
          test = "firth",
          cohort = "pheno",
          pheno = "test",
          geneSet = genesetlist,
          covar = paste0("PC", 1:4),
          verbose = FALSE
        )
      }
    )
    
    ## non-existing covariates are specified
    expect_error(
      {
        aggAssoc <- assocTest(
          aggdb(aggdb),
          gdb = gdb,
          test = "firth",
          cohort = "pheno",
          pheno = "pheno",
          geneSet = genesetlist,
          covar = "test",
          verbose = FALSE
        )
      }
    )
})


test_that("mergeAggDbs work", {
  # split varsetfiles in two
  varsetlist1 <- getVarSet(varsetfile, unit = listUnits(varsetfile)[1:6])
  varsetlist2 <- getVarSet(varsetfile, unit = listUnits(varsetfile)[7:12])
  varsetfile1 <- withr::local_tempfile()
  varsetfile2 <- withr::local_tempfile()
  write(varsetlist1, varsetfile1)
  write(varsetlist2, varsetfile2)
  aggdb1 <- withr::local_tempfile()
  aggdb2 <- withr::local_tempfile()
  aggregate(x = gdb,
            varSet = varSetFile(varsetfile1),
            maxMAF = 0.001,
            output = aggdb1,
            verbose = FALSE,
            signif = 12)
  aggregate(x = gdb,
            varSet = varSetFile(varsetfile2),
            maxMAF = 0.001,
            output = aggdb2,
            verbose = FALSE,
            signif = 12)

  # merge aggregateFiles
  aggdblist <- aggdbList(
    filelist = c(aggdb1, aggdb2)
  )
  expect_true(stringr::str_detect(capture_output({show(aggdblist)}), "aggdbList object"))
  aggdb_merged <- withr::local_tempfile()
  suppressMessages(mergeAggDbs(
    aggdblist,
    output = aggdb_merged
  ))
  test1 <- getUnit(aggdb(aggdb_merged), unit = listUnits(aggdb(aggdb_merged)))
  test2 <- getUnit(aggdb(aggdb), unit = listUnits(aggdb(aggdb_merged)))
  expect_equal(test1, test2)

  # merge aggregateFiles, collapse
  varsetlist1 <- getVarSet(varsetfile, unit = listUnits(varsetfile)[1])
  varsetlist2 <- getVarSet(varsetfile, unit = listUnits(varsetfile)[2])
  varsetfile1 <- withr::local_tempfile()
  varsetfile2 <- withr::local_tempfile()
  write(varsetlist1, varsetfile1)
  write(varsetlist2, varsetfile2)

  aggdb1 <- withr::local_tempfile()
  aggdb2 <- withr::local_tempfile()
  aggregate(x = gdb, 
            varSet = varSetFile(varsetfile1),
            maxMAF = 0.001,
            output = aggdb1,
            verbose = FALSE,
            signif = 12)
  aggregate(x = gdb, 
            varSet = varSetFile(varsetfile2),
            maxMAF = 0.001,
            output = aggdb2,
            verbose = FALSE,
            signif = 12)
  aggdblist <- aggdbList(
    filelist = c(aggdb1, aggdb2)
  )
  aggdb_merged <- withr::local_tempfile()
  suppressMessages(collapseAggDbs(
    aggdblist,
    output = aggdb_merged
  ))
  collapsed_aggregate <- readr::read_tsv(aggdb_merged, show_col_types = FALSE)

  ## compare with aggregate from GT directly
  GT <- getGT(
    gdb,
    varSet = collapseVarSetList(getVarSet(varsetfile, listUnits(varsetfile)[1:2])),
    cohort = "pheno",
    verbose = FALSE
  )
  GT <- GT[getMAF(GT) < 0.001, ]
  GT <- aggregate(recode(GT, imputeMethod = "meanImpute"))
  collapsed_aggregate <- collapsed_aggregate %>%
    dplyr::left_join(
      tibble::tibble(IID = colnames(GT), aggregate = unname(GT$aggregate)),
      by = "IID"
    )
  expect_equal(collapsed_aggregate$aggregate.x, collapsed_aggregate$aggregate.y)
})
