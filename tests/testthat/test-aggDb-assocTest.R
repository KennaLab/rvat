gdb <- create_example_gdb()
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
varsetlist <- getVarSet(
  varsetfile,
  unit = unique(listUnits(varsetfile))[1:10],
  varSetName = "Moderate"
)

genesetlist <- buildGeneSet(
  list(
    "geneset1" = listUnits(varsetlist)[c(1, 3, 4, 9)],
    "geneset2" = listUnits(varsetlist)[c(2, 4, 5, 6, 10)],
    "geneset3" = listUnits(varsetlist)[c(1, 2)],
    "geneset4" = listUnits(varsetfile)[c(7, 9, 10)],
    "geneset5" = listUnits(varsetlist)
  )
)

# generate aggregates
aggdb_file <- withr::local_tempfile(fileext = ".aggdb")
test_that("aggDb works", {
  expect_no_error(
    aggregate(
      x = gdb,
      varSet = varsetlist,
      maxMAF = 0.001,
      output = aggdb_file,
      verbose = FALSE,
      signif = 12
    )
  )
})

# compare assocTest-aggdb with assocTest-GT
test_that("assocTest-aggdb results identical to assocTest-gdb", {
  aggdb <- aggdb(aggdb_file)
  aggAssoc <- assocTest(
    aggdb,
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
  for (geneset in listGeneSets(genesetlist)) {
    test <- assocTest(
      object = gdb,
      test = c("glm", "firth"),
      varSet = collapseVarSetList(
        getVarSet(
          varsetlist,
          unit = listUnits(getGeneSet(genesetlist, geneSet = geneset))
        ),
        drop = FALSE
      ),
      cohort = "pheno",
      pheno = "pheno",
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
    as.data.frame(assoc)[, c(
      "meanCaseScore",
      "meanCtrlScore",
      "effect",
      "effectSE",
      "effectCIlower",
      "effectCIupper",
      "OR",
      "P"
    )],
    as.data.frame(aggAssoc)[, c(
      "meanCaseScore",
      "meanCtrlScore",
      "effect",
      "effectSE",
      "effectCIlower",
      "effectCIupper",
      "OR",
      "P"
    )]
  )

  # check continuous results
  aggAssoc_cont <- assocTest(
    aggdb,
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
  for (geneset in listGeneSets(genesetlist)) {
    test <- assocTest(
      object = gdb,
      test = c("lm"),
      varSet = collapseVarSetList(
        getVarSet(
          varsetlist,
          unit = listUnits(getGeneSet(genesetlist, geneSet = geneset))
        ),
        drop = FALSE
      ),
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
    as.data.frame(assoc_cont)[, c(
      "meanCaseScore",
      "meanCtrlScore",
      "effect",
      "effectSE",
      "effectCIlower",
      "effectCIupper",
      "P"
    )],
    as.data.frame(aggAssoc_cont)[, c(
      "meanCaseScore",
      "meanCtrlScore",
      "effect",
      "effectSE",
      "effectCIlower",
      "effectCIupper",
      "P"
    )]
  )

  # check dummy covariates
  aggAssoc_dummy <- assocTest(
    aggdb,
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
  for (geneset in listGeneSets(genesetlist)[1:2]) {
    test <- assocTest(
      object = gdb,
      test = c("lm"),
      varSet = collapseVarSetList(
        getVarSet(
          varsetlist,
          unit = listUnits(getGeneSet(genesetlist, geneSet = geneset))
        ),
        drop = FALSE
      ),
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
    as.data.frame(assoc_dummy)[, c(
      "meanCaseScore",
      "meanCtrlScore",
      "effect",
      "effectSE",
      "effectCIlower",
      "effectCIupper",
      "P"
    )],
    as.data.frame(aggAssoc_dummy)[, c(
      "meanCaseScore",
      "meanCtrlScore",
      "effect",
      "effectSE",
      "effectCIlower",
      "effectCIupper",
      "P"
    )]
  )

  # check dropUnits
  aggAssoc_dropunits <- assocTest(
    aggdb,
    gdb = gdb,
    test = "glm",
    cohort = "pheno",
    pheno = "pheno",
    dropUnits = listUnits(genesetlist[[1]])[1:2],
    geneSet = genesetlist[1:2],
    covar = paste0("PC", 1:4),
    verbose = FALSE,
  )
  assoc_dropunits <- list()
  for (geneset in listGeneSets(genesetlist)[1:2]) {
    units <- listUnits(getGeneSet(genesetlist, geneSet = geneset))
    units <- units[!units %in% listUnits(genesetlist[[1]])[1:2]]
    test <- assocTest(
      object = gdb,
      test = "glm",
      varSet = collapseVarSetList(
        getVarSet(
          varsetlist,
          unit = units
        ),
        drop = FALSE
      ),
      cohort = "pheno",
      pheno = "pheno",
      maxMAF = 0.001,
      covar = paste0("PC", 1:4),
      verbose = FALSE
    )
    test$unit <- geneset
    assoc_dropunits[[geneset]] <- test
  }
  assoc_dropunits <- do.call(rbind, assoc_dropunits)
  expect_equal(
    as.data.frame(assoc_dropunits)[, c(
      "meanCaseScore",
      "meanCtrlScore",
      "effect",
      "effectSE",
      "effectCIlower",
      "effectCIupper",
      "P"
    )],
    as.data.frame(aggAssoc_dropunits)[, c(
      "meanCaseScore",
      "meanCtrlScore",
      "effect",
      "effectSE",
      "effectCIlower",
      "effectCIupper",
      "P"
    )]
  )

  # check specifying multiple covariates/phenotypes
  aggAssoc_multiple_params <- assocTest(
    aggdb,
    gdb = gdb,
    test = "lm",
    cohort = "pheno",
    pheno = c("PC1", "PC2"),
    continuous = TRUE,
    geneSet = genesetlist[1:2],
    covar = list(paste0("PC", 3:4), c("PC3", "superPop")),
    verbose = FALSE,
  )

  ## compare with running separately
  for (covar in list(paste0("PC", 3:4), c("PC3", "superPop"))) {
    for (pheno in c("PC1", "PC2")) {
      aggAssoc_single_param <- assocTest(
        aggdb,
        gdb = gdb,
        test = "lm",
        cohort = "pheno",
        pheno = pheno,
        continuous = TRUE,
        geneSet = genesetlist[1:2],
        covar = covar,
        verbose = FALSE,
      )
      subset <- aggAssoc_multiple_params[
        aggAssoc_multiple_params$pheno == pheno &
          aggAssoc_multiple_params$covar == aggAssoc_single_param$covar[1],
      ]
      rownames(subset) <- NULL
      expect_equal(
        as.data.frame(subset)[, c(
          "meanCaseScore",
          "meanCtrlScore",
          "effect",
          "effectSE",
          "effectCIlower",
          "effectCIupper",
          "P"
        )],
        as.data.frame(aggAssoc_single_param)[, c(
          "meanCaseScore",
          "meanCtrlScore",
          "effect",
          "effectSE",
          "effectCIlower",
          "effectCIupper",
          "P"
        )]
      )
    }
  }

  # check consistency with weighted varsets
  varsetlist_cadd <- getVarSet(
    varsetfile,
    unit = unique(listUnits(varsetfile))[1:10],
    varSetName = "CADD"
  )
  aggdb_file_cadd_mb <- withr::local_tempfile(fileext = ".aggdb")
  suppressWarnings(aggregate(
    x = gdb,
    varSet = varsetlist_cadd,
    maxMAF = 0.001,
    MAFweights = "mb",
    output = aggdb_file_cadd_mb,
    verbose = FALSE,
    signif = 12,
    overWrite = TRUE
  ))
  aggdb_cadd_mb <- aggdb(aggdb_file_cadd_mb)
  aggAssoc_cadd_mb <- assocTest(
    aggdb_cadd_mb,
    gdb = gdb,
    test = "glm",
    cohort = "pheno",
    pheno = "pheno",
    geneSet = genesetlist[1:2],
    covar = paste0("PC", 1:4),
    verbose = FALSE,
  )

  assoc_cadd_mb <- list()
  for (geneset in listGeneSets(genesetlist)[1:2]) {
    varsetlist_cadd_geneset <- getVarSet(
      varsetlist_cadd,
      unit = listUnits(getGeneSet(genesetlist, geneSet = geneset))
    )
    VAR_id <- unlist(lapply(varsetlist_cadd_geneset@varSets, listVars))
    weight <- unlist(lapply(varsetlist_cadd_geneset@varSets, listWeights))
    weight <- ifelse(is.na(weight), ".", weight)
    varset <- varSetList(list(varSet(
      unit = "unnamed",
      varSet = "CADD",
      VAR_id = paste(VAR_id, collapse = ","),
      w = paste(weight, collapse = ",")
    )))

    expect_warning(
      {
        test <- assocTest(
          object = gdb,
          test = "glm",
          varSet = varset,
          cohort = "pheno",
          pheno = "pheno",
          maxMAF = 0.001,
          MAFweight = "mb",
          covar = paste0("PC", 1:4),
          verbose = FALSE
        )
      }
    )
    test$unit <- geneset
    assoc_cadd_mb[[geneset]] <- test
  }
  assoc_cadd_mb <- do.call(rbind, assoc_cadd_mb)
  expect_equal(
    as.data.frame(assoc_cadd_mb)[, c(
      "meanCaseScore",
      "meanCtrlScore",
      "effect",
      "effectSE",
      "effectCIlower",
      "effectCIupper",
      "P"
    )],
    as.data.frame(aggAssoc_cadd_mb)[, c(
      "meanCaseScore",
      "meanCtrlScore",
      "effect",
      "effectSE",
      "effectCIlower",
      "effectCIupper",
      "P"
    )]
  )
})

test_that("assocTest-aggdb misc", {
  # check if writing to disk results in identical results
  aggdb <- aggdb(aggdb_file)
  output <- withr::local_tempfile()
  aggAssoc <- assocTest(
    aggdb,
    gdb = gdb,
    test = c("glm", "firth"),
    cohort = "pheno",
    pheno = "pheno",
    geneSet = genesetlist[1:2],
    covar = paste0("PC", 1:4),
    verbose = FALSE
  )
  aggAssoc_output <- assocTest(
    aggdb,
    gdb = gdb,
    test = c("glm", "firth"),
    cohort = "pheno",
    pheno = "pheno",
    geneSet = genesetlist[1:2],
    covar = paste0("PC", 1:4),
    verbose = FALSE,
    output = output
  )
  aggAssoc_output <- readr::read_tsv(output, show_col_types = FALSE)
  expect_equal(
    as.data.frame(aggAssoc_output),
    as.data.frame(aggAssoc),
    tolerance = 1e-6
  )

  # check keeping subset of samples
  expect_warning(
    {
      aggAssoc <- suppressMessages(assocTest(
        aggdb,
        gdb = gdb,
        test = c("glm"),
        cohort = "pheno",
        pheno = "pheno",
        geneSet = genesetlist[1],
        covar = paste0("PC", 1:4),
        keep = listSamples(aggdb)[1:15000],
        verbose = FALSE
      ))
    },
    "15000/25000 samples in the aggdb are present in the cohort"
  )
  ## check if caseN/ctrlN match
  expect_identical(unique(aggAssoc$caseN), 5000L)
  expect_identical(unique(aggAssoc$ctrlN), 10000L)

  # check expected messages when verbose=TRUE
  messages <- capture_messages({
    aggAssoc <- assocTest(
      aggdb,
      gdb = gdb,
      test = c("glm"),
      cohort = "pheno",
      pheno = "pheno",
      geneSet = genesetlist[1],
      covar = paste0("PC", 1:4),
      verbose = TRUE
    )
  })
  expect_equal(messages[1], "Analysing geneset1\n")
  expect_equal(
    messages[2],
    "   4/4 units in the geneSet are present in the aggdb\n"
  )
  expect_equal(
    messages[3],
    "   Analysing pheno: pheno / covar: PC1,PC2,PC3,PC4\n"
  )
})

test_that("assocTest-aggdb input validation works", {
  aggdb <- aggdb(aggdb_file)

  # expect error when invalid test is specified
  expect_error(
    {
      assocTest(
        aggdb,
        gdb = gdb,
        test = "test",
        cohort = "pheno",
        pheno = "pheno",
        geneSet = genesetlist,
        covar = paste0("PC", 1:4),
        verbose = FALSE
      )
    },
    regexp = "The following tests are not valid"
  )

  # expect warning when tests are excluded (binary)
  expect_warning(
    {
      assocTest(
        aggdb,
        gdb = gdb,
        test = c("glm", "lm"),
        cohort = "pheno",
        pheno = "pheno",
        geneSet = genesetlist[1],
        covar = paste0("PC", 1:4),
        verbose = FALSE
      )
    },
    regexp = "The following tests were excluded since they are not available for binary"
  )

  # expect warning when tests are excluded (continuous)
  expect_warning(
    {
      assocTest(
        aggdb,
        gdb = gdb,
        test = c("glm", "lm"),
        cohort = "pheno",
        pheno = "PC1",
        continuous = TRUE,
        geneSet = genesetlist[1],
        covar = paste0("PC", 2:4),
        verbose = FALSE
      )
    },
    regexp = "The following tests were excluded since they are not available for continuous"
  )

  # expect error when covar is neither a character vector or a list
  expect_error(
    {
      assocTest(
        aggdb,
        gdb = gdb,
        test = "firth",
        cohort = "pheno",
        pheno = "pheno",
        geneSet = genesetlist,
        covar = 1:4,
        verbose = FALSE
      )
    },
    regexp = "`covar` must be either a character vector or a list"
  )

  # expect error when gdb is not a gdb object
  expect_error(
    {
      assocTest(
        aggdb,
        gdb = "test",
        test = "firth",
        cohort = "pheno",
        pheno = "pheno",
        geneSet = genesetlist,
        covar = paste0("PC", 1:4),
        verbose = FALSE
      )
    },
    regexp = "`gdb` must be a valid gdb object."
  )

  # expect error when non-existing phenotypes are specified
  expect_error(
    {
      assocTest(
        aggdb,
        gdb = gdb,
        test = "firth",
        cohort = "pheno",
        pheno = "test",
        geneSet = genesetlist,
        covar = paste0("PC", 1:4),
        verbose = FALSE
      )
    },
    regexp = "The `pheno` field"
  )

  # expect error when non-binary phenotype is specified and continuous != TRUE
  expect_error(
    {
      assocTest(
        aggdb,
        gdb = gdb,
        test = "firth",
        cohort = "pheno",
        pheno = "PC1",
        geneSet = genesetlist,
        covar = paste0("PC", 1:4),
        verbose = FALSE
      )
    },
    regexp = "Binary phenotypes should be coded"
  )

  # expect error when phenotype is not numeric
  expect_error(
    {
      assocTest(
        aggdb,
        gdb = gdb,
        test = "firth",
        cohort = "pheno",
        pheno = "superPop",
        geneSet = genesetlist,
        covar = paste0("PC", 1:4),
        verbose = FALSE
      )
    },
    regexp = "The following phenotypes are not numeric"
  )

  # expect error when non-existing cohort is specified
  expect_error(
    {
      assocTest(
        aggdb,
        gdb = gdb,
        test = "firth",
        cohort = "test",
        pheno = "pheno",
        geneSet = genesetlist,
        covar = paste0("PC", 1:4),
        verbose = FALSE
      )
    },
    regexp = "The specified cohort"
  )

  # expect warning when not all samples in aggdb are present in cohort
  pheno <- getCohort(gdb, "pheno")
  pheno_subset <- pheno[1:10000, ]
  uploadCohort(
    gdb,
    "pheno_subset",
    pheno_subset,
    overWrite = TRUE,
    verbose = FALSE
  )
  expect_warning(
    {
      test_subset <- assocTest(
        aggdb,
        gdb = gdb,
        test = "firth",
        cohort = "pheno_subset",
        pheno = "pheno",
        geneSet = genesetlist,
        covar = paste0("PC", 1:4),
        verbose = FALSE
      )
    },
    regexp = "10000/25000 samples in the aggdb"
  )
  expect_equal(max(test_subset$caseN), sum(pheno_subset$pheno == 1))
  expect_equal(max(test_subset$ctrlN), sum(pheno_subset$pheno == 0))

  ## this should be equal to providing a keep-list
  messages <- capture_messages({
    test_subset_keep <- suppressWarnings(assocTest(
      aggdb,
      gdb = gdb,
      test = "firth",
      cohort = "pheno",
      pheno = "pheno",
      keep = pheno_subset$IID,
      geneSet = genesetlist,
      covar = paste0("PC", 1:4),
      verbose = TRUE
    ))
  })
  expect_true(any(stringr::str_detect(messages, "Keeping 10000/25000 samples")))
  test_subset_keep$cohort = "pheno_subset"
  expect_equal(test_subset, test_subset_keep)

  # expect warning when not all in samples in cohort are present in aggdb
  aggfile_subset <- withr::local_tempfile()
  aggregate(
    x = gdb,
    varSet = varsetlist[1],
    cohort = "pheno_subset",
    maxMAF = 0.001,
    output = aggfile_subset,
    verbose = FALSE,
    signif = 12
  )
  aggdb_subset <- aggdb(aggfile_subset)
  expect_warning(
    {
      test_subset <- assocTest(
        aggdb_subset,
        gdb = gdb,
        test = "firth",
        cohort = "pheno",
        pheno = "pheno",
        geneSet = genesetlist,
        covar = paste0("PC", 1:4),
        verbose = FALSE
      )
    },
    regexp = "10000/25000 samples in the cohort"
  )

  # empty results when no samples overlap
  pheno_subset2 <- pheno[10001:20000, ]
  uploadCohort(gdb, "pheno_subset2", pheno_subset2, verbose = FALSE)
  aggdb_subset <- aggdb(aggfile_subset)
  test_subset <- suppressWarnings(assocTest(
    aggdb_subset,
    gdb = gdb,
    test = "firth",
    cohort = "pheno_subset2",
    pheno = "pheno",
    geneSet = genesetlist,
    covar = paste0("PC", 1:4),
    verbose = FALSE
  ))
  expect_true(nrow(test_subset) == 0L)

  # expect error when non-existing covariates are specified
  expect_error(
    {
      assocTest(
        aggdb,
        gdb = gdb,
        test = "firth",
        cohort = "pheno",
        pheno = "pheno",
        geneSet = genesetlist,
        covar = "test",
        verbose = FALSE
      )
    },
    regexp = "The following `covar` fields were not found"
  )

  ## same for subtractCovar
  expect_error(
    {
      assocTest(
        aggdb,
        gdb = gdb,
        test = "firth",
        cohort = "pheno",
        pheno = "pheno",
        geneSet = genesetlist,
        subtractCovar = "test",
        verbose = FALSE
      )
    },
    regexp = "The `subtractCovar` field"
  )

  # expect error when geneSet is misspecified
  expect_error(
    {
      assocTest(
        aggdb,
        gdb = gdb,
        test = "firth",
        cohort = "pheno",
        pheno = "pheno",
        geneSet = "test",
        covar = paste0("PC", 1:4),
        verbose = FALSE
      )
    },
    regexp = "`geneSet` must be a valid geneSetList or geneSetFile"
  )
})
