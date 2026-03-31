# check summariseGeno
test_that("gdb-summariseGeno output identical genoMatrix-summariseGeno", {
  # compare summariseGeno-gdb and summariseGeno-genoMatrix methods
  gdb <- create_example_gdb()
  varids <- getAnno(gdb, "var", fields = "VAR_id")$VAR_id[1:500]
  GT <- getGT(
    gdb,
    VAR_id = varids,
    cohort = "pheno",
    verbose = FALSE
  )
  sumgeno_GT <- summariseGeno(GT)
  sumgeno_gdb <- suppressMessages(summariseGeno(
    gdb,
    varSet = rvat:::.varsTovarSetList(varids),
    cohort = "pheno"
  ))
  rownames(sumgeno_GT) <- NULL
  rownames(sumgeno_gdb) <- NULL
  expect_equal(sumgeno_GT, sumgeno_gdb)
})

test_that("summariseGeno splitBy works", {
  gdb <- create_example_gdb()
  varids <- getAnno(gdb, "var", fields = "VAR_id")$VAR_id[1:500]
  GT <- getGT(
    gdb,
    VAR_id = varids,
    cohort = "pheno",
    verbose = FALSE
  )

  # check `splitBy` argument, compare to manual split
  suppressMessages({
    sumgeno_splitby <- summariseGeno(
      gdb,
      varSet = rvat:::.varsTovarSetList(varids),
      cohort = "pheno",
      splitBy = "pheno",
    )
  })
  sumgeno_manual <-
    rbind(
      summariseGeno(GT[, GT$pheno == 1]) %>%
        dplyr::mutate(pheno = 1) %>%
        dplyr::select(VAR_id, pheno, dplyr::everything()),
      summariseGeno(GT[, GT$pheno == 0]) %>%
        dplyr::mutate(pheno = 0) %>%
        dplyr::select(VAR_id, pheno, dplyr::everything())
    )
  rownames(sumgeno_splitby) <- NULL
  rownames(sumgeno_manual) <- NULL
  expect_equal(sumgeno_splitby, sumgeno_manual)
})

test_that("gdb-summariseGeno works with varSetFile", {
  gdb <- create_example_gdb()

  # run summariseGeno on varSetFile, compare to running directly on genoMatrix
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
  varsetfile_path <- withr::local_tempfile()
  ## create a smaller varSetFile for testing
  write(
    getVarSet(
      varsetfile,
      unit = c("CYP19A1", "FUS", "OPTN", "SOD1"),
      varSetName = "ModerateHigh"
    ),
    varsetfile_path
  )
  varsetfile <- varSetFile(varsetfile_path)
  GT <- getGT(
    gdb,
    varSet = collapseVarSetList(getVarSet(
      varsetfile,
      unit = listUnits(varsetfile)
    )),
    cohort = "pheno",
    verbose = FALSE
  )
  sumgeno_GT <- summariseGeno(GT)
  sumgeno_gdb <- summariseGeno(
    gdb,
    varSet = varsetfile,
    cohort = "pheno",
    verbose = FALSE
  )
  sumgeno_gdb <- sumgeno_gdb[
    match(sumgeno_GT$VAR_id, sumgeno_gdb$VAR_id),
  ]
  rownames(sumgeno_GT) <- NULL
  rownames(sumgeno_gdb) <- NULL
  expect_equal(sumgeno_GT, sumgeno_gdb)
})

test_that("gdb-summariseGeno edge cases work", {
  gdb <- create_example_gdb()

  # no remaining variants
  sumgeno <- summariseGeno(
    gdb,
    VAR_id = 1:100,
    cohort = "pheno",
    minMAF = 0.001,
    maxMAF = 0.001,
    verbose = FALSE
  )
  expect_equal(sumgeno, NULL)
})

test_that("gdb-summariseGeno writing to output works", {
  gdb <- create_example_gdb()

  sumgeno_outfile <- withr::local_tempfile(fileext = ".txt.gz")
  summariseGeno(
    gdb,
    VAR_id = 1:100,
    cohort = "pheno",
    output = sumgeno_outfile,
    verbose = FALSE
  )
  sumgeno_fromfile <- readr::read_tsv(
    sumgeno_outfile,
    show_col_types = FALSE,
    col_types = list(VAR_id = readr::col_character())
  )
  sumgeno_interactive <- summariseGeno(
    gdb,
    VAR_id = 1:100,
    cohort = "pheno",
    verbose = FALSE
  )
  expect_equal(
    as.data.frame(sumgeno_fromfile),
    sumgeno_interactive
  )

  # also check with splitBy
  sumgeno_outfile <- withr::local_tempfile(fileext = ".txt.gz")
  summariseGeno(
    gdb,
    VAR_id = 1:100,
    cohort = "pheno",
    output = sumgeno_outfile,
    splitBy = "pheno",
    verbose = FALSE
  )
  sumgeno_fromfile <- readr::read_tsv(
    sumgeno_outfile,
    show_col_types = FALSE,
    col_types = list(VAR_id = readr::col_character())
  )
  sumgeno_interactive <- summariseGeno(
    gdb,
    VAR_id = 1:100,
    cohort = "pheno",
    splitBy = "pheno",
    verbose = FALSE
  )
  expect_equal(
    as.data.frame(sumgeno_fromfile),
    sumgeno_interactive
  )
})


test_that("gdb-summariseGeno input validation works", {
  gdb <- create_example_gdb()

  # expect error when supplying varSetFile from different gdb
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))

  ## build small test gdb
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
      summariseGeno(
        test_gdb,
        varSet = varsetfile,
        cohort = "SM",
        verbose = FALSE
      )
    },
    regexp = "The varSetFile seems to be generated from a different gdb than supplied."
  )

  ## when setting strict = FALSE, no error should be raised
  expect_no_error({
    suppressWarnings(summariseGeno(
      test_gdb,
      varSet = getVarSet(varsetfile, unit = listUnits(varsetfile)[1], varSetName = listVarSets(varsetfile)[1]),
      cohort = "SM",
      verbose = FALSE,
      strict = FALSE
    ))
  })

  # expect error when missing pheno is supplied
  expect_error(
    {
      summariseGeno(
        gdb,
        VAR_id = 1:5,
        cohort = "pheno",
        pheno = "a",
        verbose = FALSE
      )
    },
    regexp = "was not found"
  )

  # expect error when missing splitBy is supplied
  expect_error(
    {
      summariseGeno(
        gdb,
        VAR_id = 1:5,
        cohort = "pheno",
        splitBy = "a",
        verbose = FALSE
      )
    },
    regexp = "not found"
  )

  # expect error when no VAR_id/varSet is specified
  expect_error(
    {
      summariseGeno(
        gdb,
        cohort = "pheno",
        verbose = FALSE
      )
    },
    regexp = "Either of"
  )

  # expect error when varSet is of wrong type
  expect_error(
    {
      summariseGeno(
        gdb,
        cohort = "pheno",
        varSet = list(1, 2, 3),
        verbose = FALSE
      )
    },
    regexp = "should be of type"
  )
})
