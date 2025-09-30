# check summariseGeno
test_that("summariseGeno works", {
  # compare summariseGeno-gdb and summariseGeno-genoMatrix methods
  gdb <- gdb(rvat_example("rvatData.gdb"))
  varids <- getAnno(gdb, "var", fields = "VAR_id")$VAR_id[1:500]
  GT <- getGT(
    gdb,
    VAR_id = varids,
    cohort = "pheno",
    verbose = FALSE
  )
  sumgeno <- summariseGeno(GT)
  sumgeno_output <- withr::local_tempfile()
  suppressMessages(summariseGeno(
    gdb,
    varSet = rvat:::.varsTovarSetList(varids),
    cohort = "pheno",
    output = sumgeno_output
  ))
  sumgeno_fromgdb <- as.data.frame(readr::read_tsv(
    sumgeno_output,
    show_col_types = FALSE
  ))
  rownames(sumgeno) <- NULL
  rownames(sumgeno_fromgdb) <- NULL
  sumgeno$VAR_id <- as.character(sumgeno$VAR_id)
  sumgeno_fromgdb$VAR_id <- as.character(sumgeno_fromgdb$VAR_id)
  expect_equal(sumgeno, sumgeno_fromgdb)

  # check `splitBy` argument
  suppressMessages(summariseGeno(
    gdb,
    varSet = rvat:::.varsTovarSetList(varids),
    cohort = "pheno",
    splitBy = "pheno",
    output = sumgeno_output
  ))
  sumgeno_fromgdb <- as.data.frame(readr::read_tsv(
    sumgeno_output,
    show_col_types = FALSE
  ))
  sumgeno <-
    rbind(
      summariseGeno(GT[, GT$pheno == 1]) %>%
        dplyr::mutate(pheno = 1) %>%
        dplyr::select(VAR_id, pheno, dplyr::everything()),
      summariseGeno(GT[, GT$pheno == 0]) %>%
        dplyr::mutate(pheno = 0) %>%
        dplyr::select(VAR_id, pheno, dplyr::everything())
    )
  rownames(sumgeno) <- NULL
  rownames(sumgeno_fromgdb) <- NULL
  sumgeno$VAR_id <- as.character(sumgeno$VAR_id)
  sumgeno_fromgdb$VAR_id <- as.character(sumgeno_fromgdb$VAR_id)
  expect_equal(sumgeno, sumgeno_fromgdb)

  # check runnin summariseGeno on a varSetFile
  varsetfile_path <- withr::local_tempfile()
  varsetfile <- buildVarSet(
    object = gdb,
    output = varsetfile_path,
    varSetName = "Moderate",
    unitTable = "varInfo",
    unitName = "gene_name",
    where = "(ModerateImpact = 1 or HighImpact = 1) and gene_name in ('CYP19A1', 'FUS', 'OPTN', 'SOD1')",
    verbose = FALSE
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
  sumgeno <- summariseGeno(GT)
  summariseGeno(
    gdb,
    varSet = varSetFile(varsetfile_path),
    cohort = "pheno",
    output = sumgeno_output,
    verbose = FALSE
  )
  sumgeno_fromgdb <- as.data.frame(readr::read_tsv(
    sumgeno_output,
    show_col_types = FALSE
  ))
  sumgeno$VAR_id <- as.character(sumgeno$VAR_id)
  sumgeno_fromgdb$VAR_id <- as.character(sumgeno_fromgdb$VAR_id)
  sumgeno_fromgdb <- sumgeno_fromgdb[
    match(sumgeno$VAR_id, sumgeno_fromgdb$VAR_id),
  ]
  rownames(sumgeno) <- NULL
  rownames(sumgeno_fromgdb) <- NULL
  expect_equal(sumgeno, sumgeno_fromgdb)

  # check whether strict testing works by supplying a different gdb
  expect_error(
    {
      summariseGeno(
        gdb(testgdb),
        varSet = varSetFile(varsetfile_path),
        cohort = "SM",
        output = sumgeno_output,
        verbose = FALSE
      )
    },
    regexp = "The varSetFile seems to be generated from a different gdb than supplied."
  )

  ## when setting strict = FALSE, no error should be raised
  expect_no_error({
    summariseGeno(
      gdb,
      varSet = varSetFile(varsetfile_path),
      cohort = "SM",
      output = sumgeno_output,
      verbose = FALSE,
      strict = FALSE
    )
  })

  # also check for assocTest -> note: should move to separate test block
  vfile <- varSetFile(varsetfile_path)
  expect_error(
    {
      assocTest(
        gdb(testgdb),
        varSet = getVarSet(vfile, unit = listUnits(vfile)[1]),
        cohort = "SM",
        pheno = "sex",
        continuous = TRUE,
        test = "lm",
        output = sumgeno_output,
        verbose = FALSE
      )
    },
    regexp = "generated from a different gdb"
  )

  ## when setting strict = FALSE, no error should be raised
  expect_no_error({
    assocTest(
      gdb(testgdb),
      varSet = getVarSet(vfile, unit = listUnits(vfile)[1]),
      cohort = "SM",
      pheno = "sex",
      continuous = TRUE,
      test = "lm",
      output = sumgeno_output,
      verbose = FALSE,
      strict = FALSE
    )
  })
})
