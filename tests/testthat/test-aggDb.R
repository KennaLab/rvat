gdb <- create_example_gdb()
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
varsetlist <- getVarSet(
  varsetfile,
  unit = listUnits(varsetfile)[1:5],
  varSetName = "Moderate"
)
aggdb_file <- withr::local_tempfile(fileext = ".aggdb")
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

test_that("aggDb input validation works", {
  # expect error when aggdb doesn't exit
  expect_error(
    {
      aggdb("nonexistent.aggdb")
    },
    regexp = "doesn't exist"
  )

  # expect error when aggdb is invalid (e.g. a directory)
  expect_error(
    {
      aggdb(tempdir())
    },
    regexp = "Invalid aggdb path"
  )
  #
})

test_that("aggDb-getUnit input validation works", {
  # expect warning when not all specified units are found
  aggdb <- aggdb(aggdb_file)

  expect_warning(
    {
      aggregates <- getUnit(aggdb, "nonexistent")
    },
    regexp = "1/1 of specified units could not be found"
  )
})

test_that("aggDbList input validation works", {
  gdb <- create_example_gdb()
  # expecterror when duplicate aggdb files are supplied
  expect_error(
    {
      aggdb_list <- aggdbList(c(aggdb_file, aggdb_file))
    },
    regexp = "The following units are duplicated"
  )

  # expect error when sample lists are not identical

  ## generate aggdb on subset of samples
  aggfile_subset <- withr::local_tempfile(fileext = ".aggdb")
  pheno <- getCohort(gdb, "pheno")
  uploadCohort(gdb, pheno[1:100, ], name = "pheno_subset", verbose = FALSE)
  varsetlist2 <- getVarSet(
    varsetfile,
    unit = listUnits(varsetfile)[10:15],
    varSetName = "Moderate"
  )
  aggregate(
    x = gdb,
    varSet = varsetlist2,
    cohort = "pheno_subset",
    maxMAF = 0.001,
    output = aggfile_subset,
    verbose = FALSE,
    signif = 12,
    overWrite = TRUE
  )

  expect_error(
    {
      aggdb_list <- aggdbList(c(aggdb_file, aggfile_subset))
    },
    regexp = "should include identical samples"
  )

  # expect error when parameters are not identical across aggdbs
  aggdb1_file <- withr::local_tempfile(fileext = ".aggdb")
  aggregate(
    x = gdb,
    varSet = varsetlist[1],
    cohort = "pheno",
    maxMAF = 0.001,
    output = aggdb1_file,
    verbose = FALSE,
    signif = 12,
    overWrite = TRUE
  )

  aggdb2_file <- withr::local_tempfile(fileext = ".aggdb")
  aggregate(
    x = gdb,
    varSet = varsetlist[2],
    cohort = "pheno",
    maxMAF = 0.005,
    output = aggdb2_file,
    verbose = FALSE,
    signif = 12,
    overWrite = TRUE
  )
  expect_error(
    {
      aggdblist <- aggdbList(c(aggdb1_file, aggdb2_file))
    },
    regexp = "Non-identical parameters were used to generate the input aggdbs."
  )

  # expect error when aggdbs were generated using different RVAT versions
  aggdb1_file <- withr::local_tempfile(fileext = ".aggdb")
  aggregate(
    x = gdb,
    varSet = varsetlist[1],
    cohort = "pheno",
    maxMAF = 0.001,
    output = aggdb1_file,
    verbose = FALSE,
    signif = 12,
    overWrite = TRUE
  )

  aggdb2_file <- withr::local_tempfile(fileext = ".aggdb")
  aggregate(
    x = gdb,
    varSet = varsetlist[2],
    cohort = "pheno",
    maxMAF = 0.001,
    output = aggdb2_file,
    verbose = FALSE,
    signif = 12,
    overWrite = TRUE
  )
  aggdb1 <- aggdb(aggdb1_file)
  aggdb2 <- aggdb(aggdb2_file)
  DBI::dbExecute(
    aggdb1,
    "UPDATE meta SET value = :new_value WHERE name = :name",
    params = list(new_value = "hello", name = "rvatVersion")
  )

  expect_error(
    {
      aggdb_list <- aggdbList(c(aggdb1_file, aggdb2_file))
    },
    regexp = "different rvat versions"
  )

  # expect error when aggdbs were generated from different gdbs
  aggdb1_file <- withr::local_tempfile(fileext = ".aggdb")
  aggregate(
    x = gdb,
    varSet = varsetlist[1],
    cohort = "pheno",
    maxMAF = 0.001,
    output = aggdb1_file,
    verbose = FALSE,
    signif = 12,
    overWrite = TRUE
  )

  aggdb2_file <- withr::local_tempfile(fileext = ".aggdb")
  aggregate(
    x = gdb,
    varSet = varsetlist[2],
    cohort = "pheno",
    maxMAF = 0.001,
    output = aggdb2_file,
    verbose = FALSE,
    signif = 12,
    overWrite = TRUE
  )
  aggdb1 <- aggdb(aggdb1_file)
  aggdb2 <- aggdb(aggdb2_file)
  DBI::dbExecute(
    aggdb1,
    "UPDATE meta SET value = :new_value WHERE name = :name",
    params = list(new_value = "hello", name = "gdbId")
  )

  expect_error(
    {
      aggdb_list <- aggdbList(c(aggdb1_file, aggdb2_file))
    },
    regexp = "different gdbs"
  )
})
