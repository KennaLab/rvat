gdb <- create_example_gdb()
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
varsetlist <- getVarSet(
  varsetfile,
  unit = unique(listUnits(varsetfile))[1:10],
  varSetName = "Moderate"
)
aggdb_file <- withr::local_tempfile(fileext = ".aggdb")
aggregate(
  x = gdb,
  varSet = varsetlist,
  maxMAF = 0.001,
  output = aggdb_file,
  verbose = FALSE,
  signif = 12
)
aggdb <- aggdb(aggdb_file)

test_that("mergeAggDbs/collapseAggDbs works", {
  # split varsetlist in two
  varsetlist1 <- getVarSet(varsetlist, unit = listUnits(varsetlist)[1:5])
  varsetlist2 <- getVarSet(varsetlist, unit = listUnits(varsetlist)[6:10])

  # write aggdbs
  aggdb1 <- withr::local_tempfile()
  aggdb2 <- withr::local_tempfile()
  aggregate(
    x = gdb,
    varSet = varsetlist1,
    maxMAF = 0.001,
    output = aggdb1,
    verbose = FALSE,
    signif = 12
  )
  aggregate(
    x = gdb,
    varSet = varsetlist2,
    maxMAF = 0.001,
    output = aggdb2,
    verbose = FALSE,
    signif = 12
  )

  # merge aggdbs
  aggdblist <- aggdbList(
    filelist = c(aggdb1, aggdb2)
  )
  expect_true(stringr::str_detect(
    capture_output({
      show(aggdblist)
    }),
    "aggdbList object"
  ))
  aggdb_merged <- withr::local_tempfile()
  suppressMessages(mergeAggDbs(
    aggdblist,
    output = aggdb_merged
  ))
  test1 <- getUnit(aggdb(aggdb_merged), unit = listUnits(aggdb(aggdb_merged)))
  test2 <- getUnit(aggdb, unit = listUnits(aggdb(aggdb_merged)))
  expect_equal(test1, test2)

  # collapse aggDbs
  varsetlist1 <- getVarSet(varsetlist, unit = listUnits(varsetlist)[1])
  varsetlist2 <- getVarSet(varsetlist, unit = listUnits(varsetlist)[2])

  aggdb1 <- withr::local_tempfile()
  aggdb2 <- withr::local_tempfile()
  aggregate(
    x = gdb,
    varSet = varsetlist1,
    maxMAF = 0.001,
    output = aggdb1,
    verbose = FALSE,
    signif = 12
  )
  aggregate(
    x = gdb,
    varSet = varsetlist2,
    maxMAF = 0.001,
    output = aggdb2,
    verbose = FALSE,
    signif = 12
  )
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
    varSet = collapseVarSetList(getVarSet(
      varsetlist,
      listUnits(varsetlist)[1:2]
    )),
    cohort = "pheno",
    verbose = FALSE
  )
  GT <- GT[getMAF(GT) < 0.001, ]
  GT <- aggregate(recode(GT, imputeMethod = "meanImpute"))
  collapsed_aggregate_check <- collapsed_aggregate %>%
    dplyr::left_join(
      tibble::tibble(IID = colnames(GT), aggregate = unname(GT$aggregate)),
      by = "IID"
    )
  expect_equal(
    collapsed_aggregate_check$aggregate.x,
    collapsed_aggregate_check$aggregate.y
  )

  ## check returning aggregates interactively
  collapsed_aggregate_interactive <- collapseAggDbs(
    aggdblist,
    verbose = FALSE
  )
  expect_equal(
    collapsed_aggregate_interactive,
    as.data.frame(collapsed_aggregate)
  )
})


test_that("mergeAggDbs input validation works", {
  # expect error when merge fails
  varsetlist1 <- getVarSet(varsetlist, unit = listUnits(varsetlist)[1])
  varsetlist2 <- getVarSet(varsetlist, unit = listUnits(varsetlist)[2])

  aggdb1 <- withr::local_tempfile()
  aggdb2 <- withr::local_tempfile()
  aggregate(
    x = gdb,
    varSet = varsetlist1,
    maxMAF = 0.001,
    output = aggdb1,
    verbose = FALSE,
    signif = 12
  )
  aggregate(
    x = gdb,
    varSet = varsetlist2,
    maxMAF = 0.001,
    output = aggdb2,
    verbose = FALSE,
    signif = 12
  )

  con <- DBI::dbConnect(RSQLite::SQLite(), aggdb1)
  ## remove the 'aggregates' table, which should lead to an error
  aggdblist <- aggdbList(
    filelist = c(aggdb1, aggdb2)
  )
  DBI::dbRemoveTable(con, "aggregates")
  aggdb_merged <- withr::local_tempfile()
  expect_error(
    {
      suppressMessages(mergeAggDbs(
        aggdblist,
        output = aggdb_merged
      ))
    },
    regexp = "Failed to merge data from"
  )

  # expect error when output is invalid
  expect_error(
    {
      suppressMessages(mergeAggDbs(
        aggdblist,
        output = TRUE
      ))
    },
    regexp = "`output` must be"
  )
})

test_that("collapseAggDbs input validation works", {
  # expect error when merge fails
  varsetlist1 <- getVarSet(varsetlist, unit = listUnits(varsetlist)[1])
  varsetlist2 <- getVarSet(varsetlist, unit = listUnits(varsetlist)[2])

  aggdb1 <- withr::local_tempfile()
  aggdb2 <- withr::local_tempfile()
  aggregate(
    x = gdb,
    varSet = varsetlist1,
    maxMAF = 0.001,
    output = aggdb1,
    verbose = FALSE,
    signif = 12
  )
  aggregate(
    x = gdb,
    varSet = varsetlist2,
    maxMAF = 0.001,
    output = aggdb2,
    verbose = FALSE,
    signif = 12
  )

  con <- DBI::dbConnect(RSQLite::SQLite(), aggdb1)
  ## remove the 'aggregates' table, which should lead to an error
  aggdblist <- aggdbList(
    filelist = c(aggdb1, aggdb2)
  )
  DBI::dbRemoveTable(con, "aggregates")
  aggdb_collapsed <- withr::local_tempfile()
  expect_error(
    {
      suppressMessages(collapseAggDbs(
        aggdblist,
        output = aggdb_collapsed
      ))
    },
    regexp = "Could not process"
  )

  # expect error when output is invalid
  expect_error(
    {
      suppressMessages(collapseAggDbs(
        aggdblist,
        output = TRUE
      ))
    },
    regexp = "`output` must be a character"
  )
})
