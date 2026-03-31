test_that("varSetFile input validation works", {
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))

  # expect error if varSetFile is called with incorrect input
  expect_error(
    {
      varSetFile("")
    },
    regexp = "must be a single, non-empty file path string"
  )

  expect_error(
    {
      varSetFile(character(0))
    },
    regexp = "must be a single, non-empty file path string"
  )

  expect_error(
    {
      varSetFile("/nonexistent/file.txt")
    },
    regexp = "File does not exist"
  )

  # expect error if varSetFile has incorrect unit count
  tmpfile <- withr::local_tempfile()
  writeLines("unit1|1,2|1,1|test", tmpfile)
  varsetfile_corrupt <- new(
    "varSetFile",
    path = tmpfile,
    units = c("unit1", "unit2"),
    metadata = list()
  )
  expect_error(
    {
      listVarSets(varsetfile_corrupt)
    },
    regexp = "does not match the number of units cached"
  )

  # expect error if varSetFile has more header lines that expected
  tmpfile <- withr::local_tempfile()
  header <- paste(rep("# header line", 15), collapse = "\n")
  writeLines(c(header, "unit1|1,2|1,1|test"), tmpfile)
  varsetfile_incorrect_header <- new(
    "varSetFile",
    path = tmpfile,
    units = character(),
    metadata = list()
  )

  expect_error(
    {
      rvat:::.varsetfile_read_metadata(
        varsetfile_incorrect_header,
        return_header = FALSE
      )
    },
    regexp = "File contains more header lines than expected"
  )
})
