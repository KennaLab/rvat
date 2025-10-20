test_that("varSetList constructor works", {
  # expect varSetList of length 0 if invoked without arguments
  varsetlist <- varSetList()
  expect_true(is(varsetlist, "varSetList") && length(varsetlist) == 0L)

  # expect error if ..

  # data.frame is supplied with incorrect fields
  expect_error(
    {
      varsetlist <- varSetList(
        data.frame(
          unit = "SOD1",
          varSetName = "moderate",
          VARid = "1,2",
          w = "1,1"
        )
      )
    },
    regexp = "Invalid input"
  )

  # input is not a data.frame or list
  expect_error(
    {
      varsetlist <- varSetList(
        as.matrix(data.frame(
          unit = "SOD1",
          varSetName = "moderate",
          VARid = "1,2",
          w = "1,1"
        ))
      )
    },
    regexp = "Invalid input"
  )
})

test_that("writing a varsetList works", {
  # write varSetList to file and read it back in
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
  varsetlist <- getVarSet(varsetfile, unit = listUnits(varsetfile))
  tmpfile <- withr::local_tempfile()
  write(varsetlist, tmpfile)
  varsetfile_read <- varSetFile(tmpfile)
  compare_varsetfile(varsetfile, varsetfile_read)

  # expect error when writing to non-writable output
  tmpdir <- withr::local_tempdir()
  expect_error(
    {
      suppressWarnings(write(varsetlist, tmpdir))
    },
    regexp = "Failed to open file"
  )
})

test_that("misc varSetList methods work", {
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
  varsetlist <- getVarSet(varsetfile, unit = listUnits(varsetfile))
  varsets <- varsetlist[which(
    listUnits(varsetlist) %in%
      c("FUS", "OPTN") &
      listVarSets(varsetlist) %in% c("Moderate")
  )]
  expect_identical(sort(listUnits(varsets)), sort(c("FUS", "OPTN")))
  expect_identical(unique(listVarSets(varsets)), "Moderate")
})
