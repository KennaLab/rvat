test_that("getVarSet basic functionality works", {
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))

  # extract some units and varset
  varsets <- getVarSet(
    varsetfile,
    unit = c("FUS", "OPTN"),
    varSetName = "Moderate"
  )

  expect_s4_class(varsets, "varSetList")
  expect_equal(length(varsets), 2L)
  expect_identical(sort(listUnits(varsets)), c("FUS", "OPTN"))
  expect_identical(unique(listVarSets(varsets)), "Moderate")
})

test_that("getVarSet varSetFile vs varSetList produces identical results", {
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))

  # extract from varSetFile
  varsets_from_file <- getVarSet(
    varsetfile,
    unit = c("FUS", "OPTN"),
    varSetName = "Moderate"
  )

  # extract from varSetList
  varsetlist <- getVarSet(varsetfile, unit = listUnits(varsetfile))
  varsets_from_list <- getVarSet(
    varsetlist,
    unit = c("FUS", "OPTN"),
    varSetName = "Moderate"
  )

  expect_identical(varsets_from_file, varsets_from_list)
})

test_that("getVarSet input validation works", {
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
  varsetlist <- getVarSet(varsetfile, unit = listUnits(varsetfile))

  # expect warning when extracting varSets that are not present in varSet
  expect_warning(
    varsets <- getVarSet(
      varsetlist,
      unit = c("FUS", "OPTN"),
      varSetName = c("Moderate", "moderate")
    ),
    "Not all specified varSets are present in the varSetList"
  )

  # expect warning when extracting units that are not present in varSet
  expect_warning(
    varsets2 <- getVarSet(
      varsetlist,
      unit = c("FUS", "OPTN", "abc"),
      varSetName = "Moderate"
    ),
    "Not all specified units are present in the varSetList"
  )

  ## same for varsetfile
  expect_warning(
    varsets2 <- getVarSet(
      varsetfile,
      unit = c("FUS", "OPTN", "abc"),
      varSetName = "Moderate"
    ),
    "Not all specified units are present in the varSetFile"
  )


  ## should still return valid results for existing units
  expect_identical(varsets, varsets2)
  
  # at least one of unit or varSetName should be specified
  expect_error({
    varsets <- getVarSet(varsetlist)
  }, regexp = "At least one of")
  
  # expect error if varSetFile contains malformed lines
  tmpfile <- withr::local_tempfile()
  ## line with three fields
  writeLines("unit1|field1|field2", tmpfile) 
  varsetfile_malformed <- new("varSetFile", 
                        path = tmpfile, 
                        units = "unit1",
                        metadata = list())
  
  expect_error({
    getVarSet(varsetfile_malformed, unit = "unit1")
  }, regexp = "Malformed line encountered")

})