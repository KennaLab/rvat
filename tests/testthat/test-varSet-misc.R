varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))

test_that("showing a varSetFile works", {
  expect_true(stringr::str_detect(
    capture_output({
      show(varsetfile)
    }),
    "varSetFile object"
  ))
})

test_that("showing a varSetList works", {
  expect_true(stringr::str_detect(
    capture_output({
      show(getVarSet(varsetfile, unit = "SOD1"))
    }),
    "varSetList"
  ))
})

test_that("showing a varSet works", {
  expect_true(stringr::str_detect(
    capture_output({
      show(getVarSet(varsetfile, unit = "SOD1")[[1]])
    }),
    "unit="
  ))
})

test_that("varSet methods are backward compatible with RVAT versions <=0.2.10", {
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))

  # write varsetfile without header
  varsetlist <- as.data.frame(getVarSet(
    varsetfile,
    unit = listUnits(varsetfile)
  ))
  path_varsetfile_noheader <- withr::local_tempfile()
  write.table(
    varsetlist,
    file = path_varsetfile_noheader,
    col.names = FALSE,
    sep = "|",
    append = FALSE,
    quote = FALSE,
    row.names = FALSE
  )

  ## connecting should work but issue a warning
  expect_warning(
    {
      varsetfile_noheader <- varSetFile(path_varsetfile_noheader)
    },
    regexp = "The following metadata fields are missing"
  )
  ## compare
  compare_varsetfile(varsetfile, varsetfile_noheader)

  ## retrieve varsets
  varsets <- getVarSet(
    varsetfile_noheader,
    unit = listUnits(varsetfile_noheader)[c(3, 1)]
  )

  # write varsetfile with unexpected metadata
  path_varsetfile_unexpected_meta <- withr::local_tempfile()
  con <- gzfile(path_varsetfile_unexpected_meta, "w")
  metadata <- list(
    test = "test"
  )
  rvat:::.write_rvat_header(
    filetype = "varSetFile",
    metadata = metadata,
    con = con
  )
  write.table(
    varsetlist,
    file = con,
    col.names = FALSE,
    sep = "|",
    append = TRUE,
    quote = FALSE,
    row.names = FALSE
  )
  close(con)

  ## expect error upon connecting
  expect_error(
    {
      vfile <- varSetFile(path_varsetfile_unexpected_meta)
    },
    regexp = "unexpected metadata fields were found: test"
  )

  # write varsetfile with unexpected metadata (2)
  con <- gzfile(path_varsetfile_unexpected_meta, "w")
  for (item in names(metadata)) {
    writeLines(
      sprintf("# %s: %s: %s", item, metadata[[item]], metadata[[item]]),
      con = con
    )
  }
  write.table(
    varsetlist,
    file = con,
    col.names = FALSE,
    sep = "|",
    append = TRUE,
    quote = FALSE,
    row.names = FALSE
  )
  close(con)

  ## expect error upon connecting
  expect_error(
    {
      vfile <- varSetFile(path_varsetfile_unexpected_meta)
    },
    regexp = "Unexpected filetype"
  )
})
