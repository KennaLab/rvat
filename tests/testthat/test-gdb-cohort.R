pheno_df <- readr::read_tsv(
  rvat_example("rvatData.pheno"),
  show_col_types = FALSE,
  progress = FALSE
)
pheno_df <- as.data.frame(pheno_df)

test_that("uploadCohort basic functionality works", {
  gdb <- create_example_gdb()

  # upload cohort from data frame, expect no error
  expect_no_error(uploadCohort(
    gdb,
    value = pheno_df,
    name = "pheno_df",
    verbose = FALSE
  ))

  # upload cohort from file, expect no error
  uploadCohort(
    gdb,
    value = rvat_example("rvatData.pheno"),
    name = "pheno_fromfile",
    verbose = FALSE
  )

  # compare the two uploaded cohorts
  pheno_retrieved <- getCohort(gdb, "pheno_df")
  pheno_fromfile_retrieved <- getCohort(gdb, "pheno_fromfile")

  expect_equal(pheno_df, pheno_retrieved)
  expect_equal(pheno_df, pheno_fromfile_retrieved)

  # check messages
  suppressMessages(expect_message(
    uploadCohort(
      gdb,
      value = pheno_df,
      name = "check_messages",
      verbose = TRUE
    ),
    regexp = "Loading cohort"
  ))

  suppressMessages(expect_message(
    uploadCohort(
      gdb,
      value = rvat_example("rvatData.pheno"),
      name = "check_messages_file",
      verbose = TRUE
    ),
    regexp = "Loading cohort"
  ))
})


test_that("uploadCohort handles subset of IIDs correctly", {
  gdb <- create_example_gdb()

  # subset pheno
  pheno_subset_df <- pheno_df[c(1:1000, 5000:6000), ]

  # upload subsetted cohort from data.frame
  uploadCohort(
    gdb,
    value = pheno_subset_df,
    name = "pheno_subset_df",
    verbose = FALSE
  )

  # upload subset cohort from file
  pheno_subset_file <- withr::local_tempfile()
  readr::write_tsv(pheno_subset_df, file = pheno_subset_file)
  uploadCohort(
    gdb,
    value = pheno_subset_file,
    name = "pheno_subset_fromfile",
    verbose = FALSE
  )

  # compare cohorts
  pheno_subset_retrieved <- getCohort(gdb, "pheno_subset_df")
  pheno_subset_fromfile_retrieved <- getCohort(gdb, "pheno_subset_fromfile")

  # remove row names for comparison
  rownames(pheno_subset_df) <- NULL
  rownames(pheno_subset_retrieved) <- NULL
  rownames(pheno_subset_fromfile_retrieved) <- NULL

  expect_equal(pheno_subset_df, pheno_subset_retrieved)
  expect_equal(pheno_subset_df, pheno_subset_fromfile_retrieved)
})


test_that("uploadCohort warns about non-matching IIDs", {
  gdb <- create_example_gdb()
  # upload a cohort that contains IIDs not present in gdb
  pheno_nonmatchingiids <- rbind(
    pheno_df[c(1:1000, 5000:6000), ],
    pheno_df[1:100, ] %>% dplyr::mutate(IID = paste0("test", 1:100))
  )

  # uploading should warn about IIDs not present in gdb
  expect_warning(
    {
      uploadCohort(
        gdb,
        value = as.data.frame(pheno_nonmatchingiids),
        name = "pheno_nonmatching",
        verbose = FALSE
      )
    },
    regexp = "100/2101 samples in the cohort"
  )

  # retrieved data should contain matching IIDs
  pheno_nonmatching_retrieved <- getCohort(gdb, "pheno_nonmatching")
  rownames(pheno_nonmatching_retrieved) <- NULL
  pheno_compare <- pheno_df[c(1:1000, 5000:6000), ]
  rownames(pheno_compare) <- NULL
  expect_equal(pheno_compare, pheno_nonmatching_retrieved)
})

test_that("getCohort field subsetting works correctly", {
  gdb <- create_example_gdb()

  # test field subsetting
  pheno_subset <- getCohort(
    gdb,
    "pheno",
    fields = c("IID", "sex", "pheno", "pop")
  )
  expect_equal(colnames(pheno_subset), c("IID", "sex", "pheno", "pop"))

  # should error for non-existing fields
  expect_error(
    getCohort(gdb, "pheno", fields = c("IID", "sex", "pheno", "pop", "hello")),
  )
})

test_that("uploadCohort metadata tracking works correctly", {
  gdb <- create_example_gdb()
  original_cohorts <- listCohort(gdb)$name

  # upload cohort from interactive session
  uploadCohort(gdb, value = pheno_df, name = "meta_test", verbose = FALSE)
  cohorts_updated <- listCohort(gdb)

  expect_equal(setdiff(cohorts_updated$name, original_cohorts), "meta_test")
  expect_equal(
    cohorts_updated[cohorts_updated$name == "meta_test", ]$value,
    "interactive_session"
  )

  # should be back to original after dropping table
  dropTable(gdb, "meta_test", verbose = FALSE)
  expect_equal(sort(original_cohorts), sort(listCohort(gdb)$name))

  # upload file-based cohort
  uploadCohort(
    gdb,
    value = rvatData::rvat_example("rvatData.pheno"),
    name = "file_meta_test",
    verbose = FALSE
  )

  cohorts_updated <- listCohort(gdb)
  expect_true(stringr::str_detect(
    cohorts_updated[cohorts_updated$name == "file_meta_test", ]$value,
    "rvatData.pheno"
  ))
})


test_that("uploadCohort input validation works correctly", {
  gdb <- create_example_gdb()
  test_df <- data.frame(IID = c("ALS1", "ALS2", "ALS3"), sex = c(1, 1, 2))

  # expect error when duplicated IID values are included in upload
  test_df_dup <- data.frame(IID = c("ALS1", "ALS1", "ALS2"), sex = c(1, 1, 2))
  expect_error(
    uploadCohort(
      gdb,
      value = test_df_dup,
      name = "pheno_dups",
      verbose = FALSE
    ),
    regexp = "contains duplicated IIDs"
  )

  # expect error when upload name is protected
  for (name in gdb_protected_tables) {
    expect_error(
      uploadCohort(gdb, value = test_df, name = name),
      "already exists as a protected gdb table"
    )
  }

  # expect error when upload name exists and `overWrite = FALSE`
  expect_error(
    uploadCohort(gdb, value = test_df, name = "pheno"),
    "already in use"
  )

  # expect error when upload name already exists as an anno table
  expect_error(
    uploadCohort(gdb, value = test_df, name = "varinfo"),
    "already in use"
  )

  # expect error when upload name is invalid
  invalid_chars <- c(".", ",", "+", "-", " ")
  for (char in invalid_chars) {
    expect_error(
      uploadCohort(
        gdb,
        value = test_df,
        name = paste0("pheno", char, "test")
      ),
      "Table name may only contain"
    )
  }

  # expect error when IID column is missing
  expect_error(
    uploadCohort(
      gdb,
      value = test_df %>% dplyr::select(-IID),
      name = "noIID",
      verbose = FALSE
    ),
    regexp = "No 'IID' field detected"
  )

  # expect error when multiple cohorts are specified
  expect_error(
    uploadCohort(
      gdb,
      value = c(rvat_example("rvatData.pheno"), rvat_example("rvatData.pheno")),
      name = "dups",
      verbose = FALSE
    ),
    regexp = "must be a single file path"
  )

  # expect error when file doesn't exist
  tempfile <- withr::local_tempfile()
  expect_error(
    uploadCohort(
      gdb,
      value = tempfile,
      name = "nonexisting",
      verbose = FALSE
    ),
    regexp = "does not exist"
  )

  # expect error when sex column is missing
  expect_error(
    uploadCohort(
      gdb,
      value = test_df %>% dplyr::select(-sex),
      name = "noSex",
      verbose = FALSE
    ),
    regexp = "No 'sex' field detected"
  )

  # expect warning when sex contains non-numeric values
  expect_warning(
    uploadCohort(
      gdb,
      value = test_df %>% dplyr::mutate(sex = c("M", "M", "F")),
      name = "nonNumericSex",
      verbose = FALSE,
      overWrite = TRUE
    ),
    regexp = "`sex` should be coded as"
  )

  ## check if sex is set to unknown (0)
  values <- getCohort(gdb, "nonNumericSex", fields = "sex")$sex
  expect_true(all(values == 0))

  # expect warning when invalid sex codings are included
  expect_warning(
    uploadCohort(
      gdb,
      value = test_df %>% dplyr::mutate(sex = c(1, 2, 3)),
      name = "invalidSex",
      verbose = FALSE,
      overWrite = TRUE
    ),
    regexp = "1 values in `sex` column not coded as 0,1,2."
  )
  ## check if sex is set to unknown (0)
  values <- getCohort(gdb, "invalidSex", fields = "sex")$sex
  expect_identical(sort(unique(values)), c(0, 1, 2))
})
