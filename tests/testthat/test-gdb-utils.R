test_that("subsetGdb and concatGdb roundtrip works", {
  # split gdb using subsetGdb and then concat using concatGdb -> original gdb and concatted gdb should be identical
  gdb <- create_example_gdb()
  vars <- getAnno(gdb, table = "var", fields = "VAR_id")$VAR_id
  vars <- split(
    vars,
    rep(1:3, each = ceiling(length(vars) / 3), length.out = length(vars))
  )
  tmpfiles <- c(
    withr::local_tempfile(),
    withr::local_tempfile(),
    withr::local_tempfile()
  )
  files <- withr::local_tempfile()
  readr::write_lines(tmpfiles, file = files)

  ## split in three parts
  for (i in 1:3) {
    suppressMessages(subsetGdb(
      gdb,
      output = tmpfiles[i],
      where = sprintf(
        "VAR_id in (%s)",
        paste(paste0("'", vars[[i]], "'"), collapse = ",")
      ),
      overWrite = TRUE
    ))
  }

  ## concat three parts
  gdbconcat <- withr::local_tempfile()
  suppressMessages(concatGdb(files, output = gdbconcat, verbose = TRUE))
  gdbconcat <- gdb(gdbconcat)

  ## compare gdbs
  compare_gdbs(gdb, gdbconcat, check_tables = FALSE)
})


test_that("subSetGdb works", {
  gdb <- create_example_gdb()

  # check if subsetting by SQL where clause works
  gdb_subset <- withr::local_tempfile()
  suppressMessages(subsetGdb(
    gdb,
    where = "CHROM = 'chr10'",
    output = gdb_subset,
    overWrite = TRUE
  ))
  gdb_subset <- gdb(gdb_subset)
  var1 <- getAnno(gdb_subset, "var")
  expect_equal(unique(var1$CHROM), "chr10")

  # check if subsetting by VAR_id works
  gdb_subset2 <- withr::local_tempfile()
  var <- getAnno(gdb, "var") %>%
    dplyr::filter(CHROM == "chr10")
  suppressMessages(subsetGdb(
    gdb,
    VAR_id = var$VAR_id,
    output = gdb_subset2,
    overWrite = TRUE
  ))
  gdb_subset2 <- gdb(gdb_subset2)
  ## check if identical to where clause subsetting above
  var2 <- getAnno(gdb_subset2, "var")
  rownames(var1) <- NULL
  var2 <- var2[match(var1$VAR_id, var2$VAR_id), ]
  rownames(var2) <- NULL
  expect_identical(var1, var2)

  # check if subsetting works
  gdb_subset3 <- withr::local_tempfile()
  var <- getAnno(gdb, "var") %>%
    dplyr::filter(CHROM == "chr10")
  uploadAnno(gdb, value = var, name = "var_chr10", verbose = FALSE)
  suppressMessages(subsetGdb(
    gdb,
    intersection = "var_chr10",
    output = gdb_subset3,
    overWrite = TRUE
  ))
  gdb_subset3 <- gdb(gdb_subset3)
  ## check if identical to where clause subsetting above
  var3 <- getAnno(gdb_subset3, "var")
  rownames(var1) <- NULL
  var3 <- var3[match(var1$VAR_id, var3$VAR_id), ]
  rownames(var3) <- NULL
  expect_equal(var1, var3)

  # check combination of where clause and VAR_id
  gdb_subset4 <- withr::local_tempfile()
  var <- getAnno(gdb, "var")
  var_chr <- var %>%
    dplyr::filter(CHROM %in% c("chr10", "chr9", "chr15"))
  suppressMessages(subsetGdb(
    gdb,
    where = "CHROM in ('chr10','chr14','chr9')",
    VAR_id = var_chr$VAR_id,
    output = gdb_subset4,
    overWrite = TRUE
  ))
  gdb_subset4 <- gdb(gdb_subset4)
  var4 <- getAnno(gdb_subset4, "var")
  expect_equal(sort(unique(var4$CHROM)), c("chr10", "chr9"))

  # check if keeping specific tables works
  gdb_subset <- withr::local_tempfile()
  suppressMessages(subsetGdb(
    gdb,
    VAR_id = 1:100,
    tables = "pheno",
    output = gdb_subset,
    overWrite = TRUE
  ))
  gdb_subset <- gdb(gdb_subset)
  expect_equal(nrow(listAnno(gdb_subset)), 0L)
  expect_equal(nrow(listCohort(gdb_subset)), 1L)
})


test_that("subsetGdb input validation works correctly", {
  gdb <- create_example_gdb()
  file <- withr::local_tempfile()
  readr::write_lines("hello", file)
  # expect error when writing to same file again when `overWrite = FALSE`
  expect_error(
    {
      subsetGdb(
        gdb,
        output = file,
        where = sprintf(
          "VAR_id in (%s)",
          paste(paste0("'", 1:10, "'"), collapse = ",")
        )
      )
    },
    "already exists"
  )
})

test_that("concatGdb works", {
  # create split gdbs
  gdb <- create_example_gdb()
  vars <- getAnno(gdb, table = "var", fields = "VAR_id")$VAR_id
  vars <- split(
    vars,
    rep(1:3, each = ceiling(length(vars) / 3), length.out = length(vars))
  )
  tmpfiles <- c(
    withr::local_tempfile(),
    withr::local_tempfile(),
    withr::local_tempfile()
  )
  files <- withr::local_tempfile()
  readr::write_lines(tmpfiles, file = files)

  ## split in three parts
  for (i in 1:3) {
    suppressMessages(subsetGdb(
      gdb,
      output = tmpfiles[i],
      where = sprintf(
        "VAR_id in (%s)",
        paste(paste0("'", vars[[i]], "'"), collapse = ",")
      ),
      overWrite = TRUE
    ))
  }

  ## concat three parts
  gdbconcat <- withr::local_tempfile()
  concatGdb(files, output = gdbconcat, verbose = FALSE)
  gdbconcat <- gdb(gdbconcat)

  # check if all tables are present
  expect_equal(
    sort(DBI::dbListTables(gdbconcat)),
    setdiff(sort(gdb_protected_tables), "tmp")
  )

  # check if indexes are present
  indexes <- DBI::dbGetQuery(
    gdbconcat,
    "SELECT * FROM sqlite_master WHERE type='index';"
  )
  expect_equal(
    sort(c("var_idx", "var_idx2", "SM_idx", "dosage_idx")),
    sort(indexes$name)
  )

  # check if indexes can be skipped
  gdbconcat_no_indexes <- withr::local_tempfile()
  concatGdb(
    files,
    output = gdbconcat_no_indexes,
    skipIndexes = TRUE,
    verbose = FALSE
  )
  gdbconcat_no_indexes <- gdb(gdbconcat_no_indexes)
  indexes <- DBI::dbGetQuery(
    gdbconcat_no_indexes,
    "SELECT * FROM sqlite_master WHERE type='index';"
  )
  expect_true(nrow(indexes) == 0L)
  compare_gdbs(gdb, gdbconcat_no_indexes, check_tables = FALSE)
})

test_that("concatGdb input validation works", {
  # expect error when duplicate files are supplied
  targets <- withr::local_tempfile()
  file <- withr::local_tempfile()
  files <- c(file, file)
  readr::write_lines(files, file = targets)
  expect_error(
    concatGdb(
      targets = targets,
      output = withr::local_tempfile()
    ),
    regexp = "Duplicate files"
  )

  # expect error when fewer than 2 files are specified
  gdb <- create_example_gdb()
  targets <- withr::local_tempfile()
  readr::write_lines(getGdbPath(gdb), file = targets)
  expect_error(
    concatGdb(
      targets = targets,
      output = withr::local_tempfile()
    ),
    regexp = "at least 2 valid"
  )

  # expect error when any of the file paths is invalid
  targets <- withr::local_tempfile()
  readr::write_lines(
    c(
      paste(sample(c(letters, 0:9), 100, replace = TRUE), collapse = ""),
      paste(sample(c(letters, 0:9), 100, replace = TRUE), collapse = "")
    ),
    file = targets
  )

  expect_error(
    concatGdb(
      targets = targets,
      output = withr::local_tempfile(),
      verbose = FALSE
    ),
    regexp = "Invalid file path"
  )

  # expect error when target file doesn't exist
  expect_error(
    concatGdb(
      targets = withr::local_tempfile(),
      output = withr::local_tempfile(),
      verbose = FALSE
    ),
    regexp = "does not exist"
  )

  # expect warning when gdbs are built with different RVAT versions
  gdb1 <- create_example_gdb()
  gdb2 <- create_example_gdb()
  meta <- RSQLite::dbGetQuery(gdb1, "select * from meta")
  meta[meta$name == "rvatVersion", ]$value <- paste(
    sample(c(letters, 0:9), 5, replace = TRUE),
    collapse = ""
  )
  RSQLite::dbWriteTable(gdb1, name = "meta", value = meta, overwrite = TRUE)
  targets <- withr::local_tempfile()
  readr::write_lines(c(getGdbPath(gdb1), getGdbPath(gdb2)), file = targets)
  expect_warning(
    concatGdb(
      targets = targets,
      output = withr::local_tempfile(),
      verbose = FALSE
    ),
    regexp = "Not all input gdbs were created"
  )

  # expect error when gdbs have different genome builds
  gdb1 <- create_example_gdb()
  gdb2 <- create_example_gdb()
  meta <- RSQLite::dbGetQuery(gdb1, "select * from meta")
  meta[meta$name == "genomeBuild", ]$value <- "GRCh37"
  RSQLite::dbWriteTable(gdb1, name = "meta", value = meta, overwrite = TRUE)
  targets <- withr::local_tempfile()
  readr::write_lines(c(getGdbPath(gdb1), getGdbPath(gdb2)), file = targets)
  expect_error(
    concatGdb(
      targets = targets,
      output = withr::local_tempfile(),
      verbose = FALSE
    ),
    regexp = "Not all input gdbs share the same genome build"
  )
  

  # expect error when gdbs contain non-matching IIDs
  gdb_subset <- create_example_gdb()
  iids_to_keep <- head(DBI::dbGetQuery(gdb_subset, "select IID from SM")$IID, 50)
  DBI::dbExecute(gdb_subset, 
    sprintf("DELETE FROM SM WHERE IID NOT IN ('%s')", 
            paste(iids_to_keep, collapse = "','")))
  targets <- withr::local_tempfile()
  readr::write_lines(c(getGdbPath(gdb), getGdbPath(gdb_subset)), file = targets)
  
  expect_error(
    concatGdb(
      targets = targets,
      output = withr::local_tempfile(),
      verbose = FALSE
    ),
    regexp = "Sample mismatch between gdbs"
  )
})
