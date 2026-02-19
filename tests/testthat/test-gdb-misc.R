# gdb validation
test_that("gdb integrity checks work", {

  # drop a table from the meta table
  gdb <- create_example_gdb()
  DBI::dbExecute(gdb, "DELETE FROM anno WHERE name = 'varInfo';")
  expect_warning({
    .gdb_check_tables(gdb, return_problems = FALSE, skip = "tmp")
  }, regexp = "not tracked in the gdb metadata")

  # add table that is not tracked in meta
  gdb <- create_example_gdb()
  DBI::dbExecute(gdb, "INSERT INTO anno (name, value, date) VALUES ('hello', 'test', 'today');")
  expect_warning({
    problems <- .gdb_check_tables(gdb, skip = "tmp")
  }, regexp = "are expected based on the gdb metadata but are not present")
})

# misc methods
test_that("misc gdb methods work", {
  gdb <- withr::local_tempfile(fileext = ".gdb")
  buildGdb(
    test_path("data/vcf_multiallelic.vcf"),
    output = gdb,
    genomeBuild = "GRCh38",
    verbose = FALSE
  )
  gdb <- gdb(gdb)
  gdb_no_meta <- create_example_gdb()
  DBI::dbExecute(gdb_no_meta, "DROP TABLE meta;")
  DBI::dbWriteTable(gdb_no_meta, 
                   "meta",
                    value = data.frame(name = character(0), value = character(0)))

  # show gdb connection
  expect_true(stringr::str_detect(
    capture_output({
      show(gdb)
    }),
    "gdb object"
  ))

  # get rvat version
  expect_equal(getRvatVersion(gdb), as.character(packageVersion("rvat")))
  expect_equal(getRvatVersion(gdb_no_meta), NA_character_)


  # get creation date
  expect_true(
    (stringr::str_detect(getCreationDate(gdb), "-") &&
      stringr::str_detect(getCreationDate(gdb), ":"))
  )
  expect_equal(getCreationDate(gdb_no_meta), NA_character_)
  
  # get gdb ID
  expect_equal(getGdbId(gdb_no_meta), NA_character_)
  
  # get genome build
  expect_equal(getGenomeBuild(gdb_no_meta), NA_character_)

  # expect error when connecting to nonexisting file
  expect_error(
    {
      gdb(paste(sample(LETTERS, 10), collapse = ""))
    },
    "doesn't exist"
  )
})

test_that("dropTable works", {
  gdb <- create_example_gdb()

  # add table and check if removed upon running `dropTable`
  all_tables <- RSQLite::dbListTables(gdb)
  uploadAnno(
    gdb,
    name = "test_table",
    value = data.frame(
      VAR_id = c(1, 1, 2, 2, 2, 3),
      transcript = c("A", "B", "C", "D", "E", "F")
    ),
    skipRemap = TRUE,
    verbose = FALSE
  )
  all_tables2 <- RSQLite::dbListTables(gdb)
  expect_equal(setdiff(all_tables2, all_tables), "test_table")

  # drop the table
  suppressMessages(dropTable(gdb, "test_table"))
  all_tables3 <- RSQLite::dbListTables(gdb)
  expect_equal(all_tables, all_tables3)

  # expect warning when trying to drop a table that doesn't exit
  expect_warning(
    {
      dropTable(gdb, "hello")
    },
    regexp = "doesn't exist"
  )
})