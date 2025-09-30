compare_gdbs <- function(gdb1, gdb2, cohort1 = "SM", cohort2 = "SM", check_var = TRUE, check_tables = TRUE) {
  # compare var tables
  if (check_var) {
    anno1 <- getAnno(gdb1, table = "var")
    anno2 <- getAnno(gdb2, table = "var")
    expect_identical(anno1, anno2)
  }
  
  # check if same tables are present
  if (check_tables) {
    tables1 <- DBI::dbListTables(gdb1)
    tables2 <- DBI::dbListTables(gdb2)
    expect_identical(sort(tables1), sort(tables2))
  }

  # compare GTs
  GT1 <- getGT(
    gdb1,
    cohort = cohort1,
    VAR_id = getAnno(gdb1, table = "var", fields = "VAR_id")$VAR_id,
    verbose = FALSE
  )
  GT2 <- getGT(
    gdb2,
    cohort = cohort2,
    VAR_id = getAnno(gdb2, table = "var", fields = "VAR_id")$VAR_id,
    verbose = FALSE
  )
  ## ignore gdb and gdbId in comparison
  metadata(GT1)$gdb <- NA_character_
  metadata(GT2)$gdb <- NA_character_
  metadata(GT1)$gdbId <- NA_character_
  metadata(GT2)$gdbId <- NA_character_

  expect_identical(GT1, GT2)
  
  invisible(NULL)
}


expect_gdb_indexes <- function(gdb, expected_indexes = c("var_idx", "var_idx2", "SM_idx", "dosage_idx")) {
  indexes <- DBI::dbGetQuery(gdb, "SELECT name FROM sqlite_master WHERE type='index'")
  expect_equal(sort(expected_indexes), sort(indexes$name))
}
