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

compare_varsetfile <- function(vsfile1, vsfile2, units = NULL) {
  if (is.null(units)) {
    units <- listUnits(vsfile1)
  }
  vs1 <- getVarSet(vsfile1, unit = units)
  vs2 <- getVarSet(vsfile2, unit = units)
  vs1@metadata <- list()
  vs2@metadata <- list()
  expect_identical(vs1, vs2)
}

compare_varsets <- function(varset1, varset2) {
  varset1@metadata <- list()
  varset2@metadata <- list()
  expect_identical(varset1, varset2)
}

create_merged_varset_file <- function() {
  gdb <- gdb(rvat_example("rvatData.gdb"))
  
  # Build varsets
  suppressMessages(buildVarSet(
    object = gdb,
    output = setup$moderate_file,
    varSetName = "Moderate",
    unitTable = "varInfo",
    unitName = "gene_name", 
    where = paste("(ModerateImpact = 1 or HighImpact = 1) and", setup$where_clause)
  ))
  
  suppressMessages(buildVarSet(
    object = gdb,
    output = setup$lof_file,
    varSetName = "HighImpact",
    unitTable = "varInfo",
    unitName = "gene_name",
    where = paste("HighImpact = 1 and", setup$where_clause),
    verbose = FALSE
  ))
  
  suppressMessages(buildVarSet(
    object = gdb,
    output = setup$cadd_file,
    varSetName = "CADD",
    unitTable = "varInfo", 
    unitName = "gene_name",
    weightName = "CADDphred",
    where = setup$where_clause
  ))
  
  # Combine varsets
  varsets <- list()
  for (i in seq_along(c(setup$moderate_file, setup$lof_file, setup$cadd_file))) {
    file <- c(setup$moderate_file, setup$lof_file, setup$cadd_file)[i]
    varset <- readr::read_delim(file, delim = "|", col_names = FALSE, col_types = "cccc", comment = "#")
    varsets[[i]] <- varset
  }
  
  varsets <- dplyr::bind_rows(varsets) |> dplyr::arrange(X1)
  metadata_obj <- metadata(varSetFile(setup$moderate_file))
  
  merged_file <- withr::local_tempfile()
  con <- gzfile(merged_file, "w") 
  rvat:::.write_rvat_header(filetype = "varSetFile", metadata = metadata_obj, con = con)
  write.table(varsets, file = con, col.names = FALSE, sep = "|", append = TRUE, quote = FALSE, row.names = FALSE)
  close(con)
  
  merged_file
}