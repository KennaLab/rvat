# check VAR_id
.gdb_check_varid <- function(VAR_id) {
  if (!is.null(VAR_id)) {
    if (
      !(is.vector(VAR_id) &&
        (is.numeric(VAR_id) || is.character(VAR_id)))
    ) {
      stop("`VAR_id` should be a numeric or character vector.", call. = FALSE)
    }
    if (anyNA(VAR_id)) {
      stop("`VAR_id` contains missing values.", call. = FALSE)
    }
    if (anyDuplicated(VAR_id) != 0L) {
      stop("`VAR_id` contains duplicated values.", call. = FALSE)
    }
  }

  invisible(NULL)
}

.gdb_check_SM <- function(SM, GT, expect_equal_dims = TRUE) {
  if (!is.data.frame(SM) && !is(SM, "DFrame")) {
    stop("`SM` must be a data.frame or DFrame.", call. = FALSE)
  }
  if (!"IID" %in% colnames(SM)) {
    stop("`SM` (sample manifest) must contain an 'IID' column.", call. = FALSE)
  }
  if (anyDuplicated(SM$IID[!is.na(SM$IID)]) != 0L) {
    stop("Duplicated IID values found in cohort.", call. = FALSE)
  }
  if (expect_equal_dims && nrow(SM) != ncol(GT)) {
    stop(
      sprintf(
        "Number of samples in `SM` (%d) doesn't match those in the genoMatrix (%d)",
        nrow(SM),
        ncol(GT)
      ),
      call. = FALSE
    )
  }
  if (!setequal(SM$IID, colnames(GT))) {
    stop(
      "IID values in SM table do not match existing genoMatrix column names.",
      call. = FALSE
    )
  }
}

.gdb_list_metadata <- function(gdb) {
  metadata <- list(
    rvatVersion = as.character(packageVersion("rvat")),
    gdbId = getGdbId(gdb),
    genomeBuild = getGenomeBuild(gdb),
    creationDate = as.character(round(Sys.time(), units = "secs"))
  )
  metadata
}

# populate metadata gdb
.gdb_populate_meta_table <- function(gdb, genomeBuild, verbose = TRUE) {
  # add rvat version to meta table
  DBI::dbExecute(
    gdb,
    "INSERT INTO meta VALUES (?, ?)",
    params = list(
      "rvatVersion",
      as.character(packageVersion("rvat"))
    )
  )

  # add random identifier
  DBI::dbExecute(
    gdb,
    "INSERT INTO meta VALUES (?, ?)",
    params = list(
      "id",
      paste(sample(c(letters, 0:9), 28L, replace = TRUE), collapse = "")
    )
  )

  # add genome build
  if (!is.null(genomeBuild)) {
    ## added this to insert into genomeBuild into the meta table
    DBI::dbExecute(
      gdb,
      "INSERT INTO meta VALUES (?, ?)",
      params = list("genomeBuild", genomeBuild)
    )
  }
  if (!genomeBuild %in% names(nonPAR)) {
    warning(
      sprintf(
        paste0(
          "The supplied genomeBuild is not supported by RVAT. ",
          "The build will be included in the gdb metadata, ",
          "but won't be used in downstream RVAT analyses, such as ",
          "correctly assigning ploidies in pseudoautosomal regions. ",
          "Supported builds include: %s."
        ),
        paste(names(nonPAR), collapse = ",")
      ),
      call. = FALSE
    )
  }
}


.gdb_create_indexes <- function(
  gdb,
  verbose
) {
  # var indexes
  if (verbose) {
    message(sprintf(
      "%s\tCreating var table indexes",
      as.character(round(Sys.time(), units = "secs"))
    ))
  }
  DBI::dbExecute(gdb, "create index IF NOT EXISTS var_idx on var (VAR_id)")
  DBI::dbExecute(
    gdb,
    "create index IF NOT EXISTS var_idx2 on var (CHROM,POS,REF,ALT)"
  )

  # SM indexes
  if (verbose) {
    message(sprintf(
      "%s\tCreating SM table index",
      as.character(round(Sys.time(), units = "secs"))
    ))
  }
  DBI::dbExecute(gdb, "create index IF NOT EXISTS SM_idx on SM (IID)")

  # dosage indexes
  if (verbose) {
    message(sprintf(
      "%s\tCreating dosage table index",
      as.character(round(Sys.time(), units = "secs"))
    ))
  }
  DBI::dbExecute(
    gdb,
    "create index IF NOT EXISTS dosage_idx on dosage (VAR_id)"
  )

  DBI::dbExecute(gdb, "ANALYZE") ##TODO: SEE IF THIS IS INTERESTING TO USE

  # return nothing
  invisible(NULL)
}

.gdb_check_tables <- function(
  gdb,
  skip = NULL,
  warning = TRUE,
  return_problems = TRUE
) {
  # expected tables based on metadata
  tables_expected <- c(
    gdb_protected_tables,
    listAnno(gdb)$name,
    listCohort(gdb)$name
  )
  if (!is.null(skip)) {
    tables_expected <- tables_expected[!tables_expected %in% skip]
  }

  # unexpected tables
  tables_found <- DBI::dbListTables(gdb)
  tables_not_expected <- setdiff(tables_found, tables_expected)
  tables_expected_not_found <- setdiff(tables_expected, tables_found)

  if (length(tables_not_expected) >= 1L && warning) {
    warning(
      sprintf(
        paste0(
          "The following table(s): '%s', are present in the gdb but are not ",
          "tracked in the gdb metadata."
        ),
        paste(tables_not_expected, collapse = ",")
      ),
      call. = FALSE
    )
  }

  if (length(tables_expected_not_found) >= 1L && warning) {
    warning(
      sprintf(
        paste0(
          "The following table(s): '%s', are expected based on the gdb metadata",
          " but are not present in the gdb."
        ),
        paste(tables_expected_not_found, collapse = ",")
      ),
      call. = FALSE
    )
  }

  if (return_problems) {
    list(
      tables_not_expected = tables_not_expected,
      tables_expected_not_found = tables_expected_not_found
    )
  } else {
    invisible(NULL)
  }
}

# validate anno/cohort name
# `uploadCohort`, `uploadAnno` and `mapVariants`
.upload_validate_name <- function(
  object,
  name,
  overWrite,
  verbose
) {
  # check if table name is protected
  if (tolower(name) %in% tolower(gdb_protected_tables)) {
    stop(
      sprintf(
        "'%s' already exists as a protected gdb table and cannot be replaced",
        name
      ),
      call. = FALSE
    )
  }

  # check if table names contains punctuation or spaces
  if (grepl("[^A-Za-z0-9_]", name)) {
    stop(
      "Table name may only contain letters, numbers and underscores.",
      call. = FALSE
    )
  }

  # check if table name already exists
  existing_tables <- tolower(c(
    listAnno(object)$name,
    listCohort(object)$name
  ))
  if (tolower(name) %in% existing_tables && !overWrite) {
    stop(
      sprintf(
        paste0(
          "Table name '%s' is already in use. ",
          "Choose a different name, set `overWrite=TRUE`, or delete the table using `dropTable`."
        ),
        name
      ),
      call. = FALSE
    )
  } else if (tolower(name) %in% existing_tables && overWrite) {
    if (verbose) {
      message(
        sprintf(
          "`%s` already exists and will be overwritten given `overWrite = TRUE`.",
          name
        )
      )
    }
    dropTable(object = object, name = name, verbose = verbose)
  }

  invisible(NULL)
}


# validation used for both `uploadCohort` and `uploadAnno`
.upload_validate_input <- function(
  object,
  name,
  value,
  overWrite,
  verbose
) {
  # validate anno/cohort name
  .upload_validate_name(
    object = object,
    name = name,
    overWrite = overWrite,
    verbose = verbose
  )

  # value must be either a data.frame of a string containing a file path
  if (!is.data.frame(value) && !is.character(value)) {
    stop(
      "`value` must be a data.frame or ",
      "a character string containing a file path.",
      call. = FALSE
    )
  }

  # if value is a character:
  # - it should be of length 1 (a single filepath)
  # - it should exist
  if (is.character(value)) {
    if (length(value) != 1L) {
      stop(
        "`value` must be a single file path, not a vector of strings.",
        call. = FALSE
      )
    }
    if (!file.exists(value)) {
      stop(
        sprintf("The specified file does not exist: '%s'", value),
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}

# update metadata on upload
.gdb_upload_update_metadata <- function(
  gdb,
  table_name,
  type,
  source_info,
  verbose = TRUE
) {
  if (type == "anno") {
    meta_table <- "anno"
    existing_tables <- listAnno(gdb)$name
  } else if (type == "cohort") {
    meta_table <- "cohort"
    existing_tables <- listCohort(gdb)$name
  } else {
    stop(
      "Invalid 'type' specified. ",
      "Must be 'anno' or 'cohort'.",
      call. = FALSE
    )
  }

  if (table_name %in% existing_tables) {
    delete_query <- sprintf(
      "DELETE FROM %s WHERE name = ?",
      meta_table
    )
    DBI::dbExecute(
      gdb,
      delete_query,
      params = list(table_name)
    )
  }

  # using :name, :value, :date as placeholders
  insert_query <- sprintf(
    "INSERT INTO %s (name, value, date) VALUES (?, ?, ?)",
    meta_table
  )
  DBI::dbExecute(
    gdb,
    insert_query,
    params = list(
      table_name,
      source_info,
      date()
    )
  )

  invisible(NULL)
}

# align uploaded cohort to SM
.align_cohort_SM <- function(gdb, cohort) {
  IID <- as.character(
    DBI::dbGetQuery(gdb, "select * from SM")$IID
  )

  if (!all(cohort$IID %in% IID)) {
    warning(
      sprintf(
        "%s/%s samples in the cohort are not present in the gdb.",
        sum(!unique(cohort$IID) %in% IID),
        length(unique(cohort$IID))
      ),
      call. = FALSE
    )
  }

  # match cohort with SM IIDs
  SM <- cohort[match(IID, cohort$IID), , drop = FALSE]

  SM
}
