#' getAnno
#' @rdname getAnno
#' @usage NULL
#' @export
setMethod(
  "getAnno",
  signature = "gdb",
  definition = function(
    object,
    table,
    fields = "*",
    left = NULL,
    inner = NULL,
    VAR_id = NULL,
    ranges = NULL,
    padding = 250L,
    where = NULL
  ) {
    # input validation
    check_character(table, allow_na = FALSE)
    check_character(fields, allow_na = FALSE)
    check_character(left, allow_na = FALSE, allow_null = TRUE)
    check_character(inner, allow_na = FALSE, allow_null = TRUE)
    check_character(where, allow_na = FALSE, allow_null = TRUE)
    check_length(table, equal = 1L, allow_na = FALSE)
    check_length(where, equal = 1L, allow_null = TRUE)
    ## padding/ranges is checked in extractRanges if specified

    # base query
    fields <- paste(fields, collapse = ",")
    query <- sprintf("select %s from %s", fields, table)

    # extract ranges, if specified
    if (!is.null(ranges)) {
      VAR_id <- extractRanges(object, ranges = ranges, padding = padding)
      if (length(VAR_id) == 0L) {
        message("No variants overlap with provided range")
      }
    }

    # add left join operation
    if (!is.null(left)) {
      for (i in left) {
        query <- sprintf("%s left join %s using (VAR_id)", query, i)
      }
    }

    # add inner join operations
    if (!is.null(inner)) {
      for (i in inner) {
        query <- sprintf("%s inner join %s using (VAR_id)", query, i)
      }
    }

    # build where statement
    if (!is.null(VAR_id)) {
      VAR_id <- sprintf(
        "VAR_id in (%s)",
        paste(as.integer(VAR_id), collapse = ",")
      )
      if (length(where) > 0L) {
        where <- sprintf("(%s) AND (%s)", VAR_id, where)
      } else {
        where <- VAR_id
      }
    }

    if (!is.null(where)) {
      query <- sprintf("%s where %s", query, where)
    }

    return(DBI::dbGetQuery(object, query))
  }
)


#' getCohort
#' @rdname getCohort
#'
#' @export
setMethod(
  "getCohort",
  signature = "gdb",
  definition = function(object, cohort, fields = "*", keepAll = FALSE) {
    # input validation
    check_character(cohort, allow_na = FALSE)
    check_length(cohort, equal = 1L, allow_na = FALSE)
    check_character(fields, allow_na = FALSE)
    check_logical(keepAll, allow_na = FALSE)
    check_length(keepAll, equal = 1L, allow_na = FALSE)

    # construct query
    fields <- paste(fields, collapse = ",")
    
    if (keepAll) {
      query <- sprintf("select %s from %s", fields, cohort)
    } else {
      query <- sprintf("select %s from %s where IID is not NULL", fields, cohort)
    }

    # extract cohort
    cohort <- DBI::dbGetQuery(object, query)

    # return cohort
    cohort
  }
)


#' @rdname uploadAnno
#' @usage NULL
#' @export
setMethod(
  "uploadAnno",
  signature = "gdb",
  definition = function(
    object,
    name,
    value,
    sep,
    skipRemap,
    skipIndexes,
    ignoreAlleles,
    keepUnmapped,
    mapRef,
    overWrite,
    verbose
  ) {
    # validate input
    .uploadAnno_validate_input(as.list(environment()))

    # validate table name and input
    .upload_validate_input(
      object = object,
      name = name,
      value = value,
      overWrite = overWrite,
      verbose = verbose
    )

    # load data and return source info
    source_info <- .uploadAnno_load_data(
      gdb = object,
      name = name,
      value = value,
      sep = sep,
      verbose = verbose
    )

    # map positions to VAR_id (if skipRemap=FALSE)
    if (!skipRemap) {
      .uploadAnno_remap_to_varid(
        gdb = object,
        name = name,
        ignoreAlleles = ignoreAlleles,
        keepUnmapped = keepUnmapped,
        mapRef = mapRef,
        verbose = verbose
      )
    }

    # create indices
    if (!skipIndexes) {
      .uploadAnno_create_indices(
        gdb = object,
        name = name,
        verbose = verbose
      )
    }

    # Update annotation meta-data table
    .gdb_upload_update_metadata(
      gdb = object,
      table_name = name,
      type = "anno",
      source_info = source_info,
      verbose = verbose
    )

    # return nothing
    invisible(NULL)
  }
)

.uploadAnno_validate_input <- function(args) {
  check_wrapper(check_character, args, "name", length_equal = 1L)
  check_wrapper(check_character, args, "sep", length_equal = 1L)
  check_wrapper(check_bool, args, "skipRemap", length_equal = 1L)
  check_wrapper(check_bool, args, "skipIndexes", length_equal = 1L)
  check_wrapper(check_bool, args, "ignoreAlleles", length_equal = 1L)
  check_wrapper(check_bool, args, "keepUnmapped", length_equal = 1L)
  check_wrapper(check_character, args, "mapRef", length_equal = 1L)
  check_wrapper(check_bool, args, "overWrite", length_equal = 1L)
  check_wrapper(check_bool, args, "verbose", length_equal = 1L)

  invisible(NULL)
}

.uploadAnno_load_data <- function(gdb, name, value, sep, verbose = TRUE) {
  if (is.data.frame(value)) {
    # if value is a data.frame upload directly
    if (verbose) {
      message(sprintf(
        "Loading table '%s' from interactive R session'",
        name
      ))
    }
    DBI::dbWriteTable(
      con = gdb,
      name = name,
      value = value,
      overwrite = TRUE
    )
    source_info <- "interactive_session"
  } else if (is.character(value)) {
    # if value is a filepath import from file
    if (verbose) {
      message(sprintf("Loading table '%s' from '%s'", name, value))
    }
    DBI::dbWriteTable(
      con = gdb,
      name = name,
      value = value,
      sep = sep,
      overwrite = TRUE
    )
    source_info <- value
  } else {
    stop(
      "`value` must be a data.frame or ",
      "a character string containing a file path.",
      call. = FALSE
    )
  }

  if (verbose) {
    fields <- DBI::dbListFields(gdb, name)
    message(sprintf(
      "%s fields detected (%s)",
      length(fields),
      paste(fields, collapse = ",")
    ))
  }

  # return source info
  source_info
}

.uploadAnno_remap_to_varid <- function(
  gdb,
  name,
  ignoreAlleles,
  keepUnmapped,
  mapRef,
  verbose
) {
  try(DBI::dbExecute(gdb, "DROP TABLE IF EXISTS tmp"), silent = TRUE)

  # determine join keys
  if (ignoreAlleles) {
    join_keys <- c("CHROM", "POS")
  } else {
    join_keys <- c("CHROM", "POS", "REF", "ALT")
  }
  join_keys_sql <- paste(join_keys, collapse = ", ")

  # check fields
  fields <- DBI::dbListFields(gdb, name)
  if (length(fields) == 0L) {
    stop(
      sprintf(
        "Table '%s' is empty or does not exist for remapping.",
        name
      ),
      call. = FALSE
    )
  }

  if (identical(fields, "VAR_id")) {
    stop(
      "Annotation table contains only a `VAR_id` field. ",
      "To upload this table, set `skipRemap = TRUE`",
      call. = FALSE
    )
  }

  if (!all(join_keys %in% fields)) {
    stop(
      sprintf(
        "Mapping not possible: %s not provided.",
        paste(join_keys[!join_keys %in% fields], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # create indices, drop if index already exists
  index_name <- paste0(name, "_idx2")

  try(
    DBI::dbExecute(gdb, sprintf("DROP INDEX IF EXISTS %s", index_name)),
    silent = TRUE
  )

  create_index_sql <- sprintf(
    "CREATE INDEX %s ON %s (%s)",
    index_name,
    name,
    join_keys_sql
  )
  DBI::dbExecute(gdb, create_index_sql)

  # join (inner or left depending on `keepUnmapped`)
  join_type <- if (keepUnmapped) "LEFT JOIN" else "INNER JOIN"
  selects <- paste(
    c(
      sprintf("%s.VAR_id", mapRef),
      sprintf(
        "%s.%s",
        name,
        fields[fields != "VAR_id"]
      )
    ),
    collapse = ","
  )

  DBI::dbExecute(
    gdb,
    sprintf(
      "CREATE TABLE tmp AS SELECT %s FROM %s %s %s USING (%s)",
      selects,
      name,
      join_type,
      mapRef,
      join_keys_sql
    )
  )

  # check number of unmapped rows
  original_count <- DBI::dbGetQuery(
    gdb,
    sprintf("SELECT COUNT(*) AS n FROM %s", name)
  )$n[1]

  if (keepUnmapped) {
    unmapped_count <- DBI::dbGetQuery(
      gdb,
      "SELECT COUNT(*) AS n FROM tmp WHERE VAR_id IS NULL"
    )$n[1]
  } else {
    remapped_count <- DBI::dbGetQuery(
      gdb,
      "SELECT COUNT(*) AS n FROM tmp"
    )$n[1]
    unmapped_count <- original_count - remapped_count
  }

  if (unmapped_count > 0L) {
    warning(
      sprintf(
        "%s rows could not be mapped to variants in the gdb",
        unmapped_count
      ),
      call. = FALSE
    )
  }

  # drop original uploaded table and rename tmp to table
  DBI::dbExecute(gdb, sprintf("DROP TABLE %s", name))
  DBI::dbExecute(gdb, sprintf("ALTER TABLE tmp RENAME TO %s", name))

  # return nothing
  invisible(NULL)
}

.uploadAnno_create_indices <- function(
  gdb,
  name,
  verbose
) {
  fields <- DBI::dbListFields(gdb, name)
  if ("VAR_id" %in% fields) {
    DBI::dbExecute(
      gdb,
      sprintf("CREATE INDEX %s_idx1 ON %s (VAR_id)", name, name)
    )
  }
  if (all(c("CHROM", "POS", "REF", "ALT") %in% fields)) {
    DBI::dbExecute(
      gdb,
      sprintf(
        "CREATE INDEX %s_idx2 ON %s (CHROM, POS, REF, ALT)",
        name,
        name
      )
    )
  } else if (all(c("CHROM", "POS") %in% fields)) {
    DBI::dbExecute(
      gdb,
      sprintf("CREATE INDEX %s_idx2 ON %s (CHROM, POS)", name, name)
    )
  }

  invisible(NULL)
}


#' @rdname uploadCohort
#' @usage NULL
#' @export
setMethod(
  "uploadCohort",
  signature = "gdb",
  definition = function(
    object,
    name,
    value,
    sep = "\t",
    overWrite = FALSE,
    verbose = TRUE
  ) {
    # input validation
    .uploadCohort_validate_input(as.list(environment()))

    # validate input
    .upload_validate_input(
      object = object,
      name = name,
      value = value,
      overWrite = overWrite,
      verbose = verbose
    )

    # load cohort
    loaded <- .uploadCohort_load_data(
      value = value,
      sep = sep,
      name = name,
      verbose = verbose
    )
    cohort <- loaded$cohort
    cohort_source <- loaded$cohort_source
    rm(loaded)

    # check and format IID,sex columns
    cohort <- .uploadCohort_validate_cohort(
      cohort = cohort,
      verbose = verbose
    )

    # match cohort with gdb IIDs
    cohort <- .align_cohort_SM(
      gdb = object,
      cohort = cohort
    )

    # upload table
    DBI::dbWriteTable(
      con = object,
      name = name,
      value = cohort,
      overwrite = TRUE
    )

    # update cohort metadata table
    .gdb_upload_update_metadata(
      gdb = object,
      table_name = name,
      type = "cohort",
      source_info = cohort_source,
      verbose = verbose
    )

    invisible(NULL)
  }
)

.uploadCohort_validate_input <- function(args) {
  check_wrapper(check_character, args, "name", length_equal = 1L)
  check_wrapper(check_character, args, "sep", length_equal = 1L)
  check_wrapper(check_bool, args, "overWrite", length_equal = 1L)
  check_wrapper(check_bool, args, "verbose", length_equal = 1L)

  invisible(NULL)
}

.uploadCohort_validate_cohort <- function(cohort, verbose) {
  # message detected fields
  if (verbose) {
    message(sprintf(
      "%s fields detected (%s)",
      ncol(cohort),
      paste(colnames(cohort), collapse = ",")
    ))
  }

  # check if cohort contains IID field
  if (!("IID" %in% colnames(cohort))) {
    stop("No 'IID' field detected, aborting upload.", call. = FALSE)
  }

  # check if cohort contains sex field
  if (!("sex" %in% colnames(cohort))) {
    stop("No 'sex' field detected, aborting upload.", call. = FALSE)
  }

  # check if cohort contains duplicated IIDs
  if (anyDuplicated(cohort$IID) != 0L) {
    stop("The cohort contains duplicated IIDs.", call. = FALSE)
  }

  # check miscoded sex values
  if (!is.numeric(cohort$sex)) {
    original <- cohort$sex
    cohort$sex <- suppressWarnings(as.numeric(as.character(cohort$sex)))
    if (any(is.na(cohort$sex) & !is.na(original))) {
      warning(
        "`sex` should be coded as 0=unknown, 1=male, 2=female.",
        "Non-numeric values in 'sex' column were coerced to NA. ",
        call. = FALSE
      )
    }
  }

  cohort$sex[is.na(cohort$sex)] <- 0
  if (!all(cohort$sex %in% c(0, 1, 2))) {
    warning(
      sprintf(
        "%s values in `sex` column not coded as 0,1,2. ",
        sum(!cohort$sex %in% c(0, 1, 2))
      ),
      "These will be set to 0 (unknown).",
      call. = FALSE
    )
    cohort$sex[!(cohort$sex %in% c(0, 1, 2))] <- 0
  }

  if (verbose) {
    message(sprintf(
      "%s males, %s females and %s unknown sex",
      sum(cohort$sex == 1, na.rm = TRUE),
      sum(cohort$sex == 2, na.rm = TRUE),
      sum(cohort$sex == 0, na.rm = TRUE)
    ))
  }

  # return formatted cohort
  cohort
}


.uploadCohort_load_data <- function(value, sep, name, verbose) {
  if (is.character(value)) {
    # if filepath is provided -> read data
    if (verbose) {
      message(sprintf("Loading cohort '%s' from '%s'\n", name, value))
    }
    cohort <- read.table(
      value,
      sep = sep,
      as.is = TRUE,
      header = TRUE,
      stringsAsFactors = FALSE
    )
    cohort_source <- value
  } else {
    # value is a data.frame
    if (verbose) {
      message(sprintf(
        "Loading cohort '%s' from interactive R session\n",
        name
      ))
    }
    cohort <- value
    cohort_source <- "interactive_session"
  }

  list(cohort = cohort, cohort_source = cohort_source)
}


#' @rdname gdb
#' @usage NULL
#' @export
setMethod(
  "dropTable",
  signature = "gdb",
  definition = function(object, name, verbose = TRUE) {
    # input validation
    check_character(name, allow_na = FALSE)
    check_length(name, equal = 1L, allow_na = FALSE)
    check_bool(verbose, allow_na = FALSE)
    check_length(verbose, equal = 1L, allow_na = FALSE)

    # check if table exists -> both as a table and in metadata
    meta_anno_exists <- tolower(name) %in% tolower(listAnno(object)$name)
    meta_cohort_exists <- tolower(name) %in% tolower(listCohort(object)$name)
    table_exists <- tolower(name) %in% tolower(DBI::dbListTables(object))

    # remove from anno meta (if exists)
    if (meta_anno_exists) {
      DBI::dbExecute(
        object,
        "DELETE FROM anno WHERE name = :table_to_remove",
        params = list(table_to_remove = name)
      )
    }

    # remove from cohort meta (if exists)
    if (meta_cohort_exists) {
      DBI::dbExecute(
        object,
        "DELETE FROM cohort WHERE name = :table_to_remove",
        params = list(table_to_remove = name)
      )
    }

    # drop table
    if (table_exists) {
      drop_sql <- sprintf("DROP TABLE IF EXISTS %s", name)
      DBI::dbExecute(object, drop_sql)
      if (verbose) message(sprintf("Table '%s' removed from gdb", name))
    } else {
      warning(
        sprintf(
          "Table '%s' doesn't exist in the gdb, and therefore could not be dropped.",
          name
        ),
        call. = FALSE
      )
    }

    # return nothing
    invisible(NULL)
  }
)
