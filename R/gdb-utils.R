#' Concatenate gdb databases
#'
#' Function to concatenate [`gdb`] databases. Only retains content of base tables (SM, var, dosage).
#'
#' @param targets File listing full paths of gdbs to concatenate.
#' @param output Output gdb file path.
#' @param skipRemap Skip resetting of VAR_id to row id after concatenation? Defaults to `FALSE`.
#' @param skipIndexes Skip generation of standard var and dosage table indexes (VAR_id;CHROM, POS,REF,ALT)?
#' Defaults to `FALSE`.
#' @param overWrite Overwrite output file if it already exists?
#' Defaults to `FALSE`.
#' @param verbose Should the function be verbose? Defaults to `TRUE`.
#' @example inst/examples/example-concatGdb.R
#'
#' @export
concatGdb <- function(
  targets,
  output,
  skipRemap = FALSE,
  skipIndexes = FALSE,
  overWrite = FALSE,
  verbose = TRUE
) {
  # validate input
  arg <- as.list(environment())
  .concatgdb_validate_input(arg)

  # check if output exists and overwrite if `overWrite = TRUE`
  .check_output(
    output = output,
    overWrite = overWrite,
    verbose = verbose
  )

  # read gdb filepaths
  gdb_list <- scan(targets, what = "character", quiet = TRUE)
  
  # validate gdb filepaths (existence, duplicates, at least 2 gdbs)
  .concatgdb_validate_gdbs(gdb_list)

  # initialize gdb and close on exit
  gdb <- gdb_init(output)
  on.exit(close(gdb), add = TRUE)

  # create tables
  if (verbose) {
    message(sprintf(
      "%s\tCreating gdb tables",
      as.character(round(Sys.time(), units = "secs"))
    ))
  }
  .concatgdb_create_schema(
    gdb = gdb,
    gdb_list = gdb_list,
    verbose = verbose
  )

  # copy var and dosage data from source db
  versions <- vector(mode = "character", length = length(gdb_list))
  builds <- vector(mode = "character", length = length(gdb_list))
  names(versions) <- gdb_list
  names(builds) <- gdb_list
  for (gdb_i in gdb_list) {
    if (verbose) message(sprintf("Merging '%s'", gdb_i))
    DBI::dbExecute(gdb, "attach :src as src", params = list(src = gdb_i))
    DBI::dbExecute(gdb, "insert into var select * from src.var")
    DBI::dbExecute(gdb, "insert into dosage select * from src.dosage")
    versions[gdb_i] <- dbGetQuery(
      gdb,
      "select * from src.meta where name = 'rvatVersion'"
    )$value
    builds[gdb_i] <- dbGetQuery(
      gdb,
      "select * from src.meta where name = 'genomeBuild'"
    )$value
    DBI::dbExecute(gdb, "detach src")
  }
  
  # check if all gdbs share the same version and genome build
  version <- unique(versions)
  build <- unique(builds)
  if (length(version) > 1L) {
    warning(
      "Not all input gdbs were created with the same RVAT version, ",
      "we recommend generating the input gdbs with the same RVAT version.",
      "The version of the first input gdb will be recorded in the output gdb.",
      call. = FALSE
    )
    version <- version[1]
  }

  if (length(build) > 1L) {
    stop("Not all input gdbs share the same genome build.", call. = FALSE)
  }

  # reset VAR_id to row id
  if (!skipRemap) {
    if (verbose) {
      message(sprintf(
        "%s\tReseting VAR_id to rowid",
        as.character(round(Sys.time(), units = "secs"))
      ))
    }
    DBI::dbExecute(gdb, "update var set VAR_id=rowid")
    DBI::dbExecute(gdb, "update dosage set VAR_id=rowid")
  }

  # generate indexes
  if (!skipIndexes) {
    .gdb_create_indexes(
      gdb = gdb,
      verbose = verbose
    )
  }

  # generate var_ranges table
  if (verbose) {
    message(sprintf(
      "%s\tCreating ranged var table",
      as.character(round(Sys.time(), units = "secs"))
    ))
  }
  addRangedVarinfo(gdb, overwrite = TRUE, verbose = verbose)

  # populate meta table
  .gdb_populate_meta_table(
    gdb = gdb,
    genomeBuild = build,
    verbose = verbose
  )

  # finished
  if (verbose) {
    message(sprintf(
      "%s\tComplete",
      as.character(round(Sys.time(), units = "secs"))
    ))
  }

  invisible(NULL)
}

.concatgdb_validate_input <- function(
  args
) {
  # targets should be a single existing filepath
  check_wrapper(check_character, args, "targets", length_equal = 1L)
  if (!file.exists(args[["targets"]])) {
    stop(
      sprintf(
        "File specified in `targets`: '%s', does not exist.",
        args[["targets"]]
      ),
      call. = FALSE
    )
  }

  # output should be a single filepath
  check_wrapper(check_character, args, "output", length_equal = 1L)

  # boolean flags
  check_wrapper(check_bool, args, "skipRemap", length_equal = 1L)
  check_wrapper(check_bool, args, "skipIndexes", length_equal = 1L)
  check_wrapper(check_bool, args, "overWrite", length_equal = 1L)
  check_wrapper(check_bool, args, "verbose", length_equal = 1L)

  invisible(NULL)
}

.concatgdb_validate_gdbs <- function(gdb_list) {
  # check if more than two files are included
  if (length(gdb_list) < 2L) {
    stop("Require at least 2 valid gdb files for merging", call. = FALSE)
  }

  # check if no duplicate files are included
  if (anyDuplicated(gdb_list) != 0L) {
    stop("Duplicate files are included in the targets file", call. = FALSE)
  }

  # check if filepaths are valid
  IID <- NULL
  for (i in gdb_list) {
    if (!file.exists(i)) {
      stop(sprintf("Invalid file path '%s'", i), call. = FALSE)
    }
    gdb_i <- gdb(i)
    if (i == gdb_list[1]) {
      IID <- DBI::dbGetQuery(gdb_i, "select * from SM")$IID 
    } else {
      IID_i <- DBI::dbGetQuery(gdb_i, "select * from SM")$IID
      if (!identical(IID, IID_i)) {
        stop(
          sprintf(
            "Sample mismatch between gdbs, samples in '%s' do not match samples in '%s'",
            i,
            gdb_list[1]
          ),
          call. = FALSE
        )
      }
    }
  }

  invisible(NULL)
}


.concatgdb_create_schema <- function(gdb, gdb_list, verbose) {
  if (verbose)
    message(sprintf(
      "%s\tCreating db tables",
      as.character(round(Sys.time(), units = "secs"))
    ))

  # create new var and dosage tables
  DBI::dbExecute(
    gdb,
    "create table var (VAR_id integer, CHROM text, POS int, ID text, REF text, ALT text, QUAL text, FILTER text, INFO text, FORMAT text);"
  )
  DBI::dbExecute(gdb, "create table dosage (VAR_id integer, GT BLOB);")

  # copy SM table and create new anno, cohort and meta metadata tables
  if (verbose) message("Creating SM, anno and cohort tables")
  DBI::dbExecute(gdb, "attach :src as src", params = list(src = gdb_list[1]))
  DBI::dbExecute(gdb, "create table SM as select * from src.SM")
  DBI::dbExecute(gdb, "detach src")
  DBI::dbExecute(gdb, "create table anno (name text,value text,date text)")
  DBI::dbExecute(gdb, "create table cohort (name text,value text,date text)")
  DBI::dbExecute(gdb, "create table meta (name text,value text)")

  invisible(NULL)
}

#' @rdname subsetGdb
#' @usage NULL
#' @export
setMethod(
  "subsetGdb",
  signature = "gdb",
  definition = function(
    object,
    output,
    intersection,
    where,
    VAR_id,
    tables,
    skipIndexes,
    overWrite,
    verbose
  ) {
    # validate input
    arg <- as.list(environment())
    .subsetgdb_validate_input(arg)

    # check if output exists and overwrite if `overWrite = TRUE`
    .check_output(
      output = output,
      overWrite = overWrite,
      verbose = verbose
    )

    # list tables
    master <- DBI::dbGetQuery(object, "select * from sqlite_master")
    tables.base <- gdb_protected_tables[
      !gdb_protected_tables %in% c("var_ranges", "tmp")
    ]
    tables.anno <- DBI::dbGetQuery(object, "select name from anno")$name
    tables.cohort <- DBI::dbGetQuery(object, "select name from cohort")$name

    # check for inconsistensies:
    # tables in metadata not present in gdb or vice versa
    .gdb_check_tables(
      gdb = object,
      skip = "tmp",
      warning = TRUE,
      return_problems = FALSE
    )

    # keep only specified tables if provided
    if (!is.null(tables)) {
      tables.anno <- tables.anno[tables.anno %in% tables]
      tables.cohort <- tables.cohort[tables.cohort %in% tables]
    }

    # create random output gdb handle
    tmp <- paste(
      c("extract", sample(letters, 28L, replace = TRUE)),
      collapse = ""
    )
    ## ensure handle name is not identical to table name
    while (TRUE) {
      if (!tmp %in% c(tables.base, tables.anno, tables.cohort)) {
        break
      }
      tmp <- paste(
        c("extract", sample(letters, 28L, replace = TRUE)),
        collapse = ""
      )
    }

    # export new var table to output gdb
    DBI::dbExecute(object, sprintf("ATTACH '%s' as %s", output, tmp))
    ## detach on exit
    on.exit(
      try(
        DBI::dbExecute(
          object,
          sprintf("DETACH DATABASE %s", DBI::dbQuoteIdentifier(object, tmp))
        ),
        silent = TRUE
      ),
      add = TRUE
    )

    # subset variants based on `intersection`,`where`,`VAR_id`
    .subsetgdb_subsetvars(
      gdb = object,
      attach_name = tmp,
      intersection = intersection,
      where = where,
      VAR_id = VAR_id
    )

    # copy tables
    .subsetgdb_copy_tables(
      gdb = object,
      attach_name = tmp,
      tables.cohort = tables.cohort,
      tables.anno = tables.anno
    )

    # copy indexes
    if (!skipIndexes) {
      .subsetgdb_copy_indexes(
        gdb = object,
        tables = c(tables.base, tables.cohort, tables.anno),
        master = master,
        attach_name = tmp
      )
    }

    # generate ranged varinfo and meta table (not copied from source)
    gdb <- gdb(output)

    ## ranged varinfo
    addRangedVarinfo(gdb, overwrite = TRUE, verbose = verbose)

    ## meta table
    DBI::dbExecute(gdb, "create table meta (name text,value text)")
    .gdb_populate_meta_table(
      gdb = gdb,
      genomeBuild = getGenomeBuild(object),
      verbose = verbose
    )

    # finished
    message(sprintf(
      "%s\tComplete",
      as.character(round(Sys.time(), units = "secs"))
    ))

    invisible(NULL)
  }
)

.subsetgdb_validate_input <- function(
  args
) {
  # output should be a single filepath
  check_wrapper(check_character, args, "output", length_equal = 1L)

  # intersection
  check_wrapper(check_character, args, "intersection", allow_null = TRUE)

  # VAR_id
  .gdb_check_varid(args[["VAR_id"]])

  # where
  check_wrapper(check_character, args, "where", length_equal = 1L, allow_null = TRUE)

  # tables
  check_wrapper(check_character, args, "tables", allow_null = TRUE)

  # boolean flags
  check_wrapper(check_bool, args, "skipIndexes", length_equal = 1L)
  check_wrapper(check_bool, args, "overWrite", length_equal = 1L)
  check_wrapper(check_bool, args, "verbose", length_equal = 1L)

  invisible(NULL)
}


.subsetgdb_subsetvars <- function(
  gdb,
  attach_name,
  intersection,
  where,
  VAR_id
) {
  # base query
  query <- sprintf(
    "create table %s.var as select distinct var.* from var",
    attach_name
  )

  # add intersections
  if (!is.null(intersection) && length(intersection) > 0L) {
    intersection <- unlist(strsplit(intersection, split = ",", fixed = TRUE))
    for (i in intersection) {
      query <- sprintf("%s inner join %s using (VAR_id)", query, i)
    }
  }

  # where clause + VAR_ids
  if (!is.null(where) && !is.null(VAR_id)) {
    query <- sprintf(
      "%s where %s and var.VAR_id in (%s)",
      query,
      where,
      paste(as.integer(VAR_id), collapse = ",")
    )
  } else if (!is.null(where)) {
    query <- sprintf("%s where %s", query, where)
  } else if (!is.null(VAR_id)) {
    query <- sprintf(
      "%s where var.VAR_id in (%s)",
      query,
      paste(as.integer(VAR_id), collapse = ",")
    )
  }

  # execute
  DBI::dbExecute(gdb, query)

  invisible(NULL)
}

.subsetgdb_copy_tables <- function(
  gdb,
  attach_name,
  tables.cohort,
  tables.anno
) {
  # copy SM table
  DBI::dbExecute(
    gdb,
    sprintf("create table %s.SM as select * from SM", attach_name)
  )

  # copy dosages
  DBI::dbExecute(
    gdb,
    sprintf(
      "create table %s.dosage as select dosage.VAR_id,dosage.GT from dosage inner join %s.var using (VAR_id)",
      attach_name,
      attach_name
    )
  )

  # cohort meta table
  DBI::dbExecute(
    gdb,
    sprintf(
      "create table %s.cohort as select * from cohort where name in (%s)",
      attach_name,
      paste(paste0("'", tables.cohort, "'"), collapse = ",")
    )
  )

  # anno meta table
  DBI::dbExecute(
    gdb,
    sprintf(
      "create table %s.anno as select * from anno where name in (%s)",
      attach_name,
      paste(paste0("'", tables.anno, "'"), collapse = ",")
    )
  )

  # copy cohort tables
  for (cohort_i in tables.cohort) {
    DBI::dbExecute(
      gdb,
      sprintf("create table %s.%s as select * from %s", attach_name, cohort_i, cohort_i)
    )
  }

  # copy anno tables
  for (anno_i in tables.anno) {
    DBI::dbExecute(
      gdb,
      sprintf(
        "create table %s.%s as select %s.* from %s inner join %s.var using (VAR_id)",
        attach_name,
        anno_i,
        anno_i,
        anno_i,
        attach_name
      )
    )
  }

  invisible(NULL)
}

.subsetgdb_copy_indexes <- function(
  gdb,
  tables,
  master,
  attach_name
) {
  for (copied in tables) {
    indexes <- master[
      master$type == "index" &
        master$tbl_name == copied,
    ]$sql

    for (index in indexes) {
      DBI::dbExecute(
        gdb,
        gsub("CREATE INDEX ", sprintf("CREATE INDEX %s.", attach_name), index)
      )
    }
  }

  invisible(NULL)
}

#' writeVcf
#' @rdname writeVcf
#' @aliases writeVcf,gdb-method
#' @usage NULL
#' @export
setMethod(
  "writeVcf",
  signature = "gdb",
  definition = function(
    object,
    output,
    VAR_id = NULL,
    IID = NULL,
    memlimit = 1000L,
    includeGeno = TRUE,
    includeVarId = FALSE,
    overWrite = FALSE,
    verbose = TRUE
  ) {
    # validate input
    arg <- as.list(environment())
    .writevcf_validate_input(arg)
    .check_output(output, overWrite = overWrite, verbose = verbose)

    # open output connection
    tryCatch(
      {
        output_con <- gzfile(output, "w")
      },
      error = function(e) {
        stop(
          sprintf("Could not write to output path '%s'", output),
          call. = FALSE
        )
      }
    )
    on.exit(close(output_con), add = TRUE)

    # sample filtering (if IID is specified)
    tmp <- .writevcf_create_query(
      gdb = object,
      IID = IID,
      VAR_id = VAR_id,
      includeGeno = includeGeno,
      includeVarId = includeVarId,
      verbose = verbose
    )
    IID_bool <- tmp$IID_bool
    SM <- tmp$SM
    query <- tmp$query

    # result handle
    query_handle <- DBI::dbSendQuery(object, query)

    # write VCF
    write("##fileformat=VCFv4.2", output_con)
    if (includeGeno) {
      write(
        sprintf(
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s",
          paste(SM$IID[IID_bool], collapse = "\t")
        ),
        output_con
      )
      while (!DBI::dbHasCompleted(query_handle)) {
        records <- DBI::dbFetch(query_handle, n = memlimit)
        records <- .writevcf_unpackGT(records, IID_bool = IID_bool)
        write.table(
          records,
          output_con,
          append = TRUE,
          sep = "\t",
          row.names = FALSE,
          col.names = FALSE,
          quote = FALSE
        )
      }
    } else {
      write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", output_con)
      while (!DBI::dbHasCompleted(query_handle)) {
        records <- DBI::dbFetch(query_handle, n = memlimit)
        write.table(
          records,
          output_con,
          append = TRUE,
          sep = "\t",
          row.names = FALSE,
          col.names = FALSE,
          quote = FALSE
        )
      }
    }

    # clear
    DBI::dbClearResult(query_handle)

    invisible(NULL)
  }
)

.writevcf_validate_input <- function(
  args
) {
  # output should be a single filepath
  check_wrapper(check_character, args, "output", length_equal = 1L)

  # VAR_id
  .gdb_check_varid(args[["VAR_id"]])

  # IID
  check_wrapper(check_character, args, "IID", allow_null = TRUE)
  
  # memlimit
  check_positive(args[["memlimit"]], arg = "memlimit", length_equal = 1L)

  # boolean flags
  check_wrapper(check_bool, args, "includeGeno", length_equal = 1L)
  check_wrapper(check_bool, args, "includeVarId", length_equal = 1L)
  check_wrapper(check_bool, args, "verbose", length_equal = 1L)

  invisible(NULL)
}

.writevcf_unpackGT <- function(records, IID_bool) {
  # decompress
  GT <- lapply(records[["GT"]], memDecompress, type = "gzip", asChar = FALSE)
  ascii_offset <- as.integer(charToRaw("0"))
  GT <- as.integer(unlist(GT)) - ascii_offset
  GT[GT > 2L] <- NA_real_
  GT <- matrix(
    GT,
    nrow = nrow(records),
    ncol = length(GT) / nrow(records),
    byrow = TRUE
  )

  # subset IIDs
  GT <- GT[, IID_bool, drop = FALSE]

  # convert to vcf format
  GT_char <- matrix(NA_character_, nrow = nrow(GT), ncol = ncol(GT))
  GT_char[is.na(GT)] <- "./."
  GT_char[GT == 0] <- "0/0"
  GT_char[GT == 1] <- "0/1"
  GT_char[GT == 2] <- "1/1"
  records$GT <- NULL
  records <- cbind(records, GT_char)

  # return
  records
}

.writevcf_create_query <- function(
  gdb,
  IID,
  VAR_id,
  includeGeno,
  includeVarId,
  verbose
) {
  id_field <- if (includeVarId) "VAR_id" else "ID"

  # sample filtering (if IID is specified)
  SM <- DBI::dbGetQuery(gdb, "select * from SM")
  if (is.null(IID)) {
    if (verbose) {
      message("No IID provided, all samples will be included in output vcf")
    }
    IID_bool <- rep(TRUE, nrow(SM))
  } else {
    IID_bool <- SM$IID %in% IID
  }
  if (verbose)
    message(sprintf(
      "%s/%s samples to be retained in output",
      sum(IID_bool),
      nrow(SM)
    ))

  # variant retrieval queries with / without genotype data
  if (includeGeno) {
    query <- sprintf(
      "select CHROM, POS, %s, REF, ALT, QUAL, FILTER, '.', FORMAT, GT from var inner join dosage using (VAR_id)",
      id_field
    )
  } else {
    query <- sprintf(
      "select CHROM, POS, %s, REF, ALT, QUAL, FILTER, '.' from var",
      id_field
    )
  }

  # send variant retrieval queries with / without VAR_id filtering
  if (is.null(VAR_id)) {
    if (verbose) {
      message(
        "No VAR_id provided, all variants will be included in output vcf"
      )
    }
    query <- query
  } else {
    if (verbose) {
      message(sprintf("%s variants to be retained in output", length(VAR_id)))
    }
    query <- paste(
      query,
      sprintf("where VAR_id in (%s)", paste(as.integer(VAR_id), collapse = ","))
    )
  }

  list(
    IID_bool = IID_bool,
    SM = SM,
    query = query
  )
}
