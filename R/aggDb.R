#' @include allClasses.R
#' @include allGenerics.R
#' @include allInternalData.R
#' @include asserts-obj-type.R
#' @include asserts-types-check.R

# aggdb constructor ----------------------------------------------------
#' @rdname aggdb
#' @usage NULL
#' @export
aggdb <- function(path) {
  if (!file.exists(path))
    stop(sprintf("'%s' doesn't exist.", path), call. = FALSE)
  tryCatch(
    {
      con <- DBI::dbConnect(DBI::dbDriver("SQLite"), path)
    },
    error = function(e) {
      stop(sprintf("Invalid aggdb path '%s'", path), call. = FALSE)
    }
  )
  new("aggdb", con)
}


# aggdb methods -------------------------------------------------------
setMethod("close", signature = "aggdb", definition = function(con) {
  DBI::dbDisconnect(con)
})

setMethod("metadata", signature = "aggdb", definition = function(x) {
  meta <- DBI::dbGetQuery(x, "select * from meta")
  meta <- as.list(setNames(meta$value, nm = meta$name))
  meta
})


setMethod("getGdbId", signature = "aggdb", definition = function(object) {
  metadata(object)$gdbId
})

setMethod("getRvatVersion", signature="aggdb",
          definition=function(object)
          {
            metadata(object)$gdbId
          })

setGeneric("listSamples", function(object) standardGeneric("listSamples"))
setMethod("listSamples", signature = "aggdb", definition = function(object) {
  IID <- DBI::dbGetQuery(object, "select * from SM")$IID
  IID
})

setMethod("listUnits", signature = "aggdb", definition = function(object) {
  unit <- DBI::dbGetQuery(object, "select unit from aggregates")$unit
  unit
})

setGeneric("listParams", function(object) standardGeneric("listParams"))
setMethod("listParams", signature = "aggdb", definition = function(object) {
  params <- DBI::dbGetQuery(object, "select * from params order by name")
  params <- as.list(setNames(params$value, nm = params$name))
  params
})

setMethod("getUnit", signature = "aggdb", definition = function(object, unit) {
  # check if unit is of type character and length > 0
  if (!is.character(unit) || length(unit) == 0L) {
    stop("'unit' must be a non-empty character vector.", call. = FALSE)
  }

  # retrieve units
  query <- sprintf(
    "SELECT unit, aggregate FROM aggregates WHERE unit IN (%s)",
    paste(sQuote(unit, "'"), collapse = ",")
  )
  aggs <- DBI::dbGetQuery(object, query)
  unit_found <- aggs$unit
  if (!all(unit %in% unit_found)) {
    warning(
      sprintf(
        "%s/%s of specified units could not be found.",
        sum(!unit %in% unit_found),
        length(unit)
      ),
      call. = FALSE
    )
  }
  if (nrow(aggs) == 0L) return(NULL)

  # decode, decompress, and unserialize blobs
  aggs <- lapply(aggs$aggregate, function(blob) {
    unserialize(memDecompress(blob, type = "gzip"))
  })

  # construct matrix
  aggs <- do.call(rbind, aggs)
  rownames(aggs) <- unit_found
  colnames(aggs) <- listSamples(object)

  # return
  aggs
})


# aggdbList constructor -----------------------------------------------------
aggdbList <- function(filelist, checkDups = TRUE) {
  # input validation
  check_character(filelist)
  check_length(filelist, min = 1L)
  check_bool(checkDups)

  # connect to aggdbs
  aggdb_list <- lapply(
    filelist,
    FUN = aggdb
  )
  on.exit(lapply(
    aggdb_list,
    function(con) if (DBI::dbIsValid(con)) try(close(con), silent = TRUE)
  ))

  # check whether duplicate units are included
  units_all <- unlist(lapply(aggdb_list, listUnits))
  if (anyDuplicated(units_all) != 0L) {
    if (checkDups) {
      stop(
        sprintf(
          "The following units are duplicated across databases: %s",
          paste(units_all[duplicated(units_all)], collapse = ",")
        ),
        call. = FALSE
      )
    } else {
      units_all <- paste0("X", seq_along(units_all))
    }
  }

  # check whether sample IDs are identical across dbs
  samples <- lapply(aggdb_list, listSamples)
  check_identical_samples <- unlist(lapply(samples, FUN = function(x) {
    identical(samples[[1L]], x)
  }))
  if (!all(check_identical_samples)) {
    stop(
      "aggdbs should include identical samples (in the same order).",
      call. = FALSE
    )
  }
  samples <- samples[[1L]]

  # check metadata
  metadata <- lapply(aggdb_list, metadata)
  rvatversion <- unique(unlist(lapply(
    metadata,
    function(x) x[["rvatVersion"]]
  )))
  gdbid <- unique(unlist(lapply(
    metadata,
    function(x) x[["gdbId"]]
  )))
  genomebuild <- unique(unlist(lapply(
    metadata,
    function(x) x[["genomeBuild"]]
  )))
  params <- lapply(
    aggdb_list,
    listParams
  )
  check_identical_params <- unlist(lapply(params, FUN = function(x) {
    identical(params[[1L]], x)
  }))
  if (!all(check_identical_params)) {
    stop(
      "Non-identical parameters were used to generate the input aggdbs.",
      call. = FALSE
    )
  }
  if (length(rvatversion) > 1L) {
    stop(
      "AggregateFiles were generated using different rvat versions.",
      call. = FALSE
    )
  }
  if (length(gdbid) > 1L) {
    stop("AggregateFiles were generated from different gdbs", call. = FALSE)
  }

  metadata <- list(
    rvatVersion = rvatversion,
    gdbId = gdbid,
    genomeBuild = if (length(genomebuild) > 0L) genomebuild[1L] else
      NA_character_,
    creationDate = as.character(round(Sys.time(), units = "secs"))
  )

  # return aggdblist
  new(
    "aggdbList",
    paths = filelist,
    units = units_all,
    samples = samples,
    params = params[[1L]],
    metadata = metadata
  )
}


# aggdbList methods -------------------------------------------------------
setMethod("show", signature = "aggdbList", definition = function(object) {
  cat(sprintf(
    "aggdbList object\naggdbs: %s\nSamples: %s\nUnits: %s\n",
    length(object@paths),
    length(object@samples),
    length(object@units)
  ))
})

#' @rdname aggdbList
#' @usage NULL
#' @export
setMethod("listUnits", signature = "aggdbList", definition = function(object) {
  object@units
})

#' @rdname aggdbList
#' @usage NULL
#' @export
setMethod(
  "listSamples",
  signature = "aggdbList",
  definition = function(object) {
    object@samples
  }
)

#' @rdname aggdbList
#' @usage NULL
#' @export
setMethod("length", signature = "aggdbList", definition = function(x) {
  length(x@paths)
})

#' @rdname aggdbList
#' @usage NULL
#' @export
setMethod("metadata", signature = "aggdbList", definition = function(x) {
  x@metadata
})

setMethod("listParams", signature = "aggdbList", definition = function(object) {
  object@params
})

setMethod(
  "mergeAggDbs",
  signature = signature(object = "aggdbList"),
  definition = function(
    object,
    output = NULL,
    overWrite = FALSE,
    verbose = TRUE
  ) {
    # input validation
    if (!is.character(output) || length(output) != 1L || !nzchar(output)) {
      stop("`output` must be a single, non-empty file path.", call. = FALSE)
    }
    check_bool(overWrite)
    check_bool(verbose)

    # initalize aggdb
    metadata <- metadata(object)
    metadata <- data.frame(
      name = names(metadata),
      value = as.character(unlist(metadata)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    params <- listParams(object)
    params <- data.frame(
      name = names(params),
      value = as.character(unlist(params, use.names = FALSE)),
      stringsAsFactors = FALSE
    )
    samples <- listSamples(object)

    aggdb <- .init_aggdb(
      output = output,
      params = params,
      samples = samples,
      metadata = metadata,
      overWrite = overWrite,
      verbose = verbose
    )
    on.exit(close(aggdb), add = TRUE)

    # merge aggdbs
    for (path in object@paths) {
      if (verbose) message(sprintf("Merging '%s'", basename(path)))

      # copy aggregates
      tryCatch(
        {
          DBI::dbExecute(aggdb, sprintf("ATTACH DATABASE '%s' AS src", path))
          DBI::dbExecute(
            aggdb,
            paste0(
              "INSERT INTO aggregates (unit, aggregate) ",
              "SELECT unit, aggregate FROM src.aggregates"
            )
          )
        },
        error = function(e) {
          warning(
            sprintf(
              "Failed to merge data from '%s': %s",
              path,
              e$message
            ),
            call. = FALSE
          )
        },
        # detach
        finally = {
          try(DBI::dbExecute(aggdb, "DETACH DATABASE src"), silent = TRUE)
        }
      )
    }

    if (verbose)
      message(sprintf(
        "Merge complete. New aggregate database created at: %s",
        output
      ))

    # check if all units are included
    check_units <- DBI::dbGetQuery(aggdb, "select unit from aggregates")$unit
    if (!all(listUnits(object) %in% check_units)) {
      stop("Not all expected units are present in merged aggdb.", call. = FALSE)
    }

    # return nothing
    invisible(NULL)
  }
)

# collapseAggDbs -------------------------------------------------------------
setMethod(
  "collapseAggDbs",
  signature = signature(object = "aggdbList"),
  definition = function(
    object,
    output = NULL,
    overWrite = FALSE,
    verbose = TRUE
  ) {
    # input validation
    check_character(output, allow_null = TRUE)
    if (!is.null(output)) {
      .check_output(output, overWrite = overWrite, verbose = verbose)
    }
    check_length(output, equal = 1L, allow_null = TRUE)
    check_bool(verbose)

    # initialize
    samples <- listSamples(object)
    agg_merged <- vector(mode = "numeric", length = length(samples))

    # loop through databases
    for (i in seq_along(object@paths)) {
      path <- object@paths[i]
      if (verbose) {
        message(sprintf(
          "Processing aggdb %d/%d: '%s'",
          i,
          length(object@paths),
          basename(path)
        ))
      }

      con <- NULL
      tryCatch(
        {
          # connect
          con <- aggdb(path)

          # retrieve aggregates
          agg_current <- getUnit(con, unit = listUnits(con))

          # sum aggregates
          agg_current_merged <- colSums(agg_current, na.rm = TRUE)

          # add to total
          agg_merged <- agg_merged + agg_current_merged
        },
        error = function(e) {
          stop(
            sprintf(
              "Could not process '%s': Error: %s",
              basename(path),
              e$message
            ),
            call. = FALSE
          )
        },
        finally = {
          if (!is.null(con) && DBI::dbIsValid(con)) {
            close(con)
          }
        }
      )
    }

    # finalize
    agg_merged <- data.frame(
      IID = samples,
      aggregate = agg_merged,
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    # return data.frame if output is not specified
    if (is.null(output)) {
      return(agg_merged)
    } else {
      if (verbose) {
        message(sprintf("Writing merged aggregates to: %s", output))
      }

      # write to output if specified
      output_con <- gzfile(output, "wb")
      on.exit(close(output_con), add = TRUE)
      write.table(
        agg_merged,
        file = output_con,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
      )

      invisible(NULL)
    }
  }
)

# helpers --------------------------------------------------------------------

.init_aggdb <- function(
  output,
  metadata,
  params,
  samples,
  overWrite,
  verbose = TRUE
) {
  # check output filepath
  .check_output(
    output = output,
    overWrite = overWrite,
    verbose = verbose
  )

  # connect to db
  tryCatch(
    {
      aggdb <- DBI::dbConnect(RSQLite::SQLite(), output)
    },
    error = function(e) {
      stop(sprintf("Invalid aggdb path '%s'", output), call. = FALSE)
    }
  )
  aggdb <- new("aggdb", aggdb)

  if (verbose) {
    message(sprintf("Initializing new aggregate database at: %s", output))
  }

  # create schema
  DBI::dbExecute(
    aggdb,
    "CREATE TABLE meta (name TEXT UNIQUE PRIMARY KEY, value TEXT)"
  )
  DBI::dbExecute(
    aggdb,
    "CREATE TABLE params (name TEXT UNIQUE PRIMARY KEY, value TEXT)"
  )
  DBI::dbExecute(
    aggdb,
    "CREATE TABLE aggregates (unit TEXT PRIMARY KEY, aggregate BLOB)"
  )
  DBI::dbExecute(
    aggdb,
    "CREATE TABLE SM (IID TEXT PRIMARY KEY)"
  )

  # write gdb metadata
  if (verbose) message("  Writing metadata...")

  # write gdb metadata
  DBI::dbWriteTable(
    conn = aggdb,
    name = "meta",
    value = metadata,
    append = TRUE,
    row.names = FALSE
  )

  # write params
  if (verbose) message("  Writing analysis parameters...")
  DBI::dbWriteTable(
    conn = aggdb,
    name = "params",
    value = params,
    append = TRUE,
    row.names = FALSE
  )

  # write sample IDs
  if (verbose) message("  Writing sample manifest (SM)...")
  DBI::dbWriteTable(
    conn = aggdb,
    name = "SM",
    value = data.frame(IID = samples, stringsAsFactors = FALSE),
    append = TRUE,
    row.names = FALSE
  )

  if (verbose) message("Aggregate database initialized successfully.")

  # return aggdb
  aggdb
}
