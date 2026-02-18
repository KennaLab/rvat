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
  if (!file.exists(path)) {
    stop(sprintf("'%s' doesn't exist.", path), call. = FALSE)
  }
  tryCatch(
    {
      con <- DBI::dbConnect(RSQLite::SQLite(), path)
    },
    error = function(e) {
      stop(sprintf("Invalid aggdb path '%s'", path), call. = FALSE)
    }
  )
  new("aggdb", con)
}


# aggdb misc methods -------------------------------------------------------
#' @rdname aggdb
#' @usage NULL
#' @export
setMethod("close", signature = "aggdb", definition = function(con) {
  DBI::dbDisconnect(con)
})

# aggdb getters/listers -------------------------------------------------------
setMethod("metadata", signature = "aggdb", definition = function(x) {
  meta <- DBI::dbGetQuery(x, "select * from meta")
  meta <- as.list(setNames(meta$value, nm = meta$name))
  meta
})

setMethod("getGdbId", signature = "aggdb", definition = function(object) {
  metadata(object)$gdbId
})

setMethod("getRvatVersion", signature = "aggdb", definition = function(object) {
  metadata(object)$rvatVersion
})

setMethod("listSamples", signature = "aggdb", definition = function(object) {
  IID <- DBI::dbGetQuery(object, "select * from SM")$IID
  IID
})

setMethod("listUnits", signature = "aggdb", definition = function(object) {
  unit <- DBI::dbGetQuery(object, "select unit from aggregates")$unit
  unit
})

setMethod("listParams", signature = "aggdb", definition = function(object) {
  params <- DBI::dbGetQuery(object, "select * from params order by name")
  params <- as.list(setNames(params$value, nm = params$name))
  params
})

# getUnit method -------------------------------------------------------------
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
  if (nrow(aggs) == 0L) {
    return(NULL)
  }

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
  if (verbose) {
    message("  Writing metadata...")
  }

  # write gdb metadata
  DBI::dbWriteTable(
    conn = aggdb,
    name = "meta",
    value = metadata,
    append = TRUE,
    row.names = FALSE
  )

  # write params
  if (verbose) {
    message("  Writing analysis parameters...")
  }
  DBI::dbWriteTable(
    conn = aggdb,
    name = "params",
    value = params,
    append = TRUE,
    row.names = FALSE
  )

  # write sample IDs
  if (verbose) {
    message("  Writing sample manifest (SM)...")
  }
  DBI::dbWriteTable(
    conn = aggdb,
    name = "SM",
    value = data.frame(IID = samples, stringsAsFactors = FALSE),
    append = TRUE,
    row.names = FALSE
  )

  if (verbose) {
    message("Aggregate database initialized successfully.")
  }

  # return aggdb
  aggdb
}
