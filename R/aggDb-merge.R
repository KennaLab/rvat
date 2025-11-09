#' @include allClasses.R
#' @include allGenerics.R
#' @include allInternalData.R

setMethod(
  "mergeAggDbs",
  signature = signature(object = "aggdbList"),
  definition = function(
    object,
    output,
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
      if (verbose) {
        message(sprintf("Merging '%s'", basename(path)))
      }

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
          stop(
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

    if (verbose) {
      message(sprintf(
        "Merge complete. New aggregate database created at: %s",
        output
      ))
    }

    # check if all units are included
    check_units <- DBI::dbGetQuery(aggdb, "select unit from aggregates")$unit
    if (!all(listUnits(object) %in% check_units)) {
      stop("Not all expected units are present in merged aggdb.", call. = FALSE)
    }

    # return nothing
    invisible(NULL)
  }
)
