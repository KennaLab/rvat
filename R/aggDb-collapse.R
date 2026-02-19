#' @include allClasses.R
#' @include allGenerics.R
#' @include allInternalData.R

#' @rdname collapseAggDbs
#' @usage NULL
#' @export
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
    check_length(output, equal = 1L, allow_null = TRUE)
    if (!is.null(output)) {
      .check_output(output, overWrite = overWrite, verbose = verbose)
    }
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