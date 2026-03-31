#' @rdname spatialClust
#' @usage NULL
#' @export
setMethod(
  "spatialClust",
  signature = "gdb",
  definition = function(
    object,
    output,
    varSetName,
    unitTable,
    unitName,
    windowSize,
    overlap,
    intersection = NULL,
    where = NULL,
    weightName = "1",
    posField = "POS",
    minTry = 5,
    memlimit = 1000L,
    warning = TRUE
  ) {
    .spatialClust_validate_input(as.list(environment()))

    # connect to output
    output <- gzfile(output, "w")
    on.exit(close(output), add = TRUE)

    # base query
    query <- sprintf(
      "select distinct %s as unit, VAR_id, %s as weight, %s as POS from %s",
      unitName,
      weightName,
      posField,
      unitTable
    )

    # add intersection(s) if specified
    if (!is.null(intersection) && length(intersection) > 0L) {
      intersection <- trimws(unlist(strsplit(
        intersection,
        split = ",",
        fixed = TRUE
      )))

      for (intersection_table in intersection) {
        query <- sprintf(
          "%s inner join %s using (VAR_id)",
          query,
          intersection_table
        )
      }
    }

    # add where clause if specified
    if (!is.null(where)) {
      query <- sprintf("%s where %s", query, where)
    }

    # add grouped concat
    query <- sprintf(
      "select unit, group_concat(VAR_id) as VAR_id, group_concat(weight) as weight, '%s' as varSetName, group_concat(POS) as POS from (%s) x group by unit",
      varSetName,
      query
    )

    # write rvat header
    metadata <- list(
      rvatVersion = as.character(packageVersion("rvat")),
      gdbId = getGdbId(object),
      genomeBuild = getGenomeBuild(object),
      creationDate = as.character(round(Sys.time(), units = "secs"))
    )
    .write_rvat_header(
      filetype = "varSetFile",
      metadata = metadata,
      con = output
    )

    handle <- RSQLite::dbSendQuery(object, query)
    on.exit(DBI::dbClearResult(handle), add = TRUE)
    while (!RSQLite::dbHasCompleted(handle)) {
      chunk <- RSQLite::dbFetch(handle, n = memlimit)

      if (nrow(chunk) == 0L) {
        break
      }

      # iterate through chunk
      for (i in seq_len(nrow(chunk))) {
        runDistanceCluster(
          varSet = chunk[i, , drop = FALSE],
          posField = posField,
          windowSize = windowSize,
          overlap = overlap,
          output = output,
          minTry = minTry,
          warning = warning
        )
      }
    }
  }
)


runDistanceCluster <- function(
  varSet,
  posField,
  windowSize,
  overlap,
  output,
  minTry = 5,
  warning = TRUE
) {
  unit <- varSet$unit
  varSetName <- varSet$varSetName

  # build varSet based on VAR_ids and weights
  varSet <- data.frame(
    VAR_id = unlist(strsplit(varSet$VAR_id, ",", fixed = TRUE)),
    weight = unlist(strsplit(varSet$weight, ",", fixed = TRUE)),
    POS = unlist(strsplit(as.character(varSet$POS), ",", fixed = TRUE)),
    stringsAsFactors = FALSE
  )

  # handle case where there are too few variants to cluster
  if (nrow(varSet) < minTry) {
    if (warning) {
      warning(
        sprintf(
          "varSet contains fewer than %s variants ('minTry' parameter), all variants are returned as a single cluster.",
          minTry
        ),
        call. = FALSE
      )
    }
    # write variants as a single cluster for each window size/overlap combination
    for (wi in seq_along(windowSize)) {
      varSeti <- paste(
        paste(unit, "0", sep = "_"),
        paste(varSet$VAR_id, collapse = ","),
        paste(varSet$weight, collapse = ","),
        paste(varSetName, windowSize[wi], overlap[wi], sep = "_"),
        sep = "|"
      )
      write(varSeti, output, append = TRUE)
    }

    # early return
    return(NULL)
  }

  # convert position field to integer based on format
  if (posField == "HGVSc") {
    # extract numeric position from HGVSc notation
    varSet$POS <- as.integer(stringr::str_replace(
      varSet$POS,
      "c.([0-9]*).*",
      "\\1"
    ))
  } else {
    varSet$POS <- as.integer(varSet$POS)
  }

  # perform clustering for each window size/overlap combination
  for (wi in seq_along(windowSize)) {
    # run clustering algorithm
    parts <- distanceCluster(
      pos = varSet$POS,
      names = varSet$VAR_id,
      windowSize = windowSize[wi],
      overlap = overlap[wi],
      return_what = "names",
      write_output = FALSE,
      path_output = NA,
      warning = warning
    )

    # write each cluster partition to output
    for (part in names(parts)) {
      # collapse VAR_ids, weights, and positions for this partition
      varSeti <- apply(
        varSet[varSet$VAR_id %in% parts[[part]], ],
        2,
        paste,
        collapse = ","
      )
      # format output line
      varSeti <- paste(
        paste(unit, part, sep = "_"),
        varSeti["VAR_id"],
        varSeti["weight"],
        paste(varSetName, windowSize[wi], overlap[wi], sep = "_"),
        sep = "|"
      )
      write(varSeti, output, append = TRUE)
    }
  }
}


# distanceCluster
# Adapted from https://github.com/heidefier/cluster_wgs_data (Fier, GenetEpidemiol, 2017).
# Partitions a supplied set of variant until the mean distance between variants
#   equals the median (as would be expected for a homogeneous poisson process).
distanceCluster <- function(
  pos,
  names,
  windowSize = 100,
  overlap = 50,
  return_what = "names",
  write_output = TRUE,
  path_output = NA,
  warning = TRUE
) {
  if (!(return_what %in% c("names", "positions"))) {
    stop(
      "Please enter 'names' OR 'positions' as argument for return_what",
      call. = FALSE
    )
  }

  if (write_output && is.na(path_output)) {
    stop("Please specify a valid output name in path_output", call. = FALSE)
  }

  b <- as.numeric(pos)
  comb <- data.frame(b = b, names = names, stringsAsFactors = FALSE)
  comb <- comb[order(b), ]
  b <- as.numeric(comb[, 1])
  names <- comb[, 2]
  bdiff <- diff(b)

  #observed mean values in the windows

  obsMeanWindow <- zoo::rollapply(
    bdiff,
    windowSize,
    mean,
    by = overlap,
    align = "left"
  )
  if (length(obsMeanWindow) == 0L) {
    if (warning) {
      warning(
        "No valid partitioning for specified parameters, all variants are returned as a single cluster.",
        call. = FALSE
      )
    }
    return(list("0" = names))
  }

  #theoretical mean value given medians in the windows
  ws <- as.numeric(windowSize)
  sh <- as.numeric(overlap)

  theoMeanWindow <- zoo::rollapply(
    bdiff,
    ws,
    function(x) 1 / (log(2) / median(x)),
    by = sh,
    align = "left"
  )

  #actual positions
  raPositions <- zoo::rollapply(
    b,
    ws + 1,
    function(x) x,
    by = sh,
    align = "left"
  )

  maximas <- vector("list", length(obsMeanWindow))
  for (i in seq_along(obsMeanWindow)) {
    if (obsMeanWindow[i] <= theoMeanWindow[i]) {
      maximas[[i]] <- NA
    } else {
      ord <- order(diff(raPositions[i, ]), decreasing = TRUE)
      j <- 1
      excl <- mean(diff(raPositions[i, ])[-c(ord[1])])
      while (theoMeanWindow[i] < excl) {
        j <- j + 1
        excl <- mean(diff(raPositions[i, ])[-c(ord[1:j])]) #The largest distances are removed until the median-informed theoretical means are larger than the mean of the edited window
      }
      j <- j - 1
      maximas[[i]] <- raPositions[i, ][c(ord[1:j]) + 1]
    }
  }
  maxima <- unlist(maximas)
  maxima <- maxima[!is.na(maxima)]
  maximaf <- sort(unique(maxima))

  if (length(maximaf) == 0L) {
    if (warning) {
      warning(
        "No valid partitioning for specified parameters, all variants are returned as a single cluster.",
        call. = FALSE
      )
    }
    return(list("0" = names))
  }

  if (!is.na(path_output)) {
    file.create(path_output)
  }
  if (return_what == "positions") {
    lo_final <- split(b, findInterval(b, maximaf))
    ll <- as.numeric(lapply(lo_final, length))
    if (return_what == "positions" && write_output) {
      lapply(lo_final, write, path_output, append = TRUE, ncolumns = max(ll))
    }
  } else {
    lo_final <- split(names, findInterval(b, maximaf))
    ll <- as.numeric(lapply(lo_final, length))
    if (return_what == "names" && write_output) {
      lapply(lo_final, write, path_output, append = TRUE, ncolumns = max(ll))
    }
  }

  lo_final
}

.spatialClust_validate_input <- function(args) {
  # type checks
  check_wrapper(check_character, args, "output", length_equal = 1L)
  check_wrapper(check_character, args, "varSetName", length_equal = 1L)
  check_wrapper(check_character, args, "unitTable", length_equal = 1L)
  check_wrapper(check_character, args, "unitName", length_equal = 1L)
  stopifnot(is.numeric(args[["windowSize"]]))
  stopifnot(is.numeric(args[["overlap"]]))
  check_wrapper(check_character, args, "intersection", allow_null = TRUE)
  check_wrapper(check_character, args, "where", allow_null = TRUE)
  check_wrapper(check_character, args, "weightName", length_equal = 1L)
  check_wrapper(check_character, args, "posField", length_equal = 1L)
  check_wrapper(check_number_whole, args, "minTry", length_equal = 1L)
  check_wrapper(check_bool, args, "warning", length_equal = 1L)

  # length of `windowSize` should match that of `overlap`
  if (length(args[["windowSize"]]) != length(args[["overlap"]])) {
    stop(
      sprintf(
        "Number of provided window sizes (%s) does not match number of provided overlaps (%s)",
        paste(args[["windowSize"]], collapse = ","),
        paste(args[["overlap"]], collapse = ",")
      ),
      call. = FALSE
    )
  }
}
