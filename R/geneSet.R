#' @rdname geneSet
#' @usage NULL
#' @export
geneSet <- function(geneSetName, units, w = NULL, metadata = "") {
  # input validation
  check_character(geneSetName)
  check_length(geneSetName, equal = 1L)
  check_character(units)
  check_length(units, equal = 1L)

  # handle weights
  n_units <- length(unlist(strsplit(units, split = ",", fixed = TRUE)))
  if (is.null(w) || is.na(w)) {
    w <- paste(
      rep(1.0, n_units),
      collapse = ","
    )
  } else {
    if (!is.character(w) || length(w) != 1L) {
      stop(
        "`w` must be a single character string (comma-separated) or NULL.",
        call. = FALSE
      )
    }

    n_w <- length(unlist(strsplit(w, split = ",", fixed = TRUE)))
    if (n_w != n_units) {
      stop(
        sprintf(
          "Number of weights (%d) does not match number of units (%d) for geneSet '%s'.",
          n_w,
          n_units,
          geneSetName
        ),
        call. = FALSE
      )
    }
  }

  if (is.null(metadata) || is.na(metadata)) {
    metadata <- ""
  }

  new(
    "geneSet",
    geneSetName = geneSetName,
    units = units,
    w = w,
    metadata = metadata
  )
}

#' @rdname geneSet
#' @usage NULL
#' @export
setMethod("metadata", signature = "geneSet", definition = function(x) {
  x@metadata
})

#' @rdname geneSet
#' @usage NULL
#' @export
setMethod("length", signature = "geneSet", definition = function(x) {
  length(listUnits(x))
})

#' @rdname geneSet
#' @usage NULL
#' @export
setMethod("as.data.frame", signature = "geneSet", definition = function(x) {
  data.frame(
    geneSetName = x@geneSetName,
    units = x@units,
    w = x@w,
    metadata = x@metadata,
    stringsAsFactors = FALSE
  )
})

#' @rdname geneSet
#' @usage NULL
#' @export
setMethod("listUnits", signature = "geneSet", definition = function(object) {
  unlist(strsplit(object@units, split = ",", fixed = TRUE))
})

#' @rdname geneSet
#' @usage NULL
#' @export
setMethod("listWeights", signature = "geneSet", definition = function(object) {
  x <- unlist(strsplit(object@w, split = ",", fixed = TRUE))
  names(x) <- listUnits(object)
  x
})


# geneSetList ------------------------------------------------------------------

#' geneSetList
#' @rdname geneSetList
#' @usage NULL
#' @export
geneSetList <- function(geneSets, metadata = list()) {
  geneSetNames <- vapply(geneSets, function(x) x@geneSetName, character(1L))
  new(
    "geneSetList",
    geneSets = geneSets,
    geneSetNames = geneSetNames,
    metadata = metadata
  )
}

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("names", signature = "geneSetList", definition = function(x) {
  x@geneSetNames
})

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("length", signature = "geneSetList", definition = function(x) {
  length(x@geneSets)
})

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("lengths", signature = "geneSetList", definition = function(x) {
  lengths(as.list(x))
})

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("show", signature = "geneSetList", definition = function(object) {
  cat(sprintf("geneSetList\nContains %s sets\n", length(object)))
})

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("metadata", signature = "geneSetList", definition = function(x) {
  x@metadata
})

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod(
  "listMetadata",
  signature = "geneSetList",
  definition = function(object) {
    unlist(lapply(object@geneSets, FUN = metadata))
  }
)

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod(
  "listGeneSets",
  signature = "geneSetList",
  definition = function(object) {
    names(object)
  }
)

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod(
  "listUnits",
  signature = "geneSetList",
  definition = function(object) {
    unique(unlist(as.list(object)))
  }
)

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("[[", c("geneSetList", "ANY", "missing"), function(x, i, j, ...) {
  x@geneSets[[i, ...]]
})

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod(
  "[",
  c("geneSetList", "ANY", "ANY"),
  function(x, i, j, ..., drop = TRUE) {
    if (1L != length(drop) || (!missing(drop) && drop)) {
      warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'", call. = FALSE)
    }

    if (missing(i) && missing(j)) {
      return(x)
    }

    initialize(
      x,
      geneSets = x@geneSets[i, ...],
      geneSetNames = x@geneSetNames[i, ...]
    )
  }
)


#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("sort", "geneSetList", function(x) {
  sorted <- stringr::str_sort(names(x))
  new(
    "geneSetList",
    geneSets = x@geneSets[match(sorted, names(x))],
    geneSetNames = sorted
  )
})


#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("as.data.frame", signature = "geneSetList", definition = function(x) {
  if (length(x) == 0L) {
    return(data.frame(
      geneSetName = character(0L),
      units = character(0L),
      w = character(0L),
      metadata = character(0L),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(x@geneSets, FUN = as.data.frame))
})


#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("as.list", signature = "geneSetList", definition = function(x) {
  names <- names(x)
  x <- strsplit(
    vapply(seq_along(x), function(i) x[[i]]@units, character(1L)),
    split = ",",
    fixed = TRUE
  )
  names(x) <- names
  x
})

#' @keywords internal
#' @noRd
setMethod(
  "mapToMatrix",
  signature = c("geneSetList", "rvbResult"),
  definition = function(object, results, ID = "unit", sparse = TRUE) {
    dat <- lapply(
      as.list(object),
      FUN = function(gene_set, IDs) {
        x <- IDs %in% gene_set
      },
      IDs = unique(as.character(results[[ID]]))
    )
    dat <- do.call(cbind, dat)
    rownames(dat) <- unique(as.character(results[[ID]]))
    colnames(dat) <- names(object)
    if (sparse) as(dat, "sparseMatrix") else dat
  }
)


#' @rdname getGeneSet
#' @usage NULL
#' @export
setMethod(
  "getGeneSet",
  signature = "geneSetList",
  definition = function(object, geneSet = NULL, unit = NULL) {
    # input validation
    check_character(geneSet, allow_null = TRUE)
    check_character(unit, allow_null = TRUE)
    if (is.null(geneSet) && is.null(unit)) {
      stop(
        "At least one of `geneSet` or `unit` should be specified.",
        call. = FALSE
      )
    }

    # filter geneSets
    if (!is.null(geneSet)) {
      if (!all(geneSet %in% names(object))) {
        warning(
          "Not all specified geneSets are present in the geneSetList, ",
          "use `listGeneSets()` to check which geneSets are available.",
          call. = FALSE
        )
      }
      object <- object[names(object) %in% geneSet]
    }

    # filter units
    if (!is.null(unit)) {
      keep_idx <- vapply(
        object@geneSets,
        function(x) {
          any(unit %in% listUnits(x))
        },
        logical(1L)
      )
      object <- object[keep_idx]
    }

    object
  }
)


#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod(
  "dropUnits",
  signature = "geneSetList",
  definition = function(object, unit) {
    # input validation
    check_character(unit)

    # remove specified units
    object@geneSets <- lapply(object@geneSets, function(x) {
      # current units and weights
      units <- listUnits(x)
      weights <- listWeights(x)

      # units to keep
      keep <- !units %in% unit

      # update
      x@units <- paste(units[keep], collapse = ",")
      x@w <- paste(weights[keep], collapse = ",")

      x
    })

    # remove empty objects
    object <- object[lengths(object) > 0L]

    object
  }
)

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod(
  "remapIDs",
  signature = "geneSetList",
  definition = function(
    object,
    dict,
    targets = NULL,
    duplicate_ids = c("keep_all", "keep_first"),
    verbose = TRUE
  ) {
    # input validation
    duplicate_ids <- match.arg(duplicate_ids)
    check_character(targets, allow_null = TRUE)
    check_bool(verbose)

    dict <- .remapIDs_validate_dict(dict, object, verbose)
    dict <- .remapIDs_handle_duplicates(dict, targets, duplicate_ids)

    # remap the IDs in each geneset based on the dictionary
    object@geneSets <- lapply(
      object@geneSets,
      FUN = function(x, dct) {
        mapping <- data.frame(
          original_id = listUnits(x),
          w = listWeights(x),
          stringsAsFactors = FALSE
        ) %>%
          dplyr::left_join(dct, by = "original_id") %>%
          dplyr::filter(!is.na(new_id))

        initialize(
          x,
          units = paste(mapping$new_id, collapse = ","),
          w = paste(mapping$w, collapse = ",")
        )
      },
      dct = dict
    )

    # return remapped object
    object
  }
)

.remapIDs_validate_dict <- function(dict, object, verbose) {
  if (ncol(dict) != 2L) {
    stop("`dict` should be a data.frame with two columns", call. = FALSE)
  }
  colnames(dict) <- c("original_id", "new_id")
  dict$original_id <- as.character(dict$original_id)
  dict$new_id <- as.character(dict$new_id)

  # number of IDs in geneSetList that are present in dictionary
  original_ids <- listUnits(object)
  if (verbose) {
    message(sprintf(
      "%s/%s IDs in the geneSetList are present in the linker file.",
      sum(original_ids %in% dict$original_id),
      length(original_ids)
    ))
  }

  # return filtered dict
  dict <- dict[dict$original_id %in% original_ids, , drop = FALSE]
  dict
}

.remapIDs_handle_duplicates <- function(dict, targets, duplicate_ids) {
  # handle IDs that map to multiple IDs
  if (!is.null(targets)) {
    dict$target <- dict$new_id %in% targets
  } else {
    dict$target <- TRUE
  }
  duplicates <- dict %>%
    dplyr::count(original_id) %>%
    dplyr::filter(n > 1) %>%
    dplyr::ungroup()

  # handle duplicates
  ## if targets are specified, only need to handle duplicates which
  ## are both in the target list.
  if (!is.null(targets)) {
    check <- dict %>%
      dplyr::filter(original_id %in% duplicates$original_id, target) %>%
      dplyr::count(original_id)
    check_1 <- check %>% dplyr::filter(n == 1)
    duplicates1 <- dict %>%
      dplyr::filter(original_id %in% check_1$original_id, target)
    check_multiple <- check %>% dplyr::filter(n > 1)
  } else {
    check_multiple <- duplicates
    duplicates1 <- data.frame(
      original_id = character(),
      new_id = character(),
      target = logical(),
      stringsAsFactors = FALSE
    )
  }

  # if duplicate_ids == "keep_all", simply keep all duplicates
  # if duplicate_ids == "keep_first", keep first
  if (duplicate_ids == "keep_all") {
    duplicates2 <- dict %>%
      dplyr::filter(original_id %in% check_multiple$original_id, target)
  } else if (duplicate_ids == "keep_first") {
    duplicates2 <- dict %>%
      dplyr::filter(original_id %in% check_multiple$original_id, target) %>%
      dplyr::group_by(original_id) %>%
      dplyr::mutate(keep = new_id[1]) %>%
      dplyr::ungroup() %>%
      dplyr::filter(new_id == keep) %>%
      dplyr::select(-keep)
  }
  duplicates_keep <- dplyr::bind_rows(duplicates1, duplicates2)

  # subset dictionary based on which duplicates to keep
  dict <- dplyr::bind_rows(
    dict %>% dplyr::filter(!original_id %in% duplicates$original_id),
    duplicates_keep
  )

  dict
}


#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("write", "geneSetList", function(x, file = "data", append = FALSE) {
  # open connection
  if (append) {
    out <- gzfile(file, "a")
  } else {
    out <- gzfile(file, "w")
  }
  on.exit(close(out), add = TRUE)

  # write header with metadata (if append = FALSE)
  if (!append || !file.exists(file)) {
    metadata <- metadata(x)
    metadata$rvatVersion <- as.character(packageVersion("rvat"))
    metadata$creationDate <- as.character(round(Sys.time(), units = "secs"))
    .write_rvat_header(filetype = "geneSetFile", metadata = metadata, con = out)
  }

  write.table(
    as.data.frame(x),
    out,
    sep = "|",
    append = append,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
})

# geneSetFile ------------------------------------------------

#' geneSetFile
#'
#' Initialize \code{\link{geneSetFile-class}} object
#' @param path Path to geneSet file.
#' @param memlimit Maximum number of records to load at a time.
#' @export
geneSetFile <- function(path, memlimit = 5000L) {
  # input validation
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a single, non-empty file path string.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("File does not exist '%s'", path), call. = FALSE)
  }
  check_number_whole(memlimit, min = 1)

  # read in metadata
  header_info <- .parse_rvat_header(
    path = path,
    expected_metadata = metadata_genesets,
    expected_filetype = "geneSetFile",
    return_skip = TRUE,
    n = length(metadata_genesets) + 1L # file description + metadata
  )
  metadata <- header_info$metadata
  skip <- header_info$skip

  # establish connection
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)

  # skip header if present
  if (skip > 0L) {
    readLines(con, n = skip)
  }

  sets <- list()
  while (length(chunk <- readLines(con, n = memlimit)) > 0L) {
    sets[[length(sets) + 1L]] <-
      stringi::stri_split_fixed(
        chunk,
        pattern = "|",
        n = 4L,
        simplify = TRUE
      )[, 1L]
  }
  sets <- unlist(sets, use.names = FALSE)
  if (is.null(sets)) {
    sets <- character(0L)
  }

  # return geneSetFile
  new("geneSetFile", path = path, sets = sets, metadata = metadata)
}

#' @rdname geneSetFile
#' @usage NULL
#' @export
setMethod("show", signature = "geneSetFile", definition = function(object) {
  cat(sprintf(
    "geneSetFile\nPath: %s\nSets: %s",
    object@path,
    length(object@sets)
  ))
})

#' @rdname geneSetFile
#' @usage NULL
#' @export
setMethod("names", signature = "geneSetFile", definition = function(x) {
  x@sets
})

#' @rdname geneSetFile
#' @usage NULL
#' @export
setMethod(
  "listGeneSets",
  signature = "geneSetFile",
  definition = function(object) {
    names(object)
  }
)

#' @rdname geneSetFile
#' @usage NULL
#' @export
setMethod("length", signature = "geneSetFile", definition = function(x) {
  length(names(x))
})

#' @rdname geneSetFile
#' @usage NULL
#' @export
setMethod("metadata", signature = "geneSetFile", definition = function(x) {
  x@metadata
})


#' @rdname getGeneSet
#' @usage NULL
#' @export
setMethod(
  "getGeneSet",
  signature = "geneSetFile",
  definition = function(object, geneSet) {
    # input validation
    check_character(geneSet)

    if (!all(geneSet %in% object@sets)) {
      warning(
        "Not all specified geneSets are present in the geneSetFile, ",
        "use `listGeneSets()` to check which geneSets are available.",
        call. = FALSE
      )
    }

    # determine indices to scan
    indices <- sort(which(object@sets %in% geneSet))

    # calculate relative skips for scanning
    relative_indices <- indices
    if (length(indices) > 1L) {
      relative_indices[2L:length(indices)] <- (dplyr::lead(indices) - indices)[
        1L:(length(indices) - 1L)
      ]
    }

    # check header and return number of lines to skip
    skip <- .genesetfile_read_metadata(object)

    # establish connection
    con <- gzfile(object@path, "rb")
    on.exit(close(con), add = TRUE)

    # skip header if present
    if (skip > 0L) {
      readLines(con, n = skip)
    }

    # retrieve genesets
    genesets <- lapply(relative_indices, FUN = function(i) {
      line_content <- scan(
        con,
        skip = i - 1L,
        nlines = 1L,
        what = "character",
        quiet = TRUE
      )
      line_content <- unlist(strsplit(line_content, split = "|", fixed = TRUE))
      if (length(line_content) < 2) {
        stop(
          sprintf(
            paste0(
              "Malformed line encountered at unit: %s. ",
              "Expected 4 fields, encountered %s."
            ),
            if (length(line_content) > 0L) line_content[1L] else "unknown",
            length(line_content)
          ),
          call. = FALSE
        )
      }
      geneSet(
        geneSetName = line_content[1],
        units = line_content[2],
        w = line_content[3],
        metadata = line_content[4]
      )
    })
    genesets <- geneSetList(genesets, metadata = metadata(object))

    # throw error if unrequested geneSets are included
    if (!all(listGeneSets(genesets) %in% geneSet)) {
      stop("Unspecified genesets retrieved, something's wrong.", call. = FALSE)
    }

    # return genesets
    genesets
  }
)

.genesetfile_read_metadata <- function(object, return_header = FALSE) {
  header <- readLines(object@path, n = length(metadata_genesets) + 10L)
  skip <- sum(startsWith(header, "#"))
  if (skip > length(metadata_genesets) + 1L) {
    stop("File contains more header lines than expected.", call. = FALSE)
  }

  if (return_header) {
    invisible(NULL)
  } else {
    skip
  }
}


#' @rdname geneSetFile
#' @usage NULL
#' @export
setMethod(
  "as.geneSetList",
  signature = "geneSetFile",
  definition = function(object) {
    getGeneSet(object, geneSet = listGeneSets(object))
  }
)


#' Build a geneSetList or geneSetFile
#'
#' Build a [`geneSetList`] or [`geneSetFile`] for use in gene set analyses
#' ([`geneSetAssoc`] or [`assocTest-aggdb`]).
#' Currently these can be built directly from GMT-files, data.frames and lists.
#'
#' @param data Can be:
#' \itemize{
#'   \item A `data.frame` where the first column includes the names of geneSets,
#'   the second column contains the genes included in the geneSet (comma-delimited).
#'   The third and fourth column are optional: the third column contains weights (comma-delimited),
#'   and the fourth column contains metadata.
#'   \item A `list`, where the names of the list represent the geneSet names and each element in the list contains
#'   a vector with the genes included in the respective geneset.
#' }
#' @param gmtpath Path to a gmt-file.
#' @param output Optional output file path (output will be gz compressed text).
#' Defaults to `NULL`, in which case a [`geneSetList`] is returned.
#' @param sep Separator used in input file. Defaults to `\t`.
#' @param verbose Should the function be verbose? (TRUE/FALSE). Defaults to `TRUE`.
#'
#' @example inst/examples/example-geneSet.R
#'
#' @export
buildGeneSet <- function(
  data = NULL,
  gmtpath = NULL,
  output = NULL,
  sep = "\t",
  verbose = TRUE
) {
  # input validation
  .buildGeneSet_validate_input(as.list(environment()))

  # build from data.frame, list or gmtpath
  if (!is.null(data)) {
    if (is.data.frame(data)) {
      genesets <- .buildGeneSet_from_df(data)
    } else {
      genesets <- .buildGeneSet_from_list(data)
    }
  } else {
    genesets <- sort(readGMT(gmtpath, sep = sep))
  }

  # set metadata
  genesets@metadata <- list(
    rvatVersion = as.character(packageVersion("rvat")),
    source = if (!is.null(gmtpath)) gmtpath else "interactive_session",
    creationDate = as.character(round(Sys.time(), units = "secs"))
  )

  # handle output
  if (!is.null(output)) {
    out <- gzfile(output, "w")

    .write_rvat_header(
      filetype = "geneSetFile",
      metadata = genesets@metadata,
      con = out
    )
    write.table(
      as.data.frame(genesets),
      out,
      sep = "|",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      append = TRUE
    )
    if (verbose) {
      message(sprintf("Generated geneSetFile: %s", output))
    }

    # return geneSetFile connection
    close(out)
    geneSetFile(output)
  } else {
    # if output is not specified, return geneSetList
    genesets
  }
}

.buildGeneSet_validate_input <- function(args) {
  # data sould be at most a four column data.frame/tibble or a list
  if (
    !is.data.frame(args[["data"]]) &&
      !is.list(args[["data"]]) &&
      !is.null(args[["data"]])
  ) {
    stop(
      "`data` should be a data.frame, list, or NULL.",
      call. = FALSE
    )
  }

  # only one of `data` or `gmtpath` should be specified
  if (!is.null(args[["data"]]) && !is.null(args[["gmtpath"]])) {
    stop(
      "Specify only one of `data` or `gmtpath`.",
      call. = FALSE
    )
  }
  #  one of `data` or `gmtpath` should be specified
  if (is.null(args[["data"]]) && is.null(args[["gmtpath"]])) {
    stop(
      "Specify one of `data` or `gmtpath`.",
      call. = FALSE
    )
  }

  # check gmtpath, if specified
  if (!is.null(args[["gmtpath"]])) {
    .check_input_path(args[["gmtpath"]])
  }
  # check output, if specified
  if (!is.null(args[["output"]])) {
    .check_output(
      args[["output"]],
      overWrite = TRUE,
      verbose = args[["verbose"]]
    )
  }

  check_wrapper(check_character, args, "sep", length_equal = 1L)
  check_wrapper(check_bool, args, "verbose", length_equal = 1L)

  invisible(NULL)
}

.buildGeneSet_from_df <- function(data) {
  # validate columns
  ncols <- ncol(data)
  if (ncols > 4L) {
    stop(
      "data.frame should have at most 4 columns (geneSet names, units, weights, metadata)",
      call. = FALSE
    )
  }
  if (ncols < 2L) {
    stop(
      "data.frame should have at least 2 columns (geneSet names and units)",
      call. = FALSE
    )
  }

  # standardize to 4 columns
  if (ncols == 2L) {
    data[, c(3, 4)] <- NA_character_
  }
  if (ncols == 3L) {
    data[, 4] <- NA_character_
  }

  # build list of genesets
  genesets <- lapply(
    seq_len(nrow(data)),
    FUN = function(i, data) {
      geneSet(
        geneSetName = as.character(data[i, 1]),
        units = as.character(data[i, 2]),
        w = as.character(data[i, 3]),
        metadata = as.character(data[i, 4])
      )
    },
    data = data
  )

  # convert list to geneSetList and return
  geneSetList(genesets)
}

.buildGeneSet_from_list <- function(data) {
  # check if each element is a character vector
  if (!all(vapply(data, is.character, logical(1)))) {
    stop(
      "Each element in the list should be a character vector.",
      call. = FALSE
    )
  }
  # build directly from list, should be all character vectors
  sets <- lapply(seq_along(data), function(i) {
    geneSet(
      geneSetName = names(data)[i],
      units = paste(data[[i]], collapse = ","),
      w = NA_character_,
      metadata = NA_character_
    )
  })

  # return geneSetList
  geneSetList(sets)
}


# readers ----------------------------------------------------------------------

#' Read a gmt-file
#'
#' Read a gmt-file and return an object of class [`geneSetList`].
#' @param path File path
#' @param sep Delimiter used
#' @export
readGMT <- function(path, sep = "\t") {
  # check if file exists
  if (!file.exists(path)) {
    stop(sprintf("GMT file '%s' does not exist.", path), call. = FALSE)
  }

  # read & split
  gmt <- readLines(path)
  gmt <- strsplit(gmt, sep, fixed = TRUE)

  # extract metadata
  metadata <- vapply(
    gmt,
    FUN = function(x) x[2],
    character(1L)
  )

  # extract gene sets
  sets <- lapply(gmt, FUN = function(x) {
    geneSet(
      geneSetName = x[1],
      units = paste(x[3:length(x)], collapse = ","),
      w = NULL,
      metadata = x[2]
    )
  })

  # return geneSetList
  geneSetList(sets)
}
