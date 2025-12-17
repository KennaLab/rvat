#' @include rvatResult.R
#' @include gsaResult.R

# geneSet --------------------------------------------------------

#' @rdname geneSet
#' @usage NULL
#' @export
geneSet = function(geneSetName, units, w = NULL, metadata = "") {
  if (is.null(w) || is.na(w)) {
    w <- paste(
      rep(1, length(unlist(strsplit(units, split = ",")))),
      collapse = ","
    )
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
  unlist(strsplit(object@units, split = ","))
})

#' @rdname geneSet
#' @usage NULL
#' @export
setMethod("listWeights", signature = "geneSet", definition = function(object) {
  x <- unlist(strsplit(object@w, split = ","))
  names(x) <- listUnits(object)
  x
})


# geneSetList ------------------------------------------------------------------

#' geneSetList
#' @rdname geneSetList
#' @usage NULL
#' @export
geneSetList <- function(geneSets, metadata = list()) {
  geneSetNames <- unlist(lapply(geneSets, function(x) {
    x@geneSetName
  }))
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
      warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
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
setMethod("sort", c("geneSetList"), function(x) {
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
  do.call(rbind, lapply(x@geneSets, FUN = as.data.frame))
})


#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("as.list", signature = "geneSetList", definition = function(x) {
  names <- names(x)
  x <- strsplit(
    unlist(lapply(1:length(x), function(i) {
      x[[i]]@units
    })),
    split = ","
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
    if (sparse) return(as(dat, "sparseMatrix")) else return(dat)
  }
)


#' @rdname getGeneSet
#' @usage NULL
#' @export
setMethod(
  "getGeneSet",
  signature = "geneSetList",
  definition = function(object, geneSet = NULL, unit = NULL) {
    if (is.null(geneSet) && is.null(unit)) {
      stop("At least one of `geneSet` or `unit` should be specified.")
    }

    if (!is.null(unit)) {
      object <- object[unlist(lapply(object@geneSets, FUN = function(x) {
        any(unit %in% listUnits(x))
      }))]
    }
    if (!is.null(geneSet)) {
      if (!all(geneSet %in% names(object))) {
        warning(
          "Not all specified geneSets are present in the geneSetList, use `listGeneSets()` to check which geneSets are available."
        )
      }
      object <- object[names(object) %in% geneSet]
    }
    return(object)
  }
)


#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod(
  "dropUnits",
  signature = "geneSetList",
  definition = function(object, unit = NULL) {
    object@geneSets <- lapply(object@geneSets, function(x) {
      x@units <- paste(listUnits(x)[!listUnits(x) %in% unit], collapse = ",")
      x@w <- paste(listWeights(x)[!listUnits(x) %in% unit], collapse = ",")
      return(x)
    })
    object <- object[lengths(object) > 0]
    return(object)
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
    duplicate_ids <- match.arg(duplicate_ids)

    # Check validity of dictionary
    if (ncol(dict) != 2) {
      stop("`dict` should be a data.frame with two columns")
    }
    colnames(dict) <- c("original_id", "new_id")
    dict$original_id <- as.character(dict$original_id)
    dict$new_id <- as.character(dict$new_id)

    # Number of IDs in geneSetList that are present in dictionary
    original_ids <- listUnits(object)
    if (verbose) {
      message(sprintf(
        "%s/%s IDs in the geneSetList are present in the linker file.",
        sum(original_ids %in% dict$original_id),
        length(original_ids)
      ))
    }
    dict <- dict[dict$original_id %in% original_ids, , drop = FALSE]

    # Handle IDs that map to multiple IDs
    if (!is.null(targets)) {
      dict$target <- ifelse(dict$new_id %in% targets, TRUE, FALSE)
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

    lengths_premapping <- lengths(object)

    # remap the IDs in each geneset based on the dictionary
    remapped <- lapply(
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
    object@geneSets <- remapped
    lengths_postmapping <- lengths(object)

    # return remapped object
    object
  }
)

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("write", "geneSetList", function(x, file = "data", append = FALSE) {
  out <- gzfile(file, "w")
  metadata <- metadata(x)
  metadata$rvatVersion <- as.character(packageVersion("rvat"))
  metadata$creationDate <- as.character(round(Sys.time(), units = "secs"))
  .write_rvat_header(filetype = "geneSetFile", metadata = metadata, con = out)
  write.table(
    as.data.frame(x),
    out,
    sep = "|",
    append = append,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  close(out)
})

# geneSetFile ------------------------------------------------------------------

#' geneSetFile
#'
#' Initialize \code{\link{geneSetFile-class}} object
#' @param path Path to geneSet file.
#' @param memlimit Maximum number of records to load at a time.
#' @export
geneSetFile = function(path, memlimit = 5000) {
  # read in metadata
  metadata <- .parse_rvat_header(
    path,
    expected_metadata = metadata_genesets,
    expected_filetype = "geneSetFile",
    n = length(metadata_genesets) + 1 # file description + metadata
  )
  header <- readLines(path, n = length(metadata_genesets) + 10)
  skip <- sum(startsWith(header, "#"))

  con = gzfile(path, "r")
  sets = c()
  counter <- 1
  while (length(i <- readLines(con, n = memlimit)) > 0) {
    if (counter == 1) {
      i <- i[(skip + 1):length(i)]
    }
    sets = c(sets, sapply(strsplit(i, split = "\\|"), "[[", 1))
    counter <- counter + 1
  }
  close(con)
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
    if (!all(geneSet %in% object@sets)) {
      warning(
        "Not all specified geneSets are present in the geneSetFile, use `listGeneSets()` to check which geneSets are available."
      )
    }

    # metadata
    header <- readLines(object@path, n = length(metadata_genesets) + 10)
    skip <- sum(startsWith(header, "#"))
    if (skip > length(metadata_genesets) + 1) {
      # metadata + filetype
      stop("File contains more header lines than expected.")
    }

    indices <- sort(which(object@sets %in% geneSet))
    set <- object@sets[indices]
    if (length(indices) > 1) {
      indices[2:length(indices)] <- (dplyr::lead(indices) - indices)[
        1:(length(indices) - 1)
      ]
    }
    con <- gzfile(object@path, "r")
    if (skip > 0) {
      skip <- readLines(con, n = skip)
    }
    x <- lapply(indices, FUN = function(i) {
      i <- scan(con, skip = i - 1, nlines = 1, what = "character", quiet = TRUE)
      i <- unlist(strsplit(i, split = "\\|"))
      geneSet(
        geneSetName = i[1],
        units = i[2],
        w = i[3],
        metadata = i[4]
      )
    })
    x <- geneSetList(x, metadata = metadata(object))
    close(con)

    # check
    if (!all(listGeneSets(x) %in% geneSet)) {
      stop("Something's wrong..")
    }
    x
  }
)


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
#' Build a [`geneSetList`] or [`geneSetFile`] for use in gene set analyses ([`geneSetAssoc`] or [`assocTest-aggdb`])
#' Currently these can be build directly from GMT-files, data.frames and lists.
#'
#' @param data Can be 1) data.frame where the first column includes the names of geneSets,
#' the second column contains the genes included in the geneSet (comma-delimited).
#' The third and fourth column are optional: the third column contains weights (comma-delimited),
#' and the fourth column contains metadata.
#' or 2) a list, where the names of the list represent the geneSet names and each element in the list contains
#' a vector with the genes included in the respective geneset.
#' @param gmtpath Path to a gmt-file
#' @param output Optional output file path (output will be gz compressed text).
#' Defaults to `NULL`, in which case a [`geneSetList`] is returned.
#' @param sep Separator used in input file.
#' @param verbose Should the function be verbose? (TRUE/FALSE), defaults to `TRUE`.
#'
#' @examples
#'
#' # build a genesetlist from a list (see ?geneSetList)
#' genesetlist <- buildGeneSet(
#'   list("geneset1" = c("SOD1", "NEK1"),
#'        "geneset2" = c("ABCA4", "SOD1", "NEK1"),
#'        "geneset3" = c("FUS", "NEK1")
#'        ))
#'
#' # specify the output parameter to write to disk in the geneSetFile format (see ?geneSetFile)
#' file <- tempfile()
#' buildGeneSet(
#'   list("geneset1" = c("SOD1", "NEK1"),
#'        "geneset2" = c("ABCA4", "SOD1", "NEK1"),
#'        "geneset3" = c("FUS", "NEK1")
#'   ),
#'   output = file
#'   )
#' genesetfile <- geneSetFile(file)
#'
#' # the `gmtpath` parameter can be used to build a geneset from a mSigDb GMT-file
#' # see the tutorials on the RVAT website for examples
#'
#' @export
buildGeneSet <- function(
  data = NULL,
  gmtpath = NULL,
  output = NULL,
  sep = "\t",
  verbose = TRUE
) {
  if (!is.null(data)) {
    if (is.data.frame(data)) {
      ncol <- ncol(data)
      if (ncol > 4) {
        stop(
          "data.frame should have at most 4 columns (geneSet names, units, weights, metadata)"
        )
      }
      if (ncol == 3) {
        data[, 4] <- NA
      }
      if (ncol == 2) {
        data[, c(3, 4)] <- NA
      }

      genesets <- lapply(
        1:nrow(data),
        FUN = function(i, data) {
          geneSet(
            geneSetName = data[i, 1],
            units = data[i, 2],
            w = data[i, 3],
            metadata = data[i, 4]
          )
        },
        data = data
      )
      genesets <- geneSetList(genesets)
    } else if (is.list(data)) {
      classes <- unlist(lapply(data, FUN = class))
      if (all(classes == "character")) {
        genesets <- lapply(
          1:length(data),
          FUN = function(i) {
            geneSet(
              geneSetName = names(data)[i],
              units = paste(data[[i]], collapse = ","),
              w = NA,
              metadata = NA
            )
          }
        )
        genesets <- geneSetList(genesets)
      } else {
        stop("Each element in the list should be a character vector")
      }
    }
  }

  if (!is.null(gmtpath)) {
    genesets <- sort(readGMT(gmtpath, sep = sep))
  }
  genesets@metadata <- list(
    rvatVersion = as.character(packageVersion("rvat")),
    source = if (!is.null(gmtpath)) gmtpath else "interactive_session",
    creationDate = as.character(round(Sys.time(), units = "secs"))
  )

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
    close(out)
    if (verbose) {
      message(sprintf("Generated geneSetFile: %s", output))
    }
    return(geneSetFile(output))
  } else {
    return(genesets)
  }
}

# readers ----------------------------------------------------------------------

#' Read a gmt-file
#'
#' Read a gmt-file and return an object of class [`geneSetList`].
#' @param path File path
#' @param sep Delimiter used
#' @export
readGMT <- function(path, sep = "\t") {
  gmt <- readLines(path)
  gmt <- strsplit(gmt, sep)
  metadata <- unlist(lapply(gmt, FUN = function(x) x[2]))

  sets <- lapply(gmt, FUN = function(x) {
    geneSet(
      geneSetName = x[1],
      units = paste(x[3:length(x)], collapse = ","),
      w = paste(rep(1, times = (length(x) - 2)), collapse = ","),
      metadata = x[2]
    )
  })
  geneSetList(sets)
}

# geneSetAssoc -----------------------------------------------------------------

setMethod(
  "checkDuplicates",
  signature = c("rvbResult"),
  definition = function(object, stop = TRUE) {
    if (length(unique(object$unit)) != nrow(object)) {
      if (stop) {
        stop("Duplicated units are present in the rvbResult object")
      } else {
        warning("Duplicated units are present in the rvbResult object")
      }
    }
  }
)
