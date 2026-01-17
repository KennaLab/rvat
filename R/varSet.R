# general methods -------------------------------------------------------------

## show

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("show", signature = "varSet", definition = function(object) {
  cat(sprintf(
    "unit=%s\nvarSetName=%s\nVAR_id=%s\nw=%s\n",
    object@unit,
    object@varSetName,
    object@VAR_id,
    object@w
  ))
})

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("metadata", signature = "varSet", definition = function(x) {
  if (!is.null(x@metadata)) x@metadata else NULL
})

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("getGdbId", signature = "varSet", definition = function(object) {
  metadata(object)$gdbId
})

#' @rdname varSet
#' @usage NULL
#' @export
setMethod(
  "getRvatVersion",
  signature = "varSet",
  definition = function(object) {
    metadata(object)$rvatVersion
  }
)

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("show", signature = "varSetList", definition = function(object) {
  cat(sprintf("varSetList\nContains %s records\n", length(object)))
  print(head(object@varSets))
})

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("show", signature = "varSetFile", definition = function(object) {
  cat(sprintf(
    "varSetFile object\nPath:%s\nUnits:%s\n",
    object@path,
    length(object@units)
  ))
})

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("metadata", signature = "varSetFile", definition = function(x) {
  x@metadata
})

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("getGdbId", signature = "varSetFile", definition = function(object) {
  metadata(object)$gdbId
})

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod(
  "getRvatVersion",
  signature = "varSetFile",
  definition = function(object) {
    metadata(object)$rvatVersion
  }
)

## length
### varSetList

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("length", signature = "varSetList", definition = function(x) {
  length(x@varSets)
})

### varSetFile
#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("length", signature = "varSetFile", definition = function(x) {
  length(listUnits(x))
})

## subsetting

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("[[", c("varSetList", "ANY", "missing"), function(x, i, j, ...) {
  y <- x@varSets[[i, ...]]
  y@metadata <- metadata(x)
  y
})

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod(
  "[",
  c("varSetList", "ANY", "ANY"),
  function(x, i, j, ..., drop = TRUE) {
    if (1L != length(drop) || (!missing(drop) && drop)) {
      warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'", call. = FALSE)
    }

    if (missing(i) && missing(j)) {
      return(x)
    }

    initialize(x, varSets = x@varSets[i, ...], units = x@units[i, ...])
  }
)

## listing items

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("listUnits", signature = "varSetFile", definition = function(object) {
  object@units
})

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("listUnits", signature = "varSetList", definition = function(object) {
  object@units
})

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod(
  "listVarSets",
  signature = "varSetFile",
  definition = function(object, memlimit = 5000L) {
    # input validation
    check_number_whole(memlimit, min = 1)

    # determine number of header lines to skip
    skip <- .varsetfile_read_metadata(object, return_header = FALSE)

    # establish connection
    con <- gzfile(object@path, "rb")
    on.exit(close(con), add = TRUE)

    # skip header if present
    if (skip > 0L) {
      readLines(con, n = skip)
    }

    # read and parse
    varsets <- list()
    while (length(chunk <- readLines(con, n = memlimit)) > 0L) {
      varsets[[length(varsets) + 1L]] <-
        stringi::stri_split_fixed(
          chunk,
          pattern = "|",
          n = 4L,
          simplify = TRUE
        )[, 4L]
    }
    varsets <- unlist(varsets, use.names = FALSE)

    # throw error if number of varsets doesn't equal length
    if (length(varsets) != length(object)) {
      stop(
        sprintf(
          paste0(
            "Error reading varSetFile: The number of data lines ",
            "read (%d) does not match the number of units cached ",
            "in the object (%d). File may be corrupt or inconsistent."
          ),
          length(varsets),
          length(object)
        ),
        call. = FALSE
      )
    }

    # return varsets
    varsets
  }
)

.varsetfile_read_metadata <- function(object, return_header = FALSE) {
  header <- readLines(object@path, n = length(metadata_varsets) + 10L)
  skip <- sum(startsWith(header, "#"))
  if (skip > length(metadata_varsets) + 1L) {
    stop("File contains more header lines than expected.", call. = FALSE)
  }

  if (return_header) {
    invisible(NULL)
  } else {
    skip
  }
}

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod(
  "listVarSets",
  signature = "varSetList",
  definition = function(object) {
    unlist(lapply(object@varSets, FUN = function(x) x@varSetName))
  }
)

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("listVars", signature = "varSet", definition = function(object) {
  unlist(strsplit(object@VAR_id, split = ",", fixed = TRUE))
})

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("listWeights", signature = "varSet", definition = function(object) {
  weights <- unlist(strsplit(object@w, split = ",", fixed = TRUE))
  weights <- ifelse(weights == ".", NA_character_, weights)
  as.numeric(weights)
})

## extract metadata
#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("getGdbId", signature = "varSetList", definition = function(object) {
  metadata(object)$gdbId
})

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod(
  "getRvatVersion",
  signature = "varSetList",
  definition = function(object) {
    metadata(object)$rvatVersion
  }
)


# Constructors -----------------------------------------------------------------

# Initialize a varSet object
varSet <- function(unit, varSetName, VAR_id, w) {
  new("varSet", unit = unit, varSetName = varSetName, VAR_id = VAR_id, w = w)
}


#' varSetList
#' @rdname varSetList
#' @usage NULL
#' @export
varSetList <- function(x, metadata = list()) {
  # return empty varSetList if called without arguments
  if (missing(x)) {
    return(new(
      "varSetList",
      varSets = list(),
      units = character(),
      metadata = list()
    ))
  }

  # x can be either a list of varSets or
  # a data.frame with fields unit,varSetName,VAR_id,w
  if (is(x, "list")) {
    units_vec <- unlist(lapply(x, function(x) {
      x@unit
    }))
    varsetlist <- new(
      "varSetList",
      varSets = x,
      units = units_vec,
      metadata = metadata
    )
  } else if (is(x, "data.frame")) {
    if (!all(c("unit", "varSetName", "VAR_id", "w") %in% colnames(x))) {
      stop(
        "Invalid input. Input must be either an object of class list that ",
        "contains one or more varSet records, or else a data frame ",
        "with fields unit, varSetName, VAR_id, w",
        call. = FALSE
      )
    }
    out <- vector(mode = "list", length = nrow(x))
    for (i in seq_len(nrow(x))) {
      out[[i]] <- varSet(
        unit = x$unit[i],
        varSetName = x$varSetName[i],
        VAR_id = x$VAR_id[i],
        w = x$w[i]
      )
    }
    varsetlist <- new(
      "varSetList",
      varSets = out,
      units = x$unit,
      metadata = metadata
    )
  } else {
    stop(
      "Invalid input. Input must be either an object of class list that ",
      "contains one or more varSet records, or else a data frame with ",
      "fields unit, varSetName, VAR_id, w",
      call. = FALSE
    )
  }

  # return varsetlist
  varsetlist
}


#' @importFrom S4Vectors setValidity2
S4Vectors::setValidity2("varSetList", function(object) {
  if (length(object) == 0L) {
    return(TRUE)
  }

  # ensure all records are of class varSet
  is_varset <-
    vapply(
      object@varSets,
      is,
      logical(1L),
      "varSet",
      USE.NAMES = FALSE
    )

  if (all(is_varset)) {
    TRUE
  } else {
    "Invalid varSetList: All elements in list must be valid 'varSet' objects."
  }
})


#' @rdname varSetFile
#' @usage NULL
#' @export
varSetFile <- function(path, memlimit = 5000L) {
  # input validation
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a single, non-empty file path string.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("File does not exist '%s'", path), call. = FALSE)
  }
  check_number_whole(memlimit, min = 1)

  # read header
  header_info <- .parse_rvat_header(
    path = path,
    expected_metadata = metadata_varsets,
    expected_filetype = "varSetFile",
    return_skip = TRUE,
    n = length(metadata_varsets) + 1L # file description + metadata
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

  # read and parse units
  units <- list()
  while (length(chunk <- readLines(con, n = memlimit)) > 0L) {
    units[[length(units) + 1L]] <-
      stringi::stri_split_fixed(
        chunk,
        pattern = "|",
        n = 4L,
        simplify = TRUE
      )[, 1L]
  }
  units <- unlist(units, use.names = FALSE)

  # return varSetFile
  new("varSetFile", path = path, units = units, metadata = metadata)
}


# Getters ----------------------------------------------------------------------
#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod(
  "getVarSet",
  signature = "varSetFile",
  definition = function(object, unit, varSetName = NULL) {
    # input validation
    check_character(unit)
    check_character(varSetName, allow_null = TRUE)

    # check if all specified units are present in varSetFile
    if (!all(unit %in% object@units)) {
      warning(
        "Not all specified units are present in the varSetFile, ",
        "use `listUnits()` to check which units are available.",
        call. = FALSE
      )
    }

    # determine indices to scan
    indices <- sort(which(object@units %in% unit))

    # calculate relative skips for scanning
    relative_indices <- indices
    if (length(indices) > 1L) {
      relative_indices[2L:length(indices)] <- (dplyr::lead(indices) - indices)[
        1L:(length(indices) - 1L)
      ]
    } else {
      relative_indices <- indices
    }

    # check header and return number of lines to skip
    skip <- .varsetfile_read_metadata(object)

    # establish connection
    con <- gzfile(object@path, "rb")
    on.exit(close(con), add = TRUE)

    # skip header if present
    if (skip > 0L) {
      readLines(con, n = skip)
    }

    # retrieve varsets
    varsets <- lapply(relative_indices, FUN = function(i) {
      line_content <- scan(
        con,
        skip = i - 1L,
        nlines = 1L,
        what = "character",
        quiet = TRUE
      )
      line_content <- unlist(strsplit(line_content, split = "|", fixed = TRUE))
      if (length(line_content) != 4L) {
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
      varSet(
        unit = line_content[1L],
        varSet = line_content[4L],
        VAR_id = line_content[2L],
        w = line_content[3L]
      )
    })
    varsets <- varSetList(varsets)

    # throw error if unrequested units are among varSet
    if (!all(listUnits(varsets) %in% unit)) {
      stop(
        "Varsets retrieved for unspecified units, something's wrong.",
        call. = FALSE
      )
    }

    # additionally subset varSets, if specified
    if (!is.null(varSetName)) {
      varsets <- getVarSet(varsets, varSetName = varSetName)
    }

    # inherit metadata from varSetFile
    varsets@metadata <- metadata(object)

    # return
    varsets
  }
)

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod(
  "getVarSet",
  signature = "varSetList",
  definition = function(object, unit = NULL, varSetName = NULL) {
    # input validation
    if (is.null(unit) && is.null(varSetName)) {
      stop(
        "At least one of `unit` or `varSetName` should be specified.",
        call. = FALSE
      )
    }
    check_character(unit, allow_null = TRUE, allow_na = FALSE)
    check_character(varSetName, allow_null = TRUE, allow_na = FALSE)

    # retrieve units
    if (!is.null(unit)) {
      # warning if not all specified units are present in varSetList
      if (!all(unit %in% listUnits(object))) {
        warning(
          "Not all specified units are present in the varSetList",
          call. = FALSE
        )
      }

      # subset
      object <- object[listUnits(object) %in% unit]
    }

    if (!is.null(varSetName)) {
      varsets <- listVarSets(object)
      # warning if not all specified varSets are present in varSetList
      if (!all(varSetName %in% varsets)) {
        warning(
          "Not all specified varSets are present in the varSetList",
          call. = FALSE
        )
      }

      # subset
      object <- object[varsets %in% varSetName]
    }

    # return subsetted varSetList
    object
  }
)


# write ------------------------------------------------------------------------
#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("write", "varSetList", function(x, file = "data") {
  # establish connection
  tryCatch(
    {
      con <- gzfile(file, "wb")
      on.exit(close(con), add = TRUE)
    },
    error = function(e) {
      stop(
        sprintf(
          "Failed to open file '%s' for writing. Error: %s",
          file,
          e$message
        ),
        call. = FALSE
      )
    }
  )

  # inherit metadata from varSetList, overwriting rvatVersion and creationDate
  metadata <- metadata(x)
  metadata$rvatVersion <- as.character(packageVersion("rvat"))
  metadata$creationDate <- as.character(round(Sys.time(), units = "secs"))

  # write header
  .write_rvat_header(filetype = "varSetFile", metadata = metadata, con = con)

  # write varSets
  write.table(
    as.data.frame(x),
    file = con,
    sep = "|",
    append = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  # return nothing
  invisible(NULL)
})


# varset utils -----------------------------------------------------------------

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("as.data.frame", signature = "varSetList", definition = function(x) {
  do.call(rbind, lapply(x@varSets, FUN = as.data.frame))
})

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("as.data.frame", signature = "varSet", definition = function(x) {
  data.frame(
    unit = x@unit,
    VAR_id = x@VAR_id,
    w = x@w,
    varSetName = x@varSetName,
    stringsAsFactors = FALSE
  )
})

#' buildVarSet-gdb
#'
#' Generate optionally weighted variant sets using annotation table(s) uploaded to the gdb.
#' See the tutorials for examples.
#'
#' @rdname buildVarSet-gdb
#' @name buildVarSet-gdb
#' @aliases buildVarSet,gdb-method
#' @param object a [`gdb`] object.
#' @param varSetName Name to assign varSet grouping.
#' This identifier column is used to allow for subsequent merging of multiple 
#' varSet files for coordinated analysis of multiple variant filtering/weighting strategies
#' @param unitTable Table containing aggregation unit mappings.
#' @param unitName Field to utilize for aggregation unit names.
#' @param output Output file name (output will be gz compressed text).
#' @param intersection Additional tables to filter through intersection 
#' (i.e. variants absent from intersection tables will not appear in output). 
#' Multiple tables should be ',' delimited.
#' @param where An SQL compliant where clause to filter output; 
#' e.g.: "CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')".
#' @param weightName Field name for desired variant weighting, 
#' must be a column within unitTable or other intersection table. 
#' Default value of 1 is equivalent to no weighting.
#' @param memlimit Chunk size used for processing rows.
#' @param verbose Should the function be verbose? Defaults to `TRUE`.
#'
#' @example inst/examples/example-buildVarSet-gdb.R
#'
#' @export
setMethod(
  "buildVarSet",
  signature = "gdb",
  definition = function(
    object,
    varSetName,
    unitTable,
    unitName,
    output = NULL,
    intersection = NULL,
    where = NULL,
    weightName = "1",
    memlimit = 1000L,
    verbose = TRUE
  ) {
    # validate input
    .buildVarSet_gdb_validate_input(as.list(environment()))

    query <- .buildVarSet_gdb_construct_query(
      unitTable = unitTable,
      unitName = unitName,
      weightName = weightName,
      intersection = intersection,
      where = where,
      varSetName = varSetName
    )

    # output to varSetFile
    metadata <- .gdb_list_metadata(object)

    if (!is.null(output)) {
      # perform query in chunks and write to disk
      .buildVarSet_gdb_write_varsetfile(
        gdb = object,
        query = query,
        output = output,
        metadata = metadata,
        memlimit = memlimit,
        verbose = verbose
      )
      if (verbose) {
        message(sprintf("Generated varSetFile: %s", output))
      }
      return(varSetFile(output))
    } else {
      # output to varSetList directly
      varsetlist <- varSetList(DBI::dbGetQuery(object, query))
      varsetlist@metadata <- metadata
      return(varsetlist)
    }
  }
)

.buildVarSet_gdb_validate_input <- function(args) {
  # type checks
  check_wrapper(check_character, args, "varSetName", length_equal = 1L)
  check_wrapper(check_character, args, "unitTable", length_equal = 1L)
  check_wrapper(check_character, args, "unitName", length_equal = 1L)
  check_wrapper(
    check_character,
    args,
    "output",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(
    check_character,
    args,
    "intersection",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(
    check_character,
    args,
    "where",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(
    check_character,
    args,
    "weightName",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(check_bool, args, "verbose", length_equal = 1L)

  # check if unitTable is present in gdb
  anno_tables <- c("var", listAnno(args[["object"]])$name)
  if (!args[["unitTable"]] %in% anno_tables) {
    stop(
      sprintf(
        "The specified `unitTable` ('%s') does not exist in the gdb.",
        args[["unitTable"]]
      ),
      call. = FALSE
    )
  }

  # check if unitName if present in unitTable
  fields_found <- DBI::dbListFields(args[["object"]], args[["unitTable"]])
  if (!args[["unitName"]] %in% fields_found) {
    stop(
      sprintf(
        "`unitName` ('%s') is not present in `unitTable` ('%s').",
        args[["unitName"]],
        args[["unitTable"]]
      ),
      call. = FALSE
    )
  }

  # check if intersection tables are present in gdb
  if (!is.null(args[["intersection"]])) {
    intersection_tables <- trimws(unlist(strsplit(
      args[["intersection"]],
      ",",
      fixed = TRUE
    )))
    if (!all(intersection_tables %in% anno_tables)) {
      stop(
        sprintf(
          "Not all intersection tables present in gdb (%s).",
          paste(
            paste0(
              "'",
              intersection_tables[!intersection_tables %in% anno_tables],
              "'"
            ),
            collapse = ","
          )
        ),
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}

.buildVarSet_gdb_construct_query <- function(
  unitTable,
  unitName,
  weightName,
  intersection,
  where,
  varSetName
) {
  # base query
  query <- sprintf(
    "select distinct %s as unit, VAR_id, %s as w from %s",
    unitName,
    weightName,
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
    "select unit, group_concat(VAR_id) as VAR_id, group_concat(w) as w, '%s' as varSetName from (%s) x group by unit",
    varSetName,
    query
  )

  query
}


.buildVarSet_gdb_write_varsetfile <- function(
  gdb,
  query,
  output,
  metadata,
  memlimit,
  verbose
) {
  tryCatch(
    {
      output_con <- gzfile(output, "wb")
      on.exit(close(output_con), add = TRUE)
    },
    error = function(e) {
      stop(
        sprintf(
          "Failed to open output file for writing: '%s'. Error: %s",
          output,
          e$message
        ),
        call. = FALSE
      )
    }
  )

  # write metadata to output
  .write_rvat_header(
    filetype = "varSetFile",
    metadata = metadata,
    con = output_con
  )

  # perform query in chunks (of size memlimit)
  handle <- DBI::dbSendQuery(gdb, query)
  on.exit(DBI::dbClearResult(handle), add = TRUE)
  while (!DBI::dbHasCompleted(handle)) {
    chunk <- DBI::dbFetch(handle, n = memlimit)
    write.table(
      chunk,
      output_con,
      append = TRUE,
      sep = "|",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }

  invisible(NULL)
}

#' buildVarSet-data.frame
#'
#' Generate optionally weighted variant sets using annotation table(s).
#' See the tutorials for examples.
#'
#' @rdname buildVarSet-data.frame
#' @name buildVarSet-data.frame
#' @aliases buildVarSet,data.frame-method
#' @param object A data.frame containing variant annotations
#' Must contain a 'VAR_id' column, and columns matching `unitName` and `fields`.
#' @param unitName Field to utilize for aggregation unit names.
#' @param fields Which fields in the data.frame to use as variant annotations.
#' These fields can be 1) an indicator (0,1 or FALSE/TRUE) that flags variants with the annotation,
#' e.g. a column named 'LOF' that indicates whether the variant is predicted to lead to loss-of-function;
#' or 2) variant weights, e.g. a column named 'CADD' that contains CADD scores.
#' Multiple fields can be specified, which will result in multiple rows per aggregation unit in the resulting varSetFile.
#' @param output Optional output file name (output will be gz compressed text).
#' @param memlimit Defaults to `NULL`, currently not implemented.
#' @param mask If the annotation field is a 0,1 or FALSE/TRUE indicator,
#' should variants with 0/FALSE be dropped? Defaults to `TRUE`.
#' If `FALSE`, all variants will be included, and will be assigned weights 0 and 1.
#'
#' @example inst/examples/example-buildVarSet-data.frame.R
#'
#'
#' @export
setMethod(
  "buildVarSet",
  signature = "data.frame",
  definition = function(
    object,
    unitName,
    fields,
    output = NULL,
    memlimit = NULL,
    mask = TRUE
  ) {
    # validate input
    .buildVarSet_df_validate_input(as.list(environment()))

    # subset relevant columns
    object <- object[, c("VAR_id", unitName, fields)]
    colnames(object)[2L] <- "unit"

    # test for each field whether it's 0,1 or weights
    is_mask <- .buildVarSet_df_get_masks(
      object,
      fields = fields,
      mask = mask
    )

    # generate varset
    varSet <- object %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(fields),
        names_to = "varSetName",
        values_to = "value"
      ) %>%
      dplyr::arrange(unit) %>%
      dplyr::left_join(
        is_mask,
        by = "varSetName"
      ) %>%
      dplyr::filter(
        (is_mask & !is.na(value) & value == 1) |
          (!is_mask & !is.na(value))
      ) %>%
      dplyr::group_by(dplyr::across(c("unit", "varSetName"))) %>%
      dplyr::summarize(
        VAR_id = paste(VAR_id, collapse = ","),
        w = paste(value, collapse = ","),
        .groups = "drop"
      ) %>%
      dplyr::select(unit, VAR_id, w, varSetName)

    # convert to varSetList
    metadata <- list(
      rvatVersion = as.character(packageVersion("rvat")),
      gdbId = NA_character_,
      genomeBuild = NA_character_,
      creationDate = as.character(round(Sys.time(), units = "secs"))
    )
    varsetlist <- varSetList(varSet, metadata = metadata)

    # return
    if (is.null(output)) {
      return(varsetlist)
    } else {
      write(varsetlist, file = output)
      return(invisible(NULL))
    }
  }
)

.buildVarSet_df_validate_input <- function(args) {
  # type checks
  check_wrapper(check_character, args, "unitName", length_equal = 1L)
  check_wrapper(check_character, args, "fields")
  check_wrapper(
    check_character,
    args,
    "output",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(
    check_number_whole,
    args,
    "memlimit",
    min = 1L,
    allow_null = TRUE
  )
  if (!is.logical(args[["mask"]]) && !is.character(args[["mask"]])) {
    stop(
      "`mask` must be logical (TRUE/FALSE) or a ",
      "character vector of field names.",
      call. = FALSE
    )
  }

  # check whether `VAR_id`, `unitName` and `fields` are present in data.frame
  if (!"VAR_id" %in% colnames(args[["object"]])) {
    stop(
      "'VAR_id' column must be present in the input data.frame.",
      call. = FALSE
    )
  }
  if (!args[["unitName"]] %in% colnames(args[["object"]])) {
    stop(
      sprintf(
        "The specified `unitName` column ('%s') was not found.",
        args[["unitName"]]
      ),
      call. = FALSE
    )
  }
  if (!all(args[["fields"]] %in% colnames(args[["object"]]))) {
    stop(
      sprintf(
        "The following specified `fields` were not found: %s.",
        paste(
          sQuote(
            args[["fields"]][!args[["fields"]] %in% colnames(args[["object"]])],
            "'"
          ),
          collapse = ","
        )
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

.buildVarSet_df_get_masks <- function(data, fields, mask) {
  if (is.logical(mask) && mask) {
    is_mask <- vapply(
      fields,
      function(field) {
        all(data[[field]] %in% c(0, 1), na.rm = TRUE)
      },
      logical(1L),
      USE.NAMES = TRUE
    )
    is_mask <- data.frame(
      varSetName = names(is_mask),
      is_mask = unname(is_mask),
      stringsAsFactors = FALSE
    )
  } else if (is.logical(mask) && !mask) {
    is_mask <- data.frame(
      varSetName = fields,
      is_mask = FALSE,
      stringsAsFactors = FALSE
    )
  } else if (is.character(mask)) {
    is_mask <- data.frame(
      varSetName = fields,
      is_mask = fields %in% mask,
      stringsAsFactors = FALSE
    )
  } else {
    stop(
      "`mask` should be either `TRUE`, `FALSE` or a character ",
      "vector of fields that should be treated as mask.",
      call. = FALSE
    )
  }

  is_mask
}

## build varSetlist from VAR_ids with 'memlimit' option
.varsTovarSetList <- function(VAR_id, chunkSize = Inf) {
  chunks <- split(VAR_id, ceiling(seq_along(VAR_id) / chunkSize))
  varset <- varSetList(lapply(seq_along(chunks), function(x) {
    varSet <- varSet(
      unit = paste0("chunk", x),
      varSet = "none",
      VAR_id = paste(chunks[[x]], collapse = ","),
      w = paste(rep("1", length(chunks[[x]])), collapse = ",")
    )
  }))
  varset
}

#' collapseVarSetList
#' @rdname varSetList
#' @usage NULL
#' @export
setMethod(
  "collapseVarSetList",
  signature = "varSetList",
  definition = function(
    object,
    unit = "unnamed",
    varSetName = "unnamed",
    drop = TRUE
  ) {
    VAR_id <- unique(unlist(lapply(object@varSets, listVars)))
    varset <- varSet(
      unit = unit,
      varSet = varSetName,
      VAR_id = paste(VAR_id, collapse = ","),
      w = paste(rep("1", length(VAR_id)), collapse = ",")
    )
    if (drop) {
      varset
    } else {
      varSetList(list(varset))
    }
  }
)

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("metadata", signature = "varSetList", definition = function(x) {
  x@metadata
})

#' @rdname getRanges
#' @usage NULL
#' @export
setMethod(
  "getRanges",
  signature = "varSetFile",
  definition = function(
    object,
    gdb,
    output = NULL,
    table = "var",
    CHROM = "CHROM",
    POS = "POS",
    where = NULL
  ) {
    # input validation
    .getRanges_validate_input(as.list(environment()))

    # read in varset data
    varset_df <- readr::read_delim(
      object@path,
      delim = "|",
      col_names = c("unit", "VAR_id", "w", "varSetName"),
      col_types = "cccc",
      comment = "#"
    )

    # get ranges
    .get_ranges(
      varset_df = varset_df,
      gdb = gdb,
      output = output,
      table = table,
      CHROM = CHROM,
      POS = POS,
      where = where
    )
  }
)

#' @rdname getRanges
#' @usage NULL
#' @export
setMethod(
  "getRanges",
  signature = "varSetList",
  definition = function(
    object,
    gdb,
    output = NULL,
    table = "var",
    CHROM = "CHROM",
    POS = "POS",
    where = NULL
  ) {
    # input validation
    .getRanges_validate_input(as.list(environment()))

    # convert varsetlist to data.frame
    varset_df <- as.data.frame(object)

    # get ranges
    .get_ranges(
      varset_df = varset_df,
      gdb = gdb,
      output = output,
      table = table,
      CHROM = CHROM,
      POS = POS,
      where = where
    )
  }
)

.get_ranges <- function(
  varset_df,
  gdb,
  output = NULL,
  table = "var",
  CHROM = "CHROM",
  POS = "POS",
  where = NULL
) {
  # check where statement,
  # can be either of length 1 or same length as varSetList/varSetFile
  if (
    !is.null(where) &&
      length(where) > 1L &&
      length(where) != nrow(varset_df)
  ) {
    stop(
      "The `where` parameter should be either of length 1, ",
      "or of the same length as the provided varSetList/varSetFile.",
      call. = FALSE
    )
  }

  # extract ranges for each record
  starts <- vector(mode = "integer", length = nrow(varset_df))
  ends <- vector(mode = "integer", length = nrow(varset_df))
  chroms <- vector(mode = "character", length = nrow(varset_df))
  not_found <- 0L

  for (i in seq_len(nrow(varset_df))) {
    # extract var ids
    varids <- unlist(stringr::str_split(
      varset_df$VAR_id[i],
      pattern = stringr::fixed(",")
    ))

    # retrieve ranges
    pos <- getAnno(
      gdb,
      table = table,
      VAR_id = varids,
      where = if (length(where) > 1L) where[i] else where
    )

    # raise error if multiple records are found for a single VAR_id
    if (anyDuplicated(pos$VAR_id) != 0L) {
      stop(
        sprintf(
          paste0(
            "More than one record per variant found in table '%s', ",
            "first encountered at unit '%s'"
          ),
          table,
          varset_df$unit[i]
        ),
        call. = FALSE
      )
    } else if (nrow(pos) > 0L) {
      # else retrieve ranges
      starts[i] <- min(pos[[POS]])
      ends[i] <- max(pos[[POS]])
      chroms[i] <- unique(pos[[CHROM]])
    } else {
      starts[i] <- NA_integer_
      ends[i] <- NA_integer_
      chroms[i] <- NA_character_
      not_found <- not_found + 1L
    }
  }

  # parse
  ranges <- data.frame(
    CHROM = chroms,
    start = starts,
    end = ends,
    unit = varset_df$unit,
    varSetName = varset_df$varSetName,
    stringsAsFactors = FALSE
  )

  # notify if no positions were found for records
  if (not_found > 0L) {
    warning(
      sprintf(
        paste0(
          "Note: for %s/%s records in the varSetFile, ",
          "no position mappings were found for any of the variants."
        ),
        not_found,
        nrow(varset_df)
      ),
      call. = FALSE
    )
  }

  # write output, if specified
  if (!is.null(output)) {
    readr::write_tsv(ranges, file = gzfile(output), col_names = TRUE)
    invisible(NULL)
  } else {
    ranges
  }
}

.getRanges_validate_input <- function(args) {
  check_wrapper(
    check_character,
    args,
    "output",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(check_character, args, "table", length_equal = 1L)
  check_wrapper(check_character, args, "CHROM", length_equal = 1L)
  check_wrapper(check_character, args, "POS", length_equal = 1L)
  check_wrapper(check_character, args, "where", allow_null = TRUE)

  # check if required fields are present in specified table
  table_fields <- DBI::dbListFields(
    args[["gdb"]],
    args[["table"]]
  )
  required_fields <- c("VAR_id", args[["CHROM"]], args[["POS"]])

  if (!all(required_fields %in% table_fields)) {
    stop(
      sprintf(
        "The following required fields are not present in table '%s': %s.",
        args[["table"]],
        paste(
          required_fields[!required_fields %in% table_fields],
          collapse = ","
        )
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}
