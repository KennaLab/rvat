# metadata ---------------------------------------------------------------------
.parse_rvat_header <- function(
  path,
  expected_metadata = NULL,
  expected_filetype = NULL,
  return_skip = FALSE,
  n = 10
) {
  lines <- readLines(path, n = n)

  # select commented lines
  header <- gsub("#\\s*", "", lines[grepl("^#", lines)])
  skip_count <- length(header)

  # if header is present, check filetype and parse metadata lines
  filetype <- NULL
  if (length(header) > 0) {
    # check if filetype is correctly specified
    if (!is.null(expected_filetype)) {
      filetype <- header[1]
      if (!filetype %in% sprintf("RVAT-%s", expected_filetype)) {
        stop(sprintf("Unexpected filetype: %s", filetype))
      }
    }

    # parse metadata
    metadata <- header[2:length(header)]
    metadata <- strsplit(metadata, ": ", fixed = TRUE)
    if (!all(lengths(metadata) == 2)) {
      stop("Unexpected metadata")
    }
    metadata <- as.list(setNames(
      unlist(lapply(metadata, "[[", 2)),
      unlist(lapply(metadata, "[[", 1))
    ))
  } else {
    metadata <- list()
  }

  # if expected metadata is specified, check if all metadata is found
  if (!is.null(expected_metadata)) {
    # check if other than expected metadata is included
    if (!all(names(metadata) %in% expected_metadata)) {
      stop(sprintf(
        "The following unexpected metadata fields were found: %s",
        paste(
          names(metadata)[!names(metadata) %in% expected_metadata],
          collapse = ","
        )
      ))
    }

    # check if metadata is missing
    if (!all(expected_metadata %in% names(metadata))) {
      warning(sprintf(
        "The following metadata fields are missing: %s",
        paste(
          expected_metadata[!expected_metadata %in% names(metadata)],
          collapse = ","
        )
      ))
    }

    # parse
    metadata_ <- list()
    for (item in expected_metadata) {
      metadata_[[item]] <- if (item %in% names(metadata)) {
        metadata[[item]]
      } else {
        NA_character_
      }
    }
    metadata <- metadata_
  }

  # return metadata
  if (return_skip) {
    list(
      metadata = metadata,
      filetype = gsub("^RVAT-", "", filetype),
      skip = skip_count
    )
  } else {
    metadata
  }
}


# write metadata
.write_rvat_filetype <- function(filetype, con) {
  writeLines(sprintf("# RVAT-%s", filetype), con = con)
}

.write_metadata <- function(metadata, con) {
  for (item in names(metadata)) {
    writeLines(sprintf("# %s: %s", item, metadata[[item]]), con = con)
  }
}

.write_rvat_header <- function(filetype, metadata, con) {
  .write_rvat_filetype(filetype, con)
  .write_metadata(metadata, con)
}

# input checks
.check_gdb_ids <- function(gdb, object, minVersion = NULL) {
  # check if gdb id in object and gdb match
  if (
    !is.null(getGdbId(object)) &&
      !is.na(getGdbId(object)) &&
      !is.null(getGdbId(gdb)) &&
      !is.na(getGdbId(gdb))
  ) {
    if (getGdbId(gdb) != getGdbId(object)) {
      stop(sprintf(
        "The %s seems to be generated from a different gdb than supplied. Please check using `getGdbId`. Set `strict` = FALSE to ignore.",
        as.character(class(object))
      ))
    }
  }

  if (is.null(minVersion)) {
    if (
      !is.null(getRvatVersion(object)) &&
        !is.na(getRvatVersion(object)) &&
        !is.null(getRvatVersion(gdb)) &&
        !is.na(getRvatVersion(gdb))
    ) {
      if (getRvatVersion(gdb) != getRvatVersion(object)) {
        warning(sprintf(
          "The gdb and %s were generated using different RVAT versions",
          as.character(class(object))
        ))
      }
    }
  } else {
    version_gdb <- getRvatVersion(gdb)
    version_object <- getRvatVersion(object)
    if (
      !is.null(version_object) &&
        !is.na(version_object) &&
        !is.null(version_gdb) &&
        !is.na(version_gdb)
    ) {
      if (
        (version_gdb < minVersion && version_object >= minVersion) ||
          (version_object < minVersion && version_gdb >= minVersion)
      ) {
        warning(sprintf(
          "The gdb was generated with RVAT version %s, while the %s was generated with RVAT version %s",
          version_gdb,
          as.character(class(object)),
          version_object
        ))
      }
    }
  }
}

.convert_chrom_to_num <- function(CHROM) {
  # return empty vector if CHROM is null or empty
  if (is.null(CHROM) || length(CHROM) == 0L) {
    return(integer(0L))
  }

  # map
  map <- c(X = 23L, Y = 24L, MT = 26L, M = 26L)
  CHROM <- toupper(gsub("^chr", "", CHROM, ignore.case = TRUE))
  idx <- CHROM %in% names(map)
  CHROM[idx] <- unname(map[CHROM[idx]])

  # convert to integer
  CHROM <- as.integer(CHROM)
  if (anyNA(CHROM)) {
    warning("Some chromosomes could not be converted.", call. = FALSE)
  }

  # return
  CHROM
}

.convert_num_to_chrom <- function(CHROM) {
  if (is.null(CHROM) || length(CHROM) == 0L) {
    return(character(0L))
  }
  map <- c(`23` = "X", `24` = "Y", `26` = "MT")
  CHROM <- as.character(CHROM)
  idx <- CHROM %in% names(map)
  CHROM[idx] <- unname(map[CHROM[idx]])
  CHROM
}

.check_output <- function(output, overWrite, verbose) {
  # output should be one filepath
  if (!is.character(output) || length(output) != 1L || !nzchar(output)) {
    stop("`output` must be a single filepath", call. = FALSE)
  }

  # check if output already exists, if so, overwrite if `overWrite = TRUE`
  if (file.exists(output)) {
    if (overWrite) {
      removed <- file.remove(output)
      if (!removed) {
        stop(
          sprintf(
            paste0(
              "Failed to remove existing output file '%s'. ",
              "Please check permissions or ensure it is not a directory."
            ),
            output
          ),
          call. = FALSE
        )
      } else if (verbose) {
        message(sprintf(
          paste0(
            "Output file '%s' already exists and ",
            "is overwritten (`overWrite = TRUE`)"
          ),
          output
        ))
      }
    } else {
      stop(
        sprintf(
          paste0(
            "Output file '%s' already exists. ",
            "Set `overWrite=TRUE` to overwrite existing files."
          ),
          output
        ),
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}

.check_input_path <- function(path) {
  # path should be a single filepath
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a single filepath", call. = FALSE)
  }

  # check if path exists
  if (!file.exists(path)) {
    stop(sprintf("The file '%s' does not exist.", path), call. = FALSE)
  }

  # check if path is a file (not a directory)
  if (dir.exists(path)) {
    stop(sprintf("'%s' is a directory, not a file.", path), call. = FALSE)
  }

  # check if file is readable
  if (file.access(path, mode = 4) != 0) {
    stop(sprintf("The file '%s' is not readable.", path), call. = FALSE)
  }

  invisible(NULL)
}
