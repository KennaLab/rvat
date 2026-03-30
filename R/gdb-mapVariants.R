#' @rdname mapVariants
#' @usage NULL
#' @export
setMethod(
  "mapVariants",
  signature = "gdb",
  definition = function(
    object,
    ranges = NULL,
    gff = NULL,
    bed = NULL,
    bedCols = character(),
    fields = NULL,
    uploadName = NULL,
    output = NULL,
    sep = "\t",
    skipIndexes = FALSE,
    overWrite = FALSE,
    verbose = TRUE
  ) {
    # validate input
    arg <- as.list(environment())
    .mapVariants_validate_input(arg)

    # read ranges from ranges/gff/bed files or
    # convert ranges to GRanges if data.frame is provided
    ranges <- .mapVariants_handle_ranges(
      ranges = ranges,
      bed = bed,
      gff = gff,
      sep = sep,
      bedCols = bedCols,
      verbose = verbose
    )

    # if uploadName is specified, check if name is valid and if already exists in gdb
    if (!is.null(uploadName)) {
      # validate uploadName
      .upload_validate_name(
        object = object,
        name = uploadName,
        overWrite = overWrite,
        verbose = verbose
      )
    }

    # if output!=NULL connect
    if (!is.null(output)) {
      .check_output(output, overWrite = overWrite, verbose = verbose)
      output_con <- gzcon(file(output, open = "wb"))
      on.exit(close(output_con), add = TRUE)
    }

    # if both output and uploadName are not specified, output will be returned
    if (is.null(uploadName) && is.null(output)) {
      container <- list()
    }

    # fields to keep from ranges object
    if (!is.null(fields)) {
      ranges <- ranges[, colnames(GenomicRanges::mcols(ranges)) %in% fields]
      # check if any of the fields are among the GRanges ranges fields
      if (any(fields %in% c("seqnames", "start", "end", "width", "strand"))) {
        fields_ranges <- fields[
          fields %in% c("seqnames", "start", "end", "width", "strand")
        ]
      } else {
        fields_ranges <- NULL
      }
    } else {
      fields_ranges <- NULL
    }

    # chromosomes that overlap between gdb and ranges
    # note: convert both to NCBI format to check overlap
    chroms <- getAnno(object, "var_ranges", fields = "CHROM")$CHROM
    chroms_ncbi <- GenomicRanges::GRanges(
      seqnames = chroms,
      ranges = IRanges::IRanges(start = 1L)
    )
    GenomeInfoDb::seqlevelsStyle(chroms_ncbi) <- "NCBI"
    chroms <- setNames(
      as.character(GenomicRanges::seqnames(chroms_ncbi)),
      nm = chroms
    )
    chroms <- names(chroms[chroms %in% GenomeInfoDb::seqlevels(ranges)])

    # map per chromosome
    for (chrom in chroms) {
      if (verbose) {
        message(sprintf("Mapping chromosome %s..", chrom))
      }

      # get genomicRanges for current chromosome
      gr <- unserialize(getAnno(
        object,
        "var_ranges",
        where = sprintf("CHROM = '%s'", chrom)
      )$ranges[[1]])
      GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"

      # generate overlaps
      overlaps <- GenomicRanges::findOverlaps(gr, ranges)
      if (!is.null(fields_ranges) && length(fields_ranges) > 0L) {
        dat <- cbind(
          VAR_id = gr[S4Vectors::queryHits(overlaps)]$VAR_id,
          as.data.frame(ranges[S4Vectors::subjectHits(overlaps)])[,
            fields_ranges,
            drop = FALSE
          ],
          as.data.frame(GenomicRanges::mcols(ranges[
            S4Vectors::subjectHits(overlaps),
          ]))
        )
      } else {
        dat <- cbind(
          VAR_id = gr[S4Vectors::queryHits(overlaps)]$VAR_id,
          as.data.frame(GenomicRanges::mcols(ranges[
            S4Vectors::subjectHits(overlaps),
          ]))
        )
      }
      dat <- dplyr::arrange(dat, VAR_id)
      dat[] <- lapply(dat, function(x) if (is.factor(x)) as.character(x) else x)

      # write to output (if specified)
      if (!is.null(output)) {
        if (chrom == chroms[1]) {
          write.table(
            dat,
            file = output_con,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
          )
        } else {
          write.table(
            dat,
            file = output_con,
            append = TRUE,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t"
          )
        }
      }

      # write to gdb (if `uploadName` is specified)
      if (!is.null(uploadName)) {
        if (chrom == chroms[1]) {
          DBI::dbCreateTable(con = object, name = uploadName, fields = dat)
          DBI::dbAppendTable(con = object, name = uploadName, value = dat)
        } else {
          DBI::dbAppendTable(con = object, name = uploadName, value = dat)
        }
      }

      # store in container (if output and uploadName are both NULL)
      if (is.null(output) && is.null(uploadName)) {
        container[[chrom]] <- dat
      }
    }

    # if output is uploaded to gdb, update metadata and index
    if (!is.null(uploadName)) {
      fields <- DBI::dbListFields(object, uploadName)
      if (verbose) {
        message(sprintf(
          "%s fields detected (%s)",
          length(fields),
          paste(fields, collapse = ",")
        ))
      }

      .gdb_upload_update_metadata(
        gdb = object,
        table_name = uploadName,
        type = "anno",
        source_info = "mapVariants",
        verbose = verbose
      )

      if (!skipIndexes) {
        .uploadAnno_create_indices(
          gdb = object,
          name = uploadName,
          verbose = verbose
        )
      }
    }

    # return if both output and uploadName are not specified
    if (is.null(output) && is.null(uploadName)) {
      container <- do.call(rbind, container)
      rownames(container) <- NULL
      return(container)
    } else {
      return(invisible(NULL))
    }
  }
)

.mapVariants_validate_input <- function(args) {
  # inputs: only one of ranges,gff,bed can be specified
  num_inputs <- sum(
    !is.null(args[["ranges"]]),
    !is.null(args[["gff"]]),
    !is.null(args[["bed"]])
  )
  if (num_inputs == 0L) {
    stop(
      "At least one of `ranges`, `gff` or `bed` should be specified.",
      call. = FALSE
    )
  }
  if (num_inputs > 1L) {
    stop(
      "Multiple region inputs specified, ",
      "specify either `ranges`, `gff`, or `bed`.",
      call. = FALSE
    )
  }

  # ranges can be either a filepath or a data.frame/GRanges object
  if (!is.null(args[["ranges"]])) {
    if (is.character(args[["ranges"]])) {
      check_wrapper(check_character, args, "ranges", length_equal = 1L)
      if (!file.exists(args[["ranges"]])) {
        stop(
          sprintf(
            "File specified for `ranges` does not exist: '%s'",
            args[["ranges"]]
          ),
          call. = FALSE
        )
      }
    } else if (
      !is.data.frame(args[["ranges"]]) && !is(args[["ranges"]], "GRanges")
    ) {
      stop(
        "`ranges` must be a GRanges object, a data.frame, ",
        "or a file path (character string).",
        call. = FALSE
      )
    }
  }

  # gff should be a filepath (or NULL)
  check_wrapper(
    check_character,
    args,
    "gff",
    length_equal = 1L,
    allow_null = TRUE
  )
  if (!is.null(args[["gff"]]) && !file.exists(args[["gff"]])) {
    stop(
      sprintf(
        "File specified for `gff` does not exist: '%s'",
        args[["gff"]]
      ),
      call. = FALSE
    )
  }

  # bed should be a filepath (or NULL)
  check_wrapper(
    check_character,
    args,
    "bed",
    length_equal = 1L,
    allow_null = TRUE
  )
  if (!is.null(args[["bed"]]) && !file.exists(args[["bed"]])) {
    stop(
      sprintf(
        "File specified for `bed` does not exist: '%s'",
        args[["bed"]]
      ),
      call. = FALSE
    )
  }

  check_wrapper(check_character, args, "bedCols")
  check_wrapper(check_character, args, "fields", allow_null = TRUE)

  check_wrapper(
    check_character,
    args,
    "uploadName",
    length_equal = 1L,
    allow_null = TRUE
  )

  check_wrapper(
    check_character,
    args,
    "output",
    length_equal = 1L,
    allow_null = TRUE
  )
  check_wrapper(check_character, args, "sep", length_equal = 1L)
  check_wrapper(check_bool, args, "skipIndexes", length_equal = 1L)
  check_wrapper(check_bool, args, "overWrite", length_equal = 1L)
  check_wrapper(check_bool, args, "verbose", length_equal = 1L)

  invisible(NULL)
}

.mapVariants_handle_ranges <- function(
  ranges,
  bed,
  gff,
  sep,
  bedCols,
  verbose
) {
  # read ranges if specified as filepath, or convert to GRanges if data.frame
  if (!is.null(ranges)) {
    if (is.character(ranges)) {
      ranges <- read.table(
        ranges,
        sep = sep,
        stringsAsFactors = FALSE,
        header = TRUE
      )
    }
    if (is.data.frame(ranges)) {
      if (
        !all(c("CHROM", "start", "end") %in% colnames(ranges)) &&
          !all(c("seqnames", "start", "end") %in% colnames(ranges))
      ) {
        stop(
          "either 'CHROM', 'start' and 'end' or 'seqnames', ",
          "'start' and 'end' should be present in ranges-file.",
          call. = FALSE
        )
      }
      gr <- GenomicRanges::makeGRangesFromDataFrame(
        ranges,
        keep.extra.columns = TRUE
      )
    } else if (is(ranges, "GRanges")) {
      gr <- ranges
    }
  } else if (!is.null(bed)) {
    # import bed-file using rtracklayer
    if (length(bedCols) > 0L && is.null(names(bedCols))) {
      bedCols <- setNames(
        object = rep("character", length(bedCols)),
        nm = bedCols
      )
    }
    gr <- rtracklayer::import(bed, extraCols = bedCols, format = "bed")
  } else if (!is.null(gff)) {
    # import gff using rtracklayer
    gr <- rtracklayer::import(gff)
  }

  # set seqlevelsStyle
  GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"

  # return
  gr
}
