# basic methods --------------------------------------------

#'
#' @rdname gdb
#' @usage NULL
#' @export
#'
setMethod("show", "gdb", function(object) {
  cat("rvat gdb object\n", "Path:", getGdbPath(object), "\n")
})

#' @rdname gdb
#' @usage NULL
#' @export
#'
setMethod("close", signature = "gdb", definition = function(con) {
  DBI::dbDisconnect(con)
})


# constructors --------------------------------------------

#' @rdname gdb
#' @usage NULL
#' @export

gdb <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("'%s' doesn't exist.", path), call. = FALSE)
  }
  tryCatch(
    {
      con <- DBI::dbConnect(duckdb::duckdb(), path)
    },
    error = function(e) {
      stop(sprintf("Invalid gdb path '%s'", path), call. = FALSE)
    }
  )
  new("gdb", con)
}


gdb_init <- function(path) {
  tryCatch(
    {
      con <- DBI::dbConnect(duckdb::duckdb(), path)
    },
    error = function(e) {
      stop(sprintf("Invalid gdb path '%s'", path), call. = FALSE)
    }
  )
  new("gdb", con)
}

# get metadata ------------------------------------------

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("listAnno", signature = "gdb", definition = function(object) {
  DBI::dbGetQuery(object, "select * from anno")
})

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("listCohort", signature = "gdb", definition = function(object) {
  DBI::dbGetQuery(object, "select * from cohort")
})

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("getRvatVersion", signature = "gdb", definition = function(object) {
  value <- DBI::dbGetQuery(
    object,
    "select * from meta where name = 'rvatVersion'"
  )$value
  if (length(value) == 0L) {
    value <- NA_character_
  }
  value
})

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("getGdbId", signature = "gdb", definition = function(object) {
  value <- DBI::dbGetQuery(object, "select * from meta where name = 'id'")$value
  if (length(value) == 0L) {
    value <- NA_character_
  }
  value
})

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("getGdbPath", signature = "gdb", definition = function(object) {
  DBI::dbGetInfo(object)$dbname
})

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("getCreationDate", signature = "gdb", definition = function(object) {
  value <- DBI::dbGetQuery(
    object,
    "select * from meta where name = 'creationDate'"
  )$value
  if (length(value) == 0L) {
    value <- NA_character_
  }
  value
})

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("getGenomeBuild", signature = "gdb", definition = function(object) {
  value <- DBI::dbGetQuery(
    object,
    "select * from meta where name = 'genomeBuild'"
  )$value
  if (length(value) == 0L) {
    value <- NA_character_
  }
  value
})


# extractRanges --------------------------------------------------

#' @rdname extractRanges
#' @usage NULL
#' @export
setMethod(
  "extractRanges",
  signature = "gdb",
  definition = function(object, ranges, padding = 250L) {
    # add padding
    if (!is.numeric(padding) || padding < 0L) {
      stop("`padding` parameter should be a positive value", call. = FALSE)
    }

    # if ranges is a data.frame, it should contain 'CHROM', 'start' and 'end' columns
    if (is.data.frame(ranges)) {
      if (!all(c("CHROM", "start", "end") %in% colnames(ranges))) {
        stop(
          "'CHROM', 'start' and 'end' should be present in ranges-file.",
          call. = FALSE
        )
      }
    } else if (is(ranges, "GRanges")) {
      # if ranges is a GRanges object, convert to a data.frame
      ranges <- data.frame(
        CHROM = as.character(seqnames(ranges)),
        start = start(ranges),
        end = end(ranges),
        stringsAsFactors = FALSE
      )
    } else {
      stop(
        "ranges should be either a data.frame or a GRanges object",
        call. = FALSE
      )
    }

    # check whether supplied ranges uses chr.. notation in CHROM field
    # update based on format used in gdb if needed
    chroms <- unique(ranges$CHROM)
    chroms_var <- unique(getAnno(object, "var_ranges", fields = "CHROM")$CHROM)
    chrom_format_with_chr <- grepl("^chr", chroms_var[1])
    if (grepl("^chr", chroms[1])) {
      if (!chrom_format_with_chr) {
        ranges$CHROM <- gsub("chr", "", ranges$CHROM, fixed = TRUE)
      }
    } else if (chrom_format_with_chr) {
      ranges$CHROM <- paste0("chr", ranges$CHROM)
    }
    chroms <- unique(ranges$CHROM)

    # if there is no chromosome overlap between ranges and gdb, raise a warning
    if (!any(chroms %in% chroms_var)) {
      warning(
        "None of the specified chromosomes are present in the gdb..",
        call. = FALSE
      )
      return(vector(mode = "character", length = 0L))
    }

    # perform queries in chunks of 500, padding is added for variants
    # where length(REF) > 1
    chunks <- split(
      seq_len(nrow(ranges)),
      ceiling(seq_len(nrow(ranges)) / 500)
    )
    run_query <- function(x, ranges, gdb, padding) {
      query <- sprintf(
        "(CHROM = '%s' and POS >= %s and POS <= %s)",
        ranges[x, ][["CHROM"]],
        ranges[x, ][["start"]] - padding,
        ranges[x, ][["end"]] + padding
      )
      query <- paste("select * from var where", paste(query, collapse = " or "))
      query <- DBI::dbGetQuery(gdb, query)
      query
    }
    queries <- do.call(
      rbind,
      lapply(
        chunks,
        FUN = run_query,
        ranges = ranges,
        gdb = object,
        padding = padding
      )
    )
    rownames(queries) <- NULL

    if (is.null(queries) || nrow(queries) == 0L) {
      return(character(0))
    }

    # use GRanges to identify overlaps taking into account variants where length(REF) > 1
    gr <- GenomicRanges::GRanges(
      seqnames = queries$CHROM,
      ranges = IRanges::IRanges(
        start = queries$POS,
        end = queries$POS + stringr::str_length(queries$REF) - 1L
      ),
      VAR_id = queries$VAR_id
    )
    names(gr) <- queries$VAR_id

    ranges <- GenomicRanges::GRanges(
      seqnames = ranges$CHROM,
      ranges = IRanges::IRanges(
        start = ranges$start,
        end = ranges$end
      )
    )
    overlaps <- GenomicRanges::findOverlaps(gr, ranges)

    # return overlapping VAR ids
    VAR_id <- as.character(gr[S4Vectors::queryHits(overlaps)]$VAR_id)
    VAR_id
  }
)
