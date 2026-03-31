#' getGT
#' @rdname getGT
#' @import GenomicRanges
#' @usage NULL
#' @export
setMethod(
  "getGT",
  signature = "gdb",
  definition = function(
    object,
    varSet = NULL,
    VAR_id = NULL,
    ranges = NULL,
    cohort = NULL,
    anno = NULL,
    annoFields = NULL,
    includeVarInfo = FALSE,
    checkPloidy = NULL,
    varSetName = "unnamed",
    unit = "unnamed",
    padding = 250L,
    verbose = TRUE,
    strict = TRUE
  ) {
    # check validity input
    arg <- as.list(environment())
    .getGT_validate_input(arg)
    rm(arg)

    # process varSet, VAR_id, or ranges to get a vector of VAR_ids
    # and their weights.
    handle_vars <- .getGT_handle_vars(
      object = object,
      varSet = varSet,
      VAR_id = VAR_id,
      ranges = ranges,
      padding = padding,
      unit = unit,
      varSetName = varSetName,
      strict = strict,
      verbose = verbose
    )
    VAR_id <- handle_vars$VAR_id
    w <- handle_vars$w
    unit <- handle_vars$unit
    varSetName <- handle_vars$varSetName
    rm(handle_vars)

    # retrieve and process cohort info
    cohort_name <- if (!is.null(cohort) && !is.character(cohort)) {
      deparse(substitute(cohort))
    } else {
      NULL
    }
    handle_cohort <- .getGT_handle_cohort(
      object = object,
      cohort = cohort,
      cohort_name = cohort_name,
      verbose = verbose
    )
    SM <- handle_cohort$SM
    cohort_name <- handle_cohort$cohort_name
    rm(handle_cohort)

    # determine ploidy for each variant,
    # checking for non-pseudoautosomal regions using either
    # user-provided build in `checkPloidy` or using ploidy stored in
    # gdb metadata
    if (is.null(checkPloidy)) {
      build_gdb <- getGenomeBuild(object)
      if (build_gdb %in% names(nonPAR)) checkPloidy <- build_gdb
    }
    ploidy <- .getGT_handle_ploidy(
      object = object,
      VAR_id = VAR_id,
      checkPloidy = checkPloidy,
      verbose = verbose
    )

    # query genotypes, decompress, and return as matrix
    GT <- .getGT_extract_geno(
      object = object,
      VAR_id = VAR_id,
      IID = SM$IID,
      verbose = verbose
    )
    VAR_id <- rownames(GT)
    w <- w[as.character(VAR_id)]
    if (is(ploidy, "data.frame")) {
      ploidy <- ploidy$ploidy[match(VAR_id, ploidy$VAR_id)]
    }

    # retrieve annotations, if specified through `anno` of `includeVarInfo`
    anno <- .getGT_handle_anno(
      object = object,
      VAR_id = VAR_id,
      includeVarInfo = includeVarInfo,
      anno = anno,
      annoFields = annoFields
    )

    # construct and return genoMatrix
    return(genoMatrix(
      GT = GT,
      SM = SM,
      anno = anno,
      VAR_id = VAR_id,
      w = w,
      ploidy = ploidy,
      varSetName = varSetName,
      unit = unit,
      cohortname = cohort_name,
      genomeBuild = getGenomeBuild(object),
      gdbpath = getGdbPath(object),
      gdbid = getGdbId(object),
      verbose = verbose
    ))
  }
)


# helpers -------------------------------------------------

.getGT_validate_input <- function(args) {
  # varSet and ranges parameters are checked in .getGT_handle_vars
  # ranges parameter is checked in extractRanges
  # Check if VAR_id is not NULL and is not a numeric or character vector
  .gdb_check_varid(args[["VAR_id"]])
  check_wrapper(
    check_character,
    args,
    "anno",
    length_equal = 1L,
    allow_null = TRUE
  )
  check_wrapper(check_character, args, "annoFields", allow_null = TRUE)
  check_wrapper(check_bool, args, "includeVarInfo", length_equal = 1L)
  check_wrapper(
    check_character,
    args,
    "checkPloidy",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(
    check_character,
    args,
    "varSetName",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(
    check_character,
    args,
    "unit",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(check_number_whole, args, "padding", length_equal = 1L)
  check_positive(args[["padding"]], arg = "padding")
  check_wrapper(check_bool, args, "verbose", length_equal = 1L)
  check_wrapper(check_bool, args, "strict", length_equal = 1L)

  # only one of varSet, VAR_id, or ranges can be specified
  sm <- sum(
    c(
      !is.null(args[["varSet"]]),
      !is.null(args[["VAR_id"]]),
      !is.null(args[["ranges"]])
    )
  )

  if (sm > 1L) {
    stop("Specify only one of `varSet`, `VAR_id`, or `ranges`.", call. = FALSE)
  }
  if (sm == 0L) {
    stop("Specify one of `varSet`, `VAR_id`, or `ranges`.", call. = FALSE)
  }

  # check if checkPloidy is either NULL or one of the supported genome builds
  if (
    !is.null(args[["checkPloidy"]]) && !args[["checkPloidy"]] %in% names(nonPAR)
  ) {
    stop(
      sprintf(
        paste0(
          "Unrecognized value for checkPloidy. Must be one of %s.\n",
          "Alternatively you can provide no value for checkPloidy ",
          "and manually reset values from the default of 'diploid' ",
          "using the resetPloidy function"
        ),
        paste(names(nonPAR), collapse = ",")
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

.getGT_handle_vars <- function(
  object,
  varSet,
  VAR_id,
  ranges,
  padding,
  varSetName,
  unit,
  strict,
  verbose
) {
  if (!is.null(varSet)) {
    # if varSet is provided, process
    processed <- .getGT_process_varSet(object, varSet, strict)
    VAR_id <- processed$VAR_id
    w <- processed$w
    unit <- processed$unit
    varSetName <- processed$varSetName

    # if ranges are supplied extract those, set weights to 1
  } else if (!is.null(ranges)) {
    VAR_id <- extractRanges(object, ranges = ranges, padding = padding)
    w <- rep(1.0, length(VAR_id))
    names(w) <- VAR_id

    ## if VAR_ids are supplied, select those, set weights to 1
  } else if (!is.null(VAR_id)) {
    VAR_id <- unique(sort(as.integer(VAR_id)))
    w <- rep(1.0, length(VAR_id))
    names(w) <- VAR_id
  }

  list(
    VAR_id = VAR_id,
    w = w,
    unit = unit,
    varSetName = varSetName
  )
}

.getGT_process_varSet <- function(object, varSet, strict) {
  # check if varSet was generated from the current gdb
  if (strict) {
    .check_gdb_ids(object, varSet, minVersion = "0.3.0")
  }

  # if a varSetList of length 1 is provided, extract the varSet
  if (is(varSet, "varSetList") && length(varSet) == 1L) {
    varSet <- varSet[[1]]
  } else if (is(varSet, "varSetList") && length(varSet) > 1L) {
    stop(
      "A varSetList with >1 varSets is supplied. ",
      "Only 1 varSet can be supplied to `getGT()`, use `collapseVarSetList` ",
      "to merge the records of a varSetList",
      call. = FALSE
    )
  } else if (!is(varSet, "varSet")) {
    stop(
      "varSet parameter should be of class `varSet` or ",
      "a `varSetList` of length 1",
      call. = FALSE
    )
  }

  # parse varSet
  VAR_id <- listVars(varSet)
  .gdb_check_varid(VAR_id)
  w <- as.numeric(listWeights(varSet))
  names(w) <- VAR_id

  list(
    VAR_id = VAR_id,
    w = w,
    unit = varSet@unit,
    varSetName = varSet@varSetName
  )
}

.getGT_handle_cohort <- function(
  object,
  cohort,
  cohort_name,
  verbose
) {
  if (is.null(cohort)) {
    # use SM table is no cohort is provided
    SM <- DBI::dbGetQuery(object, "select * from SM")
    cohort_name <- "SM"
  } else if (is.character(cohort) && length(cohort) == 1L) {
    # retrieve provided cohort
    SM <- getCohort(object, cohort, keepAll = TRUE)
    cohort_name <- cohort
  } else if (is.data.frame(cohort) || is(cohort, "DFrame")) {
    ## if cohort is supplied as either a data.frame or a DFrame,
    ## infer cohort name from name of supplied object or retrieve from
    ## metadata if class == DFrame
    if (is(cohort, "DFrame") && "name" %in% names(metadata(cohort))) {
      cohort_name <- metadata(cohort)$name
    } else {
      cohort_name <- if (!is.null(cohort_name)) cohort_name else "cohort"
    }

    # convert DFrame to data.frame
    cohort <- as.data.frame(cohort)

    # check if cohort IIDs are present in gdb and align with SM
    SM <- .align_cohort_SM(
      gdb = object,
      cohort = cohort
    )
  } else {
    # throw error if cohort is neither a cohort name or a data.frame/DFrame
    stop(
      "`cohort` should be either a cohort uploaded to the gdb ",
      "or a data.frame/DFrame containing cohort info",
      call. = FALSE
    )
  }

  list(
    SM = SM,
    cohort_name = cohort_name
  )
}

.getGT_handle_ploidy <- function(
  object,
  VAR_id,
  checkPloidy,
  verbose
) {
  # early return if checkPloidy is null
  if (is.null(checkPloidy) || is.na(checkPloidy) || checkPloidy == "NA") {
    return("diploid")
  }

  # retrieve varinfo
  varInfo <- DBI::dbGetQuery(
    object,
    sprintf(
      paste0(
        "select VAR_id, 'diploid' as ploidy, CHROM, POS, REF ",
        "from var where VAR_id in (%s);"
      ),
      paste(as.integer(VAR_id), collapse = ",")
    )
  )

  # return early if varInfo is empty
  if (nrow(varInfo) == 0L) {
    return(data.frame(
      VAR_id = character(),
      ploidy = character(),
      stringsAsFactors = FALSE
    ))
  }

  # convert varinfo to genomicRanges
  varInfo$end <- varInfo$POS + nchar(varInfo$REF) - 1L
  varInfo <- makeGRangesFromDataFrame(
    varInfo,
    seqnames.field = "CHROM",
    start.field = "POS",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  GenomeInfoDb::seqlevelsStyle(varInfo) <- "NCBI"

  # find overlap with nonPAR regions
  withCallingHandlers(
    {
      overlaps <- findOverlaps(varInfo, nonPAR[[checkPloidy]])
    },
    warning = suppressSeqinfowarning
  )

  ## update ploidy for variants that overlap nonPAR regions
  if (length(overlaps) > 0L && verbose) {
    message(sprintf(
      paste0(
        "Ploidy of non-pseudoautosomal regions of the ",
        "sex chromosomes are being set based on build %s"
      ),
      checkPloidy
    ))
  }

  varInfo$ploidy[queryHits(overlaps)] <- nonPAR[[
    checkPloidy
  ]]$ploidy[subjectHits(overlaps)]

  # return ploidy per variant
  data.frame(VAR_id = varInfo$VAR_id, ploidy = varInfo$ploidy)
}

.getGT_extract_geno <- function(
  object,
  VAR_id,
  IID,
  verbose
) {
  # query genotypes
  query <- sprintf(
    "select VAR_id, GT from dosage where VAR_id in (%s);",
    paste(as.integer(VAR_id), collapse = ",")
  )
  GT <- RSQLite::dbGetQuery(object, query)

  # check if all variants are found
  success <- sum(VAR_id %in% GT[, 1])
  if (success < length(VAR_id)) {
    warning(
      sprintf(
        "Retrieved genotypes for only %s of %s variants",
        success,
        length(VAR_id)
      ),
      call. = FALSE
    )
  }
  nvar <- nrow(GT)
  if (nvar == 0L) {
    empty_mat <- matrix(NA_real_, nrow = 0L, ncol = length(IID))
    colnames(empty_mat) <- IID
    return(empty_mat)
  }
  VAR_id <- GT[, 1]

  ## unpack genotypes and store as matrix
  GT <- lapply(GT[, 2], memDecompress, type = "gzip", asChar = FALSE)
  ascii_offset <- as.integer(charToRaw("0"))
  GT <- as.integer(unlist(GT)) - ascii_offset
  GT[GT > 2L] <- NA_real_
  GT <- matrix(GT, nrow = nvar, ncol = length(GT) / nvar, byrow = TRUE)
  rownames(GT) <- VAR_id
  if (verbose) {
    message(sprintf("Retrieved genotypes for %s variants", nvar))
  }

  # return
  GT
}

.getGT_handle_anno <- function(
  object,
  VAR_id,
  includeVarInfo,
  anno,
  annoFields
) {
  if (is.null(VAR_id)) {
    VAR_id <- integer(0)
  }
  if (includeVarInfo) {
    # extract var table
    anno <- getAnno(object, table = "var", VAR_id = VAR_id)
  } else if (!is.null(anno)) {
    # extract user-specified annotations
    anno <- getAnno(
      object,
      table = anno,
      fields = if (!is.null(annoFields)) annoFields else "*",
      VAR_id = VAR_id
    )
  }

  # stop if annotation contains duplicated VAR_ids
  if (!is.null(anno) && anyDuplicated(anno$VAR_id) != 0L) {
    stop(
      "Annotations should contain 1 row per VAR_id to be included ",
      "in a genoMatrix object.",
      call. = FALSE
    )
  }

  # return
  anno
}

suppressSeqinfowarning <- function(w) {
  if (any(grepl("sequence levels", w, fixed = TRUE))) {
    invokeRestart("muffleWarning")
  }
}
