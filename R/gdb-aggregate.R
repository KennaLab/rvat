#'  Aggregate genotypes into a single (burden) score for each individual.
#'
#'  Returns an aggregate of genotypes for each individual.
#'  Note, the [`gdb`] implementation is described here, `aggregate` can also be run directly on a
#'  [`genoMatrix`] object as described in the [`genoMatrix`] documentation.
#'  The specified genetic model, weights, MAF-weighting are taken into account when aggregating.
#'  Aggregates are written to disk in the [`aggdb`] format, which can be used as input
#'  for [`assocTest-aggdb`] to perform gene set burden analyses.
#'
#' @param x A [`gdb`] object.
#' @param cohort If a valid cohort name is provided,
#' the uploaded data for this cohort is used to filter and annotate the genotypes.
#' If not specified, all samples in the gdb will be loaded.
#' @param varSet A [`varSetList`] or [`varSetFile`] object.
#' @param VAR_id A vector of VAR_ids, alternatively the varSet parameter can be specified.
#' The `memlimit` argument controls how many variants to aggregate at a time.
#' @param pheno colData field to test as response variable, although not used within this method,
#' this can be useful to filter samples which have missing data for the response variable.
#' @param memlimit Maximum number of variants to load at once (if `VAR_id` is specified).
#' @param geneticModel Which genetic model to apply? ('allelic', 'recessive' or 'dominant').
#' Defaults to `allelic`.
#' @param imputeMethod Which imputation method to apply? ('meanImpute' or 'missingToRef').
#' Defaults to `meanImpute`.
#' @param MAFweights MAF weighting method.
#' Currently Madsen-Browning ('mb') is implemented, by default no MAF weighting is applied.
#' @param checkPloidy Version of the human genome to use when assigning variant ploidy (diploid, XnonPAR, YnonPAR).
#' Accepted inputs are 'GRCh37', 'hg19', 'GRCh38', 'hg38'.
#' If not specified, the genome build in the [`gdb`] will be used if available
#' (included if the `genomeBuild` parameter was set in [`buildGdb`]).
#' Otherwise, if the genome build is not included in the gdb metadata, and no value is provided,
#' then all variants are assigned the default ploidy of "diploid".
#' @param keep Vector of sample IDs to keep.
#' Defaults to `NULL`, in which case all samples are kept.
#' @param output Output file path for results.
#' Defaults to `NULL`, in which case results are not written.
#' @param signif Number of significant digits to store. Defaults to 6.
#' @param minCallrateVar Minimum genotype rate for variant retention.
#' @param maxCallrateVar Maximum genotype rate for variant retention.
#' @param minCallrateSM Minimum genotype rate for sample retention.
#' @param maxCallrateSM Maximum genotype rate for sample retention.
#' @param minMAF Minimum minor allele frequency for variant retention.
#' @param maxMAF Maximum minor allele frequency for variant retention.
#' @param minMAC Minimum minor allele count for variant retention.
#' @param maxMAC Maximum minor allele count for variant retention.
#' @param minCarriers Minimum carrier count for variant retention.
#' @param maxCarriers Maximum carrier count for variant retention.
#' @param minCarrierFreq Minimum carrier frequency for variant retention.
#' @param maxCarrierFreq Maximum carrier frequency for variant retention.
#' @param overWrite Should the output file be overwritten if it already exists?
#' Defaults to `FALSE`.
#' @param verbose Should the function be verbose? (TRUE/FALSE).
#' Defaults to `TRUE`.
#' @param strict Should strict checks be performed?
#' Defaults to `TRUE`. Strict tests currently includes
#' checking whether supplied varSetFile/varSetList was generated from the
#' same gdb as specified in `object`.
#'
#' @example inst/examples/example-gdb-aggregate.R
#'
#' @export
setMethod(
  "aggregate",
  signature = signature(x = "gdb"),
  definition = function(
    x,
    cohort = "SM",
    varSet = NULL,
    VAR_id = NULL,
    pheno = NULL,
    memlimit = 1000L,
    geneticModel = "allelic",
    imputeMethod = "meanImpute",
    MAFweights = "none",
    checkPloidy = NULL,
    keep = NULL,
    output = NULL,
    signif = 6L,
    minCallrateVar = 0,
    maxCallrateVar = 1,
    minCallrateSM = 0,
    maxCallrateSM = 1,
    minMAF = 0,
    maxMAF = 1,
    minMAC = 0,
    maxMAC = Inf,
    minCarriers = 0,
    maxCarriers = Inf,
    minCarrierFreq = 0,
    maxCarrierFreq = 1,
    overWrite = FALSE,
    verbose = TRUE,
    strict = TRUE
  ) {
    # input validation
    .aggregate_gdb_validate_input(as.list(environment()))

    # check gdb id
    if (!is.null(varSet) && strict) {
      .check_gdb_ids(x, varSet, minVersion = "0.3.0")
    }

    # if VAR_id is specified, generate varSet from VAR_ids
    if (!is.null(VAR_id)) {
      varSet <- .varsTovarSetList(VAR_id, chunkSize = memlimit)
    }

    # check output (aggdb) and overWrite
    # establish connection to output
    .check_output(output, overWrite = overWrite, verbose = verbose)
    inputs <- .aggregate_prepare_aggdb_inputs(as.list(environment()), gdb = x)
    aggdb <- .init_aggdb(
      output = output,
      metadata = inputs[["metadata"]],
      params = inputs[["params"]],
      samples = inputs[["samples"]],
      overWrite = overWrite,
      verbose = verbose
    )
    on.exit(close(aggdb), add = TRUE)

    # check if varSet doesn't contain duplicated units
    units <- listUnits(varSet)
    # varSet shouldn't contain multiple annotations
    # (which would result in duplicated units)
    if (anyDuplicated(units) != 0L) {
      stop(
        "Only 1 annotation per unit should be included in the varSet.",
        call. = FALSE
      )
    }

    # generate aggregates per unit and commit to aggdb
    .aggregate_gdb_generate_by_unit(
      aggdb = aggdb,
      gdb = x,
      units = units,
      varSet = varSet,
      cohort = cohort,
      pheno = pheno,
      keep = keep,
      imputeMethod = imputeMethod,
      geneticModel = geneticModel,
      MAFweights = MAFweights,
      checkPloidy = checkPloidy,
      signif = signif,
      minCallrateSM = minCallrateSM,
      maxCallrateSM = maxCallrateSM,
      minCallrateVar = minCallrateVar,
      maxCallrateVar = maxCallrateVar,
      minMAF = minMAF,
      maxMAF = maxMAF,
      minMAC = minMAC,
      maxMAC = maxMAC,
      minCarriers = minCarriers,
      maxCarriers = maxCarriers,
      minCarrierFreq = minCarrierFreq,
      maxCarrierFreq = maxCarrierFreq,
      verbose = verbose
    )

    invisible(NULL)
  }
)

.aggregate_gdb_validate_input <- function(args) {
  # type checks
  check_wrapper(check_character, args, "cohort")
  check_wrapper(check_character, args, "pheno", allow_null = TRUE)
  check_wrapper(check_number_whole, args, "memlimit")
  check_wrapper(check_character, args, "keep", allow_null = TRUE)
  check_wrapper(check_character, args, "output", allow_null = TRUE)
  check_wrapper(check_character, args, "checkPloidy", allow_null = TRUE)
  check_wrapper(check_number_whole, args, "signif")
  check_wrapper(check_bool, args, "overWrite")
  check_wrapper(check_bool, args, "verbose")
  check_wrapper(check_bool, args, "strict")

  # either varSet of VAR_id should be specified
  if (
    (!is.null(args[["varSet"]]) && !is.null(args[["VAR_id"]])) ||
      (is.null(args[["VAR_id"]]) && is.null(args[["varSet"]]))
  ) {
    stop("Either of `varSet` or `VAR_id` should be specified", call. = FALSE)
  }

  # if varSet is specified, it should be of type varSetList/varSetFile
  if (
    !is.null(args[["varSet"]]) &&
      !(is(args[["varSet"]], "varSetFile") ||
        is(args[["varSet"]], "varSetList"))
  ) {
    stop("`varSet` should be of type varSetFile or varSetList.", call. = FALSE)
  }

  # check VAR_id if specified
  .gdb_check_varid(args[["VAR_id"]])

  # use recode validator for geneticModel/MAFweights/imputeMethod validation
  .recode_validate_parameters(args)

  # sample and variant filters
  .validate_GT_filter_params(args)

  # check check lengths
  check_length(
    args[["cohort"]],
    arg = "cohort",
    equal = 1L
  )

  check_length(
    args[["pheno"]],
    arg = "pheno",
    equal = 1L,
    allow_null = TRUE
  )
  check_length(
    args[["output"]],
    arg = "output",
    equal = 1L,
    allow_null = TRUE
  )
  check_length(
    args[["checkPloidy"]],
    arg = "checkPloidy",
    equal = 1L,
    allow_null = TRUE
  )
  check_positive(args[["memlimit"]], arg = "memlimit")

  invisible(NULL)
}

.validate_GT_filter_params <- function(args) {
  # type checks
  check_wrapper(check_number_decimal, args, "minCallrateVar")
  check_wrapper(check_number_decimal, args, "maxCallrateVar")
  check_wrapper(check_number_decimal, args, "minCallrateSM")
  check_wrapper(check_number_decimal, args, "maxCallrateSM")
  check_wrapper(check_number_decimal, args, "minMAF")
  check_wrapper(check_number_decimal, args, "maxMAF")
  check_wrapper(check_number_decimal, args, "minMAC")
  check_wrapper(check_number_decimal, args, "maxMAC")
  check_wrapper(check_number_decimal, args, "minCarriers")
  check_wrapper(check_number_decimal, args, "maxCarriers")
  check_wrapper(check_number_decimal, args, "minCarrierFreq")
  check_wrapper(check_number_decimal, args, "maxCarrierFreq")

  # check values
  check_number_between(
    args[["minCallrateVar"]],
    0.0,
    1.0,
    including = TRUE,
    arg = "minCallrateVar"
  )
  check_number_between(
    args[["maxCallrateVar"]],
    0.0,
    1.0,
    including = TRUE,
    arg = "maxCallrateVar"
  )
  check_number_between(
    args[["minCallrateSM"]],
    0.0,
    1.0,
    including = TRUE,
    arg = "minCallrateSM"
  )
  check_number_between(
    args[["maxCallrateSM"]],
    0.0,
    1.0,
    including = TRUE,
    arg = "maxCallrateSM"
  )
  check_number_between(
    args[["minMAF"]],
    0.0,
    1.0,
    including = TRUE,
    arg = "minMAF"
  )
  check_number_between(
    args[["maxMAF"]],
    0.0,
    1.0,
    including = TRUE,
    arg = "maxMAF"
  )
  check_positive(args[["minMAC"]], arg = "minMAC")
  check_positive(args[["maxMAC"]], arg = "maxMAC")
  check_positive(args[["minCarriers"]], arg = "minCarriers")
  check_positive(args[["maxCarriers"]], arg = "maxCarriers")
  check_number_between(
    args[["minCarrierFreq"]],
    0.0,
    1.0,
    arg = "minCarrierFreq"
  )
  check_number_between(
    args[["maxCarrierFreq"]],
    0.0,
    1.0,
    arg = "maxCarrierFreq"
  )

  invisible(NULL)
}

.aggregate_prepare_aggdb_inputs <- function(args, gdb) {
  params <- list(
    cohort = if (is.null(args[["cohort"]])) NA_character_ else args[["cohort"]],
    pheno = if (is.null(args[["pheno"]])) NA_character_ else args[["pheno"]],
    geneticModel = args[["geneticModel"]],
    imputeMethod = args[["imputeMethod"]],
    MAFweights = args[["MAFweights"]],
    minCallrateVar = args[["minCallrateVar"]],
    maxCallrateVar = args[["maxCallrateVar"]],
    minCallrateSM = args[["minCallrateSM"]],
    maxCallrateSM = args[["maxCallrateSM"]],
    minMAF = args[["minMAF"]],
    maxMAF = args[["maxMAF"]],
    minMAC = args[["minMAC"]],
    maxMAC = args[["maxMAC"]],
    minCarriers = args[["minCarriers"]],
    maxCarriers = args[["maxCarriers"]],
    minCarrierFreq = args[["minCarrierFreq"]],
    maxCarrierFreq = args[["maxCarrierFreq"]]
  )
  params <- data.frame(
    name = names(params),
    value = as.character(unlist(params, use.names = FALSE)),
    stringsAsFactors = FALSE
  )

  metadata <- .gdb_list_metadata(gdb)
  metadata <- data.frame(
    name = names(metadata),
    value = as.character(unlist(metadata, use.names = FALSE)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  samples <- getCohort(gdb, cohort = args[["cohort"]], fields = "IID")$IID

  list(
    metadata = metadata,
    params = params,
    samples = samples
  )
}

.aggregate_gdb_generate_by_unit <- function(
  aggdb,
  gdb,
  units,
  varSet,
  cohort,
  pheno,
  keep,
  imputeMethod,
  geneticModel,
  MAFweights,
  checkPloidy,
  signif,
  minCallrateSM,
  maxCallrateSM,
  minCallrateVar,
  maxCallrateVar,
  minMAF,
  maxMAF,
  minMAC,
  maxMAC,
  minCarriers,
  maxCarriers,
  minCarrierFreq,
  maxCarrierFreq,
  verbose
) {
  lapply(units, FUN = function(unit) {
    # load genotypes
    varset <- getVarSet(varSet, unit = unit)
    if (length(varset) != 1L) {
      stop(sprintf(
        "Expected exactly 1 varSet for unit '%s', but found %s.",
        unit,
        length(varset)
      ))
    }
    varset <- varset[[1L]]
    GT <- getGT(
      gdb,
      cohort = cohort,
      varSet = varset,
      checkPloidy = checkPloidy,
      verbose = verbose,
      strict = FALSE
    )

    GT <- .filter_GT_samples_variants(
      GT = GT,
      unit = unit,
      minCallrateSM = minCallrateSM,
      maxCallrateSM = maxCallrateSM,
      minCallrateVar = minCallrateVar,
      maxCallrateVar = maxCallrateVar,
      minMAF = minMAF,
      maxMAF = maxMAF,
      minMAC = minMAC,
      maxMAC = maxMAC,
      minCarriers = minCarriers,
      maxCarriers = maxCarriers,
      minCarrierFreq = minCarrierFreq,
      maxCarrierFreq = maxCarrierFreq,
      geneticModel = geneticModel,
      pheno = pheno,
      keep = keep,
      return_value_if_empty = NULL,
      verbose = verbose
    )
    if (is.null(GT)) {
      return(NULL)
    }

    # recode
    GT <- recode(
      GT,
      imputeMethod = imputeMethod,
      geneticModel = geneticModel,
      MAFweights = MAFweights
    )

    # generate aggregates
    counts <- aggregate(
      GT,
      returnGT = FALSE,
      checkMissing = FALSE
    )
    names(counts) <- colnames(GT)
    samples <- listSamples(aggdb)
    counts <- counts[samples]

    # compress and commit to aggdb
    counts_blob <- serialize(
      signif(unname(counts), signif),
      connection = NULL
    )
    counts_blob <- list(memCompress(counts_blob, type = "gzip"))

    DBI::dbExecute(
      aggdb,
      "INSERT INTO aggregates (unit, aggregate) VALUES (?, ?)",
      params = list(unit, counts_blob)
    )
    invisible(NULL)
  })
  invisible(NULL)
}


.filter_GT_samples_variants <- function(
  GT,
  unit,
  minCallrateSM,
  maxCallrateSM,
  minCallrateVar,
  maxCallrateVar,
  minMAF,
  maxMAF,
  minMAC,
  maxMAC,
  minCarriers,
  maxCarriers,
  minCarrierFreq,
  maxCarrierFreq,
  geneticModel,
  pheno,
  keep,
  return_value_if_empty,
  verbose
) {
  # sample filtering
  keepSamples <- .filterGT_samples(
    object = GT,
    minCallrateSM = minCallrateSM,
    maxCallrateSM = maxCallrateSM,
    pheno = pheno,
    keep = keep,
    returnGT = FALSE,
    verbose = verbose
  )
  if (sum(keepSamples) == 0L) {
    if (verbose) {
      message(sprintf("No samples left after filtering for unit '%s'.", unit))
    }
    return(return_value_if_empty)
  }
  GT <- GT[, keepSamples]
  GT <- flipToMinor(GT)

  # variant filtering
  keepGeno <- .filterGT_vars(
    object = GT,
    minCallrateVar = minCallrateVar,
    maxCallrateVar = maxCallrateVar,
    minMAF = minMAF,
    maxMAF = maxMAF,
    minMAC = minMAC,
    maxMAC = maxMAC,
    minCarriers = minCarriers,
    maxCarriers = maxCarriers,
    minCarrierFreq = minCarrierFreq,
    maxCarrierFreq = maxCarrierFreq,
    geneticModel = geneticModel,
    returnGT = FALSE,
    filterWeights = TRUE,
    verbose = verbose
  )

  if (sum(keepGeno) == 0L) {
    if (verbose) {
      message(sprintf(
        "No variants left after filtering for unit '%s'.",
        unit
      ))
    }
    return(return_value_if_empty)
  }
  GT <- GT[keepGeno, ]

  # return
  GT
}
