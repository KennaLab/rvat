columns_summariseGeno <- list(
  VAR_id = "character",
  AF = "numeric",
  callRate = "numeric",
  geno0 = c("numeric", "integer"),
  geno1 = c("numeric", "integer"),
  geno2 = c("numeric", "integer"),
  hweP = "numeric"
)

setMethod(
  "summariseGeno",
  signature = signature(object = "gdb"),
  definition = function(
    object,
    cohort = "SM",
    varSet = NULL,
    VAR_id = NULL,
    pheno = NULL,
    memlimit = 1000L,
    geneticModel = "allelic",
    checkPloidy = NULL,
    keep = NULL,
    output = NULL,
    splitBy = NULL,
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
    strict = TRUE,
    verbose = TRUE
  ) {
    # input validation
    .summariseGeno_gdb_validate_input(as.list(environment()))

    # check gdb id
    if (!is.null(varSet) && strict) {
      .check_gdb_ids(object, varSet, minVersion = "0.3.0")
    }

    # if VAR_id is specified, generate varSet from VAR_ids
    if (!is.null(VAR_id)) {
      varSet <- .varsTovarSetList(VAR_id, chunkSize = memlimit)
    }

    # initialize output
    output_con <- .summariseGeno_gdb_init_output(output, splitBy = splitBy)
    if (!is.null(output_con)) {
      on.exit(close(output_con), add = TRUE)
    }

    # generate varsummary per unit
    sumgeno <- .summariseGeno_gdb_generate_by_unit(
      gdb = object,
      varSet = varSet,
      cohort = cohort,
      pheno = pheno,
      keep = keep,
      geneticModel = geneticModel,
      output_con = output_con,
      checkPloidy = checkPloidy,
      splitBy = splitBy,
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

    # return
    if (is.null(output_con)) {
      rownames(sumgeno) <- NULL
      sumgeno
    } else {
      invisible(NULL)
    }
  }
)

.summariseGeno_gdb_validate_input <- function(args) {
  # type checks
  check_wrapper(check_character, args, "cohort", length_equal = 1L)
  check_wrapper(
    check_character,
    args,
    "pheno",
    length_equal = 1L,
    allow_null = TRUE
  )
  check_wrapper(check_number_whole, args, "memlimit", length_equal = 1L)
  check_positive(args[["memlimit"]], arg = "memlimit")
  check_wrapper(
    check_character,
    args,
    "checkPloidy",
    length_equal = 1L,
    allow_null = TRUE
  )
  check_wrapper(
    check_character,
    args,
    "keep",
    allow_null = TRUE
  )
  check_wrapper(
    check_character,
    args,
    "output",
    length_equal = 1L,
    allow_null = TRUE
  )
  check_wrapper(
    check_character,
    args,
    "splitBy",
    length_equal = 1L,
    allow_null = TRUE
  )
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

  # use recode validator for geneticModel validation
  .recode_validate_parameters(args)

  # sample and variant filters
  .validate_GT_filter_params(args)

  # check if cohort fields are present in cohort
  cohort_fields <- DBI::dbListFields(args[["object"]], args[["cohort"]])
  if (!is.null(args[["pheno"]]) && !args[["pheno"]] %in% cohort_fields) {
    stop(
      sprintf(
        "The `pheno` field '%s' was not found in cohort '%s'.",
        args[["pheno"]],
        args[["cohort"]]
      ),
      call. = FALSE
    )
  }

  if (!is.null(args[["splitBy"]]) && !args[["splitBy"]] %in% cohort_fields) {
    stop(
      sprintf(
        "The `splitBy` field ('%s'), is not found in cohort '%s'.",
        paste(
          args[["splitBy"]],
          collapse = ","
        ),
        args[["cohort"]]
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

.summariseGeno_gdb_init_output <- function(output, splitBy) {
  if (is.null(output)) {
    return(NULL)
  }
  output_con <- gzcon(file(output, open = "wb"))
  if (is.null(splitBy)) {
    field_names <- names(columns_summariseGeno)
  } else {
    field_names <- c(
      "VAR_id",
      splitBy,
      names(columns_summariseGeno)[names(columns_summariseGeno) != "VAR_id"]
    )
  }
  write(
    paste(field_names, collapse = "\t"),
    file = output_con,
    append = FALSE
  )
  output_con
}

.summariseGeno_gdb_generate_by_unit <- function(
  gdb,
  varSet,
  cohort,
  pheno,
  keep,
  geneticModel,
  output_con,
  checkPloidy,
  splitBy,
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
  units_all <- unique(listUnits(varSet))
  results <- lapply(units_all, FUN = function(unit) {
    if (verbose) {
      message(sprintf("Analysing %s", unit))
    }

    # load genotypes
    varset <- getVarSet(varSet, unit = unit)[[1L]]
    GT <- getGT(
      gdb,
      cohort = cohort,
      varSet = varset,
      checkPloidy = checkPloidy,
      verbose = verbose,
      strict = FALSE
    )

    # perform variant/sample filtering
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

    # recode based on geneticModel
    GT <- recode(GT, geneticModel = geneticModel)

    if (!is.null(splitBy)) {
      # generate sumgeno per category
      splitByfct <- unique(colData(GT)[[splitBy]])
      splitByfct <- splitByfct[!is.na(splitByfct)]
      sumgeno <- lapply(splitByfct, FUN = function(fct) {
        sumgeno_fct <- summariseGeno(GT[, colData(GT)[[splitBy]] %in% fct])
        cbind(
          VAR_id = sumgeno_fct$VAR_id,
          splitBy = fct,
          sumgeno_fct[, colnames(sumgeno_fct) != "VAR_id"]
        )
      })
      sumgeno <- do.call(rbind, sumgeno)
      colnames(sumgeno)[2L] <- splitBy
      rownames(sumgeno) <- NULL
    } else {
      sumgeno <- summariseGeno(GT)
    }

    # write or return
    if (!is.null(output_con)) {
      write.table(
        sumgeno,
        file = output_con,
        append = TRUE,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = "\t"
      )
      NULL
    } else {
      sumgeno
    }
  })
  do.call(rbind, results)
}
