# validate inputs -------------------------------------------------------------
.assoctest_validate_input <- function(args) {
  # type checks ------------------------------------------------
  check_wrapper(check_character, args, "pheno")
  check_wrapper(check_character, args, "test")
  check_wrapper(check_character, args, "name")
  check_wrapper(check_bool, args, "continuous")
  check_wrapper(check_bool, args, "singlevar")
  if (!is.null(args[["covar"]])) {
    is_char <- is.character(args[["covar"]])
    is_list <- is.list(args[["covar"]]) &&
      all(vapply(args[["covar"]], is.character, logical(1L), USE.NAMES = FALSE))

    if (!is_char && !is_list) {
      stop(
        "`covar` must be either a character vector or a list ",
        "of character vectors.",
        call. = FALSE
      )
    }
  }
  check_wrapper(check_character, args, "offset", allow_null = TRUE)
  check_wrapper(check_character, args, "geneticModel")
  check_wrapper(check_character, args, "imputeMethod", allow_null = TRUE)
  check_wrapper(check_character, args, "MAFweights")
  check_wrapper(check_number_whole, args, "maxitFirth")
  check_wrapper(check_character, args, "keep", allow_null = TRUE)
  check_wrapper(check_character, args, "output", allow_null = TRUE)
  check_wrapper(check_bool, args, "append")
  check_wrapper(check_bool, args, "returnDF")
  check_wrapper(check_character, args, "methodResampling", allow_null = TRUE)
  if (
    !(is.null(args[["resamplingMatrix"]]) ||
      is.matrix(args[["resamplingMatrix"]]))
  ) {
    stop(
      "`resamplingMatrix` should either be `NULL` or a numeric matrix",
      call. = FALSE
    )
  }
  if (!is.null(args[["resamplingFile"]])) {
    if (!is(args[["resamplingFile"]], "resamplingFile")) {
      stop(
        "`resamplingFile` should be an object of class 'resamplingFile'",
        call. = FALSE
      )
    }

    if (ncol(args[["object"]]) != args[["resamplingFile"]]@nSamples) {
      stop(
        "Number of samples in resamplingFile does not match the number of samples in the GT object.",
        call. = FALSE
      )
    }
  }
  check_wrapper(check_number_whole, args, "nResampling")
  if (
    !(is.character(args[["outputResampling"]]) ||
      is.logical(args[["outputResampling"]]))
  ) {
    stop(
      "`outputResampling` should be a filepath or a boolean (TRUE or FALSE)",
      call. = FALSE
    )
  }
  check_wrapper(
    check_number_whole,
    args,
    "memlimitResampling",
    allow_null = TRUE
  )
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
  check_wrapper(check_bool, args, "verbose")

  # additional checks ----------------------------------

  # pheno
  # phenotype should be present in genoMatrix
  cohort_fields <- colnames(colData(args[["object"]]))
  if (!all(args[["pheno"]] %in% cohort_fields)) {
    stop(
      sprintf(
        "'%s' is not present in `colData(GT)`",
        paste(
          args[["pheno"]][!args[["pheno"]] %in% cohort_fields],
          collapse = ","
        )
      ),
      call. = FALSE
    )
  }

  # covar
  # covars should be present in genoMatrix
  covar_unlisted <- unlist(args[["covar"]])
  if (!all(covar_unlisted %in% cohort_fields)) {
    missing_covars <- covar_unlisted[!covar_unlisted %in% cohort_fields]
    stop(
      sprintf(
        "the following `covar` fields were not found in colData: %s",
        paste(sQuote(missing_covars), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # test
  # check if specified tests are implemented
  # check if correct test are specified for binary/continuous traits and sv/rvb tests
  args[["test"]] <- unique(args[["test"]])
  if (!all(args[["test"]] %in% assocTest_tests)) {
    stop(
      sprintf(
        "The following tests are not valid: %s",
        paste(
          args[["test"]][!args[["test"]] %in% assocTest_tests],
          collapse = ","
        )
      ),
      call. = FALSE
    )
  }
  if (args[["singlevar"]]) {
    if (!args[["continuous"]]) {
      if (!all(args[["test"]] %in% assocTest_sv_bin_tests)) {
        warning(
          "The following tests were excluded since they are not available for binary singlevar tests: ",
          paste(
            args[["test"]][!args[["test"]] %in% assocTest_sv_bin_tests],
            collapse = ","
          ),
          call. = FALSE
        )
        args[["test"]] <- args[["test"]][
          args[["test"]] %in% assocTest_sv_bin_tests
        ]
      }
    } else {
      if (!all(args[["test"]] %in% assocTest_sv_cont_tests)) {
        warning(
          "The following tests were excluded since they are not ",
          "available for continuous singlevar tests: ",
          paste(
            args[["test"]][!args[["test"]] %in% assocTest_sv_cont_tests],
            collapse = ","
          ),
          call. = FALSE
        )
        args[["test"]] <- args[["test"]][
          args[["test"]] %in% assocTest_sv_cont_tests
        ]
      }
    }
    args[["test"]] <- assocTest_sv_tests[assocTest_sv_tests %in% args[["test"]]]
  } else {
    if (!args[["continuous"]]) {
      if (!all(args[["test"]] %in% assocTest_rvb_bin_tests)) {
        warning(
          "The following tests were excluded since they are not ",
          "available for binary rvb tests: ",
          paste(
            args[["test"]][!args[["test"]] %in% assocTest_rvb_bin_tests],
            collapse = ","
          ),
          call. = FALSE
        )
        args[["test"]] <- args[["test"]][
          args[["test"]] %in% assocTest_rvb_bin_tests
        ]
      }
    } else {
      if (!all(args[["test"]] %in% assocTest_rvb_cont_tests)) {
        warning(
          "The following tests were excluded since they are not available for continuous rvb tests: ",
          paste(
            args[["test"]][!args[["test"]] %in% assocTest_rvb_cont_tests],
            collapse = ","
          ),
          call. = FALSE
        )
        args[["test"]] <- args[["test"]][
          args[["test"]] %in% assocTest_rvb_cont_tests
        ]
      }
    }
  }

  # continuous
  # if FALSE, check if phenotype is coded 0,1
  # if TRUE, check if phenotype is numeric
  for (pheno_name in args[["pheno"]]) {
    if (!is.numeric(colData(args[["object"]])[[pheno_name]])) {
      stop(
        sprintf("Phenotype '%s' should be numeric.", pheno_name),
        call. = FALSE
      )
    }

    if (!args[["continuous"]]) {
      pheno_values <- colData(args[["object"]])[, pheno_name]
      non_missing_values <- pheno_values[!is.na(pheno_values)]

      if (!all(non_missing_values %in% c(0, 1))) {
        stop(
          sprintf(
            "Binary phenotype '%s' should be coded 0,1. If the phenotype is continuous, set `continuous = TRUE`",
            pheno_name
          ),
          call. = FALSE
        )
      }
    }
  }

  # offset
  # - max length = 1
  # - should be present in genoMatrix
  # - subset tests for which offset is implemented
  if (!is.null(args[["offset"]])) {
    if (length(args[["offset"]]) > 1L) {
      stop(
        "Currently at most one variable can be specified as offset",
        call. = FALSE
      )
    }

    if (
      !args[["offset"]] %in%
        colnames(colData(args[["object"]]))
    ) {
      stop(
        sprintf("Offset variable: %s is not available ", args[["offset"]]),
        call. = FALSE
      )
    }

    if (!all(args[["test"]] %in% assocTest_offset_tests)) {
      warning(
        sprintf(
          "Including an offset is not implemented for the following tests: %s",
          paste(
            args[["test"]][!args[["test"]] %in% assocTest_offset_tests],
            collapse = ","
          )
        ),
        call. = FALSE
      )
    }
    args[["test"]] <- args[["test"]][args[["test"]] %in% assocTest_offset_tests]
  }

  # geneticModel
  # - match args
  # - MAF/MAC filters cannot be applied if geneticModel is recessive/dominant
  # - geneticModel cannot be 'allelic' if genoMatrix is set to
  #   different geneticModel
  if (!all(args[["geneticModel"]] %in% c("allelic", "recessive", "dominant"))) {
    stop(
      "`geneticModel` should be allelic,recessive or dominant",
      call. = FALSE
    )
  }

  # if geneticModel is dominant or recessive, stop if user has specified MAF or MAC filters.
  GT_geneticmodel <- metadata(args[["object"]])$geneticModel
  if (
    GT_geneticmodel != "allelic" &&
      (args[["minMAF"]] > 0 ||
        args[["maxMAF"]] < 1 ||
        args[["minMAC"]] > 0 ||
        args[["maxMAC"]] < Inf)
  ) {
    stop(
      "MAC and/or MAF filters do not apply when geneticModel != 'allelic'.",
      call. = FALSE
    )
  }

  # not possible to set geneticModel == "allelic" if it's already set
  # to 'dominant' or 'recessive'
  if (
    GT_geneticmodel != "allelic" &&
      args[["geneticModel"]] != GT_geneticmodel
  ) {
    stop(
      "Current geneticModel should be 'allelic' in order to ",
      "apply dominant or recessive models.",
      call. = FALSE
    )
  }

  if (
    GT_geneticmodel != "allelic" &&
      (length(args[["geneticModel"]]) > 1L ||
        (length(args[["geneticModel"]]) == 1L &&
          args[["geneticModel"]] != GT_geneticmodel))
  ) {
    warning(
      sprintf(
        paste0(
          "The genoMatrix provided is set to '%s', ",
          "setting the `geneticModel` parameter to '%s'"
        ),
        GT_geneticmodel,
        GT_geneticmodel
      ),
      call. = FALSE
    )
    args[["geneticModel"]] <- GT_geneticmodel
  }


  # imputeMethod
  # - max. length = 1
  # - default is 'meanImpute'
  # - check if provided method is valid
  check_length(
    args[["imputeMethod"]],
    arg = "imputeMethod",
    equal = 1L,
    allow_null = TRUE
  )
  if (is.null(args[["imputeMethod"]]) && !args[["singlevar"]]) {
    args[["imputeMethod"]] <- "meanImpute"
  } else if (
    !is.null(args[["imputeMethod"]]) &&
      !args[["imputeMethod"]] %in% c("meanImpute", "missingToRef")
  ) {
    stop(
      "`imputeMethod` should be either 'meanImpute' or 'missingToRef'",
      call. = FALSE
    )
  }

  # MAFweights
  # - check if valid
  # - check if 1 method is provided
  check_character_contains(
    args[["MAFweights"]],
    all_of = c("none", "mb"),
    arg = "MAFweights"
  )

  # maxitFirth
  # - positive number
  check_positive(args[["maxitFirth"]], arg = "maxitFirth")

  # resamplingMatrix / methodResampling
  # - stop if resamplingMatrix/methodResampling is specified for singlevar
  #   (currently not implemented)
  # - check if methodResampling is valid

  if (
    args[["singlevar"]] &&
      (!is.null(args[["resamplingMatrix"]]) ||
        !is.null(args[["resamplingFile"]]) ||
        !is.null(args[["methodResampling"]]))
  ) {
    stop(
      "Resampling is currently not implemented for singlevar tests",
      call. = FALSE
    )
  }

  if (
    (!is.null(args[["resamplingMatrix"]]) ||
      !is.null(args[["resamplingFile"]]) ||
      !is.null(args[["methodResampling"]]))
  ) {
    test_resampling <- args[["test"]][
      args[["test"]] %in% assocTest_resampling_tests
    ]
    if (length(test_resampling) == 0) {
      stop(
        sprintf(
          paste0(
            "Resampling is not implemented for any of the specified tests.\n",
            "Resampling is currently implemented for: %s"
          ),
          paste(assocTest_resampling_tests, collapse = ",")
        ),
        call. = FALSE
      )
    }

    if (length(test_resampling) < length(args[["test"]])) {
      warning(
        sprintf(
          paste0(
            "Resampling is not implemented for all of the specified tests.\n",
            "Resampling is currently implemented for: %s"
          ),
          paste(assocTest_resampling_tests, collapse = ",")
        ),
        call. = FALSE
      )
    }
  } else {
    test_resampling <- NULL
  }

  if (
    !is.null(args[["methodResampling"]]) &&
      !args[["methodResampling"]] %in% c("permutation")
  ) {
    stop(
      "Currently only 'permutation' ",
      "is an accepted option for `methodResampling`",
      call. = FALSE
    )
  }

  if (
    !is.null(args[["resamplingMatrix"]]) &&
      ncol(args[["object"]]) != nrow(args[["resamplingMatrix"]])
  ) {
    stop(
      "Number of rows in resamplingMatrix should match with the ",
      "number of samples in the GT object.",
      call. = FALSE
    )
  }
  resampling_methods <- c(
    !is.null(args[["resamplingMatrix"]]),
    !is.null(args[["resamplingFile"]]),
    !is.null(args[["methodResampling"]])
  )

  if (sum(resampling_methods) > 1) {
    stop(
      "Only one of 'resamplingMatrix', 'resamplingFile', or 'methodResampling' can be specified.",
      call. = FALSE
    )
  }


  # nResampling
  # - positive number
  check_positive(args[["nResampling"]], arg = "nResampling")

  # memlimitResampling
  # - set to `nResampling` if null
  # - should be a positive number
  check_positive(
    args[["memlimitResampling"]],
    arg = "memlimitResampling",
    allow_null = TRUE
  )
  if (is.null(args[["memlimitResampling"]])) {
    args[["memlimitResampling"]] <- args[["nResampling"]]
  }

  # variant/sample filters
  # - check boundaries
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

  # final checks

  # stop if no tests left
  if (length(args[["test"]]) == 0) {
    stop("No applicable tests specified.", call. = FALSE)
  }

  # returns
  list(
    test = args[["test"]],
    test_resampling = test_resampling,
    geneticModel = args[["geneticModel"]],
    imputeMethod = args[["imputeMethod"]],
    memlimitResampling = args[["memlimitResampling"]]
  )
}

# filter -----------------------------------------------------------------------
.filterGT_samples <- function(
  object,
  minCallrateSM = 0,
  maxCallrateSM = 1,
  covar = NULL,
  pheno = NULL,
  keep = NULL,
  mask = NULL,
  returnGT = TRUE,
  verbose = TRUE,
  verboseKeep = FALSE,
  verboseCR = FALSE
) {
  keepSamples <- rep(TRUE, ncol(object))

  if (!is.null(keep)) {
    keep_keeplist <- colnames(object) %in% keep
    if (verbose || verboseKeep) {
      message(sprintf(
        "%s/%s samples present in `keep` are kept.",
        sum(keep_keeplist),
        ncol(object)
      ))
    }
    keepSamples <- keepSamples & keep_keeplist
  }

  if (minCallrateSM > 0 || maxCallrateSM < 1) {
    cr <- getCR(object, var = FALSE)
    keepCR <- cr >= minCallrateSM &
      cr <= maxCallrateSM
    if (!all(keepCR) && (verbose || verboseCR)) {
      message(sprintf(
        "%s/%s samples don't pass callRate filters.",
        sum(!keepCR),
        ncol(object)
      ))
    }
    keepSamples <- keepSamples & keepCR
  }

  if (!is.null(covar)) {
    keep_covar <- complete.cases(colData(object)[, covar])
    if (!all(keep_covar) && (verbose || verboseCR)) {
      message(sprintf(
        "%s/%s samples have missing values for one or more covariates, these will be dropped.",
        sum(!keep_covar),
        ncol(object)
      ))
    }
    keepSamples <- keepSamples & keep_covar
  }

  if (!is.null(pheno)) {
    keep_pheno <- !is.na(colData(object)[[pheno]])
    if (!all(keep_pheno) && verbose) {
      message(sprintf(
        "%s/%s samples have missing phenotypes, these will be dropped.",
        sum(!keep_pheno),
        ncol(object)
      ))
    }
    keepSamples <- keepSamples & keep_pheno
  }

  if (verbose) {
    message(sprintf(
      "%s/%s samples are kept in total.",
      sum(keepSamples),
      ncol(object)
    ))
  }

  if (returnGT) {
    object[, keepSamples]
  } else {
    keepSamples
  }
}

.filterGT_vars <- function(
  object,
  minCallrateVar = 0,
  maxCallrateVar = 1,
  minMAF = 0,
  maxMAF = 1,
  minMAC = 0,
  maxMAC = Inf,
  minCarriers = 0,
  maxCarriers = Inf,
  minCarrierFreq = 0,
  maxCarrierFreq = 1,
  returnGT = TRUE,
  filterWeights = TRUE,
  verbose = TRUE
) {
  keepGeno <- rep(TRUE, nrow(object))

  if (minMAF > 0.0 || maxMAF < 1.0) {
    maf <- getMAF(object)
    keepGeno <- keepGeno &
      maf >= minMAF &
      maf <= maxMAF
  }

  if (minMAC > 0 || maxMAC < Inf) {
    mac <- getMAC(object)
    keepGeno <- keepGeno &
      mac >= minMAC &
      mac <= maxMAC
  }

  if (minCallrateVar > 0.0 || maxCallrateVar < 1.0) {
    cr <- getCR(object)
    keepGeno <- keepGeno &
      cr >= minCallrateVar &
      cr <= maxCallrateVar
  }

  if (
    minCarriers > 0 ||
      maxCarriers < Inf ||
      minCarrierFreq > 0.0 ||
      maxCarrierFreq < 1.0
  ) {
    carriers <- getNCarriers(object)
    carrierFreq <- carriers / ncol(object)
    keepGeno <- keepGeno &
      carriers >= minCarriers &
      carriers <= maxCarriers &
      carrierFreq >= minCarrierFreq &
      carrierFreq <= maxCarrierFreq
  }

  if (filterWeights) {
    keepGeno <- keepGeno &
      .handle_weights(object, returnGT = FALSE, verbose = verbose)
  }

  if (verbose) {
    message(sprintf(
      "%s/%s variants are retained for analysis",
      sum(keepGeno),
      nrow(object)
    ))
  }

  if (returnGT) {
    object[keepGeno, ]
  } else {
    keepGeno
  }
}

.return_empty_results <- function(singlevar, returnDF) {
  if (singlevar) {
    if (returnDF) {
      return(as.data.frame(singlevarResult()))
    } else {
      return(singlevarResult())
    }
  } else {
    if (returnDF) {
      return(as.data.frame(rvbResult()))
    } else {
      return(rvbResult())
    }
  }
}

# init results ---------------------------------------------------------------
.assoctest_init_results_rvb <- function(
  GT,
  carriers,
  callRate,
  name,
  pheno,
  covar,
  geneticModel,
  MAFweights,
  test,
  continuous
) {
  if (continuous) {
    caseCarriers <- sum(carriers)
    ctrlCarriers <- NA_real_
    meanCaseScore = mean(GT[["aggregate"]])
    meanCtrlScore <- NA_real_
    caseN <- ncol(GT)
    ctrlN <- 0
    caseCallRate <- mean(callRate)
    ctrlCallRate <- NA_real_
  } else {
    caseCarriers <- sum(carriers[colData(GT)[, pheno] == 1] >= 1)
    ctrlCarriers <- sum(carriers[colData(GT)[, pheno] == 0] >= 1)
    meanCaseScore = mean(GT[, GT[[pheno]] == 1][["aggregate"]])
    meanCtrlScore = mean(GT[, GT[[pheno]] == 0][["aggregate"]])
    caseN <- sum(colData(GT)[, pheno] == 1, na.rm = TRUE)
    ctrlN <- sum(colData(GT)[, pheno] == 0, na.rm = TRUE)
    caseCallRate <- mean(callRate[colData(GT)[, pheno] == 1])
    ctrlCallRate <- mean(callRate[colData(GT)[, pheno] == 0])
  }

  results <- data.frame(
    unit = metadata(GT)$unit,
    varSetName = metadata(GT)$varSetName,
    cohort = metadata(GT)$cohort,
    name = name,
    pheno = pheno,
    covar = paste(covar, collapse = ","),
    geneticModel = metadata(GT)$geneticModel,
    MAFweight = if (MAFweights == "none") "1" else MAFweights,
    test = test,
    nvar = nrow(GT),
    caseCarriers = caseCarriers,
    ctrlCarriers = ctrlCarriers,
    caseN = caseN,
    ctrlN = ctrlN,
    meanCaseScore = meanCaseScore,
    meanCtrlScore = meanCtrlScore,
    caseCallRate = caseCallRate,
    ctrlCallRate = ctrlCallRate,
    effect = NA_real_,
    effectSE = NA_real_,
    effectCIlower = NA_real_,
    effectCIupper = NA_real_,
    OR = NA_real_,
    P = NA_real_,
    stringsAsFactors = FALSE
  )

  results
}

.assoctest_init_results_sv <- function(
  GT,
  name,
  pheno,
  covar,
  test,
  continuous,
  caseMAF,
  ctrlMAF,
  caseN,
  ctrlN
) {
  if (continuous) {
    caseMAC <- rep(
      Matrix::rowSums(
        assays(GT)$GT,
        na.rm = TRUE
      ),
      each = length(test)
    )
    ctrlMAC <- NA_real_
    caseCallRate <- caseN / ncol(GT)
    ctrlCallRate <- NA_real_
  } else {
    caseMAC <- rep(
      Matrix::rowSums(
        assays(GT)$GT[, colData(GT)[, pheno] == 1, drop = FALSE],
        na.rm = TRUE
      ),
      each = length(test)
    )
    ctrlMAC <- rep(
      Matrix::rowSums(
        assays(GT)$GT[, colData(GT)[, pheno] == 0, drop = FALSE],
        na.rm = TRUE
      ),
      each = length(test)
    )
    caseCallRate <- caseN / sum(colData(GT)[, pheno] == 1, na.rm = TRUE)
    ctrlCallRate <- ctrlN / sum(colData(GT)[, pheno] == 0, na.rm = TRUE)
  }
  if ("effectAllele" %in% colnames(rowData(GT))) {
    effectAllele <- rep(rowData(GT)$effectAllele, each = length(test))
  } else {
    effectAllele <- NA_character_
  }
  if ("otherAllele" %in% colnames(rowData(GT))) {
    otherAllele <- rep(rowData(GT)$otherAllele, each = length(test))
  } else {
    otherAllele <- NA_character_
  }

  results <- data.frame(
    varSetName = metadata(GT)$varSetName,
    cohort = as.character(metadata(GT)$cohort),
    name = name,
    geneticModel = metadata(GT)$geneticModel,
    VAR_id = rep(rownames(GT), each = length(test)),
    pheno = pheno,
    test = rep(test, times = nrow(GT)),
    covar = paste(covar, collapse = ","),
    caseN = unname(rep(caseN, each = length(test))),
    ctrlN = unname(rep(ctrlN, each = length(test))),
    caseMAF = unname(rep(caseMAF[rownames(GT)], each = length(test))),
    caseMAC = unname(caseMAC),
    caseCallRate = unname(rep(caseCallRate, each = length(test))),
    ctrlMAF = unname(rep(ctrlMAF[rownames(GT)], each = length(test))),
    ctrlMAC = unname(ctrlMAC),
    ctrlCallRate = unname(rep(ctrlCallRate, each = length(test))),
    effectAllele = unname(effectAllele),
    otherAllele = unname(otherAllele),
    stringsAsFactors = FALSE
  )
  results
}

# generate params ---------------------------------------------------------------
.assocTest_generate_params <- function(pheno, covar, geneticModel, MAFweights) {
  if (!is.list(covar)) {
    covar <- list(covar)
  }

  params <- expand.grid(
    pheno = pheno,
    covar_list = covar,
    geneticModel = geneticModel,
    MAFweights = MAFweights,
    stringsAsFactors = FALSE
  )

  params
}


# handle covariates ------------------------------------------------------------
.handle_covar <- function(coldata, covar, pheno) {
  # early return if NULL or covar is set to 1 (intercept-only)
  if (
    is.null(covar) ||
      length(covar) == 0L ||
      (length(covar) == 1L && covar == "1")
  ) {
    return(list(coldata = coldata, covar = NULL))
  }

  # check if covar includes duplicates + check if covariates are available
  covar <- .handle_covar_keep_available(coldata = coldata, covar = covar)

  # create dummies if covariates include categorial data (factor/character fields)
  covar_types <- unlist(lapply(coldata[, covar, drop = FALSE], FUN = class))
  if (sum(covar_types %in% c("character", "factor")) > 0L) {
    output <- .handle_covar_generate_dummies(
      coldata = coldata,
      covar = covar,
      pheno = pheno,
      covar_types = covar_types
    )
    covar <- output[["covar"]]
    coldata <- output[["coldata"]]
  }

  # drop covariates with zero covariance
  covar <- .handle_covar_zero_covariance(coldata = coldata, covar = covar)

  # return covar vector and coldata
  list(coldata = coldata, covar = covar)
}

.handle_weights <- function(object, returnGT = TRUE, verbose = TRUE) {
  keepGeno <- rep(TRUE, nrow(object))

  if (sum(is.na(rowData(object)$w)) > 0) {
    warning(
      sprintf(
        "For %s variants weights are missing, these variants are excluded.",
        sum(is.na(rowData(object)$w))
      ),
      call. = FALSE
    )
    keepGeno <- keepGeno & !is.na(rowData(object)$w)
  }

  if (sum(rowData(object)$w < 0, na.rm = TRUE) > 0) {
    warning(
      sprintf(
        "%s weights are < 0, these variants are excluded.",
        sum(rowData(object)$w < 0, na.rm = TRUE)
      ),
      call. = FALSE
    )
    keepGeno <- keepGeno & rowData(object)$w >= 0
  }

  if (returnGT) {
    object[keepGeno, ]
  } else {
    keepGeno
  }
}

.generate_models <- function(covar = NULL, offset = NULL, pheno) {
  # null models
  if (is.null(covar) && is.null(offset)) {
    null <- as.formula(sprintf("%s ~ 1", pheno))
  } else if (is.null(offset)) {
    null <- as.formula(paste(
      pheno,
      paste(covar, collapse = " + "),
      sep = " ~ "
    ))
  } else {
    null <- as.formula(sprintf(
      "%s ~ %s + offset(%s)",
      pheno,
      paste(covar, collapse = " + "),
      offset
    ))
  }

  # full models
  if (is.null(offset)) {
    model <- as.formula(paste(
      pheno,
      paste(c(covar, "aggregate"), collapse = " + "),
      sep = " ~ "
    ))
    model.nbinom <- as.formula(paste(
      "aggregate",
      paste(c(covar, pheno), collapse = " + "),
      sep = " ~ "
    ))
  } else {
    model <- as.formula(sprintf(
      "%s ~ %s + offset(%s) + aggregate",
      pheno,
      paste(covar, collapse = " + "),
      offset
    ))
    model.nbinom <- as.formula(sprintf(
      "aggregate ~ %s + offset(%s) + %s",
      paste(covar, collapse = " + "),
      offset,
      pheno
    ))
  }

  list(
    null = null,
    model = model,
    model.nbinom = model.nbinom
  )
}


.handle_covar_keep_available <- function(coldata, covar) {
  # check if duplicate covariates are specified, if so issue warning and keep unique
  if (anyDuplicated(covar)) {
    warning(
      "Duplicate covariates are specified, unique covariates are kept.",
      call. = FALSE
    )
  }
  covar <- unique(covar)

  # check if all covariates are available, return error if not all are available
  if (mean(covar %in% colnames(coldata)) < 1) {
    stop(
      sprintf(
        "The following covariate(s) are not available: %s",
        paste(covar[!covar %in% colnames(coldata)], collapse = ",")
      ),
      call. = FALSE
    )
  }

  # return covar
  covar
}

.handle_covar_generate_dummies <- function(coldata, covar, pheno, covar_types) {
  modmatrix <- as.data.frame(coldata)
  modmatrix[names(covar_types[covar_types == "character"])] <- lapply(
    modmatrix[names(covar_types[covar_types == "character"])],
    factor
  )
  modmatrix <- model.matrix.lm(
    as.formula(sprintf("%s ~ %s", pheno, paste(covar, collapse = "+"))),
    data = modmatrix,
    na.action = "na.pass"
  )
  modmatrix <- modmatrix[, -1]
  colnames(modmatrix) <- make.names(colnames(modmatrix))
  covar <- colnames(modmatrix)
  coldata <- cbind(
    coldata[, setdiff(colnames(coldata), colnames(modmatrix))],
    modmatrix
  )
  list(coldata = coldata, covar = covar)
}

.handle_covar_zero_covariance <- function(coldata, covar, warning = TRUE) {
  # drop covariates with zero variance.
  covarKeep <- c()
  for (i2 in covar) {
    if (var(coldata[, i2], na.rm = TRUE) > 0) covarKeep <- c(covarKeep, i2)
  }

  if (length(covar) > length(covarKeep) && warning) {
    warning(
      sprintf(
        "The following covariate(s) have zero covariance: %s",
        paste(covar[!covar %in% covarKeep], collapse = ",")
      ),
      call. = FALSE
    )
  }

  # return covariates to keep
  covarKeep
}

# return results ------------------------
.assoctest_return_results <- function(
  object,
  results,
  output,
  returnDF,
  singlevar,
  outputResampling,
  append
) {
  if (is.character(outputResampling) || outputResampling) {
    if (is.character(outputResampling)) {
      write.table(
        as.data.frame(results),
        file = outputResampling,
        sep = "\t",
        quote = FALSE,
        append = FALSE,
        row.names = FALSE
      )
    }
    return(results)
  }

  if (!is.null(output) || !returnDF) {
    result_type <- if (singlevar) "singlevarResult" else "rvbResult"
    results <- rvatResult(results, class = result_type)
    metadata(results)$rvatVersion <- as.character(packageVersion("rvat"))
    metadata(results)$gdbId <- metadata(object)$gdbId
    metadata(results)$genomeBuild <- metadata(object)$genomeBuild
    metadata(results)$creationDate <- as.character(round(
      Sys.time(),
      units = "secs"
    ))
  } else {
    cols <- if (singlevar) {
      names(columns_singlevarResults)
    } else {
      names(columns_rvbResults)
    }
    results <- results[, cols]
  }

  if (!is.null(output)) {
    writeResult(results, file = output, append = append)
  }
  results
}

# permutations ---------------------------
.permute_rvb_wrapper <- function(
  object,
  res,
  init_results,
  test,
  test_resampling,
  pheno,
  covar,
  model,
  null,
  continuous,
  resamplingMatrix,
  resamplingFile,
  nResampling,
  methodResampling,
  memlimitResampling,
  outputResampling
) {
  if (!is.null(resamplingMatrix)) {
    # run
    permResult <- .permute_rvb(
      GT = object,
      test = test_resampling,
      pheno = pheno,
      model = model,
      null = null,
      covar = covar,
      continuous = continuous,
      perms = resamplingMatrix,
      methodResampling = methodResampling
    )
  } else if (!is.null(resamplingFile)) {
    nResampling <- resamplingFile@nResampling
    methodResampling <- resamplingFile@methodResampling
    chunks <- rep(
      memlimitResampling,
      times = nResampling %/% memlimitResampling
    )
    if (sum(chunks) - nResampling != 0) {
      chunks <- c(chunks, nResampling - sum(chunks))
    }
    permResult <- list()
    for (k in 1:length(chunks)) {
      skip <- if (k == 1) 3 else cumsum(chunks)[(k - 1)] + 3
      perms <- t(as.matrix(data.table::fread(
        resamplingFile@path,
        skip = skip,
        nrows = memlimitResampling,
        header = FALSE
      )))
      colnames(perms) <- if (k == 1) {
        paste0("perm", 1:chunks[k])
      } else {
        paste0("perm", (cumsum(chunks)[(k - 1)] + 1):(cumsum(chunks)[(k)]))
      }
      prms <- .permute_rvb(
        GT = object,
        test = test_resampling,
        pheno = pheno,
        model = model,
        null = null,
        covar = covar,
        continuous = continuous,
        perms = perms,
        methodResampling = methodResampling
      )

      permResult[[k]] <- prms
    }
    permResult <- do.call(rbind, permResult)
  } else if (!is.null(methodResampling)) {
    # define chunks based on nResampling and memlimitResampling
    chunks <- rep(
      memlimitResampling,
      times = nResampling %/% memlimitResampling
    )
    if (sum(chunks) - nResampling != 0) {
      chunks <- c(chunks, nResampling - sum(chunks))
    }

    # generate permutations per chunk
    permResult <- list()
    for (i in 1:length(chunks)) {
      perms <- buildResamplingFile(
        nSamples = ncol(object),
        methodResampling = methodResampling,
        nResampling = chunks[[i]]
      )
      colnames(perms) <- if (i == 1) {
        paste0("perm", 1:chunks[i])
      } else {
        paste0("perm", (cumsum(chunks)[(i - 1)] + 1):(cumsum(chunks)[(i)]))
      }

      # run
      prms <- .permute_rvb(
        GT = object,
        test = test_resampling,
        pheno = pheno,
        model = model,
        null = null,
        covar = covar,
        continuous = continuous,
        perms = perms,
        methodResampling = methodResampling
      )

      permResult[[i]] <- prms
    }
    permResult <- do.call(rbind, permResult)
  }

  # format results, and return/output if outputResampling = TRUE or outputResampling = <filepath>
  if (is.character(outputResampling) || outputResampling) {
    permResult <- permResult[order(permResult$test), ]
    permResult <- cbind(
      DataFrame(
        varSetName = Rle(
          init_results$varSetName[1],
          lengths = nrow(permResult)
        ),
        cohort = Rle(init_results$cohort[1], lengths = nrow(permResult)),
        name = Rle(init_results$name[1], lengths = nrow(permResult)),
        unit = Rle(init_results$unit[1], lengths = nrow(permResult)),
        covar = Rle(init_results$covar[1], lengths = nrow(permResult)),
        geneticModel = Rle(
          init_results$geneticModel[1],
          lengths = nrow(permResult)
        ),
        MAFweight = Rle(init_results$MAFweight[1], lengths = nrow(permResult))
      ),
      DataFrame(permResult)
    )
    rownames(permResult) <- NULL
    return(permResult)
  } else {
    results <- rbind(res, get_perm_pvals(res, permResult))
    return(results)
  }
}

# misc --------------------------------------------------------------
.assocTest_calc_maf_per_pheno <- function(
  GT,
  pheno,
  continuous
) {
  if (metadata(GT)$geneticModel == "allelic") {
    if (continuous) {
      caseMAF <- getAF(GT)
      ctrlMAF <- rep(NA_real_, nrow(GT))
    } else {
      caseMAF <- getAF(GT[, colData(GT)[, pheno] == 1])
      ctrlMAF <- getAF(GT[, colData(GT)[, pheno] == 0])
    }
  } else {
    caseMAF <- rep(NA_real_, nrow(GT))
    ctrlMAF <- rep(NA_real_, nrow(GT))
  }
  names(caseMAF) <- rownames(GT)
  names(ctrlMAF) <- rownames(GT)

  list(caseMAF = caseMAF, ctrlMAF = ctrlMAF)
}

.assocTest_calc_n_per_pheno <- function(GT, pheno, continuous) {
  if (continuous) {
    caseN <- Matrix::rowSums(!is.na(assays(GT)$GT), na.rm = TRUE)
    ctrlN <- rep(0L, nrow(GT))
  } else {
    caseN <- Matrix::rowSums(
      !is.na(assays(GT)$GT[, colData(GT)[, pheno] == 1, drop = FALSE]),
      na.rm = TRUE
    )
    ctrlN <- Matrix::rowSums(
      !is.na(assays(GT)$GT[, colData(GT)[, pheno] == 0, drop = FALSE]),
      na.rm = TRUE
    )
  }

  names(caseN) <- rownames(GT)
  names(ctrlN) <- rownames(GT)

  list(caseN = caseN, ctrlN = ctrlN)
}
