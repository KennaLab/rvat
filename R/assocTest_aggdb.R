columns_aggresults <- list(
  geneSetName = c("character", "Rle"),
  cohort = c("character", "Rle"),
  name = c("character", "Rle"),
  pheno = c("character", "Rle"),
  covar = c("character", "Rle"),
  test = c("character", "Rle"),
  geneSetSize = c("integer", "numeric"),
  genesObs = c("integer", "numeric"),
  caseN = c("integer", "numeric"),
  ctrlN = c("integer", "numeric"),
  meanCaseScore = c("numeric", "integer"),
  meanCtrlScore = c("numeric", "integer"),
  effect = c("numeric", "integer"),
  effectSE = c("numeric", "integer"),
  effectCIlower = c("numeric", "integer"),
  effectCIupper = c("numeric", "integer"),
  OR = c("numeric", "integer"),
  P = c("numeric", "integer")
)

#' assocTest-aggdb
#'
#' Run [`assocTest`] on a [`aggregateFile`] object.
#' See the main [`assocTest`] page for details.

#' @rdname assocTest-aggdb
#' @name assocTest-aggdb
#' @aliases assocTest,aggdb-method
#' @param object a [`aggdb`] object (generated using the [`aggregate`] method).
#' @param pheno field to test as response variable, the response variable
#' can either be binary (0/1) or continuous. If the response variable is continuous set
#' `continuous` to `TRUE`.
#' @param test Vector of statistical tests to run,
#' options include firth,glm,lm, and nbinom.
#' @param geneSet a [`geneSetList`] or [`geneSetFile`] object
#' @param gdb a [`gdb`] object, this will be used to load the cohort data.
#' @param cohort Cohort name (present in the gdb)
#' @param name Optional name of the analysis.
#' @param continuous Is the response variable continuous? (TRUE/FALSE). Defaults to `FALSE`.
#' @param covar Character vector of covariates, or list with multiple sets of covariates.
#' @param substractCovar Covariate from which aggregate should be substracted.
#' Useful when adjusting for total variant counts, by specifying the total variant count variable
#' here, the aggregate score of the gene set tested will be substracted from the total count variable.
#' @param dropUnits Optional, vector of units to exclude.
#' @param maxitFirth Maximum number of iterations to use for estimating firth confidence intervals.
#' @param keep vector of sample IDs to keep, defaults to `NULL`, in which case all samples are kept.
#' Defaults to `FALSE`.
#' @param output Output file path for results.
#' Defaults to `NULL`, in which case results are not written but returned as a `data.frame()`.
#' @param verbose Should the function be verbose? (TRUE/FALSE), defaults to `TRUE`.
#' @param strict Should strict checks be performed? Defaults to `TRUE`. Strict tests currently includes
#' checking whether supplied aggdb was generated from the same gdb as specified in `gdb`.
#' @examples
#' library(rvatData)
#' gdb <- gdb(rvat_example("rvatData.gdb"))
#'
#' # assocTest-aggdb allows for running association tests on pre-constructed aggregates.
#' # below we first generate the aggregates based on a varSetFile
#' varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
#' varset <- getVarSet(varsetfile, unit = c("NEK1", "SOD1", "ABCA4"), varSetName = "High")
#' aggfile <- tempfile()
#' aggregate(
#'   x = gdb,
#'   varSet = varset,
#'   maxMAF = 0.001,
#'   output = aggfile,
#'   verbose = FALSE
#' )
#'
#' # connect to aggdb, see ?aggdb for more details
#' aggregatefile <- aggdb(aggfile)
#'
#' # build an example genesetlist, see ?buildGeneSet for details
#' genesetlist <- buildGeneSet(
#'   list(
#'     "geneset1" = c("SOD1", "NEK1"),
#'     "geneset2" = c("ABCA4", "SOD1", "NEK1")
#'   )
#' )
#'
#' # perform association tests using assocTest
#' # note that this is very similar to running association tests on a genoMatrix (?`assocTest-genoMatrix`)
#' # or on a gdb (?`assocTest-gdb`). The main difference being that pre-constructed aggregates are used from
#' # the aggdb, and requires a genesetFile/geneSetList (?geneSetFile) to be provided to the `geneSet` argument.
#' aggAssoc <- assocTest(
#'   aggregatefile,
#'   gdb = gdb,
#'   test = c("glm", "firth"),
#'   cohort = "pheno",
#'   pheno = "pheno",
#'   geneSet = genesetlist,
#'   covar = paste0("PC", 1:4),
#'   verbose = FALSE
#' )
#'
#' @export
setMethod(
  "assocTest",
  signature = signature(object = "aggdb"),
  definition = function(
    object,
    pheno,
    test,
    geneSet,
    gdb,
    cohort = "SM",
    name = "none",
    continuous = FALSE,
    covar = NULL,
    substractCovar = NULL,
    dropUnits = NULL,
    maxitFirth = 1000L,
    keep = NULL,
    output = NULL,
    verbose = TRUE,
    strict = TRUE
  ) {
    # validate input
    arg_parsed <- .assocTest_aggdb_validate_input(as.list(environment()))
    test <- arg_parsed[["test"]]
    covar <- arg_parsed[["covar"]]
    rm(arg_parsed)

    # check gdb id
    if (strict) {
      .check_gdb_ids(gdb, object, minVersion = "0.3.0")
    }

    # initialize output
    output_con <- .assocTest_aggdb_init_output(output)
    if (!is.null(output_con)) on.exit(close(output_con), add = TRUE)

    # load and filter cohort
    cohort_df <- .assocTest_aggdb_load_cohort(
      object,
      gdb = gdb,
      cohort = cohort,
      pheno = pheno,
      continuous = continuous,
      keep = keep,
      verbose = verbose
    )
    # return empty results if no samples are left
    if (nrow(cohort_df) == 0L) {
      return(.assocTest_aggdb_return_empty_results())
    }

    # generate params data.frame
    params <- .assocTest_aggdb_generate_params(covar, pheno)

    # list genesets and units in aggdb
    genesets <- listGeneSets(geneSet)
    units_aggdb <- listUnits(object)

    # run analyses per geneset
    results <- .assocTest_aggdb_run_geneset_analysis(
      object,
      genesets = genesets,
      geneSet = geneSet,
      units_aggdb = units_aggdb,
      cohort_df = cohort_df,
      params = params,
      cohort_name = cohort,
      name = name,
      test = test,
      continuous = continuous,
      substractCovar = substractCovar,
      dropUnits = dropUnits,
      maxitFirth = maxitFirth,
      output_con = output_con,
      verbose = verbose
    )

    # return results
    if (is.null(output)) {
      results
    } else {
      invisible(NULL)
    }
  }
)

.assocTest_aggdb_validate_input <- function(args) {
  # type checks
  check_wrapper(check_character, args, "pheno", length_equal = 1L)
  check_wrapper(check_character, args, "test", length_min = 1L)
  check_wrapper(check_character, args, "cohort", length_equal = 1L)
  check_wrapper(check_character, args, "name", length_equal = 1L)
  if (!is.null(args[["covar"]])) {
    is_char <- is.character(args[["covar"]])
    is_list <- is.list(args[["covar"]]) &&
      all(vapply(args[["covar"]], is.character, logical(1), USE.NAMES = FALSE))

    if (!is_char && !is_list) {
      stop(
        "`covar` must be either a character vector or a list of character vectors.",
        call. = FALSE
      )
    }
  }
  check_wrapper(
    check_character,
    args,
    "substractCovar",
    length_equal = 1,
    allow_null = TRUE
  )
  check_wrapper(check_character, args, "dropUnits", allow_null = TRUE)
  check_wrapper(check_bool, args, "continuous")
  check_wrapper(check_bool, args, "verbose")
  check_wrapper(check_number_whole, args, "maxitFirth", length_equal = 1L)
  check_positive(args[["maxitFirth"]], arg = "maxitFirth")
  check_wrapper(check_character, args, "keep", allow_null = TRUE)
  check_wrapper(
    check_character,
    args,
    "output",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(check_bool, args, "verbose")
  check_wrapper(check_bool, args, "strict")

  # check classes
  if (!is(args[["gdb"]], "gdb")) {
    stop("`gdb` must be a valid gdb object.", call. = FALSE)
  }
  if (
    !is(args[["geneSet"]], "geneSetList") &&
      !is(args[["geneSet"]], "geneSetFile")
  ) {
    stop(
      "`geneSet` must be a valid geneSetList or geneSetFile.",
      call. = FALSE
    )
  }

  # check if tests are valid
  args[["test"]] <- unique(args[["test"]])
  if (!all(args[["test"]] %in% assocTest_aggregate_tests)) {
    stop(
      sprintf(
        "The following tests are not valid: %s",
        paste(
          args[["test"]][!args[["test"]] %in% assocTest_aggregate_tests],
          collapse = ","
        )
      ),
      call. = FALSE
    )
  }

  if (
    args[["continuous"]] && !all(args[["test"]] %in% assocTest_aggdb_cont_tests)
  ) {
    warning(
      "The following tests were excluded since they ",
      "are not available for continuous rvb tests: ",
      paste(
        args[["test"]][!args[["test"]] %in% assocTest_aggdb_cont_tests],
        collapse = ","
      ),
      call. = FALSE
    )
    args[["test"]] <- args[["test"]][
      args[["test"]] %in% assocTest_aggdb_cont_tests
    ]
  } else if (
    !args[["continuous"]] && !all(args[["test"]] %in% assocTest_aggdb_bin_tests)
  ) {
    warning(
      "The following tests were excluded since they are not ",
      "available for binary rvb tests: ",
      paste(
        args[["test"]][!args[["test"]] %in% assocTest_aggdb_bin_tests],
        collapse = ","
      ),
      call. = FALSE
    )
    args[["test"]] <- args[["test"]][
      args[["test"]] %in% assocTest_aggdb_bin_tests
    ]
  }

  # check if cohort exists in gdb
  cohort_tables <- listCohort(args[["gdb"]])$name
  if (!args[["cohort"]] %in% cohort_tables) {
    stop(
      sprintf(
        "The specified cohort ('%s') does not exist in the gdb.",
        args[["cohort"]]
      ),
      call. = FALSE
    )
  }

  # check if fields are present in cohort
  cohort_fields <- DBI::dbListFields(args[["gdb"]], args[["cohort"]])
  if (!args[["pheno"]] %in% cohort_fields) {
    stop(
      sprintf(
        "The `pheno` field '%s' was not found in cohort '%s'.",
        args[["pheno"]],
        args[["cohort"]]
      ),
      call. = FALSE
    )
  }

  covar_unlisted <- unlist(args[["covar"]])
  if (!is.null(args[["covar"]]) && !all(covar_unlisted %in% cohort_fields)) {
    missing_covars <- covar_unlisted[!covar_unlisted %in% cohort_fields]
    stop(
      sprintf(
        "The following `covar` fields were not found in cohort '%s': %s",
        args[["cohort"]],
        paste(sQuote(missing_covars), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (
    !is.null(args[["substractCovar"]]) &&
      !args[["substractCovar"]] %in% cohort_fields
  ) {
    stop(
      sprintf(
        "The `substractCovar` field '%s' was not found in cohort '%s'.",
        args[["substractCovar"]],
        args[["cohort"]]
      ),
      call. = FALSE
    )
  }

  # stop if no tests are left
  if (length(args[["test"]]) == 0L) {
    stop("No applicable tests specified.", call. = FALSE)
  }

  # return list of updated variables
  list(
    test = args[["test"]],
    covar = if (is.list(args[["covar"]])) args[["covar"]] else
      list(args[["covar"]])
  )
}

.assocTest_aggdb_init_output <- function(output) {
  if (is.null(output)) {
    return(NULL)
  }
  output_con <- gzcon(file(output, open = "wb"))
  write(
    paste(names(columns_aggresults), collapse = "\t"),
    file = output_con,
    append = FALSE
  )
  output_con
}

.assocTest_aggdb_load_cohort <- function(
  object,
  gdb,
  cohort,
  pheno,
  continuous,
  keep,
  verbose
) {
  cohort_df <- getCohort(gdb, cohort = cohort)

  # check if pheno fields are numeric
  pheno_is_numeric <- .check_pheno_numeric(cohort_df, pheno)
  if (!all(pheno_is_numeric)) {
    stop(
      sprintf(
        "The following phenotypes are not numeric: %s",
        paste(pheno[!pheno_is_numeric], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # check if binary phenotypes are coded 0,1
  if (!continuous) {
    pheno_is_binary <- .check_pheno_binary(cohort_df, pheno)
    if (!all(pheno_is_binary)) {
      stop(
        "Binary phenotypes should be coded 0,1.",
        "If the phenotype is continuous, set `continuous = TRUE`",
        call. = FALSE
      )
    }
  }

  # subset based on keep list if specified
  if (!is.null(keep)) {
    if (verbose) {
      message(sprintf(
        "Keeping %s/%s samples that are present in the keep-list.",
        sum(cohort_df[["IID"]] %in% keep),
        nrow(cohort_df)
      ))
    }

    cohort_df <- cohort_df[cohort_df[["IID"]] %in% keep, , drop = FALSE]
  }

  # early return if no samples left
  if (nrow(cohort_df) == 0L) {
    return(cohort_df)
  }

  # yield warning if not all sampples are present in aggdb
  if (!all(cohort_df$IID %in% listSamples(object))) {
    warning(
      sprintf(
        paste0(
          "%s/%s samples in the cohort are present ",
          "in the aggdb, subsetting the cohort."
        ),
        sum(cohort_df$IID %in% listSamples(object)),
        nrow(cohort_df)
      ),
      call. = FALSE
    )
    cohort_df <- cohort_df[
      cohort_df$IID %in% listSamples(object),
      ,
      drop = FALSE
    ]
  }

  if (!all(listSamples(object) %in% cohort_df$IID)) {
    warning(
      sprintf(
        paste0(
          "%s/%s samples in the aggdb are present ",
          "in the cohort, subsetting the aggdb"
        ),
        sum(listSamples(object) %in% cohort_df$IID),
        length(listSamples(object))
      ),
      call. = FALSE
    )
  }

  # return cohort
  cohort_df
}

.assocTest_aggdb_return_empty_results <- function() {
  fields <- lapply(columns_aggresults, function(type) {
    vector(mode = type[1L], length = 0L)
  })
  results <- as.data.frame(fields, stringsAsFactors = FALSE)
  results
}

.assocTest_aggdb_generate_params <- function(
  covar,
  pheno
) {
  params <- expand.grid(
    pheno = pheno,
    covar = covar,
    stringsAsFactors = FALSE
  )

  params
}

.assocTest_aggdb_run_geneset_analysis <- function(
  object,
  genesets,
  geneSet,
  units_aggdb,
  cohort_df,
  params,
  cohort_name,
  name,
  test,
  continuous,
  substractCovar,
  dropUnits,
  maxitFirth,
  output_con,
  verbose
) {
  # run analyses by geneset
  results <- lapply(genesets, FUN = function(geneset) {
    if (verbose) message(sprintf("Analysing %s", geneset))

    geneset_units <- .assocTest_aggdb_get_units(
      geneset = geneset,
      geneSet = geneSet,
      units_aggdb = units_aggdb,
      dropUnits = dropUnits,
      verbose = verbose
    )

    if (is.null(geneset_units)) {
      return(NULL)
    }

    # add aggregate to cohort
    cohort_df_agg <- .assocTest_aggdb_get_aggregates(
      object,
      cohort_df = cohort_df,
      units = geneset_units$units_keep,
      verbose = verbose
    )

    # run individual tasks (each combination of params)
    .assocTest_aggdb_run_geneset_tasks(
      cohort_df_agg = cohort_df_agg,
      cohort_name = cohort_name,
      geneset_units = geneset_units,
      params = params,
      geneset = geneset,
      name = name,
      test = test,
      continuous = continuous,
      substractCovar = substractCovar,
      maxitFirth = maxitFirth,
      output_con = output_con,
      verbose = verbose
    )
  })
  do.call(rbind, results)
}

.assocTest_aggdb_get_units <- function(
  geneset,
  geneSet,
  units_aggdb,
  dropUnits,
  verbose
) {
  # list units in geneset
  geneset_obj <- getGeneSet(geneSet, geneset)[[1L]]
  units_geneset <- listUnits(geneset_obj)

  # filter units based on overlap with aggdb and dropUnits (if specified)
  if (verbose) {
    message(sprintf(
      "   %s/%s units in the geneSet are present in the aggdb",
      sum(units_geneset %in% units_aggdb),
      length(units_geneset)
    ))
  }
  units_keep <- units_geneset[units_geneset %in% units_aggdb]

  # drop units if specified
  if (!is.null(dropUnits)) {
    units_keep <- units_keep[
      !units_keep %in% dropUnits
    ]
  }
  if (length(units_keep) == 0L) {
    if (verbose) {
      message(sprintf("No units left for %s, skipping", geneset))
    }
    return(NULL)
  }

  list(
    geneset_obj = geneset_obj,
    units_all = units_geneset,
    units_keep = units_keep
  )
}

.assocTest_aggdb_get_aggregates <- function(
  object,
  cohort_df,
  units,
  verbose
) {
  # get aggregates and collapse
  aggregates <- getUnit(object, unit = units)
  cohort_df$aggregate <- colSums(aggregates[,
    cohort_df[["IID"]],
    drop = FALSE
  ])

  # return cohort
  cohort_df
}

.assocTest_aggdb_run_geneset_tasks <- function(
  cohort_df_agg,
  cohort_name,
  geneset_units,
  params,
  geneset,
  name,
  test,
  continuous,
  substractCovar,
  maxitFirth,
  output_con,
  verbose
) {
  # run analyses per task (each combination of params)
  results <- lapply(
    seq_len(nrow(params)),
    FUN = function(i) {
      task <- params[i, ]
      task_covar <- unlist(task$covar)
      task_pheno <- task$pheno

      # add message here
      if (verbose) {
        message(sprintf(
          "   Analysing pheno: %s / covar: %s",
          task_pheno,
          paste(task_covar, collapse = ",")
        ))
      }

      # subset based on missing values in pheno,covar,substractCovar
      cohort_fields <- c(task_covar, task_pheno, substractCovar)
      task_cohort <- cohort_df_agg[
        complete.cases(cohort_df_agg[,
          cohort_fields,
          drop = FALSE
        ]),
      ]
      if (nrow(task_cohort) == 0L) {
        if (verbose) {
          message(sprintf("No units left for %s, skipping.", geneset))
        }
        return(NULL)
      } else if (verbose) {
        message(sprintf(
          "   %s/%s samples included with non-missing data for covariates and phenotype.",
          nrow(task_cohort),
          nrow(cohort_df_agg)
        ))
      }

      # handle covariates ------------
      covar_handle <- .handle_covar(
        coldata = task_cohort,
        covar = task_covar,
        pheno = task_pheno
      )
      task_covar <- covar_handle[["covar"]]
      task_cohort <- covar_handle[["coldata"]]

      if (!is.null(substractCovar)) {
        task_cohort[[substractCovar]] <- task_cohort[[substractCovar]] -
          task_cohort[["aggregate"]]
      }

      # generate models ---------
      model_handle <- .generate_models(
        covar = task_covar,
        pheno = task_pheno
      )
      null <- model_handle[["null"]]
      model <- model_handle[["model"]]
      model.nbinom <- model_handle[["model.nbinom"]]

      # init results -----------
      task_results <- .assocTest_aggdb_init_results(
        cohort = task_cohort,
        cohort_name = cohort_name,
        geneset = geneset,
        units_geneset = geneset_units$units_all,
        units_geneset_keep = geneset_units$units_keep,
        name = name,
        pheno = task_pheno,
        covar = task_covar,
        test = test,
        continuous = continuous
      )

      # generate stats
      if (continuous) {
        task_results <- .rvb_tests_rvb_cont_aggregate(
          cohort = task_cohort,
          results = task_results,
          test = test,
          pheno = task_pheno,
          model = model,
          null = null,
          covar = task_covar,
          output_con = output_con,
          append = TRUE,
          returnDF = TRUE
        )
      } else {
        task_results <- .rvb_tests_rvb_bin_aggregate(
          cohort = task_cohort,
          results = task_results,
          test = test,
          pheno = task_pheno,
          model = model,
          model.nbinom = model.nbinom,
          null = null,
          covar = task_covar,
          maxitFirth = maxitFirth,
          output_con = output_con,
          append = TRUE,
          returnDF = TRUE
        )
      }
      task_results
    }
  )
  do.call(rbind, results)
}


.assocTest_aggdb_init_results <- function(
  cohort,
  cohort_name,
  geneset,
  units_geneset,
  units_geneset_keep,
  name,
  pheno,
  covar,
  test,
  continuous
) {
  if (continuous) {
    caseN <- nrow(cohort)
    ctrlN <- 0
    meanCaseScore <- mean(cohort[, "aggregate"])
    meanCtrlScore <- NA_real_
  } else {
    caseN <- sum(cohort[, pheno] == 1, na.rm = TRUE)
    ctrlN <- sum(cohort[, pheno] == 0, na.rm = TRUE)
    meanCaseScore <- mean(cohort[cohort[, pheno] == 1, "aggregate"])
    meanCtrlScore <- mean(cohort[cohort[, pheno] == 0, "aggregate"])
  }

  results <- data.frame(
    geneSetName = geneset,
    cohort = cohort_name,
    name = name,
    pheno = pheno,
    covar = paste(covar, collapse = ","),
    test = test,
    geneSetSize = length(units_geneset),
    genesObs = length(units_geneset_keep),
    caseN = caseN,
    ctrlN = ctrlN,
    meanCaseScore = meanCaseScore,
    meanCtrlScore = meanCtrlScore,
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

.check_pheno_binary <- function(cohort, pheno) {
  check_binary <- vapply(
    pheno,
    FUN = function(phen) {
      all(cohort[[phen]] %in% c(0, 1), na.rm = TRUE)
    },
    FUN.VALUE = logical(1L),
    USE.NAMES = FALSE
  )
  check_binary
}


.check_pheno_numeric <- function(cohort, pheno) {
  check_numeric <- vapply(
    pheno,
    FUN = function(phen) {
      is.numeric(cohort[[phen]])
    },
    FUN.VALUE = logical(1L),
    USE.NAMES = FALSE
  )
  check_numeric
}


.rvb_tests_rvb_bin_aggregate <- function(
  cohort,
  results,
  test,
  pheno,
  model,
  model.nbinom,
  null,
  covar,
  maxitFirth = 1000L,
  output_con = NULL,
  append = FALSE,
  returnDF = FALSE
) {
  P <- OR <- effect <- effectSE <- effectCIupper <- effectCIlower <-
    rep(NA_real_, length(test))
  names(P) <- names(OR) <- names(effect) <- names(effectSE) <- names(
    effectCIupper
  ) <- names(effectCIlower) <-
    test

  if (
    sum(cohort[, pheno] == 1) < 2L ||
      sum(cohort[, pheno] == 0) < 2L ||
      sum(cohort$aggregate > 0) < 2L
  ) {
    test <- c()
  }

  if (sum(cohort$aggregate > 0.0) < 2L) {
    warning(
      "Less than two samples have a non-zero burden score, skipping tests.",
      call. = FALSE
    )
    test <- vector(mode = "character", length = 0L)
  }

  if (sum(cohort[, pheno] == 1) < 2L || sum(cohort[, pheno] == 0) < 2L) {
    warning("Fewer than two cases or controls, skipping tests.", call. = FALSE)
    test <- vector(mode = "character", length = 0L)
  }

  # Burden test
  if ("firth" %in% test) {
    tryCatch(
      {
        fit <- logistf::logistf(
          model,
          data = cohort,
          plconf = (which(c(covar, "aggregate") == "aggregate") + 1L),
          control = logistf::logistf.control(maxit = maxitFirth),
          plcontrol = logistf::logistpl.control(maxit = maxitFirth)
        )

        if (.check_conv_firth(fit, maxit = maxitFirth)) {
          effect["firth"] <- fit$coefficients["aggregate"]
          OR["firth"] <- exp(fit$coefficients["aggregate"])
          effectSE["firth"] <- sqrt(diag(vcov(fit)))["aggregate"]
          effectCIlower["firth"] <- fit$ci.lower["aggregate"]
          effectCIupper["firth"] <- fit$ci.upper["aggregate"]
          P["firth"] <- fit$prob["aggregate"]
        } else {
          effect["firth"] <- OR["firth"] <- effectSE["firth"] <- effectCIlower[
            "firth"
          ] <- effectCIupper["firth"] <- P["firth"] <- NA_real_
        }
      },
      error = function(e) {
        message(sprintf("Failed test '%s'\n%s", "firth", e))
      }
    )
  }

  if ("glm" %in% test) {
    tryCatch(
      {
        fit <- glm(model, data = cohort, family = "binomial")
        effect["glm"] <- summary(fit)$coef["aggregate", 1L]
        OR["glm"] <- exp(summary(fit)$coef["aggregate", 1L])
        effectSE["glm"] <- summary(fit)$coef["aggregate", 2L]
        effectCIlower["glm"] <- confint.default(fit)["aggregate", 1L]
        effectCIupper["glm"] <- confint.default(fit)["aggregate", 2L]
        P["glm"] <- summary(fit)$coef["aggregate", 4L]
      },
      error = function(e) {
        message(sprintf("Failed test '%s'\n%s", "glm", e))
      }
    )
  }

  if ("nbinom" %in% test) {
    tryCatch(
      {
        fit <- MASS::glm.nb(model.nbinom, data = cohort)
        effect["nbinom"] <- summary(fit)$coefficients[pheno, 1L]
        OR["nbinom"] <- NA_real_
        effectSE["nbinom"] <- summary(fit)$coefficients[pheno, 2L]
        effectCIlower["nbinom"] <- effect["nbinom"] -
          (1.96 * summary(fit)$coefficients[pheno, 2L])
        effectCIupper["nbinom"] <- effect["nbinom"] +
          (1.96 * summary(fit)$coefficients[pheno, 2L])
        P["nbinom"] <- summary(fit)$coefficients[pheno, 4L]
      },
      error = function(e) {
        message(sprintf("Failed test '%s'\n%s", "nbinom", e))
      }
    )
  }

  results$effect <- effect
  results$effectSE <- effectSE
  results$effectCIlower <- effectCIlower
  results$effectCIupper <- effectCIupper
  results$OR <- OR
  results$P <- P

  if (!is.null(output_con)) {
    write.table(
      results,
      sep = "\t",
      file = output_con,
      append = TRUE,
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
    return(results)
  } else {
    return(results)
  }
}


.rvb_tests_rvb_cont_aggregate <- function(
  GT,
  cohort,
  results,
  test,
  pheno,
  model,
  null,
  covar,
  output_con = NULL,
  append = FALSE,
  returnDF = FALSE
) {
  P <- OR <- effect <- effectSE <- effectCIupper <- effectCIlower <-
    rep(NA_real_, length(test))
  names(P) <- names(OR) <- names(effect) <- names(effectSE) <- names(
    effectCIupper
  ) <- names(effectCIlower) <-
    test

  if (sum(cohort[["aggregate"]] > 0) < 2) {
    warning(
      "Less than two samples have a non-zero burden score, skipping tests.",
      call. = FALSE
    )
    test <- c()
  }

  # Burden test
  if ("lm" %in% test) {
    tryCatch(
      {
        fit <- lm(model, data = cohort)
        effect["lm"] <- fit$coefficients["aggregate"]
        OR["lm"] <- NA_real_
        effectSE["lm"] <- summary(fit)$coef["aggregate", 2L]
        effectCIlower["lm"] <- confint(fit)["aggregate", 1L]
        effectCIupper["lm"] <- confint(fit)["aggregate", 2L]
        P["lm"] <- summary(fit)$coef["aggregate", 4L]
      },
      error = function(e) {
        message(sprintf("Failed test '%s'\n%s", "lm", e))
      }
    )
  }

  results$effect <- effect
  results$effectSE <- effectSE
  results$effectCIlower <- effectCIlower
  results$effectCIupper <- effectCIupper
  results$OR <- OR
  results$P <- P

  if (!is.null(output_con)) {
    write.table(
      results,
      file = output_con,
      sep = "\t",
      append = TRUE,
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE
    )
    return(results)
  } else {
    return(results)
  }
}
