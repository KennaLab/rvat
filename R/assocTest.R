#' assocTest-genoMatrix
#'
#' Run [`assocTest`] on a [`genoMatrix`] object. See the main [`assocTest`] page for details.
#'
#' @rdname assocTest-genoMatrix
#' @name assocTest-genoMatrix
#' @aliases assocTest,genoMatrix-method
#' @param object a [`genoMatrix`] object
#' @param pheno colData field to test as response variable, the response variable
#' can either be binary (0/1) or continuous. If the response variable is continuous set
#' `continuous` to `TRUE`.
#' @param test Vector of statistical tests to run,
#' options include firth,glm,lm,scoreSPA,skat,skat_burden,skato,skat_fwe,skat_burden_fwe,
#' skato_fwe,skat_robust,skato_robust,skat_burden_robust, acatv, acatvSPA. See [`assocTest`] for details.
#' @param name Optional name for the analysis, defaults to "none".
#' @param continuous Is the response variable continuous? (TRUE/FALSE). Defaults to `FALSE`.
#' @param singlevar Run single variant tests? (TRUE/FALSE).
#' Defaults to `FALSE`, in which case collapsing tests are ran.
#' @param covar Character vector of covariates. These should be present in the colData slot of the genoMatrix.
#' @param offset Optional model offset, can be used to account for regenie LOCO predictions.
#' @param overwriteAggregate In case there is already an `aggregate` column in the `colData` of the genoMatrix
#' (i.e. `aggregate` has been run on the genoMatrix), should it be overwitten? Defaults to `TRUE`.
#' @param geneticModel Which genetic model to apply? ('allelic', 'recessive' or 'dominant').
#' Defaults to `allelic`.
#' @param imputeMethod Which imputation method to apply? ('meanImpute' or 'missingToRef').
#' Defaults to `meanImpute`
#' @param MAFweights Apply MAF weighting? Currently Madsen-Browning ('mb') is implemented.
#' Defaults to 'none'.
#' @param maxitFirth Maximum number of iterations to use for estimating firth confidence intervals.
#' @param keep Vector of sample IDs to keep, defaults to `NULL`, in which case all samples are kept.
#' @param output Output file path for results.
#' Defaults to `NULL`, in which case results are not written to disk, but returned as an [`rvatResult`] object.
#' @param append Relevant if the `output` parameter is not `NULL`. Should results be appended to `output`?
#' @param returnDF Return a data.frame rather than a rvatResult. Defaults to `FALSE`.
#' @param methodResampling Which method to use for resampling? ('permutation' currently implemented)
#' Defaults to `NULL`, in which case no resampling is performed.
#' @param resamplingMatrix Pre-calculated resampling matrix (n x p), where n = number of samples, and p number of resamplings.
#' Can be generated using [`buildResamplingFile`].
#' @param resamplingFile A [`resamplingFile`] object.
#' @param nResampling Number of resamplings to perform if methodResampling is specified.
#' @param outputResampling If `TRUE` or a filepath, results for each resampling are returned (or saved to the filepath).
#' This can be useful if permutations are used to calculated to estimate correlations among genes for example.
#' Defaults to `FALSE` in which case resampling is used to calculate resampled P-values,
#' results for individual resamplings are not returned.
#' @param memlimitResampling Maximum number of resamplings to perform at a time.
#' Resampling generates a matrix of n x p, where n is the number of samples and p the number of resamplings
#' thus, for large number of resamplings it can be more efficient to split the permutations in chunks of size `memlimitResampling`.
#' Defaults to `NULL` in which case all permutations are performed.
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
#' @param verbose Should the function be verbose? (TRUE/FALSE), defaults to `TRUE`.
#' @examples
#'
#' library(rvatData)
#' data(GTsmall)
#'
#' # run a firth burden test on a binary phenotype
#' rvb <- assocTest(GTsmall,
#'                  pheno = "pheno",
#'                  covar = c("PC1", "PC2", "PC3", "PC4"),
#'                  test = "firth",
#'                  name = "example")
#'
#' # run ACAT-v and SKAT tests
#' rvb <- assocTest(GTsmall,
#'                  pheno = "pheno",
#'                  covar = c("PC1", "PC2", "PC3", "PC4"),
#'                  test = c("acatvfirth", "skat_burden_robust", "skato_robust"),
#'                  name = "example")
#'
#' # run a burden test on a continuous phenotype
#' rvb <- assocTest(GTsmall,
#'                  pheno = "age",
#'                  continuous = TRUE,
#'                  covar = c("PC1", "PC2", "PC3", "PC4"),
#'                  test = c("lm", "skat", "acatv"),
#'                  name = "example")
#'
#' # run single variant tests on a binary phenotype
#' sv <- assocTest(GTsmall,
#'                 pheno = "pheno",
#'                 singlevar = TRUE,
#'                 covar = c("PC1", "PC2", "PC3", "PC4"),
#'                 test = c("firth", "glm", "scoreSPA"),
#'                 name = "example",
#'                 minCarriers = 1
#'                 )
#'
#' # apply variant filters
#' rvb <- assocTest(GTsmall,
#'                  pheno = "pheno",
#'                  covar = c("PC1", "PC2", "PC3", "PC4"),
#'                  test = c("firth", "skat_robust", "acatv"),
#'                  name = "example",
#'                  maxMAF = 0.001,
#'                  minCarriers = 2,
#'                  minCallrateVar = 0.9,
#'                  minCallrateSM = 0.95)
#'
#' # Perform MAF-weighted burden tests (madsen-browning)
#' rvb <- assocTest(GTsmall,
#'                  pheno = "pheno",
#'                  covar = c("PC1", "PC2", "PC3", "PC4"),
#'                  test = c("firth", "skat_robust", "acatv"),
#'                  MAFweights = "mb")
#'
#' # Perform weighted burden tests with custom weights
#' gdb <- gdb(rvat_example("rvatData.gdb"))
#' # add cadd scores
#' caddscores <- getAnno(gdb, "varInfo", VAR_id = rownames(GTsmall))
#' caddscores <- caddscores[match(rownames(GTsmall), caddscores$VAR_id),]
#' # CADD-weighted burden test
#' rvb <- assocTest(recode(GTsmall, weights = caddscores$CADDphred),
#'                  pheno = "pheno",
#'                  covar = c("PC1", "PC2", "PC3", "PC4"),
#'                  test = c("firth", "skat_robust", "acatv"))
#'
#' # Perform recessive burden test
#' rvb <- assocTest(GTsmall,
#'                  pheno = "pheno",
#'                  covar = c("PC1", "PC2", "PC3", "PC4"),
#'                  test = c("firth", "skat_robust", "acatv"),
#'                  geneticModel = "recessive")
#'
#' # Resampled burden test
#' rvb <- assocTest(GTsmall,
#'                  pheno = "pheno",
#'                  covar = c("PC1", "PC2", "PC3", "PC4"),
#'                  test = c("skat", "skat_burden", "acatv"),
#'                  name = "example",
#'                  methodResampling = "permutation",
#'                  nResampling = 100)
#'
#' @export
setMethod(
  "assocTest",
  signature = signature(object = "genoMatrix"),
  definition = function(
    object,
    pheno,
    test,
    name = "none",
    continuous = FALSE,
    singlevar = FALSE,
    covar = NULL,
    offset = NULL,
    geneticModel = "allelic",
    imputeMethod = NULL,
    MAFweights = "none",
    maxitFirth = 1000L,
    keep = NULL,
    output = NULL,
    append = FALSE,
    returnDF = FALSE,
    methodResampling = NULL,
    resamplingMatrix = NULL,
    resamplingFile = NULL,
    nResampling = 1000L,
    outputResampling = FALSE,
    memlimitResampling = NULL,
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
    verbose = TRUE
  ) {
    # validity checks
    arg_parsed <- .assoctest_validate_input(as.list(environment()))
    test <- arg_parsed[["test"]]
    test_resampling <- arg_parsed[["test_resampling"]]
    geneticModel <- arg_parsed[["geneticModel"]]
    imputeMethod <- arg_parsed[["imputeMethod"]]
    memlimitResampling <- arg_parsed[["memlimitResampling"]]
    rm(arg_parsed)

    # filter samples based on callrate and/or keep list
    object <- .filterGT_samples(
      object = object,
      minCallrateSM = minCallrateSM,
      maxCallrateSM = maxCallrateSM,
      keep = keep,
      verbose = FALSE,
      verboseKeep = if (verbose) TRUE else FALSE
    )

    # generate list of parameter combinations
    params <- .assocTest_generate_params(pheno, covar, geneticModel, MAFweights)
    n_tasks <- nrow(params)

    if (verbose && nrow(params) > 1L) {
      message(sprintf(
        "Running %d parameter combinations..",
        nrow(params)
      ))
    }

    # generate sample masks for each param combination
    sample_masks <- lapply(seq_len(nrow(params)), function(i) {
      task <- params[i, ]
      cohort_fields <- c(task$pheno, unlist(task$covar))
      complete.cases(colData(object)[, cohort_fields, drop = FALSE])
    })

    # check which param combinations share identical sample lists
    mask_keys <- vapply(
      sample_masks,
      function(x) {
        paste(as.integer(x), collapse = "")
      },
      FUN.VALUE = character(1)
    )
    params$sample_mask_id <- match(mask_keys, unique(mask_keys))
    params$sample_mask <- I(sample_masks)
    params_by_sample_mask <- split(params, params$sample_mask_id)
    ## add task ids
    task_counter <- 1L
    params_by_sample_mask <- lapply(
      params_by_sample_mask,
      function(sample_mask) {
        sample_mask$task_id <- task_counter:(task_counter +
          nrow(sample_mask) -
          1L)
        task_counter <<- task_counter + nrow(sample_mask)
        sample_mask
      }
    )
    rm(task_counter)

    # generate results per sample mask
    results <- .assocTest_run_analyses_per_sample_mask(
      params_by_sample_mask = params_by_sample_mask,
      GT = object,
      n_tasks = n_tasks,
      test = test,
      name = name,
      continuous = continuous,
      singlevar = singlevar,
      offset = offset,
      append = append,
      imputeMethod = imputeMethod,
      maxitFirth = maxitFirth,
      keep = keep,
      output = output,
      returnDF = returnDF,
      methodResampling = methodResampling,
      resamplingMatrix = resamplingMatrix,
      resamplingFile = resamplingFile,
      nResampling = nResampling,
      outputResampling = outputResampling,
      memlimitResampling = memlimitResampling,
      minCallrateVar = minCallrateVar,
      maxCallrateVar = maxCallrateVar,
      minCallrateSM = minCallrateSM,
      maxCallrateSM = maxCallrateSM,
      minMAF = minMAF,
      maxMAF = maxMAF,
      minMAC = minMAC,
      maxMAC = maxMAC,
      minCarriers = minCarriers,
      maxCarriers = maxCarriers,
      minCarrierFreq = minCarrierFreq,
      maxCarrierFreq = maxCarrierFreq,
      test_resampling = test_resampling,
      verbose = verbose
    )
    rownames(results) <- NULL

    # return/write results
    .assoctest_return_results(
      object,
      results = results,
      output = output,
      returnDF = returnDF,
      singlevar = singlevar,
      outputResampling = outputResampling,
      append = append
    )
  }
)

.assocTest_run_analyses_per_sample_mask <- function(
  params_by_sample_mask,
  GT,
  n_tasks,
  test,
  name,
  continuous,
  singlevar,
  offset,
  append,
  imputeMethod,
  maxitFirth,
  keep,
  output,
  returnDF,
  methodResampling,
  resamplingMatrix,
  resamplingFile,
  nResampling,
  outputResampling,
  memlimitResampling,
  minCallrateVar,
  maxCallrateVar,
  minCallrateSM,
  maxCallrateSM,
  minMAF,
  maxMAF,
  minMAC,
  maxMAC,
  minCarriers,
  maxCarriers,
  minCarrierFreq,
  maxCarrierFreq,
  test_resampling,
  verbose
) {
  results <- lapply(
    params_by_sample_mask,
    function(task_sample_mask) {
      sample_mask <- task_sample_mask$sample_mask[[1L]]

      # sample filtering
      GT <- GT[, sample_mask]

      # subset samples and perform early return if no samples left
      if (ncol(GT) == 0L) {
        return(.return_empty_results(
          singlevar = singlevar,
          returnDF = TRUE
        ))
      }

      # flip to minor allele after sample filtering
      if (metadata(GT)$geneticModel == "allelic") {
        GT <- flipToMinor(GT)
      }

      # variant filtering
      keepGeno <- .filterGT_vars(
        GT,
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
        returnGT = FALSE,
        filterWeights = TRUE,
        verbose = FALSE
      )

      # subset variants and perform early return if no variants left
      if (sum(keepGeno) == 0L) {
        return(.return_empty_results(
          singlevar = singlevar,
          returnDF = TRUE
        ))
      }
      GT <- GT[keepGeno, ]

      # split tasks by geneticModel + MAFweights
      task_sample_mask$model_id <- paste(
        task_sample_mask$geneticModel,
        task_sample_mask$MAFweights,
        sep = " | "
      )
      params_by_model <- split(
        task_sample_mask,
        task_sample_mask$model_id
      )
      ## make sure task_ids are sequential
      task_counter <- min(task_sample_mask$task_id)
      params_by_model <- lapply(params_by_model, function(x) {
        x$task_id <- task_counter:(task_counter + nrow(x) - 1L)
        task_counter <<- task_counter + nrow(x)
        x
      })
      rm(task_counter)

      # loop through tasks per model
      results <- .assocTest_run_tasks_per_model(
        params_by_model = params_by_model,
        n_tasks = n_tasks,
        GT = GT,
        test = test,
        name = name,
        continuous = continuous,
        singlevar = singlevar,
        offset = offset,
        imputeMethod = imputeMethod,
        test_resampling = test_resampling,
        output = output,
        append = append,
        minCarriers = minCarriers,
        maxitFirth = maxitFirth,
        returnDF = returnDF,
        methodResampling = methodResampling,
        resamplingMatrix = resamplingMatrix,
        resamplingFile = resamplingFile,
        nResampling = nResampling,
        outputResampling = outputResampling,
        memlimitResampling = memlimitResampling,
        verbose = verbose
      )
      results
    }
  )
  do.call(rbind, results)
}

.assocTest_run_tasks_per_model <- function(
  params_by_model,
  n_tasks,
  GT,
  test,
  name,
  continuous,
  singlevar,
  offset,
  imputeMethod,
  test_resampling,
  output,
  append,
  minCarriers,
  maxitFirth,
  returnDF,
  methodResampling,
  resamplingMatrix,
  resamplingFile,
  nResampling,
  outputResampling,
  memlimitResampling,
  verbose
) {
  results <- lapply(seq_along(params_by_model), function(i) {
    task_model <- params_by_model[[i]]
    task_geneticModel <- task_model$geneticModel[1L]
    task_MAFweights <- task_model$MAFweights[1L]

    # single variant tests
    if (singlevar) {
      # calculate MAFs before applying geneticModel
      phenos <- unique(task_model$pheno)
      MAF_per_pheno <- lapply(
        phenos,
        function(task_pheno) {
          .assocTest_calc_maf_per_pheno(
            GT = GT,
            pheno = task_pheno,
            continuous = continuous
          )
        }
      )
      names(MAF_per_pheno) <- phenos

      # calculate sample numbers before imputation
      N_per_pheno <- lapply(
        phenos,
        function(task_pheno) {
          .assocTest_calc_n_per_pheno(
            GT = GT,
            pheno = task_pheno,
            continuous = continuous
          )
        }
      )
      names(N_per_pheno) <- phenos

      # apply genetic model
      if (task_geneticModel != metadata(GT)$geneticModel) {
        GT <- recode(GT, geneticModel = task_geneticModel)
      }

      # apply imputation, if specified
      GT <- recode(GT, imputeMethod = imputeMethod)

      # loop through tasks per covar/pheno combination
      results_per_covar_pheno <- .assocTest_run_task_per_covar_pheno_sv(
        GT,
        task_model = task_model,
        n_tasks = n_tasks,
        offset = offset,
        MAF_per_pheno = MAF_per_pheno,
        N_per_pheno = N_per_pheno,
        name = name,
        test = test,
        continuous = continuous,
        output = output,
        append = append,
        maxitFirth = maxitFirth,
        returnDF = returnDF,
        verbose = verbose
      )
      
      # rvb tests
    } else {
      
      # generate MAFweights, if specified
      if (task_MAFweights != "none") {
        GT <- recode(GT, MAFweights = task_MAFweights)
      }

      # apply genetic model
      if (task_geneticModel != metadata(GT)$geneticModel) {
        GT <- recode(GT, geneticModel = task_geneticModel)
      }
      
      # get callrate and carriers before imputation
      carriers <- as.integer(
        Matrix::colSums(assays(GT)$GT >= 1, na.rm = TRUE) >= 1
      )

      # drop variants with 0 carriers
      if (minCarriers == 0) {
        carriers_var <- getNCarriers(GT)
        if (any(carriers_var == 0)) {
          if (verbose) {
            message(sprintf(
              "  %s/%s variants have zero carriers, these are dropped.",
              sum(carriers_var == 0),
              nrow(GT)
            ))
          }
          GT <- GT[carriers_var >= 1, ]

          # return early if no variants left
          if (nrow(GT) == 0L) {
            return(.return_empty_results(
              singlevar = singlevar,
              returnDF = TRUE
            ))
          }
        }
      }

      callRate <- getCR(GT, var = FALSE)

      # apply imputation
      GT <- recode(
        GT,
        imputeMethod = imputeMethod
      )

      # aggregate
      GT <- aggregate(
        GT,
        checkMissing = FALSE
      )

      # generate results per covar/pheno combination
      results_per_covar_pheno <- .assocTest_run_task_per_covar_pheno_rvb(
        task_model = task_model,
        n_tasks = n_tasks,
        GT = GT,
        carriers = carriers,
        callRate = callRate,
        offset = offset,
        name = name,
        test = test,
        test_resampling = test_resampling,
        continuous = continuous,
        output = output,
        append = append,
        nResampling = nResampling,
        maxitFirth = maxitFirth,
        returnDF = returnDF,
        resamplingMatrix = resamplingMatrix,
        resamplingFile = resamplingFile,
        methodResampling = methodResampling,
        memlimitResampling = memlimitResampling,
        outputResampling = outputResampling,
        verbose = verbose
      )
      results_per_covar_pheno
    }
  })
  do.call(rbind, results)
}

.assocTest_run_task_per_covar_pheno_rvb <- function(
  task_model,
  n_tasks,
  GT,
  offset,
  carriers,
  callRate,
  name,
  test,
  test_resampling,
  continuous,
  output,
  append,
  nResampling,
  maxitFirth,
  returnDF,
  resamplingMatrix,
  resamplingFile,
  methodResampling,
  memlimitResampling,
  outputResampling,
  verbose
) {
  results <- lapply(seq_len(nrow(task_model)), FUN = function(i) {
    task <- task_model[i, ]
    task_covar <- unlist(task$covar)
    task_pheno <- task$pheno
    task_MAFweights <- task$MAFweights
    task_geneticModel <- task$geneticModel

    # message summarizing current task
    if (verbose) {
      str_mafweights <- if (task$MAFweights == "none") {
        ""
      } else {
        sprintf(" (%s-weighted)", task$MAFweights)
      }

      str_covar <- if (length(task_covar) == 0L) {
        "no covariates"
      } else if (length(task_covar) > 6) {
        sprintf(
          "%s + ... (%d covariate(s))",
          paste(task_covar[1:6], collapse = ","),
          length(task_covar) - 6
        )
      } else {
        paste(task_covar, collapse = ",")
      }

      message(sprintf(
        "[%d/%d] (%d samples, %d variants) %s | %s | %s%s",
        task$task_id,
        n_tasks,
        ncol(GT),
        nrow(GT),
        task_pheno,
        str_covar,
        task_geneticModel,
        str_mafweights
      ))
    }

    # handle covariates (drop covar with zero variants/duplicates + generate dummies)
    tmp <- .handle_covar(
      coldata = colData(GT),
      covar = task_covar,
      pheno = task_pheno
    )
    task_covar <- tmp[["covar"]]
    colData(GT) <- tmp[["coldata"]]
    rm(tmp)

    # generate models
    tmp <- .generate_models(
      covar = task_covar,
      offset = offset,
      pheno = task_pheno
    )
    null <- tmp[["null"]]
    model <- tmp[["model"]]
    model.nbinom <- tmp[["model.nbinom"]]
    rm(tmp)

    # initialize results
    results <- .assoctest_init_results_rvb(
      GT,
      carriers = carriers,
      callRate = callRate,
      name = name,
      pheno = task_pheno,
      covar = task_covar,
      geneticModel = task_geneticModel,
      MAFweights = task_MAFweights,
      test = test,
      continuous = continuous
    )

    # run rvb tests
    res <- .rvb_tests_rvb(
      GT = GT,
      results = results,
      test = test,
      pheno = task_pheno,
      model = model,
      model.nbinom = if (continuous) NULL else model.nbinom,
      null = null,
      covar = task_covar,
      continuous = continuous,
      output = output,
      append = append,
      nResampling = nResampling,
      maxitFirth = maxitFirth,
      returnDF = returnDF
    )

    # resampled rvb tests (if specified)
    if (!is.null(resamplingMatrix) || !is.null(resamplingFile) || !is.null(methodResampling)) {
      res <- .permute_rvb_wrapper(
        GT,
        res = res,
        init_results = results,
        test = test,
        test_resampling = test_resampling,
        pheno = task_pheno,
        covar = task_covar,
        model = model,
        null = null,
        continuous = continuous,
        resamplingMatrix = resamplingMatrix,
        resamplingFile = resamplingFile,
        nResampling = nResampling,
        methodResampling = methodResampling,
        memlimitResampling = memlimitResampling,
        outputResampling = outputResampling
      )
    }
    res
  })
  do.call(rbind, results)
}

.assocTest_run_task_per_covar_pheno_sv <- function(
  GT,
  task_model,
  n_tasks,
  offset,
  MAF_per_pheno,
  N_per_pheno,
  name,
  test,
  continuous,
  output,
  append,
  maxitFirth,
  returnDF,
  verbose
) {
  results <- lapply(seq_len(nrow(task_model)), FUN = function(i) {
    task <- task_model[i, ]
    task_covar <- unlist(task$covar)
    task_pheno <- task$pheno
    caseMAF <- MAF_per_pheno[[task_pheno]]$caseMAF
    ctrlMAF <- MAF_per_pheno[[task_pheno]]$ctrlMAF
    caseN <- N_per_pheno[[task_pheno]]$caseN
    ctrlN <- N_per_pheno[[task_pheno]]$ctrlN

    # message
    if (verbose) {
      str_mafweights <- if (task$MAFweights == "none") {
        ""
      } else {
        sprintf(" (%s-weighted)", task$MAFweights)
      }

      str_covar <- if (length(task_covar) == 0L) {
        "no covariates"
      } else if (length(task_covar) > 6) {
        sprintf(
          "%s + ... (%d covariate(s))",
          paste(task_covar[1:6], collapse = ","),
          length(task_covar) - 6
        )
      } else {
        paste(task_covar, collapse = ",")
      }

      message(sprintf(
        "[%d/%d] (%d samples, %d variants) %s | %s | %s%s",
        task$task_id,
        n_tasks,
        ncol(GT),
        nrow(GT),
        task_pheno,
        str_covar,
        task$geneticModel,
        str_mafweights
      ))
    }

    # handle covariates (drop covar with zero variants/duplicates + generate dummies)
    tmp <- .handle_covar(
      coldata = colData(GT),
      covar = task_covar,
      pheno = task_pheno
    )
    task_covar <- tmp[["covar"]]
    coldata_original <- colData(GT)
    colData(GT) <- tmp[["coldata"]]
    rm(tmp)

    # generate models
    tmp <- .generate_models(
      covar = task_covar,
      offset = offset,
      pheno = task_pheno
    )
    null <- tmp[["null"]]
    model <- tmp[["model"]]
    model.nbinom <- tmp[["model.nbinom"]]
    rm(tmp)

    # initalize results
    results <- .assoctest_init_results_sv(
      GT,
      name = name,
      pheno = task_pheno,
      covar = task_covar,
      test = test,
      continuous = continuous,
      caseMAF = caseMAF,
      ctrlMAF = ctrlMAF,
      caseN = caseN,
      ctrlN = ctrlN
    )

    # run singlevar tests
    res <- .rvb_tests_singlevar(
      GT = GT,
      results = results,
      pheno = task_pheno,
      test = test,
      model = model,
      model.nbinom = if (continuous) NULL else model.nbinom,
      null = null,
      covar = task_covar,
      continuous = continuous,
      output = output,
      append = append,
      returnDF = returnDF,
      maxitFirth = maxitFirth,
      verbose = verbose
    )
    
    # set original colData (i.e. before dummy coding)
    colData(GT) <- coldata_original

    # return
    res
  })
  do.call(rbind, results)
}