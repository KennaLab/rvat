#' assocTest-gdb
#' 
#' Run [`assocTest`] on a [`gdb`] object. See the main [`assocTest`] page for details.
#' 
#' @rdname assocTest-gdb
#' @name assocTest-gdb
#' @aliases assocTest,gdb-method
#' @param object a [`gdb`] object
#' @param pheno colData field to test as response variable, the response variable
#' can either be binary (0/1) or continuous. 
#' If the response variable is continuous set `continuous` to `TRUE`.
#' Multiple phenotypes can be specified, which will then be tested separately.
#' @param test Vector of statistical tests to run,
#' options include firth,glm,lm,nbinom,skat,skat_burden,skato,skat_fwe,skat_burden_fwe,
#' skato_fwe,skat_robust,skato_robust,skat_burden_robust, acatv, acatvSPA. 
#' See [`assocTest`] for details.
#' @param cohort If a valid cohort name is provided, then the uploaded data for 
#' this cohort is used to filter and annotate the genoMatrix object. 
#' If not specified, all samples in the gdb will be loaded.
#' @param varSet a [`varSetFile`] or [`varSetList`] object. 
#' Alternatively a vector of VAR_ids can be specified using the `VAR_id` parameter.
#' @param VAR_id a vector of VAR_ids, alternatively the `varSet` parameter can be specified.
#' If single variant tests are run, the `memlimit` argument controls how many variants to analyze at a time.
#' @param name Optional name for the analysis, defaults to "none".
#' @param continuous Is the response variable continuous? (TRUE/FALSE). Defaults to `FALSE`.
#' @param singlevar Run single variant tests? (TRUE/FALSE).
#' Defaults to `FALSE`, in which case aggregate tests are run.
#' @param covar Character vector of covariates, or a list of character vectors of covariates 
#' in which case each covariate set will be tested separately.
#' @param offset Optional model offset, can be used to account for regenie LOCO predictions.
#' @param geneticModel Which genetic model to apply? ('allelic', 'recessive' or 'dominant').
#' Defaults to `allelic`.
#' Multiple geneticModels can be specified, in which case each will be analyzed separately.
#' @param imputeMethod Which imputation method to apply? ('meanImpute' or 'missingToRef').
#' Defaults to `meanImpute`.
#' @param MAFweights MAF weighting method. Currently Madsen-Browning ('mb') is implemented, 
#' by default no MAF weighting is applied.
#' Multiple MAFweights can be specified, in which case each will be analyzed separately.
#' @param maxitFirth Maximum number of iterations to use for estimating firth confidence intervals. 
#' Defaults to 1000.
#' @param checkPloidy Version of the human genome to use when assigning variant ploidy (diploid, XnonPAR, YnonPAR). 
#' Accepted inputs are GRCh37, hg19, GRCh38, hg38.
#' If not specified, the genome build in the [`gdb`] will be used, if available 
#' (included if the `genomeBuild` parameter was set in [`buildGdb`]).
#' Otherwise, if the genome build is not included in the gdb metadata, 
#' and no value is provided, then all variants are assigned the default ploidy of "diploid"
#' @param keep Vector of sample IDs to keep, defaults to `NULL`, 
#' in which case all samples are kept.
#' @param output Output file path for results.
#' Defaults to `NULL`, in which case results are not written.
#' @param append Relevant if the `output` parameter is not `NULL`. 
#' Should results be appended to `output`?
#' @param returnDF Return a data.frame rather than a rvatResult. Defaults to `FALSE`.
#' @param methodResampling Which method to use for resampling? ('permutation' currently implemented)
#' Defaults to `NULL`, in which case no resampling is performed. 
#' @param resamplingMatrix Pre-calculated resampling matrix (n x p), 
#' where n = number of samples, and p number of resamplings.
#' Can be generated using [`buildResamplingFile`].
#' @param resamplingFile A [`resamplingFile`] object.
#' @param nResampling Number of resamplings to perform if methodResampling is specified.
#' @param outputResampling If `TRUE` or a filepath, 
#' results for each resampling are returned (or saved to the filepath).
#' This can be useful if permutations are used to estimate correlations among genes for example.
#' Defaults to `FALSE` in which case resampling is used to calculate resampled P-values, 
#' results for individual resamplings are not returned.
#' @param memlimitResampling Maximum number of resamplings to perform at a time.
#' Resampling generates a matrix of n x p, where n is the number of samples and p the number of resamplings
#' thus, for a large number of resamplings it can be more efficient to split the permutations in chunks of size `memlimitResampling`.
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
#' @param memlimit Maximum number of variants to load at once (if `VAR_id` is specified).
#' @param verbose Should the function be verbose? (TRUE/FALSE), defaults to `TRUE`.
#' @param strict Should strict checks be performed? Defaults to `TRUE`. Strict tests currently includes
#' checking whether supplied varSetFile/varSetList was generated from the same gdb as specified in `object`.
#' @example inst/examples/example-assocTest-gdb.R
#'
#' @export
#'
setMethod(
  "assocTest",
  signature = signature(object = "gdb"),
  definition = function(
    object,
    pheno,
    test,
    cohort = "SM",
    varSet = NULL,
    VAR_id = NULL,
    name = "none",
    continuous = FALSE,
    singlevar = FALSE,
    covar = NULL,
    offset = NULL,
    geneticModel = "allelic",
    imputeMethod = NULL,
    MAFweights = "none",
    maxitFirth = 1000L,
    checkPloidy = NULL,
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
    memlimit = 1000L,
    verbose = TRUE,
    strict = TRUE
  ) {
    # validity checks
    arg <- as.list(environment())
    .assoctest_gdb_validate_input(arg)
    if (!is.null(varSet) && strict) {
      .check_gdb_ids(object, varSet, minVersion = "0.3.0")
    }

    # if VAR_id is specified, generate varSet from VAR_ids
    if (!is.null(VAR_id)) {
      varSet <- .varsTovarSetList(VAR_id, chunkSize = memlimit)
    }

    # pre-load cohort
    chrt <- getCohort(
      object,
      cohort = cohort,
      fields = unique(c("IID", "sex", pheno, unlist(covar), offset))
    )
    chrt <- DataFrame(chrt[!is.na(chrt$IID), ])
    metadata(chrt)$name <- cohort

    # run
    units_all <- unique(listUnits(varSet))
    results <- lapply(units_all, FUN = function(unit) {
      if (verbose) {
        message(sprintf("Analysing %s -------------------", unit))
      }

      # load genotypes
      varsets <- getVarSet(varSet, unit = unit)
      if (singlevar) {
        GT <- getGT(
          object,
          cohort = chrt,
          varSet = if (length(varsets) > 1) {
            collapseVarSetList(varsets)
          } else {
            varsets
          },
          checkPloidy = checkPloidy,
          anno = "var",
          annoFields = c("VAR_id", "REF", "ALT"),
          verbose = verbose,
          strict = FALSE
        )
        colnames(rowData(GT))[colnames(rowData(GT)) == "REF"] <- "otherAllele"
        colnames(rowData(GT))[colnames(rowData(GT)) == "ALT"] <- "effectAllele"
      } else {
        GT <- getGT(
          object,
          cohort = chrt,
          varSet = if (length(varsets) > 1) {
            collapseVarSetList(varsets)
          } else {
            varsets
          },
          checkPloidy = checkPloidy,
          verbose = verbose,
          strict = FALSE
        )
      }
      metadata(GT)$unit <- unit
      metadata(GT)$cohort <- cohort

      # loop through varsets
      varsets_names <- listVarSets(varsets)

      results <- lapply(varsets_names, FUN = function(vs) {
        if (verbose) {
          message(sprintf("varSet: %s -------------------", vs))
        }
        varset <- getVarSet(varsets, varSetName = vs)[[1]]
        vars <- listVars(varset)
        weights <- listWeights(varset)

        metadata(GT)$varSetName <- if(is.null(VAR_id)) vs else "none"
        results <- assocTest(
          recode(GT[vars, ], weights = weights),
          pheno = pheno,
          test = test,
          name = name,
          continuous = continuous,
          singlevar = singlevar,
          covar = covar,
          offset = offset,
          geneticModel = geneticModel,
          imputeMethod = imputeMethod,
          MAFweights = MAFweights,
          maxitFirth = maxitFirth,
          keep = keep,
          output = output,
          append = append,
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
          verbose = verbose
        )
        results
      })
      results <- do.call(rbind, results)
      results
    })
    results <- do.call(rbind, results)
    results
  }
)


.assoctest_gdb_validate_input <- function(args) {
  # type checks ------------------------------------------------
  .gdb_check_varid(args[["VAR_id"]])
  check_wrapper(check_character, args, "cohort", length_equal = 1L)
  check_wrapper(check_bool, args, "strict")

  if (is.null(args[["varSet"]]) && is.null(args[["VAR_id"]])) {
    stop("Either of `varSet` or `VAR_id` should be specified", call. = FALSE)
  }

  if (!is.null(args[["varSet"]]) && !is.null(args[["VAR_id"]])) {
    stop(
      "Either of one of `varSet` or `VAR_id` should be specified, not both.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
