#' assocTest-genoMatrix
#' 
#' Run [`assocTest`] on a [`genoMatrix`] object. See the main [`assocTest`] page for details.
#' 
#' @rdname assocTest-genoMatrix
#' @name assocTest-genoMatrix
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
#' @param overwriteAggregate In case there is already an `aggregate` column in the `colData` of the genoMatrix 
#' (i.e. `aggregate` has been run on the genoMatrix), should it be overwitten? Defaults to `TRUE`.
#' @param geneticModel Which genetic model to apply? ('allelic', 'recessive' or 'dominant').
#' Defaults to `allelic`. 
#' @param imputeMethod Which imputation method to apply? ('meanImpute' or 'missingToRef').
#' Defaults to `meanImpute`
#' @param MAFweights Apply MAF weighting? Currently Madsen-Browning ('mb') is implemented.
#' Defaults to 'none'.
#' @param maxitFirth Maximum number of iterations to use for estimating firth confidence intervals.
#' @param output Output file path for results. 
#' Defaults to `NULL`, in which case results are not written to disk, but returned as an [`rvatResult`] object.
#' @param append Relevant if the `output` parameter is not `NULL`. Should results be appended to `output`?
#' Defaults to `FALSE`.
#' @param methodResampling Which method to use for resampling? ('permutation' currently implemented)
#' Defaults to `NULL`, in which case no resampling is performed. 
#' @param resamplingMatrix Pre-calculated resampling matrix (n x p), where n = number of samples, and p number of resamplings. 
#' Can be generated using [`buildResamplingFile`].
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
#' @export
setMethod("assocTest", 
          signature = signature(object="genoMatrix"),
          definition=function(object,
                              pheno,
                              test,
                              name = "none",
                              continuous = FALSE, 
                              singlevar = FALSE, 
                              covar = NULL,
                              offset = NULL,
                              overwriteAggregate = TRUE,
                              geneticModel = c("allelic", "recessive", "dominant"),
                              imputeMethod = NULL,
                              MAFweights = "none",
                              maxitFirth = 1000,
                              keep = NULL,
                              output = NULL,
                              append = FALSE,
                              returnDF = FALSE,
                              methodResampling = NULL,
                              resamplingMatrix = NULL,
                              nResampling = 1000,
                              outputResampling = FALSE,
                              memlimitResampling = NULL,
                              minCallrateVar = 0,
                              maxCallrateVar = Inf,
                              minCallrateSM = 0,
                              maxCallrateSM =Inf,
                              minMAF = 0,
                              maxMAF = 1,
                              minMAC = 0,
                              maxMAC = Inf,
                              minCarriers = 0,
                              maxCarriers = Inf,
                              minCarrierFreq = 0,
                              maxCarrierFreq = Inf,
                              verbose = TRUE
          )
          {
            # Validity checks -------------------------------------------------
            
            ## check: if overwriteAggregate=TRUE, does colData contain an 'aggregate' field?
            if (!overwriteAggregate && !"aggregate" %in% colnames(colData(object)) & !singlevar) {
              stop("`overwriteAggregate` is set to FALSE, but no 'aggregate' column is present in the colData of the genoMatrix.")
            }
            ## geneticModel
            geneticModel <- match.arg(geneticModel)
            
            ### if geneticModel is dominant or recessive, stop if user has specified MAF or MAC filters.
            if(S4Vectors::metadata(object)$geneticModel != "allelic" &
               (minMAF > 0 ||
                maxMAF < 1 ||
                minMAC > 0 ||
                maxMAC < Inf)) {
              stop("MAC and/or MAF filters do not apply when geneticModel != 'allelic'.")
            }
            
            if(S4Vectors::metadata(object)$geneticModel != "allelic" && geneticModel == "allelic") {
              warning(sprintf("The genoMatrix provided is set to '%s', setting the `geneticModel` parameter to '%s'",
                              S4Vectors::metadata(object)$geneticModel, S4Vectors::metadata(object)$geneticModel
              ))
              geneticModel <- S4Vectors::metadata(object)$geneticModel
            }
            
            ### not possible to set geneticModel == "allelic" if it's already set to 'dominant' or 'recessive'
            if(S4Vectors::metadata(object)$geneticModel != "allelic" && geneticModel != S4Vectors::metadata(object)$geneticModel) {
              stop(sprintf("Current geneticModel should be 'allelic' in order to apply dominant or recessive models."))
            }
            
            ## imputeMethod
            
            ### If imputeMethod = NULL, set to meanImpute if singlevar=FALSE
            ### + check if specified imputeMethod is available
            if(is.null(imputeMethod)) {
              if(!singlevar) imputeMethod <- "meanImpute"
            } else {
              if(!imputeMethod %in% c("meanImpute", "missingToRef")) 
                stop("`imputeMethod` should be either 'meanImpute' or 'missingToRef'")
            }
            
            ## check tests
            
            ### check if at least one of the specified tests is valid
            if(!all(test %in% assocTest_tests)) {
              stop(sprintf("The following tests are not valid: %s", 
                           paste(test[!test %in% assocTest_tests],
                                 collapse=",")
              ))
            }
            test <- unique(test)
            
            ### Check if correct tests are specified for binary/continuous traits and sv/rvb tests
            if(singlevar) {
              if(!continuous) {
            
                if(!all(test %in% assocTest_sv_bin_tests)) {
                  warning("The following tests were excluded since they are not available for binary singlevar tests: ",
                          paste(test[!test %in% assocTest_sv_bin_tests], collapse=","))
                  test <- test[test %in% assocTest_sv_bin_tests]
                }
              } else {
                if(!all(test %in% assocTest_sv_cont_tests)) {
                  warning("The following tests were excluded since they are not available for continuous singlevar tests: ",
                          paste(test[!test %in% assocTest_sv_cont_tests], collapse=","))
                  test <- test[test %in% assocTest_sv_cont_tests]
                }
              }
              test <- assocTest_sv_tests[assocTest_sv_tests %in% test]
            } else {
              if(!continuous) {
                
                if(!all(test %in% assocTest_rvb_bin_tests)) {
                  warning("The following tests were excluded since they are not available for binary rvb tests: ",
                          paste(test[!test %in% assocTest_rvb_bin_tests], collapse=","))
                  test <- test[test %in% assocTest_rvb_bin_tests]
                }
    
              } else {
                if(!all(test %in% assocTest_rvb_cont_tests)) {
                  warning("The following tests were excluded since they are not available for continuous rvb tests: ",
                          paste(test[!test %in% assocTest_rvb_cont_tests], collapse=","))
                  test <- test[test %in% assocTest_rvb_cont_tests]
                }
              }
            }
            
            ### offset
            if ( !is.null(offset) ) {
              if ( length(offset) > 1 ) {
                stop("Currently at most one variable can be specified as offset")
              }
              
              if ( !offset %in% colnames(SummarizedExperiment::colData(object)) ) {
                stop(sprintf("Offset variable: %s is not available ", offset))
              }
              
              if ( !all(test %in% assocTest_offset_tests)) {
                warning(sprintf("Including an offset is not implemented for the following tests: %s",
                                paste(test[!test %in% assocTest_offset_tests], collapse=",")))
              }
              test <- test[test %in% assocTest_offset_tests]
            }
            
            # Stop if no tests are left after filtering
            if(length(test) == 0) stop("No applicable tests specified.")
            
            ## permutation parameters
            
            ### currently resampling is only implemented for rvb tests
            if(singlevar && (!is.null(resamplingMatrix) || !is.null(methodResampling))) {
              stop("Resampling is currently not implemented for singlevar tests")
            }
            
            ### subset tests for which resampling is implemented 
            if((!is.null(resamplingMatrix) || !is.null(methodResampling))) {
              test_resampling <- test[test %in% assocTest_resampling_tests]
              if(length(test_resampling) == 0) {
                stop(sprintf("Resampling is not implemented for any of the specified tests, resampling is currently implemented for: %s",
                             paste(assocTest_resampling_tests, collapse=",")
                ))
              }
              
              if(length(test_resampling) < length(test)) {
                warning(sprintf("Resampling is not implemented for all of the specified tests, resampling is currently implemented for: %s",
                                paste(assocTest_resampling_tests, collapse=",")
                ))
              }
            }
            
            if(!(is.character(outputResampling) || is.logical(outputResampling))) {
              stop("`outputResampling` should either be a filepath or a boolean (TRUE/FALSE)")
            }
            
            if(nResampling < 1) {
              stop("`nResampling` should be a positive number")
            }
            
            if(!is.null(methodResampling) && !methodResampling %in% c("permutation")) {
              stop("Currently only 'permutation' is an accepted option for `methodResampling`")
            }
            
            if(is.null(memlimitResampling)) memlimitResampling <- nResampling
            
            
            ## phenotypes
            
            ### check if one phenotype is specified
            if( length(pheno) > 1) {stop("Only one phenotype can be specified.")}
            
            ### check if specified phenotype are present in colData
            if(!pheno %in% colnames(colData(object))) stop(sprintf("'%s' is not present in `colData(GT)`", pheno))
            
            ### keep samples for which phenotype is non-missing
            if (!is.null(keep)) {
              keepSamples <- colnames(object) %in% keep
              if(verbose) message(sprintf("%s/%s samples present in `keep` are kept.", sum(keepSamples), ncol(object)))
            } else {
              keepSamples <- rep(TRUE, ncol(object))
            }
            kp <- !is.na(colData(object)[,pheno])
            if (!all(kp)) {
              if(verbose) message(sprintf("%s samples have missing phenotype values, these will be dropped." , sum(!kp)))
              keepSamples <- keepSamples & kp
            }
            
            ### check if binary phenotype is coded as 0,1
            if (!continuous) {
              if(!all(colData(object)[,pheno][!is.na(colData(object)[,pheno])] %in% c(0,1))) stop("Binary phenotypes should be coded 0,1! If the phenotype is continuous, set `continuous = TRUE")
            }
            
            ## Add precomputed score (if specified)
            if(overwriteAggregate) {
              colData(object)$aggregate <- NA_real_
            } 
            
            ## covariates
            
            ### check if covariates are present in colData
            if(!is.null(covar) && !(length(covar) == 1 && covar == "1")) {
              if (length(unique(covar)) < length(covar))
                warning("Duplicate covariates are specified, unique covariates are kept.")
              covar <- unique(covar)
              if (mean(covar %in% colnames(colData(object))) < 1) {
                stop(sprintf("The following covariate(s) are not available: %s",
                             paste(covar[!covar %in% colnames(colData(object))], collapse=",")))
              }
              
              ### Handle character/factor covar fields 
              covar_types <- unlist(lapply(colData(object)[,covar,drop=FALSE], FUN = class))
              if(sum(covar_types %in% c("character", "factor")) > 0) {
                modmatrix <- as.data.frame(colData(object))
                modmatrix[names(covar_types[covar_types=="character"])] <- lapply( modmatrix[names(covar_types[covar_types=="character"])],
                                                                                   factor)
                modmatrix <-  model.matrix.lm(as.formula(sprintf("%s ~ %s", pheno, paste(covar, collapse="+"))), 
                                              data = modmatrix, na.action="na.pass")
                modmatrix <- modmatrix[,-1]
                colnames(modmatrix) <- make.names(colnames(modmatrix))
                covar <- colnames(modmatrix)
                colData(object) <- cbind(
                  colData(object)[,setdiff(colnames(colData(object)), colnames(modmatrix))],modmatrix)
              }
              
              ### drop covariates with zero variance.
              covarKeep <- c()
              for (i2 in covar) { 
                if(var(colData(object)[,i2], na.rm = TRUE) > 0) {covarKeep <- c(covarKeep, i2)}
              }
              
              if (length(covar) > length(covarKeep)) {
                warning(sprintf(
                  "The following covariate(s) have zero covariance: %s",
                  paste(covar[!covar %in% covarKeep], collapse = ",")
                ))
              }
              covar <- covarKeep
              
              ## Remove samples with missing covariates 
              
              kp <- complete.cases(colData(object)[,covar])
              if (!all(kp) && verbose) message(sprintf("%s samples have missing values for one or more covariates, these will be dropped.", sum(!kp)))
              keepSamples <- keepSamples & kp
            } 
            
            
            # Sample filtering --------------------------------------------------
            
            if(minCallrateSM > 0 || maxCallrateSM < Inf) {
              callRate <- Matrix::colMeans(!is.na(SummarizedExperiment::assays(object)$GT))
              kp <- callRate >= minCallrateSM & callRate <= maxCallrateSM
              if(verbose) message(sprintf("%s samples don't pass callRate filters.", sum(!kp)))
              keepSamples <- keepSamples & kp
            } 
            
            if(verbose) message(sprintf("Keeping %s/%s available samples for analysis.", sum(keepSamples),ncol(object)))
            if(!all(keepSamples)) object <- object[,keepSamples]
            
            if(sum(keepSamples) == 0) {
              if(singlevar) return(singlevarResult()) else return(rvbResult())
            }
            
            ### Flip to minor after sample filtering
            object <- flipToMinor(object)
            
            # Variant filtering --------------------------------------------------
            
            ## Check weights: variants with weight<0 or weight == 0 are removed.
            nvar <- nrow(object)
            
            if(sum(is.na(rowData(object)$w)) > 0) {
              warning(sprintf("For %s variants weights are missing, these variants are excluded.",
                              sum(is.na(rowData(object)$w))))
              object <- object[!is.na(rowData(object)$w),]
            }
            
            if(sum(rowData(object)$w < 0) > 0) {
              
              warning(sprintf("%s weights are < 0, these variants are excluded.",
                              sum(rowData(object)$w < 0)))
              
              object <- object[rowData(object)$w >= 0,]
            }
            
            ## Variant filtering, skip if filter parameters are set to the defaults.
            if(
              minCallrateVar > 0 ||
              maxCallrateVar < Inf ||
              minCallrateSM > 0 ||
              maxCallrateSM < Inf ||
              minMAF > 0 ||
              maxMAF < 1 ||
              minMAC > 0 ||
              maxMAC < Inf ||
              minCarriers > 0 ||
              maxCarriers < Inf ||
              minCarrierFreq > 0 ||
              maxCarrierFreq < Inf
            ) {
              
              sumgeno <- summariseGeno(object)
              if (geneticModel %in% c("allelic", "dominant") || 
                  S4Vectors::metadata(object)$geneticModel != "allelic") {
                carriers = (sumgeno[,"geno1"] + sumgeno[,"geno2"])
              } else if (geneticModel == "recessive") {
                carriers = (sumgeno[,"geno2"])
              }
              
              ## Callrate filtering
              keepGeno <- 
                sumgeno[,"callRate"] >= minCallrateVar & 
                sumgeno[,"callRate"] <= maxCallrateVar &
                carriers >= minCarriers & 
                carriers <= maxCarriers &
                (carriers/ncol(object)) >= minCarrierFreq & 
                (carriers/ncol(object)) <= maxCarrierFreq
              
              ## MAF/MAC filters (also if geneticModel != "allelic", but impossible if metadata(object)$geneticModel != "allelic")
              if(S4Vectors::metadata(object)$geneticModel == "allelic") {
                keepGeno <- keepGeno &
                  sumgeno[,"AF"] >= minMAF & 
                  sumgeno[,"AF"] <= maxMAF & 
                  sumgeno[,"geno1"] + 2*sumgeno[,"geno2"] >= minMAC & 
                  sumgeno[,"geno1"] + 2*sumgeno[,"geno2"] <= maxMAC 
              }
              
              if(verbose) message(sprintf("%s/%s variants are retained for analysis", sum(keepGeno), nvar))
              object <- object[keepGeno,]
              
              if(sum(keepGeno) == 0) {
                if(singlevar) return(singlevarResult()) else return(rvbResult())
              }
            } else if (verbose) {
              message(sprintf("%s/%s variants are retained for analysis", nrow(object), nvar))
            }
            
            # Specify models -----------------------------------------------------
            
            if ( (length(covar) == 0 || (length(covar) == 1 && covar == "1")) && is.null(offset) ) {
              null <- as.formula(sprintf("%s ~ 1", pheno))
              covar <- c()
            } else {
              if ( is.null(offset) ) {
                null <- as.formula(paste(pheno, paste(covar, collapse = " + "), sep = " ~ "))
              } else {
                null <- as.formula(sprintf("%s ~ %s + offset(%s)", pheno, paste(covar, collapse = " + "), offset))
              }
            }
            
            if ( is.null(offset) ) {
              model <- as.formula(paste(pheno, paste(c(covar,"aggregate"), collapse = " + "), sep = " ~ "))
              model.nbinom <- as.formula(paste("aggregate", paste(c(covar, pheno), collapse = " + "), sep = " ~ "))
              
            } else {
              model <- as.formula(sprintf("%s ~ %s + offset(%s) + aggregate", pheno, paste(covar, collapse = " + "), offset))
              model.nbinom <- as.formula(sprintf("aggregate ~ %s + offset(%s) + %s", paste(covar, collapse = " + "), offset, pheno))
            }
            
            if(singlevar) {
              
              # Calculate case/ctrl MAF before applying geneticModel
              if (S4Vectors::metadata(object)$geneticModel == "allelic") {
                caseMAF <- getAF(object[,colData(object)[,pheno] == 1])
                names(caseMAF) <- rownames(object)
                ctrlMAF = if(!continuous) getAF(object[,colData(object)[,pheno] == 0]) else rep(NA_real_, nrow(object))
                names(ctrlMAF) <- rownames(object)
              } else {
                caseMAF = rep(NA_real_, nrow(object))
                ctrlMAF = rep(NA_real_, nrow(object))
              }
              
            }
            
            ## Recode genetic model
            if(geneticModel != S4Vectors::metadata(object)$geneticModel) {
              object <- recode(object, geneticModel = geneticModel)
            }
            
            ## drop variants with 0 carriers
            if (minCarriers  == 0) {
              carriers <- getNCarriers(object)
              if (any(carriers == 0)) {
                if(verbose) message(sprintf("%s/%s variants have zero carriers, these are dropped.",
                                            sum(carriers == 0), nrow(object)))
                object <- object[carriers>=1,]
                
                if(nrow(object) == 0) {
                  if(singlevar) return(singlevarResult()) else return(rvbResult())
                }
              }
            }
            
            # rvb -------------------------------------------------------------
            
            if ( !singlevar ) {
              callRate <- Matrix::colMeans(!is.na(SummarizedExperiment::assays(object)$GT))
              results <- data.frame(
                varSetName = S4Vectors::metadata(object)$varSetName,
                cohort = S4Vectors::metadata(object)$cohort,
                name = name,
                unit = S4Vectors::metadata(object)$unit,
                pheno = pheno,
                covar = paste(covar, collapse = ","),
                geneticModel = S4Vectors::metadata(object)$geneticModel,
                MAFweight = if(MAFweights == "none") "1" else MAFweights,
                test = test,
                nvar = nrow(object),
                caseCarriers = if(!continuous) sum(Matrix::colSums(assays(object)$GT >= 1, na.rm = TRUE)[colData(object)[,pheno] == 1] >= 1) else sum(Matrix::colSums(assays(object)$GT >= 1, na.rm = TRUE) >= 1),
                ctrlCarriers = if(!continuous) sum(Matrix::colSums(assays(object)$GT >= 1, na.rm = TRUE)[colData(object)[,pheno] == 0] >= 1) else NA_real_,
                caseN = if(!continuous) sum(colData(object)[, pheno] == 1, na.rm = TRUE) else ncol(object),
                ctrlN = if(!continuous) sum(colData(object)[, pheno] == 0, na.rm = TRUE) else 0,
                meanCaseScore = NA_real_,
                meanCtrlScore = NA_real_,
                caseCallRate = if(!continuous) mean(callRate[colData(object)[,pheno] == 1]) else mean(callRate),
                ctrlCallRate = if(!continuous) mean(callRate[colData(object)[,pheno] == 0]) else NA_real_,
                effect = NA_real_,
                effectSE = NA_real_,
                effectCIlower = NA_real_,
                effectCIupper = NA_real_,
                P = NA_real_,
                stringsAsFactors = FALSE
              )
              
              if( overwriteAggregate ) {
                object <- aggregate(recode(object, 
                                             imputeMethod = imputeMethod, 
                                             MAFweights = MAFweights),checkMissing=FALSE)
              } else {
                if(S4Vectors::metadata(object)$imputeMethod != imputeMethod) {
                  object <- recode(object, 
                                   imputeMethod = imputeMethod)
                }
              }
              
              results$meanCaseScore <- if(!continuous) mean(colData(object)[colData(object)[,pheno] == 1, "aggregate"]) else mean(colData(object)[,"aggregate"])
              results$meanCtrlScore <- if(!continuous) mean(colData(object)[colData(object)[,pheno] == 0, "aggregate"]) else NA_real_
              
              
              ## permuted results
              if(!is.null(resamplingMatrix)) {
                # Check dimensions 
                if(ncol(object) != nrow(resamplingMatrix)) {
                  stop("Number of rows in resamplingMatrix should match with the number of samples in the GT object.")
                }
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
                  methodResampling = methodResampling)
                
              } else if(!is.null(methodResampling)) {
                
                ## Chunks 
                chunks <-rep(memlimitResampling, times = nResampling %/% memlimitResampling)
                if(sum(chunks)-nResampling != 0) chunks <- c(chunks, nResampling-sum(chunks))
                
                permResult <- list()
                
                for(i in 1:length(chunks))  {
                  perms <- buildResamplingFile(nSamples = ncol(object),
                                      methodResampling = methodResampling, 
                                      nResampling = chunks[[i]] )
                  colnames(perms) <- if(i==1) paste0("perm", 1:chunks[i]) else paste0("perm", 
                                                                                      (cumsum(chunks)[(i-1)]+1):(cumsum(chunks)[(i)]))
                  
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
                    methodResampling = methodResampling)
                  
                  permResult[[i]] <- prms
                }
                permResult <- do.call(rbind, permResult)
              }
              
              # format results, and return/output if outputResampling = TRUE or outputResampling = <filepath>
              if((!is.null(resamplingMatrix) || !is.null(methodResampling)) && outputResampling != FALSE) {
                permResult <- permResult[order(permResult$test),]
                permResult <- cbind(
                  S4Vectors::DataFrame(
                    varSetName = S4Vectors::Rle(results$varSetName[1], lengths = nrow(permResult)),
                    cohort = S4Vectors::Rle(results$cohort[1], lengths = nrow(permResult)),
                    name = S4Vectors::Rle(results$name[1], lengths = nrow(permResult)),
                    unit = S4Vectors::Rle(results$unit[1], lengths = nrow(permResult)),
                    covar = S4Vectors::Rle(results$covar[1], lengths = nrow(permResult)),
                    geneticModel = S4Vectors::Rle(results$geneticModel[1], lengths = nrow(permResult)),
                    MAFweight = S4Vectors::Rle(results$MAFweight[1], lengths = nrow(permResult))
                  ),
                  S4Vectors::DataFrame(permResult)
                )
                rownames(permResult) <- NULL
                if(is.character(outputResampling)) {
                  write.table(as.data.frame(permResult), 
                              file = outputResampling, 
                              sep = "\t", 
                              quote = FALSE, 
                              append = FALSE, 
                              row.names = FALSE)
                } else if (outputResampling) {
                  return(permResult)
                } 
              } else if (!is.null(resamplingMatrix) || !is.null(methodResampling)) {
                res <- .rvb_tests_rvb(
                  GT = object,
                  results = results,
                  test = test,
                  pheno = pheno,
                  model = model,
                  model.nbinom = if(!continuous) model.nbinom else NULL,
                  null = null,
                  covar = covar,
                  continuous = continuous,
                  output = output,
                  append = append,
                  maxitFirth = maxitFirth,
                  returnDF = TRUE
                )
                results <- rbind(res,get_perm_pvals(res,permResult))
                
                if(!is.null(output)) {
                  results <- rvbResult(results)
                  writeResult(results, file = output, append = append)
                  return(results)
                } else {
                  if(returnDF) return(results[,names(columns_rvbResults)]) else return(rvbResult(results))
                }
                
              } else {
                .rvb_tests_rvb(
                  GT = object,
                  results = results,
                  test = test,
                  pheno = pheno,
                  model = model,
                  model.nbinom = if(!continuous) model.nbinom else NULL,
                  null = null,
                  covar = covar,
                  continuous = continuous,
                  output = output,
                  append = append,
                  maxitFirth = maxitFirth,
                  returnDF = returnDF
                )
              }
              
            # singlevar -----------------------------------------------------------------------
            } else {
              
              results <- data.frame(
                varSetName = S4Vectors::metadata(object)$varSetName,
                cohort = as.character(S4Vectors::metadata(object)$cohort),
                name = name,
                geneticModel = S4Vectors::metadata(object)$geneticModel,
                VAR_id = rep(rownames(object), each = length(test)),
                pheno = pheno,
                test = rep(test, times = nrow(object)),
                covar = paste(covar, collapse=","),
                caseN = if(!continuous) rep(Matrix::rowSums(!is.na(assays(object)$GT[,colData(object)[,pheno] == 1,drop=FALSE]), na.rm = TRUE), each = length(test)) else rep(Matrix::rowSums(!is.na(assays(object)$GT), na.rm = TRUE), each = length(test)),
                ctrlN = if(!continuous) rep(Matrix::rowSums(!is.na(assays(object)$GT[,colData(object)[,pheno] == 0,drop=FALSE]), na.rm = TRUE), each = length(test)) else 0,
                caseMAF = rep(caseMAF[rownames(object)], each = length(test)),
                caseMAC = if(!continuous) rep(Matrix::rowSums(assays(object)$GT[,colData(object)[,pheno] == 1,drop=FALSE], na.rm = TRUE), each = length(test)) else rep(Matrix::rowSums(assays(object)$GT, na.rm = TRUE), each = length(test)),
                caseCallRate = NA_real_,
                ctrlMAF = rep(ctrlMAF[rownames(object)], each = length(test)),
                ctrlMAC = if(!continuous) rep(Matrix::rowSums(assays(object)$GT[,colData(object)[,pheno] == 0,drop=FALSE], na.rm = TRUE), each = length(test)) else NA_real_,
                ctrlCallRate =  NA_real_,
                stringsAsFactors = FALSE
              )
              if (!continuous) {
                results$caseCallRate <- results$caseN / sum(colData(object)[, pheno] == 1, na.rm = TRUE)
                results$ctrlCallRate <- results$ctrlN / sum(colData(object)[, pheno] == 0, na.rm = TRUE)
              } else {
                results$caseCallRate <- results$caseN/ncol(object)
              }
              
              
              # Impute
              if(!is.null(imputeMethod)) {
                object <- recode(object, imputeMethod = imputeMethod)
              }
              
              ## Run tests for each variant 
              .rvb_tests_singlevar(
                GT = object, 
                results = results, 
                pheno = pheno,
                test = test, 
                model = model, 
                model.nbinom = if(!continuous) model.nbinom else NULL,
                null = null,
                covar = covar, 
                continuous = continuous,
                output = output,
                append = append,
                returnDF = returnDF,
                maxitFirth = maxitFirth,
                verbose = verbose
              )
            }
          })


# rvb tests -------------------------------------------------------------------
.rvb_tests_rvb <- function(GT, 
                           results, 
                           test, 
                           pheno,
                           model, 
                           model.nbinom, 
                           null, 
                           covar, 
                           continuous,
                           output = NULL,
                           append = FALSE,
                           maxitFirth = 1000,
                           returnDF = FALSE
) {
  
  P <- effect <- effectSE <- effectCIupper <- effectCIlower <-
    rep(NA_real_, length(test))
  names(P) <- names(effect) <- names(effectSE) <- names(effectCIupper) <- names(effectCIlower) <-
    test
  
  ## skip non-sensible tests
  if( sum(colData(GT)$aggregate > 0) < 2 ) { 
    warning("Less than two samples have a non-zero burden score, skipping tests.")
    test <- c()
  }
  
  if(continuous) {
    out_type <- "C"
  } else {
    out_type <- "D"
    if (sum(colData(GT)[,pheno] == 1) < 2 | sum(colData(GT)[,pheno] == 0) < 2 ) {
      warning("Fewer than two cases or controls, skipping tests.")
      test <- c()
    }
  }
  
  # Burden test (continuous)
  if ("lm" %in% test)
  {
    tryCatch(
      {
        fit <- lm(model, data = colData(GT))
        effect["lm"] <- fit$coefficients["aggregate"]
        effectSE["lm"] <- summary(fit)$coef["aggregate",2]
        effectCIlower["lm"] <- confint(fit)["aggregate",1]
        effectCIupper["lm"] <- confint(fit)["aggregate",2]
        P["lm"] <- summary(fit)$coef["aggregate",4]
      },
      error=function(e){message(sprintf("Failed test '%s'\n%s", "lm", e))}
    )
  }
  
  # Burden test (binary)
  if ("firth" %in% test)
  {
    tryCatch(
      {
        fit <- logistf::logistf(model, data=colData(GT), 
                                plconf = (which(c(covar,"aggregate") == "aggregate")+1),
                                control = logistf::logistf.control(maxit=maxitFirth),
                                plcontrol = logistf::logistpl.control(maxit=maxitFirth)
                                )
        
        ## check if converged, if not, set to NA
        if (.check_conv_firth(fit, maxit=maxitFirth)) {
          effect["firth"] <- exp(fit$coefficients["aggregate"])
          effectSE["firth"] <- sqrt(effect["firth"] * diag(vcov(fit)))[names(fit$coefficients) == "aggregate"]
          effectCIlower["firth"] <- exp(fit$ci.lower["aggregate"])
          effectCIupper["firth"] <- exp(fit$ci.upper["aggregate"])
          P["firth"] <- fit$prob["aggregate"]
        } else {
          effect["firth"]=effectSE["firth"]=effectCIlower["firth"]=effectCIupper["firth"]=P["firth"]=NA
        }
        
      },
      error=function(e){message(sprintf("Failed test '%s'\n%s", "firth", e))}
    )
  }
  
  if ("glm" %in% test)
  {
    tryCatch(
      {
        fit <- glm(model, data=colData(GT), family = "binomial")
        effect["glm"] <- exp(summary(fit)$coef["aggregate", 1])
        effectSE["glm"] <- sqrt(effect["glm"] * diag(vcov(fit)))["aggregate"]
        effectCIlower["glm"] <- exp(confint.default(fit)["aggregate",1])
        effectCIupper["glm"] <- exp(confint.default(fit)["aggregate",2])
        P["glm"] <- summary(fit)$coef["aggregate", 4]
      },
      error=function(e){message(sprintf("Failed test '%s'\n%s", "glm", e))}
    )
  }
  
  if ("scoreSPA" %in% test)
  {
    tryCatch(
      {
        score.null <- SPAtest::ScoreTest_wSaddleApprox_NULL_Model(
          null,
          data = as.data.frame(colData(GT))
        )
        fit <- SPAtest::ScoreTest_SPA(
          genos = colData(GT)[["aggregate"]],
          obj.null = score.null,
          minmac = 0,
          beta.out = FALSE
        )
        
        P["scoreSPA"] <- fit$p.value
      },
      error=function(e){message(sprintf("Failed test '%s'\n%s", "scoreSPA", e))}
    )
  }
  
  if ("nbinom" %in% test)
  {
    tryCatch(
      {
        fit <- MASS::glm.nb(model.nbinom, data = colData(GT))
        effect["nbinom"] <- summary(fit)$coefficients[pheno,1]
        effectSE["nbinom"] <- summary(fit)$coefficients[pheno,2]
        effectCIlower["nbinom"] <- effect["nbinom"]-(1.96*summary(fit)$coefficients[pheno,2])
        effectCIupper["nbinom"] <- effect["nbinom"]+(1.96*summary(fit)$coefficients[pheno,2])
        P["nbinom"] <- summary(fit)$coefficients[pheno,4]
      },
      error=function(e){message(sprintf("Failed test '%s'\n%s", "nbinom", e))}
    )
  }
  
  # SKAT analyses
  if (sum(c("skat_burden","skat","skato", "skat_burden_robust", "skat_robust", "skato_robust") %in% test) > 0)
  {
    af <- rowData(GT)$AF
    skat.null <- SKAT::SKAT_Null_Model(null, data = colData(GT), out_type = out_type)
    
    if ("skat_burden" %in% test)
    {
      tryCatch(
        {
          if(continuous) {
            fit <- SKAT::SKAT(t(assays(GT)$GT),
                              skat.null, 
                              r.corr = 1,
                              impute.method = "fixed",
                              weights = rowData(GT)$w,
                              missing_cutoff = 1)
          } else {
            fit <- SKAT::SKATBinary(t(assays(GT)$GT),
                                    skat.null, 
                                    method="Burden",
                                    impute.method="fixed",
                                    weights=rowData(GT)$w,
                                    missing_cutoff=1)
          }
          
          P["skat_burden"] <- fit$p.value
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "skat_burden", e))}
      )
    }
    
    if ("skat" %in% test)
    {
      tryCatch(
        {
          if(continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null,
              r.corr = 0,
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          } else {
            fit <- SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "SKAT",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }
          
          P["skat"] <- fit$p.value
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "skat", e))}
      )
    }
    
    if ("skato" %in% test)
    {
      tryCatch(
        {
          if(continuous) {
            fit = SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          } else {
            fit <- SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }
          
          P["skato"] <- fit$p.value
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "skato", e))}
      )
    }
    
    if(sum(c("skat_burden_robust", "skat_robust", "skato_robust") %in% test) > 0) {
      if(S4Vectors::metadata(GT)$geneticModel != "allelic") {
        varkeep <- Matrix::rowSums(assays(GT)$GT) > 0
      } else {
        varkeep <- af > 0
      }
    }
    
    if ("skat_burden_robust" %in% test)
    {
      if (sum(varkeep) >= 1) {
        tryCatch(
          {
            fit = SKAT::SKATBinary_Robust(
              t(assays(GT)$GT)[, varkeep, drop = FALSE],
              skat.null,
              method = "Burden",
              impute.method = "fixed",
              weights = rowData(GT)$w[varkeep],
              missing_cutoff = 1
            )
            P["skat_burden_robust"] = fit$p.value
          },
          error=function(e){message(sprintf("Failed test '%s'\n%s", "skat_burden_robust", e))}
        )
      }
    }
    
    if ("skat_robust" %in% test)
    {
      if (sum(varkeep) >= 1) {
        tryCatch(
          {
            fit = SKAT::SKATBinary_Robust(
              t(assays(GT)$GT)[, varkeep, drop = FALSE],
              skat.null,
              method = "SKAT",
              impute.method = "fixed",
              weights = rowData(GT)$w[varkeep],
              missing_cutoff = 1
            )
            P["skat_robust"] = fit$p.value
          },
          error=function(e){message(sprintf("Failed test '%s'\n%s", "skat_robust", e))}
        )
      }
    }
    
    if ("skato_robust" %in% test)
    {
      if (sum(varkeep) >= 1) {
        tryCatch(
          {
            fit = SKAT::SKATBinary_Robust(
              t(assays(GT)$GT)[, varkeep, drop = FALSE],
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w[varkeep],
              missing_cutoff = 1
            )
            P["skato_robust"] <- fit$p.value
          },
          error=function(e){message(sprintf("Failed test '%s'\n%s", "skato_robust", e))}
        )
      }
    }
  }
  
  # skat with resampling for fwe estimation
  if (sum(c("skat_burden_fwe","skat_fwe","skato_fwe") %in% test) > 0)
  {
    
    skat.null = SKAT::SKAT_Null_Model(
      null,
      data = colData(GT),
      out_type = out_type,
      n.Resampling = skatResampling
    )
    
    if ("skat_burden_fwe" %in% test)
    {
      tryCatch(
        {
          if(continuous) {
            fit=SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null, 
              r.corr = 1,
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1)
          } else {
            fit = SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "Burden",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }
          
          P["skat_burden_fwe"] <- SKAT::Get_Resampling_Pvalue(fit)$p.value
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "skat_burden_fwe", e))}
      )
    }
    
    if ("skat_fwe" %in% test)
    {
      tryCatch(
        {
          if(continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null, 
              r.corr=0,
              impute.method="fixed",
              weights=rowData(GT)$w,
              missing_cutoff=1)
          } else {
            fit = SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "SKAT",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }
          
          P["skat_fwe"] <- SKAT::Get_Resampling_Pvalue(fit)$p.value
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "skat_fwe", e))}
      )
    }
    
    if ("skato_fwe" %in% test)
    {
      tryCatch(
        {
          if(continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
            
          } else {
            fit = SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }
          
          P["skato_fwe"] <- SKAT::Get_Resampling_Pvalue(fit)$p.value
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "skato_fwe", e))}
      )
    }
  }
  
  if(sum(c("acatv","acatvSPA", "acatvfirth") %in% test) > 0) {
    
    
    if("acatv" %in% test) {
      
      tryCatch(
        {
          acat_null=.acat_NULL_Model(colData(GT)[[pheno]],
                                     Z = if(!is.null(covar)) as.matrix(colData(GT)[,covar, drop = FALSE]) else NULL)
          
          # fit
          mac<-Matrix::rowSums(assays(GT)$GT)
          maf<-Matrix::rowSums(assays(GT)$GT)/(2*ncol(GT))
          fit <- .acatv_rvat(
            t(assays(GT)$GT),
            obj = acat_null, 
            method = "original",
            weights=rowData(GT)$w,
            mac.thresh=10,
            maf=maf,
            mac=mac
            )
          
          P["acatv"] <- fit
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "acatv", e))}
      )
    }
    
    if(any(c("acatvSPA", "acatvfirth") %in% test)) {
      mac<-Matrix::rowSums(assays(GT)$GT)
      maf<-Matrix::rowSums(assays(GT)$GT)/(2*ncol(GT))
    }
      
    if(any(c("acatvSPA") %in% test)) {
      tryCatch(
        {
          # null model
          obj <- SPAtest::ScoreTest_wSaddleApprox_NULL_Model(
            null,
            data = as.data.frame(colData(GT))
          )
          
          # fit
          fit <- .acatv_rvat(
            t(assays(GT)$GT),
            obj = obj, 
            method = "scoreSPA",
            weights=rowData(GT)$w,
            mac.thresh=10,
            maf=maf,
            mac=mac)
          
          P["acatvSPA"] <- fit
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "acatvSPA", e))}
      )
    }
    
    if("acatvfirth" %in% test) {
      tryCatch(
        {
          
          # fit
          fit <- .acatv_rvat(
            t(assays(GT)$GT),
            obj = NULL, 
            data = colData(GT),
            covar = covar,
            model = model, 
            method = "firth",
            weights=rowData(GT)$w,
            mac.thresh=10,
            maf=maf,
            mac=mac,
            maxitFirth=maxitFirth
            )
          
          P["acatvfirth"] <- fit
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "acatvfirth", e))}
      )
    }
    
  }
  
  results$effect <- effect
  results$effectSE <- effectSE
  results$effectCIlower <- effectCIlower
  results$effectCIupper <- effectCIupper
  results$P <- P
  
  if(!is.null(output)) {
    results <- rvbResult(results)
    writeResult(results, file = output, append = append)
    return(results)
  } else {
    if(returnDF) return(results[,names(columns_rvbResults)]) else return(rvbResult(results))
  }
}

# Singlevar tests -------------------------------------------------------------

.rvb_tests_singlevar <- function(
  GT, 
  results, 
  pheno,
  test, 
  model, 
  model.nbinom,
  null,
  covar, 
  continuous,
  output = NULL,
  append = FALSE,
  returnDF = FALSE,
  maxitFirth=1000,
  verbose = TRUE
) {
  
  testable <- which(Matrix::rowSums(results[,c("caseMAC","ctrlMAC")][!duplicated(results$VAR_id),], na.rm=TRUE) >= 2)
  if (length(testable) < nrow(GT)) {
    if (verbose) warning(sprintf("%s/%s variants have less than 2 carriers, tests will be skipped for these variants.", 
                         nrow(GT)-length(testable),
                         nrow(GT)
                         ))
  }
  Pl=effectl=effectSEl=effectCIlowerl=effectCIupperl=list()
  
  # linear model (continuous pheno)
  if ("lm" %in% test)
  {
    P=effect=effectSE=effectCIupper=effectCIlower <- rep(NA_real_, nrow(GT))
    for(i in testable) {
      tryCatch(
        {
          fit <- lm(model, data=cbind(aggregate = assays(GT)$GT[i,], data.frame(colData(GT))))
          effect[i] <- fit$coefficients["aggregate"]
          effectSE[i] <- summary(fit)$coef["aggregate",2]
          effectCIlower[i] <- confint(fit)["aggregate",1]
          effectCIupper[i] <- confint(fit)["aggregate",2]
          P[i] <- summary(fit)$coef["aggregate",4]
        },
        error=function(e){message(sprintf("Failed test '%s' for variant '%s'\n%s", 
                                          "lm", rownames(GT)[i],e))}
      )
    }
    Pl[["lm"]]=P;effectl[["lm"]]=effect;effectSEl[["lm"]]=effectSE;effectCIlowerl[["lm"]]=effectCIlower;effectCIupperl[["lm"]]=effectCIupper
  }
  
  # Firth logistic regression
  if ("firth" %in% test)
  {
    P=effect=effectSE=effectCIupper=effectCIlower <- rep(NA_real_, nrow(GT))
    for(i in testable) {
      tryCatch(
        {
          fit <- logistf::logistf(model, 
                                  data = cbind(aggregate = assays(GT)$GT[i,], data.frame(colData(GT))),
                                  plconf = (which(c(covar,"aggregate") == "aggregate")+1),
                                  control = logistf::logistf.control(maxit=maxitFirth),
                                  plcontrol = logistf::logistpl.control(maxit=maxitFirth)
                                  )
          if (.check_conv_firth(fit, maxit=maxitFirth)) {
            effect[i] <- exp(fit$coefficients["aggregate"])
            effectCIlower[i] <- exp(fit$ci.lower["aggregate"])
            effectCIupper[i] <- exp(fit$ci.upper["aggregate"])
            P[i] <- fit$prob["aggregate"]
          } 
          
        }, error=function(e){message(sprintf("Failed test '%s' for variant '%s'\n%s", 
                                             "firth", rownames(GT)[i],e))}
      )
    }
    Pl[["firth"]]=P;effectl[["firth"]]=effect;effectSEl[["firth"]]=effectSE;effectCIlowerl[["firth"]]=effectCIlower;effectCIupperl[["firth"]]=effectCIupper
  }
  
  # glm logistic regression
  if ("glm" %in% test)
  {
    P=effect=effectSE=effectCIupper=effectCIlower <- rep(NA_real_, nrow(GT))
    for(i in testable) {
      tryCatch(
        {
          fit <- glm(model, 
                     data = cbind(aggregate = assays(GT)$GT[i,], data.frame(colData(GT))), 
                     family = binomial(link = "logit"))
          effect[i] <- exp(summary(fit)$coef["aggregate", 1])
          effectCIlower[i] <- exp(confint.default(fit)["aggregate",1])
          effectCIupper[i] <- exp(confint.default(fit)["aggregate",2])
          P[i] <- summary(fit)$coef["aggregate", 4]
        },
        error=function(e){message(sprintf("Failed test '%s' for variant '%s'\n%s", 
                                          "glm", rownames(GT)[i],e))}
      )
    }
    Pl[["glm"]]=P;effectl[["glm"]]=effect;effectSEl[["glm"]]=effectSE;effectCIlowerl[["glm"]]=effectCIlower;effectCIupperl[["glm"]]=effectCIupper
  }
  if ("nbinom" %in% test) {
    P=effect=effectSE=effectCIupper=effectCIlower <- rep(NA_real_, nrow(GT))
    for(i in testable) {
      tryCatch(
        {
          fit <- MASS::glm.nb(model.nbinom, data = cbind(aggregate = assays(GT)$GT[i,], data.frame(colData(GT))))
          effect[i] <- summary(fit)$coefficients[pheno,1]
          effectSE[i] <- summary(fit)$coefficients[pheno,2]
          effectCIlower[i] <- effect[i]-(1.96*summary(fit)$coefficients[pheno,2])
          effectCIupper[i] <- effect[i]+(1.96*summary(fit)$coefficients[pheno,2])
          P[i] <- summary(fit)$coefficients[pheno,4]
        },
        error=function(e){message(sprintf("Failed test '%s' for variant '%s'\n%s", 
                                          "nbinom", rownames(GT)[i],e))}
      )
      
    }
    Pl[["nbinom"]]=P;effectl[["nbinom"]]=effect;effectSEl[["nbinom"]]=effectSE;effectCIlowerl[["nbinom"]]=effectCIlower;effectCIupperl[["nbinom"]]=effectCIupper
  }
  
  if ("scoreSPA" %in% test) {
    P=effect=effectSE=effectCIupper=effectCIlower <- rep(NA_real_, nrow(GT))
    
    if ( length(testable) > 0 ) {
      score.null <- SPAtest::ScoreTest_wSaddleApprox_NULL_Model(
        null,
        data = as.data.frame(colData(GT))
      )
      fit <- SPAtest::ScoreTest_SPA(
        genos = assays(GT[,names(score.null$y)])$GT[testable,],
        obj.null = score.null, 
        minmac = 0,
        beta.out = FALSE
      )
      P[testable] <- fit$p.value
    }
    
    Pl[["scoreSPA"]]=P;effectl[["scoreSPA"]]=effect;effectSEl[["scoreSPA"]]=effectSE;effectCIlowerl[["scoreSPA"]]=effectCIlower;effectCIupperl[["scoreSPA"]]=effectCIupper
  }
  
  if(continuous) {
    res <- cbind(results,
                 effect = effectl[["lm"]],
                 effectSE = effectSEl[["lm"]],
                 effectCIlower = effectCIlowerl[["lm"]],
                 effectCIupper= effectCIupperl[["lm"]],
                 P = Pl[["lm"]]
    )
  } else {
    res <- cbind(results,
                 P = c(do.call(rbind, Pl[test])),
                 effect = c(do.call(rbind, effectl[test])),
                 effectSE = c(do.call(rbind, effectSEl[test])),
                 effectCIlower = c(do.call(rbind, effectCIlowerl[test])),
                 effectCIupper= c(do.call(rbind, effectCIupperl[test])))
  }
  
  if(is.null(output)) {
    if(returnDF) return(res[,names(columns_singlevarResults)]) else return(singlevarResult(res))
  } else {
    writeResult(singlevarResult(res), file = output, append = append)
  }
  
  if(!is.null(output)) {
    res <- singlevarResult(res)
    writeResult(res, file = output, append = append)
    return(res)
  } else {
    rownames(res) <- NULL
    if(returnDF) return(res[,names(columns_singlevarResults)]) else return(singlevarResult(res))
  }
}

# ACAT-v ----------------------------------------------------------------------
# adapted from: https://github.com/yaowuliu/ACAT
.acatv_rvat <- function(G,
                        obj = NULL, 
                        data = NULL,
                        covar = NULL,
                        model = NULL,
                        method = c("original", "scoreSPA", "firth"),
                        weights.beta=c(1,25),
                        weights=NULL,
                        mac.thresh=10,
                        maf=NULL,
                        mac=NULL,
                        maxitFirth=1000
                        )
{
  
  method <- match.arg(method)
  
  ### check if null model is provided (if method = "original")
  if (method == "original" && is.null(obj)) {
    stop("if method == 'original', a NULL model should be provided")
  }
  
  ### check if weights match length of genotypes
  if (!is.null(weights) && length(weights) != ncol(G)){
    stop("The length of weights must equal to the number of variants!")
  }
  
  n<-nrow(G)
  if(is.null(mac)) mac <- Matrix::colSums(G,na.rm=TRUE)
  if(is.null(maf)) maf <- mac/(2*n)
  p <- length(mac)
  
  ### remove SNPs with mac=0
  if (sum(mac == 0) > 0){
    G <- G[,mac > 0, drop = FALSE]
    weights <- weights[mac > 0]
    mac <- mac[mac > 0]
    maf <- maf[mac > 0]
    if (length(mac) == 0) {
      stop("The genotype matrix does not have non-zero elements!")
    }
  }
  
  is.very.rare <- mac <= mac.thresh
  
  ## Very rare
  if (sum(is.very.rare) > 0) {
    
    if (method == "original") {
      pval.very.rare <- .acat_burden(G[,is.very.rare,drop=FALSE],
                                     obj, 
                                     weights.beta = weights.beta, 
                                     weights = weights[is.very.rare])
    } else if (method == "scoreSPA") {
      w <- if(is.null(weights)) dbeta(maf[is.very.rare],weights.beta[1],weights.beta[2]) else weights[is.very.rare]
      agg <- Matrix::rowSums(
        as(G[,is.very.rare,drop=FALSE], "sparseMatrix") %*%
          diag(w,
               ncol=length(w),
               nrow=length(w))
      )
      pval.very.rare <- SPAtest::ScoreTest_SPA(
        genos = agg,
        obj.nul = obj,
        minmac = 0
      )$p.value
    } else if (method == "firth") {
      w <- if(is.null(weights)) dbeta(maf[is.very.rare],weights.beta[1],weights.beta[2]) else weights[is.very.rare]
      agg <- Matrix::rowSums(
        as(G[,is.very.rare,drop=FALSE], "sparseMatrix") %*%
          diag(w,
               ncol=length(w),
               nrow=length(w))
      )
      fit <- logistf::logistf(model, 
                              data = cbind(aggregate = agg, as.data.frame(data[,colnames(data)!="aggregate"])),
                              plconf = (which(c(covar,"aggregate") == "aggregate")+1),
                              control = logistf::logistf.control(maxit = maxitFirth),
                              plcontrol = logistf::logistpl.control(maxit = maxitFirth))
      if (.check_conv_firth(fit, maxit=maxitFirth)) {
        pval.very.rare <- fit$prob["aggregate"]
      } else {
        pval.very.rare <- NULL
      }
    }
  }
  
  ## !Very rare
  if (sum(!is.very.rare) > 0) {
    if (is.null(weights)) {
      weights.transformed <- (dbeta(maf[!is.very.rare],weights.beta[1],weights.beta[2])/dbeta(maf[!is.very.rare],0.5,0.5))^2
    } else {
      weights.transformed <- (weights[!is.very.rare]/dbeta(maf[!is.very.rare],0.5,0.5))^2
    }
    
    if (method == "original") {
      Mpvals <- .acat_Get.marginal.pval(G[,!is.very.rare,drop=FALSE],obj)
    } else if (method == "scoreSPA") {
      Mpvals <- SPAtest::ScoreTest_SPA(
        genos = t(G[,!is.very.rare,drop=FALSE]),
        obj.null = obj,
        minmac = 0)$p.value
    } else if (method == "firth") {
      Mpvals <- rep(NA_real_, sum(!is.very.rare))
      names(Mpvals) <- names(is.very.rare)[!is.very.rare]
      for(i in which(!is.very.rare)) {
        fit <- logistf::logistf(model, 
                                data = cbind(aggregate = G[,i], as.data.frame(data[,colnames(data)!="aggregate"])),
                                plconf = (which(c(covar,"aggregate") == "aggregate")+1),
                                control = logistf::logistf.control(maxit=maxitFirth),
                                plcontrol = logistf::logistpl.control(maxit=maxitFirth))
        if (.check_conv_firth(fit, maxit=maxitFirth)) {
          Mpvals[colnames(G)[i]] <- fit$prob["aggregate"]
        } 
      }
    }
  }
  
  if(sum(is.very.rare) == 0) {
    pval <- .rvat_ACAT(Mpvals, weights.transformed)
  } else if (sum(!is.very.rare) == 0) {
    pval <- pval.very.rare
  } else {
    pvals <- c(Mpvals, pval.very.rare)
    mafs <- c(maf[!is.very.rare],mean(maf[is.very.rare]))
    
    if (is.null(weights)) {
      weights.transformed <- (dbeta(mafs,weights.beta[1],weights.beta[2])/dbeta(mafs,0.5,0.5))^2
    } else{
      weights.transformed <- (c(weights[!is.very.rare],mean(weights[is.very.rare]))/dbeta(mafs,0.5,0.5))^2
    }
    
    is.keep <- rep(TRUE,length(pvals))
    is.keep[which(pvals==1)] <- FALSE  ## remove p-values of 1.
    pval <- .rvat_ACAT(pvals[is.keep],weights.transformed[is.keep])
  }
  pval
}

# adapted from: https://github.com/yaowuliu/ACAT
.acat_NULL_Model<-function(Y,Z=NULL){
  n<-length(Y)
  #### check the type of Y
  if ((sum(Y==0)+sum(Y==1))==n){
    out_type<-"D"
  }else{
    out_type<-"C"
  }
  #### Add intercept
  Z.tilde<-cbind(rep(1,length(Y)),Z)
  if (out_type=="C"){
    #### estimate of sigma square
    Z.med<-Z.tilde%*%solve(chol(t(Z.tilde)%*%Z.tilde))   ## Z.med%*%t(Z.med) is the projection matrix of Z.tilde
    Y.res<-as.vector(Y-(Y%*%Z.med)%*%t(Z.med))
    sigma2<-sum(Y.res^2)/(n-ncol(Z.med))
    #### output
    res<-list()
    res[["out_type"]]<-out_type
    res[["Z.med"]]<-Z.med
    res[["Y.res"]]<-Y.res
    res[["sigma2"]]<-sigma2
  }else if (out_type=="D"){
    #### fit null model
    g<-glm(Y~0+Z.tilde,family = "binomial")
    prob.est<-g[["fitted.values"]]
    #### unstandarized residuals
    Y.res<-(Y-prob.est)
    ### Sigma when rho=0
    sigma2.Y<-prob.est*(1-prob.est)  ### variance of each Y_i
    ### output
    res<-list()
    res[["out_type"]]<-out_type
    res[["Z.tilde"]]<-Z.tilde
    res[["Y.res"]]<-Y.res
    res[["sigma2.Y"]]<-sigma2.Y
  }
  return(res)
}

# adapted from: https://github.com/yaowuliu/ACAT
.acat_Get.marginal.pval<-function(G,obj){
  ### check obj
  if (names(obj)[1]!="out_type"){
    stop("obj is not calculated from MOAT_NULL_MODEL!")
  }else{
    out_type<-obj[["out_type"]]
    if (out_type=="C"){
      if (!all.equal(names(obj)[2:length(obj)],c("Z.med","Y.res","sigma2"))){
        stop("obj is not calculated from MOAT_NULL_MODEL!")
      }else{
        Z.med<-obj[["Z.med"]]
        Y.res<-obj[["Y.res"]]
        n<-length(Y.res)
        SST<-obj[["sigma2"]]*(n-ncol(Z.med))
      }
    }else if (out_type=="D"){
      if (!all.equal(names(obj)[2:length(obj)],c("Z.tilde","Y.res","sigma2.Y"))){
        stop("obj is not calculated from MOAT_NULL_MODEL!")
      }else{
        Z.tilde<-obj[["Z.tilde"]]
        Y.res<-obj[["Y.res"]]
        sigma2.Y<-obj[["sigma2.Y"]]
        n<-length(Y.res)
      }
    }
  }
  
  if (!"matrix" %in% class(G) && !"dgCMatrix" %in% class(G)){
    stop("The class of G must be matrix or dgCMatrix!")
  }
  
  if (out_type=="C"){
    G_tX.med<-as.matrix(Matrix::crossprod(Z.med,G))
    ### Sigma^2 of G
    Sigma2.G<-Matrix::colSums(G^2)-Matrix::colSums(G_tX.med^2)
    SSR<-as.vector((Y.res%*%G)^2/Sigma2.G)
    SSR[Sigma2.G<=0]<-0
    df.2<-n-1-ncol(Z.med)
    t.stat<-suppressWarnings(sqrt(SSR/((SST-SSR)/df.2)))
    marginal.pvals<-2*pt(t.stat,(n-1-ncol(Z.med)),lower.tail = FALSE)
  }else if (out_type=="D"){
    Z.stat0<-as.vector(Y.res%*%G)
    ### Sigma when rho=0
    tG_X.tilde_sigma2<-as.matrix(Matrix::crossprod(G,Z.tilde*sigma2.Y))
    Sigma2.G<-Matrix::colSums(G^2*sigma2.Y)-diag(tG_X.tilde_sigma2%*%solve(t(Z.tilde)%*%(Z.tilde*sigma2.Y))%*%t(tG_X.tilde_sigma2))
    marginal.pvals<-2*pnorm(abs(Z.stat0)/sqrt(Sigma2.G),lower.tail = FALSE)
  }
  
  return(marginal.pvals)
}

# adapted from: https://github.com/yaowuliu/ACAT
.acat_burden <- function(G,obj,kernel="linear.weighted",weights.beta=c(1,25),weights=NULL){
  ### check obj
  if (names(obj)[1]!="out_type"){
    stop("obj is not calculated from NULL_MODEL!")
  }else{
    out_type<-obj[["out_type"]]
    if (out_type=="C"){
      if (!all.equal(names(obj)[2:length(obj)],c("Z.med","Y.res","sigma2"))){
        stop("obj is not calculated from NULL_MODEL!")
      }else{
        Z.med<-obj[["Z.med"]]
        Y.res<-obj[["Y.res"]]/sqrt(obj[["sigma2"]])  ## rescaled residules
        n<-length(Y.res)
      }
    }else if (out_type=="D"){
      if (!all.equal(names(obj)[2:length(obj)],c("Z.tilde","Y.res","sigma2.Y"))){
        stop("obj is not calculated from NULL_MODEL!")
      }else{
        Z.tilde<-obj[["Z.tilde"]]
        Y.res<-obj[["Y.res"]]
        sigma2.Y<-obj[["sigma2.Y"]]
        n<-length(Y.res)
      }
    }
  }
  p<-ncol(G)
  #### weights
  if (kernel=="linear.weighted"){
    if (is.null(weights)){
      MAF<-Matrix::colSums(G)/(2*dim(G)[1])
      W<-dbeta(MAF,weights.beta[1],weights.beta[2])
    }else{
      if (length(weights)==p){
        W<-weights
      }else{
        stop("The length of weights must equal to the number of variants!")
      }
    }
    
  }else if (kernel=="linear"){
    W<-rep(1,p)
  }else{
    stop("The kernel name is not valid!")
  }
  
  ###### if G is sparse or not
  if ("matrix" %in% class(G) || "dgCMatrix" %in% class(G)){
    if (out_type=="C"){
      Z.stat.sum<-as.vector((Y.res%*%G)%*%W)
      Gw<-G%*%W
      sigma.z<-sqrt(sum(Gw^2)-sum((t(Z.med)%*%(Gw))^2))
    }else if (out_type=="D"){
      Z.stat.sum<-as.vector((Y.res%*%G)%*%W)
      Gw<-as.vector(G%*%W)
      sigma.z<-sum(Gw^2*sigma2.Y)-((Gw*sigma2.Y)%*%Z.tilde)%*%solve(t(Z.tilde)%*%(Z.tilde*sigma2.Y))%*%t((Gw*sigma2.Y)%*%Z.tilde)
      sigma.z<-as.vector(sqrt(sigma.z))
    }
  }else{
    stop("The class of G must be matrix or dgCMatrix!")
  }
  
  V<-Z.stat.sum/sigma.z
  Q<-V^2   ## Q test statistic
  pval<-1-pchisq(Q,df=1)
  return(pval)
}

# adapted from: https://github.com/yaowuliu/ACAT
.rvat_ACAT <- function(Pvals,weights=NULL,is.check=TRUE){
  Pvals<-as.matrix(Pvals)
  if (is.check){
    #### check if there is NA
    if (sum(is.na(Pvals))>0){
      stop("Cannot have NAs in the p-values!")
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals<0)+sum(Pvals>1))>0){
      stop("P-values must be between 0 and 1!")
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero<-(colSums(Pvals==0)>=1)
    is.one<-(colSums(Pvals==1)>=1)
    if (sum((is.zero+is.one)==2)>0){
      stop("Cannot have both 0 and 1 p-values in the same column!")
    }
    
    if (sum(is.zero)>0){
      warning("There are p-values that are exactly 0!")
    }
    if (sum(is.one)>0){
      warning("There are p-values that are exactly 1!")
    }
    
  }
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(weights)){
    is.weights.null<-TRUE
  }else{
    is.weights.null<-FALSE
    weights<-as.matrix(weights)
    if (sum(dim(weights)!=dim(Pvals))>0){
      stop("The dimensions of weights and Pvals must be the same!")
    }else if (is.check & (sum(weights<0)>0)){
      stop("All the weights must be nonnegative!")
    }else{
      w.sum<-colSums(weights)
      if (sum(w.sum<=0)>0){
        stop("At least one weight should be positive in each column!")
      }else{
        for (j in 1:ncol(weights)){
          weights[,j]<-weights[,j]/w.sum[j]
        }
      }
    }
    
  }
  
  #### check if there are very small non-zero p values and calcuate the cauchy statistics
  is.small<-(Pvals<1e-15)
  if (is.weights.null){
    Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-1/Pvals[is.small]/pi
    cct.stat<-colMeans(Pvals)
  }else{
    Pvals[!is.small]<-weights[!is.small]*tan((0.5-Pvals[!is.small])*pi)
    Pvals[is.small]<-(weights[is.small]/Pvals[is.small])/pi
    cct.stat<-colSums(Pvals)
  }
  #### return the ACAT p value(s).
  pval<-pcauchy(cct.stat,lower.tail = FALSE)
  return(pval)
}

# misc ------------------------------------------------------------------------
.check_conv_firth <- function(fit, maxit = NULL) {
  if(max(fit$pl.iter[,"Lower"]) >= maxit ||
     max(fit$pl.iter[,"Upper"]) >= maxit ||
     max(fit$pl.iter[,"Null model"]) >= maxit) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# Permutations ----------------------------------------------------------------
.permute_rvb <- function(GT, 
                         test, 
                         pheno,
                         model, 
                         null, 
                         covar, 
                         continuous,
                         perms, 
                         methodResampling = "permutation"
) {
  
  # Currently, tests included for permutation:
  P <- list()
  
  if(sum(colData(GT)$aggregate > 0) < 2) { 
    warning("Less than two samples have a non-zero burden score, skipping tests.")
    test <- c()
  }
  
  if(continuous) {
    out_type <- "C"
  } else {
    out_type <- "D"
    if (sum(colData(GT)[,pheno] == 1) < 2 | sum(colData(GT)[,pheno] == 0) < 2 ) {
      warning("Fewer than two cases or controls, skipping tests.")
      test <- c()
    }
  }
  
  # SKAT analyses
  if (sum(c("skat_burden","skat","skato", "skat_burden_robust", "skat_robust", "skato_robust") %in% test) > 0)
  {
    af <- rowData(GT)$AF
    skat.null <- SKAT::SKAT_Null_Model(null, 
                                       data = colData(GT), 
                                       out_type = out_type)
    skat.null$n.Resampling <- ncol(perms)
    skat.null$type.Resampling <- methodResampling
    skat.null$res.out <- apply(perms, 2, function(x) skat.null$res[x])
    
    if ("skat_burden" %in% test)
    {
      tryCatch(
        {
          if(continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null, 
              r.corr = 1,
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1)
          } else {
            fit <- SKAT::SKATBinary(t(assays(GT)$GT),
                                    skat.null, 
                                    method="Burden",
                                    impute.method="fixed",
                                    weights=rowData(GT)$w,
                                    missing_cutoff=1)
          }
          
          P[["skat_burden"]] <- fit$p.value.resampling
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "skat_burden", e))}
      )
      if(length(P[["skat_burden"]]) == 0) P[["skat_burden"]] <- rep(NA_real_, ncol(perms))
    }
    
    if ("skat" %in% test)
    {
      tryCatch(
        {
          if(continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null, 
              r.corr=0,
              impute.method="fixed",
              weights=rowData(GT)$w,
              missing_cutoff=1)
          } else {
            fit <- SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "SKAT",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }
          
          P[["skat"]] <- fit$p.value.resampling
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "skat", e))}
      )
      if(length(P[["skat"]]) == 0) P[["skat"]] <- rep(NA_real_, ncol(perms))
    }
    
    if ("skato" %in% test)
    {
      tryCatch(
        {
          if(continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          } else {
            fit <- SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }
          
          P[["skato"]] <- fit$p.value.resampling
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "skato", e))}
      )
      if(length(P[["skato"]]) == 0) P[["skato"]] <- rep(NA_real_, ncol(perms))
    }
    
    if(sum(c("skat_burden_robust", "skat_robust", "skato_robust") %in% test) > 0) {
      if(S4Vectors::metadata(GT)$geneticModel != "allelic") {
        varkeep <- Matrix::rowSums(assays(GT)$GT) > 0
      } else {
        varkeep <- af > 0
      }
      
      if(sum(varkeep) > 0) {
        skat_robust_perms <- lapply(
          1:ncol(perms),
          .permute_skatrobust_skato,
          gt =  t(assays(GT)$GT)[, varkeep, drop = FALSE],
          skat_null = skat.null,
          weights = rowData(GT)$w[varkeep],
          perm = perms)
        skat_robust_perms <- matrix(unlist(skat_robust_perms), ncol=3, byrow=TRUE)
        colnames(skat_robust_perms) <- c("skat_robust", "skat_burden_robust", "skato_robust")
        P[["skat_robust"]] <- if("skat_robust" %in% test) skat_robust_perms[,"skat_robust"] else NULL
        P[["skato_robust"]] <- if("skato_robust" %in% test) skat_robust_perms[,"skato_robust"] else NULL
        P[["skat_burden_robust"]] <- if("skat_burden_robust" %in% test)  skat_robust_perms[,"skat_burden_robust"]
      }
      if(length(P[["skat_robust"]]) == 0 & "skat_robust" %in% test) P[["skat_robust"]] <- rep(NA_real_, ncol(perms))
      if(length(P[["skato_robust"]]) == 0 & "skato_robust" %in% test) P[["skato_robust"]] <- rep(NA_real_, ncol(perms))
      if(length(P[["skat_burden_robust"]] & "skat_burden_robust" %in% test) == 0) P[["skat_burden_robust"]] <- rep(NA_real_,  ncol(perms))
    }
  }
  
  if(sum(c("acatv") %in% test) > 0) {
    
    acat_null=.acat_NULL_Model(colData(GT)[[pheno]],
                               Z = if(!is.null(covar)) as.matrix(colData(GT)[,covar, drop = FALSE]) else NULL)
    
    
    if("acatv" %in% test) {
      
      tryCatch(
        {
          # run permutations 
          acatv_perms <- unlist(lapply(
            1:ncol(perms),
            .permute_acatv,
            gt = t(assays(GT)$GT),
            acat_null = acat_null,
            perm = perms,
            weights = rowData(GT)$w,
            maf = getAF(GT),
            mac = Matrix::rowSums(assays(GT)$GT)
          ))
          
          P[["acatv"]] <- acatv_perms
        },
        error=function(e){message(sprintf("Failed test '%s'\n%s", "acatv", e))}
      )
      if(length(P[["acatv"]]) == 0) P[["acatv"]] <- rep(NA_real_, ncol(perms))
    }
  }
  
  data.frame(pheno = rep(colnames(perms), 
             times = length(test)),
             test = rep(assocTest_resampling_tests[assocTest_resampling_tests %in% test], each = ncol(perms)),
             P = unname(unlist(P[assocTest_resampling_tests[assocTest_resampling_tests %in% test]])),
             stringsAsFactors = FALSE)
}

.permute_skatrobust_skato <- function(i, gt,skat_null, weights, perms) {
  perm <- perms[,i]
  skat_null_permuted <- skat_null
  skat_null_permuted$res <- skat_null_permuted$res[perm]
  skat_null_permuted$mu <- skat_null_permuted$mu[perm]
  skat_null_permuted$pi_1 <- skat_null_permuted$pi_1[perm]
  skat_null_permuted$X1 <- skat_null_permuted$X1[perm,,drop=FALSE]
  fit = SKAT::SKATBinary_Robust(
    gt,
    skat_null_permuted,
    method = "SKATO",
    impute.method = "fixed",
    weights = weights,
    missing_cutoff = 1
  )
  if(ncol(gt) > 1) {
    c(fit$p.value_each[fit$param$rho==0], 
      fit$p.value_each[fit$param$rho==1],
      fit$p.value
    )
  } else {
    c(fit$p.value, 
      fit$p.value,
      fit$p.value
    )
  }
}

.permute_acatv <- function(i, gt, acat_null, perms,weights, maf, mac) {
  perm <- perms[,i]
  acat_null_permuted <- acat_null
  acat_null_permuted$Z.tilde <- acat_null_permuted$Z.tilde[perm,,drop=FALSE]
  acat_null_permuted$Y.res <- acat_null_permuted$Y.res[perm]
  acat_null_permuted$sigma2.Y <- acat_null_permuted$sigma2.Y[perm]
  .acatv_rvat(gt,
              obj = acat_null_permuted, 
              method = "original",
              weights=weights,
              mac.thresh=10,
              maf=maf,
              mac=mac)
}




get_perm_pval <- function(pobs, pperm) {
  P <- (sum(pperm < pobs) + 1) / (length(pperm) + 1)
  P
}

get_perm_pvals <- function(x,y) {
  results <- lapply(unique(as.character(y$test)), FUN = function(tst,x,y) {
    x_ <- x[x$test == tst,]
    x_$P <- get_perm_pval(x_$P, y[y$test == tst,]$P)
    x_$test <- paste0(tst, "_resampled")
    x_
  }, x = x, y = y)
  do.call(rbind, results)
}

