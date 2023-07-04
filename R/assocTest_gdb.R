#' assocTest-gdb
#' 
#' Run [`assocTest`] on a [`gdb`] object. See the main [`assocTest`] page for details.
#' 
#' @rdname assocTest-gdb-method
#' @name assocTest-gdb-method
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
#' @param varSet a [`varSetFile`] or [`varSetList`] object, alternatively a vector of VAR_ids can be specified using the `VAR_id` parameter.
#' @param VAR_id a vector of VAR_ids, alternatively the `varSet` parameter can be specified.
#' If single variant tests are ran, the `memlimit` argument controls how many variants to analyze at a time.
#' @param name Optional name for the analysis, defaults to "none".
#' @param continuous If the response variable continuous? (TRUE/FALSE). Defaults to `FALSE`.
#' @param singlevar Run single variant tests? (TRUE/FALSE).
#' Defaults to `FALSE`, in which case aggregate tests are ran.
#' @param covar Character vector of covariates, or a list of character vectors of covariates in which case each covariate set will be tested separately.
#' @param geneticModel Which genetic model to apply? ('allelic', 'recessive' or 'dominant').
#' Defaults to `allelic`.
#' Multiple geneticModels can be specified, in which case each will be analyzed separately.
#' @param imputeMethod Which imputation method to apply? ('meanImpute' or 'missingToRef').
#' Defaults to `meanImpute`.
#' @param MAFweights MAF weighting method. Currently Madsen-Browning ('mb') is implemented, by default no MAF weighting is applied.
#' Multiple MAFweights can be specified (comma-delimited), in which case each will be analyzed separately.
#' @param maxitFirth Maximum number of iterations to use for estimating firth confidence intervals. Defaults to 1000.
#' @param checkPloidy Version of the human genome to use when assigning variant ploidy (diploid,XnonPAR,YnonPAR). 
#' Accepted inputs are GRCh37, hg19, GRCh38, hg38. 
#' If no value is provided then all variants are assigned the default ploidy of "diploid".
#' @param keep Vector of sample IDs to keep, defaults to `NULL`, in which case all samples are kept.
#' @param output Output file path for results.
#' Defaults to `NULL`, in which case results are not written.
#' @param methodResampling Which method to use for resampling? ('permutation' currently implemented)
#' Defaults to `NULL`, in which case no resampling is performed. 
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
#' @param append Relevant if the `output` parameter is not `NULL`. Should results be appended to `output`?
#' Defaults to `FALSE`.
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
#' @export

setMethod("assocTest", 
          signature = signature(object="gdb"),
          definition=function(object, 
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
                              maxitFirth = 1000,
                              checkPloidy = NULL,
                              keep = NULL,
                              output = NULL,
                              methodResampling = NULL,
                              resamplingFile = NULL,
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
                              minCarrierFreq  = 0,
                              maxCarrierFreq = Inf,
                              memlimit = Inf,
                              verbose = TRUE
                              )
          {
            
            # Validity checks -------------------------------------------------
            
            ## Either varSet (varSetList or varSetFile) or VAR_id (string of VAR_ids) should be specified.
            if (is.null(varSet) && is.null(VAR_id)) {
              stop("Either of `varSet` or `VAR_id` should be specified")
            }
            
            if (!is.null(varSet) && !is.null(VAR_id)) {
              stop("Either of one of `varSet` or `VAR_id` should be specified, not both.")
            }
            
            ## Default imputeMethod = "meanImpute" for burden tests. For singlevar the default is not  to impute.
            if(is.null(imputeMethod)) {
              if(!singlevar) imputeMethod <- "meanImpute"
            } else {
              if(!imputeMethod %in% c("meanImpute", "missingToRef"))
                stop("`imputeMethod` should be in ('meanImpute', 'missingToRef')")
            }
            
            ## geneticModel should be in allelic, recessive or dominant
            if(!all(geneticModel %in% c("allelic", "recessive", "dominant"))) {
              stop("`geneticModel` should be in ('allelic', 'recessive', 'dominant')")
            }
            
            ## Check if tests are valid
            if(!all(test %in% assocTest_tests)) {
              stop(sprintf("The following tests are not valid: %s", 
                           paste(test[!test %in% assocTest_tests],
                                 collapse=",")
              ))
            }
            test <- unique(test)
            
            ## Check if MAFweights are valid
            if(!all(MAFweights %in% c("none", "mb"))) {
              stop("MAFweights should be 'none' or 'mb'")
            }
            
            ## Check if covar is valid
            ### note: check whether covar is either a list or a character value
            ### whether the s[pecified covariates are available is checked at a later stage (when the cohort is loaded)
            if(!is.list(covar)) {
              if(is.null(covar)) {
                covar <- list(NULL) 
                } else {
                if(!is.character(covar)) stop("'covar' parameter should either be a list or a character vector.")
                covar <- list(covar)
                }
            }
            
            ## if VAR_id is specified -> create a varSetList with chunks
            if (!is.null(VAR_id)) {
              varSet <- .varsTovarSetList(VAR_id, chunkSize = memlimit)
            }
            
            ## Resampling 
            
            if(!is.null(resamplingFile)) {
              ## Connect to resampling file 
              nResampling <- resamplingFile@nResampling
              methodResampling <- resamplingFile@methodResampling
              
              # If `memlimitResampling` = NULL (the default) set it to `nResampling`
              if(is.null(memlimitResampling)) memlimitResampling <- nResampling
              
            } else if(!is.null(methodResampling)) {
              
              # If `memlimitResampling` = NULL (the default) set it to `nResampling`
              if(is.null(memlimitResampling)) memlimitResampling <- nResampling
              
              if(memlimitResampling >= nResampling && methodResampling == "permutation") {
                ## Output warning in there are a) multiple phenotypes or b) phenotype filers in place
                if(length(pheno) > 1 || minCallrateSM > 0 || maxCallrateSM < Inf) {
                  warning("...")
                }
                chrt <- getCohort(object, cohort)
                if(!is.null(keep)) {
                  chrt <- chrt[!is.na(chrt[[pheno]]) & chrt[["IID"]] %in% keep,,drop=FALSE]
                } else {
                  chrt <- chrt[!is.na(chrt[[pheno]]),,drop=FALSE]
                }
                
                resamplingMatrix <- buildResamplingFile(nSamples = nrow(chrt), nResampling = nResampling)
                colnames(resamplingMatrix) <- paste0("perm", 1:ncol(resamplingMatrix))
              } else {
                resamplingMatrix <- NULL
                resamplingFile <- NULL
              }
            } 
            
            
            ## if resampling output is to be saved, establish file connections/lists to write/store results
            if(!is.null(methodResampling)) {
              
              if ( singlevar ) {
                stop("Resampling is currently not implemented for singlevar tests")
              }
              # container to hold the resampling results
              if(!is.character(outputResampling)) {
                resamplingContainer <- 
                  vector("list", length(varSet@units) * length(pheno) * length(geneticModel) * length(MAFweights) * length(covar))
                
                # establish file connection if 'outputResampling' is a sting
              } else if (is.character(outputResampling)) {
                outputResampling <- gzcon(file(outputResampling,open='wb'))
                if(singlevar) {
                  # write(paste(c("varSetName", "cohort","name", "VAR_id", "covar","geneticModel","pheno","test", "P"), collapse="\t"), 
                  #       file = outputResampling, append = FALSE) } 
                  } else {
                          write(paste(c("varSetName", "cohort","name", "unit", "covar","geneticModel", "MAFweight", "pheno","test", "P"), collapse="\t"), 
                                file = outputResampling, append = FALSE) 
                        }
              }
            }
            
            # If specified, open file, else initialize a list to store results
            if(is.null(output)) {
              resContainer <- 
                vector("list", length(varSet) * length(pheno) * length(geneticModel) * length(MAFweights) * length(covar))
            } else {
              output <- gzcon(file(output,open='wb'))
              if(singlevar) {
                write(paste(names(columns_singlevarResults), collapse="\t"), 
                                  file = output, append = FALSE) } else {
                                    write(paste(names(columns_rvbResults), collapse="\t"), 
                                          file = output, append = FALSE) 
                                  }
            }
            
           # load cohort
           chrt <- getCohort(object, cohort = cohort, fields = unique(c("IID", "sex", pheno, unlist(covar), offset)))
           chrt <- DataFrame(chrt[!is.na(chrt$IID),])
           metadata(chrt)$name <- cohort
            
           j <- 1; reloadGT <- TRUE
           # Loop through units -----------------------------------------------
           for(unit in unique(listUnits(varSet)))
              {
            if(verbose) message(sprintf("Analysing %s", unit))
            loadGT <- TRUE
             
            # Get varset for unit
            varSets <- getVarSet(varSet, unit = unit)
                
              # Loop through phenotypes --------------------------------------
              for(phen in pheno) {
                
                if (loadGT || reloadGT) {
                  
                  # Load genotypes
                  GT <- getGT(object, 
                              cohort = chrt,
                              varSet = if( length(varSets) > 1) collapseVarSetList(varSets) else varSets,
                              checkPloidy = checkPloidy)
                  metadata(GT)$cohort <- cohort
                
                  loadGT <- FALSE
                  
                  ## keep specified samples
                  if(!is.null(keep)) {
                    if(verbose) message(sprintf("Keeping %s/%s samples that are present in the keep-list.",
                                    sum(colnames(GT) %in% keep),
                                    ncol(GT)))
                    GT <- GT[,colnames(GT) %in% keep]
                  }
                  
                  # If cohort is loaded for first time, perform some checks:
                  if(j == 1) {
                    if (length(pheno) > 1) {
                      ## check if missingness is identical across phenotypes
                      ## if so: loadGT = FALSE for subsequent iterations
                      check_pheno <- all(unlist(lapply(as.data.frame(colData(GT)[,pheno[2:length(pheno)],drop=FALSE]), 
                                                       FUN = function(x,y) identical(is.na(x), is.na(y)), 
                                                       y = is.na(colData(GT)[,pheno[1]]))))
                      if(!check_pheno) reloadGT <- TRUE
                    }
                  }
                  
                  # flipToMinor and extract variant summaries
                  GT <- flipToMinor(GT[,!is.na(colData(GT)[,phen])])
                  
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
                    sumgeno <- summariseGeno(GT)
                    rownames(sumgeno) <- rownames(GT)
                  }
                }
                
                # Loop through genetic models ----------------------------------
                for(model in geneticModel) {
                  
                  # Loop through varsets ---------------------------------------
                  for(i in 1:length(varSets)) {
                    
                    ## extract variants+weight in the current varset
                    vars <- varSets[[i]]
                    w <- listWeights(vars)
                    names(w) <- listVars(vars)
                    message(sprintf("Analysing unit %s; varSet %s", unit, vars@varSetName))
                    
                    ## checks + filtering weights
                    
                    ### remove missing weights
                    w <- w[!is.na(w)]
                    
                    ### check duplicated variants
                    if(sum(duplicated(names(w))) > 0) {
                      w <- w[!duplicated(names(w))]
                    }
                    
                    ### Check if any weights are < 0
                    if(sum(w < 0) > 0) {
                        warning(sprintf("%s weights are < 0, these are excluded.",
                                        sum(w < 0)))
                      w <- w[w > 0]
                    }
                    
                    ## check if same variants are included as previous varSet
                    if(i == 1) {
                      varIndex <- (rownames(GT) %in% names(w))
                      varIndexIdentical <- FALSE
                    } else {
                      varIndex_ <- (rownames(GT) %in% names(w))
                      if(identical(varIndex, varIndex_)) {
                        varIndexIdentical <- TRUE
                      } else {
                        varIndexIdentical <- FALSE
                        varIndex <- varIndex_
                      }
                    }
                    
                    ## subset and match 
                    w <- w[rownames(GT)[varIndex]]
                    
                    ## if varIndexIdentical=FALSE, perform filtering
                    if(!varIndexIdentical) {
                      GT_ <- recode(GT[varIndex,], weights=w)
                      
                      if(minCallrateSM > 0 || maxCallrateSM < Inf) {
                        callRate <- colMeans(!is.na(SummarizedExperiment::assays(GT_)$GT))
                        sampleKeep <- callRate >= minCallrateSM & callRate <= maxCallrateSM
                        if (verbose && sum(!sampleKeep) > 0 ) message(sprintf("%s samples don't pass callRate filters.", sum(!sampleKeep)))
                        
                      } else {
                        sampleKeep <- rep(TRUE, ncol(GT_))
                      }
                      
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
                        
                        # Recalculate variant summary if samples have been excluded
                        if(!all(sampleKeep)) {
                          sumgeno_ <- summariseGeno(flipToMinor(GT_[,sampleKeep]))
                        } else {
                          sumgeno_ <- sumgeno[rownames(GT_),,drop=FALSE]
                        }
                        
                        # calcualte number of carriers per variant
                        if (model %in% c("allelic", "dominant")) {
                          carriers = (sumgeno_[,"geno1"] + sumgeno_[,"geno2"])
                        } else if (model == "recessive") {
                          carriers = (sumgeno_[,"geno2"])
                        }
                        
                        # variant filtering
                        varKeep <- 
                          sumgeno_[,"callRate"] >= minCallrateVar & 
                          sumgeno_[,"callRate"] <= maxCallrateVar & 
                          sumgeno_[,"AF"] >= minMAF & 
                          sumgeno_[,"AF"] <= maxMAF & 
                          sumgeno_[,"geno1"] + 2*sumgeno_[,"geno2"] >= minMAC & 
                          sumgeno_[,"geno1"] + 2*sumgeno_[,"geno2"] <= maxMAC &
                          carriers >= minCarriers & 
                          carriers <= maxCarriers &
                          (carriers/sum(sampleKeep)) >= minCarrierFreq &
                          (carriers/sum(sampleKeep)) <= maxCarrierFreq
                        
                      } else {
                        varKeep <- rep(TRUE, nrow(GT_))
                      }
                      
                      ## print warning if no samples are left after filtering
                      if ( sum(sampleKeep) == 0 ) {
                        warning("No samples left after filtering!")
                        varKeep <- rep(FALSE, nrow(GT_))
                        noSamplesLeft <- TRUE
                      } else {
                        noSamplesLeft <- FALSE
                      }
                      ## print warning if no variants are left after filtering
                      if(sum(varKeep) == 0) warning("No variants left after filtering!")
                      if(sum(varKeep) == 0 || sum(sampleKeep) == 0) {
                        varIndex <- rep(FALSE,nrow(GT))
                        next
                      }
                      GT_ <- GT_[varKeep, sampleKeep]
                      
                      ## calculate call-rates before imputing
                      if(!singlevar) {
                        callRate <- colMeans(!is.na(assays(GT_)$GT))
                        caseCallRate = if(!continuous) mean(callRate[colData(GT_)[,phen] == 1]) else mean(callRate)
                        ctrlCallRate = if(!continuous) mean(callRate[colData(GT_)[,phen] == 0]) else NA_real_
                      }
                      
                      ## impute GT
                      GT_ <- recode(GT_,
                                    imputeMethod = imputeMethod,
                                    geneticModel = model)
                    } 
                    metadata(GT_)$unit <- vars@unit
                    metadata(GT_)$varSetName <- if(is.null(VAR_id)) vars@varSetName else "none"
                    weight <- w[rownames(GT_)]
                    
                    # Loop through MAFweights ----------------------------------
                    for(MAFweight in MAFweights) {
                      
                      ## calculate aggregate if mode = rvb and a test is included that is based on aggregates 
                      if(!singlevar) {
                        GT_ <- aggregate(recode(GT_, weights = weight, MAFweights = MAFweight), returnGT=TRUE,checkMissing=FALSE)
                      } else {
                        GT_ <- recode(GT_, weights = weight, MAFweights = MAFweight)
                      }
                       
                      # Loop through covariates --------------------------------
                      for(cov in covar) {
                        
                        ## non-resampled assocTest

                        if(is.null(methodResampling) || (!is.null(resamplingFile) && outputResampling == FALSE)) {
                          
                          res <- assocTest(
                            object = GT_,
                            pheno = phen,
                            test = test,
                            name = name,
                            continuous = continuous, 
                            singlevar = singlevar, 
                            covar = cov,
                            offset = offset,
                            overwriteAggregate = FALSE,
                            geneticModel = model,
                            imputeMethod = imputeMethod,
                            MAFweights = MAFweight,
                            output = NULL,
                            append = FALSE,
                            returnDF = TRUE,
                            maxitFirth = maxitFirth,
                            minCallrateVar = 0,
                            maxCallrateVar = Inf,
                            minCallrateSM = 0,
                            maxCallrateSM =Inf,
                            minMAF = 0,
                            maxMAF = 1,
                            minMAC = 0,
                            maxMAC = Inf,
                            minCarriers = 0,
                            maxCarriers = Inf
                          )
                          
                          # Replace callrate (calculated before imputation)
                          if(!singlevar && nrow(res) > 0) {
                            res$caseCallRate <- caseCallRate
                            res$ctrlCallRate <- ctrlCallRate
                          }
                          
                          if(noSamplesLeft && nrow(res) > 0) {
                            res$caseN <- 0
                            res$ctrlN <- 0
                          }
                          
                          # write file to connection (if specified), otherwise store result
                          if(is.null(resamplingFile)) {
                            if(is.null(output)) {
                              resContainer[[j]] <- res
                            } else {
                              write.table(res, file = output, append = TRUE, sep = "\t", 
                                          col.names = FALSE, row.names = FALSE, quote = FALSE)
                            }
                          }
                        }
                        
                        ## resampled assocTest
                        if(!is.null(methodResampling)) {
                          
                          # if resamplingFile is specified, loop through the resamplingfile 
                          if(!is.null(resamplingFile)) {
                            chunks <-rep(memlimitResampling, times = nResampling %/% memlimitResampling)
                            if(sum(chunks)-nResampling != 0) chunks <- c(chunks, nResampling-sum(chunks))
                            listResampling <- list()
                            
                            # loop through chunks 
                            for(k in 1:length(chunks))  {
                              skip <- if(k==1) 3 else cumsum(chunks)[(k-1)]+3
                              perms <- t(as.matrix(data.table::fread(resamplingFile@path, skip = skip, nrows = memlimitResampling, header = FALSE)))
                              colnames(perms) <- if(k==1) paste0("perm", 1:chunks[k]) else paste0("perm", 
                                                                                                  (cumsum(chunks)[(k-1)]+1):(cumsum(chunks)[(k)]))
                              resResampling <- assocTest(
                                object = GT_,
                                pheno = phen,
                                test = test,
                                name = name,
                                continuous = continuous, 
                                singlevar = singlevar, 
                                covar = cov,
                                overwriteAggregate = FALSE,
                                geneticModel = model,
                                imputeMethod = imputeMethod,
                                MAFweights = MAFweight,
                                output = NULL,
                                append = FALSE,
                                returnDF = TRUE,
                                maxitFirth = maxitFirth,
                                resamplingMatrix = perms,
                                nResampling = nResampling,
                                outputResampling = TRUE,
                                methodResampling = methodResampling,
                                memlimitResampling = memlimitResampling,
                                minCallrateVar = 0,
                                maxCallrateVar = Inf,
                                minCallrateSM = 0,
                                maxCallrateSM =Inf,
                                minMAF = 0,
                                maxMAF = 1,
                                minMAC = 0,
                                maxMAC = Inf,
                                minCarriers = 0,
                                maxCarriers = Inf
                              )
                              
                              ## if output is a connection, if output = TRUE store in listResampling, otherwise, 
                              if(is(outputResampling, "connection")) {
                                write.table(resResampling, file = outputResampling, append = TRUE, sep = "\t", 
                                            col.names = FALSE, row.names = FALSE, quote = FALSE)
                              } else {
                                listResampling[[k]] <- resResampling
                              }
                            }
                            
                            ## collect results if not written to output
                            if(!is(outputResampling, "connection"))  {
                              resResampling <- do.call(rbind, listResampling)
                              resResampling <- resResampling[order(as.character(resResampling$test)),]
                              
                              # if outputResampling is true, store to be returned at the end of the function
                              if(outputResampling == TRUE) {
                                resamplingContainer[[j]] <- resResampling
                                
                                # otherwise calculate resampled P-values
                              } else {
                                ## calculate resampled P-values 
                                results <- rbind(res,get_perm_pvals(res,resResampling))
                                if(is.null(output)) {
                                  resContainer[[j]] <- results
                                } else {
                                  #writeResult(res, file = output, append = TRUE)
                                  write.table(results, file = output, append = TRUE, sep = "\t", 
                                              col.names = FALSE, row.names = FALSE, quote = FALSE)
                                }
                              }
                            }
                            
                          } else {
                            resResampling <- assocTest(
                              object = GT_,
                              pheno = phen,
                              test = test,
                              name = name,
                              continuous = continuous, 
                              singlevar = singlevar, 
                              covar = cov,
                              overwriteAggregate = FALSE,
                              geneticModel = model,
                              imputeMethod = imputeMethod,
                              MAFweights = MAFweight,
                              output = NULL,
                              append = FALSE,
                              returnDF = TRUE,
                              maxitFirth = maxitFirth,
                              resamplingMatrix = resamplingMatrix,
                              nResampling = nResampling,
                              outputResampling = if(outputResampling == FALSE) FALSE else TRUE,
                              methodResampling = methodResampling,
                              memlimitResampling = memlimitResampling,
                              minCallrateVar = 0,
                              maxCallrateVar = Inf,
                              minCallrateSM = 0,
                              maxCallrateSM =Inf,
                              minMAF = 0,
                              maxMAF = 1,
                              minMAC = 0,
                              maxMAC = Inf,
                              minCarriers = 0,
                              maxCarriers = Inf
                            )
                            if(is(outputResampling, "connection")) {
                              write.table(resResampling, file = outputResampling, append = TRUE, sep = "\t", 
                                          col.names = FALSE, row.names = FALSE, quote = FALSE)
                            } else if (outputResampling == TRUE)  {
                              resamplingContainer[[j]] <- resResampling
                            } else {
                              if (is.null(output)) {
                                resContainer[[j]] <- resResampling
                              } else {
                                write.table(resResampling, file = output, append = TRUE, sep = "\t", 
                                            col.names = FALSE, row.names = FALSE, quote = FALSE)
                              }
                            }
                          }
                        }
                        
                        j <- j + 1
                      }
                    }
                  }
                }
              }
           }
            
            if(!is.null(methodResampling)) {
              if(!is.null(outputResampling) && outputResampling == TRUE) {
                resamplingContainer <- do.call(rbind, resamplingContainer)
                if(!is.null(output)) close(output)
                return(resamplingContainer)
              } 
              if(is(outputResampling, "connection")) {
                close(outputResampling)
                }
              } 
            
            # Return results or close file
            if(is.null(output) && (is.null(outputResampling) || !outputResampling)) {
              resContainer <- do.call(rbind, resContainer)
              rownames(resContainer) <- NULL
              return(if(singlevar) singlevarResult(resContainer) else rvbResult(resContainer))
            } else {
              if(!is.null(output) && is(output, "connection")) close(output)
            }
          }
)