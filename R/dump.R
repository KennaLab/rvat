# dump -------------------------------------------------------------------------
#' @rdname summariseGeno
#' @aliases summariseGeno,gdb-method
#' @export
setMethod("summariseGeno", 
          signature = signature(object="gdb"),
          definition=function(object, 
                              cohort = "SM", 
                              varSet = NULL, 
                              VAR_id = NULL,
                              pheno = NULL, 
                              memlimit = 1000,
                              geneticModel = "allelic", 
                              checkPloidy = NULL,
                              keep = NULL,
                              output = NULL,
                              splitBy = NULL,
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
                              strict = TRUE,
                              verbose = TRUE
                              )
          {
            .dump(
              object = object, 
              cohort = cohort, 
              what = "varSummary",
              varSet = varSet, 
              VAR_id = VAR_id,
              pheno = pheno, 
              memlimit = memlimit,
              geneticModel = geneticModel, 
              checkPloidy = checkPloidy,
              keep = keep,
              output = output,
              splitBy = splitBy,
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
              minCarrierFreq  = minCarrierFreq,
              maxCarrierFreq = maxCarrierFreq,
              verbose = verbose,
              strict = strict
            )
          }
)

#'  Aggregate genotypes into a single (burden) score for each individual
#' 
#'  Returns an aggregate of genotypes for each individual. 
#'  Note, the [`gdb`] implementation is described here, `aggregate` can also be run directly on a 
#'  [`genoMatrix`] object as described in the [`genoMatrix`] documentation.
#'  Specified genetic model, weights, MAF-weighting are taken into account when aggregating.
#'  Aggregates are written to disk in the [`aggregateFile`] format, which can be used as input
#'  for [`assocTest-aggregateFile`] to perform gene set burden analyses.
#' 
#' @param x a [`gdb`] object
#' @param cohort If a valid cohort name is provided, then the uploaded data for this cohort is used to filter and annotate the genotypes 
#' If not specified, all samples in the gdb will be loaded.
#' @param varSet a [`varSetList`] or [`varSetFile`] object.
#' @param VAR_id A vector of VAR_ids, alternatively the varSet parameter can be specified.
#' The `memlimit` argument controls how many variants to aggregate at a time.
#' @param pheno colData field to test as response variable, although not used within this method,
#' this can be useful to filter samples which have missing data for the response variable.
#' @param memlimit Maximum number of variants to load at once (if `VAR_id` is specified).
#' @param geneticModel Which genetic model to apply? ('allelic', 'recessive' or 'dominant').
#' Defaults to `allelic`.
#' @param imputeMethod Which imputation method to apply? ('meanImpute' or 'missingToRef').
#' Defaults to `meanImpute`.
#' @param MAFweights MAF weighting method. Currently Madsen-Browning ('mb') is implemented, by default no MAF weighting is applied.
#' @param checkPloidy Version of the human genome to use when assigning variant ploidy (diploid, XnonPAR, YnonPAR). 
#' Accepted inputs are GRCh37, hg19, GRCh38, hg38.
#' If not specified, the genome build in the [`gdb`] will be used, if available (included in the `genomeBuild` parameter was set in [`buildGdb`]).
#' Otherwise, if the genome build is not included in the gdb metadata, and no value is provided, then all variants are assigned the default ploidy of "diploid"
#' @param keep vector of sample IDs to keep, defaults to `NULL`, in which case all samples are kept.
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
#' @param verbose Should the function be verbose? (TRUE/FALSE), defaults to `TRUE`.
#' @param strict Should strict checks be performed? Defaults to `TRUE`. Strict tests currently includes
#' checking whether supplied varSetFile/varSetList was generated from the same gdb as specified in `object`.
#' 
#' @export
setMethod("aggregate", 
          signature = signature(x="gdb"),
          definition=function(x, 
                              cohort = "SM", 
                              varSet = NULL, 
                              VAR_id = NULL,
                              pheno = NULL, 
                              memlimit = 1000,
                              geneticModel = "allelic", 
                              imputeMethod = "meanImpute",
                              MAFweights = "none", 
                              checkPloidy = NULL,
                              keep = NULL,
                              output = NULL,
                              signif = 6,
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
                              verbose = TRUE,
                              strict = TRUE
                              )
          {
            .dump(object=x, 
              cohort = cohort,
              what = "aggregate",
              varSet = varSet, 
              VAR_id = VAR_id,
              pheno = pheno, 
              memlimit = memlimit,
              geneticModel = geneticModel,
              imputeMethod = imputeMethod,
              MAFweights = MAFweights,
              checkPloidy = checkPloidy,
              keep = keep,
              output = output,
              signif = signif,
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
              minCarrierFreq  = minCarrierFreq,
              maxCarrierFreq = maxCarrierFreq,
              verbose = verbose,
              strict = strict
              )
          }
)

setMethod(".dump", 
          signature = signature(object="gdb"),
          definition=function(object, 
                              cohort = "SM", 
                              what = c("aggregate", "varSummary"),
                              varSet = NULL, 
                              VAR_id = NULL,
                              pheno = NULL, 
                              memlimit = 1000,
                              geneticModel = "allelic", 
                              imputeMethod = "meanImpute",
                              MAFweights = "none", 
                              checkPloidy = NULL,
                              keep = NULL,
                              output = NULL,
                              signif = 3,
                              splitBy = NULL,
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
                              verbose = TRUE,
                              strict = TRUE
                              )
          {
            what <- match.arg(what)
            if(!what %in% c("aggregate","varSummary")) stop("Available options are 'aggregate' and 'varSummary', please check the documentation.")
            
            # Validity checks -------------------------------------------------- 
            
            if((!is.null(varSet) & !is.null(VAR_id)) || (is.null(VAR_id) && is.null(varSet))) {
              stop("Either of `varSet` or `VAR_id` should be specified")
            }
            
            # check if varSet was generated from the current gdb
            if (!is.null(varSet) && strict) {
              .check_gdb_ids(object, varSet)
            }
            
            if (!is.null(VAR_id)) {
              varSet <- .varsTovarSetList(VAR_id, chunkSize = memlimit)
            }
            
            ## Default imputeMethod = "meanImpute"
            if(!imputeMethod %in% c("meanImpute", "missingToRef"))
              stop("`imputeMethod` should be in ('meanImpute', 'missingToRef')")
            
            ## check geneticModels
            if(length(geneticModel) > 1) {
              stop("Only 1 `geneticModel` can be specified.")
            }
            
            if(!all(geneticModel %in% c("allelic", "recessive", "dominant"))) {
              stop("`geneticModel` should be one of ('allelic', 'recessive', 'dominant')")
            } 
            
            # Check if MAFweights are valid
            if(!all(MAFweights %in% c("none", "mb"))) {
              stop("MAFweights should be 'none' or 'mb'")
            }
            
            # open output
            output <- gzcon(file(output,open='wb'))
            j <- 1
            n <-  length(varSet)
            
            for(i in 1:n)
            {
              # Get varset, and retrieve all VAR_ids
              unit <- unique(listUnits(varSet))[i]
              if (verbose) message(sprintf("Analysing %s", unit))
              varset <- getVarSet(varSet, unit = unit)
              if(length(varset) > 1) {
                stop("Only 1 annotation per unit should be included in the varSet")
              }
              varset <- varset[[1]]
              
              # Load genotypes
              GT <- getGT(object, 
                          cohort = cohort,
                          varSet = varset,
                          checkPloidy = checkPloidy,
                          verbose = verbose,
                          strict = strict
                          )
              
              ## subset samples based on keep list
              if(!is.null(keep)) {
                if (verbose) {
                  message(sprintf("Keeping %s/%s samples that are present in the keep-list.",
                                sum(colnames(GT) %in% keep),
                                ncol(GT)))
                }
                GT <- GT[,colnames(GT) %in% keep]
              }
              
              ## subset samples based on `pheno`
              if(!is.null(pheno)) {
                if(!pheno %in% colnames(colData(GT))) stop(sprintf("%s field not present in cohort!", pheno))
                GT <- GT[,!is.na(colData(GT)[,pheno])]
              } 
              
              ## check splitBy parameter
              if(!is.null(splitBy)) {
                # Check if field is present
                if(!splitBy %in% colnames(colData(GT))) stop(sprintf("field '%s' not present in cohort!", splitBy))
                # Check levels of cohort
                splitByfct <- unique(colData(GT)[[splitBy]])
                splitByfct <- splitByfct[!is.na(splitByfct)]
              }
              
              ## write metadata/samples/units to aggregateFile
              if(j == 1 && what == "aggregate") {
                  metadata <- list(
                    rvatVersion = as.character(packageVersion("rvat")),
                    gdbId = getGdbId(object),
                    genomeBuild = getGenomeBuild(object),
                    creationDate = as.character(round(Sys.time(), units = "secs"))
                  )
                  .write_rvat_header(filetype = "aggregateFile", 
                                     metadata = metadata, 
                                     con = output)
                  samples <- colnames(GT)
                  write(paste(samples, collapse = ","), 
                        file = output, append = FALSE) 
                  write(paste(unique(varSet@units), collapse=","), 
                        file = output, append = TRUE)
                }
                
                ## pull weights and remove variants with missing weights
                w <- listWeights(varset)
                names(w) <- listVars(varset)
                if (verbose) message(sprintf("Analysing unit %s; varSet %s", unit, varset@varSetName))
                w <- w[!is.na(w)]
                
                ## Check duplicated variants
                if( sum(duplicated(names(w))) > 0 ) {
                  w <- w[!duplicated(names(w))]
                }
                
                ### Check if any weights are < 0
                if(sum(w < 0) > 0) {
                  warning(sprintf("%s weights are < 0, these are excluded.",
                                  sum(w < 0)))
                  w <- w[w > 0]
                }
                
                GT <- GT[rownames(GT) %in% names(w),]
                w <- w[rownames(GT)]
                GT <- recode(GT, weights=w)
              
              ## Sample filtering  
              if(minCallrateSM > 0 || maxCallrateSM < Inf) {
                callRate <- Matrix::colMeans(!is.na(SummarizedExperiment::assays(GT)$GT))
                sampleKeep <- callRate >= minCallrateSM & callRate <= maxCallrateSM
                if(verbose) message(sprintf("%s samples don't pass callRate filters.", sum(!sampleKeep)))
              } else {
                sampleKeep <- rep(TRUE, ncol(GT))
              }
              
              ## Site filtering
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
                
                # Calculate variant summary
                sumgeno <- summariseGeno(flipToMinor(GT[,sampleKeep]))
                
                if (geneticModel %in% c("allelic", "dominant")) {
                  carriers = (sumgeno[,"geno1"] + sumgeno[,"geno2"])
                } else if (model == "recessive") {
                  carriers = (sumgeno[,"geno2"])
                }
                
                # variant filtering
                varKeep <- 
                  sumgeno[,"callRate"] >= minCallrateVar & 
                  sumgeno[,"callRate"] <= maxCallrateVar & 
                  sumgeno[,"AF"] >= minMAF & 
                  sumgeno[,"AF"] <= maxMAF & 
                  sumgeno[,"geno1"] + 2*sumgeno[,"geno2"] >= minMAC & 
                  sumgeno[,"geno1"] + 2*sumgeno[,"geno2"] <= maxMAC &
                  carriers >= minCarriers & 
                  carriers <= maxCarriers &
                  (carriers/sum(sampleKeep)) >= minCarrierFreq &
                  (carriers/sum(sampleKeep)) <= maxCarrierFreq
                
              } else {
                varKeep <- rep(TRUE, nrow(GT))
              }
              if(sum(sampleKeep) == 0) {
                message("No samples left after filtering!")
                varKeep <- rep(FALSE, nrow(GT))
                sampleKeep <- rep(TRUE, ncol(GT))
                noSamplesLeft <- TRUE
              } else {
                noSamplesLeft <- FALSE
              }
              if(sum(varKeep) == 0) message("No variants left after filtering!")
              
              GT <- GT[varKeep, sampleKeep]
              
              if(what == "varSummary") {
                GT <- recode(GT,
                             geneticModel = geneticModel)
                
                ## splitBy, if specified
                if(!is.null(splitBy)) {
                  sumgeno <- list()
                  for(fct in splitByfct) {
                    sumgeno_ <- summariseGeno(GT[,colData(GT)[[splitBy]] %in% fct])
                    sumgeno[[as.character(fct)]] <- cbind(VAR_id = sumgeno_$VAR_id, 
                                                          splitBy = fct, 
                                                          sumgeno_[,colnames(sumgeno_) != "VAR_id"])
                  }
                  sumgeno <- do.call(rbind, sumgeno)
                  colnames(sumgeno)[2] <- splitBy
                  rownames(sumgeno) <- NULL
                  
                } else {
                  sumgeno <- summariseGeno(GT)
                }
                
                if(i == 1) {
                  write.table(sumgeno, file = output, row.names = FALSE, quote = FALSE, sep = "\t")
                } else{
                  write.table(sumgeno, file = output, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
                }
              } else if(what == "aggregate") {
                GT <- recode(GT,
                             imputeMethod = imputeMethod,
                             geneticModel = geneticModel)
                weight <- w[rownames(GT)]
                counts <- aggregate(flipToMinor(recode(GT, weights = weight, MAFweights = MAFweights)), returnGT=FALSE,checkMissing=FALSE)
                names(counts) <- colnames(GT)
                counts <- counts[samples]
                
                write(paste(as.character(memCompress(as.character(signif(counts, signif)), type='gzip')), collapse = ","), 
                      file = output, append = TRUE) 
              }
              
              j <- j+1
            }
            close(output)
          }
)