#' assocTest-aggregateFile
#' 
#' Run [`assocTest`] on a [`aggregateFile`] object. 
#' See the main [`assocTest`] page for details.

#' @rdname assocTest-aggregateFile
#' @name assocTest-aggregateFile
#' @aliases assocTest,aggregateFile-method
#' @param object a [`aggregateFile`] object (generated using the [`aggregate`] method).
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
#' checking whether supplied aggregateFile was generated from the same gdb as specified in `gdb`.
#' @examples
#' library(rvatData)
#' gdb <- gdb(rvat_example("rvatData.gdb"))
#' 
#' # assocTest-aggregateFile allows for running association tests on pre-constructed aggregates.
#' # below we first generate the aggregates based on a varSetFile
#' varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
#' varset <- getVarSet(varsetfile, unit = c("NEK1", "SOD1", "ABCA4"), varSetName = "High")
#' aggfile <- tempfile()
#' aggregate(x = gdb,
#'           varSet = varset,
#'           maxMAF = 0.001,
#'           output = aggfile,
#'           verbose = FALSE)
#' 
#' # connect to aggregateFile, see ?aggregateFile for more details
#' aggregatefile <- aggregateFile(aggfile)
#' 
#' # build an example genesetlist, see ?buildGeneSet for details
#' genesetlist <- buildGeneSet(
#'   list("geneset1" = c("SOD1", "NEK1"),
#'        "geneset2" = c("ABCA4", "SOD1", "NEK1")))
#' 
#' # perform association tests using assocTest
#' # note that this is very similar to running association tests on a genoMatrix (?`assocTest-genoMatrix`)
#' # or on a gdb (?`assocTest-gdb`). The main difference being that pre-constructed aggregates are used from 
#' # the aggregateFile, and requires a genesetFile/geneSetList (?geneSetFile) to be provided to the `geneSet` argument.
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
setMethod("assocTest", 
          signature = signature(object="aggregateFile"),
          definition=function(object,
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
                              maxitFirth = 1000,
                              keep = NULL,
                              output = NULL,
                              verbose = TRUE,
                              strict = TRUE
          )
          {
            
            # Check if tests are valid
            if(!all(test %in% assocTest_aggregate_tests)) {
              stop(sprintf("The following tests are not valid: %s", 
                           paste(test[!test %in% assocTest_aggregate_tests],
                                 collapse=",")
              ))
            }
            test <- unique(test)
            
            if (strict) {
              .check_gdb_ids(gdb, object, minVersion = "0.3.0")
            }
            
            ## Check if covar is valid
            ### note: check whether covar is either a list or a character value
            ### whether the specified covariates are available is checked at a later stage (when the cohort is loaded)
            if(!is.list(covar)) {
              if(is.null(covar)) {
                covar <- list(NULL) 
              } else {
                if(!is.character(covar)) stop("'covar' parameter should either be a list or a character vector.")
                covar <- list(covar)
              }
            }
            
            ## initialize output
            if(is.null(output)) {
              resContainer <- 
                vector("list", length(geneSet) * length(pheno) * length(covar))
            } else {
              output <- gzcon(file(output,open='wb'))
              write(paste(c("geneSetName", "cohort", "name", "pheno", "covar", "test",
                            "geneSetSize", "genesObs", "caseN", "ctrlN", "meanCaseScore", "meanCtrlScore",
                            "effect", "effectSE", "effectCIlower", "effectCIupper", "OR", "P"
                            ), collapse="\t"), 
                     file = output, append = FALSE) 
            }
           
            j <- 1
            
            ## Load cohort 
            cohort_name <- cohort
            cohort <- getCohort(object = gdb, cohort = cohort)
            cohort <- cohort[!is.na(cohort$IID),]
            
            ## Subset if vector samples to keep is provided
            if(!is.null(keep)) {
              if (verbose) {
                message(sprintf("Keeping %s/%s samples that are present in the keep-list.",
                                sum(cohort[["IID"]] %in% keep),
                                nrow(cohort)))
              }
              
              cohort <- cohort[cohort[["IID"]] %in% keep,,drop = FALSE]
            }
            
            # Check overlap between samples in aggregateFile and cohort
            if(length(intersect(cohort$IID, listSamples(object))) == 0) {
              stop("No overlap between cohort and aggregateFile!")
            }
            if (!all(cohort$IID %in% listSamples(object))) {
              warning(sprintf("%s/%s samples in the cohort are present in the aggregateFile, subsetting the cohort.",
                              sum(cohort$IID %in% listSamples(object)), nrow(cohort)
                              ))
              cohort <- cohort[cohort$IID %in% listSamples(object),,drop=FALSE]
            } 
            
            if(!all(listSamples(object) %in% cohort$IID)) {
              warning(sprintf("%s/%s samples in the aggregateFile are present in the cohort, subsetting the aggregateFile.",
                              sum(listSamples(object) %in% cohort$IID), length(listSamples(object))
                              ))
            }
            
            if(!all(pheno %in% colnames(cohort))) 
              stop(sprintf("The following phenotypes are not present in cohort '%s': %s",
                           cohort_name, paste(pheno[!pheno %in% colnames(cohort)], collapse = ",")))
            
            ## Check if covariates are present in cohort
            if(!all(unlist(covar) %in% colnames(cohort))) 
              stop(sprintf("The following covariates are not present in the cohort: %s",
                           paste(unique(unlist(covar)[!unlist(covar) %in% colnames(cohort)]), collapse = ",")
              ))
            
            if (!is.null(substractCovar)) {
              if (!substractCovar %in% colnames(cohort)) {
                stop(sprintf("substractCovar variable: %s is not available ", substractCovar))
              }
            }
            
            cohort_all <- cohort
            
            for(phen in pheno) {
              if (verbose) message(sprintf("Analysing phenotype: %s", phen))
              if(sum(is.na(cohort_all[[phen]])) != 0) {
                if (verbose) message(sprintf("'%s' is available for %s/%s samples", phen, sum(!is.na(cohort_all[[phen]])), nrow(cohort_all)))
                cohort <- cohort_all[!is.na(cohort[[phen]]),,drop=FALSE]
              } else {
                cohort <- cohort_all
              }
              
              # Loop through units
              for(geneset in unique(names(geneSet)))
              {
                if (verbose) message(sprintf("Analysing %s", geneset))
                
                # Get varset, and retrieve all VAR_ids
                geneSets <- getGeneSet(geneSet, geneSet = geneset)[[1]]
                
                # Check how many units are present in aggregateFile, throw a warning if not all are present
                
                if (verbose) {
                  message(sprintf("%s/%s units in the geneSet are present in the aggregateFile", 
                                  sum(listUnits(geneSets) %in% listUnits(object)),
                                  length(geneSets)))
                }
                
                if(sum(listUnits(geneSets) %in% listUnits(object)) == 0) {
                  next()
                }
                
                # Extract units 
                units <- listUnits(geneSets)[listUnits(geneSets) %in% listUnits(object)]
                if (!is.null(dropUnits) ) {
                  if (!is(dropUnits, "character")) {
                    stop("The `dropUnits` parameter should be a character vector")
                  }
                  units <- units[!units %in% dropUnits]
                } 
                dat <- getUnit(object, unit = units)
                
                # aggregate
                cohort[["aggregate"]] <- colSums(dat[,cohort[["IID"]],drop=FALSE])
                
                for(cov in covar) {
                  
                  cohort_ <- cohort
                  
                  ## parse covariates
                  if(!is.null(cov) && !(length(cov) == 1 && cov == "1")) {
                    ## remove duplicate covariates + check if covariates are present in colData(object)
                    ## generate dummies for categorical covariates
                    ## remove covariates with zero covariance
                    tmp <- .handle_covar(coldata = cohort_, covar = cov, pheno = phen)
                    cov <- tmp[["covar"]]
                    cohort_ <- tmp[["coldata"]]
                  }
                  
                  
                 if (!is.null(substractCovar)) {
                   cohort_[[substractCovar]] <- cohort_[[substractCovar]] - cohort_[["aggregate"]]
                  }

                # Models
                if (length(cov) == 0 || (length(cov) == 1 && cov == "1")) {
                  null <- as.formula(sprintf("%s ~ 1", phen))
                  cov <- c()
                } else {
                  null <- as.formula(paste(phen, paste(cov, collapse = " + "), sep = " ~ "))
                }
                model <- as.formula(paste(phen, paste(c(cov,"aggregate"), collapse = " + "), sep = " ~ "))
                model.nbinom <- as.formula(paste("aggregate", paste(c(cov, phen), collapse = " + "), sep = " ~ "))
                
                
                results <- data.frame(
                  geneSetName = geneset,
                  cohort = cohort_name,
                  name = name,
                  pheno = phen,
                  covar = paste(cov, collapse = ","),
                  test = test,
                  geneSetSize = length(geneSets),
                  genesObs = nrow(dat),
                  caseN = if(!continuous) sum(cohort_[, phen] == 1, na.rm = TRUE) else nrow(cohort_),
                  ctrlN = if(!continuous) sum(cohort_[, phen] == 0, na.rm = TRUE) else 0,
                  meanCaseScore = if(!continuous) mean(cohort_[cohort_[,phen] == 1, "aggregate"]) else mean(cohort_[,"aggregate"]),
                  meanCtrlScore = if(!continuous) mean(cohort_[cohort_[,phen] == 0, "aggregate"]) else NA_real_,
                  effect = NA_real_,
                  effectSE = NA_real_,
                  effectCIlower = NA_real_,
                  effectCIupper = NA_real_,
                  OR = NA_real_,
                  P = NA_real_,
                  stringsAsFactors = FALSE
                )
                
                if(!continuous) {
                  res <- .rvb_tests_rvb_bin_aggregate(
                    cohort = cohort_,
                    results = results,
                    test = test,
                    pheno = phen,
                    model = model,
                    model.nbinom = model.nbinom,
                    null = null,
                    covar = cov,
                    maxitFirth = maxitFirth,
                    output = output,
                    append = TRUE,
                    returnDF = TRUE
                  )
                } else {
                  res <- .rvb_tests_rvb_cont_aggregate(
                    cohort = cohort_,
                    results = results,
                    test = test,
                    pheno = phen,
                    model = model,
                    null = null,
                    covar = cov,
                    output = output,
                    append = TRUE,
                    returnDF = TRUE
                  )
                }
                
                if(is.null(output)) {
                  resContainer[[j]] <- res
                } 
                j <- j + 1
              }
             }
            }
            
            # Return results or close file
            if(is.null(output)) {
              resContainer <- do.call(rbind, resContainer)
              resContainer
            } else {
              close(output)
            }
          }
)

.rvb_tests_rvb_bin_aggregate <- function(
                                  cohort,
                                  results, 
                                  test, 
                                  pheno,
                                  model, 
                                  model.nbinom, 
                                  null, 
                                  covar, 
                                  maxitFirth = 1000,
                                  output = NULL,
                                  append = FALSE,
                                  returnDF = FALSE
) {
  
  P <- OR <- effect <- effectSE <- effectCIupper <- effectCIlower <-
    rep(NA_real_, length(test))
  names(P) <- names(OR) <- names(effect) <- names(effectSE) <- names(effectCIupper) <- names(effectCIlower) <-
    test
  
  if(sum(cohort[,pheno] == 1) < 2 | sum(cohort[,pheno] == 0) < 2 | sum(cohort$aggregate > 0) < 2) {
    test <- c()
  }
  
  if( sum(cohort$aggregate > 0) < 2 ) { 
    warning("Less than two samples have a non-zero burden score, skipping tests.")
    test <- c()
  }

  if (sum(cohort[,pheno] == 1) < 2 | sum(cohort[,pheno] == 0) < 2 ) {
      warning("Fewer than two cases or controls, skipping tests.")
      test <- c()
  }
  
  # Burden test
  if ("firth" %in% test)
  {
    tryCatch(
      {
        fit <- logistf::logistf(model, 
                                data=cohort, 
                                plconf = (which(c(covar,"aggregate") == "aggregate")+1),
                                control = logistf::logistf.control(maxit=maxitFirth),
                                plcontrol = logistf::logistpl.control(maxit=maxitFirth)
                                )
        
        if (.check_conv_firth(fit, maxit=maxitFirth)) {
          effect["firth"] <- fit$coefficients["aggregate"]
          OR["firth"] <- exp(fit$coefficients["aggregate"])
          effectSE["firth"] <- sqrt(diag(vcov(fit)))["aggregate"]
          effectCIlower["firth"] <- fit$ci.lower["aggregate"]
          effectCIupper["firth"] <- fit$ci.upper["aggregate"]
          P["firth"] <- fit$prob["aggregate"]
        } else {
          effect["firth"]=OR["firth"]=effectSE["firth"]=effectCIlower["firth"]=effectCIupper["firth"]=P["firth"]=NA
        }
      },
      error=function(e){message(sprintf("Failed test '%s'\n%s", "firth", e))}
    )
  }
  
  if ("glm" %in% test)
  {
    tryCatch(
      {
        fit <- glm(model, data=cohort, family = "binomial")
        effect["glm"] <- summary(fit)$coef["aggregate", 1]
        OR["glm"] <- exp(summary(fit)$coef["aggregate", 1])
        effectSE["glm"] <- summary(fit)$coef["aggregate", 2]
        effectCIlower["glm"] <- confint.default(fit)["aggregate",1]
        effectCIupper["glm"] <- confint.default(fit)["aggregate",2]
        P["glm"] <- summary(fit)$coef["aggregate", 4]
      },
      error=function(e){message(sprintf("Failed test '%s'\n%s", "glm", e))}
    )
  }
  
  if ("nbinom" %in% test)
  {
    tryCatch(
      {
        fit <- MASS::glm.nb(model.nbinom, data = cohort)
        effect["nbinom"] <- summary(fit)$coefficients[pheno,1]
        OR["nbinom"] <- NA_real_
        effectSE["nbinom"] <- summary(fit)$coefficients[pheno,2]
        effectCIlower["nbinom"] <- effect["nbinom"]-(1.96*summary(fit)$coefficients[pheno,2])
        effectCIupper["nbinom"] <- effect["nbinom"]+(1.96*summary(fit)$coefficients[pheno,2])
        P["nbinom"] <- summary(fit)$coefficients[pheno,4]
      },
      error=function(e){message(sprintf("Failed test '%s'\n%s", "nbinom", e))}
    )
  }
  
  results$effect <- effect
  results$effectSE <- effectSE
  results$effectCIlower <- effectCIlower
  results$effectCIupper <- effectCIupper
  results$OR <- OR
  results$P <- P
  
  if(!is.null(output)) {
    write.table(results, sep = "\t", file = output, append = TRUE, col.names = FALSE, row.names = FALSE,quote=FALSE)
    return(results)
  } else {
    if(returnDF) return(results) else return(results)
  }
 
}

  
.rvb_tests_rvb_cont_aggregate <- function(GT, 
                                          cohort,
                                          results, 
                                          test, 
                                          pheno,
                                          model, 
                                          null, 
                                          covar, 
                                          output = NULL,
                                          append = FALSE,
                                          returnDF = FALSE
) {
  
  P <- OR <- effect <- effectSE <- effectCIupper <- effectCIlower <-
    rep(NA_real_, length(test))
  names(P) <- names(OR) <- names(effect) <- names(effectSE) <- names(effectCIupper) <- names(effectCIlower) <-
    test
  
  if(sum(cohort[["aggregate"]] > 0) < 2) { 
    warning("Less than two samples have a non-zero burden score, skipping tests.")
    test <- c()
  }
  
  # Burden test
  if ("lm" %in% test)
  {
    tryCatch(
      {
        fit <- lm(model, data = cohort)
        effect["lm"] <- fit$coefficients["aggregate"]
        OR["lm"] <- NA_real_
        effectSE["lm"] <- summary(fit)$coef["aggregate",2]
        effectCIlower["lm"] <- confint(fit)["aggregate",1]
        effectCIupper["lm"] <- confint(fit)["aggregate",2]
        P["lm"] <- summary(fit)$coef["aggregate",4]
      },
      error=function(e){message(sprintf("Failed test '%s'\n%s", "lm", e))}
    )
  }
  
  results$effect <- effect
  results$effectSE <- effectSE
  results$effectCIlower <- effectCIlower
  results$effectCIupper <- effectCIupper
  results$OR <- OR
  results$P <- P
  
  if(!is.null(output)) {
    write.table(results, file = output, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE,quote=FALSE)
    return(results)
  } else {
    if(returnDF) return(results) else return(results)
  }
}