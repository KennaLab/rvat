rvat_cli_methods <- list(
  
  ## available methods/functions
  optparse::make_option(c("--buildGdb"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--concatGdb"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--uploadAnno"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--mapVariants"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--uploadCohort"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--listAnno"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--listCohort"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--dropTable"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--subsetGdb"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--buildVarSet"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--spatialClust"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--summariseGeno"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--aggregate"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--mergeAggregateFiles"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--collapseAggregateFiles"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--vcfInfo2Table"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--assocTest"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--buildResamplingFile"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--buildGeneSet"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--buildCorMatrix"), 
                        action = "store_true", 
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--geneSetAssoc"), 
                        action = "store_true", 
                        default = NULL,
                        help = "")
)

# parameter list
rvat_cli_parameters <- list(
  
  ## parameters used across multiple methods/functions
  optparse::make_option(c("--vcf"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--gdb"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--output"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--skipIndexes"), 
              action = "store_true", 
              type = "logical",
              default = FALSE,
              help = ""),
  
  optparse::make_option(c("--overWrite"), 
                        action = "store_true", 
                        type = "logical",
                        default = FALSE,
                        help = ""),
  
  optparse::make_option(c("--memlimit"), 
                        action = "store", 
                        type = "integer",
                        default = 1000L,
                        help = ""),
  
  optparse::make_option(c("--skipRemap"), 
                        action = "store_true", 
                        type = "logical",
                        default = FALSE,
                        help = ""),
  
  optparse::make_option(c("--tables"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--sep"), 
              action = "store", 
              type = "character",
              default = "\t",
              help = ""),
  
  optparse::make_option(c("--name"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--value"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--verbose"), 
              action = "store_true", 
              default = TRUE,
              help = ""),
  
  optparse::make_option(c("--quiet"), 
                        action = "store_true", 
                        default = FALSE,
                        help = ""),
  
  optparse::make_option(c("--strict"), 
                        action = "store_true", 
                        default = TRUE,
                        help = ""),
  
  optparse::make_option(c("--not-strict"), 
                        action = "store_true", 
                        default = FALSE,
                        help = ""),
  
  optparse::make_option(c("--where"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--intersection"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--varSetName"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--unitTable"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--unitName"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--weightName"), 
              action = "store", 
              type = "character",
              default = "1",
              help = ""),
  
  optparse::make_option(c("--cohort"), 
              action = "store", 
              type = "character",
              default = "SM",
              help = ""),
  
  optparse::make_option(c("--varSet"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--VAR_id"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--keep"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--pheno"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--geneticModel"), 
              action = "store", 
              type = "character",
              default = "allelic",
              help = ""),
  
  optparse::make_option(c("--checkPloidy"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--minCallrateVar"), 
              action = "store", 
              type = "numeric",
              default = 0,
              help = ""),
  
  optparse::make_option(c("--maxCallrateVar"), 
              action = "store", 
              type = "numeric",
              default = Inf,
              help = ""),
  
  optparse::make_option(c("--minCallrateSM"), 
              action = "store", 
              type = "numeric",
              default = 0,
              help = ""),
  
  optparse::make_option(c("--maxCallrateSM"), 
              action = "store", 
              type = "numeric",
              default = Inf,
              help = ""),
  
  optparse::make_option(c("--minMAF"), 
              action = "store", 
              type = "numeric",
              default = 0,
              help = ""),
  
  optparse::make_option(c("--maxMAF"), 
              action = "store", 
              type = "numeric",
              default = 1,
              help = ""),
  
  optparse::make_option(c("--minMAC"), 
              action = "store", 
              type = "numeric",
              default = 0,
              help = ""),
  
  optparse::make_option(c("--maxMAC"), 
              action = "store", 
              type = "numeric",
              default = Inf,
              help = ""),
  
  optparse::make_option(c("--minCarriers"), 
              action = "store", 
              type = "numeric",
              default = 0,
              help = ""),
  
  optparse::make_option(c("--maxCarriers"), 
              action = "store", 
              type = "numeric",
              default = Inf,
              help = ""),
  
  optparse::make_option(c("--minCarrierFreq"), 
              action = "store", 
              type = "numeric",
              default = 0,
              help = ""),
  
  optparse::make_option(c("--maxCarrierFreq"), 
              action = "store", 
              type = "numeric",
              default = Inf,
              help = ""),
  
  optparse::make_option(c("--maxitFirth"), 
              action = "store", 
              type = "numeric",
              default = 1000,
              help = ""),
  
  optparse::make_option(c("--nResampling"), 
              action = "store", 
              type = "numeric",
              default = 1000,
              help = ""),
  
  optparse::make_option(c("--methodResampling"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--outputResampling"), 
              action = "store", 
              type = "character",
              default = "FALSE",
              help = ""),
  
  optparse::make_option(c("--test"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--imputeMethod"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--MAFweights"), 
              action = "store", 
              type = "character",
              default = "none",
              help = ""),
  
  optparse::make_option(c("--covar"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--seed"), 
                        action = "store", 
                        type = "integer",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--aggregateFile"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--geneSet"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  ## buildGdb ------------------------------------
  
  optparse::make_option(c("--skipVarRanges"), 
              action = "store_true", 
              type = "logical",
              default = FALSE,
              help = ""),
  
  
  optparse::make_option(c("--genomeBuild"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  ## concatGdb ---------------------------------
  
  optparse::make_option(c("--targets"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  
  ## uploadAnno --------------------------------
  optparse::make_option(c("--ignoreAlleles"), 
              action = "store_true", 
              type = "logical",
              default = FALSE,
              help = ""),
  
  optparse::make_option(c("--keepUnmapped"), 
                        action = "store_true", 
                        type = "logical",
                        default = FALSE,
                        help = ""),
  
  optparse::make_option(c("--mapRef"), 
              action = "store", 
              type = "character",
              default = "var",
              help = ""),
  
  ## mapVariants -------------------------------
  optparse::make_option(c("--ranges"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--gff"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--bed"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--bedCols"), 
              action = "store", 
              type = "character",
              default = character(),
              help = ""),
  
  optparse::make_option(c("--fields"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--uploadName"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  ## spatialClust --------------------------------------------------------------
  optparse::make_option(c("--windowSize"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--overlap"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--posField"), 
              action = "store", 
              type = "character",
              default = "POS",
              help = ""),
  
  optparse::make_option(c("--minTry"), 
              action = "store", 
              type = "integer",
              default = 5L,
              help = ""),
  
  ## summariseGeno -------------------------------------------------------------
  
  optparse::make_option(c("--splitBy"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  ## aggregate -----------------------------------------------------------------
  
  optparse::make_option(c("--signif"), 
              action = "store", 
              type = "numeric",
              default = 6,
              help = ""),
  
  ## mergeAggregateFiles -------------------------------------------------------
  
  optparse::make_option(c("--filelist"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--checkDups"), 
                        action = "store_true", 
                        type = "logical",
                        default = TRUE,
                        help = ""),
  
  optparse::make_option(c("--not-checkDups"), 
                        action = "store_true", 
                        type = "logical",
                        default = FALSE,
                        help = ""),
  
  
  ## assocTest-gdb -----------------------------------------------------------------
  
  optparse::make_option(c("--resamplingFile"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--offset"), 
              action = "store", 
              type = "character",
              default = NULL,
              help = ""),
  
  optparse::make_option(c("--continuous"), 
              action = "store_true", 
              type = "logical",
              default = FALSE,
              help = ""),
  
  optparse::make_option(c("--singlevar"), 
              action = "store_true", 
              type = "logical",
              default = FALSE,
              help = ""),
  
  optparse::make_option(c("--memlimitResampling"), 
                        action = "store", 
                        type = "numeric",
                        default = NULL,
                        help = ""),
  
  ## assocTest-aggregateFile ----------------------------------------------------
  optparse::make_option(c("--substractCovar"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--dropUnits"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  ## buildResamplingFile -------------------------------------------------------
  optparse::make_option(c("--nSamples"), 
                        action = "store", 
                        type = "numeric",
                        default = NULL,
                        help = ""),
  
  ## buildGeneSet -------------------------------------------------------------
  optparse::make_option(c("--gmtpath"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  ## buildCorMatrix -----------------------------------------------------------
  optparse::make_option(c("--minR2"), 
                        action = "store", 
                        type = "numeric",
                        default = 1e-04,
                        help = ""),
  
  optparse::make_option(c("--maxDist"), 
                        action = "store", 
                        type = "numeric",
                        default = 2.5e6,
                        help = ""),
  
  optparse::make_option(c("--makePD"), 
                        action = "store_true", 
                        type = "logical",
                        default = TRUE,
                        help = ""),
  
  optparse::make_option(c("--not-makePD"), 
                        action = "store_true", 
                        type = "logical",
                        default = FALSE,
                        help = ""),
  
  optparse::make_option(c("--absolute"), 
                        action = "store_true", 
                        type = "logical",
                        default = TRUE,
                        help = ""),
  
  optparse::make_option(c("--not-absolute"), 
                        action = "store_true", 
                        type = "logical",
                        default = FALSE,
                        help = ""),

  ## geneSetAssoc -------------------------------------------------------------
  optparse::make_option(c("--rvbResult"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),

  optparse::make_option(c("--scoreMatrix"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--cormatrix"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--condition"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--condition_type"), 
                        action = "store", 
                        type = "character",
                        default = "geneSet",
                        help = ""),
  
  optparse::make_option(c("--threshold"), 
                        action = "store", 
                        type = "numeric",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--Zcutoffs"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--minSetSize"), 
                        action = "store", 
                        type = "numeric",
                        default = 1,
                        help = ""),
  
  optparse::make_option(c("--maxSetSize"), 
                        action = "store", 
                        type = "numeric",
                        default = Inf,
                        help = ""),
  
  optparse::make_option(c("--oneSided"), 
                        action = "store_true", 
                        type = "logical",
                        default = TRUE,
                        help = ""),
  
  optparse::make_option(c("--twoSided"), 
                        action = "store_true", 
                        type = "logical",
                        default = FALSE,
                        help = ""),
  
  optparse::make_option(c("--INT"), 
                        action = "store_true", 
                        type = "logical",
                        default = FALSE,
                        help = ""),
  
  optparse::make_option(c("--scoreCutoffs"), 
                        action = "store", 
                        type = "character",
                        default = NULL,
                        help = ""),
  
  optparse::make_option(c("--ID"), 
                        action = "store", 
                        type = "character",
                        default = "unit",
                        help = ""),
  
  ## vcfInfo2Table -------------------------------------------------------------
  optparse::make_option(c("--splitMultiallelic"), 
                        action = "store_true", 
                        type = "logical",
                        default = TRUE,
                        help = ""),
  
  optparse::make_option(c("--not-splitMultiallelic"), 
                        action = "store_true", 
                        type = "logical",
                        default = FALSE,
                        help = ""),
  
  ## help
  optparse::make_option(c("--help"), 
                        action = "store_true", 
                        type = "logical",
                        default = FALSE,
                        help = "")
)

# list all available options
rvat_cli_options <- c(rvat_cli_methods, rvat_cli_parameters)

# list all valid flags
rvat_cli_methods_flags <- unlist(lapply(rvat_cli_methods, 
                                FUN = function(x) {
                                  if(x@action == "store_true" && !is.na(x@long_flag)) {
                                    x@long_flag
                                  } else if (x@action == "store_true" && !is.na(x@short_flag)) {
                                    x@short_flag
                                  }
                                }))
rvat_cli_valid_calls <- paste(rvat_cli_methods_flags, collapse = "\n")
rvat_cli_methods_flags <- gsub("^-{1,2}", "", rvat_cli_methods_flags)

rvat_cli_flags <- unlist(lapply(rvat_cli_options, 
                       FUN = function(x) {
                         if(x@action == "store_true" && !is.na(x@long_flag)) {
                           x@long_flag
                         } else if (x@action == "store_true" && !is.na(x@short_flag)) {
                           x@short_flag
                         }
                       }))
rvat_cli_flags <- gsub("^-{1,2}", "", rvat_cli_flags)