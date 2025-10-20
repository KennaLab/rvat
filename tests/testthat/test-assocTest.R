data(GT)
data(GTsmall)

## generate params 
generate_param <- function(
  pheno_bin = "pheno", 
  pheno_quant = "age", 
  covar = paste0("PC",1:4),
  covarCat = "superPop",
  tests_bin_rvb = c("firth", "scoreSPA", "glm", "skat","skat_burden", "skato", "skat_robust", "skat_burden_robust", "skato_robust", "acatv", "acatvSPA"),
  tests_cont_rvb = c("lm", "skat", "skat_burden","skato","acatv"),
  tests_bin_sv =  c("firth", "glm", "scoreSPA"),
  tests_cont_sv = c("lm")
) {
  params <- list()
  ## rvb
  for (mode in c("binary", "continuous")) {
    if(mode == "binary") {
      pheno <- pheno_bin
      tests <- tests_bin_rvb
    } else {
      pheno <- pheno_quant
      tests <- tests_cont_rvb
    }
    param <- list(
      # unfiltered, all tests
      list(
        pheno = pheno,
        test = tests,
        covar = covar
      ),
      # couple of filters, all tests+acatvfirth
      list(
        pheno = pheno,
        test = if(mode=="continuous") tests else c(tests, "acatvfirth"),
        covar = covar,
        maxMAF = 0.001,
        minCarriers = 5,
        minCallrateVar = 0.8,
        minCallrateSM = 0.9
      ),
      
      # couple of filters, all tests, dominant model
      list(
        pheno = pheno,
        test = tests,
        covar = covar,
        geneticModel = "dominant",
        maxMAF = 0.001,
        minCarriers = 5,
        minCallrateVar = 0.8,
        minCallrateSM = 0.9
      ),
      # couple of filters, all tests, missingToRef imputation
      list(
        pheno = pheno,
        test = tests,
        covar = covar,
        imputeMethod = "missingToRef",
        maxMAF = 0.001,
        minCarriers = 5,
        minCallrateVar = 0.8,
        minCallrateSM = 0.9
      ),
      # unfiltered, all tests, recessive model
      list(
        pheno = pheno,
        test = tests,
        covar = covar,
        geneticModel = "recessive"
      ),
      # couple of (different) filters, all tests, MAF-weighted 
      list(
        pheno = pheno,
        test = tests,
        covar = covar,
        MAFweights = "mb",
        minMAF = 0.0001,
        minCarriers = 10
      ),
      # couple of (different) filters, all tests, MAF-weighted 
      list(
        pheno = pheno,
        test = tests,
        covar = covar,
        MAFweights = "mb",
        minMAC=2,
        maxMAC=50,
        minCarrierFreq = 0.00001
      ),
      # unfiltered, dominant, MAF-weighted 
      # list(
      #   pheno = pheno,
      #   test = tests,
      #   covar = covar,
      #   geneticModel = "dominant",
      #   MAFweights = "mb"
      # ),
      # couple of filters, no covar
      list(
        pheno = pheno,
        test = tests,
        minMAC = 2,
        maxMAC = 50,
        minCarrierFreq = 0.00001
      ),
      # couple of filters, include factor covar
      list(
        pheno = pheno,
        test = tests,
        maxCarrierFreq = 0.001,
        covar = c(covar,covarCat)
      )
    )
    params[["rvb"]][[mode]] <- param
  }
  
  
  ## sv
  for (mode in c("binary", "continuous")) {
    if(mode == "binary") {
      pheno <- pheno_bin
      tests <- tests_bin_sv
    } else {
      pheno <- pheno_quant
      tests <- tests_cont_sv
    }
    param <- list(
      # unfiltered, all tests
      list(
        pheno = pheno,
        test = tests,
        covar = covar
      ),
      # couple of filters, all tests
      list(
        pheno = pheno,
        test = tests,
        covar = covar,
        maxMAF = 0.001,
        minCarriers = 5,
        minCallrateVar = 0.8,
        minCallrateSM = 0.9
      ),
      # couple of filters, all tests, dominant model
      list(
        pheno = pheno,
        test = tests,
        covar = covar,
        geneticModel = "dominant",
        maxMAF = 0.001,
        minCarriers = 5,
        minCallrateVar = 0.8,
        minCallrateSM = 0.9
      ),
      
      # unfiltered, all tests, recessive model
      list(
        pheno = pheno,
        test = tests,
        covar = covar,
        geneticModel = "recessive"
      ),
      
      # couple of filters, no covar
      list(
        pheno = pheno,
        test = tests,
        minMAC=2,
        maxMAC=50,
        minCarrierFreq = 0.00001
      ),
      # couple of filters, include factor covar
      list(
        pheno = pheno,
        test = tests,
        maxCarrierFreq = 0.001,
        covar = c(covar,covarCat)
      )
    )
    params[["sv"]][[mode]] <- param
  }
  params
}

run_tests_gdb <- function(gdb, params, VAR_id, var_sv = NULL,cohort="pheno") {
  results <- list()
  
  ## binary burden tests
  results_bin <- list()
  for (i in 1:length(params[["rvb"]][["binary"]])) {
   # print(sprintf("rvb binary: %s/%s", i, length(params[["rvb"]][["binary"]])))
    param_i <- params[["rvb"]][["binary"]][[i]]
    assoc <- assocTest(
      object = gdb,
      VAR_id=VAR_id,
      cohort=cohort,
      pheno = param_i$pheno,
      test = param_i$test,
      covar = if(!is.null(param_i$covar)) param_i$covar else NULL,
      imputeMethod = if(!is.null(param_i$imputeMethod)) param_i$imputeMethod else NULL,
      geneticModel = if(!is.null(param_i$geneticModel)) param_i$geneticModel else "allelic",
      MAFweights = if(!is.null(param_i$MAFweights)) param_i$MAFweights else "none",
      minCallrateVar = if(!is.null(param_i$minCallrateVar)) param_i$minCallrateVar else 0,
      minCallrateSM = if(!is.null(param_i$minCallrateSM)) param_i$minCallrateSM else 0,
      minMAF = if(!is.null(param_i$minMAF)) param_i$minMAF else 0,
      maxMAF = if(!is.null(param_i$maxMAF)) param_i$maxMAF else 1,
      minMAC = if(!is.null(param_i$minMAC)) param_i$minMAC else 0,
      maxMAC = if(!is.null(param_i$maxMAC)) param_i$maxMAC else Inf,
      minCarriers = if(!is.null(param_i$minCarriers)) param_i$minCarriers else 0,
      maxCarriers = if(!is.null(param_i$maxCarriers)) param_i$maxCarriers else Inf,
      minCarrierFreq = if(!is.null(param_i$minCarrierFreq)) param_i$minCarrierFreq else 0,
      maxCarrierFreq = if(!is.null(param_i$maxCarrierFreq)) param_i$maxCarrierFreq else 1,
      verbose = FALSE
    )
    assoc$i <- as.character(i)
    results_bin[[as.character(i)]] <- assoc
  }
  results[["rvb"]][["binary"]] <- as.data.frame(do.call(rbind,results_bin))
  
  ## quantitative burden tests
  results_cont <- list()
  for (i in 1:length(params[["rvb"]][["continuous"]])) {
    #print(sprintf("rvb continuous: %s/%s", i, length(params[["rvb"]][["continuous"]])))
    param_i <- params[["rvb"]][["continuous"]][[i]]
    assoc <- assocTest(
      object = gdb,
      VAR_id=VAR_id,
      cohort=cohort,
      pheno = param_i$pheno,
      test = param_i$test,
      continuous = TRUE,
      covar = if(!is.null(param_i$covar)) param_i$covar else NULL,
      imputeMethod = if(!is.null(param_i$imputeMethod)) param_i$imputeMethod else NULL,
      geneticModel = if(!is.null(param_i$geneticModel)) param_i$geneticModel else "allelic",
      MAFweights = if(!is.null(param_i$MAFweights)) param_i$MAFweights else "none",
      minCallrateVar = if(!is.null(param_i$minCallrateVar)) param_i$minCallrateVar else 0,
      minCallrateSM = if(!is.null(param_i$minCallrateSM)) param_i$minCallrateSM else 0,
      minMAF = if(!is.null(param_i$minMAF)) param_i$minMAF else 0,
      maxMAF = if(!is.null(param_i$maxMAF)) param_i$maxMAF else 1,
      minMAC = if(!is.null(param_i$minMAC)) param_i$minMAC else 0,
      maxMAC = if(!is.null(param_i$maxMAC)) param_i$maxMAC else Inf,
      minCarriers = if(!is.null(param_i$minCarriers)) param_i$minCarriers else 0,
      maxCarriers = if(!is.null(param_i$maxCarriers)) param_i$maxCarriers else Inf,
      minCarrierFreq = if(!is.null(param_i$minCarrierFreq)) param_i$minCarrierFreq else 0,
      maxCarrierFreq = if(!is.null(param_i$maxCarrierFreq)) param_i$maxCarrierFreq else 1,
      verbose = FALSE
    )
    assoc$i <- as.character(i)
    results_cont[[as.character(i)]] <- assoc
  }
  results[["rvb"]][["continuous"]] <- as.data.frame(do.call(rbind,results_cont))
  
  ## binary sv tests
  results_bin <- list()
  for (i in 1:length(params[["sv"]][["binary"]])) {
    #print(sprintf("sv binary: %s/%s", i, length(params[["sv"]][["binary"]])))
    param_i <- params[["sv"]][["binary"]][[i]]
    assoc <- assocTest(
      object = gdb,
      VAR_id=var_sv,
      cohort=cohort,
      pheno = param_i$pheno,
      test = param_i$test,
      covar = if(!is.null(param_i$covar)) param_i$covar else NULL,
      geneticModel = if(!is.null(param_i$geneticModel)) param_i$geneticModel else "allelic",
      minCallrateVar = if(!is.null(param_i$minCallrateVar)) param_i$minCallrateVar else 0,
      minCallrateSM = if(!is.null(param_i$minCallrateSM)) param_i$minCallrateSM else 0,
      minMAF = if(!is.null(param_i$minMAF)) param_i$minMAF else 0,
      maxMAF = if(!is.null(param_i$maxMAF)) param_i$maxMAF else 1,
      minMAC = if(!is.null(param_i$minMAC)) param_i$minMAC else 0,
      maxMAC = if(!is.null(param_i$maxMAC)) param_i$maxMAC else Inf,
      minCarriers = if(!is.null(param_i$minCarriers)) param_i$minCarriers else 0,
      maxCarriers = if(!is.null(param_i$maxCarriers)) param_i$maxCarriers else Inf,
      minCarrierFreq = if(!is.null(param_i$minCarrierFreq)) param_i$minCarrierFreq else 0,
      maxCarrierFreq = if(!is.null(param_i$maxCarrierFreq)) param_i$maxCarrierFreq else 1,
      verbose = FALSE,
      singlevar=TRUE
    )
    
    if(nrow(assoc) > 0) {
      assoc$i <- as.character(i)
      results_bin[[as.character(i)]] <- assoc
    } else {
      results_bin[[as.character(i)]] <- NULL
    }
    results_bin[[as.character(i)]] <- assoc
  }
  results[["sv"]][["binary"]] <- as.data.frame(do.call(rbind,results_bin))
  
  ## continuous sv tests
  results_cont <- list()
  for (i in 1:length(params[["sv"]][["continuous"]])) {
    #print(sprintf("sv continuous: %s/%s", i, length(params[["sv"]][["continuous"]])))
    param_i <- params[["sv"]][["continuous"]][[i]]
    assoc <- assocTest(
      object = gdb,
      VAR_id=var_sv,
      cohort=cohort,
      pheno = param_i$pheno,
      test = param_i$test,
      covar = if(!is.null(param_i$covar)) param_i$covar else NULL,
      geneticModel = if(!is.null(param_i$geneticModel)) param_i$geneticModel else "allelic",
      minCallrateVar = if(!is.null(param_i$minCallrateVar)) param_i$minCallrateVar else 0,
      minCallrateSM = if(!is.null(param_i$minCallrateSM)) param_i$minCallrateSM else 0,
      minMAF = if(!is.null(param_i$minMAF)) param_i$minMAF else 0,
      maxMAF = if(!is.null(param_i$maxMAF)) param_i$maxMAF else 1,
      minMAC = if(!is.null(param_i$minMAC)) param_i$minMAC else 0,
      maxMAC = if(!is.null(param_i$maxMAC)) param_i$maxMAC else Inf,
      minCarriers = if(!is.null(param_i$minCarriers)) param_i$minCarriers else 0,
      maxCarriers = if(!is.null(param_i$maxCarriers)) param_i$maxCarriers else Inf,
      minCarrierFreq = if(!is.null(param_i$minCarrierFreq)) param_i$minCarrierFreq else 0,
      maxCarrierFreq = if(!is.null(param_i$maxCarrierFreq)) param_i$maxCarrierFreq else 1,
      verbose = FALSE,
      singlevar=TRUE,
      continuous=TRUE
    )
    
    if(nrow(assoc) > 0) {
      assoc$i <- as.character(i)
      results_cont[[as.character(i)]] <- assoc
    } else {
      results_cont[[as.character(i)]] <- NULL
    }
    
  }
  results[["sv"]][["continuous"]] <- as.data.frame(do.call(rbind,results_cont))
  results
}

run_tests_GT <- function(GT, params,var_sv = NULL) {
  results <- list()
  
  ## 
  if(is.null(var_sv)) {
    var_sv <- rownames(GT)
  }
  
  ## binary burden tests
  results_bin <- list()
  for (i in 1:length(params[["rvb"]][["binary"]])) {
    #print(sprintf("rvb binary: %s/%s", i, length(params[["rvb"]][["binary"]])))
    param_i <- params[["rvb"]][["binary"]][[i]]
    assoc <- assocTest(
      object = GT,
      pheno = param_i$pheno,
      test = param_i$test,
      covar = if(!is.null(param_i$covar)) param_i$covar else NULL,
      imputeMethod = if(!is.null(param_i$imputeMethod)) param_i$imputeMethod else NULL,
      geneticModel = if(!is.null(param_i$geneticModel)) param_i$geneticModel else "allelic",
      MAFweights = if(!is.null(param_i$MAFweights)) param_i$MAFweights else "none",
      minCallrateVar = if(!is.null(param_i$minCallrateVar)) param_i$minCallrateVar else 0,
      minCallrateSM = if(!is.null(param_i$minCallrateSM)) param_i$minCallrateSM else 0,
      minMAF = if(!is.null(param_i$minMAF)) param_i$minMAF else 0,
      maxMAF = if(!is.null(param_i$maxMAF)) param_i$maxMAF else 1,
      minMAC = if(!is.null(param_i$minMAC)) param_i$minMAC else 0,
      maxMAC = if(!is.null(param_i$maxMAC)) param_i$maxMAC else Inf,
      minCarriers = if(!is.null(param_i$minCarriers)) param_i$minCarriers else 0,
      maxCarriers = if(!is.null(param_i$maxCarriers)) param_i$maxCarriers else Inf,
      minCarrierFreq = if(!is.null(param_i$minCarrierFreq)) param_i$minCarrierFreq else 0,
      maxCarrierFreq = if(!is.null(param_i$maxCarrierFreq)) param_i$maxCarrierFreq else 1,
      verbose = FALSE
    )
    assoc$i <- as.character(i)
    metadata(assoc)$creationDate <- NA_character_
    results_bin[[as.character(i)]] <- assoc
  }
  results[["rvb"]][["binary"]] <- as.data.frame(do.call(rbind,results_bin))
  
  ## quantitative burden tests
  results_cont <- list()
  for (i in 1:length(params[["rvb"]][["continuous"]])) {
    #print(sprintf("rvb continuous: %s/%s", i, length(params[["rvb"]][["continuous"]])))
    param_i <- params[["rvb"]][["continuous"]][[i]]
    assoc <- assocTest(
      object = GT,
      pheno = param_i$pheno,
      test = param_i$test,
      continuous = TRUE,
      covar = if(!is.null(param_i$covar)) param_i$covar else NULL,
      imputeMethod = if(!is.null(param_i$imputeMethod)) param_i$imputeMethod else NULL,
      geneticModel = if(!is.null(param_i$geneticModel)) param_i$geneticModel else "allelic",
      MAFweights = if(!is.null(param_i$MAFweights)) param_i$MAFweights else "none",
      minCallrateVar = if(!is.null(param_i$minCallrateVar)) param_i$minCallrateVar else 0,
      minCallrateSM = if(!is.null(param_i$minCallrateSM)) param_i$minCallrateSM else 0,
      minMAF = if(!is.null(param_i$minMAF)) param_i$minMAF else 0,
      maxMAF = if(!is.null(param_i$maxMAF)) param_i$maxMAF else 1,
      minMAC = if(!is.null(param_i$minMAC)) param_i$minMAC else 0,
      maxMAC = if(!is.null(param_i$maxMAC)) param_i$maxMAC else Inf,
      minCarriers = if(!is.null(param_i$minCarriers)) param_i$minCarriers else 0,
      maxCarriers = if(!is.null(param_i$maxCarriers)) param_i$maxCarriers else Inf,
      minCarrierFreq = if(!is.null(param_i$minCarrierFreq)) param_i$minCarrierFreq else 0,
      maxCarrierFreq = if(!is.null(param_i$maxCarrierFreq)) param_i$maxCarrierFreq else 1,
      verbose = FALSE
    )
    assoc$i <- as.character(i)
    metadata(assoc)$creationDate <- NA_character_
    results_cont[[as.character(i)]] <- assoc
  }
  results[["rvb"]][["continuous"]] <- as.data.frame(do.call(rbind,results_cont))
  
  ## binary sv tests
  results_bin <- list()
  for (i in 1:length(params[["sv"]][["binary"]])) {
    #print(sprintf("sv binary: %s/%s", i, length(params[["sv"]][["binary"]])))
    param_i <- params[["sv"]][["binary"]][[i]]
    assoc <- assocTest(
      object = GT[var_sv,],
      pheno = param_i$pheno,
      test = param_i$test,
      covar = if(!is.null(param_i$covar)) param_i$covar else NULL,
      geneticModel = if(!is.null(param_i$geneticModel)) param_i$geneticModel else "allelic",
      minCallrateVar = if(!is.null(param_i$minCallrateVar)) param_i$minCallrateVar else 0,
      minCallrateSM = if(!is.null(param_i$minCallrateSM)) param_i$minCallrateSM else 0,
      minMAF = if(!is.null(param_i$minMAF)) param_i$minMAF else 0,
      maxMAF = if(!is.null(param_i$maxMAF)) param_i$maxMAF else 1,
      minMAC = if(!is.null(param_i$minMAC)) param_i$minMAC else 0,
      maxMAC = if(!is.null(param_i$maxMAC)) param_i$maxMAC else Inf,
      minCarriers = if(!is.null(param_i$minCarriers)) param_i$minCarriers else 0,
      maxCarriers = if(!is.null(param_i$maxCarriers)) param_i$maxCarriers else Inf,
      minCarrierFreq = if(!is.null(param_i$minCarrierFreq)) param_i$minCarrierFreq else 0,
      maxCarrierFreq = if(!is.null(param_i$maxCarrierFreq)) param_i$maxCarrierFreq else 1,
      verbose = FALSE,
      singlevar=TRUE
    )
    
    if(nrow(assoc) > 0) {
      assoc$i <- as.character(i)
      metadata(assoc)$creationDate <- NA_character_
      results_bin[[as.character(i)]] <- assoc
    } else {
      results_bin[[as.character(i)]] <- NULL
    }
    results_bin[[as.character(i)]] <- assoc
  }
  results[["sv"]][["binary"]] <- as.data.frame(do.call(rbind,results_bin))
  
  ## continuous sv tests
  results_cont <- list()
  for (i in 1:length(params[["sv"]][["continuous"]])) {
    #print(sprintf("sv continuous: %s/%s", i, length(params[["sv"]][["continuous"]])))
    param_i <- params[["sv"]][["continuous"]][[i]]
    assoc <- assocTest(
      object = GT[var_sv,],
      pheno = param_i$pheno,
      test = param_i$test,
      covar = if(!is.null(param_i$covar)) param_i$covar else NULL,
      geneticModel = if(!is.null(param_i$geneticModel)) param_i$geneticModel else "allelic",
      minCallrateVar = if(!is.null(param_i$minCallrateVar)) param_i$minCallrateVar else 0,
      minCallrateSM = if(!is.null(param_i$minCallrateSM)) param_i$minCallrateSM else 0,
      minMAF = if(!is.null(param_i$minMAF)) param_i$minMAF else 0,
      maxMAF = if(!is.null(param_i$maxMAF)) param_i$maxMAF else 1,
      minMAC = if(!is.null(param_i$minMAC)) param_i$minMAC else 0,
      maxMAC = if(!is.null(param_i$maxMAC)) param_i$maxMAC else Inf,
      minCarriers = if(!is.null(param_i$minCarriers)) param_i$minCarriers else 0,
      maxCarriers = if(!is.null(param_i$maxCarriers)) param_i$maxCarriers else Inf,
      minCarrierFreq = if(!is.null(param_i$minCarrierFreq)) param_i$minCarrierFreq else 0,
      maxCarrierFreq = if(!is.null(param_i$maxCarrierFreq)) param_i$maxCarrierFreq else 1,
      verbose = FALSE,
      singlevar=TRUE,
      continuous=TRUE
    )
    
    if(nrow(assoc) > 0) {
      assoc$i <- as.character(i)
      metadata(assoc)$creationDate <- NA_character_
      results_cont[[as.character(i)]] <- assoc
    } else {
      results_cont[[as.character(i)]] <- NULL
    }
    
  }
  results[["sv"]][["continuous"]] <- as.data.frame(do.call(rbind,results_cont))
  results
}

test_that("assocTest-GT works" ,{
  params <- generate_param()
  warnings <- capture_warnings(run_test <- run_tests_GT(GT,params,var_sv=rownames(GT)[1:25]))
  expect_true(all(warnings %in% c("There are p-values that are exactly 0!", "no non-missing arguments to min; returning Inf")))
  expect_snapshot_value(run_test, style = "serialize")
})

test_that("assocTest-gdb and assocTest-GT result in identical output" ,{
  params <- generate_param()
  warnings <- capture_warnings(run_test_gdb <- run_tests_gdb(gdb(rvat_example("rvatData.gdb")), params, VAR_id=rownames(GT), var_sv=rownames(GT)[1:25]))
  expect_true(all(warnings %in% c("There are p-values that are exactly 0!", "no non-missing arguments to min; returning Inf")))
  warnings <- capture_warnings(run_test_GT <- run_tests_GT(GT, params, var_sv=rownames(GT)[1:25]))
  expect_true(all(warnings %in% c("There are p-values that are exactly 0!", "no non-missing arguments to min; returning Inf")))
  run_test_gdb$rvb$binary$varSetName="unnamed"
  run_test_gdb$rvb$binary$unit="unnamed"
  run_test_gdb$rvb$continuous$varSetName="unnamed"
  run_test_gdb$rvb$continuous$unit="unnamed"
  run_test_gdb$sv$binary$varSetName = "unnamed"
  run_test_gdb$sv$continuous$varSetName = "unnamed"
  run_test_gdb$sv$binary$caseMAF = NA_real_
  run_test_gdb$sv$continuous$caseMAF = NA_real_
  run_test_gdb$sv$binary$ctrlMAF = NA_real_
  run_test_gdb$sv$continuous$ctrlMAF = NA_real_
  run_test_GT$sv$binary$caseMAF = NA_real_
  run_test_GT$sv$continuous$caseMAF = NA_real_
  run_test_GT$sv$binary$ctrlMAF = NA_real_
  run_test_GT$sv$continuous$ctrlMAF = NA_real_
  rownames( run_test_GT$sv$continuous) <- NULL
  rownames( run_test_GT$sv$binary) <- NULL
  rownames( run_test_gdb$sv$continuous) <- NULL
  rownames( run_test_gdb$sv$binary) <- NULL
  expect_equal(run_test_gdb,run_test_GT,tolerance=1.490116e-08)
})


gdb <- gdb(rvat_example("rvatData.gdb"))
moderate <- withr::local_tempfile()
LOF <- withr::local_tempfile()
CADD <- withr::local_tempfile()
buildVarSet(object = gdb,
            output = moderate,
            varSetName = "Moderate",
            unitTable = "varInfo",
            unitName = "gene_name",
            where = "(ModerateImpact = 1 or HighImpact = 1) and gene_name in ('CYP19A1', 'FUS', 'OPTN')",
            verbose = FALSE
            )

##
buildVarSet(object = gdb,
            output = LOF,
            varSetName = "HighImpact",
            unitTable = "varInfo",
            unitName = "gene_name",
            where = "HighImpact = 1 and gene_name in ('CYP19A1', 'FUS', 'OPTN')",
            verbose = FALSE
            )

# Build a varset containing CADD weights
buildVarSet(object = gdb,
            output = CADD,
            varSetName = "CADD",
            unitTable = "varInfo",
            unitName = "gene_name",
            weightName = "CADDphred",
            where = "gene_name in ('CYP19A1', 'FUS', 'OPTN')",
            verbose = FALSE
            )

## combined varsets
varsets <- list()
i <- 1
for(file in c(moderate,LOF,CADD)) {
  varset <- readr::read_delim(file, delim = "|", col_names = FALSE, col_types = "cccc", comment = "#", progress = FALSE)
  varsets[[i]] <- varset
  i <- i+1
}
varsets <- dplyr::bind_rows(varsets)
varsets <- varsets %>%
  dplyr::arrange(X1)
metadata <- metadata(varSetFile(moderate))


merged <- withr::local_tempfile()
con <- gzfile(merged, "w")
rvat:::.write_rvat_header(filetype = "varSetFile", metadata = metadata, con = con)
write.table(varsets, file = con, col.names = FALSE, sep="|", append=TRUE, quote = FALSE, row.names = FALSE)
close(con)

test_that("assocTest-gdb works" ,{
  ## looping through either a varSetList or a varSetFile should yield identical result
  varset <- varSetFile(moderate)
  test1 <- assocTest(
    gdb,
    cohort = "pheno",
    varSet = varset,
    pheno = "pheno",
    covar = paste0("PC",1:4),
    test = "scoreSPA",
    verbose = FALSE
  )
  metadata(test1)$creationDate <- NA_character_

  test2 <- assocTest(
    gdb,
    cohort = "pheno",
    varSet = getVarSet(varset, unit = listUnits(varset)),
    pheno = "pheno",
    covar = paste0("PC",1:4),
    test = "scoreSPA",
    verbose = FALSE
  )
  metadata(test2)$creationDate <- NA_character_
  expect_identical(test1,test2)

  ## looping through varSet should give same results as separate varSets (to ensure re-using pheno cohort etc. works as it should)
  varset <- varSetFile(merged)
  warnings <- capture_warnings({test1 <- assocTest(
    gdb,
    cohort = "pheno",
    varSet = varset,
    pheno = "pheno",
    geneticModel = c("allelic", "dominant", "recessive"),
    covar = paste0("PC",1:4),
    maxMAF = 0.0001,
    test = "scoreSPA",
    verbose = FALSE
  )})
  expect_true(all(warnings %in% c("NAs introduced by coercion", "There are p-values that are exactly 0!",
   "no non-missing arguments to min; returning Inf",
    "Less than two samples have a non-zero burden score, skipping tests.") | grepl("weights are missing", warnings)))
  
  test2 <- list()
  for(i in 1:length(varset)) {
    warnings <- capture_warnings({test <- assocTest(
      gdb,
      cohort = "pheno",
      varSet = getVarSet(varset,unit=listUnits(varset)[i], varSetName = listVarSets(varset)[i]),
      pheno = "pheno",
      geneticModel = c("allelic", "dominant", "recessive"),
      covar = paste0("PC",1:4),
      maxMAF = 0.0001,
      test = "scoreSPA",
      verbose = FALSE
    )})
    expect_true(all(warnings %in% c("NAs introduced by coercion", "There are p-values that are exactly 0!",
   "no non-missing arguments to min; returning Inf",
    "Less than two samples have a non-zero burden score, skipping tests.") | grepl("weights are missing", warnings)))
    test2[[i]] <- test
  }
  test2 <- do.call(rbind,test2)
  test1$ID <- paste(test1$varSetName,test1$unit,test1$geneticModel)
  test2$ID <- paste(test2$varSetName,test2$unit,test2$geneticModel)
  test2 <- test2[match(test1$ID,test2$ID),]
  metadata(test1)$creationDate <- NA_character_
  metadata(test2)$creationDate <- NA_character_
  expect_identical(test1,test2)

  # also compare with running directly on genoMatrix
  test3 <- list() 
  j <- 1
  for (i in 1:length(varset)) {
    
    varSet <- getVarSet(varset,
      unit = listUnits(varset)[i], varSetName = listVarSets(varset)[i])
    gt <- suppressWarnings(getGT(
      gdb,
      varSet = varSet,
      cohort = "pheno",
      verbose = FALSE
    ))

    for (geneticmodel in c("allelic", "dominant", "recessive")) {

     test <- suppressWarnings(assocTest(
        gt,
        pheno = "pheno",
        geneticModel = geneticmodel,
        covar = paste0("PC",1:4),
        maxMAF = 0.0001,
        test = "scoreSPA",
        verbose = FALSE
  ))
      test3[[j]] <- test
      j <- j + 1
    }
  }
  test3 <- do.call(rbind,test3)
  test3$ID <- paste(test3$varSetName,test3$unit,test3$geneticModel)
  test3 <- test3[match(test2$ID,test3$ID),]
  metadata(test3)$creationDate <- NA_character_
  expect_identical(test2,test3)
})


test_that("effectAllele is assigned correctly in sv tests" ,{

  GT1 <- GT

  # manually flip some alleles
  flip <- rep(FALSE, nrow(GT1))
  names(flip) <- rownames(GT1)
  flip[c(1, 7, 21)] <- rep(TRUE, 3)
  assays(GT1)$GT <- abs(assays(GT1)$GT - 2 * matrix(rep(flip, each = metadata(GT1)$m), nrow = metadata(GT1)$nvar, byrow = TRUE))
  rowData(GT1)$AF <- getAF(GT1)
  sv <- assocTest(
    GT1,
    covar = paste0("PC", 1:4),
    pheno = "pheno",
    singlevar = TRUE,
    test = c("firth", "scoreSPA", "glm"),
    verbose = FALSE
  )
  var <- as.data.frame(rowData(GT1))
  var$VAR_id <- as.character(rownames(var))
  var <- var %>%
    dplyr::left_join(as.data.frame(sv), by = "VAR_id")
  var_check <- var %>% dplyr::filter(effectAllele != ALT)
  expect_identical(sort(unique(var_check$VAR_id)), as.character(sort(names(which(flip)))))

})


# test resampling
test_that("resampled assocTests work" ,{
  set.seed(10)
  warnings <- capture_warnings({test <- assocTest(
    GT,
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    test = c("skat", "skat_robust", "skat_burden", "skat_burden_robust", "skato_robust", "acatv"),
    nResampling = 10,
    methodResampling = "permutation",
    verbose = FALSE
  )})
  expect_true(all(warnings %in% c("There are p-values that are exactly 0!")))
  metadata(test)$creationDate <- NA_character_
  metadata(test)$gdbPath <- NA_character_
  metadata(test)$rvatVersion <- NA_character_
  expect_snapshot_value(test, style = "serialize")

  # gdb <- gdb(rvat_example("rvatData.gdb"))
  # var <- getAnno(gdb, "varInfo", where = "gene_name = 'CYP19A1'")
  # GT1 <- getGT(
  #   gdb,
  #   VAR_id = var$VAR_id,
  #   anno = "var",
  #   verbose = FALSE
  # )
  # GT1 <- flipToMinor(GT1)

  # set.seed(10)
  # resamplingfile_path <- tempfile(fileext = ".gz")
  # buildResamplingFile(nSamples = ncol(GT1), nResampling = 10, methodResampling = "permutation", output = resamplingfile_path)
  # resamplingfile <- resamplingFile(resamplingfile_path)
  # test_gdb <- assocTest(
  #   gdb,
  #   VAR_id = rownames(GT1),
  #   cohort = "pheno",
  #   pheno = "pheno",
  #   covar = paste0("PC", 1:4),
  #   test = c("skat", "skat_robust", "skat_burden", "skat_burden_robust", "skato_robust", "acatv"),
  #   resamplingFile = resamplingfile,
  #   verbose = FALSE,
  #   outputResampling = TRUE
  # )

  # set.seed(10)
  # test_GT <- assocTest(
  #   gdb,
  #   VAR_id = rownames(GT1),
  #   cohort = "pheno",
  #   pheno = "pheno",
  #   covar = paste0("PC", 1:4),
  #   test = c("skat", "skat_robust", "skat_burden", "skat_burden_robust", "skato_robust", "acatv"),
  #   methodResampling = "permutation",
  #   nResampling = 10,
  #   verbose = FALSE,
  #   outputResampling = TRUE
  # )
  # test_GT$ID = paste(test_GT$unit,test_GT$test, test_GT$geneticModel, test_GT$MAFweight, test_GT$pheno, test_GT$varSetName, sep="_")
  # test_gdb$ID = paste(test_gdb$unit,test_gdb$test, test_gdb$geneticModel, test_gdb$MAFweight, test_gdb$pheno, test_gdb$varSetName, sep="_")
  # test_GT <- test_GT[match(test_gdb$ID, test_GT$ID),]
  # metadata(test_GT)$creationDate <- metadata(test_gdb)$creationDate <- NULL
  # expect_equal(test_gdb, test_GT, tolerance = 1e-10)

  # reSamplingFile tests
  out <- withr::local_tempfile()
  buildResamplingFile(nSamples=1000,
                      nResampling=1000,
                      methodResampling="permutation",
                      output=out)
  resamplingfile <- resamplingFile(out)
  expect_true(stringr::str_detect(capture_output({show(resamplingfile)}), "resamplingFile object"))

  # output resampling
  # compare returning nad writing
  outputresampling <- tempfile()
  set.seed(10)
  test <- assocTest(
    GTsmall,
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    test = c("skat"),
    methodResampling = "permutation",
    nResampling = 10,
    outputResampling = outputresampling,
    verbose = FALSE
  )
  outputresampling <- readr::read_tsv(outputresampling,
    show_col_types = FALSE,
    col_types = list(MAFweight = "character")
  )

  set.seed(10)
  test <- assocTest(
    GTsmall,
    pheno = "pheno",
    covar = paste0("PC", 1:4),
    test = c("skat"),
    methodResampling = "permutation",
    nResampling = 10,
    outputResampling = TRUE,
    verbose = FALSE
  )
  expect_equal(as.data.frame(outputresampling), as.data.frame(test))
})



# input checks

test_that("assocTest input checks work" ,{
  GT_test <- GTsmall
  
  # expect a warning when a covariate with zero covariance is included
  GT_test$covar_novar <- 1
  expect_warning({test <- assocTest(
    GT_test,
    test = "glm",
    covar = c(paste0("PC", 1:4), "covar_novar"),
    pheno = "pheno",
    verbose = FALSE
  )}, regexp = "The following covariate\\(s\\) have zero covariance: covar_novar")
  ## also: the covariate should not be present in the output
  expect_identical(as.character(test$covar), "PC1,PC2,PC3,PC4")

  # expect error when binary phenotype is specified for non-binary phenotype
  expect_error({test <- assocTest(
    GT_test,
    test = "glm",
    covar = paste0("PC", 1:4),
    pheno = "PC1",
    verbose = FALSE
  )}, regexp = "should be coded")

  # only one phenotype can be specified
  # expect_error({test <- assocTest(
  #   GT_test,
  #   test = "glm",
  #   covar = paste0("PC", 1:4),
  #   pheno = c("PC1", "PC2"),
  #   verbose = FALSE
  # )}, regexp = "Only one phenotype")

  # keep list
  keep <- colnames(GT_test)[1:100]
  expect_true(any(stringr::str_detect( capture_messages(
  {
     assocTest(
      GT_test,
      keep = keep,
      test = "glm",
      covar = paste0("PC", 1:4),
      pheno = "pheno",
      verbose = TRUE
    )
    }
  ), sprintf("100/%s", ncol(GT_test)))))

  # expect_error({test <- assocTest(
  #   GT_test,
  #   test = "glm",
  #   covar = paste0("PC", 1:4),
  #   pheno = c("PC1", "PC2"),
  #   verbose = FALSE
  # )}, regexp = "Only one phenotype")


  # expect message when phenotype contains missing values (and verbose = TRUE)
  # tmp <- GT_test
  # tmp$pheno[1:5]<- NA
  # expect_true(any(stringr::str_detect(capture_messages(
  # {
  #   test <- assocTest(
  #       tmp,
  #       test = "glm",
  #       covar = c(paste0("PC", 1:4)),
  #       pheno = "pheno",
  #       verbose = TRUE
  #     )
  #   }
  # ), "have missing phenotype")))
  # rm(tmp)

  # expect warning when missing weights are included
  tmp <- GT_test
  rowData(tmp)$w[1:5] <- NA_real_
  expect_warning({test <- assocTest(
    tmp,
    test = c("glm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    verbose = FALSE
  )}, regexp = "are missing")

  # expect warning when negative weights are included
  tmp <- GT_test
  rowData(tmp)$w[1:5] <- -1
  expect_warning({test <- assocTest(
    tmp,
    test = c("glm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    verbose = FALSE
  )}, regexp = "are < 0")
  rm(tmp)

  # expect empty results when no variants pass thresholds
  test <- assocTest(
    GT_test,
    test = c("glm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    minMAF = 0.011,
    maxMAF = 0.011,
    verbose = FALSE
  )
  expect_true(nrow(test) == 0 && is(test, "rvatResult"))
  
  # expect empty results when no samples pass thresholds
  test <- assocTest(
    GT_test,
    test = c("glm"),
    covar = c(paste0("PC", 1:4)),
    minCallrateSM = 0.9,
    maxCallrateSM = 0.9,
    pheno = "pheno",
    verbose = FALSE
  )
  expect_true(nrow(test) == 0 && is(test, "rvatResult"))

  # expect an error when a non-implemented test is specified
  expect_error({test <- assocTest(
    GT_test,
    test = c("glm", "test"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    verbose = FALSE
  )}, regexp = "The following tests are not valid: test")
  
  # expect a warning when tests are specified that are only valid for binary tests or vice versa
  expect_warning({test <- assocTest(
    GT_test,
    test = c("glm", "lm"),
    covar = c(paste0("PC", 2:4)),
    pheno = "PC1",
    verbose = FALSE,
    continuous = TRUE
  )}, regexp = "The following tests were excluded since they are not available for continuous rvb tests: glm")
  
  # expect an error when a non-existing phenotype is specified
  expect_error({test <- assocTest(
    GT_test,
    test = c("glm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "test",
    verbose = FALSE
  )}, regexp = "'test' is not present in `colData\\(GT\\)`")

  # expect an error when overwriteAggregate is FALSE, but no aggregate field is present
  # expect_error({test <- assocTest(
  #   GT_test,
  #   test = c("glm"),
  #   covar = c(paste0("PC", 1:4)),
  #   pheno = "pheno",
  #   overwriteAggregate = FALSE,
  #   verbose = FALSE
  # )}, regexp = "no 'aggregate' column")

  # expect error the genoMatrix is set to dominant/recessive 
  # and another geneticModel is specified
  expect_error({test <- assocTest(
    recode(GT_test, geneticModel = "recessive"),
    test = c("glm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    geneticModel = "dominant",
    verbose = FALSE
  )}, regexp = "in order to apply")

  # expect warning the genoMatrix is set to dominant/recessive 
  # and 'allelic' geneticModel is specified
  expect_error({test <- assocTest(
    recode(GT_test, geneticModel = "recessive"),
    test = c("glm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    geneticModel = "allelic",
    verbose = FALSE
  )}, regexp = "Current geneticModel should be")


  # expect error when MAC/MAF filters are set while geneticModel != allelic
  expect_error({test <- assocTest(
    recode(GT_test, geneticModel = "recessive"),
    test = c("glm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    geneticModel = "dominant",
    maxMAC = 10,
    verbose = FALSE
  )}, regexp = "do not apply")

  # export error when non-implemented imputeMethod is specified
  expect_error({test <- assocTest(
    GT_test,
    test = c("glm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    imputeMethod = "hello!",
    verbose = FALSE
  )}, regexp = "should be either") 

  # warning on non-implemented tests
  expect_warning({test <- assocTest(
    GT_test,
    test = c("glm", "lm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    verbose = FALSE
  )}, regexp = "The following tests were excluded")

  expect_warning({test <- assocTest(
    GT_test,
    test = c("glm", "lm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    singlevar = TRUE,
    verbose = FALSE
  )}, regexp = "The following tests were excluded")

  expect_warning({test <- assocTest(
    GT_test,
    test = c("glm", "lm"),
    covar = c(paste0("PC", 2:4)),
    pheno = "PC1",
    continuous = TRUE,
    verbose = FALSE
  )}, regexp = "The following tests were excluded")
  
  expect_warning({test <- assocTest(
    GT_test,
    test = c("glm", "lm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    singlevar = TRUE,
    continuous = TRUE,
    verbose = FALSE
  )}, regexp = "The following tests were excluded")

  # errror if no applicable tests are left
  expect_error({suppressWarnings(test <- assocTest(
    GT_test,
    test = c("lm"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    verbose = FALSE
  ))}, regexp = "No applicable tests")

  # only one offset can be specified
  expect_error({test <- assocTest(
    GT_test,
    test = "firth",
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    offset = c("PC", 5:6),
    verbose = FALSE
  )}, regexp = "Currently at most")
  
  # offset should be available
  expect_error({test <- assocTest(
    GT_test,
    test = "firth",
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    offset = "hello!",
    verbose = FALSE
  )}, regexp = "not available")

  # warning 
  expect_error({test <- assocTest(
    GT_test,
    test = "firth",
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    offset = "hello!",
    verbose = FALSE
  )}, regexp = "not available")

  # resampling

  ## not implemented for singlevar tetss
  expect_error({test <- assocTest(
    GT_test,
    test = "firth",
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    methodResampling = "permutation",
    verbose = FALSE,
    singlevar = TRUE
  )}, regexp = "not implemented")

  ## throw warning if some tests are not valid 
  expect_warning({test <- assocTest(
    GT_test,
    test = c("skat", "firth"),
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    methodResampling = "permutation",
    nResampling = 2,
    verbose = FALSE
  )}, regexp = "not implemented")

  ## throw error if no valid tests are implemented
  expect_error({test <- assocTest(
    GT_test,
    test = "firth",
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    methodResampling = "permutation",
    verbose = FALSE
  )}, regexp = "not implemented")

  ## throw error if no valid permutation methods are specified
  expect_error({test <- assocTest(
    GT_test,
    test = "skat",
    covar = c(paste0("PC", 1:4)),
    pheno = "pheno",
    methodResampling = "hello!",
    verbose = FALSE
  )}, regexp = "accepted option")

}) 
