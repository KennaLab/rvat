compare_gdbs <- function(
  gdb1,
  gdb2,
  cohort1 = "SM",
  cohort2 = "SM",
  check_var = TRUE,
  check_tables = TRUE
) {
  # compare var tables
  if (check_var) {
    anno1 <- getAnno(gdb1, table = "var")
    anno2 <- getAnno(gdb2, table = "var")
    expect_identical(anno1, anno2)
  }

  # check if same tables are present
  if (check_tables) {
    tables1 <- DBI::dbListTables(gdb1)
    tables2 <- DBI::dbListTables(gdb2)
    expect_identical(sort(tables1), sort(tables2))
  }

  # compare GTs
  GT1 <- getGT(
    gdb1,
    cohort = cohort1,
    VAR_id = getAnno(gdb1, table = "var", fields = "VAR_id")$VAR_id,
    verbose = FALSE
  )
  GT2 <- getGT(
    gdb2,
    cohort = cohort2,
    VAR_id = getAnno(gdb2, table = "var", fields = "VAR_id")$VAR_id,
    verbose = FALSE
  )
  ## ignore gdb and gdbId in comparison
  metadata(GT1)$gdb <- NA_character_
  metadata(GT2)$gdb <- NA_character_
  metadata(GT1)$gdbId <- NA_character_
  metadata(GT2)$gdbId <- NA_character_

  expect_identical(GT1, GT2)

  invisible(NULL)
}


expect_gdb_indexes <- function(
  gdb,
  expected_indexes = c("var_idx", "var_idx2", "SM_idx", "dosage_idx")
) {
  indexes <- DBI::dbGetQuery(
    gdb,
    "SELECT name FROM sqlite_master WHERE type='index'"
  )
  expect_equal(sort(expected_indexes), sort(indexes$name))
}

compare_varsetfile <- function(vsfile1, vsfile2, units = NULL) {
  if (is.null(units)) {
    units <- listUnits(vsfile1)
  }
  vs1 <- getVarSet(vsfile1, unit = units)
  vs2 <- getVarSet(vsfile2, unit = units)
  vs1@metadata <- list()
  vs2@metadata <- list()
  expect_identical(vs1, vs2)
}

compare_varsets <- function(varset1, varset2) {
  varset1@metadata <- list()
  varset2@metadata <- list()
  expect_identical(varset1, varset2)
}

create_merged_varset_file <- function() {
  gdb <- gdb(rvat_example("rvatData.gdb"))

  # Build varsets
  suppressMessages(buildVarSet(
    object = gdb,
    output = setup$moderate_file,
    varSetName = "Moderate",
    unitTable = "varInfo",
    unitName = "gene_name",
    where = paste(
      "(ModerateImpact = 1 or HighImpact = 1) and",
      setup$where_clause
    )
  ))

  suppressMessages(buildVarSet(
    object = gdb,
    output = setup$lof_file,
    varSetName = "HighImpact",
    unitTable = "varInfo",
    unitName = "gene_name",
    where = paste("HighImpact = 1 and", setup$where_clause),
    verbose = FALSE
  ))

  suppressMessages(buildVarSet(
    object = gdb,
    output = setup$cadd_file,
    varSetName = "CADD",
    unitTable = "varInfo",
    unitName = "gene_name",
    weightName = "CADDphred",
    where = setup$where_clause
  ))

  # Combine varsets
  varsets <- list()
  for (i in seq_along(c(
    setup$moderate_file,
    setup$lof_file,
    setup$cadd_file
  ))) {
    file <- c(setup$moderate_file, setup$lof_file, setup$cadd_file)[i]
    varset <- readr::read_delim(
      file,
      delim = "|",
      col_names = FALSE,
      col_types = "cccc",
      comment = "#"
    )
    varsets[[i]] <- varset
  }

  varsets <- dplyr::bind_rows(varsets) |> dplyr::arrange(X1)
  metadata_obj <- metadata(varSetFile(setup$moderate_file))

  merged_file <- withr::local_tempfile()
  con <- gzfile(merged_file, "w")
  rvat:::.write_rvat_header(
    filetype = "varSetFile",
    metadata = metadata_obj,
    con = con
  )
  write.table(
    varsets,
    file = con,
    col.names = FALSE,
    sep = "|",
    append = TRUE,
    quote = FALSE,
    row.names = FALSE
  )
  close(con)

  merged_file
}

.get_param <- function(param_list, name, default = NULL) {
  if (!is.null(param_list[[name]])) param_list[[name]] else default
}

generate_param <- function(
  pheno_bin = "pheno",
  pheno_quant = "age",
  covar = paste0("PC", 1:4),
  covarCat = "superPop",
  tests_bin_rvb = c(
    "firth",
    "scoreSPA",
    "glm",
    "skat",
    "skat_burden",
    "skato",
    "skat_robust",
    "skat_burden_robust",
    "skato_robust",
    "acatv",
    "acatvSPA"
  ),
  tests_cont_rvb = c("lm", "skat", "skat_burden", "skato", "acatv"),
  tests_bin_sv = c("firth", "glm", "scoreSPA"),
  tests_cont_sv = c("lm")
) {
  params <- list()
  ## rvb
  for (mode in c("binary", "continuous")) {
    if (mode == "binary") {
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
        test = if (mode == "continuous") tests else c(tests, "acatvfirth"),
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
        minMAC = 2,
        maxMAC = 50,
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
        covar = c(covar, covarCat)
      )
    )
    params[["rvb"]][[mode]] <- param
  }

  ## sv
  for (mode in c("binary", "continuous")) {
    if (mode == "binary") {
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
        minMAC = 2,
        maxMAC = 50,
        minCarrierFreq = 0.00001
      ),
      # couple of filters, include factor covar
      list(
        pheno = pheno,
        test = tests,
        maxCarrierFreq = 0.001,
        covar = c(covar, covarCat)
      )
    )
    params[["sv"]][[mode]] <- param
  }
  params
}

run_tests_gdb <- function(
  gdb,
  params,
  VAR_id,
  var_sv = NULL,
  cohort = "pheno"
) {
  results <- list()

  ## binary burden tests
  results_bin <- list()
  for (i in 1:length(params[["rvb"]][["binary"]])) {
    # print(sprintf("rvb binary: %s/%s", i, length(params[["rvb"]][["binary"]])))
    param_i <- params[["rvb"]][["binary"]][[i]]
    assoc <- suppressMessages(assocTest(
      object = gdb,
      VAR_id = VAR_id,
      cohort = cohort,
      pheno = param_i$pheno,
      test = param_i$test,
      covar = .get_param(param_i, "covar"),
      imputeMethod = .get_param(param_i, "imputeMethod"),
      geneticModel = .get_param(param_i, "geneticModel", "allelic"),
      MAFweights = .get_param(param_i, "MAFweights", "none"),
      minCallrateVar = .get_param(param_i, "minCallrateVar", 0),
      minCallrateSM = .get_param(param_i, "minCallrateSM", 0),
      minMAF = .get_param(param_i, "minMAF", 0),
      maxMAF = .get_param(param_i, "maxMAF", 1),
      minMAC = .get_param(param_i, "minMAC", 0),
      maxMAC = .get_param(param_i, "maxMAC", Inf),
      minCarriers = .get_param(param_i, "minCarriers", 0),
      maxCarriers = .get_param(param_i, "maxCarriers", Inf),
      minCarrierFreq = .get_param(param_i, "minCarrierFreq", 0),
      maxCarrierFreq = .get_param(param_i, "maxCarrierFreq", 1),
      verbose = TRUE
    ))
    assoc$i <- as.character(i)
    results_bin[[as.character(i)]] <- assoc
  }
  results[["rvb"]][["binary"]] <- as.data.frame(do.call(rbind, results_bin))

  ## quantitative burden tests
  results_cont <- list()
  for (i in 1:length(params[["rvb"]][["continuous"]])) {
    #print(sprintf("rvb continuous: %s/%s", i, length(params[["rvb"]][["continuous"]])))
    param_i <- params[["rvb"]][["continuous"]][[i]]
    assoc <- suppressMessages(assocTest(
      object = gdb,
      VAR_id = VAR_id,
      cohort = cohort,
      pheno = param_i$pheno,
      test = param_i$test,
      continuous = TRUE,
      covar = .get_param(param_i, "covar"),
      imputeMethod = .get_param(param_i, "imputeMethod"),
      geneticModel = .get_param(param_i, "geneticModel", "allelic"),
      MAFweights = .get_param(param_i, "MAFweights", "none"),
      minCallrateVar = .get_param(param_i, "minCallrateVar", 0),
      minCallrateSM = .get_param(param_i, "minCallrateSM", 0),
      minMAF = .get_param(param_i, "minMAF", 0),
      maxMAF = .get_param(param_i, "maxMAF", 1),
      minMAC = .get_param(param_i, "minMAC", 0),
      maxMAC = .get_param(param_i, "maxMAC", Inf),
      minCarriers = .get_param(param_i, "minCarriers", 0),
      maxCarriers = .get_param(param_i, "maxCarriers", Inf),
      minCarrierFreq = .get_param(param_i, "minCarrierFreq", 0),
      maxCarrierFreq = .get_param(param_i, "maxCarrierFreq", 1),
      verbose = TRUE
    ))
    assoc$i <- as.character(i)
    results_cont[[as.character(i)]] <- assoc
  }
  results[["rvb"]][["continuous"]] <- as.data.frame(do.call(
    rbind,
    results_cont
  ))

  ## binary sv tests
  results_bin <- list()
  for (i in 1:length(params[["sv"]][["binary"]])) {
    #print(sprintf("sv binary: %s/%s", i, length(params[["sv"]][["binary"]])))
    param_i <- params[["sv"]][["binary"]][[i]]
    assoc <- suppressMessages(assocTest(
      object = gdb,
      VAR_id = var_sv,
      cohort = cohort,
      pheno = param_i$pheno,
      test = param_i$test,
      covar = .get_param(param_i, "covar"),
      geneticModel = .get_param(param_i, "geneticModel", "allelic"),
      minCallrateVar = .get_param(param_i, "minCallrateVar", 0),
      minCallrateSM = .get_param(param_i, "minCallrateSM", 0),
      minMAF = .get_param(param_i, "minMAF", 0),
      maxMAF = .get_param(param_i, "maxMAF", 1),
      minMAC = .get_param(param_i, "minMAC", 0),
      maxMAC = .get_param(param_i, "maxMAC", Inf),
      minCarriers = .get_param(param_i, "minCarriers", 0),
      maxCarriers = .get_param(param_i, "maxCarriers", Inf),
      minCarrierFreq = .get_param(param_i, "minCarrierFreq", 0),
      maxCarrierFreq = .get_param(param_i, "maxCarrierFreq", 1),
      verbose = TRUE,
      singlevar = TRUE
    ))

    if (nrow(assoc) > 0) {
      assoc$i <- as.character(i)
      results_bin[[as.character(i)]] <- assoc
    } else {
      results_bin[[as.character(i)]] <- NULL
    }
    results_bin[[as.character(i)]] <- assoc
  }
  results[["sv"]][["binary"]] <- as.data.frame(do.call(rbind, results_bin))

  ## continuous sv tests
  results_cont <- list()
  for (i in 1:length(params[["sv"]][["continuous"]])) {
    #print(sprintf("sv continuous: %s/%s", i, length(params[["sv"]][["continuous"]])))
    param_i <- params[["sv"]][["continuous"]][[i]]
    assoc <- suppressMessages(assocTest(
      object = gdb,
      VAR_id = var_sv,
      cohort = cohort,
      pheno = param_i$pheno,
      test = param_i$test,
      covar = .get_param(param_i, "covar"),
      geneticModel = .get_param(param_i, "geneticModel", "allelic"),
      minCallrateVar = .get_param(param_i, "minCallrateVar", 0),
      minCallrateSM = .get_param(param_i, "minCallrateSM", 0),
      minMAF = .get_param(param_i, "minMAF", 0),
      maxMAF = .get_param(param_i, "maxMAF", 1),
      minMAC = .get_param(param_i, "minMAC", 0),
      maxMAC = .get_param(param_i, "maxMAC", Inf),
      minCarriers = .get_param(param_i, "minCarriers", 0),
      maxCarriers = .get_param(param_i, "maxCarriers", Inf),
      minCarrierFreq = .get_param(param_i, "minCarrierFreq", 0),
      maxCarrierFreq = .get_param(param_i, "maxCarrierFreq", 1),
      verbose = TRUE,
      singlevar = TRUE,
      continuous = TRUE
    ))

    if (nrow(assoc) > 0) {
      assoc$i <- as.character(i)
      results_cont[[as.character(i)]] <- assoc
    } else {
      results_cont[[as.character(i)]] <- NULL
    }
  }
  results[["sv"]][["continuous"]] <- as.data.frame(do.call(rbind, results_cont))
  results
}

run_tests_GT <- function(GT, params, var_sv = NULL) {
  results <- list()

  ##
  if (is.null(var_sv)) {
    var_sv <- rownames(GT)
  }

  ## binary burden tests
  results_bin <- list()
  for (i in 1:length(params[["rvb"]][["binary"]])) {
    #print(sprintf("rvb binary: %s/%s", i, length(params[["rvb"]][["binary"]])))
    param_i <- params[["rvb"]][["binary"]][[i]]
    assoc <- suppressMessages(assocTest(
      object = GT,
      pheno = param_i$pheno,
      test = param_i$test,
      covar = .get_param(param_i, "covar"),
      imputeMethod = .get_param(param_i, "imputeMethod"),
      geneticModel = .get_param(param_i, "geneticModel", "allelic"),
      MAFweights = .get_param(param_i, "MAFweights", "none"),
      minCallrateVar = .get_param(param_i, "minCallrateVar", 0),
      minCallrateSM = .get_param(param_i, "minCallrateSM", 0),
      minMAF = .get_param(param_i, "minMAF", 0),
      maxMAF = .get_param(param_i, "maxMAF", 1),
      minMAC = .get_param(param_i, "minMAC", 0),
      maxMAC = .get_param(param_i, "maxMAC", Inf),
      minCarriers = .get_param(param_i, "minCarriers", 0),
      maxCarriers = .get_param(param_i, "maxCarriers", Inf),
      minCarrierFreq = .get_param(param_i, "minCarrierFreq", 0),
      maxCarrierFreq = .get_param(param_i, "maxCarrierFreq", 1),
      verbose = TRUE
    ))
    assoc$i <- as.character(i)
    metadata(assoc)$creationDate <- NA_character_
    results_bin[[as.character(i)]] <- assoc
  }
  results[["rvb"]][["binary"]] <- as.data.frame(do.call(rbind, results_bin))

  ## quantitative burden tests
  results_cont <- list()
  for (i in 1:length(params[["rvb"]][["continuous"]])) {
    #print(sprintf("rvb continuous: %s/%s", i, length(params[["rvb"]][["continuous"]])))
    param_i <- params[["rvb"]][["continuous"]][[i]]
    assoc <- suppressMessages(assocTest(
      object = GT,
      pheno = param_i$pheno,
      test = param_i$test,
      continuous = TRUE,
      covar = .get_param(param_i, "covar"),
      imputeMethod = .get_param(param_i, "imputeMethod"),
      geneticModel = .get_param(param_i, "geneticModel", "allelic"),
      MAFweights = .get_param(param_i, "MAFweights", "none"),
      minCallrateVar = .get_param(param_i, "minCallrateVar", 0),
      minCallrateSM = .get_param(param_i, "minCallrateSM", 0),
      minMAF = .get_param(param_i, "minMAF", 0),
      maxMAF = .get_param(param_i, "maxMAF", 1),
      minMAC = .get_param(param_i, "minMAC", 0),
      maxMAC = .get_param(param_i, "maxMAC", Inf),
      minCarriers = .get_param(param_i, "minCarriers", 0),
      maxCarriers = .get_param(param_i, "maxCarriers", Inf),
      minCarrierFreq = .get_param(param_i, "minCarrierFreq", 0),
      maxCarrierFreq = .get_param(param_i, "maxCarrierFreq", 1),
      verbose = TRUE
    ))
    assoc$i <- as.character(i)
    metadata(assoc)$creationDate <- NA_character_
    results_cont[[as.character(i)]] <- assoc
  }
  results[["rvb"]][["continuous"]] <- as.data.frame(do.call(
    rbind,
    results_cont
  ))

  ## binary sv tests
  results_bin <- list()
  for (i in 1:length(params[["sv"]][["binary"]])) {
    #print(sprintf("sv binary: %s/%s", i, length(params[["sv"]][["binary"]])))
    param_i <- params[["sv"]][["binary"]][[i]]
    assoc <- suppressMessages(assocTest(
      object = GT[var_sv, ],
      pheno = param_i$pheno,
      test = param_i$test,
      covar = .get_param(param_i, "covar"),
      geneticModel = .get_param(param_i, "geneticModel", "allelic"),
      minCallrateVar = .get_param(param_i, "minCallrateVar", 0),
      minCallrateSM = .get_param(param_i, "minCallrateSM", 0),
      minMAF = .get_param(param_i, "minMAF", 0),
      maxMAF = .get_param(param_i, "maxMAF", 1),
      minMAC = .get_param(param_i, "minMAC", 0),
      maxMAC = .get_param(param_i, "maxMAC", Inf),
      minCarriers = .get_param(param_i, "minCarriers", 0),
      maxCarriers = .get_param(param_i, "maxCarriers", Inf),
      minCarrierFreq = .get_param(param_i, "minCarrierFreq", 0),
      maxCarrierFreq = .get_param(param_i, "maxCarrierFreq", 1),
      verbose = TRUE,
      singlevar = TRUE
    ))

    if (nrow(assoc) > 0) {
      assoc$i <- as.character(i)
      metadata(assoc)$creationDate <- NA_character_
      results_bin[[as.character(i)]] <- assoc
    } else {
      results_bin[[as.character(i)]] <- NULL
    }
    results_bin[[as.character(i)]] <- assoc
  }
  results[["sv"]][["binary"]] <- as.data.frame(do.call(rbind, results_bin))

  ## continuous sv tests
  results_cont <- list()
  for (i in 1:length(params[["sv"]][["continuous"]])) {
    #print(sprintf("sv continuous: %s/%s", i, length(params[["sv"]][["continuous"]])))
    param_i <- params[["sv"]][["continuous"]][[i]]
    assoc <- suppressMessages(assocTest(
      object = GT[var_sv, ],
      pheno = param_i$pheno,
      test = param_i$test,
      covar = .get_param(param_i, "covar"),
      geneticModel = .get_param(param_i, "geneticModel", "allelic"),
      minCallrateVar = .get_param(param_i, "minCallrateVar", 0),
      minCallrateSM = .get_param(param_i, "minCallrateSM", 0),
      minMAF = .get_param(param_i, "minMAF", 0),
      maxMAF = .get_param(param_i, "maxMAF", 1),
      minMAC = .get_param(param_i, "minMAC", 0),
      maxMAC = .get_param(param_i, "maxMAC", Inf),
      minCarriers = .get_param(param_i, "minCarriers", 0),
      maxCarriers = .get_param(param_i, "maxCarriers", Inf),
      minCarrierFreq = .get_param(param_i, "minCarrierFreq", 0),
      maxCarrierFreq = .get_param(param_i, "maxCarrierFreq", 1),
      verbose = TRUE,
      singlevar = TRUE,
      continuous = TRUE
    ))

    if (nrow(assoc) > 0) {
      assoc$i <- as.character(i)
      metadata(assoc)$creationDate <- NA_character_
      results_cont[[as.character(i)]] <- assoc
    } else {
      results_cont[[as.character(i)]] <- NULL
    }
  }
  results[["sv"]][["continuous"]] <- as.data.frame(do.call(rbind, results_cont))
  results
}
