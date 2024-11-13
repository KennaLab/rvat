# rationale
# - circumvents need of rerunning all functions for CLI testing, while they're all ready being tested separately
# - directly checks if arguments are passed correctly, rather than indirectly by checking outputs etc.
# - tests will directly fail if parameters or defaults ihave not been properly updated in CLI

# buildGdb
test_that("--buildGdb works", {
  # generate mocked function/method that returns the environment
  forms <- formals(buildGdb)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    buildGdb = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--buildGdb",
        "--vcf=example.vcf.gz",
        "--output=example.gdb"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    vcf = "example.vcf.gz",
    output = "example.gdb"
  )
  expected_args <- c(expected_args, formals(buildGdb)[!names(formals(buildGdb)) %in% names(expected_args)])
  expect_equal(expected_args, mock_args)
  
  # test non-defaults
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--buildGdb",
        "--vcf=example.vcf.gz",
        "--output=example.gdb",
        "--skipIndexes",
        "--skipVarRanges",
        "--overWrite",
        "--genomeBuild=GRCh38",
        "--memlimit=1000",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    vcf = "example.vcf.gz",
    output = "example.gdb",
    skipIndexes = TRUE,
    skipVarRanges = TRUE,
    overWrite = TRUE,
    genomeBuild = "GRCh38",
    memlimit = 1000,
    verbose = FALSE
  )
  expected_args <- c(expected_args, formals(buildGdb)[!names(formals(buildGdb)) %in% names(expected_args)])
  expect_equal(expected_args, mock_args)
  
})


# concatGdb
test_that("--concatGdb works", {
  # generate mocked function/method that returns the environment
  forms <- formals(concatGdb)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    concatGdb = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--concatGdb",
        "--targets=targets.txt",
        "--output=example.gdb"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    targets = "targets.txt",
    output = "example.gdb"
  )
  expected_args <- c(expected_args, formals(concatGdb)[!names(formals(concatGdb)) %in% names(expected_args)])
  expect_equal(expected_args, mock_args)
  
  # test non-defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--concatGdb",
        "--targets=targets.txt",
        "--output=example.gdb",
        "--skipRemap",
        "--skipIndexes",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    targets = "targets.txt",
    output = "example.gdb",
    skipRemap = TRUE,
    skipIndexes = TRUE,
    verbose = FALSE
  )
  expected_args <- c(expected_args, formals(concatGdb)[!names(formals(concatGdb)) %in% names(expected_args)])
  expect_equal(expected_args, mock_args)
})


# subsetGdb
test_that("--subsetGdb works", {
  # generate mocked function/method that returns the environment
  forms <- formals(subsetGdb)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    subsetGdb = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--subsetGdb",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--output=example_subset.gdb"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    output = "example_subset.gdb"
  )
  expected_args <- c(expected_args, formals(subsetGdb)[!names(formals(subsetGdb)) %in% names(expected_args)])
  expect_s4_class(mock_args[[1]], "gdb")
  expect_equal(expected_args[2:length(expected_args)], mock_args[2:length(mock_args)])
  
  # test non-defaults 
  tmpfile <- withr::local_tempfile()
  readr::write_lines(1:10, tmpfile)
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--subsetGdb",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--output=example.gdb",
        "--intersection=varInfo",
        "--where='ModerateImpact = 1'",
        sprintf("--VAR_id=%s", tmpfile),
        "--tables=varInfo",
        "--skipIndexes",
        "--overWrite",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    output = "example.gdb",
    intersection = "varInfo",
    where = "'ModerateImpact = 1'",
    VAR_id = as.character(1:10),
    tables = "varInfo",
    skipIndexes = TRUE,
    overWrite = TRUE,
    verbose = FALSE
  )
  expect_s4_class(mock_args[[1]], "gdb")
  expected_args <- c(expected_args, formals(subsetGdb)[!names(formals(subsetGdb)) %in% names(expected_args)])
  expect_equal(expected_args[2:length(expected_args)], mock_args[2:length(mock_args)])
  
})

# uploadAnno
test_that("--uploadAnno works", {
  # generate mocked function/method that returns the environment
  forms <- formals(uploadAnno)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    uploadAnno = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--uploadAnno",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--name=data",
        "--value=data.txt"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    name = "data",
    value = "data.txt"
  )
  expected_args <- c(expected_args, formals(uploadAnno)[!names(formals(uploadAnno)) %in% names(expected_args)])
  expect_s4_class(mock_args[[1]], "gdb")
  expect_equal(expected_args[2:length(expected_args)], mock_args[2:length(mock_args)])
  
  # test non-defaults
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--uploadAnno",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--name=data",
        "--value=data.txt",
        "--sep='|'",
        "--skipRemap",
        "--skipIndexes",
        "--ignoreAlleles",
        "--keepUnmapped",
        "--mapRef='var2'",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    name = "data",
    value = "data.txt",
    sep = "'|'",
    skipRemap = TRUE,
    skipIndexes = TRUE,
    ignoreAlleles = TRUE,
    keepUnmapped = TRUE,
    mapRef = "'var2'",
    verbose = FALSE
  )
  expect_s4_class(mock_args[[1]], "gdb")
  expected_args <- c(expected_args, formals(uploadAnno)[!names(formals(uploadAnno)) %in% names(expected_args)])
  expect_equal(expected_args[2:length(expected_args)], mock_args[2:length(mock_args)])
  
})

# uploadCohort
test_that("--uploadCohort works", {
  # generate mocked function/method that returns the environment
  forms <- formals(uploadCohort)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    uploadCohort = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--uploadCohort",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--name=data",
        "--value=data.txt"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    name = "data",
    value = "data.txt"
  )
  expected_args <- c(expected_args, formals(uploadCohort)[!names(formals(uploadCohort)) %in% names(expected_args)])
  expect_s4_class(mock_args[[1]], "gdb")
  expect_equal(expected_args[2:length(expected_args)], mock_args[2:length(mock_args)])
  
  # test non-defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--uploadCohort",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--name=data",
        "--value=data.txt",
        "--sep='|'",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    name = "data",
    value = "data.txt",
    sep = "'|'",
    verbose = FALSE
  )
  expect_s4_class(mock_args[[1]], "gdb")
  expected_args <- c(expected_args, formals(uploadCohort)[!names(formals(uploadCohort)) %in% names(expected_args)])
  expect_equal(expected_args[2:length(expected_args)], mock_args[2:length(mock_args)])
  
})


# listAnno
test_that("--listAnno works", {
  # generate mocked function/method that returns the environment
  forms <- formals(listAnno)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    listAnno = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--listAnno",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb"))
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb"))
  )
  expected_args <- c(expected_args, formals(listAnno)[!names(formals(listAnno)) %in% names(expected_args)])
  expect_s4_class(mock_args[[1]], "gdb")
})



# listCohort
test_that("--listCohort works", {
  # generate mocked function/method that returns the environment
  forms <- formals(listCohort)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    listCohort = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--listCohort",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb"))
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb"))
  )
  expected_args <- c(expected_args, formals(listCohort)[!names(formals(listCohort)) %in% names(expected_args)])
  expect_s4_class(mock_args[[1]], "gdb")
})

# mapVariants
test_that("--mapVariants works", {
  # generate mocked function/method that returns the environment
  forms <- formals(mapVariants)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    mapVariants = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--mapVariants",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb"))
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb"))
  )
  expected_args <- c(expected_args, formals(mapVariants)[!names(formals(mapVariants)) %in% names(expected_args)])
  expect_s4_class(mock_args[[1]], "gdb")
  
  # skip bedCols
  expect_equal(expected_args[!names(expected_args) %in% c("object", "bedCols")], 
               mock_args[!names(mock_args) %in% c("object", "bedCols")])
  
  
  # test non-defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--mapVariants",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--ranges=ranges.txt",
        "--gff=ensembl.gtf",
        "--bed=bedfile.bed",
        "--bedCols=gene,gene_name",
        "--fields=gene,gene_name",
        "--uploadName=ensembl",
        "--output=output.txt.gz",
        "--sep='|'",
        "--skipIndexes",
        "--overWrite",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    ranges = "ranges.txt",
    gff = "ensembl.gtf",
    bed = "bedfile.bed",
    bedCols = c("gene", "gene_name"),
    fields = c("gene", "gene_name"),
    uploadName = "ensembl",
    output = "output.txt.gz",
    sep = "'|'",
    skipIndexes = TRUE,
    overWrite = TRUE,
    verbose = FALSE
  )
  expect_s4_class(mock_args[[1]], "gdb")
  expected_args <- c(expected_args, formals(mapVariants)[!names(formals(mapVariants)) %in% names(expected_args)])
  expect_equal(expected_args[2:length(expected_args)], mock_args[2:length(mock_args)])
  
})


# buildVarSet
test_that("--buildVarSet works", {
  # generate mocked function/method that returns the environment
  forms <- rvat:::.retrieve_formals_dispatch("buildVarSet", "gdb")
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    buildVarSet = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--buildVarSet",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--varSetName=LOF",
        "--unitTable=gene",
        "--unitName=gene_id",
        "--output=varsetfile.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    varSetName = "LOF",
    unitTable = "gene",
    unitName = "gene_id",
    output = "varsetfile.txt.gz"
  )
  expected_args <- c(expected_args, formals(buildVarSet)[!names(formals(buildVarSet)) %in% names(expected_args)])
  expect_s4_class(mock_args[[1]], "gdb")
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("object")], 
               mock_args[!names(mock_args) %in% c("object")])
  
  
  # test non-defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--buildVarSet",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--varSetName=LOF",
        "--unitTable=gene",
        "--unitName=gene_id",
        "--output=varsetfile.txt.gz",
        "--intersection=QCpass",
        "--where='QC=1'",
        "--weightName=CADD",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    varSetName = "LOF",
    unitTable = "gene",
    unitName = "gene_id",
    output = "varsetfile.txt.gz",
    intersection = "QCpass",
    where = "'QC=1'",
    weightName = "CADD",
    verbose = FALSE
  )
  expect_s4_class(mock_args[[1]], "gdb")
  expected_args <- c(expected_args, formals(buildVarSet)[!names(formals(buildVarSet)) %in% names(expected_args)])
  expect_equal(expected_args[!names(expected_args) %in% c("object")], 
               mock_args[!names(mock_args) %in% c("object")])
  
})


# spatialClust
test_that("--spatialClust works", {
  # generate mocked function/method that returns the environment
  forms <- formals(spatialClust)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    spatialClust = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--spatialClust",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--output=varsetfile.txt.gz",
        "--varSetName=LOF",
        "--unitTable=gene",
        "--unitName=gene_id",
        "--windowSize=60",
        "--overlap=30"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    varSetName = "LOF",
    unitTable = "gene",
    unitName = "gene_id",
    output = "varsetfile.txt.gz",
    windowSize = 60,
    overlap = 30
  )
  expected_args <- c(expected_args, formals(spatialClust)[!names(formals(spatialClust)) %in% names(expected_args)])
  expect_s4_class(mock_args[["object"]], "gdb")
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("object")], 
               mock_args[!names(mock_args) %in% c("object")])
  
  
  # test non-defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--spatialClust",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--varSetName=LOF",
        "--unitTable=gene",
        "--unitName=gene_id",
        "--output=varsetfile.txt.gz",
        "--intersection=QCpass",
        "--where='QC=1'",
        "--weightName=CADD",
        "--windowSize=60,100",
        "--overlap=30,50",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    varSetName = "LOF",
    unitTable = "gene",
    unitName = "gene_id",
    output = "varsetfile.txt.gz",
    intersection = "QCpass",
    where = "'QC=1'",
    weightName = "CADD",
    windowSize = c(60, 100),
    overlap = c(30, 50)
  )
  expect_s4_class(mock_args[["object"]], "gdb")
  expected_args <- c(expected_args, formals(spatialClust)[!names(formals(spatialClust)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("object")], 
               mock_args[!names(mock_args) %in% c("object")])
  
})


# summariseGeno
test_that("--summariseGeno works", {
  # generate mocked function/method that returns the environment
  forms <- rvat:::.retrieve_formals_dispatch("summariseGeno", "gdb")
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    summariseGeno = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--summariseGeno",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--output=sumgeno.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    output = "sumgeno.txt.gz"
  )
  expected_args <- c(expected_args, formals(summariseGeno)[!names(formals(summariseGeno)) %in% names(expected_args)])
  expect_s4_class(mock_args[[1]], "gdb")
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  
  expect_equal(expected_args[!names(expected_args) %in% c("object")], 
               mock_args[!names(mock_args) %in% c("object")])
  
  
  # test non-defaults 
  ## save varids
  varids <- withr::local_tempfile()
  readr::write_lines(1:10, varids)
  
  ## save varsetfile
  varsetfile <- withr::local_tempfile()
  buildVarSet(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    varSetName = "Moderate",
    unitTable = "varInfo",
    unitName = "gene_name",
    where = "ModerateImpact = 1",
    output = varsetfile,
    verbose = FALSE
  )
  
  ## save keep-list
  keeplist <- withr::local_tempfile()
  readr::write_lines(paste0("ALS", 1:100),
                     file = keeplist)
  
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--summariseGeno",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--output=sumgeno.txt.gz",
        "--cohort=pheno",
        sprintf("--varSet=%s", varsetfile),
        sprintf("--VAR_id=%s", varids),
        "--pheno=pheno",
        "--memlimit=5000",
        "--geneticModel=allelic,dominant",
        "--checkPloidy=GRCh38",
        sprintf("--keep=%s", keeplist),
        "--splitBy=sex",
        "--minCallrateVar=0.9",
        "--maxCallrateVar=1",
        "--minCallrateSM=0.9",
        "--maxCallrateSM=1",
        "--minMAF=0.00001",
        "--maxMAF=0.05",
        "--minMAC=1",
        "--maxMAC=1000",
        "--minCarriers=1",
        "--maxCarriers=1000",
        "--minCarrierFreq=0.00001",
        "--maxCarrierFreq=0.05",
        "--not-strict",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    cohort = "pheno",
    VAR_id = as.character(1:10),
    pheno = "pheno",
    memlimit = 5000,
    geneticModel = c("allelic", "dominant"),
    checkPloidy = "GRCh38",
    keep = paste0("ALS", 1:100),
    output = "sumgeno.txt.gz",
    splitBy = "sex",
    minCallrateVar = 0.9,
    maxCallrateVar = 1,
    minCallrateSM = 0.9,
    maxCallrateSM = 1,
    minMAF = 0.00001,
    maxMAF = 0.05,
    minMAC = 1,
    maxMAC = 1000,
    minCarriers = 1,
    maxCarriers = 1000,
    minCarrierFreq = 0.00001,
    maxCarrierFreq = 0.05,
    strict = FALSE,
    verbose = FALSE
  )
  expect_s4_class(mock_args[["object"]], "gdb")
  expect_s4_class(mock_args[["varSet"]], "varSetFile")
  expected_args <- c(expected_args, formals(summariseGeno)[!names(formals(summariseGeno)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("object", "varSet")], 
               mock_args[!names(mock_args) %in% c("object", "varSet")])
  
})


# aggregate
test_that("--aggregate works", {
  # generate mocked function/method that returns the environment
  forms <- rvat:::.retrieve_formals_dispatch("aggregate", "gdb")
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    aggregate = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  varids <- withr::local_tempfile()
  readr::write_lines(1:10, varids)
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--aggregate",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        sprintf("--VAR_id=%s", varids),
        "--output=aggregate.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    x = gdb(rvatData::rvat_example("rvatData.gdb")),
    VAR_id = as.character(1:10),
    output = "aggregate.txt.gz"
  )
  expected_args <- c(expected_args, formals(aggregate)[!names(formals(aggregate)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_s4_class(mock_args[["x"]], "gdb")
  expect_equal(expected_args[!names(expected_args) %in% c("x")], 
               mock_args[!names(mock_args) %in% c("x")])
  
  
  # test non-defaults 
  ## save varids
  varids <- withr::local_tempfile()
  readr::write_lines(1:10, varids)
  
  ## save varsetfile
  varsetfile <- withr::local_tempfile()
  buildVarSet(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    varSetName = "Moderate",
    unitTable = "varInfo",
    unitName = "gene_name",
    where = "ModerateImpact = 1",
    output = varsetfile,
    verbose = FALSE
  )
  
  ## save keep-list
  keeplist <- withr::local_tempfile()
  readr::write_lines(paste0("ALS", 1:100),
                     file = keeplist)
  
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--aggregate",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--cohort=pheno",
        sprintf("--varSet=%s", varsetfile),
        sprintf("--VAR_id=%s", varids),
        "--pheno=pheno",
        "--memlimit=5000",
        "--geneticModel=allelic,dominant", 
        "--imputeMethod=missingToRef",
        "--MAFweights=mb",
        "--checkPloidy=GRCh38",
        sprintf("--keep=%s", keeplist),
        "--output=aggregate.txt.gz",
        "--signif=3",
        "--minCallrateVar=0.9",
        "--maxCallrateVar=1",
        "--minCallrateSM=0.9",
        "--maxCallrateSM=1",
        "--minMAF=0.00001",
        "--maxMAF=0.05",
        "--minMAC=1",
        "--maxMAC=1000",
        "--minCarriers=1",
        "--maxCarriers=1000",
        "--minCarrierFreq=0.00001",
        "--maxCarrierFreq=0.05",
        "--not-strict",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    x = gdb(rvatData::rvat_example("rvatData.gdb")),
    cohort = "pheno",
    VAR_id = as.character(1:10),
    pheno = "pheno",
    memlimit = 5000,
    geneticModel = c("allelic", "dominant"),
    imputeMethod = "missingToRef",
    MAFweights = "mb",
    checkPloidy = "GRCh38",
    keep = paste0("ALS", 1:100),
    output = "aggregate.txt.gz",
    signif = 3,
    minCallrateVar = 0.9,
    maxCallrateVar = 1,
    minCallrateSM = 0.9,
    maxCallrateSM = 1,
    minMAF = 0.00001,
    maxMAF = 0.05,
    minMAC = 1,
    maxMAC = 1000,
    minCarriers = 1,
    maxCarriers = 1000,
    minCarrierFreq = 0.00001,
    maxCarrierFreq = 0.05,
    verbose = FALSE,
    strict = FALSE
  )
  expect_s4_class(mock_args[["x"]], "gdb")
  expect_s4_class(mock_args[["varSet"]], "varSetFile")
  expected_args <- c(expected_args, formals(aggregate)[!names(formals(aggregate)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("x", "varSet")], 
               mock_args[!names(mock_args) %in% c("x", "varSet")])
  
})

# mergeAggregateFiles
test_that("--mergeAggregateFiles works", {
  # generate mocked function/method that returns the environment
  forms <- formals(mergeAggregateFiles)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    mergeAggregateFiles = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  aggregatefile1 <- withr::local_tempfile()
  aggregatefile2 <- withr::local_tempfile()
  gdb <- gdb(rvatData::rvat_example("rvatData.gdb"))
  agg1 <- aggregate(gdb, VAR_id = 1:10, cohort = "pheno", output = aggregatefile1, verbose = FALSE)
  agg2 <- aggregate(gdb, VAR_id = 11:20, cohort = "pheno", output = aggregatefile2, verbose = FALSE)
  filelist <- withr::local_tempfile()
  readr::write_lines(c(aggregatefile1, aggregatefile2), 
                     file = filelist
                     )
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--mergeAggregateFiles",
        sprintf("--filelist=%s", filelist),
        "--not-checkDups",
        "--output=mergeAggregateFiles.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    output = "mergeAggregateFiles.txt.gz"
  )
  expected_args <- c(expected_args, formals(mergeAggregateFiles)[!names(formals(mergeAggregateFiles)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_s4_class(mock_args[["object"]], "aggregateFileList")
  expect_equal(expected_args[!names(expected_args) %in% c("object")], 
               mock_args[!names(mock_args) %in% c("object")])
  
  
  # test non-defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--mergeAggregateFiles",
        sprintf("--filelist=%s", filelist),
        "--not-checkDups",
        "--output=mergeAggregateFiles.txt.gz",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    output = "mergeAggregateFiles.txt.gz",
    verbose = FALSE
  )
  expected_args <- c(expected_args, formals(mergeAggregateFiles)[!names(formals(mergeAggregateFiles)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_s4_class(mock_args[["object"]], "aggregateFileList")
  expect_equal(expected_args[!names(expected_args) %in% c("object")], 
               mock_args[!names(mock_args) %in% c("object")])
  
})

# collapseAggregateFiles
test_that("--collapseAggregateFiles works", {
  # generate mocked function/method that returns the environment
  forms <- formals(collapseAggregateFiles)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    collapseAggregateFiles = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  aggregatefile1 <- withr::local_tempfile()
  aggregatefile2 <- withr::local_tempfile()
  gdb <- gdb(rvatData::rvat_example("rvatData.gdb"))
  agg1 <- aggregate(gdb, VAR_id = 1:10, cohort = "pheno", output = aggregatefile1, verbose = FALSE)
  agg2 <- aggregate(gdb, VAR_id = 11:20, cohort = "pheno", output = aggregatefile2, verbose = FALSE)
  filelist <- withr::local_tempfile()
  readr::write_lines(c(aggregatefile1, aggregatefile2), 
                     file = filelist
  )
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--collapseAggregateFiles",
        sprintf("--filelist=%s", filelist),
        "--not-checkDups",
        "--output=collapseAggregateFiles.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    output = "collapseAggregateFiles.txt.gz"
  )
  expected_args <- c(expected_args, formals(collapseAggregateFiles)[!names(formals(collapseAggregateFiles)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_s4_class(mock_args[["object"]], "aggregateFileList")
  expect_equal(expected_args[!names(expected_args) %in% c("object")], 
               mock_args[!names(mock_args) %in% c("object")])
  
  
  # test non-defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--collapseAggregateFiles",
        sprintf("--filelist=%s", filelist),
        "--not-checkDups",
        "--output=collapseAggregateFiles.txt.gz",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    output = "collapseAggregateFiles.txt.gz",
    verbose = FALSE
  )
  expected_args <- c(expected_args, formals(collapseAggregateFiles)[!names(formals(collapseAggregateFiles)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_s4_class(mock_args[["object"]], "aggregateFileList")
  expect_equal(expected_args[!names(expected_args) %in% c("object")], 
               mock_args[!names(mock_args) %in% c("object")])
  
})


# assocTest-gdb
test_that("--assocTest works", {
  # generate mocked function/method that returns the environment
  forms <- rvat:::.retrieve_formals_dispatch("assocTest", "gdb")
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    assocTest = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--assocTest",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--pheno=pheno",
        "--test=firth",
        "--output=assocTest.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    pheno = "pheno",
    test = "firth",
    output = "assocTest.txt.gz"
  )
  expected_args <- c(expected_args, formals(assocTest)[!names(formals(assocTest)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_s4_class(mock_args[["object"]], "gdb")
  expect_equal(expected_args[names(mock_args)][!names(expected_args) %in% c("object")], 
               mock_args[!names(mock_args) %in% c("object")])
  
  
  # test non-defaults 
  ## save varids
  varids <- withr::local_tempfile()
  readr::write_lines(1:10, varids)
  
  ## save varsetfile
  varsetfile <- withr::local_tempfile()
  null <- buildVarSet(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    varSetName = "Moderate",
    unitTable = "varInfo",
    unitName = "gene_name",
    where = "ModerateImpact = 1",
    output = varsetfile,
    verbose = FALSE
  )
  
  ## save resamplingfile
  resamplingfile <- withr::local_tempfile()
  null <- buildResamplingFile(nSamples = 25000, nResampling = 100, output = resamplingfile)
  
  ## save keep-list
  keeplist <- withr::local_tempfile()
  readr::write_lines(paste0("ALS", 1:100),
                     file = keeplist)
  
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--assocTest",
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--pheno=pheno,sex",
        "--test=firth,glm,scoreSPA",
        "--cohort=pheno",
        sprintf("--varSet=%s", varsetfile),
        sprintf("--VAR_id=%s", varids),
        "--name=test",
        "--continuous",
        "--singlevar",
        "--covar=PC1/PC1,PC2,PC3,PC4",
        "--offset=offset",
        "--geneticModel=allelic,dominant", 
        "--imputeMethod=missingToRef",
        "--MAFweights=mb",
        "--maxitFirth=5000",
        "--checkPloidy=GRCh38",
        sprintf("--keep=%s", keeplist),
        "--output=assocTest.txt.gz",
        "--methodResampling=permutation",
        sprintf("--resamplingFile=%s",resamplingfile),
        "--nResampling=1000",
        "--outputResampling=outputresampling.txt.gz",
        "--memlimitResampling=5000",
        "--minCallrateVar=0.9",
        "--maxCallrateVar=1",
        "--minCallrateSM=0.9",
        "--maxCallrateSM=1",
        "--minMAF=0.00001",
        "--maxMAF=0.05",
        "--minMAC=1",
        "--maxMAC=1000",
        "--minCarriers=1",
        "--maxCarriers=1000",
        "--minCarrierFreq=0.00001",
        "--maxCarrierFreq=0.05",
        "--memlimit=5000",
        "--not-strict",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    pheno = c("pheno", "sex"),
    test = c("firth", "glm", "scoreSPA"),
    cohort = "pheno",
    VAR_id = as.character(1:10),
    name = "test",
    continuous = TRUE,
    singlevar = TRUE,
    covar = list("PC1", paste0("PC", 1:4)),
    offset = "offset",
    geneticModel = c("allelic", "dominant"),
    imputeMethod = "missingToRef",
    MAFweights = "mb",
    maxitFirth = 5000,
    checkPloidy = "GRCh38",
    keep = paste0("ALS", 1:100),
    output = "assocTest.txt.gz",
    methodResampling = "permutation",
    nResampling = 1000,
    outputResampling = "outputresampling.txt.gz",
    memlimitResampling = 5000,
    minCallrateVar = 0.9,
    maxCallrateVar = 1,
    minCallrateSM = 0.9,
    maxCallrateSM = 1,
    minMAF = 0.00001,
    maxMAF = 0.05,
    minMAC = 1,
    maxMAC = 1000,
    minCarriers = 1,
    maxCarriers = 1000,
    minCarrierFreq = 0.00001,
    maxCarrierFreq = 0.05,
    memlimit = 5000,
    verbose = FALSE,
    strict = FALSE
  )
  expect_s4_class(mock_args[["object"]], "gdb")
  expect_s4_class(mock_args[["varSet"]], "varSetFile")
  expect_s4_class(mock_args[["resamplingFile"]], "resamplingFile")
  expected_args <- c(expected_args, formals(assocTest)[!names(formals(assocTest)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("object", "varSet", "resamplingFile")], 
               mock_args[!names(mock_args) %in% c("object", "varSet", "resamplingFile")])
  
})



# assocTest-aggregateFile
test_that("--assocTest works", {
  # generate mocked function/method that returns the environment
  forms <- rvat:::.retrieve_formals_dispatch("assocTest", "aggregateFile")
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    assocTest = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  varsetfile <- withr::local_tempfile()
  aggregatefile <- withr::local_tempfile()
  genesetfile <- withr::local_tempfile()
  
  null <- buildVarSet(
    object = gdb(rvatData::rvat_example("rvatData.gdb")),
    varSetName = "Moderate",
    unitTable = "varInfo",
    unitName = "gene_name",
    where = "ModerateImpact = 1",
    output = varsetfile,
    verbose = FALSE
  )
  
  null <- aggregate(
    x = gdb(rvatData::rvat_example("rvatData.gdb")),
    varSet = getVarSet(varSetFile(varsetfile), unit = c("SOD1", "ABCA4")),
    output = aggregatefile,
    verbose = FALSE
  )
 
  null <- buildGeneSet(
    data = list(
      "A" = c("SOD1"),
      "B" = c("SOD1", "ABCA4")
    ),
    output = genesetfile,
    verbose = FALSE
  )
  
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--assocTest",
        "--pheno=pheno",
        "--test=firth",
        sprintf("--aggregateFile=%s",aggregatefile),
        sprintf("--geneSet=%s",genesetfile),
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--output=assocTest.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    pheno = "pheno",
    test = "firth",
    output = "assocTest.txt.gz"
  )
  expected_args <- c(expected_args, formals(assocTest)[!names(formals(assocTest)) %in% names(expected_args)])
  expect_s4_class(mock_args[["object"]], "aggregateFile")
  expect_s4_class(mock_args[["gdb"]], "gdb")
  expect_s4_class(mock_args[["geneSet"]], "geneSetFile")
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("object", "gdb", "geneSet")], 
               mock_args[!names(mock_args) %in% c("object", "gdb", "geneSet")])
  
  
  # test non-defaults 
  
  ## save keep-list
  keeplist <- withr::local_tempfile()
  readr::write_lines(paste0("ALS", 1:100),
                     file = keeplist)
  
  ## save keep-list
  dropunits <- withr::local_tempfile()
  readr::write_lines("SOD1",
                     file = dropunits)
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--assocTest",
        sprintf("--aggregateFile=%s",aggregatefile),
        "--pheno=pheno",
        "--test=firth,glm",
        sprintf("--geneSet=%s",genesetfile),
        sprintf("--gdb=%s", rvatData::rvat_example("rvatData.gdb")),
        "--cohort=pheno",
        "--name=test",
        "--continuous",
        "--covar=PC1/PC1,PC2,PC3,PC4,varcount",
        "--substractCovar=varcount",
        sprintf("--dropUnits=%s", dropunits),
        "--maxitFirth=5000",
        sprintf("--keep=%s", keeplist),
        "--output=assocTest.txt.gz",
        "--quiet",
        "--not-strict"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    pheno = "pheno",
    test = c("firth", "glm"),
    cohort = "pheno",
    name = "test",
    continuous = TRUE,
    covar = list(c("PC1"), c(paste0("PC", 1:4), "varcount")),
    substractCovar = "varcount",
    dropUnits = "SOD1",
    maxitFirth = 5000,
    keep = paste0("ALS", 1:100),
    output = "assocTest.txt.gz",
    verbose = FALSE,
    strict = FALSE
  )
  
  expect_s4_class(mock_args[["object"]], "aggregateFile")
  expect_s4_class(mock_args[["gdb"]], "gdb")
  expect_s4_class(mock_args[["geneSet"]], "geneSetFile")
  expected_args <- c(expected_args, formals(assocTest)[!names(formals(assocTest)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("object", "gdb", "geneSet")], 
               mock_args[!names(mock_args) %in% c("object", "gdb", "geneSet")])
  
})

# buildResamplingFile
test_that("--buildResamplingFile works", {
  # generate mocked function/method that returns the environment
  forms <- formals(buildResamplingFile)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    buildResamplingFile = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--buildResamplingFile",
        "--nSamples=10000",
        "--output=buildResamplingFile.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    nSamples = 10000,
    output = "buildResamplingFile.txt.gz"
  )
  expected_args <- c(expected_args, formals(buildResamplingFile)[!names(formals(buildResamplingFile)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args,
               mock_args)
  
  # test non-defaults 
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--buildResamplingFile",
        "--nSamples=10000",
        "--nResampling=100000",
        "--memlimit=500",
        "--methodResampling=permutation",
        "--output=buildResamplingFile.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    nSamples = 10000,
    nResampling = 100000,
    memlimit = 500,
    methodResampling = "permutation",
    output = "buildResamplingFile.txt.gz"
  )
  
  expected_args <- c(expected_args, formals(buildResamplingFile)[!names(formals(buildResamplingFile)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args, 
               mock_args)
  
})

# buildGeneSet
test_that("--buildGeneSet works", {
  # generate mocked function/method that returns the environment
  forms <- formals(buildGeneSet)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    buildGeneSet = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--buildGeneSet",
        "--gmtpath=data.gmt",
        "--output=buildGeneSet.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    gmtpath = "data.gmt",
    output = "buildGeneSet.txt.gz"
  )
  expected_args <- c(expected_args, formals(buildGeneSet)[!names(formals(buildGeneSet)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args,
               mock_args)
  
  # test non-defaults 
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--buildGeneSet",
        "--gmtpath=data.gmt",
        "--output=buildGeneSet.txt.gz",
        "--sep=|",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    gmtpath = "data.gmt",
    output = "buildGeneSet.txt.gz",
    sep = "|",
    verbose = FALSE
  )
  
  expected_args <- c(expected_args, formals(buildGeneSet)[!names(formals(buildGeneSet)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args, 
               mock_args)
  
})

# buildCorMatrix
test_that("--buildCorMatrix works", {
  # generate mocked function/method that returns the environment
  forms <- formals(buildCorMatrix)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    buildCorMatrix = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  aggregatefile <- withr::local_tempfile()
  gdb <- gdb(rvatData::rvat_example("rvatData.gdb"))
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
  varsetlist <- getVarSet(varsetfile, unit = listUnits(varsetfile)[1:4], varSetName = "Moderate")
  agg <- aggregate(gdb, 
                   varSet = varsetlist,
                   cohort = "pheno", 
                   output = aggregatefile, 
                   verbose = FALSE)
  
  results <- withr::local_tempfile()
  data(rvbresults)
  writeResult(rvbresults[rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",], 
              file = results)
  output <- withr::local_tempfile()
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--buildCorMatrix",
        sprintf("--rvbResult=%s", results),
        sprintf("--aggregateFile=%s", aggregatefile),
        sprintf("--output=%s", output)
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list()
  expect_s4_class(mock_args[["object"]], "rvbResult")
  expect_s4_class(mock_args[["aggregateFile"]], "aggregateFile")
  
  expected_args <- c(expected_args, formals(buildCorMatrix)[!names(formals(buildCorMatrix)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("object", "aggregateFile")], 
               mock_args[!names(mock_args) %in% c("object", "aggregateFile")])
  
  # test non-defaults 
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--buildCorMatrix",
        sprintf("--rvbResult=%s", results),
        sprintf("--aggregateFile=%s", aggregatefile),
        "--memlimit=10000",
        "--minR2=1e-3",
        "--not-makePD",
        "--not-absolute",
        "--maxDist=1000000",
        "--quiet",
        sprintf("--output=%s", output)
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()

  # run cli
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))

  # check against expected args
  expected_args <- list(
    memlimit = 10000,
    minR2 = 1e-3,
    makePD = FALSE,
    maxDist = 1000000,
    absolute = FALSE,
    verbose = FALSE
  )

  expected_args <- c(expected_args, formals(buildCorMatrix)[!names(formals(buildCorMatrix)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("object", "aggregateFile")], 
               mock_args[!names(mock_args) %in% c("object", "aggregateFile")])
  
})


# geneSetAssoc
test_that("--geneSetAssoc works", {
  # generate mocked function/method that returns the environment
  forms <- formals(geneSetAssoc)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    geneSetAssoc = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  
  ## geneSetFile
  genesetfile <- withr::local_tempfile()
  null <- buildGeneSet(
    gmtpath = "../data/c5.go.mf.v2023.2.Hs.symbols.gmt",
    output = genesetfile,
    verbose = FALSE
  )
  
  results <- withr::local_tempfile()
  data(rvbresults)
  writeResult(rvbresults[rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",], file = results)
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--geneSetAssoc",
        sprintf("--rvbResult=%s", results),
        sprintf("--geneSet=%s", genesetfile),
        "--test=lm",
        "--output=geneSetAssoc.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    test = "lm",
    output = "geneSetAssoc.txt.gz"
  )
  expect_s4_class(mock_args[["object"]], "rvbResult")
  expect_s4_class(mock_args[["geneSet"]], "geneSetFile")
  expected_args <- c(expected_args, formals(geneSetAssoc)[!names(formals(geneSetAssoc)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("object", "geneSet")], 
               mock_args[!names(mock_args) %in% c("object", "geneSet")])
  
  # test non-defaults 
  
  ## scorematrix
  scorematrix <- withr::local_tempfile()
  data <- readr::read_tsv("../data/GSE67835_Human_Cortex.txt.gz", 
                          show_col_types = FALSE, 
                          progress = FALSE)
  genes <- data$GENE
  data$GENE <- NULL
  data <- as.matrix(data)
  rownames(data) <- genes
  saveRDS(data, file = scorematrix)

  ## cormatrix
  aggregatefile <- withr::local_tempfile()
  cormatrix <-  withr::local_tempfile()
  gdb <- gdb(rvatData::rvat_example("rvatData.gdb"))
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
  varsetlist <- getVarSet(varsetfile, unit = listUnits(varsetfile)[1:4], varSetName = "Moderate")
  agg <- aggregate(gdb, 
                   varSet = varsetlist,
                   cohort = "pheno", 
                   output = aggregatefile, 
                   verbose = FALSE)
  
  
  mat <- suppressMessages(buildCorMatrix(
    rvbresults[rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",],
    aggregateFile = aggregateFile(aggregatefile),
    verbose = FALSE
  ))
  saveRDS(mat, file = cormatrix)
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--geneSetAssoc",
        sprintf("--rvbResult=%s", results),
        sprintf("--geneSet=%s", genesetfile),
        sprintf("--scoreMatrix=%s", scorematrix),
        sprintf("--cormatrix=%s", cormatrix),
        sprintf("--condition=%s", genesetfile),
        sprintf("--covar=%s", "nvar,nCarriers"),
        "--test=lm",
        "--threshold=1e-5",
        "--Zcutoffs=-3,3",
        "--INT",
        "--scoreCutoffs=-3,3",
        "--minSetSize=5",
        "--maxSetSize=1000",
        "--twoSided",
        "--memlimit=5000",
        "--ID=gene",
        "--output=geneSetAssoc.txt.gz",
        "--quiet"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))

  # check against expected args
  expected_args <- list(
    covar = c("nvar", "nCarriers"),
    test = "lm",
    threshold = 1e-5,
    Zcutoffs = c(-3, 3),
    INT = TRUE,
    scoreCutoffs = c(-3, 3),
    minSetSize = 5,
    maxSetSize = 1000,
    oneSided = FALSE,
    memlimit = 5000,
    ID = "gene",
    output = "geneSetAssoc.txt.gz",
    verbose = FALSE
  )

  expect_s4_class(mock_args[["object"]], "rvbResult")
  expect_s4_class(mock_args[["geneSet"]], "geneSetFile")
  expect_s4_class(mock_args[["cormatrix"]], "Matrix")
  expect_s4_class(mock_args[["condition"]], "geneSetFile")
  expect_true(is(mock_args[["scoreMatrix"]], "matrix"))
  
  expected_args <- c(expected_args, formals(geneSetAssoc)[!names(formals(geneSetAssoc)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args[!names(expected_args) %in% c("object", "geneSet", "cormatrix", "condition", "scoreMatrix")],
               mock_args[!names(mock_args) %in% c("object", "geneSet", "cormatrix", "condition", "scoreMatrix")])

})


# vcfInfo2Table
test_that("--vcfInfo2Table works", {
  # generate mocked function/method that returns the environment
  forms <- formals(vcfInfo2Table)
  mock_function <- function() {}
  formals(mock_function) <- forms
  body(mock_function) <- quote({
    mock_args <<- as.list(environment())
    invisible()
  })
  local_mocked_bindings(
    vcfInfo2Table = mock_function,
    .package = "rvat"
  )
  
  # test defaults 
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--vcfInfo2Table",
        sprintf("--vcf=%s", rvatData::rvat_example("rvatData.vcf.gz")),
        "--output=vcfInfo2Table.txt.gz"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    vcf = rvatData::rvat_example("rvatData.vcf.gz"),
    output = "vcfInfo2Table.txt.gz"
  )
  expected_args <- c(expected_args, formals(vcfInfo2Table)[!names(formals(vcfInfo2Table)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args,
               mock_args)
  
  # test non-defaults 
  
  local_mocked_bindings(
    commandArgs = function(trailingOnly = TRUE) {
      c("--vcfInfo2Table",
        sprintf("--vcf=%s", rvatData::rvat_example("rvatData.vcf.gz")),
        "--output=vcfInfo2Table.txt.gz",
        "--not-splitMultiallelic"
      )
    },
    .package = "base"
  )
  collected_args <- rvat:::collect_args()
  
  # run cli 
  suppressMessages(rvat_cli(
    args = collected_args[["args"]],
    args_raw = collected_args[["args_raw"]]
  ))
  
  # check against expected args
  expected_args <- list(
    vcf = rvatData::rvat_example("rvatData.vcf.gz"),
    output = "vcfInfo2Table.txt.gz",
    splitMultiallelic = FALSE
  )
  
  expected_args <- c(expected_args, formals(vcfInfo2Table)[!names(formals(vcfInfo2Table)) %in% names(expected_args)])
  expected_args <- expected_args[sort(names(expected_args))]
  mock_args <- mock_args[sort(names(mock_args))]
  expect_equal(expected_args, 
               mock_args)
  
})

