rvat_cli <- function(args, args_raw) {
  
  # check number of  
  if (sum(args_raw %in% rvat_cli_methods_flags) > 1) {
    message(sprintf("Rare variant analysis toolkit (rvat)

Error: More than one valid function call provided: %s

Provide a valid function call for detailed help.\n%s", 
                    paste(rvat_cli_methods_flags[rvat_cli_methods_flags %in% args_raw], collapse = ","),
                    rvat_cli_valid_calls))
    stopQuietly()
  }
  
  if (sum(args_raw %in% rvat_cli_methods_flags) == 0) {
    message(sprintf("Rare variant analysis toolkit (rvat)

Error: Please provide a valid function call.

Provide a valid function call for detailed help.\n%s", 
                    rvat_cli_valid_calls))
    
    stopQuietly()
  }
  
  # return help message if --help argument is specified or a function is called
  # without arguments
  check_help(args = args_raw)
  
  # check verbosity
  if ( !args[["quiet"]] || "verbose" %in% args_raw ) verbose <- TRUE else verbose <- FALSE
  
  # check strictness
  if ( !args[["not-strict"]] || "strict" %in% args_raw ) strict <- TRUE else strict <- FALSE
  
  # print options
  if (verbose) message(sprintf("Analysis started at: %s\n\n", as.character(round(Sys.time(), units = "secs"))))
  if (verbose) print_args(args, args_raw, rvat_cli_flags)
  
  # buildGdb
  if ( "buildGdb" %in% args_raw )
  {
    required <- c("output", "vcf")
    expected <- c("buildGdb", "quiet", names(formals(buildGdb)))
    check_args(func_name = "buildGdb", args = args_raw, help = help, required = required, expected = expected)
    
    buildGdb(vcf = args[["vcf"]],
             output = args[["output"]],
             skipIndexes = args[["skipIndexes"]],
             skipVarRanges = args[["skipVarRanges"]],
             overWrite = args[["overWrite"]],
             genomeBuild = args[["genomeBuild"]],
             memlimit = args[["memlimit"]],
             verbose = verbose)
  }
  
  # concatGdb
  if ( "concatGdb" %in% args_raw )
  {
    required <- c("targets", "output")
    expected <- c("concatGdb", "quiet", names(formals(concatGdb)))
    check_args(func_name = "concatGdb", args = args_raw, help = help, required = required, expected = expected)
    
    concatGdb(targets = args[["targets"]],
              output = args[["output"]],
              skipRemap = args[["skipRemap"]],
              skipIndexes  = args[["skipIndexes"]],
              verbose = verbose)
  }
  
  # subsetGdb
  if ( "subsetGdb" %in% args_raw )
  {
    required <- c("gdb", "output")
    expected <- c("subsetGdb", "quiet", "gdb", names(formals(subsetGdb)))
    check_args(func_name = "subsetGdb", args = args_raw, help = help, required = required, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    
    if (is.null(args[["VAR_id"]])) {VAR_id <- NULL} else {VAR_id <- readLines(args[["VAR_id"]])}
    if (is.null(args[["tables"]])) {tables <- NULL} else {tables <- unlist(strsplit(args[["tables"]], split = ","))}
    
    subsetGdb(object = gdb,
              output = args[["output"]],
              intersection = args[["intersection"]],
              where = args[["where"]],
              VAR_id = VAR_id,
              tables = tables,
              skipIndexes = args[["skipIndexes"]],
              overWrite = args[["overWrite"]],
              verbose = verbose
    )
  }
  
  # uploadAnno
  if ("uploadAnno" %in% args_raw ) {
    
    required <- c("gdb", "name", "value")
    expected <- c("uploadAnno", "gdb", "quiet", names(formals(uploadAnno)))
    check_args(func_name = "uploadAnno", args = args_raw, help = help, required = required, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    uploadAnno(object = gdb,
               name = args[["name"]],
               value = args[["value"]],
               sep = args[["sep"]],
               skipIndexes = args[["skipIndexes"]],
               skipRemap = args[["skipRemap"]],
               ignoreAlleles = args[["ignoreAlleles"]],
               keepUnmapped = args[["keepUnmapped"]],
               mapRef = args[["mapRef"]],
               verbose = verbose
               )
  }
  
  # uploadCohort
  if ( "uploadCohort" %in% args_raw )
  {
    required <- c("gdb", "name", "value")
    expected <- c("uploadCohort", "gdb", "quiet", names(formals(uploadCohort)))
    check_args(func_name = "uploadCohort", args = args_raw, help = help, required = required, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    uploadCohort(gdb,
                 name = args[["name"]],
                 value = args[["value"]],
                 sep = args[["sep"]],
                 verbose = verbose
                 )
  }
  
  # listAnno
  if ( "listAnno" %in% args_raw )
  {
    required <- c("gdb")
    expected <- c("listAnno", "gdb",  "quiet", names(formals(listAnno)))
    check_args(func_name = "listAnno", args = args_raw, help = help, required = required, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    listAnno(gdb)
  }
  
  # listCohort
  if ( "listCohort" %in% args_raw )
  {
    required <- c("gdb")
    expected <- c("listCohort", "gdb", "quiet", names(formals(listCohort)))
    check_args(func_name = "listCohort", args = args_raw, help = help, required = required, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    listCohort(gdb)
  }
  
  # dropTable
  if ( "dropTable" %in% args_raw )
  {
    required <- c("gdb", "name")
    expected <- c("dropTable", "gdb", "quiet", names(formals(dropTable)))
    check_args(func_name = "dropTable", args = args_raw, help = help, required = required, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    dropTable(
      object = gdb,
      name = args[["name"]],
      verbose = verbose
    )
  }
  
  # mapVariants
  if ( "mapVariants" %in% args_raw )
  {
    required <- c("gdb")
    expected <- c("mapVariants", "gdb",  "quiet", names(formals(mapVariants)))    
    check_args(func_name = "mapVariants", args = args_raw, help = help, required = required, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    if (!is.null(args[["fields"]])) fields <- unlist(strsplit(args[["fields"]], split = ",")) else {fields <- NULL}
    if (length(args[["bedCols"]]) > 0) {bedCols <- unlist(strsplit(args[["bedCols"]], split = ","))} else {bedCols <- args[["bedCols"]]}
    
    mapVariants(gdb,
                ranges = args[["ranges"]],
                gff = args[["gff"]],
                bed = args[["bed"]],
                bedCols = bedCols,
                fields = fields,
                uploadName = args[["uploadName"]],
                output = args[["output"]],
                sep = args[["sep"]],
                overWrite = args[["overWrite"]],
                skipIndexes = args[["skipIndexes"]],
                verbose = verbose
    )
  }
  
  # buildVarSet
  if ( "buildVarSet" %in% args_raw )
  {
    
    required <- c("gdb", "varSetName", "unitTable", "unitName", "output")
    expected <- c("buildVarSet", "gdb", "quiet", names(.retrieve_formals_dispatch("buildVarSet", "gdb")))    
    check_args(func_name = "buildVarSet", args = args_raw, help = help, required = required, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    buildVarSet(
      object = gdb,
      varSetName = args[["varSetName"]],
      unitTable = args[["unitTable"]],
      unitName = args[["unitName"]],
      output = args[["output"]],
      intersection = args[["intersection"]],
      where = args[["where"]],
      weightName = args[["weightName"]],
      verbose = verbose
      )
  }
  
  # spatialClust
  if ( "spatialClust" %in% args_raw )
  {
    required <- c("gdb", "varSetName", "unitTable", "unitName", "output", "windowSize", "overlap")
    expected <- c("spatialClust", "gdb", "quiet", names(formals(spatialClust)))
    check_args(func_name = "spatialClust", args = args_raw, help = help, required = required, expected = expected)

    gdb <- gdb(args[["gdb"]])
    windowSize <- as.numeric(unlist(strsplit(args[["windowSize"]], split = ",")))
    overlap <- as.numeric(unlist(strsplit(args[["overlap"]], split = ",")))

    spatialClust(
      object = gdb,
      output = args[["output"]],
      varSetName = args[["varSetName"]],
      unitTable = args[["unitTable"]],
      unitName = args[["unitName"]],
      windowSize = windowSize,
      overlap = overlap,
      intersection = args[["intersection"]],
      where = args[["where"]],
      weightName = args[["weightName"]],
      posField = args[["posField"]],
      minTry = args[["minTry"]]
    )
  }

  # summariseGeno
  if ( "summariseGeno" %in% args_raw )
  {
    required <- c("gdb", "output")
    expected <- c("summariseGeno", "gdb", "quiet", "not-strict", names(.retrieve_formals_dispatch("summariseGeno", "gdb")))
    check_args(func_name = "summariseGeno", args = args_raw, help = help, required = required, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    if (!is.null(args[["varSet"]])) varSet <- varSetFile(args[["varSet"]]) else varSet <- NULL
    if (!is.null(args[["VAR_id"]])) VAR_id <- readLines(args[["VAR_id"]]) else VAR_id <- NULL
    if (!is.null(args[["keep"]])) keep <- readLines(args[["keep"]]) else keep <- NULL
    
    summariseGeno(
      object = gdb, 
      cohort = args[["cohort"]], 
      varSet = varSet,
      VAR_id = VAR_id,
      pheno = args[["pheno"]],
      memlimit = args[["memlimit"]],
      geneticModel = unlist(strsplit(args[["geneticModel"]], split = ",")),
      checkPloidy = args[["checkPloidy"]],
      keep = keep,
      output = args[["output"]],
      splitBy = args[["splitBy"]],
      minCallrateVar = args[["minCallrateVar"]],
      maxCallrateVar = args[["maxCallrateVar"]],
      minCallrateSM = args[["minCallrateSM"]],
      maxCallrateSM = args[["maxCallrateSM"]],
      minMAF = args[["minMAF"]],
      maxMAF = args[["maxMAF"]],
      minMAC = args[["minMAC"]],
      maxMAC = args[["maxMAC"]],
      minCarriers = args[["minCarriers"]],
      maxCarriers = args[["maxCarriers"]],
      minCarrierFreq = args[["minCarrierFreq"]],
      maxCarrierFreq = args[["maxCarrierFreq"]],
      strict = strict,
      verbose = verbose
    )
  }
  
  # aggregate 
  if ( "aggregate" %in% args_raw )
  {
    required <- c("gdb", "output")
    required_one_of <- c("varSet", "VAR_id")
    expected <- c("aggregate", "gdb", "quiet", "not-strict", names(.retrieve_formals_dispatch("aggregate", "gdb")))
    check_args(func_name = "aggregate", args = args_raw, help = help, required = required, required_one_of = required_one_of, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    if (!is.null(args[["varSet"]])) varSet <- varSetFile(args[["varSet"]]) else varSet <- NULL
    if (!is.null(args[["VAR_id"]])) VAR_id <- readLines(args[["VAR_id"]]) else VAR_id <- NULL
    if (!is.null(args[["keep"]])) keep <- readLines(args[["keep"]]) else keep <- NULL
    if (!is.null(args[["imputeMethod"]])) imputeMethod <- args[["imputeMethod"]] else imputeMethod <- "meanImpute"
    
    aggregate(
      x = gdb, 
      cohort = args[["cohort"]],
      varSet = varSet,
      VAR_id = VAR_id,
      pheno = args[["pheno"]],
      memlimit = args[["memlimit"]],
      geneticModel = unlist(strsplit(args[["geneticModel"]], split = ",")),
      imputeMethod = imputeMethod,
      MAFweights = unlist(strsplit(args[["MAFweights"]],split=",")),
      checkPloidy = args[["checkPloidy"]],
      keep = keep,
      output = args[["output"]],
      signif = args[["signif"]],
      minCallrateVar = args[["minCallrateVar"]],
      maxCallrateVar = args[["maxCallrateVar"]],
      minCallrateSM = args[["minCallrateSM"]],
      maxCallrateSM = args[["maxCallrateSM"]],
      minMAF = args[["minMAF"]],
      maxMAF = args[["maxMAF"]],
      minMAC = args[["minMAC"]],
      maxMAC = args[["maxMAC"]],
      minCarriers = args[["minCarriers"]],
      maxCarriers = args[["maxCarriers"]],
      minCarrierFreq = args[["minCarrierFreq"]],
      maxCarrierFreq = args[["maxCarrierFreq"]],
      verbose = verbose,
      strict = strict
    )
  }
  
  # mergeAggregateFiles 
  if ( "mergeAggregateFiles" %in% args_raw )
  {
    required <- c("filelist", "output")
    expected <- c("mergeAggregateFiles", "quiet", "not-checkDups", names(formals(mergeAggregateFiles)), names(formals(aggregateFileList)))
    check_args(func_name = "mergeAggregateFiles", args = args_raw, help = help, required = required, expected = expected)

    files <- readLines(args[["filelist"]])
    if ( !args[["not-checkDups"]] || "checkDups" %in% args_raw ) checkDups <- TRUE else checkDups <- FALSE
    aggregateFileList <- aggregateFileList(files,
                                           checkDups = checkDups)
    

    mergeAggregateFiles(
      object = aggregateFileList,
      output = args[["output"]],
      verbose = verbose
    )
  }
  
  # collapseAggregateFiles 
  if ( "collapseAggregateFiles" %in% args_raw )
  {
    required <- c("filelist", "output")
    expected <- c("collapseAggregateFiles", "quiet", "not-checkDups", names(formals(collapseAggregateFiles)), names(formals(aggregateFileList)))
    check_args(func_name = "collapseAggregateFiles", args = args_raw, help = help, required = required, expected = expected)
    
    files <- readLines(args[["filelist"]])
    if ( !args[["not-checkDups"]] || "checkDups" %in% args_raw ) checkDups <- TRUE else checkDups <- FALSE
    aggregateFileList <- aggregateFileList(files,
                                           checkDups = checkDups)
    
    collapseAggregateFiles(
      object = aggregateFileList,
      output = args[["output"]],
      verbose = verbose
    )
  }

  # assocTest-gdb
  if ( "assocTest" %in% args_raw && !"aggregateFile" %in% args_raw )
  {
    required <- c("gdb", "pheno", "test")
    required_one_of <- c("output", "outputResampling")
    expected <- c("assocTest", "gdb", "quiet", "not-strict", "seed", names(.retrieve_formals_dispatch("assocTest", "gdb")))
    check_args(func_name = "assocTest", args = args_raw, help = help, required = required, required_one_of = required_one_of, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    if (!is.null(args[["varSet"]])) varSet <- varSetFile(args[["varSet"]]) else varSet <- NULL
    if (!is.null(args[["VAR_id"]])) VAR_id <- readLines(args[["VAR_id"]]) else VAR_id <- NULL
    if (!is.null(args[["name"]])) name <- args[["name"]] else name <- "none"
    if (!is.null(args[["covar"]])) {
      covar <- lapply(unlist(strsplit(args[["covar"]], split = "/")), 
                      FUN = function(x) unlist(strsplit(x, split = ",")))
    } else {
      covar <- NULL
    }
    if (!is.null(args[["keep"]])) keep <- readLines(args[["keep"]]) else keep <- NULL
    if (!is.null(args[["resamplingFile"]])) resamplingFile <- resamplingFile(args[["resamplingFile"]]) else resamplingFile <- NULL
    if (!is.null(args[["seed"]])) set.seed(args[["seed"]])
    if (args[["outputResampling"]] == "FALSE") outputResampling <- FALSE else outputResampling <- args[["outputResampling"]]
    
    assocTest(
      object = gdb, 
      pheno = unlist(strsplit(args[["pheno"]], split = ",")), 
      test = unlist(strsplit(args[["test"]], split = ",")), 
      cohort = args[["cohort"]],
      varSet = varSet, 
      VAR_id = VAR_id,
      name = name,
      continuous = args[["continuous"]],
      singlevar = args[["singlevar"]],
      covar = covar,
      offset = args[["offset"]],
      geneticModel = unlist(strsplit(args[["geneticModel"]], split = ",")),
      imputeMethod = args[["imputeMethod"]],
      MAFweights = unlist(strsplit(args[["MAFweights"]],split = ",")),
      maxitFirth = args[["maxitFirth"]],
      checkPloidy = args[["checkPloidy"]],
      keep = keep,
      output = args[["output"]],
      methodResampling = args[["methodResampling"]],
      resamplingFile = resamplingFile,
      nResampling = args[["nResampling"]],
      outputResampling = outputResampling,
      memlimitResampling = args[["memlimitResampling"]],
      minCallrateVar = args[["minCallrateVar"]],
      maxCallrateVar = args[["maxCallrateVar"]],
      minCallrateSM = args[["minCallrateSM"]],
      maxCallrateSM = args[["maxCallrateSM"]],
      minMAF = args[["minMAF"]],
      maxMAF = args[["maxMAF"]],
      minMAC = args[["minMAC"]],
      maxMAC = args[["maxMAC"]],
      minCarriers = args[["minCarriers"]],
      maxCarriers = args[["maxCarriers"]],
      minCarrierFreq = args[["minCarrierFreq"]],
      maxCarrierFreq = args[["maxCarrierFreq"]],
      memlimit = args[["memlimit"]],
      verbose = verbose,
      strict = strict
    )
  }
  
  
  # assocTest-aggregateFile
  if ( "assocTest" %in% args_raw && "aggregateFile" %in% args_raw )
    {
    
    required <- c("aggregateFile", "pheno", "test", "geneSet", "gdb", "output")
    expected <- c("assocTest", "aggregateFile", "quiet", "not-strict", names(.retrieve_formals_dispatch("assocTest", "aggregateFile")))
    check_args(func_name = "assocTest", args = args_raw, help = help, required = required, expected = expected)
    
    gdb <- gdb(args[["gdb"]])
    genesetfile <- geneSetFile(args[["geneSet"]])
    aggregatefile <- aggregateFile(args[["aggregateFile"]])
    if (!is.null(args[["covar"]])) {
      covar <- lapply(unlist(strsplit(args[["covar"]], split = "/")), 
                      FUN = function(x) unlist(strsplit(x, split = ",")))
    } else {
      covar <- NULL
    }
    if (!is.null(args[["name"]])) name <- args[["name"]] else name <- "none"
    if (!is.null(args[["keep"]])) keep <- readLines(args[["keep"]]) else keep <- NULL
    if (!is.null(args[["dropUnits"]])) dropUnits <- readLines(args[["dropUnits"]]) else dropUnits <- NULL
    
    assocTest(
      object = aggregatefile,
      pheno = unlist(strsplit(args[["pheno"]], split = ",")), 
      test = unlist(strsplit(args[["test"]], split = ",")), 
      geneSet = genesetfile,
      gdb = gdb, 
      cohort = args[["cohort"]], 
      name = name,
      continuous = args[["continuous"]],
      covar = covar,
      substractCovar = args[["substractCovar"]],
      dropUnits = dropUnits,
      maxitFirth = args[["maxitFirth"]],
      keep = keep,
      output = args[["output"]],
      verbose = verbose,
      strict = strict
    )
  }
  
  # buildResamplingFile
  if ( "buildResamplingFile" %in% args_raw )
  {
    
    required <- c("nSamples", "output")
    expected <- c("buildResamplingFile", names(formals(buildResamplingFile)))
    check_args(func_name = "buildResamplingFile", args = args_raw, help = help, required = required, expected = expected)
    
    if (!is.null(args[["seed"]])) set.seed(args[["seed"]])
    if (!is.null(args[["methodResampling"]])) methodResampling <- args[["methodResampling"]] else  methodResampling <- "permutation"
    
    buildResamplingFile(
      nSamples = args[["nSamples"]],
      nResampling = args[["nResampling"]],
      memlimit = args[["memlimit"]],
      methodResampling = methodResampling,
      output = args[["output"]]
    )
  }
  
  # buildGeneSet ----------------------------------------------------------
  if ( "buildGeneSet" %in% args_raw ) {
    required <- c("buildGeneSet", "gmtpath", "output")
    expected <- c("buildGeneSet", "quiet", names(formals(buildGeneSet)))
    check_args(func_name = "buildGeneSet", args = args_raw, help = help, required = required, expected = expected)
    
    buildGeneSet(
      gmtpath = args[["gmtpath"]],
      output = args[["output"]],
      sep = args[["sep"]],
      verbose = verbose
    )
  }
  
  # buildCorMatrix ----------------------------------------------------------
  if ( "buildCorMatrix" %in% args_raw ) {
    required <- c("buildCorMatrix", "rvbResult", "aggregateFile", "output")
    expected <- c("buildCorMatrix", "rvbResult", "output", "quiet", "not-makePD", "not-absolute", names(formals(buildCorMatrix)))
    check_args(func_name = "buildCorMatrix", args = args_raw, help = help, required = required, expected = expected)
    
    result <- rvbResult(args[["rvbResult"]])
    aggregatefile <- aggregateFile(args[["aggregateFile"]])
    if ( !args[["not-makePD"]] || "makePD" %in% args_raw ) makePD <- TRUE else makePD <- FALSE
    if ( !args[["not-absolute"]] || "absolute" %in% args_raw ) absolute <- TRUE else absolute <- FALSE
    
    cormatrix <- buildCorMatrix(
      object = result,
      aggregateFile = aggregatefile,
      memlimit = args[["memlimit"]],
      minR2 = args[["minR2"]],
      makePD = makePD,
      absolute = absolute,
      maxDist = args[["maxDist"]],
      verbose = verbose
    )
    saveRDS(cormatrix, file = args[["output"]])
  }
  
  # geneSetAssoc --------------------------------------------------------
  if ( "geneSetAssoc" %in% args_raw )
  {
    required <- c("geneSet", "rvbResult", "test", "output")
    expected <- c("geneSetAssoc", "geneSet", "rvbResult", "twoSided", "quiet", names(formals(geneSetAssoc)))
    check_args(func_name = "geneSetAssoc", args = args_raw, help = help, required = required, expected = expected)
    
    # load results
    result <- rvbResult(args[["rvbResult"]])
    
    # connect to genesetfile
    genesetfile <- geneSetFile(args[["geneSet"]])
    
    # load scoreMatrix, if specified
    if (!is.null(args[["scoreMatrix"]])) scoreMatrix <- readRDS(args[["scoreMatrix"]]) else scoreMatrix <- NULL
    
    # load cormatrix, if specified
    if (!is.null(args[["cormatrix"]])) cormatrix <- readRDS(args[["cormatrix"]]) else cormatrix <- NULL
    
    # condition
    if (!is.null(args[["condition"]])) {
      if (args[["condition_type"]] == "geneSet") {
        condition <- geneSetFile(args[["condition"]])
      } else if (args[["condition_type"]] == "matrix") {
        condition <- readRDS(args[["condition"]])
      } else if (args[["condition_type"]] == "vector") {
        condition <- unlist(strsplit(args[["condition"]], split = ","))
      } else {
        stop("--condition_type should be one of 'geneSet', 'matrix', or 'vector'")
      }
    } else {
      condition <- NULL
    }
    
    if (!is.null(args[["covar"]])) covar <- unlist(strsplit(args[["covar"]], split = ",")) else covar <- NULL
    if (!is.null(args[["Zcutoffs"]])) Zcutoffs <- as.numeric(unlist(strsplit(args[["Zcutoffs"]], split = ","))) else Zcutoffs <- NULL
    if (!is.null(args[["scoreCutoffs"]])) scoreCutoffs <- as.numeric(unlist(strsplit(args[["scoreCutoffs"]], split = ","))) else scoreCutoffs <- NULL
    if ( !args[["twoSided"]] || "oneSided" %in% args_raw ) oneSided <- TRUE else oneSided <- FALSE
    
    geneSetAssoc(
      object = result,
      geneSet = genesetfile,
      scoreMatrix = scoreMatrix,
      cormatrix = cormatrix,
      condition = condition,
      covar = covar,
      test = unlist(strsplit(args[["test"]], split = ",")), 
      threshold = args[["threshold"]],
      Zcutoffs = Zcutoffs,
      INT = args[["INT"]],
      scoreCutoffs = scoreCutoffs,
      minSetSize = args[["minSetSize"]],
      maxSetSize = args[["maxSetSize"]],
      oneSided = oneSided,
      memlimit = args[["memlimit"]],
      ID = args[["ID"]],
      output = args[["output"]],
      verbose = verbose
    )
  }
  
  # vcfInfo2Table
  if ( "vcfInfo2Table" %in% args_raw )
  {
    required <- c("vcfInfo2Table", "vcf", "output")
    expected <- c("vcfInfo2Table", "not-splitMultiallelic", names(formals(vcfInfo2Table)))
    check_args(func_name = "vcfInfo2Table", args = args_raw, help = help, required = required, expected = expected)
    
    if ( !args[["not-splitMultiallelic"]] || "splitMultiAllelic" %in% args_raw ) splitMultiallelic <- TRUE else splitMultiallelic <- FALSE
    
    vcfInfo2Table(
      vcf = args[["vcf"]],
      output = args[["output"]],
      splitMultiallelic = splitMultiallelic
    )
  } 

  if(length(warnings()) > 0) {
    message("Overview of warnings:")
    summary(warnings())
  }
  
  if (verbose) {
    message("\nFinished!")
    message(sprintf("End time: %s", as.character(round(Sys.time(), units = "secs"))))
  }
}