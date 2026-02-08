library(rvatData)
aggfile <- tempfile()
gdb <- create_example_gdb()

# generate aggregates for varSets
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
varsets <- getVarSet(varsetfile, unit = c("SOD1", "FUS"), varSetName = "High")
aggregate(x = gdb,
          varSet = varsets,
          maxMAF = 0.001,
          output = aggfile,
          verbose = FALSE)

# generate for aggregates for list of variants
aggregate(x = gdb,
          VAR_id = 1:100,
          maxMAF = 0.001,
          output = aggfile,
          verbose = FALSE)

# use recessive model
aggregate(x = gdb,
          varSet = varsets,
          maxMAF = 0.001,
          geneticModel = "recessive",
          output = aggfile,
          verbose = FALSE)

# apply MAF weighting 
aggregate(x = gdb,
          varSet = varsets,
          maxMAF = 0.001,
          MAFweights = "mb",
          output = aggfile,
          verbose = FALSE)