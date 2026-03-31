library(rvatData)
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
varset <- getVarSet(varsetfile, unit = c("NEK1"), varSetName = "High")[[1]]

# list variants and weights included in the varSet
listVars(varset)
listWeights(varset)

# note that usually you'll work with varSets in either a varSetList or
# a varSetFile (see ?varSetList and ?varSetFile)
