library(rvatData)
gdb <- create_example_gdb()

# generate the aggregates based on a varSetFile
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
varset <- getVarSet(varsetfile, unit = c("NEK1", "SOD1", "ABCA4"), varSetName = "High")
aggfile <- tempfile()
aggregate(x = gdb,
          varSet = varset,
          maxMAF = 0.001,
          output = aggfile,
          verbose = FALSE)

# connect to aggdb, see ?aggdb for more details
aggdb <- aggdb(aggfile)

# list units and samples in aggdb 
head(listUnits(aggdb))
head(listSamples(aggdb))

# retrieve aggregates 
aggregates <- getUnit(aggdb, unit = "SOD1")

# see ?`assocTest-aggdb` for details on running association tests on an aggdb
