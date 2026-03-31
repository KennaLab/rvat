library(rvatData)
data(rvbresults)
gdb <- create_example_gdb()

# generate the aggregates based on a varSetFile
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
varset <- getVarSet(
  varsetfile,
  unit = c("NEK1", "SOD1", "ABCA4"),
  varSetName = "High"
)
aggfile <- tempfile()
aggregate(
  x = gdb,
  varSet = varset,
  maxMAF = 0.001,
  output = aggfile,
  verbose = FALSE
)

# build a block-wise correlation matrix
cormatrix <- buildCorMatrix(
  rvbresults,
  aggdb = aggdb(aggfile)
)
