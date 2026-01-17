library(rvatData)
gdb <- create_example_gdb()

# assocTest-aggdb allows for running association tests on pre-constructed aggregates.
# below we first generate the aggregates based on a varSetFile
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
varset <- getVarSet(varsetfile, unit = c("NEK1", "SOD1", "ABCA4"), varSetName = "High")
aggfile <- tempfile()
aggregate(
  x = gdb,
  varSet = varset,
  maxMAF = 0.001,
  output = aggfile,
  verbose = FALSE
)

# connect to aggdb, see ?aggdb for more details
aggdb <- aggdb(aggfile)

# build an example genesetlist, see ?buildGeneSet for details
genesetlist <- buildGeneSet(
  list(
    "geneset1" = c("SOD1", "NEK1"),
    "geneset2" = c("ABCA4", "SOD1", "NEK1")
  )
)

# perform association tests using assocTest
# note that this is very similar to running association tests on a genoMatrix (?`assocTest-genoMatrix`)
# or on a gdb (?`assocTest-gdb`). The main difference being that pre-constructed aggregates are used from
# the aggdb, and requires a genesetFile/geneSetList (?geneSetFile) to be provided to the `geneSet` argument.
aggAssoc <- assocTest(
  aggdb,
  gdb = gdb,
  test = c("glm", "firth"),
  cohort = "pheno",
  pheno = "pheno",
  geneSet = genesetlist,
  covar = paste0("PC", 1:4),
  verbose = FALSE
)
