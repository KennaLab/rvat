library(rvatData)
gdb <- create_example_gdb()

# generate two aggregate files
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
aggdb1 <- tempfile()
aggregate(x = gdb,
          varSet = getVarSet(varsetfile, unit = c("SOD1", "FUS"), varSetName = "High"),
          maxMAF = 0.001,
          output = aggdb1,
          verbose = FALSE)

aggdb2 <- tempfile()
aggregate(x = gdb,
          varSet = getVarSet(varsetfile, unit = c("NEK1"), varSetName = "High"),
          maxMAF = 0.001,
          output = aggdb2,
          verbose = FALSE)

# merge using mergeAggDbs
aggdb <- tempfile()
agglist <- aggdbList(c(aggdb1, aggdb2))
mergeAggDbs(
  agglist,
  output = aggdb
  )
