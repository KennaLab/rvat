# build a genesetlist from a list (see ?geneSetList)
genesetlist <- buildGeneSet(
  list("geneset1" = c("SOD1", "NEK1"),
       "geneset2" = c("ABCA4", "SOD1", "NEK1"),
       "geneset3" = c("FUS", "NEK1")
       ))

# specify the output parameter to write to disk in the geneSetFile format (see ?geneSetFile)
file <- tempfile()
buildGeneSet(
  list("geneset1" = c("SOD1", "NEK1"),
       "geneset2" = c("ABCA4", "SOD1", "NEK1"),
       "geneset3" = c("FUS", "NEK1")
  ),
  output = file
  )
genesetfile <- geneSetFile(file)

# the `gmtpath` parameter can be used to build a geneset from a mSigDb GMT-file
# see the tutorials on the RVAT website for examples
