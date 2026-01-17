genesetlist <- buildGeneSet(
  list("geneset1" = c("SOD1", "NEK1"),
       "geneset2" = c("ABCA4", "SOD1", "NEK1"),
       "geneset3" = c("FUS", "NEK1")
       ))
geneset <- genesetlist[[1]]

# list units and weights included in the geneSet
listUnits(geneset)
listWeights(geneset)

# note that usually you'll work with geneSets in either a geneSetList or 
# a geneSetFile (see ?geneSetList and ?geneSetFile)
