library(rvatData)

# build a geneSetList
# can also be build based on a GMT-file (see ?buildGeneSet)
file <- tempfile()
genesetfile <- buildGeneSet(
  list(
    "geneset1" = c("SOD1", "NEK1"),
    "geneset2" = c("ABCA4", "SOD1", "NEK1"),
    "geneset3" = c("FUS", "NEK1")
  ),
  output = file
)

# connect to a geneSetFile
genesetfile <- geneSetFile(file)

# extract a couple of gene sets from the geneSetFile, which will return a geneSetList (see ?geneSetList)
getGeneSet(genesetfile, c("geneset1", "geneset2"))

# list included gene sets
genesets <- listGeneSets(genesetfile)
head(genesets)

# see ?geneSetAssoc and ?`assocTest-aggdb` to run gene set analyses
