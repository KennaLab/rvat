library(rvatData)

# build a geneSetList
# can also be build based on a GMT-file (see ?buildGeneSet)
genesetlist <- buildGeneSet(
  list(
    "geneset1" = c("SOD1", "NEK1"),
    "geneset2" = c("ABCA4", "SOD1", "NEK1"),
    "geneset3" = c("FUS", "NEK1")
  )
)

# extract a couple of gene sets from the geneSetList, which will return a new geneSetList
getGeneSet(genesetlist, c("geneset1", "geneset2"))

# list included gene sets and units
genesets <- listGeneSets(genesetlist)
head(genesets)
units <- listUnits(genesetlist)
head(units)

# several basic list operations work on a geneSetList
length(genesetlist)
genesetlist[1:2]
genesetlist[[1]]

# write a geneset list to a geneSetFile on disk (see ?geneSetFile)
file <- tempfile()
write(genesetlist, file)
genesetfile <- geneSetFile(file)

# exclude units from all genesets included in a geneSetList
dropUnits(genesetlist, unit = "NEK1")

# remap IDs
linker <- data.frame(
  gene_name = c("SOD1", "NEK1", "FUS", "ABCA4"),
  gene_id = c(
    "ENSG00000142168",
    "ENSG00000137601",
    "ENSG00000089280",
    "ENSG00000198691"
  )
)
genesetlist_remapped <- remapIDs(genesetlist, linker)
listUnits(genesetlist_remapped)

# see ?geneSetAssoc and ?`assocTest-aggdb` to run gene set analyses
