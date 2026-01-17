library(rvatData)

# connect to varsetfile from rvatData package
# to build a varsetfile, see ?buildVarSet
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))

# extract a couple of genes from the varSetFile, which whill return a varSetList
varsetlist <- getVarSet(varsetfile, c("SOD1", "FUS"))
varsetlist

# many of the varSetFile methods are also available for a varSetList
# for example, getVarSet can also be used to extract specific genes or varSets from a varSetList
getVarSet(varsetlist, unit = "SOD1")

# list included units and varSets 
units <- listUnits(varsetlist)
head(units)
varsets <- listVarSets(varsetlist)
head(varsets)

# several basic list operations work on a varSetList
length(varsetlist)
varsetlist[1:2]
varsetlist[[1]]

# extract metadata
metadata(varsetlist)
getRvatVersion(varsetlist)
getGdbId(varsetlist)

# all varsets in in a varsetlist can be collapsed into one varSet using the
# collapseVarSetList method
collapseVarSetList(varsetlist)

# varSetLists can be written to a varSetFile on disk using the write method
output <- tempfile() 
write(varsetlist, output)
varsetfile <- varSetFile(output)
varsetfile

# see e.g., ?assocTest and ?aggregate for downstream methods that can loop through varSets included in a varSetlist
