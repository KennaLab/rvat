library(rvatData)

# connect to varsetfile from rvatData package
# to build a varsetfile, see ?buildVarSet
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))

# retrieve specific genes, this will return a varSetList (see ?varSetList)
varsets <- getVarSet(varsetfile, unit = c("SOD1", "FUS"))
head(varsets)

# the varsetfile contains multiple records per gene, (High impact, CADD scores etc.)
unique(listVarSets(varsetfile))

# specific varSets can be selected using the `varSetName` parameter
varsets <- getVarSet(varsetfile, unit = c("SOD1", "FUS"), varSetName = "CADD")
head(varsets)

# see ?varSetFile and ?varSetList for more details on connecting and handling varsetfiles.
# see e.g., ?assocTest and ?aggregate for downstream methods that can loop through varsetfiles and varsetlists.
