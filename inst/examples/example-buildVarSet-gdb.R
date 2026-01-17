library(rvatData)

# Build a varset including variants with a moderate predicted impact
gdb <- create_example_gdb()
varsetfile_moderate <- tempfile()
buildVarSet(object = gdb,
            output = varsetfile_moderate,
            varSetName = "Moderate",
            unitTable = "varInfo",
            unitName = "gene_name",
            where = "ModerateImpact = 1")

# Build a varset that contains CADD scores
varsetfile_cadd <- tempfile()
buildVarSet(object = gdb,
            output = varsetfile_cadd,
            varSetName = "CADD",
            unitTable = "varInfo",
            unitName = "gene_name",
            weightName = "CADDphred")

# connect to varsetfile and retrieve variant sets
varsetfile <- varSetFile(varsetfile_moderate)
varsets <- getVarSet(varsetfile, unit = c("SOD1", "FUS"))

# see ?getVarSet, ?varSetFile and ?varSetList for more details on connecting and handling varsetfiles.
# see e.g., ?assocTest and ?aggregate for downstream methods that can loop through varsetfiles and varsetlists.
