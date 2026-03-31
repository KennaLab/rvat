library(rvatData)

gdb <- create_example_gdb()
anno <- getAnno(gdb, "varinfo", where = "gene_name in ('SOD1', 'FUS')")
varsetfile_from_df <- tempfile()
buildVarSet(
  anno,
  unitName = "gene_name",
  fields = c("HighImpact"),
  output = varsetfile_from_df
)

# connect to varsetfile and retrieve variant sets
varsetfile <- varSetFile(varsetfile_from_df)
varsets <- getVarSet(varsetfile, unit = c("SOD1", "FUS"))

# see ?getVarSet, ?varSetFile and ?varSetList for more details on connecting and handling varsetfiles.
# see e.g., ?assocTest and ?aggregate for downstream methods that can loop through varsetfiles and varsetlists.
