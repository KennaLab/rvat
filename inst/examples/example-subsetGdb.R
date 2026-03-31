library(rvatData)
gdb <- create_example_gdb()

# Make a gdb subset that includes only variants annotated to SOD1
output <- tempfile()
subsetGdb(
  gdb,
  intersection = "varInfo",
  where = "gene_name = 'SOD1'",
  output = output
)
gdb_subset <- gdb(output)

# Specific tables can be selected to include.
# all other user-uploaded annotation and cohort tables will be excluded
subsetGdb(
  gdb,
  intersection = "varInfo",
  where = "gene_name = 'SOD1'",
  tables = "varInfo",
  output = output,
  overWrite = TRUE
)
gdb_subset <- gdb(output)

# subset gdbs based on list of VAR ids
anno <- getAnno(
  gdb,
  "var",
  range = data.frame(CHROM = "chr16", start = 31191399, end = 31191605)
)

subsetGdb(
  gdb,
  VAR_id = anno$VAR_id,
  output = output,
  overWrite = TRUE
)
gdb_subset <- gdb(output)
