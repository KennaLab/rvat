library(rvatData)
gdb <- create_example_gdb()

# generate for variant summaries for list of variants
sumgeno <- tempfile()
summariseGeno(gdb, cohort = "pheno", VAR_id = 1:100, output = sumgeno)

# generate for variant summaries for varSetFile
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
varsets <- getVarSet(varsetfile, unit = c("SOD1", "FUS"), varSetName = "High")
summariseGeno(gdb, cohort = "pheno", varSet = varsets, output = sumgeno)

# variant summaries can be generated for subgroups using the `splitBy` parameter.
# this will result in an additional column in the output for the subgroups
summariseGeno(
  gdb,
  cohort = "pheno",
  VAR_id = 1:100,
  splitBy = "pheno",
  output = sumgeno
)
data <- read.table(sumgeno, header = TRUE)
# contains 'pheno' column
head(data)

# summariseGeno can be ran directly on a genoMatrix
data(GT)
sumgeno <- summariseGeno(GT)
