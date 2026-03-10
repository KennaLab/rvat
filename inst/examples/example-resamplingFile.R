library(rvatData)
# build and connect to a resamplingFile
file <- tempfile(fileext = ".gz")
buildResamplingFile(nSamples = 25000, nResampling = 100, output = file)
resamplingfile <- resamplingFile(file)

# perform resampled association tests
gdb <- create_example_gdb()
assoc <- assocTest(
  gdb,
  VAR_id = 1:10,
  cohort = "pheno",
  pheno = "pheno",
  covar = paste0("PC", 1:4),
  test = c("skat"),
  resamplingFile = resamplingfile,
  verbose = FALSE
)
