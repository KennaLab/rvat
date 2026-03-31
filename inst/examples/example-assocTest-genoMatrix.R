library(rvatData)
data(GTsmall)

# run a firth burden test on a binary phenotype
rvb <- assocTest(
  GTsmall,
  pheno = "pheno",
  covar = c("PC1", "PC2", "PC3", "PC4"),
  test = "firth",
  name = "example"
)

# run ACAT-v and SKAT tests
rvb <- assocTest(
  GTsmall,
  pheno = "pheno",
  covar = c("PC1", "PC2", "PC3", "PC4"),
  test = c("acatvfirth", "skat_burden_robust", "skato_robust"),
  name = "example"
)

# run a burden test on a continuous phenotype
rvb <- assocTest(
  GTsmall,
  pheno = "age",
  continuous = TRUE,
  covar = c("PC1", "PC2", "PC3", "PC4"),
  test = c("lm", "skat", "acatv"),
  name = "example"
)

# run single variant tests on a binary phenotype
sv <- assocTest(
  GTsmall,
  pheno = "pheno",
  singlevar = TRUE,
  covar = c("PC1", "PC2", "PC3", "PC4"),
  test = c("firth", "glm", "scoreSPA"),
  name = "example",
  minCarriers = 1
)

# apply variant filters
rvb <- assocTest(
  GTsmall,
  pheno = "pheno",
  covar = c("PC1", "PC2", "PC3", "PC4"),
  test = c("firth", "skat_robust", "acatv"),
  name = "example",
  maxMAF = 0.001,
  minCarriers = 2,
  minCallrateVar = 0.9,
  minCallrateSM = 0.95
)

# Perform MAF-weighted burden tests (madsen-browning)
rvb <- assocTest(
  GTsmall,
  pheno = "pheno",
  covar = c("PC1", "PC2", "PC3", "PC4"),
  test = c("firth", "skat_robust", "acatv"),
  MAFweights = "mb"
)

# Perform weighted burden tests with custom weights
gdb <- create_example_gdb()
# add cadd scores
caddscores <- getAnno(gdb, "varInfo", VAR_id = rownames(GTsmall))
caddscores <- caddscores[match(rownames(GTsmall), caddscores$VAR_id), ]
# CADD-weighted burden test
rvb <- assocTest(
  recode(GTsmall, weights = caddscores$CADDphred),
  pheno = "pheno",
  covar = c("PC1", "PC2", "PC3", "PC4"),
  test = c("firth", "skat_robust", "acatv")
)

# Perform recessive burden test
rvb <- assocTest(
  GTsmall,
  pheno = "pheno",
  covar = c("PC1", "PC2", "PC3", "PC4"),
  test = c("firth", "skat_robust", "acatv"),
  geneticModel = "recessive"
)

# Resampled burden test
rvb <- assocTest(
  GTsmall,
  pheno = "pheno",
  covar = c("PC1", "PC2", "PC3", "PC4"),
  test = c("skat", "skat_burden", "acatv"),
  name = "example",
  methodResampling = "permutation",
  nResampling = 100
)
