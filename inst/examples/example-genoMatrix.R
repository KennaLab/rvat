library(rvatData)
data(GT)

# Basic operations -------------------------------------------------------------

# retrieve rowData (i.e. variant info)
rowData(GT)

# retrieve colData (i.e. sample info)
colData(GT)

# retrieve sample IDs
samples <- colnames(GT)
head(samples)

# retrieve VAR ids
vars <- rownames(GT)
head(vars)

# check dimensions
nrow(GT)
ncol(GT)
dim(GT)

# A genoMatrix object can be subsetted similarly as a data.frame:
GT[1:5, 1:5]

# Subset samples based on sample info
GT[, GT$pheno == 1]

# Subset first two variants
GT[1:2, ]

# Extract variant/sample summaries --------------------------------------------------

# calculate allele frequencies, allele counts etc.
af <- getAF(GT)
maf <- getMAF(GT)
ac <- getAC(GT)
mac <- getMAC(GT)
carriers <- getNCarriers(GT)

# generate call-rates
var_cr <- getCR(GT) # variant call-rates
sample_cr <- getCR(GT, var = FALSE) # sample call-rates

# generate variant summaries
varsummary <- summariseGeno(GT)

# Recode genotypes  --------------------------------------------------

# flip variants with AF > 0.5 to the minor allele
GT <- flipToMinor(GT)

# recode genotypes to domiant/recessive models
recode(GT, geneticModel = "dominant")
recode(GT, geneticModel = "recessive")
# see ?recode for details

# recode missing genotypes
recode(GT, imputeMethod = "meanImpute")
recode(GT, imputeMethod = "missingToRef")
# see ?recode for details

# generate aggregate (burden) scores

# by default, scores will be added to the colData
aggregate(recode(GT, imputeMethod = "meanImpute"))

# set `returnGT = FALSE` to return aggregates as a vector instead
aggregate <- aggregate(
  recode(GT, imputeMethod = "meanImpute"),
  returnGT = FALSE
)
head(aggregate)

# Update cohort and variant info in genoMatrix --------------------------
gdb <- create_example_gdb()
anno <- getAnno(gdb, "varInfo", fields = c("VAR_id", "CADDphred", "PolyPhen"))
updateGT(GT, anno = anno)
pheno <- colData(GT)
updateGT(GT, SM = colData(GT)[, 1:4])

# Get variant carriers

# retrieve a data.frame that lists the samples carrying each respective
# variant in the genoMatrix. Additional variant and sample info can be included
# using the `rowDataFields` and `colDataFields` respectively.
carriers <- getCarriers(
  GT,
  rowDataFields = c("REF", "ALT"),
  colDataFields = c("superPop")
)
head(carriers)

# Perform rare variant tests:

# burden test (firth)
rvb <- assocTest(
  GT,
  pheno = "pheno",
  covar = c("sex", "PC1", "PC2", "PC3", "PC4"),
  test = "firth"
)

# single variant tests
sv <- assocTest(
  GT,
  pheno = "pheno",
  covar = c("sex", "PC1", "PC2", "PC3", "PC4"),
  test = "scoreSPA",
  singlevar = TRUE
)
# see ?assocTest for details
