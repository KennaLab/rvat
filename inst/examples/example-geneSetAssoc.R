library(rvatData)
data(rvbresults)
res <- rvbresults[
  rvbresults$test == "firth" &
    rvbresults$varSetName == "ModerateImpact",
]

# example genesetlist used in examples below (see ?buildGeneSet on build geneSetLists/geneSetFiles)
genesetlist <- buildGeneSet(
  list(
    "geneset1" = c("SOD1", "NEK1"),
    "geneset2" = c("ABCA4", "SOD1", "NEK1"),
    "geneset3" = c("FUS", "NEK1")
  )
)

# Perform competitive gene set analysis using a linear model
GSAresults <- geneSetAssoc(
  res,
  genesetlist,
  covar = c("nvar"),
  test = c("lm")
)

# Outlying gene association scores can be remedied by either setting Z-score cutoffs (i.e. all Z-scores exceeding these values will be set to the respective cutoff),
# or inverse normal transforming the Z-scores:
GSAresults <- geneSetAssoc(
  res,
  genesetlist,
  covar = c("nvar"),
  test = c("lm"),
  maxSetSize = 500,
  Zcutoffs = c(-4, 4) # lower and upper bounds
)

GSAresults <- geneSetAssoc(
  res,
  genesetlist,
  covar = c("nvar"),
  test = c("lm"),
  maxSetSize = 500,
  INT = TRUE # perform inverse normal transformation
)


# Conditional gene set analyses can be performed to test whether gene sets are associated independently with the phenotype of interest.
# In the example below we test whether gene sets are independent of geneset1
GSAresults <- geneSetAssoc(
  res,
  condition = getGeneSet(genesetlist, "geneset1"),
  geneSet = genesetlist,
  covar = c("nvar"),
  test = c("lm"),
  maxSetSize = 500
)

# perform two-sided tests by setting `oneSided = FALSE`
GSAresults <- geneSetAssoc(
  res,
  genesetlist,
  covar = c("nvar"),
  test = c("lm"),
  maxSetSize = 500,
  oneSided = TRUE
)

# Test whether the proportion of P-values below a specified threshold is greater than the proportion outside of it.
# Using Fisher's exact test
# The `threshold` parameter specifies the P-value cutoff to define significant genes:
GSAresults <- geneSetAssoc(
  res,
  genesetlist,
  test = c("fisher"),
  threshold = 1e-4,
  maxSetSize = 500
)
