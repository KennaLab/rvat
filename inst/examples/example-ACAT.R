library(rvatData)

# this will first ACAT P-values across statistical tests
# and then ACAT these P-values across varSets
data(rvbresults)
rvbresults <- rvbresults[1:1000, ]
ACAT(
  rvbresults,
  aggregate = list("test", "varSetName")
)

# alternatively, by providing a vector P-values across aggregates will be combined in one stage
ACAT(
  rvbresults,
  aggregate = c("test", "varSetName")
)

# Input P-value shouldn't be exactly 0 or 1 (see: https://github.com/yaowuliu/ACAT).
# By default, P-values that are exactly 0 or 1 are reset (`fixpval = TRUE`) to the
# minimum P-value (>0) and maximum P-value (<1) in the results.
# Alternatives include:
# manual minmax values:
ACAT(
  rvbresults,
  aggregate = list("test", "varSetName"),
  fixpval = TRUE,
  fixpval_method = "manual",
  fixpval_maxP = 0.9999,
  fixpval_minP = 1e-32
)

# Liu method (see FAQ on: https://github.com/yaowuliu/ACAT)
ACAT(
  rvbresults,
  aggregate = list("test", "varSetName"),
  fixpval = TRUE,
  fixpval_method = "Liu"
)
