library(rvatData)
data(rvbresults)

# rvatResult inherits from the DataFrame classes
# standard methods also work on an rvatResult
head(rvbresults)
nrow(rvbresults)
ncol(rvbresults)
dim(rvbresults)
rvbresults_df <- as.data.frame(rvbresults)
head(rvbresults_df)
rvbresults[1:3,]
rvbresults[,1:10]
rvbresults[["unit"]][1:5]

# write results
file <- tempfile()
writeResult(rvbresults, file = file)
rvbresults <- rvbResult(file)

# similarly, for single variant results
data(GTsmall)
sv <- assocTest(
  GTsmall,
  covar = c("sex", paste0("PC", 1:4)),
  pheno = "pheno",
  test = "scoreSPA",
  singlevar = TRUE,
  verbose = FALSE
)
svresultfile <- tempfile()
writeResult(sv, file = svresultfile)
sv <- singlevarResult(sv)
head(sv)

# merge 
merge <- merge(rvbresults[1:23], as.data.frame(rvbresults[,c(1, 23:28)]), by = "unit")

# show summary
summary(rvbresults)

# show top results
topResult(rvbresults, n = 10)

# qqplot
man <- qqplot(rvbresults[rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",],
                 label = "unit")

# generate a manhattan plot
man <- manhattan(rvbresults[rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",],
                 label = "unit", 
                 contigs = "GRCh38")


# see ?ACAT and ?geneSetAssoc for downstream analyses on rvatResults
