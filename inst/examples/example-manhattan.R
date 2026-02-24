library(rvatData)
data(rvbresults)

# generate manhatan plot
man <- manhattan(
  rvbresults[
    rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",
  ],
  label = "unit",
  contigs = "GRCh38"
)

# if many overlapping gene label, try setting `labelRepel = TRUE`
man <- manhattan(
  rvbresults[
    rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",
  ],
  label = "unit",
  labelRepel = TRUE,
  contigs = "GRCh38"
)

# alter the significane threshold using the `threshold` parameter
man <- manhattan(
  rvbresults[
    rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",
  ],
  label = "unit",
  labelRepel = TRUE,
  threshold = 1e-5,
  contigs = "GRCh38"
)

# the threshold for displaying labels can be set using the `labelTrheshold` parameter
man <- manhattan(
  rvbresults[
    rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",
  ],
  label = "unit",
  labelRepel = TRUE,
  labelThreshold = 1e-12,
  contigs = "GRCh38"
)
