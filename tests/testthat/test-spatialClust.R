gdb <- create_example_gdb()

# # check mapToCds
test_that("spatialClust works",{

  gtf <- rtracklayer::import("../data/Homo_sapiens.GRCh38.105.rvatData.gtf")
  transcripts <- unique(gtf$transcript_id)
  output <- withr::local_tempfile()
  expect_no_error({suppressMessages(mapToCDS(
    gdb,
    gff = gtf,
    transcript_id = transcripts,
    exonPadding = 12,
    output = output
  ))})
  output2 <- withr::local_tempfile()
  R.utils::gunzip(output, destname = output2, remove = TRUE, overwrite = TRUE)
  # upload
  uploadAnno(gdb, value = output2, name = "cdsPOS", skipRemap = TRUE, verbose = FALSE)
}
)

# check snapshot
test_that("spatialClust snapshot is identical",{
  output <- withr::local_tempfile()
  suppressWarnings(spatialClust(
    gdb,
    unitTable = "varInfo",
    varSetName = "moderate",
    unitName = "transcript_id",
    where = "ModerateImpact = 1",
    windowSize = 30,
    overlap = 15,
    intersection = "cdsPOS",
    posField = "cdsPOS",
    output = output
  ))

  varsetfile <- varSetFile(output)
  assoc <- assocTest(gdb,
                     varSet = getVarSet(varsetfile, unit = listUnits(varsetfile)[1:10]),
                     cohort="pheno",
                     pheno = "pheno",
                     covar = c("sex", "PC1", "PC2", "PC3", "PC4"),
                     test = "glm",
                     name = "spatial_clustering",
                     verbose = FALSE)
  metadata(assoc)$creationDate <- NA_character_
  metadata(assoc)$gdbPath <- NA_character_
  expect_snapshot_value(assoc, style = "serialize")
    }
  )