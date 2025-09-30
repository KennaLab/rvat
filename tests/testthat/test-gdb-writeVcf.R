# write vcf from gdb and build gdb from that vcf
test_that("writeVcf roundtrip works", {
  gdb1 <- create_example_gdb()

  ## write gdb to vcf
  output_writevcf <- withr::local_tempfile()
  suppressMessages(writeVcf(
    gdb1,
    output_writevcf
  ))

  ## build a gdb again from vcf written from gdb
  output_writevcf_buildgdb <- withr::local_tempfile()
  suppressMessages(buildGdb(
    vcf = output_writevcf,
    output = output_writevcf_buildgdb,
    overWrite = TRUE,
    genomeBuild = "GRCh38"
  ))
  gdb2 <- gdb(output_writevcf_buildgdb)
  SM1 <- getCohort(gdb1, "SM")
  SM2 <- getCohort(gdb2, "SM")
  expect_identical(SM1$IID, SM2$IID)

  # include SM with sex info
  DBI::dbWriteTable(gdb2, "SM", SM1, overwrite = TRUE)

  # compare gdbs
  compare_gdbs(gdb1, gdb2, check_tables = FALSE, check_var = FALSE)

  ## var tables won't be identical because writeVcf drops some INFO fields
  var1 <- getAnno(gdb1, "var")
  var1$INFO <- "."
  var2 <- getAnno(gdb2, "var")
  expect_identical(var1, var2)
})

# same as above but without genotypes (includeGeno = FALSE)
test_that("writeVcf w/o genotypes roundtrip works ", {
  gdb1 <- create_example_gdb()
  output_writevcf <- withr::local_tempfile()
  suppressMessages(writeVcf(
    object = gdb1,
    includeGeno = FALSE,
    output_writevcf
  ))

  var1 <- getAnno(gdb1, "var")
  var2 <- readr::read_tsv(
    output_writevcf,
    show_col_types = FALSE,
    col_names = FALSE,
    comment = "#"
  )
  colnames(var2) <- colnames(var1)[2:9]
  var1$QUAL <- as.numeric(var1$QUAL)
  var1$INFO <- "."
  expect_equal(var1[, 2:9], as.data.frame(var2))
})


test_that("writeVcf with subset of IIDs and VAR_ids works", {
  # include subset of IIDs/VAR_ids
  gdb1 <- create_example_gdb()
  output_writevcf <- withr::local_tempfile()
  output_writevcf_buildgdb <- withr::local_tempfile()
  IIDs <- sample(getCohort(gdb1, "SM")$IID, 100)
  writeVcf(
    object = gdb1,
    IID = IIDs,
    VAR_id = 1:100,
    output_writevcf,
    verbose = FALSE
  )
  suppressMessages(buildGdb(
    vcf = output_writevcf,
    output = output_writevcf_buildgdb,
    overWrite = TRUE,
    genomeBuild = "GRCh38"
  ))
  gdb2 <- gdb(output_writevcf_buildgdb)
  GT1 <- getGT(gdb1, VAR_id = 1:100, cohort = "SM", verbose = FALSE)
  GT2 <- getGT(
    gdb2,
    VAR_id = getAnno(gdb2, table = "var", fields = "VAR_id")$VAR_id,
    cohort = "SM",
    verbose = FALSE
  )

  ### ignore metadata in comparison
  metadata(GT1)$gdb <- NA_character_
  metadata(GT1)$gdbId <- NA_character_
  colData(GT1)$sex <- 0
  metadata(GT2)$gdb <- NA_character_
  metadata(GT2)$gdbId <- NA_character_
  expect_identical(GT1[, colnames(GT1) %in% IIDs], GT2)
})
