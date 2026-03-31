test_that("getRanges works", {
  # run getRanges on a varsetfile
  gdb <- create_example_gdb()
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
  ranges <- expect_no_error(getRanges(varsetfile, gdb = gdb))

  # check whether running getRanges on varSetList vs. varSetFile
  # results in identical output
  varsetlist <- getVarSet(varsetfile, unit = listUnits(varsetfile))
  ranges_from_varsetlist <- getRanges(varsetlist, gdb = gdb)
  expect_identical(ranges, ranges_from_varsetlist)

  # compare with manually retrieving ranges
  var <- getAnno(gdb, "varInfo")
  manual <- var %>%
    dplyr::filter(ModerateImpact == 1) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarize(
      CHROM = CHROM[1],
      start = min(as.numeric(POS)),
      end = max(as.numeric(POS))
    ) %>%
    dplyr::mutate(varSetName = "Moderate") %>%
    dplyr::rename(unit = gene_name) %>%
    as.data.frame()
  ranges_moderate <- ranges %>%
    dplyr::filter(varSetName == "Moderate")
  expect_equal(ranges_moderate, manual[, colnames(ranges_moderate)])

  # check if writing output works
  output <- withr::local_tempfile()
  expect_no_error(getRanges(varsetfile, gdb = gdb, output = output))
  ranges_from_file <- read.table(
    output,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  expect_equal(ranges, ranges_from_file)

  # check if running getRanges with a different table works
  uploadAnno(
    gdb,
    name = "var_test",
    value = var %>% dplyr::mutate(POS = as.numeric(POS)),
    verbose = FALSE,
    skipRemap = TRUE
  )
  ranges_var_test <- getRanges(varsetfile, gdb = gdb, table = "var_test")
  expect_equal(ranges, ranges_var_test)

  # check if using different CHROM, POS fields works
  uploadAnno(
    gdb,
    name = "var_CHR_POSITION",
    value = var %>%
      dplyr::rename(CHR = CHROM, POSITION = POS) %>%
      dplyr::mutate(POSITION = as.numeric(POSITION)),
    skipRemap = TRUE,
    verbose = FALSE,
    overWrite = TRUE
  )
  ranges_chr_position <- getRanges(
    varsetfile,
    gdb = gdb,
    table = "var_CHR_POSITION",
    CHROM = "CHR",
    POS = "POSITION"
  )
  expect_equal(ranges, ranges_chr_position)

  # test `where` clause
  var_keep <- rbind(
    var %>% dplyr::mutate(keep = "n"),
    var %>% dplyr::mutate(keep = "y")
  )
  uploadAnno(
    gdb,
    name = "var_keep",
    value = var_keep %>% dplyr::mutate(POS = as.numeric(POS)),
    verbose = FALSE,
    skipRemap = TRUE,
    overWrite = TRUE
  )
  ranges_var_keep <- getRanges(
    varsetfile,
    gdb = gdb,
    table = "var_keep",
    where = "keep = 'y'"
  )
  expect_equal(ranges, ranges_var_keep)

  # test a where clause per varSet
  var_gene_name <- rbind(var, var %>% dplyr::mutate(gene_name = ""))
  uploadAnno(
    gdb,
    name = "var_gene_name",
    value = var_gene_name %>% dplyr::mutate(POS = as.numeric(POS)),
    verbose = FALSE,
    skipRemap = TRUE,
    overWrite = TRUE
  )
  ## should work if filtered and equal ranges
  ranges_gene_name <- getRanges(
    varsetfile,
    gdb = gdb,
    table = "var_gene_name",
    where = sprintf("gene_name = '%s'", listUnits(varsetfile))
  )
  expect_equal(ranges, ranges_gene_name)
})

test_that("getRanges input validation works", {
  gdb <- create_example_gdb()
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))

  # expect an error if `CHROM` and `POS` fields are not specified
  var_no_chrom_pos <- getAnno(gdb, "varInfo") %>%
    dplyr::select(-CHROM, -POS)
  expect_error(
    {
      ranges <- getRanges(varsetfile, gdb = gdb, table = "var_no_chrom_pos")
    }
  )

  # expect an error with multiple records per variant
  var_multiple_records <- getAnno(gdb, "varInfo") %>%
    dplyr::slice(rep(1:dplyr::n(), each = 2))

  uploadAnno(
    gdb,
    name = "var_multiple_records",
    value = var_multiple_records %>% dplyr::mutate(POS = as.numeric(POS)),
    skipRemap = TRUE,
    overWrite = TRUE,
    verbose = FALSE
  )
  expect_error(
    {
      ranges <- getRanges(varsetfile, gdb = gdb, table = "var_multiple_records")
    },
    regexp = "More than one record per variant found"
  )

  # expect error when where clause if of incorrect length
  expect_error(
    {
      getRanges(varsetfile, gdb = gdb, where = c("condition1", "condition2"))
    },
    regexp = "should be either of length 1"
  )

  # expect error when not all required fields are present in ranges table
  var <- getAnno(gdb, "varInfo")
  uploadAnno(
    gdb,
    name = "var_noCHROM",
    value = var %>% dplyr::select(-CHROM),
    verbose = FALSE,
    skipRemap = TRUE,
    overWrite = TRUE
  )
  expect_error(
    {
      getRanges(varsetfile, gdb = gdb, table = "var_noCHROM")
    },
    regexp = "The following required fields are not present in table 'var_noCHROM': CHROM."
  )

  # check if warning is thrown if no variants are found for a varSet
  var <- getAnno(gdb, "varInfo")
  uploadAnno(
    gdb,
    name = "var_noSOD1",
    value = var %>% dplyr::filter(gene_name != "SOD1"),
    verbose = FALSE,
    skipRemap = TRUE,
    overWrite = TRUE
  )
  varsetlist_SOD1 <- getVarSet(varsetfile, unit = "SOD1")
  expect_warning(
    {
      ranges <- getRanges(varsetlist_SOD1, gdb = gdb, table = "var_noSOD1")
    },
    regexp = "no position mappings were found for any of the variants."
  )
  ## ranges should be NA
  expect_true(
    all(is.na(ranges$CHROM)) &&
      all(is.na(ranges$start)) &&
      all(is.na(ranges$end))
  )
})
