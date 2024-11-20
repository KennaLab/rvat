# tests
test_that("getRanges works", {
  
  # run getRanges on a varsetfile
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
  gdb <- create_example_gdb()
  ranges <- getRanges(varsetfile, gdb = gdb)
  
  # check whether running on varSetList vs. varSetFile results in identical output
  varsetlist <- getVarSet(varsetfile, unit = listUnits(varsetfile))
  ranges_from_varsetlist <- getRanges(varsetlist, gdb = gdb)
  
  expect_equal(ranges, ranges_from_varsetlist)
  
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
  expect_equal(ranges_moderate, manual[,colnames(ranges_moderate)])
  
  # check if using different table works 
  uploadAnno(gdb, name = "var_test", value = var %>% dplyr::mutate(POS = as.numeric(POS)), verbose = FALSE, skipRemap = TRUE)
  ranges2 <- getRanges(varsetfile, gdb = gdb, table = "var_test")
  expect_equal(ranges, ranges2)
  
  # check if using different CHROM, POS fields works 
  uploadAnno(gdb, 
             name = "var_test",
             value = var %>% dplyr::rename(CHR = CHROM, POSITION = POS) %>% dplyr::mutate(POSITION = as.numeric(POSITION)), 
             skipRemap = TRUE, 
             verbose = FALSE) 
  
  ## expect an error if `CHROM` and `POS` fields are not specified
  expect_error(
    {
      ranges2 <- getRanges(varsetfile, gdb = gdb, table = "var_test")
    }
  )
 
  ## should work when `CHROM` and `POS` fields are set
  ranges2 <- getRanges(varsetfile, gdb = gdb, table = "var_test", CHROM = "CHR", POS = "POSITION")
  expect_equal(ranges, ranges2)
  
  # test `where` clause
  var_ <- rbind(var %>% dplyr::mutate(keep = "n"), var %>% dplyr::mutate(keep = "y"))
  uploadAnno(gdb, name = "var_test", value = var_ %>% dplyr::mutate(POS = as.numeric(POS)), verbose = FALSE, skipRemap = TRUE)
  
  ## should return an error if no filter is specified
  expect_error(
    {
      ranges2 <- getRanges(varsetfile, gdb = gdb, table = "var_test")
    }, regexp = "More than one record per variant found"
  )

  ## should work if filtered and equal ranges
  ranges2 <- getRanges(varsetfile, gdb = gdb, table = "var_test", where = "keep = 'y'")
  expect_equal(ranges, ranges2)
  
  # test a where clause per varSet
  var_ <- rbind(var, var %>% dplyr::mutate(gene_name = ""))
  uploadAnno(gdb, name = "var_test", value = var_ %>% dplyr::mutate(POS = as.numeric(POS)), verbose = FALSE, skipRemap = TRUE)
  
  ## expect an error without filters
  expect_error(
    {
      ranges2 <- getRanges(varsetfile, gdb = gdb, table = "var_test")
    }, regexp = "More than one record per variant found"
  )
  
  ## should work if filtered and equal ranges
  ranges2 <- getRanges(varsetfile, gdb = gdb, table = "var_test", where = sprintf("gene_name = '%s'", listUnits(varsetfile)))
  expect_equal(ranges, ranges2)
  
  # check if no overlaps is not appropriately handled
  varsetlist <- getVarSet(varsetfile, unit = "SOD1")
  uploadAnno(gdb, name = "var_test", value = var %>% dplyr::filter(gene_name != "SOD1"), verbose = FALSE, skipRemap = TRUE)
  expect_warning(
    {
      ranges2 <- getRanges(varsetlist, gdb = gdb, table = "var_test")
    }, regexp = "no position mappings where found for any of the variants."
    )
  expect_true(all(is.na(ranges2$CHROM)) && all(is.na(ranges2$start)) && all(is.na(ranges2$end)))
  }
)
