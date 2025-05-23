gdb <- gdb(rvat_example("rvatData.gdb"))
moderate <- withr::local_tempfile()
LOF <- withr::local_tempfile()
CADD <- withr::local_tempfile()
tmpfile <- withr::local_tempfile()
tmpgdb <- withr::local_tempfile()
file.copy(rvatData::rvat_example("rvatData.gdb"), tmpgdb, overwrite = TRUE)


test_that("buildVarSet works",{
  # build varsets
  ## moderate
  expect_no_error(
    suppressMessages(buildVarSet(object = gdb, 
                output = moderate,
                varSetName = "Moderate", 
                unitTable = "varInfo", 
                unitName = "gene_name",
                where = "(ModerateImpact = 1 or HighImpact = 1) and gene_name in ('CYP19A1', 'FUS', 'OPTN', 'SOD1')")
   )
  )
  
  ## LOF
  expect_no_error(
    buildVarSet(object = gdb, 
                output = LOF,
                varSetName = "HighImpact", 
                unitTable = "varInfo", 
                unitName = "gene_name",
                where = "HighImpact = 1 and gene_name in ('CYP19A1', 'FUS', 'OPTN', 'SOD1')",
                verbose = FALSE
                )
  )
  
  ## CADD weights
  expect_no_error(
    suppressMessages(buildVarSet(object = gdb, 
                output = CADD,
                varSetName = "CADD", 
                unitTable = "varInfo", 
                unitName = "gene_name",
                weightName = "CADDphred",
                where = "gene_name in ('CYP19A1', 'FUS', 'OPTN', 'SOD1')"))
  )
  
  ## build varSets based on a data.frame
  anno <- getAnno(gdb, "varinfo")
  anno$Moderate=ifelse(anno$ModerateImpact == 1 | anno$HighImpact == 1, 1, 0)
  anno$CADD=anno$CADDphred
  
  # retrieve varsets
  moderate_df <- buildVarSet(
    anno,
    unitName = "gene_name",
    fields = c("Moderate")
  )
  LOF_df <- buildVarSet(
    anno,
    unitName = "gene_name",
    fields = c("HighImpact")
  )
  cadd_df <- buildVarSet(
    anno,
    unitName = "gene_name",
    fields = c("CADD")
  )
  
  varset1 <- getVarSet(moderate_df, unit = c("CYP19A1", "FUS", "OPTN", "SOD1"))
  varset2 <- getVarSet(varSetFile(moderate), unit = c("CYP19A1", "FUS", "OPTN", "SOD1"))
  varset1@metadata = list()
  varset2@metadata = list()
  expect_identical(
    varset1,
    varset2
  )
  varset1 <- getVarSet(LOF_df, unit = c("CYP19A1", "FUS", "OPTN", "SOD1"))
  varset2 <- getVarSet(varSetFile(LOF), unit = c("CYP19A1", "FUS", "OPTN", "SOD1"))
  varset1@metadata = list()
  varset2@metadata = list()
  expect_identical(
    varset1,
    varset2
  )
  
  varset1 <- getVarSet(cadd_df, unit = c("CYP19A1", "FUS", "OPTN", "SOD1"))
  varset2 <- getVarSet(varSetFile(CADD), unit = c("CYP19A1", "FUS", "OPTN", "SOD1"))
  varset1@metadata = list()
  varset2@metadata = list()
  expect_identical(
    varset1,
    varset2
  )
  
  # test intersect
  gdb_tmp <- gdb(tmpgdb)
  vars <- getAnno(gdb_tmp, "var", fields = "VAR_id")$VAR_id
  vars <- sort(sample(vars, size = 200))
  uploadAnno(gdb_tmp, 
             value = data.frame(VAR_id = vars), 
             name = "random_vars", 
             skipRemap = TRUE,
             verbose = FALSE
             )
  
  # 
  moderate_intersect <- tmpfile
  moderate_intersect <- buildVarSet(object = gdb_tmp, 
              varSetName = "Moderate", 
              unitTable = "varInfo", 
              unitName = "gene_name",
              intersect = "random_vars",
              where = "(ModerateImpact = 1 or HighImpact = 1) and gene_name in ('CYP19A1', 'FUS', 'OPTN', 'SOD1')")
  vars_check <- listVars(collapseVarSetList(moderate_intersect))
  expect_true(all(vars_check %in% vars))
})


## combine 
varsets <- list()
i <- 1
for(file in c(moderate,LOF,CADD)) {
  varset <- readr::read_delim(file, delim="|",col_names=FALSE,col_types="cccc", comment = "#")
  varsets[[i]] <- varset
  i <- i+1
}
varsets <- dplyr::bind_rows(varsets)
varsets <- varsets %>% 
  dplyr::arrange(X1)
metadata <- metadata(varSetFile(moderate))
merged <- tmpfile
con <- gzfile(merged, "w") 
rvat:::.write_rvat_header(filetype = "varSetFile", metadata = metadata, con = con)
write.table(varsets, file = con, col.names = FALSE, sep="|", append=TRUE, quote = FALSE, row.names = FALSE)
close(con)


test_that("getVarSet works",{
  # build varsets
  ## get some varsets from varsetfile
  varsetfile <- varSetFile(merged)
  varsets <- getVarSet(
    varsetfile,
    unit = c("FUS", "OPTN"),
    varSetName = "Moderate"
  )
  varsetfile_moderate <- varSetFile(moderate)
  varsets_moderate <- getVarSet(
    varsetfile_moderate,
    unit = c("FUS", "OPTN")
  )
  expect_identical(varsets, varsets_moderate)
  
  ## extract from varsetlist
  varsetlist <- getVarSet(varsetfile, unit = listUnits(varsetfile))
  varsets2 <- getVarSet(varsetlist,
                        unit = c("FUS", "OPTN"),
                        varSetName = "Moderate"
                        )
  expect_identical(varsets, varsets2)
  
  
  ## try extracting varSets that are not present in varSet
  expect_warning({
    varsets2 <- getVarSet(varsetlist,
                          unit = c("FUS", "OPTN"),
                          varSetName = "moderate")
  }, regexp = "Not all specified varSets are present in the varSetList"
  )
  
  ## try extracting units that are not present in varSet
  expect_warning({
   varsets2 <- getVarSet(varsetlist,
              unit = c("FUS", "OPTN", "abc"),
              varSetName = "Moderate")
  }, regexp = "Not all specified units are present in the varSetList"
  )
  expect_identical(varsets, varsets2)
  
  ## extract in different order
  varsets2 <- getVarSet(varsetlist,
                        unit = c("OPTN", "FUS"),
                        varSetName = "Moderate")
  expect_identical(varsets, varsets2)
  
  ## compare subsetting 
  varsets3=varsetlist[which(listUnits(varsetlist) %in% c("FUS", "OPTN") & listVarSets(varsetlist) %in% c("Moderate"))]
  expect_identical(varsets2, varsets3)
  
})

test_that("writing a varsetList works",{
  varsetfile <- varSetFile(merged)
  varsetlist <- getVarSet(varsetfile, unit = listUnits(varsetfile))
  write(varsetlist, tmpfile)
  varsetfile2 <- varSetFile(tmpfile)
  varsetlist1 <- getVarSet(varsetfile,unit=listUnits(varsetfile))
  varsetlist2 <- getVarSet(varsetfile2,unit=listUnits(varsetfile2))
  ## ignore creation dates, since they'll differ
  varsetlist1@metadata$creationDate <- NA_character_
  varsetlist2@metadata$creationDate <- NA_character_
  expect_identical(varsetlist1, varsetlist2)
})

test_that("showing a varSetFile works", {
  expect_true(stringr::str_detect(capture_output({show(varSetFile(merged))}), "varSetFile object"))
}
)

test_that("showing a varSetList works", {
  expect_true(stringr::str_detect(capture_output({show(getVarSet(varSetFile(merged), unit = "SOD1"))}), "varSetList"))
  }
)

test_that("showing a varSet works", {
  expect_true(stringr::str_detect(capture_output({show(getVarSet(varSetFile(merged), unit = "SOD1")[[1]])}), "unit="))
}
)

test_that("varSet methods are backward compatible with RVAT versions <=0.2.10", {
  # write varsetfile without header
  varsetlist <- as.data.frame(getVarSet(varSetFile(moderate), unit = listUnits(varSetFile(moderate))))
  write.table(varsetlist, file = gzfile(tmpfile), col.names = FALSE, sep="|", append=FALSE, quote = FALSE, row.names = FALSE)
  
  ## connect
  expect_warning({vfile <- varSetFile(tmpfile)}, regexp = "The following metadata fields are missing")
  
  ## retrieve varsets 
  varsets <- getVarSet(vfile, unit = listUnits(vfile)[c(3,1)])
  
  # write varsetfile with unexpected metadata
  con <- gzfile(tmpfile, "w") 
  metadata <- list(
    test = "test"
  )
  rvat:::.write_rvat_header(filetype = "varSetFile", metadata = metadata, con = con)
  write.table(varsetlist, file = con, col.names = FALSE, sep="|", append=TRUE, quote = FALSE, row.names = FALSE)
  close(con)
  
  ## connect
  expect_error({vfile <- varSetFile(tmpfile)}, regexp = "unexpected")
  
  # write varsetfile with unexpected metadata (2)
  con <- gzfile(tmpfile, "w") 
  for (item in names(metadata)) {
    writeLines(sprintf("# %s: %s: %s", item, metadata[[item]], metadata[[item]]), con = con)
  }
  write.table(varsetlist, file = con, col.names = FALSE, sep="|", append=TRUE, quote = FALSE, row.names = FALSE)
  close(con)
  
  ## connect
  expect_error({vfile <- varSetFile(tmpfile)}, regexp = "Unexpected")
}
)

# getRanges
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

