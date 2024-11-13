# simple test: build a gdb and expect no errors
test_that("buildGdb works", {
  expect_no_error(suppressMessages(buildGdb(
    vcf = rvat_example("rvatData.vcf.gz"),
    output = testgdb,
    genomeBuild = "GRCh38"
  )))
  expect_no_error({gdb <- gdb(testgdb)})
})

# write gdb to vcf, then build a gdb from this vcf, which should result in an identical gdb
test_that("writeVcf works",{
  
  ## write gdb to vcf
  output_writevcf <- withr::local_tempfile()
  suppressMessages(writeVcf(
    gdb(rvat_example("rvatData.gdb")),
    output_writevcf
  ))

  ## build a gdb again from vcf written from gdb
  output_writevcf_buildgdb <- withr::local_tempfile()
  suppressMessages(gdb2 <- buildGdb(
    vcf = output_writevcf,
    output = output_writevcf_buildgdb,
    overWrite=TRUE,
    genomeBuild = "GRCh38"
  ))

  gdb1 <- gdb(rvat_example("rvatData.gdb"))
  gdb2 <- gdb(output_writevcf_buildgdb)
  cohort <- getCohort(gdb1, "pheno")
  cohort <- cohort[!is.na(cohort$IID),]
  suppressMessages(upload <- uploadCohort(gdb2, name = "pheno", value = cohort))
  GT1 <- getGT(gdb1,
               VAR_id = getAnno(gdb1, table = "var", fields = "VAR_id")$VAR_id,
               cohort = "pheno",
               verbose = FALSE
  )
  GT2 <- getGT(gdb2,
               VAR_id = getAnno(gdb2, table = "var", fields = "VAR_id")$VAR_id,
               cohort = "pheno",
               verbose = FALSE
  )
  ### ignore gdb and gdbId in comparison
  metadata(GT1)$gdb <- NA_character_
  metadata(GT2)$gdb <- NA_character_
  metadata(GT1)$gdbId <- NA_character_
  metadata(GT2)$gdbId <- NA_character_
  expect_identical(GT1,GT2)
  
  ### also expect varinfo to be identical
  var1 <- getAnno(gdb1, "var")
  var1$INFO <- "."
  var2 <- getAnno(gdb2, "var")
  expect_identical(var1, var2) 
  
  # also check w/o including genotypes
  suppressMessages(writeVcf(object = gdb(rvat_example("rvatData.gdb")),
           includeGeno = FALSE,
           output_writevcf
           ))
  var2 <- readr::read_tsv(output_writevcf, show_col_types = FALSE, col_names = FALSE, comment = "#")
  colnames(var2) <- colnames(var1)[2:9]
  var1$QUAL <- as.numeric(var1$QUAL)
  expect_equal(var1[,2:9], as.data.frame(var2))

})

# test whether subsetgdb works by splitting gdb in three subsets, and then glueing them back together with concatGdb
test_that("subsetGdb and concatGdb work",{
  # split gdb using subsetGdb and then concat using concatGdb -> original gdb and concatted gdb should be identical
  gdb <- gdb(rvat_example("rvatData.gdb"))
  vars <- getAnno(gdb, table = "var", fields = "VAR_id")$VAR_id
  vars <- split(vars, rep(1:3, each = ceiling(length(vars)/3), length.out = length(vars)))
  tmpfiles <- c(
    withr::local_tempfile(),
    withr::local_tempfile(),
    withr::local_tempfile()
  )
  files <-  withr::local_tempfile()
  readr::write_lines(tmpfiles, file = files)
  for(i in 1:3) {
    suppressMessages(subsetGdb(
      gdb,
      output = tmpfiles[i],
      where = sprintf("VAR_id in (%s)", paste(paste0("'", vars[[i]], "'"), collapse = ",")),
      overWrite = TRUE
    ))
  }
  ## expect error when writing to same file again when `overWrite = FALSE`
  expect_error({subsetGdb(
    gdb,
    output = tmpfiles[i],
    where = sprintf("VAR_id in (%s)",paste(paste0("'", vars[[i]], "'"),collapse=","))
  )}, "already exists")
  
  ## concat three parts
  gdbconcat <- withr::local_tempfile()
  suppressMessages(concatGdb(files, output = gdbconcat))
  gdbconcat <- gdb(gdbconcat)
  
  ## compare gdbs
  GT1 <- getGT(gdb,
               VAR_id = getAnno(gdb, table = "var", fields = "VAR_id")$VAR_id,
               verbose = FALSE
  )
  GT2 <- getGT(gdbconcat,
               VAR_id = getAnno(gdbconcat, table = "var", fields = "VAR_id")$VAR_id,
               verbose = FALSE
  )
  ### ignore gdb and gdbId in comparison
  metadata(GT1)$gdb <- NA_character_
  metadata(GT2)$gdb <- NA_character_
  metadata(GT1)$gdbId <- NA_character_
  metadata(GT2)$gdbId <- NA_character_
  
  expect_identical(GT1,GT2)
  
  # additional tests subsetGdb
  gdb_subset <- withr::local_tempfile()
  suppressMessages(subsetGdb(gdb,
            where = "CHROM = 'chr10'", 
            output = gdb_subset, 
            overWrite = TRUE))
  gdb_subset <- gdb(gdb_subset)
  var1 <- getAnno(gdb_subset, "var")
  expect_equal(unique(var1$CHROM), "chr10")
  
  
  # do the same using `VAR_id` parameter
  gdb_subset2 <- withr::local_tempfile()
  var <- getAnno(gdb, "var") %>%
    dplyr::filter(CHROM == "chr10")
  suppressMessages(subsetGdb(gdb,
            VAR_id = var$VAR_id,
            output = gdb_subset2, 
            overWrite = TRUE))
  gdb_subset2 <- gdb(gdb_subset2)
  var2 <- getAnno(gdb_subset2, "var")
  rownames(var1) = NULL
  var2 <- var2[match(var1$VAR_id,var2$VAR_id),]
  rownames(var2) = NULL
  expect_identical(var1, var2)
  
  # check combination of where and VAR_id
  gdb_subset3 <- withr::local_tempfile()
  var <- getAnno(gdb, "var")
  var_chr <- var %>% 
    dplyr::filter(CHROM %in% c("chr10", "chr9", "chr15"))
  suppressMessages(subsetGdb(gdb,
            where = "CHROM in ('chr10','chr14','chr9')", 
            VAR_id = var_chr$VAR_id,
            output = gdb_subset3, 
            overWrite = TRUE))
  gdb_subset3 <- gdb(gdb_subset3)
  var3 <- getAnno(gdb_subset3, "var")
  expect_equal(sort(unique(var3$CHROM)), c("chr10", "chr9"))
  
  # additional tests concatGdb
  files <- withr::local_tempfile()
  readr::write_lines(c(tmpfiles[1], tmpfiles), file = files)
  expect_error(concatGdb(
    targets = files,
    output = withr::local_tempfile()
  ))
})


test_that("uploading and retrieving annotations works",{
  # subset gdb in three parts
  gdb <- gdb(testgdb)
  suppressMessages(uploadAnno(gdb, 
             name = "varinfo", 
             value = rvat_example("rvatData.varinfo")))
  
  a <- data.frame(
    VAR_id = 1:10,
    gene = LETTERS[10:1],
    width = c(5,9,59,65,39,10,12,94,52,11)
  )
  
  b <- data.frame(
    VAR_id = 5:10,
    gene = LETTERS[6:1],
    transcript = LETTERS[16:21]
  )

  suppressMessages(uploadAnno(gdb, name = "a", value = a, skipRemap = TRUE))
  suppressMessages(uploadAnno(gdb, name = "b", value = b, skipRemap = TRUE))

  query <- getAnno(gdb, "a", inner = "b")
  expect_equal(query$VAR_id, 5:10)
  expect_equal(query[,2], query[,4])

  # check with/without excluding unmapped variants
  var <- getAnno(gdb,"var", VAR_id = c(5,6,100,193,900))
  a <- data.frame(
    CHROM = c(var$CHROM, c("chrX", "chrY")),
    POS = c(var$POS, 1000,10000),
    REF = c(var$REF, "A", "C"),
    ALT = c(var$ALT,"G", "T"),
    gene = LETTERS[7:1],
    width = c(5,9,59,65,39,10,12)
  )
  expect_warning(uploadAnno(gdb, name = "c", value = a, verbose = FALSE))
  expect_warning(uploadAnno(gdb, name = "d", value = a, keepUnmapped = TRUE, verbose = FALSE))
  c <- getAnno(gdb, "c")
  d <- getAnno(gdb, "d")
  d <- d[!is.na(d$VAR_id),]
  expect_equal(c,d)

  # check with/without ignoring alleles
  expect_warning(suppressMessages(uploadAnno(gdb, name = "c", value = a, ignoreAlleles = TRUE)))
  expect_warning(suppressMessages(uploadAnno(gdb, name = "d", value = a, ignoreAlleles = TRUE, keepUnmapped = TRUE)))
  c <- getAnno(gdb, "c")
  d <- getAnno(gdb, "d")
  
  # expect errors, 
  
  ## when name is protected
  expect_error({uploadAnno(gdb, 
                           value = data.frame(VAR_id = c("1", "2"),
                                              anno_field = c("A", "B")),
                           name = "SM",
                           skipRemap = TRUE
                           )},
               regexp = "already exists as a protected table")
  expect_error({uploadAnno(gdb, 
                           value = data.frame(VAR_id = c("1", "2"),
                                              anno_field = c("A", "B")),
                           name = "meta",
                           skipRemap = TRUE
                           )},
               regexp = "already exists as a protected table")
  expect_error({uploadAnno(gdb, 
                           value = data.frame(VAR_id = c("1", "2"),
                                                anno_field = c("A", "B")),
                           name = "anno",
                           skipRemap = TRUE
                           )},
               regexp = "already exists as a protected table")
  expect_error({uploadAnno(gdb, 
                           value = data.frame(VAR_id = c("1", "2"),
                                              anno_field = c("A", "B")),
                           name = "var_ranges",
                           skipRemap = TRUE
                           )},
               regexp = "already exists as a protected table")
  expect_error({uploadAnno(gdb, 
                           value = data.frame(VAR_id = c("1", "2"),
                                              anno_field = c("A", "B")),
                           name = "dosage",
                           skipRemap = TRUE
                           )},
               regexp = "already exists as a protected table")
  
  ## when name is already included as pheno table
  uploadCohort(gdb, 
               value = data.frame(IID = c("ALS1", "ALS2"),
                                  sex = c(1, 2)),
               name = "pheno_tmp",
               verbose = FALSE
               )
  expect_error({uploadAnno(gdb, 
                          value = data.frame(VAR_id = c("1", "2"),
                                                anno_field = c("A", "B")),
                          name = "pheno_tmp",
                          skipRemap = TRUE
                          )},
               regexp = "already exists as a cohort table")
  
  ## when upload name includes punctuation marks
  for (mark in c(".", ",", "+", "-", " ")) {
    expect_error({uploadAnno(gdb, 
                             value = data.frame(VAR_id = c("1", "2"),
                                                anno_field = c("A", "B")),
                             name = sprintf("anno%stest", mark),
                             skipRemap = TRUE
                             )},
                 regexp = "Table name shouldn't contain")
    }
  }
)

test_that("getGT works",{
  ## GT with 100 variants for testing
  gdb <- gdb(rvat_example("rvatData.gdb"))
  GT <- getGT(gdb,
              VAR_id = 1:100,
              cohort = "pheno",
              verbose = FALSE
              )
  
  ## order of VAR_ids shouldn't matter
  GT_check <- getGT(gdb,
              VAR_id = sample(1:100,size=100),
              cohort = "pheno",
              verbose = FALSE
              )
  expect_identical(GT,GT_check)
  
  ## loading VAR_ids that are not present in gdb should yield a warning
  expect_warning(
    {
    GT_check <- getGT(gdb,
                      VAR_id = c(1:100, 6000, 6001),
                      cohort = "pheno",
                      verbose = FALSE
                      )
    }
  )
  expect_identical(GT,GT_check)
  
  
  ## check weights
  CADD <- withr::local_tempfile(fileext = ".txt.gz")
  
  # Build a varset containing CADD weights
  varset <- buildVarSet(object = gdb, 
              output = CADD,
              varSetName = "CADD", 
              unitTable = "varInfo", 
              unitName = "gene_name",
              weightName = "CADDphred",
              where = "gene_name in ('CYP19A1', 'FUS', 'OPTN')",
              verbose = FALSE
              )
  varset <- varSetFile(CADD)
  GT_check <- suppressWarnings(getGT(gdb, varSet = getVarSet(varset, unit = "FUS"), cohort = "pheno", verbose = FALSE))
  
  # Compare with manually extracting and adding CADD scores
  varinfo <- getAnno(gdb, "varinfo", VAR_id=listVars(getVarSet(varset, unit = "FUS")[[1]]))
  GT_check2 <- getGT(gdb,cohort="pheno", VAR_id=listVars(getVarSet(varset, unit = "FUS")[[1]]),varSetName="CADD",unit="FUS",verbose=FALSE)
  caddscores <- varinfo$CADDphred
  GT_check2 <- suppressWarnings(recode(GT_check2, weights = as.numeric(varinfo$CADDphred)))
  
  expect_identical(GT_check,GT_check2)
  
  ## check different ways of specifying a cohort
  GT <- getGT(gdb,
              VAR_id = 1:100,
              cohort = "pheno",
              verbose = FALSE
              )
  chrt <- getCohort(gdb, "pheno")
  GT2 <- getGT(gdb,
               VAR_id = 1:100,
               cohort = chrt,
               verbose = FALSE
               )
  chrt <- S4Vectors::DataFrame(chrt)
  S4Vectors::metadata(chrt)$name <- "pheno"
  GT3 <- getGT(gdb,
               VAR_id = 1:100,
               cohort = chrt,
               verbose = FALSE
               )
  metadata(GT2)$cohort <- "pheno"
  expect_equal(GT,GT2)
  expect_equal(GT,GT3)
  
  ## check specifying anno 
  GT1 <- getGT(gdb,
              VAR_id = 1:100,
              cohort = "pheno",
              verbose = FALSE
              )
  GT2 <- getGT(gdb,
              VAR_id = 1:100,
              cohort = "pheno",
              anno = "var",
              verbose = FALSE
              )
  expect_gt(ncol(rowData(GT2)), ncol(rowData(GT1)))
  
  GT3 <- getGT(gdb,
               VAR_id = 1:100,
               cohort = "pheno",
               anno = "var",
               annoFields = c("VAR_id", "CHROM", "POS", "ID", "REF", "ALT"),
               verbose = FALSE
               )
  expect_gt(ncol(rowData(GT2)), ncol(rowData(GT3)))
  
  ## expect error when anno table contains more than one line per VAR_id
  gdb_tmp <- withr::local_tempfile()
  file.copy(rvat_example("rvatData.gdb"), gdb_tmp, overwrite = TRUE)
  gdb_tmp <- gdb(gdb_tmp)
  uploadAnno(gdb_tmp, 
             name = "transcriptInfo", 
             value = data.frame(VAR_id = c(1, 1, 2, 2, 2, 3),
                                transcript = c("A", "B", "C", "D", "E", "F")),
             skipRemap = TRUE,
             verbose = FALSE
             )
  expect_error({
    GT <- getGT(gdb_tmp,
              VAR_id = 1:100,
              cohort = "pheno",
              anno = "transcriptInfo")
  })
  
  ## expect error when a varsetlist with >1 varSets is provided
  expect_error({
    GT <- getGT(gdb_tmp,
                VAR_id = 1:100,
                cohort = "pheno",
                varSet = getVarSet(varset, listUnits(varset))
                )
    }, regexp = "Only 1 varSet can be supplied")
  }
)

# check uploadCohort/getCohort

test_that("uploading and retrieving cohorts works",{
  gdb <- gdb(testgdb)
  pheno <- readr::read_tsv(rvat_example("rvatData.pheno"), show_col_types = FALSE)
  
  # upload a cohort that includes all IIDs in gdb
  uploadCohort(gdb, value = as.data.frame(pheno), name = "pheno", verbose = FALSE)
  
  # upload directly from file
  uploadCohort(gdb, value = rvat_example("rvatData.pheno"), name = "pheno_fromfile", verbose = FALSE)
  
  ## compare 
  pheno_retrieved <- getCohort(gdb, "pheno")
  expect_equal(as.data.frame(pheno), pheno_retrieved)
  pheno_fromfile_retrieved <- getCohort(gdb, "pheno_fromfile")
  expect_equal(as.data.frame(pheno), pheno_fromfile_retrieved, tolerance = 1e-8)
  
  # upload a cohort that contains subset of IIDs
  pheno_subset <- pheno[c(1:1000, 5000:6000),]
  uploadCohort(gdb, value = as.data.frame(pheno_subset), name = "pheno_subset", verbose = FALSE)
  pheno_subset_file <- withr::local_tempfile()
  readr::write_tsv(pheno_subset, file = pheno_subset_file)
  uploadCohort(gdb, value = pheno_subset_file, name = "pheno_subset_fromfile", verbose = FALSE)
  
  ## compare 
  pheno_subset_retrieved <- getCohort(gdb, "pheno_subset")
  pheno_subset_retrieved <- pheno_subset_retrieved[!is.na(pheno_subset_retrieved$IID),]
  pheno_subset_fromfile_retrieved <- getCohort(gdb, "pheno_subset_fromfile")
  pheno_subset_fromfile_retrieved <- pheno_subset_fromfile_retrieved[!is.na(pheno_subset_fromfile_retrieved$IID),]
  pheno_subset_fromfile_retrieved$pheno <- as.numeric(pheno_subset_fromfile_retrieved$pheno) # integer/double differences
  rownames(pheno_subset) <- NULL
  rownames(pheno_subset_retrieved) <- NULL
  rownames(pheno_subset_fromfile_retrieved) <- NULL
  expect_equal(as.data.frame(pheno_subset), pheno_subset_retrieved)
  expect_equal(as.data.frame(pheno_subset), pheno_subset_fromfile_retrieved)
  
  # upload a cohort that contains IIDs not present in gdb
  pheno_nonmatchingiids <- dplyr::bind_rows(
    pheno[c(1:1000, 5000:6000),],
    pheno[1:100,] %>% dplyr::mutate(IID = paste0("test", 1:100))
  )
  suppressMessages(uploadCohort(gdb, value = as.data.frame(pheno_nonmatchingiids), name = "pheno_subset2"))
  
  ## compare 
  pheno_nonmatchingiids_retrieved <- getCohort(gdb, "pheno_subset2")
  pheno_nonmatchingiids_retrieved <- pheno_nonmatchingiids_retrieved[!is.na(pheno_nonmatchingiids_retrieved$IID),]
  rownames(pheno_nonmatchingiids_retrieved) <- NULL
  expect_failure(expect_equal(pheno_nonmatchingiids, pheno_nonmatchingiids_retrieved))
  expect_equal(as.data.frame(pheno_subset), pheno_nonmatchingiids_retrieved)

  # check subsetting fields
  pheno <- getCohort(gdb, "pheno", fields = c("IID", "sex", "pheno", "pop"))
  expect_equal(colnames(pheno), c("IID", "sex", "pheno", "pop"))
  
  # expect errors, 
  
  ## when subsetting non-existing fields
  expect_error(pheno <- getCohort(gdb, "pheno", fields = c("IID", "sex", "pheno", "pop", "random")))
  
  ## when duplicated IID values are included in upload
  expect_error(uploadCohort(gdb, 
               value = data.frame(IID = c("ALS1", "ALS1", "ALS2"),
                                  sex = c(1, 1, 2)),
               name = "pheno_dups",
               verbose = FALSE
               ))
  
  ## when name is protected
  expect_error({uploadCohort(gdb, 
                             value = data.frame(IID = c("ALS1", "ALS2"),
                                                sex = c(1, 2)),
                             name = "SM")},
               regexp = "already exists as a protected table")
  expect_error({uploadCohort(gdb, 
                             value = data.frame(IID = c("ALS1", "ALS2"),
                                                sex = c(1, 2)),
                             name = "meta")},
               regexp = "already exists as a protected table")
  expect_error({uploadCohort(gdb, 
                             value = data.frame(IID = c("ALS1", "ALS2"),
                                                sex = c(1, 2)),
                             name = "anno")},
               regexp = "already exists as a protected table")
  expect_error({uploadCohort(gdb, 
                             value = data.frame(IID = c("ALS1", "ALS2"),
                                                sex = c(1, 2)),
                             name = "var_ranges")},
               regexp = "already exists as a protected table")
  expect_error({uploadCohort(gdb, 
                             value = data.frame(IID = c("ALS1", "ALS2"),
                                                sex = c(1, 2)),
                             name = "dosage")},
               regexp = "already exists as a protected table")
  
  ## when name is already included as anno table
  expect_error({uploadCohort(gdb, 
                             value = data.frame(IID = c("ALS1", "ALS2"),
                                                sex = c(1, 2)),
                             name = "dosage")},
               regexp = "already exists as a protected table")
  
  ## when upload name includes punctuation marks
  expect_error({uploadCohort(gdb, 
                            value = data.frame(IID = c("ALS1", "ALS2"),
                                               sex = c(1, 2)),
                            name = "pheno.test")},
               regexp = "Table name shouldn't contain")
  expect_error({uploadCohort(gdb, 
                             value = data.frame(IID = c("ALS1", "ALS2"),
                                                sex = c(1, 2)),
                             name = "pheno,test")},
               regexp = "Table name shouldn't contain")
  expect_error({uploadCohort(gdb, 
                             value = data.frame(IID = c("ALS1", "ALS2"),
                                                sex = c(1, 2)),
                             name = "pheno+test")},
               regexp = "Table name shouldn't contain")
  expect_error({uploadCohort(gdb, 
                             value = data.frame(IID = c("ALS1", "ALS2"),
                                                sex = c(1, 2)),
                             name = "pheno-test")},
               regexp = "Table name shouldn't contain")
}
)

test_that("dropTable works",{
  # upload a table
  gdb <- gdb(testgdb)
  all_tables <- RSQLite::dbListTables(gdb)
  uploadAnno(gdb, 
             name = "test_table", 
             value = data.frame(VAR_id = c(1, 1, 2, 2, 2, 3),
                                transcript = c("A", "B", "C", "D", "E", "F")),
             skipRemap = TRUE,
             verbose = FALSE
  )
  all_tables2 <- RSQLite::dbListTables(gdb)
  expect_equal(setdiff(all_tables2, all_tables), "test_table")
  
  # drop the table
  suppressMessages(dropTable(gdb, "test_table"))
  all_tables3 <- RSQLite::dbListTables(gdb)
  expect_equal(all_tables, all_tables3)
 
}
)


# misc methods
test_that("misc gdb methods work",{
  gdb <- gdb(testgdb)
  
  # show gdb connection
  expect_true(stringr::str_detect(capture_output({show(gdb)}), "gdb object"))
  
  # get rvat version
  expect_equal(getRvatVersion(gdb), as.character(packageVersion("rvat")))
  
  # get creation date
  expect_true((stringr::str_detect(getCreationDate(gdb), "-") && stringr::str_detect(getCreationDate(gdb), ":")))
  
  # expect error when connecting to nonexisting file
  expect_error({gdb(paste(sample(LETTERS, 10), collapse = ""))}, "doesn't exist")
  }
)

# check summariseGeno
test_that("summariseGeno works" ,{
  # compare summariseGeno-gdb and summariseGeno-genoMatrix methods
  gdb <- gdb(rvat_example("rvatData.gdb"))
  varids <- getAnno(gdb, "var", fields = "VAR_id")$VAR_id[1:500]
  GT <- getGT(
    gdb,
    VAR_id = varids,
    cohort = "pheno",
    verbose = FALSE
  )
  sumgeno <- summariseGeno(GT)
  sumgeno_output <- withr::local_tempfile()
  suppressMessages(summariseGeno(
    gdb,
    varSet = .varsTovarSetList(varids),
    cohort = "pheno",
    output = sumgeno_output
  ))
  sumgeno_fromgdb <- as.data.frame(readr::read_tsv(sumgeno_output, show_col_types = FALSE))
  rownames(sumgeno) <- NULL
  rownames(sumgeno_fromgdb) <- NULL
  sumgeno$VAR_id <- as.character(sumgeno$VAR_id)
  sumgeno_fromgdb$VAR_id <- as.character(sumgeno_fromgdb$VAR_id)
  expect_equal(sumgeno, sumgeno_fromgdb)
  
  # check `splitBy` argument
  suppressMessages(summariseGeno(
    gdb,
    varSet = .varsTovarSetList(varids),
    cohort = "pheno",
    splitBy = "pheno",
    output = sumgeno_output
  ))
  sumgeno_fromgdb <- as.data.frame(readr::read_tsv(sumgeno_output, show_col_types = FALSE))
  sumgeno <- 
    rbind(summariseGeno(GT[,GT$pheno == 1]) %>% dplyr::mutate(pheno = 1) %>% dplyr::select(VAR_id,pheno,dplyr::everything()),
          summariseGeno(GT[,GT$pheno == 0]) %>% dplyr::mutate(pheno = 0) %>% dplyr::select(VAR_id,pheno,dplyr::everything()))
  rownames(sumgeno) <- NULL
  rownames(sumgeno_fromgdb) <- NULL
  sumgeno$VAR_id <- as.character(sumgeno$VAR_id)
  sumgeno_fromgdb$VAR_id <- as.character(sumgeno_fromgdb$VAR_id)
  expect_equal(sumgeno, sumgeno_fromgdb) 
  
  # check runnin summariseGeno on a varSetFile
  varsetfile_path <- withr::local_tempfile()
  varsetfile <- buildVarSet(object = gdb, 
              output = varsetfile_path,
              varSetName = "Moderate", 
              unitTable = "varInfo", 
              unitName = "gene_name",
              where = "(ModerateImpact = 1 or HighImpact = 1) and gene_name in ('CYP19A1', 'FUS', 'OPTN', 'SOD1')",
              verbose = FALSE
              )
  varsetfile <- varSetFile(varsetfile_path)
  GT <- getGT(
    gdb,
    varSet = collapseVarSetList(getVarSet(varsetfile, unit = listUnits(varsetfile))),
    cohort = "pheno",
    verbose = FALSE
  )
  sumgeno <- summariseGeno(GT)
  summariseGeno(
    gdb,
    varSet = varSetFile(varsetfile_path),
    cohort = "pheno",
    output = sumgeno_output,
    verbose = FALSE
  )
  sumgeno_fromgdb <- as.data.frame(readr::read_tsv(sumgeno_output, show_col_types = FALSE))
  sumgeno$VAR_id <- as.character(sumgeno$VAR_id)
  sumgeno_fromgdb$VAR_id <- as.character(sumgeno_fromgdb$VAR_id)
  sumgeno_fromgdb <- sumgeno_fromgdb[match(sumgeno$VAR_id, sumgeno_fromgdb$VAR_id),]
  rownames(sumgeno) <- NULL
  rownames(sumgeno_fromgdb) <- NULL
  expect_equal(sumgeno, sumgeno_fromgdb)
  
  # check whether strict testing works by supplying a different gdb
  expect_error({summariseGeno(
    gdb(testgdb),
    varSet = varSetFile(varsetfile_path),
    cohort = "pheno",
    output = sumgeno_output,
    verbose = FALSE
  )}, regexp = "The varSetFile seems to be generated from a different gdb than supplied.")
  
  ## when setting strict = FALSE, no error should be raised
  expect_no_error({summariseGeno(
    gdb(testgdb),
    varSet = varSetFile(varsetfile_path),
    cohort = "SM",
    output = sumgeno_output,
    verbose = FALSE,
    strict = FALSE
  )})
  
  # also check for assocTest -> note: should move to separate test block
  expect_error({assocTest(
    gdb(testgdb),
    varSet = varSetFile(varsetfile_path),
    cohort = "SM",
    pheno = "sex",
    test = "glm",
    output = sumgeno_output,
    verbose = FALSE
  )}, regexp = "The varSetFile seems to be generated from a different gdb than supplied.")
  
  ## when setting strict = FALSE, no error should be raised
  varsetfile <- varSetFile(varsetfile_path)
  expect_no_error({suppressWarnings(assocTest(
    gdb(testgdb),
    varSet = getVarSet(varsetfile, listUnits(varsetfile)[1]),
    cohort = "SM",
    pheno = "sex",
    test = "glm",
    output = sumgeno_output,
    verbose = FALSE,
    strict = FALSE
  ))})
  }
)

# check extractRanges
test_that("extractRanges works" ,{
  # extract a couple of genes
  gdb <- gdb(rvatData::rvat_example("rvatData.gdb"))
  gtf <- readr::read_tsv("../data/Homo_sapiens.GRCh38.105.gene.txt.gz", 
                         show_col_types = FALSE)
  colnames(gtf)[1] <- "CHROM"
  
  ## extract three genes, check if matches varinfo
  gtf <- gtf[gtf$gene_name %in% c("RIN3", "FUS", "ZNF483"),]
  check <- extractRanges(gdb, ranges = gtf)
  check_varinfo <- getAnno(gdb, "varInfo", where = "gene_name in ('RIN3', 'FUS', 'ZNF483')")
  expect_identical(sort(check), sort(as.character(check_varinfo$VAR_id) ))
  
  ## data.frame and GRanges should result in identical output
  check2 <- extractRanges(gdb, ranges = GenomicRanges::makeGRangesFromDataFrame(gtf) )
  expect_identical(check, check2)
  
  ## 'chr' prefix shouldn't matter
  check3 <- extractRanges(gdb, ranges = gtf %>% dplyr::mutate(CHROM = paste0("chr", CHROM)) )
  expect_identical(check, check3)
  
  # check whether INDELs are extracted correctly
  ## includes three deletions that should overlap the start coordinates of ranges
  ranges <- data.frame(
    CHROM = c("chr1", "chr4", "chr16", "chrX"),
    start = c(13112464 + 1, 169433612 + 1,31184283 + 1, 1382393 + 1),
    end = c(13112464 + 1, 169433612 + 1, 31184283 + 1, 1382393 + 1),
    VAR_id = c(644, 1362, 1206, 1765)
  )
  
  ## without padding, none should be extracted
  check <- extractRanges(gdb, 
                       ranges = ranges,
                       padding = 0)
  expect_identical(check, vector(mode = "character", length = 0))
  check <- extractRanges(gdb, 
                         ranges = ranges,
                         padding = 10)
  expect_identical(check, c("1362", "1206", "1765"))
  
  # check couple of regions by comparing to running manually with GenomicRanges
  ranges <- data.frame(
    CHROM = c("chr21", "chr21",  "chr4", "chrX"),
    start = c(31659820, 31668483, 169433614, 1356275),
    end = c(31666490, 31668554, 169438160, 1365198)
  )
  
  check1 <- extractRanges(gdb, ranges = ranges)
  var <- getAnno(gdb, "var", fields = c("VAR_id", "CHROM", "POS", "REF", "ALT"))
  var_granges <- GenomicRanges::GRanges(
    seqnames = var$CHROM,
    ranges = IRanges::IRanges(
      start = var$POS,
      end = var$POS + stringr::str_length(var$REF) - 1
    ),
    VAR_id = var$VAR_id
  )
  check2 <- subsetByOverlaps(var_granges,  GenomicRanges::makeGRangesFromDataFrame(ranges))$VAR_id
  expect_identical(sort(as.numeric(check1)), sort(as.numeric(check2)))
  }
)
