test_that("getGT basic functionality works", {
  # GT with 100 variants for testing
  gdb <- gdb(rvatData::rvat_example("rvatData.gdb"))
  GT <- expect_no_error(getGT(
    gdb,
    VAR_id = 1:100,
    cohort = "pheno",
    verbose = FALSE
  ))

  # check basic components
  expect_s4_class(GT, "genoMatrix")
  expect_true(is(GT, "genoMatrix"))
  expect_equal(nrow(GT), 100L)
  expect_equal(ncol(GT), 25000L)
  expect_equal(colnames(rowData(GT)), c("ploidy", "w"))

  # order of VAR_ids shouldn't matter
  GT_shuffled <- getGT(
    gdb,
    VAR_id = sample(1:100, size = 100),
    cohort = "pheno",
    verbose = FALSE
  )
  expect_identical(GT, GT_shuffled)
})


test_that("getGT handles invalid VAR_ids correctly", {
  gdb <- gdb(rvatData::rvat_example("rvatData.gdb"))

  # loading VAR_ids that are not present in gdb should yield a warning
  GT_reference <- getGT(gdb, VAR_id = 1:100, cohort = "pheno", verbose = FALSE)
  expect_warning(
    GT_with_invalid <- getGT(
      gdb,
      VAR_id = c(1:100, 6000, 6001),
      cohort = "pheno",
      verbose = FALSE
    ),
    "Retrieved genotypes for only 100 of 102 variants"
  )

  # should return same result as valid VAR_ids only
  expect_identical(GT_reference, GT_with_invalid)

  # provide VAR_ids that are all absent in gdb
  expect_warning({
    GT <- getGT(
      gdb,
      VAR_id = 5e6 + 1:10
    )
  })
  ## output should be an empty genoMatrix
  expect_identical(nrow(GT), 0L)
  sm <- getCohort(gdb, "SM")
  expect_identical(ncol(GT), nrow(sm))
})

test_that("getGT works with varSet and weights", {
  gdb <- gdb(rvatData::rvat_example("rvatData.gdb"))

  # build a varset containing CADD weights
  CADD_file <- withr::local_tempfile(fileext = ".txt.gz")
  varset <- buildVarSet(
    object = gdb,
    output = CADD_file,
    varSetName = "CADD",
    unitTable = "varInfo",
    unitName = "gene_name",
    weightName = "CADDphred",
    where = "gene_name in ('CYP19A1', 'FUS', 'OPTN')",
    verbose = FALSE
  )
  varset <- varSetFile(CADD_file)

  # extract genotypes with varSet
  GT_varset <- suppressWarnings(getGT(
    gdb,
    varSet = getVarSet(varset, unit = "FUS"),
    cohort = "pheno",
    verbose = FALSE
  ))

  # compare with manually extracting and adding CADD scores
  varinfo <- getAnno(
    gdb,
    "varinfo",
    VAR_id = listVars(getVarSet(varset, unit = "FUS")[[1]])
  )
  GT_manual <- getGT(
    gdb,
    cohort = "pheno",
    VAR_id = listVars(getVarSet(varset, unit = "FUS")[[1]]),
    varSetName = "CADD",
    unit = "FUS",
    verbose = FALSE
  )
  GT_manual <- suppressWarnings(recode(
    GT_manual,
    weights = as.numeric(varinfo$CADDphred)
  ))

  expect_equal(GT_varset, GT_manual)
})

test_that("getGT handles different cohort specifications", {
  gdb <- gdb(rvatData::rvat_example("rvatData.gdb"))

  # test specifying gdb cohort name
  GT1 <- getGT(gdb, VAR_id = 1:100, cohort = "pheno", verbose = FALSE)

  # test cohort as data.frame
  chrt <- getCohort(gdb, "pheno")
  GT2 <- getGT(gdb, VAR_id = 1:100, cohort = chrt, verbose = FALSE)

  # test cohort as DataFrame with metadata
  chrt_df <- S4Vectors::DataFrame(chrt)
  metadata(chrt_df)$name <- "pheno"
  GT3 <- getGT(gdb, VAR_id = 1:100, cohort = chrt_df, verbose = FALSE)

  metadata(GT2)$cohort <- "pheno"
  expect_equal(GT1, GT2)
  expect_equal(GT1, GT3)
})

test_that("getGT annotation functionality works", {
  gdb <- gdb(rvatData::rvat_example("rvatData.gdb"))

  # GT without annotations
  GT1 <- getGT(gdb, VAR_id = 1:100, cohort = "pheno", verbose = FALSE)

  # GT with var table anno
  GT2 <- getGT(
    gdb,
    VAR_id = 1:100,
    cohort = "pheno",
    anno = "var",
    verbose = FALSE
  )
  var <- getAnno(gdb, "var", VAR_id = 1)
  expect_gt(ncol(rowData(GT2)), ncol(rowData(GT1)))
  expect_true(all(setdiff(colnames(var), "VAR_id") %in% colnames(rowData(GT2))))

  # GT with specific annotation fields
  GT3 <- getGT(
    gdb,
    VAR_id = 1:100,
    cohort = "pheno",
    anno = "var",
    annoFields = c("VAR_id", "CHROM", "POS", "ID", "REF", "ALT"),
    verbose = FALSE
  )
  expect_gt(
    ncol(SummarizedExperiment::rowData(GT2)),
    ncol(SummarizedExperiment::rowData(GT3))
  )
  expect_equal(ncol(SummarizedExperiment::rowData(GT3)), 7L)
  expect_true(all(
    c("CHROM", "POS", "ID", "REF", "ALT") %in% colnames(rowData(GT2))
  ))
  expect_true(
    !all(setdiff(colnames(var), "VAR_id") %in% colnames(rowData(GT3)))
  )

  # includeVarInfo should produce same output as anno = "var"
  GT_includevarinfo <- getGT(
    gdb,
    VAR_id = 1:10,
    cohort = "pheno",
    includeVarInfo = TRUE,
    verbose = FALSE
  )
  GT_includevar <- getGT(
    gdb,
    VAR_id = 1:10,
    cohort = "pheno",
    anno = "var",
    verbose = FALSE
  )
  expect_equal(GT_includevarinfo, GT_includevar)
})

test_that("getGT ranges works", {
  gdb <- create_example_gdb()

  # check if ranges are extracted correctly
  var <- getAnno(
    gdb,
    "var",
    VAR_id = 100:200,
    fields = c("VAR_id", "CHROM", "POS")
  )
  GT_var <- getGT(gdb, VAR_id = 100:200, verbose = FALSE)
  GT_ranges <- getGT(
    gdb,
    ranges = data.frame(
      CHROM = var$CHROM,
      start = var$POS,
      end = var$POS
    ),
    verbose = FALSE
  )
  expect_equal(GT_var, GT_ranges)
})

test_that("ploidy is set correctly", {
  gdb <- create_example_gdb()

  # get VAR_ids on chrX
  vars <- getAnno(
    gdb,
    "var",
    fields = "VAR_id",
    where = "CHROM = 'chrX'"
  )$VAR_id
  GT_chrX <- getGT(gdb, cohort = "pheno", VAR_id = vars, verbose = FALSE)
  # should contain both 'diploid' and 'XnonPAR'
  expect_identical(
    sort(unique(rowData(GT_chrX)$ploidy)),
    sort(c("diploid", "XnonPAR"))
  )
  # none of the male genotypes should be > 1 for non-PAR chrX variants
  count_XnonPAR <- sum(
    assays(GT_chrX[
      rowData(GT_chrX)$ploidy == "XnonPAR",
      GT_chrX$sex == 1
    ])$GT >
      1,
    na.rm = TRUE
  )
  expect_identical(count_XnonPAR, 0L)
  # counts of >1 do occur for diploid variants on chrX
  count_diploid <-
    sum(
      assays(GT_chrX[
        rowData(GT_chrX)$ploidy == "diploid",
        GT_chrX$sex == 1
      ])$GT >
        1,
      na.rm = TRUE
    )
  expect_gt(count_diploid, 0L)

  # when genomeBuild is missing, all variants should be treated as diploid
  meta <- DBI::dbGetQuery(gdb, "SELECT * FROM meta")
  meta[meta$name == "genomeBuild", "value"] <- NA_character_
  DBI::dbWriteTable(gdb, "meta", meta, overwrite = TRUE)
  GT_chrX <- getGT(gdb, cohort = "pheno", VAR_id = vars, verbose = FALSE)
  # should contain both 'diploid' and 'XnonPAR'
  expect_identical(
    sort(unique(rowData(GT_chrX)$ploidy)),
    "diploid"
  )
})

test_that("getGT validation works", {
  gdb <- create_example_gdb()

  # expect error when VAR_id is not numeric
  expect_error(
    getGT(
      gdb,
      VAR_id = c(TRUE, FALSE),
      cohort = "pheno",
      verbose = FALSE
    ),
    regexp = "`VAR_id` should be"
  )

  # expect error when VAR_id contains NAs
  expect_error(
    getGT(
      gdb,
      VAR_id = c(1, NA, 3),
      cohort = "pheno",
      verbose = FALSE
    ),
    regexp = "`VAR_id` contains missing values"
  )

  # expect error when VAR_id contains duplicates
  expect_error(
    getGT(
      gdb,
      VAR_id = c(1, 1, 3),
      cohort = "pheno",
      verbose = FALSE
    ),
    regexp = "`VAR_id` contains duplicated values"
  )
  
  # expect error
  

  # expect error when annotation has multiple rows per VAR_id

  ## upload annotation with multiple rows per VAR_id
  uploadAnno(
    gdb,
    name = "transcriptInfo",
    value = data.frame(
      VAR_id = c(1, 1, 2, 2, 2, 3),
      transcript = c("A", "B", "C", "D", "E", "F")
    ),
    skipRemap = TRUE,
    verbose = FALSE
  )

  expect_error(
    getGT(
      gdb,
      VAR_id = 1:100,
      cohort = "pheno",
      anno = "transcriptInfo",
      verbose = FALSE
    ),
    regexp = "should contain 1 row per VAR_id"
  )

  # expect error when varSetFile is supplied
  varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
  expect_error(
    getGT(
      gdb,
      varSet = varsetfile,
      cohort = "pheno",
      verbose = FALSE
    ),
    regexp = "should be of class"
  )

  # expect error when varSetList has multiple rows
  varsetlist <- getVarSet(varsetfile, unit = c("CYP19A1", "FUS"))
  expect_error(
    getGT(
      gdb,
      varSet = varsetlist,
      cohort = "pheno",
      verbose = FALSE
    ),
    regexp = "A varSetList with >1 varSets is supplied."
  )

  # expect error when unrecognized value of checkPloidy is provided
  expect_error(
    getGT(
      gdb,
      VAR_id = 1:10,
      cohort = "pheno",
      verbose = FALSE,
      checkPloidy = "GRCh40"
    ),
    regexp = "Unrecognized value for checkPloidy"
  )

  # expect error when none of VAR_id, varSet of ranges is provided
  expect_error(
    getGT(
      gdb,
      cohort = "pheno",
      verbose = FALSE
    ),
    regexp = "Specify one of"
  )

  # expect error when invalid cohort is provided
  pheno <- getCohort(gdb, "pheno")
  expect_error(
    getGT(
      gdb,
      VAR_id = 1:10,
      cohort = as.matrix(pheno),
      verbose = FALSE
    ),
    regexp = "should be either"
  )

  # expect message when non-PAR chrX variants are loaded
  vars <- getAnno(
    gdb,
    "var",
    fields = "VAR_id",
    where = "CHROM = 'chrX'"
  )$VAR_id
  expect_message(
    expect_message(
      test <- getGT(gdb, VAR_id = vars),
      regexp = "Retrieved genotypes"
    ),
    regexp = "Ploidy of.*non-pseudoautosomal"
  )

  # should error when varSet has multiple units
  CADD_file <- withr::local_tempfile(fileext = ".txt.gz")
  varset <- buildVarSet(
    object = gdb,
    output = CADD_file,
    varSetName = "CADD",
    unitTable = "varInfo",
    unitName = "gene_name",
    weightName = "CADDphred",
    where = "gene_name in ('CYP19A1', 'FUS', 'OPTN')",
    verbose = FALSE
  )
  varset <- varSetFile(CADD_file)
  expect_error(
    getGT(
      gdb,
      VAR_id = 1:100,
      cohort = "pheno",
      varSet = getVarSet(varset, listUnits(varset))
    ),
    regexp = "Specify only one of"
  )
})
