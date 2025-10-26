data(GT)
gdb <- create_example_gdb()
varids <- getAnno(gdb, "var", fields = "VAR_id")$VAR_id
GT_all <- getGT(
  gdb,
  VAR_id = varids,
  cohort = "pheno",
  verbose = FALSE
)

test_that("updateGT annotation updating work", {
  # test annotation data
  rowdata <- data.frame(
    VAR_id = rownames(GT),
    A = sample(LETTERS, size = nrow(GT), replace = TRUE),
    B = sample(LETTERS, size = nrow(GT), replace = TRUE),
    stringsAsFactors = FALSE
  )

  # check if the input order doesn't matter
  GT1 <- GT
  GT1 <- updateGT(GT1, anno = rowdata)
  GT2 <- GT
  GT2 <- updateGT(GT2, anno = rowdata[sample(1:nrow(GT)), ])
  expect_identical(GT1, GT2)
})


test_that("updateGT sample metadata updating works", {
  # add additional column to sample metadata
  GT_addcolumn <- GT
  sm <- colData(GT_addcolumn)
  sm$random_column <- sample(LETTERS, size = ncol(GT_addcolumn), replace = TRUE)
  GT_addcolumn <- updateGT(GT_addcolumn, SM = sm)
  expect_identical(sm, colData(GT_addcolumn))

  # check if order matters (it shouldn't)
  GT_addcolumn_randomorder <- GT
  sm <- colData(GT_addcolumn_randomorder)
  sm$random_column <- sample(
    LETTERS,
    size = ncol(GT_addcolumn_randomorder),
    replace = TRUE
  )
  GT_addcolumn_randomorder <- updateGT(
    GT_addcolumn_randomorder,
    SM = sm[sample(1:nrow(sm), size = nrow(sm)), ]
  )
  expect_identical(sm, colData(GT_addcolumn_randomorder))
})

test_that("updateGT input validation works", {
  # expect error if input anno contains duplicate VAR_ids
  anno_dup <- data.frame(
    VAR_id = c(rownames(GT)[1], rownames(GT)[1]),
    A = c("X", "Y"),
    stringsAsFactors = FALSE
  )
  expect_error(
    {
      updateGT(GT, anno = anno_dup)
    },
    regexp = "shouldn't contain duplicated VAR_ids"
  )

  # expect warning when annotations are added to empty genoMatrix
  expect_warning(
    {
      updateGT(GT[getMAF(GT) > 0.5, ], anno = anno_dup[1, ])
    },
    regexp = "Annotations can't be added"
  )

  # expect error when input anno misses VAR_id column
  expect_error(
    {
      updateGT(GT, anno = anno_dup[1, ] %>% dplyr::select(-VAR_id))
    },
    regexp = "anno should contain a `VAR_id`"
  )

  # expect error when anno contains protected rowData columns
  anno <- anno_dup[1, ]
  for (protected_col in c("ploidy", "w", "AF")) {
    expect_error(
      {
        updateGT(GT, anno = anno %>% dplyr::mutate(!!protected_col := "test"))
      },
      regexp = "are protected"
    )
  }

  # expect warning when anno contains subset of variants
  expect_warning(
    {
      updateGT(GT, anno = anno)
    },
    regexp = "Not all variants present in the genoMatrix"
  )

  # expect error if input SM contains duplicate IIDs
  sm_dupIIDs <- colData(GT)
  sm_dupIIDs <- rbind(
    sm_dupIIDs[, c("IID", "sex")],
    DataFrame(
      IID = c("sample1", "sample2"),
      sex = c(0, 0),
      row.names = c("sample1", "sample2")
    )
  )
  expect_error({
    updateGT(GT, SM = sm_dupIIDs)
  })

  # expect error if input SM is missing IIDs
  sm_subset <- colData(GT)[2000:3000, ]
  expect_error(
    {
      updateGT(GT, SM = sm_subset)
    },
    regexp = "Number of samples"
  )

  # expect error when SM has non-existent IIDs
  sm <- colData(GT)
  sm$IID[1:5] <- paste0("sample", 1:5)
  expect_error(
    {
      updateGT(GT, SM = sm)
    },
    regexp = "IID values in SM table do not match"
  )

  # expect error when SM has invalid type
  sm <- colData(GT)
  expect_error(
    {
      updateGT(GT, SM = as.matrix(sm))
    },
    regexp = "must be a data.frame or DFrame"
  )

  # expect error when genoMatrix includes non-autosomal variants
  # and updateGT is called with SM
  expect_error(
    {
      updateGT(GT_all, SM = colData(GT_all))
    },
    regexp = "Genotype dosage values at non-diploid"
  )
})

test_that("flipToMinor works", {
  # load test data
  gdb <- create_example_gdb()
  var <- getAnno(gdb, "varInfo", where = "gene_name = 'ABCA4'")
  GT <- getGT(
    gdb,
    VAR_id = var$VAR_id,
    anno = "var",
    verbose = FALSE
  )

  # before flipping, MAC should equal AC for all variants
  expect_equal(getMAC(GT), getAC(GT))
  GT <- flipToMinor(GT)

  # manually flip some alleles
  flip <- rep(FALSE, nrow(GT))
  flip[c(47, 52, 71)] <- rep(TRUE, 3)
  assays(GT)$GT <- abs(
    assays(GT)$GT -
      2 *
        matrix(
          rep(flip, each = metadata(GT)$m),
          nrow = metadata(GT)$nvar,
          byrow = TRUE
        )
  )
  rowData(GT)$AF <- getAF(GT)

  # now MAC and AC should differ for flipped variants
  expect_failure(expect_equal(getMAC(GT), getAC(GT)))

  # and flip to minor again
  GT_flipped <- flipToMinor(GT)
  flip <- getAF(GT) > 0.5
  expect_equal((1 - getAF(GT_flipped)[flip]), getAF(GT)[flip])
  expect_equal((rowData(GT_flipped)$effectAllele != rowData(GT)$ALT), flip)
  expect_equal((rowData(GT_flipped)$otherAllele != rowData(GT)$REF), flip)
})

test_that("flipToMinor works with X-chromosomal variants", {
  # test data
  gdb <- create_example_gdb()
  var <- getAnno(gdb, "varInfo", where = "gene_name = 'UBQLN2'")
  GT <- getGT(
    gdb,
    VAR_id = var$VAR_id,
    anno = "var",
    verbose = FALSE
  )

  # before flipping, MAC should equal AC for all variants
  expect_equal(getMAC(GT), getAC(GT))
  GT <- flipToMinor(GT)

  # manually flip some X-chromosome alleles
  carriers <- getCarriers(GT, colDataFields = "sex")

  # select variants that have carriers in both sexes for testing
  vars <- carriers %>%
    dplyr::group_by(VAR_id) %>%
    dplyr::filter(all(c(1, 2) %in% sex)) %>%
    dplyr::pull(VAR_id) %>%
    unique() %>%
    head(5)

  # get carriers for selected variants
  vars_carriers <- carriers %>%
    dplyr::filter(VAR_id %in% vars)

  # manually flip alleles
  for (varid in unique(vars_carriers$VAR_id)) {
    vars_male <- vars_carriers[
      vars_carriers$VAR_id == varid & vars_carriers$sex == 1,
    ]
    vars_female <- vars_carriers[
      vars_carriers$VAR_id == varid & vars_carriers$sex == 2,
    ]

    # flip male genotypes
    assays(GT)$GT[varid, vars_male$IID] <- 0
    assays(GT)$GT[varid, !colnames(GT) %in% vars_male$IID] <- 1

    # flip female genotypes
    assays(GT)$GT[varid, vars_female$IID] <- abs(
      assays(GT)$GT[varid, vars_female$IID] - 2
    )
  }

  # recalculate AFs
  rowData(GT)$AF <- getAF(GT)

  # MAC != AC
  expect_failure(expect_equal(getMAC(GT), getAC(GT)))

  # flip again
  GT_flipped <- flipToMinor(GT)
  flip <- getAF(GT) > 0.5

  # verify
  expect_equal((1 - getAF(GT_flipped)[flip]), getAF(GT)[flip])
  expect_equal((rowData(GT_flipped)$effectAllele != rowData(GT)$ALT), flip)
  expect_equal((rowData(GT_flipped)$otherAllele != rowData(GT)$REF), flip)
})

test_that("flipToMinor input validation works", {
  # expect warning when geneticModel is not allelic
  expect_warning(
    {
      flipToMinor(recode(GT, geneticModel = "dominant"))
    },
    regexp = "flipToMinor only applies when geneticModel"
  )
})
