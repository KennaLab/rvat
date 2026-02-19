data(GT)
gdb <- create_example_gdb()
varids <- getAnno(gdb, "var", fields = "VAR_id")$VAR_id
GT_all <- getGT(
  gdb,
  VAR_id = varids,
  cohort = "pheno",
  verbose = FALSE
)

test_that("genoMatrix getters work", {
  gdb <- create_example_gdb()

  # load freq/counts/hwe/cr generated with plink
  plink_freq <- readr::read_table(
    test_path("data/rvatData.frq.gz"),
    show_col_types = FALSE
  )
  plink_counts <- readr::read_table(
    test_path("data/rvatData.frq.counts.gz"),
    show_col_types = FALSE
  )
  plink_hwe <- suppressWarnings(readr::read_table(
    test_path("data/rvatData.hwe.gz"),
    show_col_types = FALSE
  ))
  plink_cr <- readr::read_table(
    test_path("data/rvatData.lmiss.gz"),
    show_col_types = FALSE
  )

  ## format
  plink <- plink_freq %>%
    dplyr::select(SNP, MAF) %>%
    dplyr::left_join(
      plink_counts %>% dplyr::select(SNP, C1, C2),
      by = "SNP"
    ) %>%
    dplyr::left_join(plink_hwe %>% dplyr::select(SNP, P), by = "SNP") %>%
    dplyr::left_join(plink_cr %>% dplyr::select(SNP, F_MISS), by = "SNP") %>%
    dplyr::mutate(
      SNP = stringr::str_replace(SNP, "^25:", "X:"),
      SNP = stringr::str_replace(SNP, "^23:", "X:")
    )

  ## generate stats with RVAT
  var <- getAnno(gdb, "var", fields = c("VAR_id", "CHROM", "POS", "REF", "ALT"))
  var$VAR_id <- as.character(var$VAR_id)
  var$ID <- sprintf(
    "%s:%s:%s:%s",
    stringr::str_replace(var$CHROM, "chr", ""),
    var$POS,
    var$REF,
    var$ALT
  )
  sumgeno <- summariseGeno(GT_all)
  var <- var[match(rownames(GT_all), as.character(var$VAR_id)), ]
  var$MAF <- getMAF(GT_all)
  var$AF <- getAF(GT_all)
  var$AC <- getAC(GT_all)
  var$CR <- getCR(GT_all)
  var$ncarriers <- getNCarriers(GT_all)

  ## combine RVAT and plink results
  var <- var %>%
    dplyr::left_join(
      sumgeno[, c("VAR_id", "geno0", "geno1", "geno2", "hweP")],
      by = "VAR_id"
    ) %>%
    dplyr::left_join(plink, by = c("ID" = "SNP"))

  # compare results
  expect_equal(var$MAF.x, var$MAF.y, tolerance = 1e-3)
  expect_equal(var$AC, var$C1)
  expect_equal(var$CR, (1 - var$F_MISS), tolerance = 1e-3)
  expect_equal(sumgeno$AF, var$AF)

  var <- var %>%
    dplyr::mutate(check_carriers = AC - geno2)
  expect_equal(var$ncarriers, var$check_carriers)

  # check sample call rates
  plink_cr <- readr::read_table(
    test_path("data/rvatData.imiss.gz"),
    show_col_types = FALSE
  )
  rvat_cr <- getCR(GT_all, var = FALSE)
  expect_equal(unname(rvat_cr), (1 - plink_cr$F_MISS), tolerance = 1e-3)
})


test_that("genoMatrix getters work for random subset", {
  # compare for a random subset

  # load freq/counts/hwe/cr generated with plink
  plink_freq <- readr::read_table(
    test_path("data/rvatData.randomsubset.frq.gz"),
    show_col_types = FALSE
  )
  plink_counts <- readr::read_table(
    test_path("data/rvatData.randomsubset.frq.counts.gz"),
    show_col_types = FALSE
  )
  plink_hwe <- suppressWarnings(readr::read_table(
    test_path("data/rvatData.randomsubset.hwe.gz"),
    show_col_types = FALSE
  ))
  plink_cr <- readr::read_table(
    test_path("data/rvatData.randomsubset.lmiss"),
    show_col_types = FALSE
  )

  ## format
  plink <- plink_freq %>%
    dplyr::select(SNP, MAF) %>%
    dplyr::left_join(
      plink_counts %>% dplyr::select(SNP, C1, C2),
      by = "SNP"
    ) %>%
    dplyr::left_join(plink_hwe %>% dplyr::select(SNP, P), by = "SNP") %>%
    dplyr::left_join(plink_cr %>% dplyr::select(SNP, F_MISS), by = "SNP") %>%
    dplyr::mutate(
      SNP = stringr::str_replace(SNP, "^25:", "X:"),
      SNP = stringr::str_replace(SNP, "^23:", "X:")
    )

  ## generate stats with RVAT
  var <- getAnno(gdb, "var", fields = c("VAR_id", "CHROM", "POS", "REF", "ALT"))
  var$VAR_id <- as.character(var$VAR_id)
  var$ID <- sprintf(
    "%s:%s:%s:%s",
    stringr::str_replace(var$CHROM, "chr", ""),
    var$POS,
    var$REF,
    var$ALT
  )
  keep <- readr::read_table(
    test_path("data/rvatData.random.keep"),
    col_names = FALSE,
    show_col_types = FALSE
  )
  GT_subset <- GT_all[, keep$X1]
  sumgeno <- summariseGeno(GT_subset)
  var <- var[match(rownames(GT_subset), as.character(var$VAR_id)), ]
  var$MAF <- getMAF(GT_subset)
  var$AF <- getAF(GT_subset)
  var$AC <- getAC(GT_subset)
  var$CR <- getCR(GT_subset)
  var$ncarriers <- getNCarriers(GT_subset)

  ## combine RVAT and plink results
  var <- var %>%
    dplyr::left_join(
      sumgeno[, c("VAR_id", "geno0", "geno1", "geno2", "hweP")],
      by = "VAR_id"
    ) %>%
    dplyr::left_join(plink, by = c("ID" = "SNP"))

  # compare results
  expect_equal(var$MAF.x, var$MAF.y, tolerance = 1e-3)
  expect_equal(var$AC, var$C1)
  expect_equal(var$CR, (1 - var$F_MISS), tolerance = 1e-3)
  expect_equal(sumgeno$AF, var$AF)

  var <- var %>%
    dplyr::mutate(check_carriers = AC - geno2)
  expect_equal(var$ncarriers, var$check_carriers)

  # Also check sample callrates for subset
  plink_cr <- readr::read_table(
    test_path("data/rvatData.randomsubset.imiss.gz"),
    show_col_types = FALSE
  )
  plink_cr <- plink_cr[match(colnames(GT_subset), plink_cr$IID), ]
  rvat_cr <- getCR(GT_subset, var = FALSE)
  expect_equal(unname(rvat_cr), (1 - plink_cr$F_MISS), tolerance = 1e-3)
})


test_that("getAF works with sex chromosome variants", {
  # test YnonPAR variants (currently not in rvatData, might add later)
  # YnonPAR variants should only be present in males, so AF calculation should
  # be based on male samples only, ignoring females entirely
  # create test case by converting XnonPAR variants to YnonPAR
  GT_y <- GT_all
  rowData(GT_y)$ploidy <- dplyr::case_when(
    rowData(GT_y)$ploidy == "XnonPAR" ~ "YnonPAR",
    TRUE ~ rowData(GT_y)$ploidy
  )
  metadata(GT_y)$ploidyLevels <- unique(rowData(GT_y)$ploidy)
  GT_y <- rvat:::.resetSexChromDosage(GT_y)

  # calculate allele frequencies
  AF <- getAF(GT_y)
  AF_males <- getAF(GT_y[, GT_y$sex == 1])

  # for YnonPAR variants, overall AF should equal
  # AF calculated from males only.
  expect_equal(
    AF[rowData(GT_y)$ploidy == "YnonPAR"],
    AF_males[rowData(GT_y)$ploidy == "YnonPAR"]
  )
  expect_equal(
    getAF(GT_all[, GT_all$sex == 1])[rowData(GT_all)$ploidy == "XnonPAR"],
    AF_males[rowData(GT_y)$ploidy == "YnonPAR"]
  )

  # AF calculated on YnonPAR-only subset should match full matrix
  AF_Y_only <- getAF(GT_y[rowData(GT_y)$ploidy == "YnonPAR", ])
  AF_Y_only_males <- getAF(GT_y[
    rowData(GT_y)$ploidy == "YnonPAR",
    GT_y$sex == 1
  ])

  # subsetting should not change AF calculations
  expect_equal(
    AF[rowData(GT_y)$ploidy == "YnonPAR"],
    AF_Y_only
  )

  # from full matrix should equal AF from males-only subset
  expect_equal(
    AF[rowData(GT_y)$ploidy == "YnonPAR"],
    AF_Y_only_males
  )

  # test consistency with summariseGeno
  sumgeno_y <- summariseGeno(GT_y)
  expect_equal(
    unname(AF[rowData(GT_y)$ploidy == "YnonPAR"]),
    sumgeno_y[rowData(GT_y)$ploidy == "YnonPAR", ]$AF
  )
  expect_equal(
    sumgeno_y[rowData(GT_y)$ploidy == "YnonPAR", ]$hweP,
    rep(1.0, sum(rowData(GT_y)$ploidy == "YnonPAR"))
  )

  # expect NAs when sex is missing
  GT_y$sex <- 0
  expect_warning(
    {
      AF <- getAF(GT_y)
    },
    regexp = "Note that sex is missing"
  )
  expect_true(all(is.na(AF[rowData(GT_y)$ploidy == "YnonPAR"])))

  # test mixed X/Y chromosome scenaratios
  GT_xy <- GT_all
  rowData(GT_xy)$ploidy[rowData(GT_xy)$ploidy == "XnonPAR"][1:15] <- "YnonPAR"
  metadata(GT_xy)$ploidyLevels <- unique(rowData(GT_xy)$ploidy)
  GT_xy <- rvat:::.resetSexChromDosage(GT_xy)

  ## calculate sex-stratified allele frequencies
  AF <- getAF(GT_xy)
  AF_males <- getAF(GT_xy[, GT_xy$sex == 1])
  AF_females <- getAF(GT_xy[, GT_xy$sex == 2])

  ## for YnonPAR variants, overall AF should equal male-only AF
  expect_equal(
    AF[rowData(GT_xy)$ploidy == "YnonPAR"],
    AF_males[rowData(GT_xy)$ploidy == "YnonPAR"]
  )

  # verify that the YnonPAR variants have same AF
  # as the original XnonPAR variants in males
  expect_equal(
    getAF(GT_all[, GT_all$sex == 1])[rowData(GT_all)$ploidy == "XnonPAR"][1:15],
    AF_males[rowData(GT_xy)$ploidy == "YnonPAR"]
  )

  # for remaining XnonPAR variants, AF should match female AF from original data
  expect_equal(
    getAF(GT_all[, GT_all$sex == 2])[rowData(GT_all)$ploidy == "XnonPAR"][
      16:sum(rowData(GT_all)$ploidy == "XnonPAR")
    ],
    AF_females[rowData(GT_xy)$ploidy == "XnonPAR"]
  )

  ## return NA when sex is missing
  GT_xy$sex <- 0
  expect_warning(
    {
      AF <- getAF(GT_xy)
    },
    regexp = "Note that sex is missing"
  )
  expect_true(all(is.na(AF[
    rowData(GT_xy)$ploidy %in% c("YnonPAR", "XnonPAR")
  ])))
})

test_that("genoMatrix getters work on recoded genoMatrices", {
  ac <- getAC(GT_all)
  ac_imputed <- getAC(recode(GT_all, imputeMethod = "meanImpute"))
  expect_equal(ac, ac_imputed)
})


test_that("genoMatrix getters input validation works", {
  expect_warning(
    {
      af_recessive <- getAF(recode(GT_all, geneticModel = "recessive"))
    },
    regexp = "are set to NA"
  )
  expect_true(length(af_recessive) == nrow(GT_all) && all(is.na(af_recessive)))

  # ploidy == XnonPAR and sex contains unknowns
  GT_check <- GT_all
  colData(GT_check)$sex[sample(1:ncol(GT_check), size = 100)] <- 0
  expect_warning(
    {
      af <- getAF(GT_check)
    },
    regexp = "contains variants with ploidy"
  )

  # expect warning when getAC is called on a non-allelic genoMatrix
  GT_all_recessive <- recode(GT_all, geneticModel = "recessive")
  expect_warning(
    {
      ac <- getAC(GT_all_recessive)
    },
    regexp = "Allele counts can't be calculated"
  )
})

test_that("genoMatrix-summarisegeno input validation works", {
  # expect error when summariseGeno is called on inputed genoMatrix
  GT_imputed <- recode(GT_all, imputeMethod = "meanImpute")
  expect_error(
    {
      summariseGeno(GT_imputed)
    },
    regexp = "provide a non-imputed genoMatrix"
  )

  # expact warning when XnonPAR variants are included and sex is missing
  GT_missing_sex <- GT_all[rowData(GT_all)$ploidy == "XnonPAR", ]
  GT_missing_sex$sex[1:5] <- 0
  ## this will actually return two warnings (AF and hweP)
  expect_warning(
    {
      expect_warning(
        {
          summariseGeno(GT_missing_sex)
        },
        regexp = "contains variants with ploidy='XnonPAR'"
      )
    },
    regexp = "contains variants with ploidy='XnonPAR'"
  )
})
