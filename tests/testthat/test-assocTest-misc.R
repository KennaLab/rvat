data(GT)
data(GTsmall)

test_that("assocTest covar/pheno missingness is handled correctly", {
  # compare missingness in covariates with keep-list
  GT_missing_covar <- GTsmall
  
  ## introduce some missing values in a covariate
  GT_missing_covar$covar_with_missings <- GT_missing_covar$PC1
  GT_missing_covar_check <- GT_missing_covar
  GT_missing_covar$covar_with_missings[c(1:10, 3500:3600)] <- NA
  keeplist <- colnames(GT_missing_covar)[
    !is.na(GT_missing_covar$covar_with_missings)
  ]
  ## run assocTest
  assoc_missing_covar <- assocTest(
    object = GT_missing_covar,
    pheno = "pheno",
    test = "scoreSPA",
    covar = c("PC2", "covar_with_missings"),
    verbose = TRUE
  )
  
  ## run assocTest with keep-list based on non-missing covariate values
  assoc_missing_covar_check <- suppressMessages(assocTest(
    object = GT_missing_covar_check,
    pheno = "pheno",
    test = "scoreSPA",
    covar = c("PC2", "covar_with_missings"),
    keep = keeplist,
    verbose = TRUE
  ))
  metadata(assoc_missing_covar)$creationDate <- NA_character_
  metadata(assoc_missing_covar_check)$creationDate <- NA_character_
  
  ## compare
  expect_equal(assoc_missing_covar, assoc_missing_covar_check)

  # compare missingness in pheno with keep-list
  GT_missing_pheno <- GTsmall

  ## introduce some missing values in a phenotype
  GT_missing_pheno$pheno_with_missings <- GT_missing_pheno$PC1
  GT_missing_pheno_check <- GT_missing_pheno
  GT_missing_pheno$pheno_with_missings[c(1:10, 3500:3600)] <- NA
  keeplist <- colnames(GT_missing_pheno)[
    !is.na(GT_missing_pheno$pheno_with_missings)
  ]

  assoc_missing_pheno <- suppressMessages(assocTest(
    object = GT_missing_pheno,
    pheno = "pheno_with_missings",
    test = "lm",
    covar = paste0("PC", 2:4),
    continuous = TRUE,
    verbose = TRUE
  ))

  ## run assocTest with keep-list based on non-missing pheno values
  assoc_missing_pheno_check <- suppressMessages(assocTest(
    object = GT_missing_pheno_check,
    pheno = "pheno_with_missings",
    test = "lm",
    covar = paste0("PC", 2:4),
    continuous = TRUE,
    keep = keeplist,
    verbose = TRUE
  ))
  metadata(assoc_missing_pheno)$creationDate <- NA_character_
  metadata(assoc_missing_pheno_check)$creationDate <- NA_character_
  expect_equal(assoc_missing_pheno, assoc_missing_pheno_check)
})

test_that("effectAllele is assigned correctly in sv tests", {
  GT_check <- GT

  # manually flip some alleles to test effectAllele assignment
  flip <- rep(FALSE, nrow(GT_check))
  names(flip) <- rownames(GT_check)
  flip[c(1, 7, 21)] <- rep(TRUE, 3)
  assays(GT_check)$GT <- abs(
    assays(GT_check)$GT -
      2 *
        matrix(
          rep(flip, each = metadata(GT_check)$m),
          nrow = metadata(GT_check)$nvar,
          byrow = TRUE
        )
  )
  rowData(GT_check)$AF <- getAF(GT_check)
  sv <- assocTest(
    GT_check,
    covar = paste0("PC", 1:4),
    pheno = "pheno",
    singlevar = TRUE,
    test = c("firth", "scoreSPA", "glm"),
    verbose = FALSE
  )
  var <- as.data.frame(rowData(GT_check))
  var$VAR_id <- as.character(rownames(var))
  var <- var %>%
    dplyr::left_join(as.data.frame(sv), by = "VAR_id")
  var_check <- var %>% dplyr::filter(effectAllele != ALT)
  expect_identical(
    sort(unique(var_check$VAR_id)),
    as.character(sort(names(which(flip))))
  )
})