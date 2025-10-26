data(GT)

test_that("recoding and aggregate work", {
  # test varying parameters for recode and aggregate using snapshots
  expect_snapshot_value(
    aggregate(recode(GT, imputeMethod = "meanImpute"), returnGT = FALSE),
    style = "serialize"
  )
  expect_snapshot_value(
    aggregate(recode(GT, imputeMethod = "missingToRef"), returnGT = FALSE),
    style = "serialize"
  )
  expect_snapshot_value(
    aggregate(
      recode(GT, imputeMethod = "meanImpute", geneticModel = "dominant"),
      returnGT = FALSE
    ),
    style = "serialize"
  )
  expect_snapshot_value(
    aggregate(
      recode(GT, imputeMethod = "meanImpute", geneticModel = "recessive"),
      returnGT = FALSE
    ),
    style = "serialize"
  )
  expect_snapshot_value(
    aggregate(
      recode(GT, imputeMethod = "meanImpute", MAFweights = "mb"),
      returnGT = FALSE
    ),
    style = "serialize"
  )
  expect_snapshot_value(
    aggregate(
      recode(
        GT,
        imputeMethod = "meanImpute",
        MAFweights = "mb",
        geneticModel = "dominant"
      ),
      returnGT = FALSE
    ),
    style = "serialize"
  )
  expect_snapshot_value(
    aggregate(
      recode(
        GT,
        imputeMethod = "meanImpute",
        MAFweights = "mb",
        geneticModel = "recessive"
      ),
      returnGT = FALSE
    ),
    style = "serialize"
  )
})

test_that("genoMatrix recode input validation works", {
  # expect error when MAFweights is invalid
  expect_error(
    recode(GT, MAFweights = "invalid"),
    regexp = "should be either"
  )
  expect_error(
    .calc_maf_weights(c(1, 2, 3), af = NULL, method = "invalid"),
    regexp = "should be either"
  )

  # expect error when calc_maf_weights af/w don't match
  expect_error(
    .calc_maf_weights(c(1, 2, 3), af = c(0.01, 0.02), method = "mb"),
    regexp = "Unequal lengths"
  )

  # expect error when calc_maf_weights
  expect_error(
    .calc_maf_weights(c(1, 2, 3), af = NULL, method = "mb"),
    regexp = "AF should be specified for 'mb' weighting"
  )

  # expect error when geneticModel is invalid
  expect_error(
    recode(GT, geneticModel = "invalid"),
    regexp = "does not represent"
  )

  # expect error when geneticModel is invalid
  expect_error(
    recode(GT, imputeMethod = "invalid"),
    regexp = "not a recognized"
  )

  # expect error when weights has incorrect length
  expect_error(
    recode(GT, weight = c(1, 2)),
    regexp = "Length of `weights` should equal"
  )

  # expect error when geneticModel is specified
  # and current geneticModel != allelic
  GT_recessive <- recode(GT, geneticModel = "recessive")
  GT_dominant <- recode(GT, geneticModel = "dominant")
  expect_error(
    recode(GT_recessive, geneticModel = "dominant"),
    regexp = "should be 'allelic'"
  )
  expect_error(
    recode(GT_dominant, geneticModel = "allelic"),
    regexp = "should be 'allelic'"
  )

  # expect non-imputed genoMatrix when genoMatrix is specified
  GT_imputed <- recode(GT, imputeMethod = "meanImpute")
  expect_error(
    recode(GT_imputed, geneticModel = "dominant"),
    regexp = "Provide a non-imputed genoMatrix"
  )
  
  # expect message when GT is already imputed
  expect_message(
    recode(GT_imputed, imputeMethod = "meanImpute"),
    regexp = "GT is already imputed"
  )
  
  # expect aggregates to be rest
  GT_agg <- aggregate(recode(GT, imputeMethod = "meanImpute"), returnGT = TRUE)
  GT_agg_recoded <- recode(GT_agg, weights = rowData(GT_agg)$w * 2)
  expect_true(all(is.na(GT_agg_recoded$aggregate)))
})

test_that("genoMatrix aggregate input validation works", {
  # expect error when GT contains missing values
  expect_error(
    aggregate(GT),
    regexp = "contains missing values"
  )

  # expect error when all weights are missing or infinite
  GT_missing_weights <- GT
  GT_missing_weights <- recode(GT_missing_weights, imputeMethod = "meanImpute")
  rowData(GT_missing_weights)$w <- NA_real_
  expect_error(
    aggregate(GT_missing_weights),
    regexp = "All weights are missing or infinite"
  )

  # expact warning when some weights are missing
  rowData(GT_missing_weights)$w[1:2] <- c(1, 2)
  expect_warning(
    aggregate(GT_missing_weights),
    regexp = "contain missing or infinite weights"
  )
})
