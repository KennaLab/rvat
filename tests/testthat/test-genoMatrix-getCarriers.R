data(GT)

test_that("getCarriers works", {
  # compare getCarriers to aggregae
  carriers <- getCarriers(GT)
  carriers_agg <- carriers %>%
    dplyr::group_by(IID) %>%
    dplyr::summarize(aggregate = sum(genotype), .groups = "drop")

  agg <- aggregate(recode(GT, imputeMethod = "meanImpute"), returnGT = FALSE)
  agg <- agg[agg >= 1]
  agg <- data.frame(
    IID = names(agg),
    aggregate = round(agg),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(carriers_agg, by = "IID")

  expect_identical(agg$aggregate.x, agg$aggregate.y)
})


test_that("getCarriers with row/colDataFields", {
  # test adding colDataFields
  carriers <- getCarriers(GT, colDataFields = c("sex", "age"))
  expect_true(all(c("sex", "age") %in% colnames(carriers)))

  # test adding rowDataFields
  carriers <- getCarriers(GT, rowDataFields = c("REF", "ALT"))
  expect_true(all(c("REF", "ALT") %in% colnames(carriers)))
})

test_that("getCarriers groupBy parameter works", {
  # test groupBy
  carriers <- getCarriers(GT)
  carriers_grouped <- getCarriers(GT, groupBy = "superPop")
  expect_true(sum(carriers_grouped$carrierN) == nrow(carriers))
  expect_true("superPop" %in% colnames(carriers_grouped))
})

test_that("getCarriers aggregate parameter works", {
  # test aggregate param by comparing to directly running aggregate per group
  carriers <- getCarriers(GT, aggregate = TRUE, groupBy = "superPop")
  GT_imputed <- recode(GT, imputeMethod = "meanImpute")
  check <- unlist(lapply(unique(GT_imputed$superPop), FUN = function(x) {
    mean(aggregate(GT_imputed[, GT_imputed$superPop == x])$aggregate)
  }))
  names(check) <- unique(GT_imputed$superPop)

  carriers <- carriers %>%
    dplyr::left_join(
      tibble::tibble(superPop = names(check), aggregate = unname(check)),
      by = "superPop"
    )
  expect_identical(carriers$meanBurdenScore, carriers$aggregate)
})

test_that("getCarriers input validation works", {
  # expect error when not all specified fields are present
  expect_error(
    {
      carriers <- getCarriers(GT, colDataFields = c("A"))
    },
    regexp = "Not all specified colDataFields are present in"
  )
  expect_error(
    {
      carriers <- getCarriers(GT, rowDataFields = c("A"))
    },
    regexp = "Not all specified rowDataFields are present in"
  )
  
  # expect error when non-existent VAR_ids are provided
  expect_error(
    {
      getCarriers(GT, VAR_id = c("a", "b"))
      
    },
    regexp = "Not all specified VAR_ids are present in the genotype matrix."
  )
  
  # expect error when groupBy field is not present
  expect_error(
  {
      getCarriers(GT, aggregate = TRUE, groupBy = "a")
    },
    regexp = "Not all specified `groupBy`"
  ) 
  
  # expect error when aggregate = TRUE, but no groupBy provided
  expect_error(
  {
      getCarriers(GT, aggregate = TRUE)
    },
    regexp = "The `groupBy` argument should be specified"
  ) 
  
  # expect empty data.frame when no carriers found
  carriers <- getCarriers(GT)
  varid <- unique(carriers$VAR_id)[1]
  drop <- carriers$IID[carriers$VAR_id == varid]
  carriers_empty <- getCarriers(
    GT[, !colnames(GT) %in% drop],
    VAR_id = varid
  )
  expect_true(nrow(carriers_empty) == 0)
  expect_equal(colnames(carriers_empty), colnames(carriers))
  
})
