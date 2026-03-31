data(GT)

test_that("show works", {
  expect_true(stringr::str_detect(
    capture_output({
      show(GT)
    }),
    "rvat genoMatrix"
  ))
})

test_that("genoMatrix validation works", {
  expect_true(validObject(GT))

  GT_subset <- GT[1:10, 1:100]
  expect_true(validObject(GT_subset))
  expect_equal(nrow(GT_subset), 10)
  expect_equal(ncol(GT_subset), 100)
})

test_that("genoMatric constructor validation works", {
  gt <- assays(GT)$GT
  sm <- colData(GT)
  anno <- rowData(GT) %>% as.data.frame()
  VAR_id <- rownames(GT)

  # expect error when GT is not a matrix
  expect_error(
    {
      genoMatrix(
        GT = as.data.frame(gt),
        SM = sm,
        VAR_id = VAR_id
      )
    },
    regexp = "`GT` must be a matrix"
  )

  # expect error when SM is not a data.frame/DFrame
  expect_error(
    {
      genoMatrix(
        GT = gt,
        SM = as.matrix(sm),
        VAR_id = VAR_id
      )
    },
    regexp = "`SM` must be a data.frame or DFrame"
  )

  # expect error when IID field is missing from SM
  expect_error(
    {
      genoMatrix(
        GT = gt,
        SM = as.data.frame(sm) %>% dplyr::select(-IID),
        VAR_id = VAR_id
      )
    },
    regexp = "must contain an 'IID' column"
  )

  # expect error when duplicated IIDs are present in SM
  expect_error(
    {
      genoMatrix(
        GT = gt,
        SM = as.data.frame(sm) %>% dplyr::bind_rows(as.data.frame(sm)[1, ]),
        VAR_id = VAR_id
      )
    },
    regexp = "Duplicated IID values found in cohort"
  )

  # expect error when number of samples in SM doesn't match
  expect_error(
    {
      genoMatrix(
        GT = gt,
        SM = sm[1:10, ],
        VAR_id = VAR_id
      )
    },
    regexp = "Number of samples in `SM`"
  )

  # expect error when length of VAR_id doesn't match number of rows in GT
  expect_error(
    {
      genoMatrix(
        GT = gt,
        SM = sm,
        VAR_id = VAR_id[1:2]
      )
    },
    regexp = "The number of rows in `GT`"
  )

  # expect error when anno is not a data.frame/DFrame
  expect_error(
    {
      genoMatrix(
        GT = gt,
        SM = sm,
        VAR_id = VAR_id,
        anno = as.matrix(rowData(GT))
      )
    },
    regexp = "`anno` must be a data.frame"
  )

  # expect error when anno missed VAR_id column
  expect_error(
    {
      genoMatrix(
        GT = gt,
        SM = sm,
        VAR_id = VAR_id,
        anno = as.data.frame(rowData(GT))
      )
    },
    regexp = "`anno` must contain a `VAR_id` column"
  )

  # expect error when weights don't match number of variants
  expect_error(
    {
      genoMatrix(
        GT = gt,
        SM = sm,
        VAR_id = VAR_id,
        w = c(1, 2)
      )
    },
    regexp = "Length of `w`"
  )

  # expect errors when weights are non-numeric
  expect_error(
    {
      genoMatrix(
        GT = gt,
        SM = sm,
        VAR_id = VAR_id,
        w = as.character(c(1, 2))
      )
    },
    regexp = "must be a numeric vector or scalar"
  )

  # expect error when ploidy don't match number of variants
  expect_error(
    {
      genoMatrix(
        GT = gt,
        SM = sm,
        VAR_id = VAR_id,
        ploidy = c("diploid", "diploid")
      )
    },
    regexp = "`ploidy` vector doesn't match `VAR_id` vector"
  )

  # expect message when some sample IIDs are missing
  expect_message(
    {
      genoMatrix(
        GT = gt,
        SM = as.data.frame(sm) %>% dplyr::mutate(IID = replace(IID, 1:5, NA)),
        VAR_id = VAR_id
      )
    },
    regexp = "samples in the gdb are present"
  )
})
