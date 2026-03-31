test_that("ACAT snapshots are consistent", {
  data(rvbresults, envir = environment())
  rvbresults_ACAT <- suppressMessages(ACAT(
    rvbresults,
    aggregate = list("test", "varSetName"),
    fixpval = TRUE,
    fixpval_method = "Liu",
    fixpval_maxP = 0.99,
    fixpval_minP = 1e-16
  ))
  metadata(rvbresults_ACAT)$gdbPath <- NA_character_
  metadata(rvbresults_ACAT)$creationDate <- NA_character_
  expect_snapshot_value(rvbresults_ACAT, style = "serialize")
})

test_that("ACAT input validation works", {
  data(rvbresults, envir = environment())
  rvbresults_P0 <- rvbresults[1:50, ]
  rvbresults_P0$P[1:5] <- 0
  expect_warning(
    {
      rvbresults_ACAT <- suppressMessages(ACAT(
        rvbresults_P0,
        aggregate = list("test", "varSetName"),
        fixpval = TRUE,
        fixpval_method = "Liu",
        fixpval_maxP = 0.99,
        fixpval_minP = 1e-16
      ))
    },
    "are exactly 0"
  )
})
