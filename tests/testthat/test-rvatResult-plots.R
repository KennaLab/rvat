# helper to filter down rvb results to one row per gene
get_filtered_rvb <- function() {
  data(rvbresults, envir = environment())
  rvbresults[
    rvbresults$varSetName == "ModerateImpact" &
      rvbresults$test == "firth",
  ]
}

test_that("qqplot works", {
  rvb_filtered <- get_filtered_rvb()
  ## with label
  expect_no_error(qqplot(
    rvb_filtered,
    label = "unit",
    lambda1000 = TRUE
  ))

  ## w/o label
  expect_no_error(qqplot(
    rvb_filtered,
    lambda1000 = TRUE
  ))
})


test_that("manhattan works", {
  rvb_filtered <- get_filtered_rvb()

  expect_no_error(manhattan(
    rvb_filtered,
    label = "unit",
    contigs = "GRCh38"
  ))

  ## labelRepel
  expect_no_error(manhattan(
    rvb_filtered,
    label = "unit",
    labelRepel = TRUE,
    contigs = "GRCh38"
  ))
})

test_that("densityPlot works", {
  rvb_filtered <- get_filtered_rvb()
  genesetlist <- buildGeneSet(list("genesetA" = sample(rvb_filtered$unit, 50)))
  expect_no_error({
    suppressMessages(densityPlot(
      rvb_filtered,
      "genesetA",
      genesetlist,
      showMeans = TRUE,
      title = "test"
    ))
  })
})

test_that("manhattan input validation works", {
  rvbresults <- get_filtered_rvb()
  rvbresults$P <- NULL
  expect_error(
    {
      manhattan(
        rvbresults,
        contigs = "GRCh38"
      )
    },
    regexp = "Results should contain a 'P' column"
  )

  ## expect warning when contigs are not specified
  rvbresults <- get_filtered_rvb()
  metadata(rvbresults)$genomeBuild <- NA_character_

  expect_warning(
    {
      man <- manhattan(
        rvbresults
      )
    },
    regexp = "Contigs not specified, defaulting to"
  )
})

test_that("qqplot input validation works", {
  rvbresults <- get_filtered_rvb()

  ## expect warning when misisng P-values are included
  rvbresults$P[1:10] <- NA_real_
  expect_message(
    {
      qqplot(
        rvbresults,
        lambda1000 = TRUE
      )
    },
    regexp = "are excluded"
  )
})
