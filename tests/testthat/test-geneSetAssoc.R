data(rvbresults, envir = environment())
genesetfile_c5 <- withr::local_tempfile()
suppressMessages(
  genesetfile <- buildGeneSet(
    gmtpath = "../data/c5.go.mf.v2023.2.Hs.symbols.gmt",
    output = genesetfile_c5
  )
)
genesetfile <- geneSetFile(genesetfile_c5)
genesetlist <- as.geneSetList(genesetfile)
rvbresults_moderate <- rvbresults[
  rvbresults$test == "firth" & rvbresults$varSetName == "ModerateImpact",
]

test_that("prepareStatsGSA works", {
  # check if Z-cutoffs are messaged correctly
  zscores <- qnorm(1 - rvbresults_moderate$P)
  messages <- capture_messages(
    {
      rvbresults_zcutoffs <- rvat:::.prepare_stats_GSA(
        rvbresults_moderate,
        Zcutoffs = c(-3, 4),
        INT = FALSE,
        covar = NULL
      )
    }
  )
  expect_equal(
    messages[1],
    sprintf("%s Z-scores <-3 are set to -3\n", sum(zscores < -3, na.rm = TRUE))
  )
  expect_equal(
    messages[2],
    sprintf("%s Z-scores >4 are set to 4\n", sum(zscores > 4, na.rm = TRUE))
  )

  # check if Z-cutoffs are applied correctly
  rvbresults_moderate_zcutoffs <- rvat:::.prepare_stats_GSA(
    rvbresults_moderate,
    Zcutoffs = c(-3, 4),
    INT = FALSE,
    covar = NULL,
    verbose = FALSE
  )
  expect_true(all(rvbresults_moderate_zcutoffs$Z[zscores > 4] == 4))
  expect_true(all(rvbresults_moderate_zcutoffs$Z[zscores < -3] == -3))

  # check if inverse normal transformation works as expected
  rvbresults_moderate_INT <- rvat:::.prepare_stats_GSA(
    rvbresults_moderate,
    Zcutoffs = NULL,
    INT = TRUE,
    covar = NULL,
    verbose = FALSE
  )
  expect_equal(mean(rvbresults_moderate_INT$Z), 0, tolerance = 1e-4)
  expect_equal(var(rvbresults_moderate_INT$Z), 1, tolerance = 1e-4)

  # expect a warning when both INT and Z-cutoffs are supplied
  expect_warning(
    rvbresults_moderate_INT2 <- rvat:::.prepare_stats_GSA(
      rvbresults_moderate,
      Zcutoffs = c(-3, 4),
      INT = TRUE,
      covar = NULL
    )
  )
  expect_equal(rvbresults_moderate_INT, rvbresults_moderate_INT2)

  # check if infinite Z-scores are handled correctly
  zscores <- qnorm(1 - rvbresults_moderate$P)
  rvbresults_moderate_notrans <- suppressMessages(rvat:::.prepare_stats_GSA(
    rvbresults_moderate,
    Zcutoffs = NULL,
    INT = FALSE,
    covar = NULL
  ))
  expect_true({
    all(dplyr::near(
      rvbresults_moderate_notrans$Z[is.infinite(zscores) & zscores > 0],
      max(zscores[!is.infinite(zscores)])
    ))
  })
  expect_true({
    all(dplyr::near(
      rvbresults_moderate_notrans$Z[is.infinite(zscores) & zscores < 0],
      min(zscores[!is.infinite(zscores)])
    ))
  })

  ## check if missing covariate values are handled correctly
  rvbresults_moderate$random_covar <- sample(
    c("a", "b", NA_character_),
    size = nrow(rvbresults_moderate),
    replace = TRUE
  )
  suppressMessages(
    rvbresults_moderate_prep <- rvat:::.prepare_stats_GSA(
      rvbresults_moderate,
      Zcutoffs = NULL,
      INT = FALSE,
      covar = c("nvar", "random_covar")
    )
  )
  expect_equal(
    nrow(rvbresults_moderate_prep),
    (nrow(rvbresults_moderate) - sum(is.na(rvbresults_moderate$random_covar)))
  )

  # expect error when duplicated units are included
  expect_error({
    rvbresults_zcutoffs <- rvat:::.prepare_stats_GSA(
      rvbresults,
      Zcutoffs = c(-3, 4),
      INT = FALSE,
      covar = NULL
    )
  })
})

test_that("geneSetAssoc snapshots are equal", {
  res <- rvbresults[
    rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",
  ]
  # snapshot tests for various param combinations
  results_list <- list()
  results_list[["test1"]] <- suppressMessages(geneSetAssoc(
    res,
    genesetlist[1:100],
    covar = c("nvar"),
    test = c("fisher", "lm", "ttest"),
    minSetSize = 10,
    maxSetSize = 500,
    threshold = 1e-4
  ))

  results_list[["test2"]] <- suppressMessages(geneSetAssoc(
    res,
    genesetlist[1:100],
    covar = c("nvar"),
    test = c("fisher", "lm", "ttest"),
    minSetSize = 10,
    maxSetSize = 500,
    Zcutoffs = c(-3, 4),
    threshold = 1e-4
  ))

  results_list[["test3"]] <- suppressMessages(geneSetAssoc(
    res,
    genesetlist[1:100],
    covar = c("nvar"),
    test = c("fisher", "lm", "ttest"),
    minSetSize = 10,
    maxSetSize = 500,
    INT = TRUE,
    threshold = 1e-4
  ))

  results_list[["test4"]] <- suppressMessages(geneSetAssoc(
    res,
    genesetlist[1:10],
    condition = genesetlist[11:20],
    covar = c("nvar"),
    test = c("fisher", "lm", "ttest"),
    minSetSize = 10,
    maxSetSize = 500,
    INT = TRUE,
    threshold = 1e-4
  ))

  results_list[["test5"]] <- suppressMessages(geneSetAssoc(
    res,
    genesetlist[1:10],
    covar = NULL,
    test = c("fisher", "lm", "ttest"),
    minSetSize = 10,
    maxSetSize = 2000,
    threshold = 1e-4,
    Zcutoffs = c(-3, 2),
  ))

  res$carriers <- (res$caseCarriers + res$ctrlCarriers)
  results_list[["test6"]] <- suppressMessages(geneSetAssoc(
    res,
    genesetlist[1000:1020],
    covar = c("nvar", "carriers"),
    test = c("fisher", "lm", "ttest"),
    minSetSize = 10,
    maxSetSize = 2000,
    threshold = 1e-4,
    Zcutoffs = c(-3, 2),
    oneSided = FALSE
  ))

  for (i in 1:length(results_list)) {
    metadata(results_list[[i]])$creationDate <- NA_character_
    metadata(results_list[[i]])$gdbPath <- NA_character_
    metadata(results_list[[i]])$rvatVersion <- NA_character_
  }
  expect_snapshot_value(results_list, style = "serialize")
})

test_that("geneSetAssoc results match manual calculations", {
  # reference geneSetAssoc results
  res <- rvbresults[
    rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",
  ]
  res$carriers <- (res$caseCarriers + res$ctrlCarriers)
  gsa_genesetlist <- geneSetAssoc(
    res,
    genesetlist,
    covar = c("nvar", "carriers"),
    test = c("fisher", "lm"),
    minSetSize = 10,
    maxSetSize = 500,
    threshold = 1e-4,
    verbose = FALSE
  )
  genesets <- sample(unique(gsa_genesetlist$geneSetName), size = 10)
  gsa <- geneSetAssoc(
    res,
    getGeneSet(genesetlist, geneSet = genesets),
    covar = c("nvar", "carriers"),
    test = c("fisher", "lm"),
    minSetSize = 10,
    maxSetSize = 500,
    threshold = 1e-4,
    oneSided = FALSE,
    verbose = FALSE
  )

  ## prepare stats for manual calculations
  res_prep <- rvat:::.prepare_stats_GSA(
    res,
    covar = c("nvar", "carriers"),
    Zcutoffs = NULL,
    INT = FALSE,
    verbose = FALSE
  )

  ## run linear model
  run_linear_model <- function(x, gslist, results) {
    units <- listUnits(getGeneSet(gslist, geneSet = x))
    results$in_geneset <- results$unit %in% units
    model <- lm(Z ~ in_geneset + nvar + carriers, data = as.data.frame(results))
    tibble(
      geneSetName = x,
      effect = summary(model)$coef["in_genesetTRUE", "Estimate"],
      P = summary(model)$coef["in_genesetTRUE", "Pr(>|t|)"]
    )
  }
  lm_results <- dplyr::bind_rows(lapply(
    genesets,
    FUN = run_linear_model,
    gslist = genesetlist,
    results = res_prep
  ))
  lm_results <- lm_results %>%
    dplyr::left_join(
      as.data.frame(gsa) %>%
        dplyr::filter(test == "lm") %>%
        dplyr::select(geneSetName, effect, P),
      by = "geneSetName"
    )
  expect_equal(lm_results$effect.x, lm_results$effect.y, tolerance = 1e-4)
  expect_equal(lm_results$P.x, lm_results$P.y, tolerance = 1e-4)

  ## run fisher tests
  run_fisher <- function(x, gslist, results, threshold) {
    units <- listUnits(getGeneSet(gslist, geneSet = x))
    results$in_geneset <- results$unit %in% units
    results$sig <- results$P < threshold
    test <- fisher.test(
      matrix(
        c(
          sum(results$sig & results$in_geneset, na.rm = TRUE),
          sum(results$sig & !results$in_geneset, na.rm = TRUE),
          sum(!results$sig & results$in_geneset, na.rm = TRUE),
          sum(!results$sig & !results$in_geneset, na.rm = TRUE)
        ),
        nrow = 2,
        byrow = TRUE
      )
    )
    tibble(geneSetName = x, effect = test$estimate, P = test$p.value)
  }
  fisher_results <- dplyr::bind_rows(lapply(
    genesets,
    FUN = run_fisher,
    gslist = genesetlist,
    results = res_prep,
    threshold = 1e-4
  ))
  fisher_results <- fisher_results %>%
    dplyr::left_join(
      as.data.frame(gsa) %>%
        dplyr::filter(test == "fisher") %>%
        dplyr::select(geneSetName, effect, P),
      by = "geneSetName"
    )
  expect_equal(lm_results$effect.x, lm_results$effect.y, tolerance = 1e-4)
  expect_equal(lm_results$P.x, lm_results$P.y, tolerance = 1e-4)
})

test_that("geneSetAssoc misc", {
  res <- rvbresults[
    rvbresults$varsetname == "moderateimpact" & rvbresults$test == "firth",
  ]
  res$carriers <- (res$casecarriers + res$ctrlcarriers)

  # expect identical output geneSetList vs. geneSetFile
  gsa_genesetlist <- geneSetAssoc(
    res,
    genesetlist,
    covar = c("nvar", "carriers"),
    test = c("fisher", "lm"),
    minSetSize = 10,
    maxSetSize = 500,
    threshold = 1e-4,
    verbose = FALSE
  )
  gsa_genesetfile <- geneSetAssoc(
    res,
    genesetfile,
    covar = c("nvar", "carriers"),
    test = c("fisher", "lm"),
    minSetSize = 10,
    maxSetSize = 500,
    threshold = 1e-4,
    verbose = FALSE
  )
  metadata(gsa_genesetlist)$creationDate <- NA_character_
  metadata(gsa_genesetfile)$creationDate <- NA_character_
  expect_equal(gsa_genesetlist, gsa_genesetfile)

  # expect identical output with different chunk size
  gsa_genesetfile_check_chunksize <- geneSetAssoc(
    res,
    genesetfile,
    covar = c("nvar", "carriers"),
    test = c("fisher", "lm"),
    minSetSize = 10,
    maxSetSize = 500,
    memlimit = 100,
    threshold = 1e-4,
    verbose = FALSE
  )
  gsa_genesetfile$ID <- paste(gsa_genesetfile$geneSetName, gsa_genesetfile$test)
  gsa_genesetfile_check_chunksize$ID <- paste(
    gsa_genesetfile_check_chunksize$geneSetName,
    gsa_genesetfile_check_chunksize$test
  )
  metadata(gsa_genesetfile_check_chunksize)$creationDate <- NA_character_
  expect_equal(
    gsa_genesetfile,
    gsa_genesetfile_check_chunksize[
      match(
        as.character(gsa_genesetfile$ID),
        as.character(gsa_genesetfile_check_chunksize$ID)
      ),
    ]
  )

  # check summary gsaResult
  expect_true(stringr::str_detect(
    capture_output({
      summary(gsa_genesetlist)
    }),
    "n gene sets = "
  ))

  # check one-sided vs. two-sided tests
  gsa_onesided <- geneSetAssoc(
    res,
    genesetlist[1:100],
    covar = c("nvar", "carriers"),
    test = c("fisher", "lm", "ttest"),
    minSetSize = 10,
    maxSetSize = 500,
    threshold = 1e-2,
    oneSided = TRUE,
    verbose = FALSE
  )

  gsa_twosided <- geneSetAssoc(
    res,
    genesetlist[1:100],
    covar = c("nvar", "carriers"),
    test = c("fisher", "lm", "ttest"),
    minSetSize = 10,
    maxSetSize = 500,
    threshold = 1e-2,
    oneSided = FALSE,
    verbose = FALSE
  )
  expect_equal(
    gsa_onesided[gsa_onesided$effect > 0 & gsa_onesided$test == "lm", ]$P * 2,
    gsa_twosided[gsa_twosided$effect > 0 & gsa_twosided$test == "lm", ]$P
  )
  expect_equal(
    gsa_onesided[gsa_onesided$effect > 0 & gsa_onesided$test == "lm", ]$effect,
    gsa_twosided[gsa_twosided$effect > 0 & gsa_twosided$test == "lm", ]$effect
  )
  expect_equal(
    gsa_onesided[gsa_onesided$effect < 0 & gsa_onesided$test == "lm", ]$effect,
    gsa_twosided[gsa_twosided$effect < 0 & gsa_twosided$test == "lm", ]$effect
  )
})
