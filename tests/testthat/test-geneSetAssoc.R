data(rvbresults)
genesetfile_c5 <- withr::local_tempfile()
suppressMessages(genesetfile <- buildGeneSet(gmtpath = "../data/c5.go.mf.v2023.2.Hs.symbols.gmt",
                              output = genesetfile_c5))
genesetfile <- geneSetFile(genesetfile_c5)
genesetlist <- as.geneSetList(genesetfile)

test_that("prepareStatsGSA works",{

  # expect error when duplicated units are included
  expect_error({
    rvbresults_zcutoffs <- .prepare_stats_GSA(rvbresults, Zcutoffs = c(-3, 4), INT = FALSE, covar = NULL)
  })

  # check Z-cutoffs
  rvbresults_moderate <- rvbresults[rvbresults$test == "firth" & rvbresults$varSetName == "ModerateImpact", ]
  zscores <- qnorm(1-rvbresults_moderate$P)
  messages <- testthat::capture_messages(
    {rvbresults_zcutoffs <- .prepare_stats_GSA(rvbresults_moderate, Zcutoffs = c(-3, 4), INT = FALSE, covar = NULL)}
  )
  expect_equal(messages[1], sprintf("%s Z-scores <-3 are set to -3\n", sum(zscores < -3, na.rm = TRUE)))
  expect_equal(messages[2], sprintf("%s Z-scores >4 are set to 4\n", sum(zscores > 4, na.rm = TRUE)))
  
  rvbresults_moderate_zcutoffs <- .prepare_stats_GSA(rvbresults_moderate, Zcutoffs = c(-3, 4), INT = FALSE, covar = NULL, verbose = FALSE)
  expect_true(all(rvbresults_moderate_zcutoffs$Z[zscores > 4] == 4))
  expect_true(all(rvbresults_moderate_zcutoffs$Z[zscores < -3] == -3))

  # check INT
  rvbresults_moderate_INT <- .prepare_stats_GSA(rvbresults_moderate, Zcutoffs = NULL, INT = TRUE, covar = NULL, verbose = FALSE)
  expect_equal(mean(rvbresults_moderate_INT$Z), 0, tolerance = 1e-4)
  expect_equal(var(rvbresults_moderate_INT$Z), 1, tolerance = 1e-4)

  # expect a warning when both INT and Z-cutoffs are supplied
  expect_warning(rvbresults_moderate_INT2 <- .prepare_stats_GSA(rvbresults_moderate, Zcutoffs = c(-3, 4), INT = TRUE, covar = NULL))
  expect_equal(rvbresults_moderate_INT, rvbresults_moderate_INT2)

  # check fixing infinites
  rvbresults_moderate <- rvbresults[rvbresults$test == "skat_burden_robust" & rvbresults$varSetName == "ModerateImpact", ]
  zscores <- qnorm(1-rvbresults_moderate$P)

  rvbresults_moderate_notrans <- suppressMessages(.prepare_stats_GSA(rvbresults_moderate, Zcutoffs = NULL, INT = FALSE, covar = NULL))
  expect_true({
    all(dplyr::near(rvbresults_moderate_notrans$Z[is.infinite(zscores) & zscores > 0], max(zscores[!is.infinite(zscores)])))
  })
  expect_true({
    all(dplyr::near(rvbresults_moderate_notrans$Z[is.infinite(zscores) & zscores < 0], min(zscores[!is.infinite(zscores)])))
  })

  ## check adding covars
  rvbresults_moderate$random_covar <- sample(c("a", "b", NA_character_), size = nrow(rvbresults_moderate), replace = TRUE)
 suppressMessages(rvbresults_moderate_prep <- .prepare_stats_GSA(rvbresults_moderate, Zcutoffs = NULL, INT = FALSE, covar = c("nvar", "random_covar")))
  expect_equal(nrow(rvbresults_moderate_prep), (nrow(rvbresults_moderate) - sum(is.na(rvbresults_moderate$random_covar))))

}
)

test_that("geneSetAssoc snapshots are equal",{

  # results
  res <- rvbresults[rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth", ]

  # snapshot run
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
  
  for(i in 1:length(results_list)) {
    metadata(results_list[[i]])$creationDate <- NA_character_
    metadata(results_list[[i]])$gdbPath <- NA_character_
  }
  expect_snapshot_value(results_list, style = "serialize")
  }
)


test_that("geneSetAssoc works",{

  res <- rvbresults[rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth", ]
  res$carriers <- (res$caseCarriers + res$ctrlCarriers)

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
  gsa_genesetfile2 <- geneSetAssoc(
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
  gsa_genesetfile2$ID <- paste(gsa_genesetfile2$geneSetName, gsa_genesetfile2$test)
  metadata(gsa_genesetfile2)$creationDate <- NA_character_
  expect_equal(gsa_genesetfile, gsa_genesetfile2[match(as.character(gsa_genesetfile$ID), as.character(gsa_genesetfile2$ID)),])

  # check summary gsaResult
  expect_true(stringr::str_detect(capture_output({summary(gsa_genesetlist)}), "n gene sets = "))

  # perform tests manually for a couple of genesets
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
  res <- .prepare_stats_GSA(res, covar = c("nvar", "carriers"), Zcutoffs = NULL, INT = FALSE, verbose = FALSE)

  ## linear model
  run_linear_model <- function(x, gslist, results) {
    units <- listUnits(getGeneSet(gslist, geneSet = x))
    results$in_geneset <- results$unit %in% units
    model <- lm(Z ~ in_geneset + nvar + carriers, data = as.data.frame(results))
    tibble(geneSetName = x, effect = summary(model)$coef["in_genesetTRUE","Estimate"], P = summary(model)$coef["in_genesetTRUE","Pr(>|t|)"])
  }
  lm_results <- dplyr::bind_rows(lapply(genesets, FUN = run_linear_model, gslist = genesetlist, results = res))
  lm_results <- lm_results %>%
    dplyr::left_join(as.data.frame(gsa) %>% dplyr::filter(test=="lm") %>% dplyr::select(geneSetName, effect, P), by = "geneSetName")
  expect_equal(lm_results$effect.x, lm_results$effect.y, tolerance = 1e-4)
  expect_equal(lm_results$P.x, lm_results$P.y, tolerance = 1e-4)

  ## fisher
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
        ), nrow = 2, byrow = TRUE
      )
    )
    tibble(geneSetName = x, effect = test$estimate, P = test$p.value)
  }
  fisher_results <- dplyr::bind_rows(lapply(genesets, FUN = run_fisher, gslist = genesetlist, results = res, threshold = 1e-4))
  fisher_results <- fisher_results %>%
    dplyr::left_join(as.data.frame(gsa) %>% dplyr::filter(test=="fisher") %>% dplyr::select(geneSetName, effect, P), by = "geneSetName")
  expect_equal(lm_results$effect.x, lm_results$effect.y, tolerance = 1e-4)
  expect_equal(lm_results$P.x, lm_results$P.y, tolerance = 1e-4)

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
  expect_equal(gsa_onesided[gsa_onesided$effect > 0 & gsa_onesided$test == "lm",]$P*2,
               gsa_twosided[gsa_twosided$effect > 0 & gsa_twosided$test == "lm",]$P
               )
  expect_equal(gsa_onesided[gsa_onesided$effect > 0 & gsa_onesided$test == "lm",]$effect,
               gsa_twosided[gsa_twosided$effect > 0 & gsa_twosided$test == "lm",]$effect
  )
  expect_equal(gsa_onesided[gsa_onesided$effect < 0 & gsa_onesided$test == "lm",]$effect,
               gsa_twosided[gsa_twosided$effect < 0 & gsa_twosided$test == "lm",]$effect
  )
}
)


# Cell-type enrichments
#
test_that("cell-type enrichment snapshots are equal",{

  # prepare data
  res <- rvbresults[rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth", ]
  gene_anno <- readr::read_tsv("../data/Homo_sapiens.GRCh38.105.gene.txt.gz", show_col_types = FALSE, progress = FALSE)

  ## from https://github.com/Kyoko-wtnb/FUMA_scRNA_data/tree/master/processed_data
  data <- readr::read_tsv("../data/GSE67835_Human_Cortex.txt.gz", show_col_types = FALSE, progress = FALSE)
  genes <- data$GENE
  data$GENE <- NULL
  data <- as.matrix(data)
  rownames(data) <- genes

  res <- merge(res, gene_anno[,c("gene_id", "gene_name")] %>% dplyr::filter(!duplicated(gene_name)), by = c("unit" = "gene_name"))
  res$unit <- res$gene_id
  res <- res[res$unit %in% rownames(data),]
  data <- data[rownames(data) %in% res$unit,]
  data <- data[as.character(res$unit),]
  res$average <- data[,"Average"]
  data <- data[,colnames(data) != "Average"]

  # run tests
  results_list <- list()
  results_list[["test1"]] <- suppressMessages(geneSetAssoc(
    res,
    scoreMatrix = data,
    covar = c("nvar", "average"),
    test = "lm" ))

  results_list[["test2"]] <- suppressMessages(geneSetAssoc(
    res,
    scoreMatrix = data,
    covar = c("nvar", "average"),
    INT = TRUE,
    test = "lm" ))

  results_list[["test3"]] <- suppressMessages(geneSetAssoc(
    res,
    scoreMatrix = data,
    covar = c("nvar", "average"),
    test = "lm",
    Zcutoffs = c(-3, 4)
    ))

  expect_snapshot_value(results_list, style = "serialize")


  ## manual
  ce <- geneSetAssoc(
    res,
    scoreMatrix = data,
    covar = c("nvar", "average"),
    test = "lm",
    oneSided = FALSE,
    verbose = FALSE
    )
  res <- .prepare_stats_GSA(res, covar = c("nvar", "average"), Zcutoffs = NULL, INT = FALSE, verbose = FALSE)

  ## linear model
  run_linear_model <- function(x, mat, results) {
    row <- data.frame(unit = rownames(data), expr = data[, x,drop=TRUE])
    results <- as.data.frame(results) %>%
      dplyr::left_join(row, by = "unit")
    model <- lm(Z ~ expr + nvar + average, data = results)
    tibble(name = x, effect = summary(model)$coef["expr","Estimate"], P = summary(model)$coef["expr","Pr(>|t|)"])
  }
  lm_results <- dplyr::bind_rows(lapply(colnames(data), FUN = run_linear_model, mat = data, results = res))
  lm_results <- lm_results %>%
    dplyr::left_join(as.data.frame(ce) %>% dplyr::filter(test == "lm") %>% dplyr::select(name, effect, P), by = "name")
  expect_equal(lm_results$effect.x, lm_results$effect.y, tolerance = 1e-4)
  expect_equal(lm_results$P.x, lm_results$P.y, tolerance = 1e-4)
  }
)

