data(rvbresults, envir = environment())

test_that("cell-type enrichment snapshots are equal", {
  # prepare data
  res <- rvbresults[
    rvbresults$varSetName == "ModerateImpact" & rvbresults$test == "firth",
  ]
  gene_anno <- readr::read_tsv(
    "../data/Homo_sapiens.GRCh38.105.gene.txt.gz",
    show_col_types = FALSE,
    progress = FALSE
  )

  ## from https://github.com/Kyoko-wtnb/FUMA_scRNA_data/tree/master/processed_data
  data <- readr::read_tsv(
    "../data/GSE67835_Human_Cortex.txt.gz",
    show_col_types = FALSE,
    progress = FALSE
  )
  genes <- data$GENE
  data$GENE <- NULL
  data <- as.matrix(data)
  rownames(data) <- genes

  res <- merge(
    res,
    gene_anno[, c("gene_id", "gene_name")] %>%
      dplyr::filter(!duplicated(gene_name)),
    by = c("unit" = "gene_name")
  )
  res$unit <- res$gene_id
  res <- res[res$unit %in% rownames(data), ]
  data <- data[rownames(data) %in% res$unit, ]
  data <- data[as.character(res$unit), ]
  res$average <- data[, "Average"]
  data <- data[, colnames(data) != "Average"]

  # run tests
  results_list <- list()
  results_list[["test1"]] <- suppressMessages(geneSetAssoc(
    res,
    scoreMatrix = data,
    covar = c("nvar", "average"),
    test = "lm"
  ))

  results_list[["test2"]] <- suppressMessages(geneSetAssoc(
    res,
    scoreMatrix = data,
    covar = c("nvar", "average"),
    INT = TRUE,
    test = "lm"
  ))

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
  res <- rvat:::.prepare_stats_GSA(
    res,
    covar = c("nvar", "average"),
    Zcutoffs = NULL,
    INT = FALSE,
    verbose = FALSE
  )

  ## linear model
  run_linear_model <- function(x, mat, results) {
    row <- data.frame(unit = rownames(data), expr = data[, x, drop = TRUE])
    results <- as.data.frame(results) %>%
      dplyr::left_join(row, by = "unit")
    model <- lm(Z ~ expr + nvar + average, data = results)
    tibble(
      name = x,
      effect = summary(model)$coef["expr", "Estimate"],
      P = summary(model)$coef["expr", "Pr(>|t|)"]
    )
  }
  lm_results <- dplyr::bind_rows(lapply(
    colnames(data),
    FUN = run_linear_model,
    mat = data,
    results = res
  ))
  lm_results <- lm_results %>%
    dplyr::left_join(
      as.data.frame(ce) %>%
        dplyr::filter(test == "lm") %>%
        dplyr::select(name, effect, P),
      by = "name"
    )
  expect_equal(lm_results$effect.x, lm_results$effect.y, tolerance = 1e-4)
  expect_equal(lm_results$P.x, lm_results$P.y, tolerance = 1e-4)
})

