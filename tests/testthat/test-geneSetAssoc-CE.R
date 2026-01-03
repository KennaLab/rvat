data(rvbresults, envir = environment())
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

test_that("cell-type enrichment snapshots are equal", {
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

  # check application of scoreCutoffs
  ce <- geneSetAssoc(
    res,
    scoreMatrix = data,
    covar = c("nvar", "average"),
    test = "lm",
    scoreCutoffs = c(1, 1),
    oneSided = FALSE,
    verbose = FALSE
  )

  # compare writing to output with interactive
  ce_interactive <- geneSetAssoc(
    res,
    scoreMatrix = data,
    covar = c("nvar", "average"),
    test = "lm",
    oneSided = FALSE,
    verbose = FALSE
  )
  output <- withr::local_tempfile()
  geneSetAssoc(
    res,
    scoreMatrix = data,
    covar = c("nvar", "average"),
    test = "lm",
    oneSided = FALSE,
    verbose = FALSE,
    output = output
  )
  ce_from_file <- readr::read_tsv(output, progress = FALSE, show_col_types = FALSE)
  rownames(ce_interactive) <- NULL
  rownames(ce_from_file) <- NULL
  expect_equal(ce_interactive, as.data.frame(ce_from_file))
})


test_that("Cell-type enrichment input validation works", {
  # expect messages when some units are not present in scoreMatrix or vice versa
  messages <- capture_messages({
    geneSetAssoc(
      res,
      scoreMatrix = data[1:100, ],
      covar = c("nvar", "average"),
      test = "lm",
      oneSided = FALSE,
      verbose = TRUE
    )
  })
  expect_true(any(stringr::str_detect(
    messages,
    "units in the results are not present in the scoreMatrix"
  )))

  data_nonoverlap <- data
  rownames(data_nonoverlap)[1:10] <- LETTERS[1:10]
  messages <- capture_messages({
    geneSetAssoc(
      res,
      scoreMatrix = data_nonoverlap[1:100, ],
      covar = c("nvar", "average"),
      test = "lm",
      oneSided = FALSE,
      verbose = TRUE
    )
  })
  expect_true(any(stringr::str_detect(
    messages,
    "are not present in the results"
  )))

  # expect error when scoreCutoffs is invalid
  expect_error(
    {
      geneSetAssoc(
        res,
        scoreMatrix = data,
        covar = c("nvar", "average"),
        test = "lm",
        scoreCutoffs = 1,
        oneSided = FALSE,
        verbose = TRUE
      )
    },
    regexp = "`scoreCutoffs` should be a vector of length 2"
  )

  expect_error(
    {
      geneSetAssoc(
        res,
        scoreMatrix = data,
        covar = c("nvar", "average"),
        test = "lm",
        scoreCutoffs = c(-1, -1),
        oneSided = FALSE,
        verbose = TRUE
      )
    },
    regexp = "scoreCutoffs should be a positive vector"
  )
})
