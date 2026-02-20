#' @include rvatResult.R
#' @include gsaResult.R
#' @include geneSet.R

#' @rdname geneSetAssoc
#' @usage NULL
#' @export
setMethod(
  "geneSetAssoc",
  signature = "rvbResult",
  definition = function(
    object,
    geneSet = NULL,
    scoreMatrix = NULL,
    cormatrix = NULL,
    condition = NULL,
    covar = NULL,
    test = c("lm", "mlm", "fisher", "ttest", "ztest", "ACAT"),
    threshold = NULL,
    Zcutoffs = NULL,
    INT = FALSE,
    scoreCutoffs = NULL,
    minSetSize = 1L,
    maxSetSize = Inf,
    oneSided = TRUE,
    memlimit = 1000L,
    ID = "unit",
    output = NULL,
    verbose = TRUE
  ) {
    # validate input
    .geneSetAssoc_validate_input(as.list(environment()))

    # check methods and available tests
    ## only 'lm' and 'mlm' are implemented for `scoreMatrix` object
    if (!is.null(scoreMatrix)) {
      test <- test[test %in% geneSetAssoc_tests_score]
    }

    # prepare test-stats (Z-scores, covariates, etc.)
    object <- .prepare_stats_GSA(
      object,
      covar,
      Zcutoffs,
      INT,
      verbose = verbose
    )

    # handle MLM null model
    nullmodel <- NULL
    if ("mlm" %in% test) {
      if (oneSided) {
        warning(
          "Note that currently mlm will return two-sided P-values, ",
          "regardless of the `oneSided` argument.",
          call. = FALSE
        )
      }
      if (is.null(cormatrix)) {
        stop(
          "The mlm test requires specifying the 'cormatrix' argument,",
          "see the `buildCorMatrix` method.",
          call. = FALSE
        )
      }

      nullmodel <- fitNullModelGSA(
        object = object,
        cormatrix = cormatrix,
        covar = covar,
        Zcutoffs = Zcutoffs
      )
      object <- getResults(nullmodel)
    }

    # handle condition parameter, if specified
    cond_data <- NULL
    if (!is.null(condition)) {
      cond_data <- .geneSetAssoc_handle_condition_param(
        condition = condition,
        object = object,
        geneSet = geneSet,
        scoreMatrix = scoreMatrix,
        verbose = verbose
      )
      object <- cond_data$object
    }

    ## scoreMatrix -----------------------------------------------------

    if (!is.null(scoreMatrix)) {
      result_scorematrix <- .geneSetAssoc_scoreMatrix(
        object,
        scoreMatrix = scoreMatrix,
        scoreCutoffs = scoreCutoffs,
        nullmodel = nullmodel,
        covar = covar,
        test = test,
        oneSided = oneSided,
        ID = ID,
        memlimit = memlimit,
        verbose = verbose
      )
      if (!is.null(output)) {
        write.table(
          result_scorematrix,
          sep = "\t",
          quote = FALSE,
          file = gzfile(output),
          row.names = FALSE
        )
        return(invisible(NULL))
      } else {
        return(result_scorematrix)
      }
    } else {
      result_gsa <- .geneSetAssoc_gsa(
        object,
        geneSet = geneSet,
        cond_data = cond_data,
        nullmodel = nullmodel,
        covar = covar,
        test = test,
        threshold = threshold,
        oneSided = oneSided,
        ID = ID,
        memlimit = memlimit,
        minSetSize = minSetSize,
        maxSetSize = maxSetSize,
        verbose = verbose
      )
      result_gsa <- gsaResult(result_gsa)
      metadata(result_gsa)$rvatVersion <- as.character(packageVersion("rvat"))
      metadata(result_gsa)$gdbId <- getGdbId(object)
      metadata(result_gsa)$genomeBuild <- getGenomeBuild(object)
      metadata(result_gsa)$creationDate <- as.character(round(
        Sys.time(),
        units = "secs"
      ))

      if (verbose) {
        message(sprintf(
          "%s out of %s sets are kept.",
          length(unique(result_gsa$geneSetName)),
          length(geneSet)
        ))
      }

      if (!is.null(output)) {
        writeResult(result_gsa, file = output)
        return(invisible(NULL))
      } else {
        return(result_gsa)
      }
    }
  }
)

.geneSetAssoc_validate_input <- function(args) {
  check_wrapper(check_character, args, "covar", allow_null = TRUE)
  check_wrapper(check_character, args, "test")
  check_wrapper(check_number_whole, args, "minSetSize")
  check_wrapper(check_number_whole, args, "maxSetSize", allow_infinite = TRUE)
  check_wrapper(check_bool, args, "oneSided")
  check_wrapper(check_bool, args, "INT")
  check_wrapper(check_number_whole, args, "memlimit")
  if (!is.null(args[["threshold"]])) {
    lapply(args[["threshold"]], check_number_decimal)
  }
  if (!is.null(args[["Zcutoffs"]])) {
    lapply(args[["Zcutoffs"]], check_number_decimal)
  }
  if (!is.null(args[["scoreCutoffs"]])) {
    lapply(args[["scoreCutoffs"]], check_number_decimal)
  }
  check_wrapper(check_character, args, "ID", length_equal = 1L)
  check_wrapper(check_bool, args, "verbose")

  # cormatrix is tested within fitNullModelGSA
  if (
    !is.null(args[["geneSet"]]) &&
      !is(args[["geneSet"]], "geneSetList") &&
      !is(args[["geneSet"]], "geneSetFile")
  ) {
    stop("`geneSet` should be a geneSetList or geneSetFile.", call. = FALSE)
  }
  # scoreMatrix should be of type ..
  if (
    !is.null(args[["scoreMatrix"]]) &&
      !is(args[["scoreMatrix"]], "Matrix") &&
      !is(args[["scoreMatrix"]], "matrix")
  ) {
    stop("`scoreMatrix` should be a (sparse) matrix.", call. = FALSE)
  }

  if (!is.null(args[["scoreCutoffs"]])) {
    if (length(args[["scoreCutoffs"]]) != 2L) {
      stop(
        "`scoreCutoffs` should be a vector of length 2 (minimum and maximum sd)",
        call. = FALSE
      )
    }

    if (any(args[["scoreCutoffs"]] < 0L)) {
      stop(
        "scoreCutoffs should be a positive vector (the number of standard deviations below and above the mean)",
        call. = FALSE
      )
    }
  }

  if (!is.null(args[["output"]])) {
    .check_output(
      args[["output"]],
      overWrite = TRUE,
      verbose = args[["verbose"]]
    )
  }

  invisible(NULL)
}

.geneSetAssoc_handle_condition_param <- function(
  condition,
  object,
  geneSet,
  scoreMatrix,
  verbose
) {
  type <- "vector"
  if (is(condition, "geneSetList") || is(condition, "geneSetFile")) {
    type <- "geneSet"
  } else if (is(condition, "matrix")) {
    type <- "matrix"
  }

  if (type == "matrix") {
    nonoverlap <- sum(!object$unit %in% rownames(condition))
    if (nonoverlap > 0L) {
      if (verbose) {
        message(sprintf(
          "%s/%s units in the results are not present in the condition matrix, these are excluded.",
          nonoverlap,
          nrow(object)
        ))
      }
      object <- object[object$unit %in% rownames(condition), ]
    }

    nonoverlap <- sum(!rownames(condition) %in% object$unit)
    if (nonoverlap > 0L) {
      if (verbose) {
        message(sprintf(
          "%s/%s units in the condition matrix are not present in the results, these are excluded",
          nonoverlap,
          nrow(object)
        ))
      }
    }

    # make sure the order is correct
    condition <- condition[as.character(object$unit), , drop = FALSE]
    names_vec <- colnames(condition)

    if (nrow(condition) == 0L || nrow(object) == 0L) {
      stop("No overlapping units left to test!", call. = FALSE)
    }
  } else if (type == "vector") {
    if (!is.null(geneSet)) {
      type <- "geneSet"
      names_vec <- condition[condition %in% names(geneSet)]
      condition <- getGeneSet(geneSet, names_vec)
    } else if (!is.null(scoreMatrix)) {
      type <- "matrix"
      names_vec <- condition[condition %in% colnames(scoreMatrix)]
      condition <- scoreMatrix[, names_vec, drop = FALSE]
    }
  } else if (type == "geneSet") {
    names_vec <- names(condition)
  }
  if (length(names_vec) == 0L) {
    stop("Condition invalid or empty overlap.", call. = FALSE)
  }

  # return
  list(type = type, data = condition, names = names_vec, object = object)
}

.geneSetAssoc_gsa <- function(
  object,
  geneSet,
  cond_data,
  nullmodel,
  covar,
  test,
  threshold,
  oneSided,
  ID,
  memlimit,
  minSetSize,
  maxSetSize,
  verbose
) {
  # handle thresholds
  if (any(geneSetAssoc_tests_competitive_threshold %in% test)) {
    if (is.null(threshold)) {
      threshold <- 0.05 / nrow(object)
      if (verbose) {
        message(sprintf(
          "`threshold` not specified for defining significant genes, using a bonferroni threshold: %s",
          signif(threshold, 4)
        ))
      }
    }
  } else {
    threshold <- NULL
  }

  # split geneset into chunks based on memlimit
  chunks <- split(
    seq_along(geneSet),
    ceiling(seq_along(geneSet) / memlimit)
  )

  results <- lapply(chunks, function(chunk) {
    # extract current chunk
    geneSetList_chunk <- getGeneSet(
      geneSet,
      geneSet = listGeneSets(geneSet)[chunk]
    )
    mappedMatrix <- mapToMatrix(
      geneSetList_chunk,
      object,
      ID = ID,
      sparse = TRUE
    )

    # filter based on set size
    set_sizes <- Matrix::colSums(mappedMatrix)
    sets_keep <- names(set_sizes[
      set_sizes >= minSetSize & set_sizes <= maxSetSize
    ])
    ## return early if no sets left
    if (length(sets_keep) == 0L) {
      empty_result <- gsaResult()
      if (!is.null(cond_data)) {
        empty_result$condition <- character(0)
      }
      return(empty_result)
    }
    ## subset mapped matrix and genesetlist
    mappedMatrix <- mappedMatrix[, sets_keep, drop = FALSE]
    current_sets <- getGeneSet(geneSetList_chunk, geneSet = sets_keep)

    # run GSA (conditional or not)
    if (is.null(cond_data)) {
      gsa(
        as.data.frame(object),
        genesetlist = current_sets,
        mappedMatrix = as.matrix(mappedMatrix),
        nullmodel = nullmodel,
        covar = covar,
        test = test,
        threshold = threshold,
        oneSided = oneSided,
        ID = ID
      )
    } else {
      gsa_conditional(
        object,
        condition = cond_data[["data"]],
        condition.type = cond_data[["type"]],
        genesetlist = current_sets,
        mappedMatrix = as.matrix(mappedMatrix),
        geneSetFile = geneSetFile,
        geneSetList = geneSetList,
        nullmodel = nullmodel,
        covar = covar,
        test = test,
        threshold = threshold,
        oneSided = oneSided,
        ID = ID,
        memlimit = memlimit
      )
    }
  })
  results <- do.call(rbind, results)
  rownames(results) <- NULL
  results
}

.geneSetAssoc_scoreMatrix <- function(
  object,
  scoreMatrix,
  scoreCutoffs,
  nullmodel,
  covar,
  test,
  oneSided,
  ID,
  memlimit,
  verbose
) {
  # check overlap between scoreMatrix and results file
  nonoverlap <- sum(!object$unit %in% rownames(scoreMatrix))
  if (nonoverlap > 0L) {
    if (verbose) {
      message(sprintf(
        "%s/%s units in the results are not present in the scoreMatrix, these are excluded.",
        nonoverlap,
        nrow(object)
      ))
    }
    object <- object[object$unit %in% rownames(scoreMatrix), ]
  }
  nonoverlap <- sum(!rownames(scoreMatrix) %in% object$unit)
  if (nonoverlap > 0L && verbose) {
    message(sprintf(
      "%s/%s units in the scoreMatrix are not present in the results, these are excluded",
      nonoverlap,
      nrow(object)
    ))
  }

  # make sure the order is correct
  scoreMatrix <- scoreMatrix[as.character(object$unit), , drop = FALSE]

  if (nrow(scoreMatrix) == 0L || nrow(object) == 0L) {
    stop("No overlapping units left to test!", call. = FALSE)
  }

  # apply score cutoffs
  if (!is.null(scoreCutoffs)) {
    scoreMatrix <- apply(scoreMatrix, 2, function(x) {
      sdev <- sd(x)
      mn <- mean(x)
      x[x < (mn - (scoreCutoffs[1] * sdev))] <- (mn -
        (scoreCutoffs[1] * sdev))
      x[x > (mn + (scoreCutoffs[2] * sdev))] <- (mn +
        (scoreCutoffs[2] * sdev))
      x
    })
  }

  # split scoreMatrix into chunks based on memlimit
  chunks <- split(
    seq_len(ncol(scoreMatrix)),
    ceiling(seq_len(ncol(scoreMatrix)) / memlimit)
  )

  results <- lapply(chunks, function(chunk) {
    scoreMatrix_chunk <- scoreMatrix[, chunk, drop = FALSE]
    enrich_test(
      as.data.frame(object),
      scorematrix = scoreMatrix_chunk,
      nullmodel = nullmodel,
      covar = covar,
      test = test,
      oneSided = oneSided,
      ID = ID
    )
  })
  results <- do.call(rbind, results)
  rownames(results) <- NULL
  results
}


gsa_conditional <- function(
  object,
  condition,
  condition.type,
  genesetlist,
  mappedMatrix,
  geneSetFile = NULL,
  geneSetList = NULL,
  nullmodel = NULL,
  covar = NULL,
  test = "lm",
  threshold = NULL,
  oneSided = TRUE,
  ID = "unit",
  memlimit = 1000
) {
  Pl <- effectl <- effectSEl <- effectCIlowerl <- effectCIupperl <- list()
  Ngenes_available <- Matrix::colSums(mappedMatrix)

  if ("lm" %in% test) {
    Plt <- effectlt <- effectSElt <- effectCIlowerlt <- effectCIupperlt <- list()
    if (condition.type == "geneSet") {
      chunksCond <- split(
        seq_along(condition),
        ceiling(seq_along(condition) / memlimit)
      )
      condNames <- names(condition)
    } else {
      chunksCond <- list(seq_len(ncol(condition)))
      condNames <- colnames(condition)
    }

    for (chunk in chunksCond) {
      if (condition.type == "geneSet") {
        geneSetList_chunk <- getGeneSet(
          condition,
          geneSet = names(condition)[chunk]
        )
        mappedMatrixCond <- mapToMatrix(
          geneSetList_chunk,
          object,
          ID = ID,
          sparse = TRUE
        )
      } else {
        mappedMatrixCond <- condition
      }

      for (geneset in colnames(mappedMatrixCond)) {
        stats <- cbind(
          as.data.frame(object),
          as.matrix(mappedMatrixCond)[, geneset, drop = FALSE]
        )

        cov <- if (length(covar) == 0L) geneset else c(covar, geneset)
        if (var(stats[, geneset]) == 0) {
          cov <- cov[cov != geneset]
        }

        if (length(cov) > 0L) {
          X <- cbind(1, as.matrix(stats[, c(cov), drop = FALSE]))
        } else {
          X <- cbind(rep(1, nrow(stats)))
        }
        U1 <- crossprod(X, stats$Z)
        U2 <- solve(crossprod(X), U1)
        ytr <- stats$Z - X %*% U2
        U3 <- crossprod(X, as.matrix(mappedMatrix))
        U4 <- solve(crossprod(X), U3)
        Str <- mappedMatrix - X %*% U4
        effect <- as.vector(crossprod(ytr, Str)) / colSums(Str^2)
        Str2 <- colSums(Str^2)
        sig <- (sum(ytr^2) - effect^2 * Str2) / (nrow(stats) - ncol(X) - 2)
        effectSE <- sqrt(sig * (1 / Str2))

        if (oneSided) {
          P <- pt(
            (effect / effectSE),
            nrow(stats) - ncol(X) - 1L,
            lower.tail = FALSE
          )
          fac <- qt(0.95, df = nrow(stats) - ncol(X) - 1L)
          effectCIlower <- effect - (fac * effectSE)
          effectCIupper <- rep(Inf, length(P))
        } else {
          P <- 2 *
            pt(
              abs(effect / effectSE),
              nrow(stats) - ncol(X) - 1L,
              lower.tail = FALSE
            )
          fac <- qt(0.975, df = nrow(stats) - ncol(X) - 1L)
          effectCIlower <- effect - (fac * effectSE)
          effectCIupper <- effect + (fac * effectSE)
        }
        Plt[[geneset]] <- P
        effectlt[[geneset]] <- effect
        effectSElt[[geneset]] <- effectSE
        effectCIlowerlt[[geneset]] <- effectCIlower
        effectCIupperlt[[geneset]] <- effectCIupper
      }
    }
    Plt <- Plt[condNames]
    effectlt <- effectlt[condNames]
    effectSElt <- effectSElt[condNames]
    effectCIlowerlt <- effectCIlowerlt[condNames]
    effectCIupperlt <- effectCIupperlt[condNames]
    Pl[["lm"]] <- unlist(Plt)
    effectl[["lm"]] <- unlist(effectlt)
    effectSEl[["lm"]] <- unlist(effectSElt)
    effectCIlowerl[["lm"]] <- unlist(effectCIlowerlt)
    effectCIupperl[["lm"]] <- unlist(effectCIupperlt)
  }

  Ngenes_available <- Matrix::colSums(mappedMatrix)
  results <- data.frame(
    geneSetName = c(rep(
      names(genesetlist),
      times = (length(test[
        test %in% geneSetAssoc_tests_competitive_condition
      ]) *
        length(condNames))
    )),
    test = c(rep(
      geneSetAssoc_tests_competitive_condition[
        geneSetAssoc_tests_competitive_condition %in% test
      ],
      each = (length(genesetlist) * length(condNames))
    )),
    covar = c(rep(
      paste(covar, collapse = ","),
      times = (length(genesetlist) *
        length(condNames) *
        length(test[test %in% geneSetAssoc_tests_competitive_condition]))
    )),
    condition = c(rep(
      condNames,
      each = (length(test[test %in% geneSetAssoc_tests_competitive_condition]) *
        length(genesetlist))
    )),
    threshold = c(rep(rep(
      NA_real_,
      times = (length(genesetlist) * length(condNames))
    ))),
    geneSetSize = rep(
      lengths(genesetlist),
      times = (length(test[
        test %in% geneSetAssoc_tests_competitive_condition
      ]) *
        length(condNames))
    ),
    genesObs = c(rep(
      Ngenes_available,
      times = (length(test[
        test %in% geneSetAssoc_tests_competitive_condition
      ]) *
        length(condNames))
    )),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  results <- cbind(
    results,
    effect = unname(c(effectl[["lm"]])),
    effectSE = unname(c(effectSEl[["lm"]])),
    effectCIlower = unname(c(effectCIlowerl[["lm"]])),
    effectCIupper = unname(c(effectCIupperl[["lm"]])),
    P = unname(c(Pl[["lm"]]))
  )
  results
}


gsa <- function(
  stats,
  genesetlist,
  mappedMatrix,
  nullmodel = NULL,
  covar = NULL,
  test = c("lm", "mlm", "fisher", "ttest", "ztest", "ACAT"),
  threshold = NULL,
  oneSided = TRUE,
  ID = "unit"
) {
  Pl <- effectl <- effectSEl <- effectCIlowerl <- effectCIupperl <- list()
  results_list <- list()
  Ngenes_available <- Matrix::colSums(mappedMatrix)
  Ngenes <- lengths(genesetlist)

  if (any(test %in% geneSetAssoc_tests_competitive)) {
    if ("lm" %in% test) {
      res <- .geneSetAssoc_compute_lm_stats(
        object,
        mat = mappedMatrix,
        covar = covar,
        oneSided = oneSided
      )
      results_list[["lm"]] <- .geneSetAssoc_format_row(
        names = names(genesetlist),
        test = "lm",
        res = res,
        Ngenes = Ngenes,
        Ngenes_available = Ngenes_available,
        covar = covar
      )
    }

    if ("mlm" %in% test) {
      P <- effect <- effectSE <- effectCIupper <- effectCIlower <- rep(
        NA_real_,
        length(genesetlist)
      )

      if (
        is(getNullModel(nullmodel), "GENESIS.nullMixedModel") ||
          is(getNullModel(nullmodel), "GENESIS.nullModel")
      ) {
        mat <- GWASTools::MatrixGenotypeReader(
          genotype = t(as.matrix(mappedMatrix)),
          snpID = as.integer(1:ncol(mappedMatrix)),
          chromosome = rep(1L, ncol(mappedMatrix)),
          position = seq_len(ncol(mappedMatrix)),
          scanID = seq_len(nrow(mappedMatrix))
        )
        mat <- GWASTools::GenotypeBlockIterator(
          GWASTools::GenotypeData(
            mat,
            scanAnnot = GWASTools::ScanAnnotationDataFrame(
              data.frame(scanID = stats$scanID, stringsAsFactors = FALSE)
            )
          ),
          snpBlock = ncol(mappedMatrix)
        )
        res <- GENESIS::assocTestSingle(
          gdsobj = mat,
          null.model = getNullModel(nullmodel),
          test = "Score"
        )
        ## double check
        if (!identical(res$variant.id, as.integer(1:ncol(mappedMatrix)))) {
          stop("", call. = FALSE)
        }
        P <- res$Score.pval
        effect <- res$Est
        effectSE <- res$Est.SE
        effectCIlower <- NA_real_
        effectCIlupper <- NA_real_
      }

      Pl[["mlm"]] <- P
      effectl[["mlm"]] <- effect
      effectSEl[["mlm"]] <- effectSE
      effectCIlowerl[["mlm"]] <- effectCIlower
      effectCIupperl[["mlm"]] <- effectCIupper
    }

    if ("fisher" %in% test) {
      res <- .geneSetAssoc_compute_fisher_stats(
        stats = stats,
        mappedMatrix = mappedMatrix,
        thresholds = threshold,
        oneSided = oneSided,
        ID = ID
      )
      results_list[["fisher"]] <- .geneSetAssoc_format_row(
        names = names(genesetlist),
        test = "lm",
        res = res,
        Ngenes = Ngenes,
        Ngenes_available = Ngenes_available,
        covar = covar
      )
    }

    results <- data.frame(
      geneSetName = c(
        rep(
          names(genesetlist),
          times = length(test[
            test %in% geneSetAssoc_tests_competitive_nothreshold
          ])
        ),
        rep(
          names(genesetlist),
          times = (length(test[
            test %in% geneSetAssoc_tests_competitive_threshold
          ]) *
            length(threshold))
        )
      ),
      test = c(
        rep(
          geneSetAssoc_tests_competitive_nothreshold[
            geneSetAssoc_tests_competitive_nothreshold %in% test
          ],
          each = length(genesetlist)
        ),
        rep(
          geneSetAssoc_tests_competitive_threshold[
            geneSetAssoc_tests_competitive_threshold %in% test
          ],
          each = (length(genesetlist) * length(threshold))
        )
      ),
      covar = c(
        rep(
          paste(covar, collapse = ","),
          times = (length(genesetlist) *
            length(test[test %in% geneSetAssoc_tests_competitive_nothreshold]))
        ),
        rep(
          NA_character_,
          times = (length(genesetlist) *
            length(threshold) *
            length(test[test %in% "fisher"]))
        )
      ),
      threshold = c(
        rep(
          rep(NA_real_, times = length(genesetlist)),
          times = length(test[
            test %in% geneSetAssoc_tests_competitive_nothreshold
          ])
        ),
        rep(
          rep(as.character(threshold), each = length(genesetlist)),
          times = length(test[
            test %in% geneSetAssoc_tests_competitive_threshold
          ])
        )
      ),
      geneSetSize = c(
        rep(
          lengths(genesetlist),
          times = length(test[
            test %in% geneSetAssoc_tests_competitive_nothreshold
          ])
        ),
        rep(
          lengths(genesetlist),
          times = (length(test[
            test %in% geneSetAssoc_tests_competitive_threshold
          ]) *
            length(threshold))
        )
      ),
      genesObs = c(
        rep(
          Ngenes_available,
          times = length(test[
            test %in% geneSetAssoc_tests_competitive_nothreshold
          ])
        ),
        rep(
          Ngenes_available,
          times = (length(test[
            test %in% geneSetAssoc_tests_competitive_threshold
          ]) *
            length(threshold))
        )
      ),
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    tst <- c(
      geneSetAssoc_tests_competitive_nothreshold,
      geneSetAssoc_tests_competitive_threshold
    )[
      c(
        geneSetAssoc_tests_competitive_nothreshold,
        geneSetAssoc_tests_competitive_threshold
      ) %in%
        test
    ]
    results_competitive <- cbind(
      results,
      effect = unlist(effectl[tst]),
      effectSE = unlist(effectSEl[tst]),
      effectCIlower = unlist(effectCIlowerl[tst]),
      effectCIupper = unlist(effectCIupperl[tst]),
      P = unlist(Pl[tst])
    )
    rownames(results_competitive) <- NULL
  }

  if (any(test %in% geneSetAssoc_tests_selfcontained)) {
    if ("ttest" %in% test) {
      P <- effect <- effectSE <- effectCIupper <- effectCIlower <- rep(
        NA_real_,
        length(genesetlist)
      )
      Z <- stats$Z
      for (i in seq_along(genesetlist)) {
        geneset <- names(genesetlist)[i]
        Z.tmp <- Z[as.matrix(mappedMatrix[, geneset, drop = FALSE])[, 1]]
        tryCatch(
          {
            if (oneSided) {
              fit <- t.test(Z.tmp, alternative = "greater")
            } else {
              fit <- t.test(Z.tmp, alternative = "two.sided")
            }

            effect[i] <- fit$estimate
            effectSE[i] <- fit$estimate / fit$statistic
            effectCIlower[i] <- fit$conf.int[1]
            effectCIupper[i] <- fit$conf.int[2]
            P[i] <- fit$p.value
          }
        )
      }
      Pl[["ttest"]] <- P
      effectl[["ttest"]] <- effect
      effectSEl[["ttest"]] <- effectSE
      effectCIlowerl[["ttest"]] <- effectCIlower
      effectCIupperl[["ttest"]] <- effectCIupper
    }

    if ("ztest" %in% test) {
      P <- effect <- effectSE <- effectCIupper <- effectCIlower <- rep(
        NA_real_,
        length(genesetlist)
      )
      Z <- stats$Z
      for (i in 1:length(genesetlist)) {
        geneset <- names(genesetlist)[i]
        Z.tmp <- Z[as.matrix(mappedMatrix[, geneset, drop = FALSE])[, 1]]
        tryCatch(
          {
            if (oneSided) {
              fit <- z_test(Z.tmp, alternative = "greater")
            } else {
              fit <- z_test(Z.tmp, alternative = "two.sided")
            }

            effect[i] <- fit$mean
            effectSE[i] <- fit$SE
            effectCIlower[i] <- fit$CIlower
            effectCIupper[i] <- fit$CIupper
            P[i] <- fit$P
          }
        )
      }
      Pl[["ztest"]] <- P
      effectl[["ztest"]] <- effect
      effectSEl[["ztest"]] <- effectSE
      effectCIlowerl[["ztest"]] <- effectCIlower
      effectCIupperl[["ztest"]] <- effectCIupper
    }

    if ("ACAT" %in% test) {
      P <- effect <- effectSE <- effectCIupper <- effectCIlower <- rep(
        NA_real_,
        length(genesetlist)
      )
      Pvals <- stats$P
      for (i in seq_along(genesetlist)) {
        geneset <- names(genesetlist)[i]
        P.tmp <- Pvals[as.matrix(mappedMatrix[, geneset, drop = FALSE])[, 1]]
        tryCatch(
          {
            Pval = .rvat_ACAT(P.tmp)
            P[i] <- Pval
          }
        )
      }
      Pl[["ACAT"]] <- P
      effectl[["ACAT"]] <- effect
      effectSEl[["ACAT"]] <- effectSE
      effectCIlowerl[["ACAT"]] <- effectCIlower
      effectCIupperl[["ACAT"]] <- effectCIupper
    }

    results <- data.frame(
      geneSetName = c(rep(
        names(genesetlist),
        times = length(test[test %in% geneSetAssoc_tests_selfcontained])
      )),
      test = c(rep(
        geneSetAssoc_tests_selfcontained[
          geneSetAssoc_tests_selfcontained %in% test
        ],
        each = length(genesetlist)
      )),
      covar = c(rep(
        NA_character_,
        times = (length(genesetlist) *
          length(test[test %in% geneSetAssoc_tests_selfcontained]))
      )),
      threshold = c(rep(
        rep(NA_real_, times = length(genesetlist)),
        times = length(test[test %in% geneSetAssoc_tests_selfcontained])
      )),
      geneSetSize = c(rep(
        lengths(genesetlist),
        times = length(test[test %in% geneSetAssoc_tests_selfcontained])
      )),
      genesObs = c(rep(
        Ngenes_available,
        times = length(test[test %in% geneSetAssoc_tests_selfcontained])
      )),
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    tst <- geneSetAssoc_tests_selfcontained[
      geneSetAssoc_tests_selfcontained %in% test
    ]
    results_selfcontained <- cbind(
      results,
      effect = unlist(effectl[tst]),
      effectSE = unlist(effectSEl[tst]),
      effectCIlower = unlist(effectCIlowerl[tst]),
      effectCIupper = unlist(effectCIupperl[tst]),
      P = unlist(Pl[tst])
    )
    rownames(results_selfcontained) <- NULL
  }

  if (
    any(test %in% geneSetAssoc_tests_selfcontained) &&
      any(test %in% geneSetAssoc_tests_competitive)
  ) {
    rbind(results_competitive, results_selfcontained)
  } else if (any(test %in% geneSetAssoc_tests_competitive)) {
    results_competitive
  } else if (any(test %in% geneSetAssoc_tests_selfcontained)) {
    results_selfcontained
  }
}

enrich_test <- function(
  stats,
  scorematrix,
  nullmodel = NULL,
  covar = NULL,
  test = c("lm", "mlm"),
  oneSided = TRUE,
  ID = "unit"
) {
  Pl <- effectl <- effectSEl <- effectCIlowerl <- effectCIupperl <- list()

  if ("lm" %in% test) {
    if (length(covar) > 0L) {
      X <- cbind(1, as.matrix(stats[, c(covar), drop = FALSE]))
    } else {
      X <- cbind(rep(1, nrow(stats)))
    }
    U1 <- crossprod(X, stats$Z)
    U2 <- solve(crossprod(X), U1)
    ytr <- stats$Z - X %*% U2
    U3 <- crossprod(X, scorematrix)
    U4 <- solve(crossprod(X), U3)
    Str <- scorematrix - X %*% U4
    effect <- as.vector(crossprod(ytr, Str)) / colSums(Str^2)
    Str2 <- colSums(Str^2)
    sig <- (sum(ytr^2) - effect^2 * Str2) / (nrow(stats) - ncol(X) - 2L)
    effectSE <- sqrt(sig * (1 / Str2))

    if (oneSided) {
      P <- pt(
        (effect / effectSE),
        nrow(stats) - ncol(X) - 1L,
        lower.tail = FALSE
      )
      fac <- qt(0.95, df = nrow(stats) - ncol(X) - 1L)
      effectCIlower <- effect - (fac * effectSE)
      effectCIupper <- rep(Inf, length(P))
    } else {
      P <- 2.0 *
        pt(
          abs(effect / effectSE),
          nrow(stats) - ncol(X) - 1L,
          lower.tail = FALSE
        )
      fac <- qt(0.975, df = nrow(stats) - ncol(X) - 1L)
      effectCIlower <- effect - (fac * effectSE)
      effectCIupper <- effect + (fac * effectSE)
    }

    Pl[["lm"]] <- P
    effectl[["lm"]] <- effect
    effectSEl[["lm"]] <- effectSE
    effectCIlowerl[["lm"]] <- effectCIlower
    effectCIupperl[["lm"]] <- effectCIupper
  }

  if ("mlm" %in% test) {
    # ..
  }

  results <- data.frame(
    name = rep(
      colnames(scorematrix),
      times = length(test[test %in% geneSetAssoc_tests_score])
    ),
    test = test[test %in% c("lm", "mlm")],
    covar = rep(
      paste(covar, collapse = ","),
      times = (ncol(scorematrix) *
        length(test[test %in% geneSetAssoc_tests_score]))
    ),
    genesObs = nrow(scorematrix),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  results <- cbind(
    results,
    effect = c(effectl[["lm"]], effectl[["mlm"]]),
    effectSE = c(effectSEl[["lm"]], effectSEl[["mlm"]]),
    effectCIlower = c(effectCIlowerl[["lm"]], effectCIlowerl[["mlm"]]),
    effectCIupper = c(effectCIupperl[["lm"]], effectCIupperl[["mlm"]]),
    P = c(Pl[["lm"]], Pl[["mlm"]])
  )
  results
}


.prepare_stats_GSA <- function(object, covar, Zcutoffs, INT, verbose = TRUE) {
  if (!is.null(Zcutoffs) && length(Zcutoffs) != 2L) {
    stop(
      "`Zcutoffs should be a vector of length 2 (minimum and maximum)",
      call. = FALSE
    )
  }

  if (!is.null(Zcutoffs) && INT) {
    warning(
      "Z-score cutoffs are specified while INT = TRUE. Z-score cutoffs won't be applied.",
      call. = FALSE
    )
  }

  if (!all(covar %in% colnames(object))) {
    stop(
      sprintf(
        "The following covariates are not present in the rvbResult: %s",
        paste(covar[!covar %in% colnames(object)], collapse = ",")
      ),
      call. = FALSE
    )
  }

  # Check if there are duplicate units
  if (anyDuplicated(object[["unit"]]) != 0L) {
    stop(
      "Units are duplicated; first filter the input so each row has a unique unit",
      call. = FALSE
    )
  }
  # Exclude rows with missing covariate values
  check <- complete.cases(as.data.frame(object)[, covar, drop = FALSE])
  if (sum(!check) > 0L) {
    if (verbose) {
      message(sprintf(
        "%s row(s) are excluded because of missing covariate values.",
        sum(!check)
      ))
    }
    object <- object[check, ]
  }

  # Exclude missing P-values
  if (any(is.na(object[["P"]]))) {
    if (verbose) {
      message(sprintf(
        "%s P-values are missing, these are excluded.",
        sum(is.na(object$P))
      ))
    }
    object <- object[!is.na(object$P), ]
  }

  # Add Z-score
  object$Z <- qnorm(1 - object$P)

  # Apply Z score cutoffs
  if (INT) {
    object$Z <- qnorm(
      (rank(object$Z, na.last = "keep") - 0.5) / sum(!is.na(object$Z))
    )
  } else if (!is.null(Zcutoffs)) {
    if (length(Zcutoffs) != 2L) {
      stop(
        "The length of `Zcutoffs` should be 2 (minimum and maximum).",
        call. = FALSE
      )
    }
    if (verbose) {
      message(sprintf(
        "%s Z-scores <%s are set to %s",
        sum(object$Z < Zcutoffs[1]),
        Zcutoffs[1],
        Zcutoffs[1]
      ))
    }
    if (verbose) {
      message(sprintf(
        "%s Z-scores >%s are set to %s",
        sum(object$Z > Zcutoffs[2]),
        Zcutoffs[2],
        Zcutoffs[2]
      ))
    }
    object$Z <- ifelse(object$Z < Zcutoffs[1], Zcutoffs[1], object$Z)
    object$Z <- ifelse(object$Z > Zcutoffs[2], Zcutoffs[2], object$Z)
  } else if (any(is.infinite(object$Z))) {
    # If Z-score cutoffs are not specified, check if any Z-scores are ±infinite
    if (any(is.infinite(object$Z) & object$Z < 0)) {
      minZ <- min(object$Z[!is.infinite(object$Z)])
      if (verbose) {
        message(sprintf(
          "%s Z-scores are -Inf, these are set to the minimum observed Z-score: %s.",
          sum(is.infinite(object$Z) & object$Z < 0),
          signif(minZ, 4)
        ))
      }
      object$Z[is.infinite(object$Z) & object$Z < 0] <- minZ
    }
    if (any(is.infinite(object$Z) & object$Z > 0)) {
      maxZ <- max(object$Z[!is.infinite(object$Z)])
      if (verbose) {
        message(sprintf(
          "%s Z-scores are +Inf, these are set to the maximum observed Z-score: %s.",
          sum(is.infinite(object$Z) & object$Z > 0),
          signif(maxZ, 4)
        ))
      }
      object$Z[is.infinite(object$Z) & object$Z > 0] <- maxZ
    }
  }
  object
}


z_test <- function(x, mu = 0, var = 1, alternative = "two.sided") {
  se <- var / sqrt(length(x))
  b <- mean(x)
  list(
    mean = b - mu,
    SE = se,
    P = if (alternative == "two.sided") {
      2.0 * pnorm(abs((mean(x) - mu) / se), lower.tail = FALSE)
    } else {
      pnorm(((mean(x) - mu) / se), lower.tail = FALSE)
    },
    CIlower = if (alternative == "two.sided") {
      b - qnorm(0.975) * se
    } else {
      b - qnorm(0.95) * se
    },
    CIupper = if (alternative == "two.sided") b + qnorm(0.975) * se else Inf
  )
}

.geneSetAssoc_compute_lm_stats <- function(stats, mat, covar, oneSided) {
  # based on: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-166
  if (length(covar) > 0L) {
    X <- cbind(1, as.matrix(stats[, covar, drop = FALSE]))
  } else {
    X <- cbind(rep(1, nrow(stats)))
  }
  U1 <- crossprod(X, stats$Z)
  U2 <- solve(crossprod(X), U1)
  ytr <- stats$Z - X %*% U2
  U3 <- crossprod(X, mat)
  U4 <- solve(crossprod(X), U3)
  Str <- mat - X %*% U4
  effect <- as.vector(crossprod(ytr, Str)) / colSums(Str^2)
  Str2 <- colSums(Str^2)
  sig <- (sum(ytr^2) - effect^2 * Str2) / (nrow(stats) - ncol(X) - 2L)
  effectSE <- sqrt(sig * (1 / Str2))

  if (oneSided) {
    P <- pt(
      (effect / effectSE),
      nrow(stats) - ncol(X) - 1L,
      lower.tail = FALSE
    )
    fac <- qt(0.95, df = nrow(stats) - ncol(X) - 1L)
    effectCIlower <- effect - (fac * effectSE)
    effectCIupper <- rep(Inf, length(P))
  } else {
    P <- 2.0 *
      pt(
        abs(effect / effectSE),
        nrow(stats) - ncol(X) - 1L,
        lower.tail = FALSE
      )
    fac <- qt(0.975, df = nrow(stats) - ncol(X) - 1L)
    effectCIlower <- effect - (fac * effectSE)
    effectCIupper <- effect + (fac * effectSE)
  }

  list(
    effect = effect,
    effectSE = effectSE,
    effectCIlower = effectCIlower,
    effectCIupper = effectCIupper,
    P = P
  )
}

.geneSetAssoc_compute_fisher_stats <- function(
  stats,
  mappedMatrix,
  thresholds,
  oneSided,
  ID
) {
  out_list <- list()
  set_sizes <- Matrix::colSums(mappedMatrix)
  n_units <- nrow(stats)
  for (thresh in thresholds) {
    is_sig <- as.integer(stats$P < thresh)
    n_sig <- sum(is_sig)
    n_not_sig <- n_units - n_sig

    sig_in_set <- as.vector(Matrix::crossprod(mappedMatrix, is_sig))
    sig_not_in_set <- n_sig - sig_in_set
    not_sig_in_set <- set_sizes - sig_in_set
    not_sig_not_in_set <- n_not_sig - sig_not_in_set

    results <- vapply(
      seq_along(sig_in_set),
      function(i) {
        mat <- matrix(
          c(
            sig_in_set[i],
            sig_not_in_set[i],
            not_sig_in_set[i],
            not_sig_not_in_set[i]
          ),
          nrow = 2,
          byrow = TRUE
        )

        tryCatch(
          {
            tst <- fisher.test(
              mat,
              alternative = if (oneSided) "greater" else "two.sided"
            )
            c(
              effect = unname(tst$estimate),
              effectCIlower = unname(tst$conf.int[1]),
              effectCIupper = unname(tst$conf.int[2]),
              P = unname(tst$p.value)
            )
          },
          error = function(e) {
            c(
              effect = NA_real_,
              effectCIlower = NA_real_,
              effectCIupper = NA_real_,
              P = NA_real_
            )
          }
        )
      },
      FUN.VALUE = numeric(4)
    )
    results <- data.frame(
      threshold = thresh,
      effect = results["effect", ],
      effectSE = rep(NA_real_, ncol(results)),
      effectCIlower = results["effectCIlower", ],
      effectCIupper = results["effectCIupper", ],
      P = results["P", ],
      stringsAsFactors = FALSE
    )

    out_list[[as.character(thresh)]] <- results
  }
  do.call(rbind, out_list)
}

.geneSetAssoc_format_row <- function(
  names,
  test,
  res,
  Ngenes,
  Ngenes_available,
  covar
) {
  results <- data.frame(
    geneSetName = names,
    test = test,
    covar = if (is.na(covar)) {
      NA_character_
    } else {
      paste(covar, collapse = ",")
    },
    threshold = if ("threshold" %in% colnames(res)) {
      res$threshold
    } else {
      NA_real_
    },
    geneSetSize = Ngenes,
    genesObs = Ngenes_available,
    effect = res$effect,
    effectSE = res$effectSE,
    effectCIlower = res$effectCIlower,
    effectCIupper = res$effectCIupper,
    P = res$P,
    stringsAsFactors = FALSE
  )
  rownames(results) <- NULL
  results
}


.geneSetAssoc_compute_selfcontained <- function(
  stats,
  mappedMatrix,
  test_selfcontained,
  oneSided
) {
  Z <- stats$Z
  Pvals <- stats$P

  # Iterate over columns (genesets)
  # mappedMatrix is sparse, so we assume iterating cols is relatively efficient
  # But direct subsetting `mappedMatrix[, i]` can be slow if not careful.
  # We use the fact that it is sparse to get indices efficiently.

  results <- vapply(
    seq_len(ncol(mappedMatrix)),
    function(i) {
      # Get indices of units in this set
      # For a sparse matrix, `which` or subsetting returns the indices
      mappedMatrix_i <- mappedMatrix[, i, drop = TRUE]

      subset_Z <- Z[mappedMatrix_i]
      subset_P <- Pvals[mappedMatrix_i]

      out_template <- c(
        effect = NA_real_,
        effectSE = NA_real_,
        effectCIlower = NA_real_,
        effectCIupper = NA_real_,
        P = NA_real_
      )

      result_list <- list()
      if (test_selfcontained == "ttest") {
        out_ttest <- out_template
        tryCatch(
          {
            fit <- t.test(
              subset_Z,
              alternative = if (oneSided) "greater" else "two.sided"
            )
            out_ttest["effect"] <- unname(fit$estimate)
            out_ttest["effectSE"] <- unname(fit$estimate / fit$statistic)
            out_ttest["effectCIlower"] <- fit$conf.int[1]
            out_ttest["effectCIupper"] <- fit$conf.int[2]
            out_ttest["P"] <- fit$p.value
          }
        )
        result_list[["ttest"]] <- out_ttest
      } else if (test_selfcontained == "ztest") {
        out_ztest <- out_template
        tryCatch(
          {
            fit <- z_test(
              subset_Z,
              alternative = if (oneSided) "greater" else "two.sided"
            )
            out_ztest["effect"] <- fit$mean
            out_ztest["effectSE"] <- fit$SE
            out_ztest["effectCIlower"] <- fit$CIlower
            out_ztest["effectCIupper"] <- fit$CIupper
            out_ztest["P"] <- fit$P
          }
        )
        result_list[["ztest"]] <- out_ztest
      } else if (test_type == "ACAT") {
        out_acat <- out_template
        tryCatch(
          {
            out_acat["P"] <- .rvat_ACAT(subset_P)
          }
        )
        result_list[["ACAT"]] <- out_acat
      }

      result_list
    },
    FUN.VALUE = numeric(5)
  )

  list(
    effect = vapply(res, function(x) x$eff, numeric(1)),
    effectSE = vapply(res, function(x) x$se, numeric(1)),
    effectCIlower = vapply(res, function(x) x$low, numeric(1)),
    effectCIupper = vapply(res, function(x) x$up, numeric(1)),
    P = vapply(res, function(x) x$p, numeric(1))
  )
}
