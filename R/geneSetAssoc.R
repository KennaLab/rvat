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
    ## only 'lm' is implemented for `scoreMatrix` object
    if (!is.null(scoreMatrix) && !is.null(geneSet)) {
      stop(
        "Specify either `geneSet` or `scoreMatrix`, not both.",
        call. = FALSE
      )
    }

    if (!is.null(scoreMatrix)) {
      test <- test[test %in% geneSetAssoc_tests_score]
    } else if (!is.null(geneSet)) {
      if (!is.null(condition)) {
        test <- test[test %in% geneSetAssoc_tests_competitive_condition]
      } else {
        test <- test[test %in% geneSetAssoc_tests]
      }
    } else {
      stop(
        "Either `geneSet` or `scoreMatrix` must be specified.",
        call. = FALSE
      )
    }

    if (length(test) == 0L) {
      stop("No valid test specified.", call. = FALSE)
    }

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
          "The mlm test requires specifying the 'cormatrix' argument, ",
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
    } else {
      # prepare test-stats (Z-scores, covariates, etc.)
      object <- .prepare_stats_GSA(
        object,
        covar,
        Zcutoffs,
        INT,
        verbose = verbose
      )
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

    # run analyses based on input type
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
        out_con <- gzfile(output, "w")
        on.exit(close(out_con), add = TRUE)
        write.table(
          result_scorematrix,
          sep = "\t",
          quote = FALSE,
          file = out_con,
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

      # format results
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

      # write or return
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
  # scoreMatrix should be of type matrix/Matrix
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
  # handle thresholds (defaults to bonferroni if not specified)
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
  nullmodel = NULL,
  covar = NULL,
  test = "lm",
  threshold = NULL,
  oneSided = TRUE,
  ID = "unit",
  memlimit = 1000L
) {
  if (condition.type == "geneSet") {
    chunksCond <- split(
      seq_along(condition),
      ceiling(seq_along(condition) / memlimit)
    )
    condNames <- names(condition)
  } else if (condition.type == "matrix") {
    chunksCond <- split(
      seq_len(ncol(condition)),
      ceiling(seq_len(ncol(condition)) / memlimit)
    )
    condNames <- colnames(condition)
  }

  results_list <- lapply(chunksCond, function(chunk) {
    # retrieve condition matrix for current chunk
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
    } else if (condition.type == "matrix") {
      mappedMatrixCond <- condition[, chunk, drop = FALSE]
    }

    # iterate over each condition in this chunk
    cond_results <- lapply(colnames(mappedMatrixCond), function(cond_name) {
      stats <- as.data.frame(object)
      stats[[cond_name]] <- as.matrix(mappedMatrixCond)[, cond_name]

      # handle covariates
      current_covar <- if (length(covar) == 0L) {
        cond_name
      } else {
        c(covar, cond_name)
      }
      if (var(stats[[cond_name]], na.rm = TRUE) == 0) {
        current_covar <- current_covar[current_covar != cond_name]
      }

      # compute lm stats
      res <- .geneSetAssoc_compute_lm_stats(
        stats = stats,
        mat = mappedMatrix,
        covar = current_covar,
        oneSided = oneSided
      )

      # format output
      res_formatted <- .geneSetAssoc_format_row(
        names = names(genesetlist),
        test = "lm",
        res = res,
        Ngenes = lengths(genesetlist),
        Ngenes_available = Matrix::colSums(mappedMatrix),
        covar = covar
      )
      res_formatted$condition <- cond_name
      res_formatted
    })
    do.call(rbind, cond_results)
  })
  # merge chunks
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL

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
  results_list <- list()
  Ngenes_available <- Matrix::colSums(mappedMatrix)
  Ngenes <- lengths(genesetlist)

  if ("lm" %in% test) {
    res <- .geneSetAssoc_compute_lm_stats(
      stats,
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
        snpID = as.integer(seq_len(ncol(mappedMatrix))),
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
      if (!identical(res$variant.id, as.integer(seq_len(ncol(mappedMatrix))))) {
        stop("MLM results do not align with gene sets.", call. = FALSE)
      }
      P <- res$Score.pval
      effect <- res$Est
      effectSE <- res$Est.SE
      effectCIlower <- NA_real_
      effectCIupper <- NA_real_
    }

    results_list[["mlm"]] <- .geneSetAssoc_format_row(
      names = names(genesetlist),
      test = "mlm",
      res = list(
        P = P,
        effect = effect,
        effectSE = effectSE,
        effectCIlower = effectCIlower,
        effectCIupper = effectCIupper
      ),
      Ngenes = Ngenes,
      Ngenes_available = Ngenes_available,
      covar = covar
    )
  }

  if ("fisher" %in% test) {
    results_fisher <- lapply(
      threshold,
      FUN = function(thresh) {
        res <- .geneSetAssoc_compute_fisher_stats(
          stats = stats,
          mappedMatrix = mappedMatrix,
          threshold = thresh,
          oneSided = oneSided,
          ID = ID
        )
        res_formatted <- .geneSetAssoc_format_row(
          names = names(genesetlist),
          test = "fisher",
          res = res,
          Ngenes = Ngenes,
          Ngenes_available = Ngenes_available,
          covar = NA_character_
        )
      }
    )
    results_list[["fisher"]] <- do.call(rbind, results_fisher)
  }

  if (any(test %in% geneSetAssoc_tests_selfcontained)) {
    results_selfcontained <- .geneSetAssoc_compute_selfcontained(
      stats = stats,
      mappedMatrix = mappedMatrix,
      test_selfcontained = test[
        test %in% geneSetAssoc_tests_selfcontained
      ],
      oneSided = oneSided
    )
    for (tst in names(results_selfcontained)) {
      results_list[[tst]] <- .geneSetAssoc_format_row(
        names = names(genesetlist),
        test = tst,
        res = results_selfcontained[[tst]],
        Ngenes = Ngenes,
        Ngenes_available = Ngenes_available,
        covar = NA_character_
      )
    }
  }

  # merge results and return
  do.call(rbind, results_list)
}

enrich_test <- function(
  stats,
  scorematrix,
  nullmodel = NULL,
  covar = NULL,
  test = "lm",
  oneSided = TRUE,
  ID = "unit"
) {
  results_list <- list()
  Ngenes_available <- nrow(scorematrix)
  Ngenes <- rep(NA_integer_, ncol(scorematrix))

  if ("lm" %in% test) {
    res <- .geneSetAssoc_compute_lm_stats(
      stats,
      mat = scorematrix,
      covar = covar,
      oneSided = oneSided
    )

    res <- data.frame(
      name = colnames(scorematrix),
      test = "lm",
      covar = if (length(covar) == 0L) {
        NA_character_
      } else {
        paste(covar, collapse = ",")
      },
      genesObs = nrow(scorematrix),
      effect = res$effect,
      effectSE = res$effectSE,
      effectCIlower = res$effectCIlower,
      effectCIupper = res$effectCIupper,
      P = res$P,
      stringsAsFactors = FALSE
    )
    results_list[["lm"]] <- res
  }

  do.call(rbind, results_list)
}


.prepare_stats_GSA <- function(object, covar, Zcutoffs, INT, verbose = TRUE) {
  if (!is.null(Zcutoffs) && length(Zcutoffs) != 2L) {
    stop(
      "`Zcutoffs` should be a vector of length 2 (minimum and maximum)",
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

  # check if there are duplicate units
  if (anyDuplicated(object[["unit"]]) != 0L) {
    stop(
      "Units are duplicated; first filter the input so each row has a unique unit",
      call. = FALSE
    )
  }
  # exclude rows with missing covariate values
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

  # exclude missing P-values
  if (any(is.na(object[["P"]]))) {
    if (verbose) {
      message(sprintf(
        "%s P-values are missing, these are excluded.",
        sum(is.na(object$P))
      ))
    }
    object <- object[!is.na(object$P), ]
  }

  # add Z-score
  object$Z <- qnorm(1 - object$P)

  # apply Z score cutoffs
  if (INT) {
    object$Z <- qnorm(
      (rank(object$Z, na.last = "keep") - 0.5) / sum(!is.na(object$Z))
    )
  } else if (!is.null(Zcutoffs)) {
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

z_test <- function(x, mu = 0, sd = 1, alternative = "two.sided") {
  se <- sd / sqrt(length(x))
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
  threshold,
  oneSided,
  ID
) {
  out_list <- list()
  set_sizes <- Matrix::colSums(mappedMatrix)
  n_units <- nrow(stats)

  is_sig <- as.integer(stats$P < threshold)
  n_sig <- sum(is_sig)
  n_not_sig <- n_units - n_sig

  sig_in_set <- as.vector(Matrix::crossprod(mappedMatrix, is_sig))
  sig_not_in_set <- n_sig - sig_in_set
  not_sig_in_set <- set_sizes - sig_in_set
  not_sig_not_in_set <- n_not_sig - not_sig_in_set

  results <- vapply(
    seq_along(sig_in_set),
    function(i) {
      mat <- round(matrix(
        c(
          sig_in_set[i],
          sig_not_in_set[i],
          not_sig_in_set[i],
          not_sig_not_in_set[i]
        ),
        nrow = 2,
        byrow = TRUE
      ))

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
    threshold = threshold,
    effect = results["effect", ],
    effectSE = rep(NA_real_, ncol(results)),
    effectCIlower = results["effectCIlower", ],
    effectCIupper = results["effectCIupper", ],
    P = results["P", ],
    stringsAsFactors = FALSE
  )
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
  n_sets <- ncol(mappedMatrix)

  # initialize output
  results <- lapply(test_selfcontained, function(x) {
    list(
      effect = rep(NA_real_, n_sets),
      effectSE = rep(NA_real_, n_sets),
      effectCIlower = rep(NA_real_, n_sets),
      effectCIupper = rep(NA_real_, n_sets),
      P = rep(NA_real_, n_sets)
    )
  })
  names(results) <- test_selfcontained

  for (i in seq_len(n_sets)) {
    mappedMatrix_i <- mappedMatrix[, i, drop = TRUE]

    # skip if empty
    if (!any(mappedMatrix_i)) {
      next
    }

    # Z-scores for current set
    subset_Z <- Z[mappedMatrix_i]

    # loop through tests
    for (tst in test_selfcontained) {
      tryCatch(
        {
          if (tst == "ttest") {
            fit <- t.test(
              subset_Z,
              alternative = if (oneSided) "greater" else "two.sided"
            )

            results[[tst]]$effect[i] <- fit$estimate
            results[[tst]]$effectSE[i] <- fit$estimate / fit$statistic
            results[[tst]]$effectCIlower[i] <- fit$conf.int[1]
            results[[tst]]$effectCIupper[i] <- fit$conf.int[2]
            results[[tst]]$P[i] <- fit$p.value
          } else if (tst == "ztest") {
            fit <- z_test(
              subset_Z,
              alternative = if (oneSided) "greater" else "two.sided"
            )

            results[[tst]]$effect[i] <- fit$mean
            results[[tst]]$effectSE[i] <- fit$SE
            results[[tst]]$effectCIlower[i] <- fit$CIlower
            results[[tst]]$effectCIupper[i] <- fit$CIupper
            results[[tst]]$P[i] <- fit$P
          } else if (tst == "ACAT") {
            subset_P <- Pvals[mappedMatrix_i]
            results[[tst]]$P[i] <- .rvat_ACAT(subset_P)
          }
        },
        error = function(e) {}
      )
    }
  }

  results
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
    covar = if (is.null(covar) || all(is.na(covar))) {
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
