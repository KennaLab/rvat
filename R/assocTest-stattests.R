# rvb tests -------------------------------------------------------------------
.rvb_tests_rvb <- function(
  GT,
  results,
  test,
  pheno,
  model,
  model.nbinom,
  null,
  covar,
  continuous,
  maxitFirth = 1000L,
  nResampling = 1000L,
  returnDF = FALSE
) {
  P <- OR <- effect <- effectSE <- effectCIupper <- effectCIlower <-
    rep(NA_real_, length(test))
  names(P) <- names(OR) <- names(effect) <- names(effectSE) <- names(
    effectCIupper
  ) <- names(effectCIlower) <-
    test

  ## skip non-sensible tests
  if (sum(colData(GT)$aggregate > 0) < 2) {
    warning(
      "Less than two samples have a non-zero burden score, skipping tests.",
      call. = FALSE
    )
    test <- c()
  }

  if (continuous) {
    out_type <- "C"
  } else {
    out_type <- "D"
    if (
      sum(colData(GT)[, pheno] == 1) < 2 || sum(colData(GT)[, pheno] == 0) < 2
    ) {
      warning(
        "Fewer than two cases or controls, skipping tests.",
        call. = FALSE
      )
      test <- c()
    }
  }

  # Burden test (continuous)
  if ("lm" %in% test) {
    tryCatch(
      {
        fit <- lm(model, data = colData(GT))
        effect["lm"] <- fit$coefficients["aggregate"]
        effectSE["lm"] <- summary(fit)$coef["aggregate", 2]
        effectCIlower["lm"] <- confint(fit)["aggregate", 1]
        effectCIupper["lm"] <- confint(fit)["aggregate", 2]
        P["lm"] <- summary(fit)$coef["aggregate", 4]
      },
      error = function(e) {
        message(sprintf("Failed test '%s'\n%s", "lm", e))
      }
    )
  }

  # Burden test (binary)
  if ("firth" %in% test) {
    tryCatch(
      {
        fit <- logistf::logistf(
          model,
          data = colData(GT),
          plconf = (which(c(covar, "aggregate") == "aggregate") + 1),
          control = logistf::logistf.control(maxit = maxitFirth),
          plcontrol = logistf::logistpl.control(maxit = maxitFirth)
        )

        ## check if converged, if not, set to NA
        if (.check_conv_firth(fit, maxit = maxitFirth)) {
          effect["firth"] <- fit$coefficients["aggregate"]
          OR["firth"] <- exp(fit$coefficients["aggregate"])
          effectSE["firth"] <- sqrt(diag(vcov(fit)))["aggregate"]
          effectCIlower["firth"] <- fit$ci.lower["aggregate"]
          effectCIupper["firth"] <- fit$ci.upper["aggregate"]
          P["firth"] <- fit$prob["aggregate"]
        } else {
          effect["firth"] <- OR["firth"] <- effectSE["firth"] <- effectCIlower[
            "firth"
          ] <- effectCIupper["firth"] <- P["firth"] <- NA
        }
      },
      error = function(e) {
        message(sprintf("Failed test '%s'\n%s", "firth", e))
      }
    )
  }

  if ("glm" %in% test) {
    tryCatch(
      {
        fit <- glm(model, data = colData(GT), family = "binomial")
        effect["glm"] <- summary(fit)$coef["aggregate", 1]
        OR["glm"] <- exp(summary(fit)$coef["aggregate", 1])
        effectSE["glm"] <- summary(fit)$coef["aggregate", 2]
        effectCIlower["glm"] <- confint.default(fit)["aggregate", 1]
        effectCIupper["glm"] <- confint.default(fit)["aggregate", 2]
        P["glm"] <- summary(fit)$coef["aggregate", 4]
      },
      error = function(e) {
        message(sprintf("Failed test '%s'\n%s", "glm", e))
      }
    )
  }

  if ("scoreSPA" %in% test) {
    tryCatch(
      {
        score.null <- SPAtest::ScoreTest_wSaddleApprox_NULL_Model(
          null,
          data = as.data.frame(colData(GT))
        )
        fit <- SPAtest::ScoreTest_SPA(
          genos = colData(GT)[["aggregate"]],
          obj.null = score.null,
          minmac = 0,
          beta.out = FALSE
        )

        P["scoreSPA"] <- fit$p.value
      },
      error = function(e) {
        message(sprintf("Failed test '%s'\n%s", "scoreSPA", e))
      }
    )
  }

  if ("nbinom" %in% test) {
    tryCatch(
      {
        fit <- MASS::glm.nb(model.nbinom, data = colData(GT))
        effect["nbinom"] <- summary(fit)$coefficients[pheno, 1]
        effectSE["nbinom"] <- summary(fit)$coefficients[pheno, 2]
        effectCIlower["nbinom"] <- effect["nbinom"] -
          (1.96 * summary(fit)$coefficients[pheno, 2])
        effectCIupper["nbinom"] <- effect["nbinom"] +
          (1.96 * summary(fit)$coefficients[pheno, 2])
        P["nbinom"] <- summary(fit)$coefficients[pheno, 4]
      },
      error = function(e) {
        message(sprintf("Failed test '%s'\n%s", "nbinom", e))
      }
    )
  }

  # SKAT analyses
  if (
    sum(
      c(
        "skat_burden",
        "skat",
        "skato",
        "skat_burden_robust",
        "skat_robust",
        "skato_robust"
      ) %in%
        test
    ) >
      0
  ) {
    skat.null <- SKAT::SKAT_Null_Model(
      null,
      data = colData(GT),
      out_type = out_type
    )

    if ("skat_burden" %in% test) {
      tryCatch(
        {
          if (continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null,
              r.corr = 1,
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          } else {
            fit <- SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "Burden",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }

          P["skat_burden"] <- fit$p.value
        },
        error = function(e) {
          message(sprintf("Failed test '%s'\n%s", "skat_burden", e))
        }
      )
    }

    if ("skat" %in% test) {
      tryCatch(
        {
          if (continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null,
              r.corr = 0,
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          } else {
            fit <- SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "SKAT",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }

          P["skat"] <- fit$p.value
        },
        error = function(e) {
          message(sprintf("Failed test '%s'\n%s", "skat", e))
        }
      )
    }

    if ("skato" %in% test) {
      tryCatch(
        {
          if (continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          } else {
            fit <- SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }

          P["skato"] <- fit$p.value
        },
        error = function(e) {
          message(sprintf("Failed test '%s'\n%s", "skato", e))
        }
      )
    }

    if (
      sum(c("skat_burden_robust", "skat_robust", "skato_robust") %in% test) > 0
    ) {
      varkeep <- Matrix::rowSums(assays(GT)$GT) > 0
    }

    if ("skat_burden_robust" %in% test) {
      if (sum(varkeep) >= 1) {
        tryCatch(
          {
            fit <- SKAT::SKATBinary_Robust(
              t(assays(GT)$GT)[, varkeep, drop = FALSE],
              skat.null,
              method = "Burden",
              impute.method = "fixed",
              weights = rowData(GT)$w[varkeep],
              missing_cutoff = 1
            )
            P["skat_burden_robust"] <- fit$p.value
          },
          error = function(e) {
            message(sprintf("Failed test '%s'\n%s", "skat_burden_robust", e))
          }
        )
      }
    }

    if ("skat_robust" %in% test) {
      if (sum(varkeep) >= 1) {
        tryCatch(
          {
            fit <- SKAT::SKATBinary_Robust(
              t(assays(GT)$GT)[, varkeep, drop = FALSE],
              skat.null,
              method = "SKAT",
              impute.method = "fixed",
              weights = rowData(GT)$w[varkeep],
              missing_cutoff = 1
            )
            P["skat_robust"] <- fit$p.value
          },
          error = function(e) {
            message(sprintf("Failed test '%s'\n%s", "skat_robust", e))
          }
        )
      }
    }

    if ("skato_robust" %in% test) {
      if (sum(varkeep) >= 1) {
        tryCatch(
          {
            fit <- SKAT::SKATBinary_Robust(
              t(assays(GT)$GT)[, varkeep, drop = FALSE],
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w[varkeep],
              missing_cutoff = 1
            )
            P["skato_robust"] <- fit$p.value
          },
          error = function(e) {
            message(sprintf("Failed test '%s'\n%s", "skato_robust", e))
          }
        )
      }
    }
  }

  if (sum(c("acatv", "acatvSPA", "acatvfirth") %in% test) > 0) {
    if ("acatv" %in% test) {
      tryCatch(
        {
          acat_null <- .acat_NULL_Model(
            colData(GT)[[pheno]],
            Z = if (!is.null(covar)) {
              as.matrix(colData(GT)[, covar, drop = FALSE])
            } else {
              NULL
            }
          )

          # fit
          mac <- Matrix::rowSums(assays(GT)$GT)
          maf <- Matrix::rowSums(assays(GT)$GT) / (2 * ncol(GT))
          fit <- .acatv_rvat(
            t(assays(GT)$GT),
            obj = acat_null,
            method = "original",
            weights = rowData(GT)$w,
            mac.thresh = 10,
            maf = maf,
            mac = mac
          )

          P["acatv"] <- fit
        },
        error = function(e) {
          message(sprintf("Failed test '%s'\n%s", "acatv", e))
        }
      )
    }

    if (any(c("acatvSPA", "acatvfirth") %in% test)) {
      mac <- Matrix::rowSums(assays(GT)$GT)
      maf <- Matrix::rowSums(assays(GT)$GT) / (2 * ncol(GT))
    }

    if (any(c("acatvSPA") %in% test)) {
      tryCatch(
        {
          # null model
          obj <- SPAtest::ScoreTest_wSaddleApprox_NULL_Model(
            null,
            data = as.data.frame(colData(GT))
          )

          # fit
          fit <- .acatv_rvat(
            t(assays(GT)$GT),
            obj = obj,
            method = "scoreSPA",
            weights = rowData(GT)$w,
            mac.thresh = 10,
            maf = maf,
            mac = mac
          )

          P["acatvSPA"] <- fit
        },
        error = function(e) {
          message(sprintf("Failed test '%s'\n%s", "acatvSPA", e))
        }
      )
    }

    if ("acatvfirth" %in% test) {
      tryCatch(
        {
          # fit
          fit <- .acatv_rvat(
            t(assays(GT)$GT),
            obj = NULL,
            data = colData(GT),
            covar = covar,
            model = model,
            method = "firth",
            weights = rowData(GT)$w,
            mac.thresh = 10,
            maf = maf,
            mac = mac,
            maxitFirth = maxitFirth
          )

          P["acatvfirth"] <- fit
        },
        error = function(e) {
          message(sprintf("Failed test '%s'\n%s", "acatvfirth", e))
        }
      )
    }
  }

  results$effect <- effect
  results$effectSE <- effectSE
  results$effectCIlower <- effectCIlower
  results$effectCIupper <- effectCIupper
  results$OR <- OR
  results$P <- P
  results
}

# Singlevar tests -------------------------------------------------------------

.rvb_tests_singlevar <- function(
  GT,
  results,
  pheno,
  test,
  model,
  model.nbinom,
  null,
  covar,
  continuous,
  returnDF = FALSE,
  maxitFirth = 1000L,
  verbose = TRUE
) {
  testable <- which(
    Matrix::rowSums(
      results[, c("caseMAC", "ctrlMAC")][!duplicated(results$VAR_id), ],
      na.rm = TRUE
    ) >=
      2
  )
  if (length(testable) < nrow(GT)) {
    if (verbose) {
      warning(
        sprintf(
          "%s/%s variants have less than 2 carriers, tests will be skipped for these variants.",
          nrow(GT) - length(testable),
          nrow(GT)
        ),
        call. = FALSE
      )
    }
  }
  Pl <- ORl <- effectl <- effectSEl <- effectCIlowerl <- effectCIupperl <- list()

  # linear model (continuous pheno)
  if ("lm" %in% test) {
    P <- OR <- effect <- effectSE <- effectCIupper <- effectCIlower <- rep(
      NA_real_,
      nrow(GT)
    )
    for (i in testable) {
      tryCatch(
        {
          fit <- lm(
            model,
            data = cbind(
              aggregate = assays(GT)$GT[i, ],
              data.frame(colData(GT))
            )
          )
          effect[i] <- fit$coefficients["aggregate"]
          effectSE[i] <- summary(fit)$coef["aggregate", 2]
          effectCIlower[i] <- confint(fit)["aggregate", 1]
          effectCIupper[i] <- confint(fit)["aggregate", 2]
          P[i] <- summary(fit)$coef["aggregate", 4]
        },
        error = function(e) {
          message(sprintf(
            "Failed test '%s' for variant '%s'\n%s",
            "lm",
            rownames(GT)[i],
            e
          ))
        }
      )
    }
    Pl[["lm"]] <- P
    ORl[["lm"]] <- OR
    effectl[["lm"]] <- effect
    effectSEl[["lm"]] <- effectSE
    effectCIlowerl[["lm"]] <- effectCIlower
    effectCIupperl[["lm"]] <- effectCIupper
  }

  # Firth logistic regression
  if ("firth" %in% test) {
    P <- OR <- effect <- effectSE <- effectCIupper <- effectCIlower <- rep(
      NA_real_,
      nrow(GT)
    )
    for (i in testable) {
      tryCatch(
        {
          fit <- logistf::logistf(
            model,
            data = cbind(
              aggregate = assays(GT)$GT[i, ],
              data.frame(colData(GT))
            ),
            plconf = (which(c(covar, "aggregate") == "aggregate") + 1),
            control = logistf::logistf.control(maxit = maxitFirth),
            plcontrol = logistf::logistpl.control(maxit = maxitFirth)
          )
          if (.check_conv_firth(fit, maxit = maxitFirth)) {
            effect[i] <- fit$coefficients["aggregate"]
            effectSE[i] <- sqrt(diag(vcov(fit)))["aggregate"]
            effectCIlower[i] <- fit$ci.lower["aggregate"]
            effectCIupper[i] <- fit$ci.upper["aggregate"]
            P[i] <- fit$prob["aggregate"]
          }
        },
        error = function(e) {
          message(sprintf(
            "Failed test '%s' for variant '%s'\n%s",
            "firth",
            rownames(GT)[i],
            e
          ))
        }
      )
    }
    Pl[["firth"]] <- P
    ORl[["firth"]] <- exp(effect)
    effectl[["firth"]] <- effect
    effectSEl[["firth"]] <- effectSE
    effectCIlowerl[["firth"]] <- effectCIlower
    effectCIupperl[["firth"]] <- effectCIupper
  }

  # glm logistic regression
  if ("glm" %in% test) {
    P <- OR <- effect <- effectSE <- effectCIupper <- effectCIlower <- rep(
      NA_real_,
      nrow(GT)
    )
    for (i in testable) {
      tryCatch(
        {
          fit <- glm(
            model,
            data = cbind(
              aggregate = assays(GT)$GT[i, ],
              data.frame(colData(GT))
            ),
            family = binomial(link = "logit")
          )
          effect[i] <- summary(fit)$coef["aggregate", 1]
          effectSE[i] <- summary(fit)$coef["aggregate", 2]
          effectCIlower[i] <- confint.default(fit)["aggregate", 1]
          effectCIupper[i] <- confint.default(fit)["aggregate", 2]
          P[i] <- summary(fit)$coef["aggregate", 4]
        },
        error = function(e) {
          message(sprintf(
            "Failed test '%s' for variant '%s'\n%s",
            "glm",
            rownames(GT)[i],
            e
          ))
        }
      )
    }
    Pl[["glm"]] <- P
    ORl[["glm"]] <- exp(effect)
    effectl[["glm"]] <- effect
    effectSEl[["glm"]] <- effectSE
    effectCIlowerl[["glm"]] <- effectCIlower
    effectCIupperl[["glm"]] <- effectCIupper
  }
  if ("nbinom" %in% test) {
    P <- OR <- effect <- effectSE <- effectCIupper <- effectCIlower <- rep(
      NA_real_,
      nrow(GT)
    )
    for (i in testable) {
      tryCatch(
        {
          fit <- MASS::glm.nb(
            model.nbinom,
            data = cbind(
              aggregate = assays(GT)$GT[i, ],
              data.frame(colData(GT))
            )
          )
          effect[i] <- summary(fit)$coefficients[pheno, 1]
          effectSE[i] <- summary(fit)$coefficients[pheno, 2]
          effectCIlower[i] <- effect[i] -
            (1.96 * summary(fit)$coefficients[pheno, 2])
          effectCIupper[i] <- effect[i] +
            (1.96 * summary(fit)$coefficients[pheno, 2])
          P[i] <- summary(fit)$coefficients[pheno, 4]
        },
        error = function(e) {
          message(sprintf(
            "Failed test '%s' for variant '%s'\n%s",
            "nbinom",
            rownames(GT)[i],
            e
          ))
        }
      )
    }
    Pl[["nbinom"]] <- P
    ORl[["nbinom"]] <- OR
    effectl[["nbinom"]] <- effect
    effectSEl[["nbinom"]] <- effectSE
    effectCIlowerl[["nbinom"]] <- effectCIlower
    effectCIupperl[["nbinom"]] <- effectCIupper
  }

  if ("scoreSPA" %in% test) {
    P <- OR <- effect <- effectSE <- effectCIupper <- effectCIlower <- rep(
      NA_real_,
      nrow(GT)
    )

    if (length(testable) > 0) {
      score.null <- SPAtest::ScoreTest_wSaddleApprox_NULL_Model(
        null,
        data = as.data.frame(colData(GT))
      )
      fit <- SPAtest::ScoreTest_SPA(
        genos = assays(GT[, names(score.null$y)])$GT[testable, ],
        obj.null = score.null,
        minmac = 0,
        beta.out = FALSE
      )
      P[testable] <- fit$p.value
    }

    Pl[["scoreSPA"]] <- P
    ORl[["scoreSPA"]] <- OR
    effectl[["scoreSPA"]] <- effect
    effectSEl[["scoreSPA"]] <- effectSE
    effectCIlowerl[["scoreSPA"]] <- effectCIlower
    effectCIupperl[["scoreSPA"]] <- effectCIupper
  }

  if (continuous) {
    res <- cbind(
      results,
      effect = effectl[["lm"]],
      effectSE = effectSEl[["lm"]],
      effectCIlower = effectCIlowerl[["lm"]],
      effectCIupper = effectCIupperl[["lm"]],
      OR = ORl[["lm"]],
      P = Pl[["lm"]]
    )
  } else {
    res <- cbind(
      results,
      effect = c(do.call(rbind, effectl[test])),
      effectSE = c(do.call(rbind, effectSEl[test])),
      effectCIlower = c(do.call(rbind, effectCIlowerl[test])),
      effectCIupper = c(do.call(rbind, effectCIupperl[test])),
      OR = c(do.call(rbind, ORl[test])),
      P = c(do.call(rbind, Pl[test]))
    )
  }
  rownames(res) <- NULL
  res
}

# ACAT-v ----------------------------------------------------------------------
# adapted from: https://github.com/yaowuliu/ACAT
.acatv_rvat <- function(
  G,
  obj = NULL,
  data = NULL,
  covar = NULL,
  model = NULL,
  method = c("original", "scoreSPA", "firth"),
  weights.beta = c(1, 25),
  weights = NULL,
  mac.thresh = 10,
  maf = NULL,
  mac = NULL,
  maxitFirth = 1000L
) {
  method <- match.arg(method)

  ### check if null model is provided (if method = "original")
  if (method == "original" && is.null(obj)) {
    stop(
      "if method == 'original', a NULL model should be provided",
      call. = FALSE
    )
  }

  ### check if weights match length of genotypes
  if (!is.null(weights) && length(weights) != ncol(G)) {
    stop(
      "The length of weights must equal to the number of variants!",
      call. = FALSE
    )
  }

  n <- nrow(G)
  if (is.null(mac)) {
    mac <- Matrix::colSums(G, na.rm = TRUE)
  }
  if (is.null(maf)) {
    maf <- mac / (2 * n)
  }
  p <- length(mac)

  ### remove SNPs with mac=0
  if (sum(mac == 0) > 0) {
    G <- G[, mac > 0, drop = FALSE]
    weights <- weights[mac > 0]
    maf <- maf[mac > 0]
    mac <- mac[mac > 0]
    if (length(mac) == 0) {
      stop(
        "The genotype matrix does not have non-zero elements!",
        call. = FALSE
      )
    }
  }

  is.very.rare <- mac <= mac.thresh

  ## Very rare
  if (sum(is.very.rare) > 0) {
    if (method == "original") {
      pval.very.rare <- .acat_burden(
        G[, is.very.rare, drop = FALSE],
        obj,
        weights.beta = weights.beta,
        weights = weights[is.very.rare]
      )
    } else if (method == "scoreSPA") {
      w <- if (is.null(weights)) {
        dbeta(maf[is.very.rare], weights.beta[1], weights.beta[2])
      } else {
        weights[is.very.rare]
      }
      agg <- Matrix::rowSums(
        as(G[, is.very.rare, drop = FALSE], "sparseMatrix") %*%
          diag(w, ncol = length(w), nrow = length(w))
      )
      pval.very.rare <- SPAtest::ScoreTest_SPA(
        genos = agg,
        obj.nul = obj,
        minmac = 0
      )$p.value
    } else if (method == "firth") {
      w <- if (is.null(weights)) {
        dbeta(maf[is.very.rare], weights.beta[1], weights.beta[2])
      } else {
        weights[is.very.rare]
      }
      agg <- Matrix::rowSums(
        as(G[, is.very.rare, drop = FALSE], "sparseMatrix") %*%
          diag(w, ncol = length(w), nrow = length(w))
      )
      fit <- logistf::logistf(
        model,
        data = cbind(
          aggregate = agg,
          as.data.frame(data[, colnames(data) != "aggregate"])
        ),
        plconf = (which(c(covar, "aggregate") == "aggregate") + 1),
        control = logistf::logistf.control(maxit = maxitFirth),
        plcontrol = logistf::logistpl.control(maxit = maxitFirth)
      )
      if (.check_conv_firth(fit, maxit = maxitFirth)) {
        pval.very.rare <- fit$prob["aggregate"]
      } else {
        pval.very.rare <- NULL
      }
    }
  }

  ## !Very rare
  if (sum(!is.very.rare) > 0) {
    if (is.null(weights)) {
      weights.transformed <- (dbeta(
        maf[!is.very.rare],
        weights.beta[1],
        weights.beta[2]
      ) /
        dbeta(maf[!is.very.rare], 0.5, 0.5))^2
    } else {
      weights.transformed <- (weights[!is.very.rare] /
        dbeta(maf[!is.very.rare], 0.5, 0.5))^2
    }

    if (method == "original") {
      Mpvals <- .acat_Get.marginal.pval(G[, !is.very.rare, drop = FALSE], obj)
    } else if (method == "scoreSPA") {
      Mpvals <- SPAtest::ScoreTest_SPA(
        genos = t(G[, !is.very.rare, drop = FALSE]),
        obj.null = obj,
        minmac = 0
      )$p.value
    } else if (method == "firth") {
      Mpvals <- rep(NA_real_, sum(!is.very.rare))
      names(Mpvals) <- names(is.very.rare)[!is.very.rare]
      for (i in which(!is.very.rare)) {
        fit <- logistf::logistf(
          model,
          data = cbind(
            aggregate = G[, i],
            as.data.frame(data[, colnames(data) != "aggregate"])
          ),
          plconf = (which(c(covar, "aggregate") == "aggregate") + 1),
          control = logistf::logistf.control(maxit = maxitFirth),
          plcontrol = logistf::logistpl.control(maxit = maxitFirth)
        )
        if (.check_conv_firth(fit, maxit = maxitFirth)) {
          Mpvals[colnames(G)[i]] <- fit$prob["aggregate"]
        }
      }
    }
  }

  if (sum(is.very.rare) == 0) {
    pval <- .rvat_ACAT(Mpvals, weights.transformed)
  } else if (sum(!is.very.rare) == 0) {
    pval <- pval.very.rare
  } else {
    pvals <- c(Mpvals, pval.very.rare)
    mafs <- c(maf[!is.very.rare], mean(maf[is.very.rare]))

    if (is.null(weights)) {
      weights.transformed <- (dbeta(mafs, weights.beta[1], weights.beta[2]) /
        dbeta(mafs, 0.5, 0.5))^2
    } else {
      weights.transformed <- (c(
        weights[!is.very.rare],
        mean(weights[is.very.rare])
      ) /
        dbeta(mafs, 0.5, 0.5))^2
    }

    is.keep <- rep(TRUE, length(pvals))
    is.keep[which(pvals == 1)] <- FALSE ## remove p-values of 1.
    pval <- .rvat_ACAT(pvals[is.keep], weights.transformed[is.keep])
  }
  pval
}

# adapted from: https://github.com/yaowuliu/ACAT
.acat_NULL_Model <- function(Y, Z = NULL) {
  n <- length(Y)
  #### check the type of Y
  if ((sum(Y == 0) + sum(Y == 1)) == n) {
    out_type <- "D"
  } else {
    out_type <- "C"
  }
  #### Add intercept
  Z.tilde <- cbind(rep(1, length(Y)), Z)
  if (out_type == "C") {
    #### estimate of sigma square
    Z.med <- Z.tilde %*% solve(chol(t(Z.tilde) %*% Z.tilde)) ## Z.med%*%t(Z.med) is the projection matrix of Z.tilde
    Y.res <- as.vector(Y - (Y %*% Z.med) %*% t(Z.med))
    sigma2 <- sum(Y.res^2) / (n - ncol(Z.med))
    #### output
    res <- list()
    res[["out_type"]] <- out_type
    res[["Z.med"]] <- Z.med
    res[["Y.res"]] <- Y.res
    res[["sigma2"]] <- sigma2
  } else if (out_type == "D") {
    #### fit null model
    g <- glm(Y ~ 0 + Z.tilde, family = "binomial")
    prob.est <- g[["fitted.values"]]
    #### unstandarized residuals
    Y.res <- (Y - prob.est)
    ### Sigma when rho=0
    sigma2.Y <- prob.est * (1 - prob.est) ### variance of each Y_i
    ### output
    res <- list()
    res[["out_type"]] <- out_type
    res[["Z.tilde"]] <- Z.tilde
    res[["Y.res"]] <- Y.res
    res[["sigma2.Y"]] <- sigma2.Y
  }
  return(res)
}

# adapted from: https://github.com/yaowuliu/ACAT
.acat_Get.marginal.pval <- function(G, obj) {
  ### check obj
  if (names(obj)[1] != "out_type") {
    stop("obj is not calculated from MOAT_NULL_MODEL!", call. = FALSE)
  } else {
    out_type <- obj[["out_type"]]
    if (out_type == "C") {
      if (
        !all.equal(names(obj)[2:length(obj)], c("Z.med", "Y.res", "sigma2"))
      ) {
        stop("obj is not calculated from MOAT_NULL_MODEL!", call. = FALSE)
      } else {
        Z.med <- obj[["Z.med"]]
        Y.res <- obj[["Y.res"]]
        n <- length(Y.res)
        SST <- obj[["sigma2"]] * (n - ncol(Z.med))
      }
    } else if (out_type == "D") {
      if (
        !all.equal(names(obj)[2:length(obj)], c("Z.tilde", "Y.res", "sigma2.Y"))
      ) {
        stop("obj is not calculated from MOAT_NULL_MODEL!", call. = FALSE)
      } else {
        Z.tilde <- obj[["Z.tilde"]]
        Y.res <- obj[["Y.res"]]
        sigma2.Y <- obj[["sigma2.Y"]]
        n <- length(Y.res)
      }
    }
  }

  if (!"matrix" %in% class(G) && !"dgCMatrix" %in% class(G)) {
    stop("The class of G must be matrix or dgCMatrix!", call. = FALSE)
  }

  if (out_type == "C") {
    G_tX.med <- as.matrix(Matrix::crossprod(Z.med, G))
    ### Sigma^2 of G
    Sigma2.G <- Matrix::colSums(G^2) - Matrix::colSums(G_tX.med^2)
    SSR <- as.vector((Y.res %*% G)^2 / Sigma2.G)
    SSR[Sigma2.G <= 0] <- 0
    df.2 <- n - 1 - ncol(Z.med)
    t.stat <- suppressWarnings(sqrt(SSR / ((SST - SSR) / df.2)))
    marginal.pvals <- 2 * pt(t.stat, (n - 1 - ncol(Z.med)), lower.tail = FALSE)
  } else if (out_type == "D") {
    Z.stat0 <- as.vector(Y.res %*% G)
    ### Sigma when rho=0
    tG_X.tilde_sigma2 <- as.matrix(Matrix::crossprod(G, Z.tilde * sigma2.Y))
    Sigma2.G <- Matrix::colSums(G^2 * sigma2.Y) -
      diag(
        tG_X.tilde_sigma2 %*%
          solve(t(Z.tilde) %*% (Z.tilde * sigma2.Y)) %*%
          t(tG_X.tilde_sigma2)
      )
    marginal.pvals <- 2 *
      pnorm(abs(Z.stat0) / sqrt(Sigma2.G), lower.tail = FALSE)
  }

  return(marginal.pvals)
}

# adapted from: https://github.com/yaowuliu/ACAT
.acat_burden <- function(
  G,
  obj,
  kernel = "linear.weighted",
  weights.beta = c(1, 25),
  weights = NULL
) {
  ### check obj
  if (names(obj)[1] != "out_type") {
    stop("obj is not calculated from NULL_MODEL!", call. = FALSE)
  } else {
    out_type <- obj[["out_type"]]
    if (out_type == "C") {
      if (
        !all.equal(names(obj)[2:length(obj)], c("Z.med", "Y.res", "sigma2"))
      ) {
        stop("obj is not calculated from NULL_MODEL!", call. = FALSE)
      } else {
        Z.med <- obj[["Z.med"]]
        Y.res <- obj[["Y.res"]] / sqrt(obj[["sigma2"]]) ## rescaled residules
        n <- length(Y.res)
      }
    } else if (out_type == "D") {
      if (
        !all.equal(names(obj)[2:length(obj)], c("Z.tilde", "Y.res", "sigma2.Y"))
      ) {
        stop("obj is not calculated from NULL_MODEL!", call. = FALSE)
      } else {
        Z.tilde <- obj[["Z.tilde"]]
        Y.res <- obj[["Y.res"]]
        sigma2.Y <- obj[["sigma2.Y"]]
        n <- length(Y.res)
      }
    }
  }
  p <- ncol(G)
  #### weights
  if (kernel == "linear.weighted") {
    if (is.null(weights)) {
      MAF <- Matrix::colSums(G) / (2 * dim(G)[1])
      W <- dbeta(MAF, weights.beta[1], weights.beta[2])
    } else {
      if (length(weights) == p) {
        W <- weights
      } else {
        stop(
          "The length of weights must equal to the number of variants!",
          call. = FALSE
        )
      }
    }
  } else if (kernel == "linear") {
    W <- rep(1, p)
  } else {
    stop("The kernel name is not valid!", call. = FALSE)
  }

  ###### if G is sparse or not
  if ("matrix" %in% class(G) || "dgCMatrix" %in% class(G)) {
    if (out_type == "C") {
      Z.stat.sum <- as.vector((Y.res %*% G) %*% W)
      Gw <- G %*% W
      sigma.z <- sqrt(sum(Gw^2) - sum((t(Z.med) %*% (Gw))^2))
    } else if (out_type == "D") {
      Z.stat.sum <- as.vector((Y.res %*% G) %*% W)
      Gw <- as.vector(G %*% W)
      sigma.z <- sum(Gw^2 * sigma2.Y) -
        ((Gw * sigma2.Y) %*% Z.tilde) %*%
          solve(t(Z.tilde) %*% (Z.tilde * sigma2.Y)) %*%
          t((Gw * sigma2.Y) %*% Z.tilde)
      sigma.z <- as.vector(sqrt(sigma.z))
    }
  } else {
    stop("The class of G must be matrix or dgCMatrix!", call. = FALSE)
  }

  V <- Z.stat.sum / sigma.z
  Q <- V^2 ## Q test statistic
  pval <- 1 - pchisq(Q, df = 1)
  return(pval)
}

# adapted from: https://github.com/yaowuliu/ACAT
.rvat_ACAT <- function(Pvals, weights = NULL, is.check = TRUE) {
  Pvals <- as.matrix(Pvals)
  if (is.check) {
    #### check if there is NA
    if (sum(is.na(Pvals)) > 0) {
      stop("Cannot have NAs in the p-values!", call. = FALSE)
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals < 0) + sum(Pvals > 1)) > 0) {
      stop("P-values must be between 0 and 1!", call. = FALSE)
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero <- (colSums(Pvals == 0) >= 1)
    is.one <- (colSums(Pvals == 1) >= 1)
    if (sum((is.zero + is.one) == 2) > 0) {
      stop(
        "Cannot have both 0 and 1 p-values in the same column!",
        call. = FALSE
      )
    }

    if (sum(is.zero) > 0) {
      warning("There are p-values that are exactly 0!", call. = FALSE)
    }
    if (sum(is.one) > 0) {
      warning("There are p-values that are exactly 1!", call. = FALSE)
    }
  }
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(weights)) {
    is.weights.null <- TRUE
  } else {
    is.weights.null <- FALSE
    weights <- as.matrix(weights)
    if (sum(dim(weights) != dim(Pvals)) > 0) {
      stop(
        "The dimensions of weights and Pvals must be the same!",
        call. = FALSE
      )
    } else if (is.check && (sum(weights < 0) > 0)) {
      stop("All the weights must be nonnegative!", call. = FALSE)
    } else {
      w.sum <- colSums(weights)
      if (sum(w.sum <= 0) > 0) {
        stop(
          "At least one weight should be positive in each column!",
          call. = FALSE
        )
      } else {
        for (j in 1:ncol(weights)) {
          weights[, j] <- weights[, j] / w.sum[j]
        }
      }
    }
  }

  # check if there are very small non-zero p values and calcuate the cauchy statistics
  is.small <- (Pvals < 1e-15)
  if (is.weights.null) {
    Pvals[!is.small] <- tan((0.5 - Pvals[!is.small]) * pi)
    Pvals[is.small] <- 1 / Pvals[is.small] / pi
    cct.stat <- colMeans(Pvals)
  } else {
    Pvals[!is.small] <- weights[!is.small] * tan((0.5 - Pvals[!is.small]) * pi)
    Pvals[is.small] <- (weights[is.small] / Pvals[is.small]) / pi
    cct.stat <- colSums(Pvals)
  }
  #### return the ACAT p value(s).
  pval <- pcauchy(cct.stat, lower.tail = FALSE)
  return(pval)
}

# misc ------------------------------------------------------------------------
.check_conv_firth <- function(fit, maxit = NULL) {
  if (
    max(fit$pl.iter[, "Lower"]) >= maxit ||
      max(fit$pl.iter[, "Upper"]) >= maxit ||
      max(fit$pl.iter[, "Null model"]) >= maxit
  ) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}


# Permutations ----------------------------------------------------------------
.permute_rvb <- function(
  GT,
  test,
  pheno,
  model,
  null,
  covar,
  continuous,
  perms,
  methodResampling = "permutation"
) {
  # Currently, tests included for permutation:
  P <- list()

  if (sum(colData(GT)$aggregate > 0) < 2) {
    warning(
      "Less than two samples have a non-zero burden score, skipping tests.",
      call. = FALSE
    )
    test <- c()
  }

  if (continuous) {
    out_type <- "C"
  } else {
    out_type <- "D"
    if (
      sum(colData(GT)[, pheno] == 1) < 2 || sum(colData(GT)[, pheno] == 0) < 2
    ) {
      warning(
        "Fewer than two cases or controls, skipping tests.",
        call. = FALSE
      )
      test <- c()
    }
  }

  # SKAT analyses
  if (
    sum(
      c(
        "skat_burden",
        "skat",
        "skato",
        "skat_burden_robust",
        "skat_robust",
        "skato_robust"
      ) %in%
        test
    ) >
      0
  ) {
    skat.null <- SKAT::SKAT_Null_Model(
      null,
      data = colData(GT),
      out_type = out_type
    )
    skat.null$n.Resampling <- ncol(perms)
    skat.null$type.Resampling <- methodResampling
    skat.null$res.out <- apply(perms, 2, function(x) skat.null$res[x])

    if ("skat_burden" %in% test) {
      tryCatch(
        {
          if (continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null,
              r.corr = 1,
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          } else {
            fit <- SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "Burden",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }

          P[["skat_burden"]] <- fit$p.value.resampling
        },
        error = function(e) {
          message(sprintf("Failed test '%s'\n%s", "skat_burden", e))
        }
      )
      if (length(P[["skat_burden"]]) == 0) {
        P[["skat_burden"]] <- rep(NA_real_, ncol(perms))
      }
    }

    if ("skat" %in% test) {
      tryCatch(
        {
          if (continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null,
              r.corr = 0,
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          } else {
            fit <- SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "SKAT",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }

          P[["skat"]] <- fit$p.value.resampling
        },
        error = function(e) {
          message(sprintf("Failed test '%s'\n%s", "skat", e))
        }
      )
      if (length(P[["skat"]]) == 0) P[["skat"]] <- rep(NA_real_, ncol(perms))
    }

    if ("skato" %in% test) {
      tryCatch(
        {
          if (continuous) {
            fit <- SKAT::SKAT(
              t(assays(GT)$GT),
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          } else {
            fit <- SKAT::SKATBinary(
              t(assays(GT)$GT),
              skat.null,
              method = "SKATO",
              impute.method = "fixed",
              weights = rowData(GT)$w,
              missing_cutoff = 1
            )
          }

          P[["skato"]] <- fit$p.value.resampling
        },
        error = function(e) {
          message(sprintf("Failed test '%s'\n%s", "skato", e))
        }
      )
      if (length(P[["skato"]]) == 0) P[["skato"]] <- rep(NA_real_, ncol(perms))
    }

    if (
      sum(c("skat_burden_robust", "skat_robust", "skato_robust") %in% test) > 0
    ) {
      varkeep <- Matrix::rowSums(assays(GT)$GT) > 0

      if (sum(varkeep) > 0) {
        skat_robust_perms <- lapply(
          1:ncol(perms),
          .permute_skatrobust_skato,
          gt = t(assays(GT)$GT)[, varkeep, drop = FALSE],
          skat_null = skat.null,
          weights = rowData(GT)$w[varkeep],
          perm = perms
        )
        skat_robust_perms <- matrix(
          unlist(skat_robust_perms),
          ncol = 3,
          byrow = TRUE
        )
        colnames(skat_robust_perms) <- c(
          "skat_robust",
          "skat_burden_robust",
          "skato_robust"
        )
        P[["skat_robust"]] <- if ("skat_robust" %in% test) {
          skat_robust_perms[, "skat_robust"]
        } else {
          NULL
        }
        P[["skato_robust"]] <- if ("skato_robust" %in% test) {
          skat_robust_perms[, "skato_robust"]
        } else {
          NULL
        }
        P[["skat_burden_robust"]] <- if ("skat_burden_robust" %in% test) {
          skat_robust_perms[, "skat_burden_robust"]
        }
      }
      if (length(P[["skat_robust"]]) == 0 & "skat_robust" %in% test) {
        P[["skat_robust"]] <- rep(NA_real_, ncol(perms))
      }
      if (length(P[["skato_robust"]]) == 0 & "skato_robust" %in% test) {
        P[["skato_robust"]] <- rep(NA_real_, ncol(perms))
      }
      if (
        length(P[["skat_burden_robust"]]) == 0 & "skat_burden_robust" %in% test
      ) {
        P[["skat_burden_robust"]] <- rep(NA_real_, ncol(perms))
      }
    }
  }

  if (sum(c("acatv") %in% test) > 0) {
    acat_null = .acat_NULL_Model(
      colData(GT)[[pheno]],
      Z = if (!is.null(covar)) {
        as.matrix(colData(GT)[, covar, drop = FALSE])
      } else {
        NULL
      }
    )

    if ("acatv" %in% test) {
      tryCatch(
        {
          # run permutations
          acatv_perms <- unlist(lapply(
            1:ncol(perms),
            .permute_acatv,
            gt = t(assays(GT)$GT),
            acat_null = acat_null,
            perm = perms,
            weights = rowData(GT)$w,
            maf = getAF(GT),
            mac = Matrix::rowSums(assays(GT)$GT)
          ))

          P[["acatv"]] <- acatv_perms
        },
        error = function(e) {
          message(sprintf("Failed test '%s'\n%s", "acatv", e))
        }
      )
      if (length(P[["acatv"]]) == 0) P[["acatv"]] <- rep(NA_real_, ncol(perms))
    }
  }

  data.frame(
    pheno = rep(colnames(perms), times = length(test)),
    test = rep(
      assocTest_resampling_tests[assocTest_resampling_tests %in% test],
      each = ncol(perms)
    ),
    P = unname(unlist(P[assocTest_resampling_tests[
      assocTest_resampling_tests %in% test
    ]])),
    stringsAsFactors = FALSE
  )
}

.permute_skatrobust_skato <- function(i, gt, skat_null, weights, perms) {
  perm <- perms[, i]
  skat_null_permuted <- skat_null
  skat_null_permuted$res <- skat_null_permuted$res[perm]
  skat_null_permuted$mu <- skat_null_permuted$mu[perm]
  skat_null_permuted$pi_1 <- skat_null_permuted$pi_1[perm]
  skat_null_permuted$X1 <- skat_null_permuted$X1[perm, , drop = FALSE]
  fit = SKAT::SKATBinary_Robust(
    gt,
    skat_null_permuted,
    method = "SKATO",
    impute.method = "fixed",
    weights = weights,
    missing_cutoff = 1
  )
  if (ncol(gt) > 1) {
    c(
      fit$p.value_each[fit$param$rho == 0],
      fit$p.value_each[fit$param$rho == 1],
      fit$p.value
    )
  } else {
    c(fit$p.value, fit$p.value, fit$p.value)
  }
}

.permute_acatv <- function(i, gt, acat_null, perms, weights, maf, mac) {
  perm <- perms[, i]
  acat_null_permuted <- acat_null
  acat_null_permuted$Z.tilde <- acat_null_permuted$Z.tilde[perm, , drop = FALSE]
  acat_null_permuted$Y.res <- acat_null_permuted$Y.res[perm]
  acat_null_permuted$sigma2.Y <- acat_null_permuted$sigma2.Y[perm]
  .acatv_rvat(
    gt,
    obj = acat_null_permuted,
    method = "original",
    weights = weights,
    mac.thresh = 10,
    maf = maf,
    mac = mac
  )
}


get_perm_pval <- function(pobs, pperm) {
  P <- (sum(pperm < pobs) + 1) / (length(pperm) + 1)
  P
}

get_perm_pvals <- function(x, y) {
  results <- lapply(
    unique(as.character(y$test)),
    FUN = function(tst, x, y) {
      x_ <- x[x$test == tst, ]
      x_$P <- get_perm_pval(x_$P, y[y$test == tst, ]$P)
      x_$test <- paste0(tst, "_resampled")
      x_
    },
    x = x,
    y = y
  )
  do.call(rbind, results)
}
