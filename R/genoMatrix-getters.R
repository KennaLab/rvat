#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod("getAF", signature = "genoMatrix", definition = function(object) {
  meta <- metadata(object)
  gt <- assays(object)$GT

  # return vector of NAs if geneticModel is not 'allelic'
  if (meta$geneticModel != "allelic") {
    warning(
      sprintf(
        "Allele frequencies are set to NA for geneticModel '%s'.",
        meta$geneticModel
      ),
      call. = FALSE
    )
    af <- rep(NA_real_, nrow(object))
    names(af) <- rownames(object)
    return(af)
  }

  # return numeric(0) if no variants are present in genoMatrix
  if (nrow(object) == 0L) {
    return(numeric(0))
  }

  # diploid
  if (all(meta$ploidyLevels == "diploid")) {
    AF <- Matrix::rowMeans(gt, na.rm = TRUE) / 2.0
  } else {
    # non-diploid variants
    rowdata <- rowData(object)
    coldata <- colData(object)

    # initialize AF vector with NAs
    AF <- rep(NA_real_, nrow(object))
    names(AF) <- rownames(object)

    # diploid
    if ("diploid" %in% meta$ploidyLevels) {
      AF[rowdata$ploidy == "diploid"] <- Matrix::rowMeans(
        gt[rowdata$ploidy == "diploid", , drop = FALSE],
        na.rm = TRUE
      ) /
        2.0
    }

    # XnonPAR and YnonPAR

    ## issue warning if sex is unknown/missing for (a subset of) samples
    if (
      any(c("XnonPAR", "YnonPAR") %in% meta$ploidyLevels) &&
        any(object$sex == 0)
    ) {
      warning(
        sprintf(
          paste0(
            "GT contains variants with ploidy=%s, ",
            "MAF is calculated in samples with non-missing sex info. ",
            "Note that sex is missing for %s samples"
          ),
          paste(
            paste0("'", meta$ploidyLevels[meta$ploidyLevels != "diploid"], "'"),
            collapse = ","
          ),
          sum(object$sex == 0)
        ),
        call. = FALSE
      )
    }

    # XnonPAR
    if ("XnonPAR" %in% meta$ploidyLevels) {
      AF[rowdata$ploidy == "XnonPAR"] <- .getAF_XnonPAR(
        gt[rowdata$ploidy == "XnonPAR", , drop = FALSE],
        coldata = coldata,
        meta = meta
      )
    }

    # YnonPAR
    if ("YnonPAR" %in% meta$ploidyLevels) {
      AF[rowdata$ploidy == "YnonPAR"] <- .getAF_YnonPAR(
        gt[rowdata$ploidy == "YnonPAR", , drop = FALSE],
        coldata = coldata
      )
    }
  }

  # return AF
  AF
})


.getAF_XnonPAR <- function(
  GT,
  coldata,
  meta
) {
  is_male <- coldata$sex == 1
  is_female <- coldata$sex == 2
  if ((sum(is_male) + sum(is_female)) == 0L) {
    return(rep(NA_real_, nrow(GT)))
  }

  if (meta$imputeMethod != "none") {
    ac_males <- rowSums(
      floor(GT[, is_male, drop = FALSE]),
      na.rm = TRUE
    )
    ac_females <- rowSums(
      floor(GT[, is_female, drop = FALSE]),
      na.rm = TRUE
    )
  } else {
    ac_males <- rowSums(
      GT[, is_male, drop = FALSE],
      na.rm = TRUE
    )
    ac_females <- rowSums(
      GT[, is_female, drop = FALSE],
      na.rm = TRUE
    )
  }

  non_missing_males <- Matrix::rowSums(
    !is.na(GT[, is_male, drop = FALSE])
  )
  non_missing_females <- Matrix::rowSums(
    !is.na(GT[, is_female, drop = FALSE])
  )
  AF_XnonPAR <- (ac_males + ac_females) /
    (non_missing_males + (2 * non_missing_females))

  AF_XnonPAR
}

.getAF_YnonPAR <- function(GT, coldata) {
  is_male <- coldata$sex == 1
  if (sum(is_male) == 0L) {
    return(rep(NA_real_, nrow(GT)))
  }
  AF_YnonPAR <- Matrix::rowMeans(
    GT[, is_male, drop = FALSE],
    na.rm = TRUE
  )

  AF_YnonPAR
}


#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod("getMAF", signature = "genoMatrix", definition = function(object) {
  object <- flipToMinor(object)
  getAF(object)
})

#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod("getAC", signature = "genoMatrix", definition = function(object) {
  # calculate allele counts if geneticModel == 'allelic'
  if (metadata(object)$geneticModel == "allelic") {
    # floor if imputed
    if (metadata(object)$imputeMethod != "none") {
      ac <- rowSums(
        floor(assays(object)$GT),
        na.rm = TRUE
      )
    } else {
      ac <- rowSums(assays(object)$GT, na.rm = TRUE)
    }
  } else {
    # AC can only be calculated when geneticModel is 'allelic'
    warning(
      sprintf(
        "Allele counts can't be calculated for geneticModel '%s'. Returning NA vector.",
        metadata(object)$geneticModel
      ),
      call. = FALSE
    )
    ac <- rep(NA_real_, nrow(object))
    if (nrow(object) > 0L) {
      names(ac) <- rownames(object)
    }
  }

  ac
})

#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod("getMAC", signature = "genoMatrix", definition = function(object) {
  object <- flipToMinor(object)
  getAC(object)
})

#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod(
  "getNCarriers",
  signature = "genoMatrix",
  definition = function(object) {
    carriers <- rowSums(
      assays(object)$GT >= 1,
      na.rm = TRUE
    )
    carriers
  }
)

#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod(
  "getCR",
  signature = "genoMatrix",
  definition = function(object, var = TRUE) {
    if (var) {
      Matrix::rowMeans(!is.na(assays(object)$GT))
    } else {
      Matrix::colMeans(!is.na(assays(object)$GT))
    }
  }
)
