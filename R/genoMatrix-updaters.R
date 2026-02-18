#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod(
  "flipToMinor",
  signature = "genoMatrix",
  definition = function(object) {
    # flipToMinor works for allelic model only
    if (metadata(object)$geneticModel != "allelic") {
      warning(
        "flipToMinor only applies when geneticModel == 'allelic', ",
        "genoMatrix is returned unchanged.",
        call. = FALSE
      )
      return(object)
    }

    # variants to flip
    flip <- getAF(object) > 0.5
    flip[is.na(flip)] <- FALSE
    
    # flip effect alleles in variant info
    if (all(c("effectAllele", "otherAllele") %in% colnames(rowData(object)))) {
      effectAllele <- ifelse(
        flip,
        rowData(object)[["otherAllele"]],
        rowData(object)[["effectAllele"]]
      )
      rowData(object)[["otherAllele"]] <- ifelse(
        flip,
        rowData(object)[["effectAllele"]],
        rowData(object)[["otherAllele"]]
      )
      rowData(object)[["effectAllele"]] <- effectAllele
    } else if (all(c("REF", "ALT") %in% colnames(rowData(object)))) {
      rowData(object)[["effectAllele"]] <- ifelse(
        flip,
        rowData(object)[["REF"]],
        rowData(object)[["ALT"]]
      )
      rowData(object)[["otherAllele"]] <- ifelse(
        flip,
        rowData(object)[["ALT"]],
        rowData(object)[["REF"]]
      )
    }

    # swap effect allele for flipped variants
    if (sum(flip) > 0L) {
      if (all(metadata(object)$ploidyLevels == "diploid")) {
        assays(object)$GT <- abs(
          assays(object)$GT -
            2 *
              matrix(
                rep(flip, each = metadata(object)$m),
                nrow = metadata(object)$nvar,
                byrow = TRUE
              )
        )
      } else {
        # if genoMatrix contains non-diploid variants,
        # apply appropriate flipping based on ploidy
        GT <- assays(object)$GT
        for (i in which(flip)) {
          if (rowData(object)$ploidy[i] == "XnonPAR") {
            # flip 0s and 1s for males
            for (i2 in which(colData(object)$sex == 1)) {
              GT[i, i2] <- abs(GT[i, i2] - 1)
            }
            # flip 0s and 2s for females
            for (i2 in which(colData(object)$sex == 2)) {
              GT[i, i2] <- abs(GT[i, i2] - 2)
            }
          }
          if (rowData(object)$ploidy[i] == "YnonPAR") {
            # flip 0s and 1s
            GT[i, ] <- abs(GT[i, ] - 1)
          }

          if (rowData(object)$ploidy[i] == "diploid") {
            # flip 0s and 2s
            GT[i, ] <- abs(GT[i, ] - 2)
          }
        }
        assays(object)$GT <- GT
      }
    }

    object
  }
)

.calc_maf_weights <- function(w, af = NULL, method = "none") {
  if (!method %in% c("mb", "none")) {
    stop("method should be either 'none' or 'mb'.", call. = FALSE)
  }

  if (!is.null(af) && length(w) != length(af)) {
    stop("Unequal lengths", call. = FALSE)
  }

  if (method == "none") {
    w
  } else if (method == "mb") {
    if (is.null(af)) {
      stop("AF should be specified for 'mb' weighting", call. = FALSE)
    }
    w / (sqrt(af * (1 - af)))
  }
}

#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod(
  "recode",
  signature = "genoMatrix",
  definition = function(
    object,
    geneticModel = NULL,
    imputeMethod = NULL,
    weights = NULL,
    MAFweights = NULL
  ) {
    # validate input
    .recode_validate_input(as.list(environment()))

    # return early if genoMatrix is empty
    if (nrow(object) == 0L) {
      if (!is.null(imputeMethod)) {
        metadata(object)$imputeMethod <- imputeMethod
      }
      return(object)
    }

    # set custom weights, if weights are of length 1 they are recycled
    if (!is.null(weights)) {
      if (length(weights) == 1L) {
        rowData(object)$w <- as.numeric(rep(weights, nrow(object)))
      } else {
        rowData(object)$w <- as.numeric(weights)
      }
    }

    # set MAF weights (only 'mb' currently)
    if (!is.null(MAFweights) && MAFweights == "mb") {
      rowData(object)$w <- .calc_maf_weights(
        w = rowData(object)$w,
        af = getAF(object),
        method = "mb"
      )
    }

    # set geneticModel
    if (
      !is.null(geneticModel) &&
        metadata(object)$geneticModel != geneticModel
    ) {
      object <- flipToMinor(object)
      if (geneticModel == "allelic") {
        assays(object)$GT <- assays(object)$GT
      } else if (geneticModel == "dominant") {
        assays(object)$GT <- (assays(object)$GT > 0) * 1
      } else if (geneticModel == "recessive") {
        assays(object)$GT <- (assays(object)$GT == 2) * 1
      }
      metadata(object)$geneticModel <- geneticModel
    }

    # impute
    if (!is.null(imputeMethod)) {
      # skip if already imputed
      if (metadata(object)$imputeMethod != "none") {
        message(sprintf(
          "GT is already imputed (method='%s'), skipping imputation",
          metadata(object)$imputeMethod
        ))
      } else {
        # mean impute missing genotype dosages
        if (imputeMethod == "meanImpute") {
          missings <- which(is.na(assays(object)$GT), arr.ind = TRUE)
          if (nrow(missings) > 0L) {
            assays(object)$GT[missings] <- rowMeans(
              assays(object)$GT,
              na.rm = TRUE
            )[missings[, 1]]
          }
        }

        # reset missing genotypes to 0
        if (imputeMethod == "missingToRef") {
          assays(object)$GT[is.na(assays(object)$GT)] <- 0
        }

        # update metadata
        metadata(object)$imputeMethod <- imputeMethod
      }
    }

    # reset aggregates to missing if present
    if ("aggregate" %in% colnames(colData(object))) {
      colData(object)$aggregate <- NA_real_
    }

    # check validity and return
    validObject(object)
    object
  }
)

.recode_validate_parameters <- function(args) {
  # MAFweights should be of type character and length 1, and implemented
  check_wrapper(
    check_character,
    args,
    "MAFweights",
    length_equal = 1L,
    allow_null = TRUE
  )

  if (
    !is.null(args[["MAFweights"]]) &&
      !all(args[["MAFweights"]] %in% c("mb", "none"))
  ) {
    stop(
      "`MAFweights` parameter should be either 'none' or 'mb'.",
      call. = FALSE
    )
  }

  # geneticModel should be of type character, length 1, and implemented
  check_wrapper(
    check_character,
    args,
    "geneticModel",
    length_equal = 1L,
    allow_null = TRUE
  )
  if (
    !is.null(args[["geneticModel"]]) &&
      !any(args[["geneticModel"]] %in% c("allelic", "recessive", "dominant"))
  ) {
    stop(
      sprintf(
        "%s does not represent a valid genetic model",
        args[["geneticModel"]]
      ),
      call. = FALSE
    )
  }

  # imputeMethod should be of type character, length 1, and implemented
  check_wrapper(
    check_character,
    args,
    "imputeMethod",
    length_equal = 1L,
    allow_null = TRUE
  )
  if (
    !is.null(args[["imputeMethod"]]) &&
      (!args[["imputeMethod"]] %in% c("meanImpute", "missingToRef"))
  ) {
    stop(
      sprintf(
        "'%s' is not a recognized imputation method",
        args[["imputeMethod"]]
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

.recode_validate_input <- function(args) {
  # validate parameters
  .recode_validate_parameters(args)

  # check if parameters are consistent with genoMatrix state

  ## weights should be either of length 1 or equal to number of variants in genoMatrix
  if (
    !is.null(args[["weights"]]) &&
      length(args[["weights"]]) != 1L &&
      length(args[["weights"]]) != nrow(args[["object"]])
  ) {
    stop(
      "Length of `weights` should equal the number ",
      "of variants in the genoMatrix.",
      call. = FALSE
    )
  }

  # current model should be allelic in order to apply dominant/recessive models
  # can only be applied to a non-imputed genoMatrix
  if (!is.null(args[["geneticModel"]])) {
    if (
      metadata(args[["object"]])$geneticModel != args[["geneticModel"]] &&
        metadata(args[["object"]])$geneticModel != "allelic"
    ) {
      stop(
        "Current geneticModel should be 'allelic' in ",
        "order to apply dominant or recessive models.",
        call. = FALSE
      )
    }

    if (
      metadata(args[["object"]])$geneticModel != args[["geneticModel"]] &&
        metadata(args[["object"]])$imputeMethod != "none"
    ) {
      stop(
        "Provide a non-imputed genoMatrix to ",
        "perform dominant/recessive recoding",
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}


#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod(
  "updateGT",
  signature = "genoMatrix",
  definition = function(object, SM = NULL, anno = NULL) {
    # validate input
    .updateGT_validate_input(as.list(environment()))

    # process SM is provided
    if (!is.null(SM)) {
      SM <- DataFrame(SM)
      rownames(SM) <- SM[["IID"]]

      # match IIDs
      SM <- SM[match(colnames(object), SM$IID), , drop = FALSE]

      # reset colData
      colData(object) <- SM
    }

    # process anno is provided
    if (!is.null(anno)) {
      # early return if genoMatrix is empty
      if (nrow(object) == 0L) {
        warning(
          "Annotations can't be added to an empty genoMatrix. ",
          "The input genoMatrix is returned unchanged.",
          call. = FALSE
        )
        return(object)
      }

      # merge
      rowdata <- rowData(object)
      rowdata$VAR_id <- rownames(object)
      rowdata <- merge(
        rowdata[, c("VAR_id", "ploidy", "w")],
        anno,
        all.x = TRUE,
        by = "VAR_id"
      )
      rownames(rowdata) <- rowdata$VAR_id
      rowdata$VAR_id <- NULL

      # update rowData
      rowData(object) <- rowdata[rownames(object), ]
    }

    # validate and return
    validObject(object)
    return(object)
  }
)

.updateGT_validate_input <- function(args) {
  if (!is.null(args[["SM"]])) {
    .gdb_check_SM(
      SM = args[["SM"]],
      GT = args[["object"]]
    )
    if (!all(metadata(args[["object"]])$ploidyLevels == "diploid")) {
      stop(
        "Genotype dosage values at non-diploid sites were previously ",
        "reset according to the original sample sex values.",
        call. = FALSE
      )
    }
  }

  if (!is.null(args[["anno"]])) {
    if (!"VAR_id" %in% colnames(args[["anno"]])) {
      stop("anno should contain a `VAR_id` column", call. = FALSE)
    }
    if (any(c("ploidy", "w", "AF") %in% colnames(args[["anno"]]))) {
      stop(
        "`ploidy`, `w` and `AF` are protected rowData column names",
        call. = FALSE
      )
    }
    if (anyDuplicated(args[["anno"]]$VAR_id) != 0L) {
      stop("`anno` shouldn't contain duplicated VAR_ids.", call. = FALSE)
    }
    if (!all(rownames(args[["object"]]) %in% args[["anno"]]$VAR_id)) {
      warning(
        "Not all variants present in the genoMatrix are present in the ",
        "anno table. Fields for missing variants will be filled with NAs.",
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}
