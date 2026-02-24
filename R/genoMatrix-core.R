#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod("show", signature = "genoMatrix", function(object) {
  cat("rvat genoMatrix Object\n")
  callNextMethod()
})


#' @importFrom S4Vectors setValidity2
S4Vectors::setValidity2("genoMatrix", function(object) {
  msg <- character()

  # assays -----------------------------
  if (!"GT" %in% assayNames(object)) {
    msg <- c(msg, "Assay 'GT' must be present.")
  } else if (!is.matrix(assays(object)$GT)) {
    msg <- c(msg, "Assay 'GT' must be a matrix.")
  }

  # rowData ---------------------------
  if (!"ploidy" %in% colnames(rowData(object))) {
    msg <- c(msg, "'ploidy' column must be present in rowData.")
  } else if (!is.character(rowData(object)$ploidy)) {
    msg <- c(msg, "'ploidy' column in rowData must be of type character.")
  } else if (
    !is.null(metadata(object)$ploidyLevels) &&
      !all(rowData(object)$ploidy %in% metadata(object)$ploidyLevels)
  ) {
    msg <- c(
      msg,
      "All values in rowData$ploidy must be among those in metadata$ploidyLevels."
    )
  }

  if (!"w" %in% colnames(rowData(object))) {
    msg <- c(msg, "'w' column must be present in rowData.")
  } else if (!is.numeric(rowData(object)$w)) {
    msg <- c(msg, "'w' column in rowData must be numeric.")
  }

  # colData ---------------------------
  # sample IDs
  if (!"IID" %in% colnames(colData(object))) {
    msg <- c(msg, "'IID' column must be present in colData.")
  } else if (!identical(colData(object)$IID, colnames(object))) {
    msg <- c(msg, "'IID' column in colData should match colnames.")
  }

  # sex
  if (!"sex" %in% colnames(colData(object))) {
    msg <- c(
      msg,
      "'sex' column is mandatory in colData."
    )
  } else if (anyNA(colData(object)$sex)) {
    msg <- c(
      msg,
      "'sex' column in colData must not contain NA values."
    )
  } else if (!all(colData(object)$sex %in% c(0, 1, 2))) {
    msg <- c(
      msg,
      "'sex' column values in colData must be 0, 1, or 2."
    )
  }

  # metadata ---------------------------

  # ploidy levels
  if (
    nrow(object) > 0L &&
      !all(
        metadata(object)$ploidyLevels %in% c("diploid", "XnonPAR", "YnonPAR")
      )
  ) {
    msg <- c(msg, "Ploidy must be set to one of diploid, XnonPAR or YnonPAR")
  }

  # return potential messages, TRUE otherwise
  if (length(msg) == 0L) {
    TRUE
  } else {
    msg
  }
})


setMethod(
  ".resetSexChromDosage",
  signature = "genoMatrix",
  definition = function(object) {
    if (
      nrow(object) == 0L ||
        all(metadata(object)$ploidyLevels == "diploid")
    ) {
      return(object)
    }

    # XnonPAR
    # recode males to 0,1
    # set to NA if sex is unknown/missing
    assays(object)$GT[
      rowData(object)$ploidy == "XnonPAR",
      object$sex == 1
    ] <- as.numeric(
      assays(object)$GT[rowData(object)$ploidy == "XnonPAR", object$sex == 1] >=
        1
    )
    assays(object)$GT[
      rowData(object)$ploidy == "XnonPAR",
      object$sex == 0
    ] <- NA_real_

    # YnonPAR
    # recode males to 0,1
    # set to NA if female or sex is unknown/missing
    assays(object)$GT[
      rowData(object)$ploidy == "YnonPAR",
      object$sex == 1
    ] <- as.numeric(
      assays(object)$GT[rowData(object)$ploidy == "YnonPAR", object$sex == 1] >=
        1
    )
    assays(object)$GT[
      rowData(object)$ploidy == "YnonPAR",
      object$sex %in% c(0, 2)
    ] <- NA_real_

    # return
    object
  }
)

#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod("[", signature = "genoMatrix", function(x, i, j, drop = TRUE) {
  out <- callNextMethod()

  # update variant counts and ploidyLevels in event of variant filtering
  if (!missing(i)) {
    metadata(out)$nvar <- nrow(out)
    metadata(out)$ploidyLevels <- unique(rowData(out)$ploidy)
  }

  # update sample counts in the event of sample filtering
  if (!missing(j)) {
    metadata(out)$m <- ncol(out)
  }

  validObject(out)
  out
})
