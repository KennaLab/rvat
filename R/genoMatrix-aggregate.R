#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod(
  "aggregate",
  signature = "genoMatrix",
  definition = function(x, returnGT = TRUE, checkMissing = TRUE) {
    # input validation
    check_bool(returnGT)
    check_bool(checkMissing)

    # stop if genoMatrix contains missing values
    if (checkMissing) {
      callRate <- getCR(x)
      if (any(callRate < 1.0, na.rm = TRUE)) {
        stop(
          "genoMatrix contains missing values, ",
          "impute missing values using the `recode` method.",
          call. = FALSE
        )
      }
    }

    # skip non-usable weights
    skip <- is.na(rowData(x)$w) | is.infinite(rowData(x)$w)
    if (all(skip)) {
      stop(
        "All weights are missing or infinite, ",
        "cannot compute aggregate.",
        call. = FALSE
      )
    } else if (any(skip)) {
      warning(
        sprintf(
          paste0(
            "%s/%s variants in the genoMatrix contain missing ",
            "or infinite weights, these are excluded."
          ),
          sum(skip),
          length(skip)
        ),
        call. = FALSE
      )
      x <- x[!skip, ]
    }

    # generate aggregate counts
    colData(x)$aggregate <- Matrix::crossprod(assays(x)$GT, rowData(x)$w)[, 1L]

    # return either genoMatrix or aggregate vector
    if (returnGT) {
      return(x)
    } else {
      return(colData(x)$aggregate)
    }
  }
)
