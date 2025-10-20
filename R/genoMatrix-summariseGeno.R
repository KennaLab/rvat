#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod(
  "summariseGeno",
  signature = "genoMatrix",
  definition = function(object) {
    if (metadata(object)$imputeMethod != "none") {
      stop("Please provide a non-imputed genoMatrix.", call. = FALSE)
    }

    # callRate
    callRate <- getCR(object)

    # ref,het,hom
    ref <- Matrix::rowSums(assays(object)$GT == 0, na.rm = TRUE)
    het <- Matrix::rowSums(assays(object)$GT == 1, na.rm = TRUE)
    hom <- Matrix::rowSums(assays(object)$GT == 2, na.rm = TRUE)

    # allele frequencies
    af <- getAF(object)

    # HWE
    hweP <- rep(NA_real_, nrow(object))
    meta <- metadata(object)

    ## diploid
    hweP[rowData(object)$ploidy == "diploid"] <- hweTest(
      ref[rowData(object)$ploidy == "diploid"],
      het[rowData(object)$ploidy == "diploid"],
      hom[rowData(object)$ploidy == "diploid"],
      af[rowData(object)$ploidy == "diploid"]
    )

    ## XnonPAR
    if ("XnonPAR" %in% meta$ploidyLevels) {
      if (any(object$sex == 0)) {
        warning(
          sprintf(
            paste0(
              "GT contains variants with ploidy='XnonPAR', ",
              "hweP is calculated within females, ",
              "note that sex is missing for %s samples."
            ),
            sum(object$sex == 0)
          ),
          call. = FALSE
        )
      }
      object_fem_XnonPAR <- object[
        rowData(object)$ploidy == "XnonPAR",
        colData(object)$sex == 2
      ]
      ref_X <- Matrix::rowSums(assays(object_fem_XnonPAR)$GT == 0, na.rm = TRUE)
      het_X <- Matrix::rowSums(assays(object_fem_XnonPAR)$GT == 1, na.rm = TRUE)
      hom_X <- Matrix::rowSums(assays(object_fem_XnonPAR)$GT == 2, na.rm = TRUE)
      hweP[rowData(object)$ploidy == "XnonPAR"] <-
        hweTest(
          ref = ref_X,
          het = het_X,
          hom = hom_X,
          af = getAF(object_fem_XnonPAR)
        )
    }

    ## YnonPAR
    if ("YnonPAR" %in% meta$ploidyLevels) {
      hweP[rowData(object)$ploidy == "YnonPAR"] <- 1.0
    }

    # return
    output <- data.frame(
      VAR_id = rownames(object),
      AF = af,
      callRate = callRate,
      geno0 = ref,
      geno1 = het,
      geno2 = hom,
      hweP = hweP,
      stringsAsFactors = FALSE
    )
    output
  }
)

hweTest <- function(ref, het, hom, af) {
  n <- ref + het + hom
  ref0 <- round(n * (1 - af)^2)
  het0 <- round(n * 2 * af * (1 - af))
  hom0 <- round(n * (af)^2)
  X <- ((ref - ref0)^2) /
    ref0 +
    ((het - het0)^2) / het0 +
    ((hom - hom0)^2) / hom0
  X[!(is.finite(X))] <- 0
  return(1 - pchisq(X, df = 1))
}
