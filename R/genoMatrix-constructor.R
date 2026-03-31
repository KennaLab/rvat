#' @rdname genoMatrix
#' @usage NULL
genoMatrix <- function(
  GT,
  SM,
  VAR_id,
  anno = NULL,
  w = 1.0,
  ploidy = "diploid",
  varSetName = "unnamed",
  unit = "unnamed",
  cohortname = "unnamed",
  genomeBuild = NA_character_,
  gdbpath = NA_character_,
  gdbid = NA_character_,
  verbose = TRUE
) {
  # validate input
  .genomatrix_validate_input(as.list(environment()))

  # set weights to 1 if not set
  if (length(w) == 1L) {
    w <- rep(w, length(VAR_id))
  }

  # set ploidy for all variants if not set
  if (length(ploidy) == 1L) {
    ploidy <- rep(ploidy, length(VAR_id))
  }

  # label GT matrix
  colnames(GT) <- SM$IID
  rownames(GT) <- VAR_id

  # remove samples with missing IID values
  missing_samples <- is.na(SM$IID)
  if (sum(missing_samples) > 0L) {
    if (verbose) {
      message(sprintf(
        "%s/%s samples in the gdb are present in cohort '%s'",
        sum(!missing_samples),
        nrow(SM),
        cohortname
      ))
    }
    GT <- GT[, !missing_samples, drop = FALSE]
    SM <- SM[!missing_samples, , drop = FALSE]
    rownames(SM) <- as.character(SM$IID)
  }

  # construct SummarizedExperiment
  GT <- SummarizedExperiment(
    assays = list(GT = GT),
    colData = DataFrame(SM),
    rowData = DataFrame(ploidy = ploidy, w = unname(w))
  )

  # set metadata
  metadata(GT)$gdb <- gdbpath
  metadata(GT)$gdbId <- gdbid
  metadata(GT)$genomeBuild <- genomeBuild
  metadata(GT)$ploidyLevels <- unique(as.character(ploidy))
  metadata(GT)$m <- ncol(GT)
  metadata(GT)$nvar <- nrow(GT)
  metadata(GT)$varSetName <- varSetName
  metadata(GT)$unit <- unit
  metadata(GT)$cohort <- cohortname
  metadata(GT)$geneticModel <- "allelic"
  metadata(GT)$imputeMethod <- "none"

  # create genoMatrix object
  GT <- new("genoMatrix", GT)

  # reset genotype dosages in accordance with male / female /
  # missing sex at sites with XnonPAR or YnonPAR ploidy
  GT <- .resetSexChromDosage(GT)

  # add optional row annotations
  if (!is.null(anno)) {
    GT <- updateGT(GT, anno = anno)
  }

  # return
  GT
}

.genomatrix_validate_input <- function(args) {
  # GT should be a matrix and same number of rows as `VAR_id`
  if (!is.matrix(args[["GT"]])) {
    stop("`GT` must be a matrix.", call. = FALSE)
  }
  # check if length of VAR_id equals number of rows in GT
  if (length(args[["VAR_id"]]) != nrow(args[["GT"]])) {
    stop(
      "The number of rows in `GT` doesn't equal the length of `VAR_id`.",
      call. = FALSE
    )
  }
  # VAR_id checks
  .gdb_check_varid(args[["VAR_id"]])

  # sample manifest should be a data.frame/DFrame with a non-duplicated IID column.
  # the number of samples should be equal to the number of columns in GT
  if (!is.data.frame(args[["SM"]]) && !is(args[["SM"]], "DFrame")) {
    stop("`SM` must be a data.frame or DFrame.", call. = FALSE)
  }
  if (!"IID" %in% colnames(args[["SM"]])) {
    stop("`SM` (sample manifest) must contain an 'IID' column.", call. = FALSE)
  }
  if (anyDuplicated(args[["SM"]]$IID[!is.na(args[["SM"]]$IID)]) != 0L) {
    stop("Duplicated IID values found in cohort.", call. = FALSE)
  }
  if (nrow(args[["SM"]]) != ncol(args[["GT"]])) {
    stop(
      sprintf(
        "Number of samples in `SM` (%d) doesn't match those in the genoMatrix (%d)",
        nrow(args[["SM"]]),
        ncol(args[["GT"]])
      ),
      call. = FALSE
    )
  }

  # anno should be a data.frame with a VAR_id column
  if (!is.null(args[["anno"]]) && !is.data.frame(args[["anno"]])) {
    stop("`anno` must be a data.frame", call. = FALSE)
  }
  if (!is.null(args[["anno"]]) && !"VAR_id" %in% colnames(args[["anno"]])) {
    stop("`anno` must contain a `VAR_id` column", call. = FALSE)
  }

  # weight vector should have same length as VAR_id vector
  if (!is.numeric(args[["w"]])) {
    stop("`w` (weights) must be a numeric vector or scalar.", call. = FALSE)
  }
  if (
    length(args[["w"]]) != 1L && length(args[["w"]]) != length(args[["VAR_id"]])
  ) {
    stop(
      sprintf(
        "Length of `w` (%d) must be 1 or match length of `VAR_id` (%d).",
        length(args[["w"]]),
        length(args[["VAR_id"]])
      ),
      call. = FALSE
    )
  }

  # ploidy vector should be of type character and have same length as VAR_id vector
  check_wrapper(check_character, args, "ploidy")
  if (
    length(args[["ploidy"]]) != 1L &&
      length(args[["ploidy"]]) != length(args[["VAR_id"]])
  ) {
    stop("`ploidy` vector doesn't match `VAR_id` vector.", call. = FALSE)
  }

  # type checks for remaining arguments
  check_wrapper(check_character, args, "varSetName", length_equal = 1L)
  check_wrapper(check_character, args, "unit", length_equal = 1L)
  check_wrapper(check_character, args, "cohortname", length_equal = 1L)
  check_wrapper(
    check_character,
    args,
    "genomeBuild",
    allow_na = TRUE,
    length_equal = 1L
  )
  check_wrapper(
    check_character,
    args,
    "gdbpath",
    allow_na = TRUE,
    length_equal = 1L
  )
  check_wrapper(
    check_character,
    args,
    "gdbid",
    allow_na = TRUE,
    length_equal = 1L
  )
  check_wrapper(check_bool, args, "verbose", length_equal = 1L)

  invisible(NULL)
}
