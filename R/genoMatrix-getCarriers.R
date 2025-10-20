#' @rdname genoMatrix
#' @usage NULL
#' @export
setMethod(
  "getCarriers",
  signature = "genoMatrix",
  definition = function(
    object,
    VAR_id = NULL,
    colDataFields = NULL,
    rowDataFields = NULL,
    groupBy = NULL,
    aggregate = FALSE,
    imputeMethod = "meanImpute"
  ) {
    # validate input
    .getCarriers_validate_input(as.list(environment()))

    # subset VAR_ids (if specified)
    if (is.null(VAR_id)) {
      VAR_id <- rownames(object)
    } else {
      VAR_id <- as.character(VAR_id)
    }
    object <- object[VAR_id, ]

    # generate data.frame of individual carriers
    if (is.null(groupBy)) {
      # carrier indices
      indices <- which(
        !is.na(assays(object)$GT) & assays(object)$GT >= 1L,
        arr.ind = TRUE
      )

      if (nrow(indices) > 0L) {
        indices <- indices[
          order(indices[, "row"], indices[, "col"]),
          ,
          drop = FALSE
        ]
        carriers <- data.frame(
          VAR_id = rownames(object)[indices[, "row"]],
          IID = colnames(object)[indices[, "col"]],
          genotype = assays(object)$GT[indices],
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      } else {
        # return empty data.frame if zero carriers are found
        carriers <- data.frame(
          VAR_id = character(0),
          IID = character(0),
          genotype = numeric(0),
          stringsAsFactors = FALSE
        )
      }

      # add optional colDataFields
      if (!is.null(colDataFields)) {
        carriers <- dplyr::left_join(
          carriers,
          cbind(
            IID = colnames(object),
            as.data.frame(colData(object))[,
              setdiff(colDataFields, "IID"),
              drop = FALSE
            ]
          ),
          by = "IID"
        )
      }

      # add optional rowDataFields
      if (!is.null(rowDataFields)) {
        carriers <- dplyr::left_join(
          carriers,
          cbind(
            VAR_id = rownames(object),
            as.data.frame(rowData(object))[,
              setdiff(rowDataFields, "VAR_id"),
              drop = FALSE
            ]
          ),
          by = "VAR_id"
        )
      }

      # return carrier frequencies per group
    } else if (!is.null(groupBy) && !aggregate) {
      groupBy <- unique(c("VAR_id", groupBy))
      carriers <-
        cbind(
          tibble::as_tibble(colData(object)[,
            setdiff(groupBy, "VAR_id"),
            drop = FALSE
          ]),
          tibble::as_tibble(t(assays(object)$GT[VAR_id, , drop = FALSE]))
        ) %>%
        tidyr::pivot_longer(
          dplyr::all_of(VAR_id),
          names_to = "VAR_id",
          values_to = "genotype"
        ) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(groupBy))) %>%
        dplyr::summarize(
          carrierN = sum(genotype >= 1, na.rm = TRUE),
          carrierFreq = mean(genotype >= 1, na.rm = TRUE),
          carrierFreqSE = sd(genotype >= 1, na.rm = TRUE) / sqrt(dplyr::n()),
          .groups = "drop"
        )

      # return mean burden scores per group
    } else if (!is.null(groupBy) && aggregate) {
      object <- flipToMinor(object)
      object <- aggregate(
        recode(object, imputeMethod = imputeMethod),
        checkMissing = FALSE
      )

      carriers <- as.data.frame(colData(object))[,
        c("aggregate", groupBy),
        drop = FALSE
      ] %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(groupBy))) %>%
        dplyr::summarize(
          meanBurdenScore = mean(aggregate),
          .groups = "drop"
        )
    }

    carriers
  }
)

.getCarriers_validate_input <- function(args) {
  # general checks VAR_id and check if present in genoMatrix
  .gdb_check_varid(args[["VAR_id"]])
  if (
    !is.null(args[["VAR_id"]]) &&
      !all(args[["VAR_id"]] %in% rownames(args[["object"]]))
  ) {
    stop(
      "Not all specified VAR_ids are present in the genotype matrix.",
      call. = FALSE
    )
  }

  # type checks
  check_wrapper(check_character, args, "colDataFields", allow_null = TRUE)
  check_wrapper(check_character, args, "rowDataFields", allow_null = TRUE)
  check_wrapper(check_character, args, "groupBy", allow_null = TRUE)
  check_wrapper(check_bool, args, "aggregate", length_equal = 1L)
  check_wrapper(check_character, args, "imputeMethod", length_equal = 1L)

  # check if colDataFields are present in genoMatrix colData
  if (
    !is.null(args[["colDataFields"]]) &&
      !all(args[["colDataFields"]] %in% colnames(colData(args[["object"]])))
  ) {
    stop(
      "Not all specified colDataFields are present in colData(object)!",
      call. = FALSE
    )
  }

  # check if rowDataFields are present in genoMatrix rowData
  if (
    !is.null(args[["rowDataFields"]]) &&
      !all(args[["rowDataFields"]] %in% colnames(rowData(args[["object"]])))
  ) {
    stop(
      "Not all specified rowDataFields are present in rowData(object)!",
      call. = FALSE
    )
  }

  # check if groupBy fields are present in genoMatrix colData
  if (
    !is.null(args[["groupBy"]]) &&
      !all(args[["groupBy"]] %in% colnames(colData(args[["object"]])))
  ) {
    stop(
      "Not all specified `groupBy` fields are present in colData(object)!",
      call. = FALSE
    )
  }

  # if aggregate = TRUE, `groupBy` should be specified
  if (args[["aggregate"]] && is.null(args[["groupBy"]])) {
    stop(
      "The `groupBy` argument should be specified if `aggregate = TRUE`.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
