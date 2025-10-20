# Methods ---------------------------------------------------------------------

## reading & writing ----------------------------------------------------------

#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod(
  "writeResult",
  "rvatResult",
  function(
    object,
    file = "",
    append = FALSE,
    quote = FALSE,
    sep = "\t",
    eol = "\n",
    na = "NA",
    dec = ".",
    col.names = !append,
    qmethod = c("escape", "double"),
    fileEncoding = ""
  ) {
    if (append) {
      write.table(
        object,
        file = file,
        append = append,
        quote = quote,
        sep = sep,
        eol = eol,
        na = na,
        dec = dec,
        row.names = FALSE,
        col.names = col.names,
        qmethod = qmethod,
        fileEncoding = fileEncoding
      )
    } else {
      file <- file(file, "w")
      on.exit(close(file), add = TRUE)

      # write metadata
      metadata <- metadata(object)
      metadata <- metadata[names(metadata) %in% metadata_rvatresult]
      .write_rvat_header(
        filetype = as.character(class(object)[1L]),
        metadata = metadata,
        con = file
      )

      # write results
      write.table(
        object,
        file = file,
        append = FALSE,
        quote = quote,
        sep = sep,
        eol = eol,
        na = na,
        dec = dec,
        row.names = FALSE,
        col.names = col.names,
        qmethod = qmethod,
        fileEncoding = fileEncoding
      )
    }
  }
)


#' Read association results
#'
#' Read results generated using the \code{\link{assocTest}} method.

#' @param path File path to data
#' @param type Result type ('singlevarResult', 'rvbResult', 'gsaResult').
#' Defaults to `NULL` in which case the result type is inferred from the header.
#' @param sep The field separator. Defaults to `\\t`, which is the default separator using in \code{\link{assocTest}}.
#' @return An object of type \code{\link{rvbResult}} or \code{\link{singlevarResult}}.
#' @export
readResults <- function(path, type = NULL, sep = "\t") {
  # input validation
  check_character(path, allow_null = TRUE)
  check_length(path, equal = 1L, allow_null = TRUE)
  check_character(type, allow_null = TRUE)
  check_length(type, equal = 1L, allow_null = TRUE)
  check_character(sep)
  check_length(sep, equal = 1L)

  # parse header
  header_info <- .parse_rvat_header(
    path,
    expected_metadata = metadata_rvatresult,
    expected_filetype = if (!is.null(type)) type else names(columns_rvatResult),
    return_skip = TRUE,
    n = length(metadata_rvatresult) + 1L
  )

  # determine filetype, either from header or as provided by user
  header_filetype <- header_info$filetype
  if (is.null(header_filetype) && is.null(type)) {
    # try to infer from colnames
    type <- checkClassrvatResult(path, sep = sep)
  }
  filetype <- if (!is.null(type)) type else header_filetype

  # read results
  dat <- read.table(
    file = path,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    comment.char = "#"
  )

  # check if all expected fields are present
  expected_fields <- names(columns_rvatResult[[filetype]])
  if (!all(expected_fields %in% colnames(dat))) {
    missing_cols <- setdiff(expected_fields, colnames(dat))
    stop(
      sprintf(
        paste0(
          "The results file is missing the following ",
          "required columns for type '%s': %s"
        ),
        filetype,
        paste(sQuote(missing_cols), collapse = ",")
      ),
      call. = FALSE
    )
  }

  # convert to rvatResult and add metadata
  result <- rvatResult(dat, class = filetype)
  metadata(result) <- header_info$metadata

  # return rvatResult
  result
}

## getters ----------------------------------------------------------
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("getGdbId", signature = "rvatResult", definition = function(object) {
  value <- metadata(object)$gdbId
  if (length(value) == 0L) value <- NA_character_
  value
})

setMethod("getIdCol", signature = "rvatResult", definition = function(object) {
  result_class <- checkClassrvatResult(object)
  id_col <- switch(
    result_class,
    rvbResult = "unit",
    singlevarResult = "VAR_id",
    gsaResult = "geneSetName",
    stop(
      "Unknown result type.",
      call. = FALSE
    )
  )
  id_col
})

#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod(
  "getGenomeBuild",
  signature = "rvatResult",
  definition = function(object) {
    value <- metadata(object)$genomeBuild
    if (length(value) == 0L) value <- NA_character_
    value
  }
)

## combine -------------------------------------------------------------------
#' merge.rvatResult.data.frame
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("merge", c("rvatResult", "data.frame"), function(x, y, by, ...) {
  result_class <- checkClassrvatResult(x)
  nrows <- nrow(x)
  metadata <- metadata(x)
  x <- dplyr::left_join(as.data.frame(x), y, by = by)

  if (nrow(x) != nrows) {
    message("Rows have been added to the rvatResult object")
  }

  x <- rvatResult(x, class = result_class)
  metadata(x) <- metadata
  x
})

#' merge.rvatResult.DataFrame
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("merge", c("rvatResult", "DataFrame"), function(x, y, by, ...) {
  result_class <- checkClassrvatResult(x)
  nrows <- nrow(x)
  metadata <- metadata(x)
  x <- dplyr::left_join(as.data.frame(x), y, by = by)

  if (nrow(x) != nrows) {
    message("Rows have been added to the rvatResult object")
  }

  x <- rvatResult(x, class = result_class)
  metadata(x) <- metadata
  x
})

## misc -----------------------------------------------------------------------

# topResult
### topResults --------------------------------------------------------------

#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("topResult", "rvatResult", function(object, n = 10L) {
  head(object[order(object$P), ], n)
})

# Check child class, for internal use
checkClassrvatResult <- function(object, type = NULL, sep = "\t") {
  # input validation
  if (!is.null(type) && !type %in% names(columns_rvatResult)) {
    stop(
      sprintf(
        "Type should be one of the following: %s.",
        paste(paste0("'", names(columns_rvatResult), "'"), collapse = ",")
      ),
      call. = FALSE
    )
  }
  check_character(sep)
  check_length(sep, equal = 1L)

  # if rvatResult is provided, directly return type
  if (is(object, "rvatResult")) {
    type_check <- vapply(
      names(columns_rvatResult),
      FUN = function(x) {
        is(object, x)
      },
      FUN.VALUE = logical(1L)
    )
    return(names(type_check)[type_check])
  }

  # if object is not an rvatResult, it should be a filepath
  if (!is(object, "character") || length(object) != 1L) {
    stop(
      "Input must be an rvatResult object or a single file path.",
      call. = FALSE
    )
  }

  # read header and check overlap with rvatResult headers
  cols <- colnames(read.table(
    object,
    header = TRUE,
    nrows = 1L,
    sep = sep,
    comment.char = "#"
  ))
  overlap <- vapply(
    names(columns_rvatResult),
    FUN = function(result_type) {
      all(names(columns_rvatResult[[result_type]]) %in% cols)
    },
    FUN.VALUE = logical(1L)
  )

  # determine result class
  result_class <- if (!is.null(type)) {
    # if provided by user
    type
  } else if (sum(overlap) == 1L) {
    # otherwise infer it from overlap
    names(overlap)[overlap]
  } else {
    NA_character_
  }

  # check validity of result type
  if (is.na(result_class) && sum(overlap) > 1L) {
    stop(
      sprintf(
        paste0(
          "Cannot distinguish result type. Columns match multiple types: ",
          "%s. Please specify the `type` argument."
        ),
        paste(sQuote(names(overlap)[overlap]), collapse = ", ")
      ),
      call. = FALSE
    )
  } else if (is.na(result_class)) {
    stop(
      "Could not infer result type from file columns. ",
      "Please specify the `type` argument.",
      call. = FALSE
    )
  } else if (!result_class %in% names(columns_rvatResult)) {
    stop(
      sprintf(
        "Specified `type` ('%s') is not a recognized result type.",
        result_class
      ),
      call. = FALSE
    )
  }

  # check if expected names
  expected_names <- names(columns_rvatResult[[result_class]])
  if (!all(expected_names %in% cols)) {
    missing_cols <- setdiff(expected_names, cols)
    stop(
      sprintf(
        paste0(
          "The file is missing required columns for the ",
          "specified or inferred type '%s': %s"
        ),
        result_class,
        paste(sQuote(missing_cols), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  result_class
}

# rvatResult constructor, for internal use
rvatResult <- function(object, class) {
  switch(
    class,
    rvbResult = rvbResult(object),
    singlevarResult = singlevarResult(object),
    gsaResult = gsaResult(object),
    stop(
      sprintf(
        "Class '%s' is not a known rvatResult type.",
        class
      ),
      call. = FALSE
    )
  )
}

## plots -----------------------------------------------------------------------
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod(
  "qqplot",
  "rvatResult",
  function(
    object,
    title = "",
    label = NULL,
    threshold = NULL,
    showThreshold = TRUE,
    labelThreshold = NULL,
    cex = 16,
    lambda1000 = FALSE,
    case = NULL,
    control = NULL,
    verbose = TRUE
  ) {
    # validate input
    .qqplot_validate_inputs(as.list(environment()))

    # get caseN/ctrlN from results if not specified
    if (lambda1000 && (is.null(case) || is.null(control))) {
      case <- max(object$caseN, na.rm = TRUE)
      control <- max(object$ctrlN, na.rm = TRUE)
    }

    # subset required fields
    fields <- c("P", label)
    qq_data <- as.data.frame(object)[, fields, drop = FALSE]

    # generate qqplot
    .qqplot(
      qq_data = qq_data,
      title = title,
      label = label,
      threshold = threshold,
      showThreshold = showThreshold,
      labelThreshold = labelThreshold,
      cex = cex,
      lambda1000 = lambda1000,
      case = case,
      control = control
    )
  }
)

.qqplot_validate_inputs <- function(args) {
  # type checks
  check_wrapper(check_character, args, "title", length_equal = 1L)
  check_wrapper(
    check_character,
    args,
    "label",
    length_equal = 1L,
    allow_null = TRUE
  )
  check_wrapper(
    check_number_decimal,
    args,
    "threshold",
    length_equal = 1L,
    allow_null = TRUE
  )
  check_wrapper(check_bool, args, "showThreshold", length_equal = 1L)
  check_wrapper(
    check_number_decimal,
    args,
    "labelThreshold",
    length_equal = 1L,
    allow_null = TRUE
  )
  check_wrapper(check_number_decimal, args, "cex", length_equal = 1L)
  check_wrapper(check_bool, args, "lambda1000", length_equal = 1L)
  check_wrapper(
    check_number_whole,
    args,
    "case",
    length_equal = 1L,
    allow_null = TRUE
  )
  check_wrapper(
    check_number_whole,
    args,
    "control",
    length_equal = 1L,
    allow_null = TRUE
  )

  # check if label exists
  if (
    !is.null(args[["label"]]) &&
      !args[["label"]] %in% colnames(args[["object"]])
  ) {
    stop(
      sprintf("Label '%s' not found in object.", args[["label"]]),
      call. = FALSE
    )
  }

  # lambda1000 doesn't apply to gsaResult
  if (
    args[["lambda1000"]] &&
      (checkClassrvatResult(args[["object"]]) == "gsaResult")
  ) {
    stop(
      "lambda1000 does not apply to a gsaResult object.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

.qqplot <- function(
  qq_data,
  title = "",
  label = "label",
  threshold = NULL,
  showThreshold = TRUE,
  labelThreshold = NULL,
  cex = 16,
  lambda1000 = FALSE,
  case = NULL,
  control = NULL,
  verbose = TRUE
) {
  # set significance threshold to bonferroni if not specified
  if (is.null(threshold)) {
    threshold <- -log10(0.05 / nrow(qq_data))
  } else {
    threshold <- -log10(threshold)
  }

  # threshold for showing labels is identical to significant threshold,
  # unless specified otherwise
  if (is.null(labelThreshold)) {
    labelThreshold <- threshold
  } else {
    labelThreshold <- -log10(labelThreshold)
  }

  # subset records with missing P-values
  qq_data <- .qqplot_prepare_data(qq_data, verbose = verbose)
  lambda <- .get_lambda(qq_data$P)
  lambda1000_value <- if (lambda1000)
    .get_lambda1000(lambda, case, control) else NULL
  anno <- .qqplot_annotation(
    lambda_value = lambda,
    lambda1000_value = lambda1000_value,
    show_lambda1000 = lambda1000
  )

  # main plot
  mx <- max(qq_data$null)
  qq <- ggplot2::ggplot(qq_data, ggplot2::aes(x = null, y = logp)) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(expression("-log"[10] * "(Null)")) +
    ggplot2::ylab(expression("-log"[10] * "(Obs)")) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "black"),
      axis.line.y = ggplot2::element_line(color = "black"),
      axis.text.x = ggplot2::element_text(size = cex),
      axis.text.y = ggplot2::element_text(size = cex),
      axis.title.x = ggplot2::element_text(size = cex + 2),
      axis.title.y = ggplot2::element_text(size = cex + 2)
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      col = "slateblue"
    ) +
    ggplot2::geom_text(
      data = anno,
      ggplot2::aes(
        x = xpos,
        y = ypos,
        hjust = hjust,
        vjust = vjust,
        label = annotateText
      ),
      parse = FALSE,
      size = cex / 2.5
    )

  # add significance threshold
  if (showThreshold) {
    qq <- qq +
      ggplot2::geom_segment(
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        col = "red",
        linetype = 2,
        data = data.frame(x = 0, xend = mx, y = threshold, yend = threshold)
      )
  }

  # add label
  if (!is.null(label)) {
    qq <- qq +
      ggrepel::geom_text_repel(
        data = (qq_data[
          !is.na(qq_data$logp) & qq_data$logp > labelThreshold,
          ,
          drop = FALSE
        ]),
        ggplot2::aes(label = .data[[label]]),
        min.segment.length = 0,
        box.padding = 0.5
      )
  }

  qq
}

.qqplot_prepare_data <- function(qq_data, verbose) {
  if (anyNA(qq_data$P)) {
    if (verbose) {
      message(sprintf(
        "%d rows with missing P-values are excluded.",
        sum(is.na(qq_data$P))
      ))
    }
    qq_data <- qq_data[!is.na(qq_data$P), , drop = FALSE]
  }

  qq_data$logp <- -log10(qq_data$P)
  qq_data <- qq_data[order(qq_data$logp), ]
  qq_data$null <- sort(-log10(runif(nrow(qq_data))))

  qq_data
}

.qqplot_annotation <- function(
  lambda_value,
  lambda1000_value,
  show_lambda1000
) {
  anno <- data.frame(
    xpos = Inf,
    ypos = -Inf,
    annotateText = if (show_lambda1000) {
      c(
        sprintf("\u03BB = %s", signif(lambda_value, 5L)),
        sprintf("\u03BB1000 = %s", signif(lambda1000_value, 5L))
      )
    } else {
      sprintf("\u03BB = %s", signif(lambda_value, 5L))
    },
    hjust = 1,
    vjust = if (show_lambda1000) c(-3, -1) else -1
  )
  anno
}

.get_lambda <- function(P) {
  chisq <- qchisq(P, 1L, lower.tail = FALSE)
  lambda <- median(chisq) / qchisq(0.5, 1L)
  lambda
}

.get_lambda1000 <- function(lambda, case, control) {
  lambda1000 <- 1.0 + (lambda - 1.0) * (1.0 / case + 1.0 / control) * 500.0
  lambda1000
}

#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod(
  "manhattan",
  "rvatResult",
  function(
    object,
    highlight = NULL,
    label = NULL,
    threshold = NULL,
    labelThreshold = NULL,
    labelRepel = FALSE,
    labelSize = 3.88,
    contigs = NULL,
    title = "",
    verbose = TRUE
  ) {
    # validate input
    .manhattan_validate_input(as.list(environment()))

    # contigs
    contigs_df <- .manhattan_handle_contigs(object, contigs = contigs)

    # set label
    if (is.null(label) && "label" %in% colnames(object)) {
      label <- "label"
    }

    # configuration
    man_data <- .manhattan_prepare_data(
      object,
      contigs_df = contigs_df,
      label = label,
      highlight = highlight,
      verbose = verbose
    )

    if (is.null(threshold)) {
      threshold <- 0.05 / nrow(man_data)
    }
    if (is.null(labelThreshold)) {
      labelThreshold <- threshold
    }

    # Make + filter dataframe with points that need to be highlighted
    mplot <- ggplot2::ggplot(
      man_data,
      ggplot2::aes(x = POS_cumulative, y = logp, color = code)
    ) +
      ggplot2::xlab("") +
      ggplot2::ylab(expression("-Log"[10] * "(p)")) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line(color = "black"),
        axis.line.y = ggplot2::element_line(color = "black"),
        axis.text.x = ggplot2::element_text(size = 14, angle = 90),
        axis.text.y = ggplot2::element_text(size = 16),
        axis.title.y = ggplot2::element_text(size = 18)
      ) +
      ggplot2::scale_x_continuous(
        limits = c(min(man_data$POS_cumulative), max(man_data$POS_cumulative)),
        expand = c(0, 0),
        breaks = (contigs_df$increment + contigs_df$Length / 2),
        labels = c(1:22, "X", "Y")
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, max(man_data$logp) + 0.75),
        expand = c(0, 0)
      ) +
      ggplot2::guides(colour = "none") +
      ggplot2::geom_segment(
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        col = "grey",
        linetype = 2,
        data = data.frame(
          x = min(man_data$POS_cumulative),
          xend = max(man_data$POS_cumulative),
          y = -log10(threshold),
          yend = -log10(threshold)
        )
      ) +
      ggplot2::geom_point()

    # add labels if specified
    if (!is.null(label)) {
      if (labelRepel) {
        mplot <- mplot +
          ggrepel::geom_text_repel(
            ggplot2::aes(x = POS_cumulative, y = logp, label = .data[[label]]),
            data = (man_data[
              !is.na(man_data$logp) & man_data$P < labelThreshold,
            ]),
            size = labelSize
          )
      } else {
        mplot <- mplot +
          ggplot2::geom_text(
            ggplot2::aes(x = POS_cumulative, y = logp, label = .data[[label]]),
            data = (man_data[
              !is.na(man_data$logp) & man_data$P < labelThreshold,
            ]),
            nudge_y = 0.7,
            angle = 0,
            vjust = "inward",
            hjust = "inward",
            size = labelSize
          )
      }
    }

    # add highlight if specified
    if (!is.null(highlight)) {
      mplot <- mplot +
        ggplot2::geom_point(
          data = man_data[man_data$highlight, , drop = FALSE],
          ggplot2::aes(x = POS_cumulative, y = logp),
          color = "red"
        )
    }

    # add title if specified
    if (title != "") {
      mplot <- mplot + ggplot2::ggtitle(title)
    }

    # return plot
    mplot
  }
)

.manhattan_validate_input <- function(args) {
  # type checks
  check_wrapper(check_character, args, "highlight", allow_null = TRUE)
  check_wrapper(
    check_character,
    args,
    "label",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(
    check_number_decimal,
    args,
    "threshold",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(
    check_number_decimal,
    args,
    "labelThreshold",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(check_bool, args, "labelRepel", length_equal = 1L)
  check_wrapper(check_number_decimal, args, "labelSize", length_equal = 1L)
  check_wrapper(check_character, args, "title", length_equal = 1L)

  # check contigs
  if (
    !is.null(args[["contigs"]]) &&
      !is.character(args[["contigs"]]) &&
      !is.data.frame(args[["contigs"]])
  ) {
    stop(
      "`contigs` must be a character string (e.g., 'GRCh38') or a data.frame.",
      call. = FALSE
    )
  }
  if (
    is.data.frame(args[["contigs"]]) &&
      !all(c("CHROM", "Length") %in% colnames(args[["contigs"]]))
  ) {
    stop(
      "A custom `contigs` data.frame must contain ",
      "'CHROM' and 'Length' columns.",
      call. = FALSE
    )
  }

  # validate results
  cols <- colnames(args[["object"]])
  if (!"P" %in% cols) {
    stop(
      "Results should contain a 'P' column for P-values.",
      call. = FALSE
    )
  }

  if (!"CHROM" %in% cols) {
    stop("Results must contain a 'CHROM' column.", call. = FALSE)
  }

  # check chrom,pos fields
  if (!"POS" %in% cols && !all(c("start", "end") %in% cols)) {
    stop(
      "Results should contain either a 'POS' column or both ",
      "'start' and 'end' columns.",
      call. = FALSE
    )
  }

  # check if the specified label column exists
  if (!is.null(args[["label"]]) && !args[["label"]] %in% cols) {
    stop(
      sprintf(
        "Specified `label` column '%s' not found in results.",
        args[["label"]]
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

.manhattan_handle_contigs <- function(object, contigs) {
  # if contigs not specified, try to get from results metadata
  # otherwise set to default (GRCh37) with a warning.
  if (is.null(contigs)) {
    results_build <- getGenomeBuild(object)
    if (
      is.na(results_build) || !nzchar(results_build) || is.null(results_build)
    ) {
      warning("Contigs not specified, defaulting to GRCh37.", call. = FALSE)
      contigs <- "GRCh37"
    } else {
      contigs <- results_build
    }
  }

  # get data.frame of contigs based on specified contigs
  if (is.character(contigs) || length(contigs) == 1L) {
    contigs_df <- .build_contigs[[contigs]]
  } else if (is.data.frame(contigs)) {
    contigs_df <- contigs
  }

  # add increment
  contigs_df$increment <- c(
    0L,
    cumsum(as.numeric(contigs_df$Length))[-nrow(contigs_df)]
  )

  # recode chromosomes
  contigs_df$CHROM <- .convert_chrom_to_num(contigs_df$CHROM)

  # return
  contigs_df
}

.manhattan_prepare_data <- function(
  object,
  contigs_df,
  label,
  highlight,
  verbose
) {
  # data.frame
  plot_data <- as.data.frame(object)

  # generate POS based on start/end
  if (!"POS" %in% colnames(plot_data)) {
    plot_data$POS <- (plot_data$start + plot_data$end) / 2.0
  }

  # standardize chromosome names
  plot_data$CHROM <- .convert_chrom_to_num(plot_data$CHROM)
  contigs_df$CHROM <- .convert_chrom_to_num(contigs_df$CHROM)

  # remove records with missing P-values or which are not present in contigs
  missing_P <- is.na(plot_data$P)
  not_in_contig <- !(plot_data$CHROM %in% contigs_df$CHROM)

  if (any(missing_P) && verbose) {
    message(sprintf(
      "  Removing %d row(s) with missing P-values.",
      sum(missing_P)
    ))
  }
  if (any(not_in_contig) && verbose) {
    message(sprintf(
      "  Removing %d row(s) with non-standard chromosomes.",
      sum(not_in_contig & !missing_P)
    ))
  }
  plot_data <- plot_data[!missing_P & !not_in_contig, ]

  if (nrow(plot_data) == 0L) return(plot_data) # Return empty if no data remains

  # add contigs
  plot_data <- dplyr::left_join(
    plot_data,
    contigs_df,
    by = "CHROM"
  )

  # add cumulative position
  plot_data$POS_cumulative <- plot_data$POS + plot_data$increment

  # add -log10 P-value
  plot_data$logp <- -log10(plot_data$P)

  # add code for alternating chromosome colors
  plot_data$code <- plot_data$CHROM %% 2L

  # add highlight field if specified
  if (!is.null(highlight)) {
    id_col <- getIdCol(object)
    plot_data$highlight <- plot_data[[id_col]] %in% highlight
    highlight_col <- "highlight"
  } else {
    highlight_col <- NULL
  }

  # keep only necessary columns
  cols_keep <- c("POS_cumulative", "logp", "P", "code", label, highlight_col)
  plot_data[, cols_keep, drop = FALSE]
}

#' @rdname densityPlot
#' @usage NULL
#' @export
setMethod(
  "densityPlot",
  signature = signature(object = "rvatResult"),
  function(
    object,
    geneSet,
    geneSetList,
    showMeans = FALSE,
    INT = FALSE,
    Zcutoffs = NULL,
    title = ""
  ) {
    if (checkClassrvatResult(object) == "gsaResult") {
      stop(
        "gsaResult objects do not have per-test P-values, only per gene set.",
        "For 'object', select an rvbResult or singlevarResult object.",
        call. = FALSE
      )
    }

    object <- .prepare_stats_GSA(
      object,
      INT = INT,
      Zcutoffs = Zcutoffs,
      covar = NULL
    )
    units <- listUnits(getGeneSet(geneSetList, geneSet = geneSet))
    object$in_geneset <- object$unit %in% units

    if (showMeans) {
      means <- as.data.frame(object) %>%
        dplyr::group_by(in_geneset) %>%
        dplyr::summarize(mean = mean(Z, na.rm = TRUE))
    }

    # https://github.com/clauswilke/colorblindr/blob/master/R/palettes.R
    colorBlindFriendly <- c(
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7",
      "#999999"
    )

    ggplot2::ggplot(
      as.data.frame(object),
      ggplot2::aes(x = Z, color = in_geneset, fill = in_geneset)
    ) +
      ggplot2::ggtitle(title) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::geom_density(alpha = 0.6) +
      {
        if (showMeans) {
          ggplot2::geom_vline(
            data = means,
            mapping = ggplot2::aes(xintercept = mean)
          )
        }
      } +
      ggplot2::theme(text = ggplot2::element_text(size = 12)) +
      ggplot2::scale_colour_manual(values = colorBlindFriendly) +
      ggplot2::scale_fill_manual(values = colorBlindFriendly) +
      ggplot2::theme_classic()
  }
)

# rvbResult --------------------------------------------------------------------

# Constructor
proto_rvbResult <- function(n = 0L) {
  object <- S4Vectors::DataFrame(
    unit = rep(NA_character_, n),
    cohort = Rle(rep(NA_character_, n)),
    varSetName = Rle(rep(NA_character_, n)),
    name = Rle(rep(NA_character_, n)),
    pheno = Rle(rep(NA_character_, n)),
    covar = Rle(rep(NA_character_, n)),
    geneticModel = Rle(rep(NA_character_, n)),
    MAFweight = Rle(rep(NA_character_, n)),
    test = Rle(rep(NA_character_, n)),
    nvar = rep(NA_real_, n),
    caseCarriers = rep(NA_real_, n),
    ctrlCarriers = rep(NA_real_, n),
    meanCaseScore = rep(NA_real_, n),
    meanCtrlScore = rep(NA_real_, n),
    caseN = rep(NA_real_, n),
    ctrlN = rep(NA_real_, n),
    caseCallRate = rep(NA_real_, n),
    ctrlCallRate = rep(NA_real_, n),
    effect = rep(NA_real_, n),
    effectSE = rep(NA_real_, n),
    effectCIlower = rep(NA_real_, n),
    effectCIupper = rep(NA_real_, n),
    OR = rep(NA_real_, n),
    P = rep(NA_real_, n)
  )
  object
}


#' rvbResult
#' @rdname rvatResult
#' @usage NULL
#' @export
rvbResult <- function(object) {
  # return an empty rvbResult if input is empty
  if (missing(object) || is.null(object)) {
    return(new("rvbResult", proto_rvbResult()))
  }

  # object should be a filepath or a data.frame/DFrame
  if (
    !is.character(object) && !is.data.frame(object) && !is(object, "DFrame")
  ) {
    stop(
      "`object` should be either a data.frame/DataFrame or a filepath.",
      call. = FALSE
    )
  }

  # assume filepath if object if of type character
  if (is.character(object)) {
    return(readResults(path = object, type = "rvbResult", sep = "\t"))
  }

  # format and fill in missing fields if object is a data.frame/DFrame
  if (is.data.frame(object) || is(object, "DFrame")) {
    object <- .rvatResult_format(
      object,
      id_col = "unit",
      columns_info = columns_rvbResults,
      proto_func = proto_rvbResult,
      rle_cols = columns_rvbResults_rle,
      num_cols = columns_rvbResults_numeric
    )

    # return rvbResult object
    new("rvbResult", object)
  }
}

.rvatResult_format <- function(
  object,
  id_col,
  columns_info,
  proto_func,
  rle_cols,
  num_cols
) {
  # otherwise it should be a data.frame/DFrame with at least unit and P columns
  if (!all(c(id_col, "P") %in% colnames(object))) {
    stop(
      sprintf("Object must have at least '%s' and 'P' columns.", id_col),
      call. = FALSE
    )
  }

  if (all(names(columns_info) %in% colnames(object))) {
    object <- cbind(
      object[, names(columns_info)],
      object[, !colnames(object) %in% names(columns_info), drop = FALSE]
    )
  } else {
    # fill in missing fields
    proto <- proto_func(n = nrow(object))
    proto[, intersect(
      colnames(proto),
      colnames(object)
    )] <- object[, intersect(colnames(proto), colnames(object))]
    object <- cbind(
      proto,
      object[, !colnames(object) %in% colnames(proto), drop = FALSE]
    )
  }

  # set correct types
  if (is.data.frame(object)) object <- S4Vectors::DataFrame(object)
  object[, rle_cols] <- lapply(
    object[, rle_cols],
    S4Vectors::Rle
  )
  object[, num_cols] <- lapply(
    object[, num_cols],
    as.numeric
  )
  object[[id_col]] <- as.character(object[[id_col]])

  object
}

# validity
setValidity2("rvbResult", function(object) {
  .rvatResult_setvalidity(object, columns_info = columns_rvbResults)
})

.rvatResult_setvalidity <- function(object, columns_info) {
  missing_cols <- names(columns_info)[
    !names(columns_info) %in% colnames(object)
  ]
  if (ncol(object) != 1L && length(missing_cols) > 0L) {
    msg <-
      sprintf(
        "The following columns are missing: %s",
        paste(missing_cols, collapse = ",")
      )
    return(msg)
  }

  if (ncol(object) != 1L) {
    columns <- names(columns_info)[
      names(columns_info) %in% colnames(object)
    ]
    types <- unlist(lapply(object@listData, class))
    type_check <- unlist(lapply(
      columns,
      FUN = function(x) types[x] %in% columns_info[[x]]
    ))
    names(type_check) <- columns
    if (!all(type_check)) {
      type_msg <- lapply(
        names(type_check)[!type_check],
        FUN = function(x) {
          sprintf(
            "%s: %s",
            x,
            paste(columns_info[[x]], collapse = " or ")
          )
        }
      )

      msg <- sprintf(
        "The following columns should be of type:\n%s",
        paste(type_msg, collapse = "\n")
      )
      return(msg)
    }
  }

  # return TRUE if all checks are passed
  TRUE
}


## Methods ---------------------------------------------------------------------

### summary --------------------------------------------------------------------

#' summarize rvbResult
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("summary", "rvbResult", function(object, asList = FALSE, ...) {
  # generate summary
  info <- list()
  info[["Nunits"]] <- length(unique(object[["unit"]]))

  info[["ctrlNmin"]] <- if ((all(is.na(object[["ctrlN"]])))) {
    NA_real_
  } else if (nrow(object) > 0L) {
    min(object[["ctrlN"]], na.rm = TRUE)
  } else {
    0L
  }
  info[["ctrlNmax"]] <- if ((all(is.na(object[["ctrlN"]])))) {
    NA_real_
  } else if (nrow(object) > 0L) {
    max(object[["ctrlN"]], na.rm = TRUE)
  } else {
    0L
  }
  info[["caseNmin"]] <- if ((all(is.na(object[["caseN"]])))) {
    NA_real_
  } else if (nrow(object) > 0L) {
    min(object[["caseN"]], na.rm = TRUE)
  } else {
    0L
  }
  info[["caseNmax"]] <- if ((all(is.na(object[["caseN"]])))) {
    NA_real_
  } else if (nrow(object) > 0L) {
    max(object[["caseN"]], na.rm = TRUE)
  } else {
    0L
  }
  info[["name"]] <- unique(object[["name"]])
  info[["varSetName"]] <- unique(object[["varSetName"]])
  info[["test"]] <- unique(object[["test"]])
  info[["covar"]] <- unique(object[["covar"]])
  info[["geneticModel"]] <- unique(object[["geneticModel"]])
  info[["MAFweight"]] <- unique(object[["MAFweight"]])
  info[["pheno"]] <- unique(object[["pheno"]])
  info[["cohort"]] <- unique(object[["cohort"]])

  # return as is if asList=TRUE
  if (asList) {
    return(info)
  }

  # otherwise print summary
  cat(sprintf("%s object\n-------------------\n", class(object)[1L]))
  cat(sprintf(
    paste(
      "n units = %s",
      "N cases = %s",
      "N controls = %s",
      "names = %s",
      "varSets = %s",
      "tests = %s",
      "covars = %s",
      "geneticModels = %s",
      "MAFweights = %s",
      "pheno = %s",
      "cohorts = %s",
      sep = "\n"
    ),
    info$Nunits,
    info$caseNmax,
    info$ctrlNmax,
    .summary_truncate(info$name),
    .summary_truncate(info$varSetName),
    .summary_truncate(info$test),
    .summary_truncate(info$covar),
    .summary_truncate(info$geneticModel),
    .summary_truncate(info$MAFweight),
    .summary_truncate(info$pheno),
    .summary_truncate(info$cohort)
  ))
  cat("\n")
})

.summary_truncate <- function(string, max_items = 5L, sep = "; ") {
  # return none if empty
  if (length(string) == 0L || all(is.na(string))) {
    return("none")
  }

  # if the doesn't exceed max items, simply collapse and return
  if (length(string) <= max_items) {
    return(paste(string, collapse = sep))
  } else {
    # if does exceed max items, trunacate
    return(sprintf(
      "%s ...+ %d more",
      paste(string[1L:max_items], collapse = sep),
      length(string) - max_items
    ))
  }
}

### ACAT -----------------------------------------------------------------------

#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod(
  "ACAT",
  "rvatResult",
  function(
    object,
    aggregate = "test",
    group = c(
      "unit",
      "cohort",
      "varSetName",
      "name",
      "pheno",
      "covar",
      "geneticModel",
      "MAFweight",
      "test"
    ),
    fixpval = TRUE,
    fixpval_method = c("minmax", "manual", "Liu"),
    fixpval_maxP = 0.99,
    fixpval_minP = 1e-32,
    warning = TRUE
  ) {
    # validate input
    .ACAT_validate_input(as.list(environment()))
    fixpval_method <- match.arg(fixpval_method)
    if (is.character(aggregate)) aggregate <- list(aggregate)
    aggregate_all <- unlist(aggregate)

    # group should contain all aggregate columns
    if (!all(aggregate_all %in% group)) {
      group <- c(group, aggregate_all[!aggregate_all %in% group])
    }

    # remove missing P-values
    if (sum(is.na(object[["P"]])) > 0L) {
      if (warning) {
        message(sprintf(
          "Removing %s missing P-values.",
          sum(is.na(object[["P"]]))
        ))
      }
      object <- object[!is.na(object[["P"]]), ]
    }

    # perform ACAT per aggregate group
    metadata <- metadata(object)
    for (agg in aggregate) {
      group <- group[!group %in% agg]

      # add this round of ACAT to metadata
      metadata <- .ACAT_update_metadata(
        object,
        metadata = metadata,
        agg = agg,
        group = group,
        fixpval = fixpval,
        fixpval_method = fixpval_method,
        fixpval_maxP = fixpval_maxP,
        fixpval_minP = fixpval_minP
      )

      # fix pvalues that are exactly 0 or 1
      dat <- .ACAT_fix_pvals(
        object,
        group = group,
        fixpval = fixpval,
        fixpval_method = fixpval_method,
        fixpval_maxP = fixpval_maxP,
        fixpval_minP = fixpval_minP,
        warning = warning
      )

      # check if all groups are of length 1,
      # which usually indicates a misspecification
      count <- dat %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>%
        dplyr::count()
      if (all(count$n == 1L)) {
        warning("All groupings have length 1.", call. = FALSE)
      }

      # perform ACAT
      acat <- dat %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>%
        dplyr::summarize(P = .rvat_ACAT(P), .groups = "drop")

      # add ACAT P-values + summarise other columns
      object <- .ACAT_parse_results(
        object,
        acat = acat,
        group = group,
        agg = agg
      )
    }
    metadata(object) <- metadata
    object
  }
)

.ACAT_validate_input <- function(args) {
  # type checks
  check_wrapper(check_character, args, "group")
  check_wrapper(check_bool, args, "fixpval", length_equal = 1L)
  check_wrapper(check_character, args, "fixpval_method")
  check_wrapper(check_number_decimal, args, "fixpval_maxP", length_equal = 1L)
  check_wrapper(check_number_decimal, args, "fixpval_minP", length_equal = 1L)
  check_wrapper(check_bool, args, "warning", length_equal = 1L)
  if (!is.character(args[["aggregate"]]) && !is.list(args[["aggregate"]])) {
    stop(
      "`aggregate` should be either a list or a character vector.",
      call. = FALSE
    )
  }

  # number ranges
  check_number_between(
    args[["fixpval_maxP"]],
    arg = "fixpval_maxP",
    min = 0.0,
    max = 1.0
  )
  check_number_between(
    args[["fixpval_minP"]],
    arg = "fixpval_minP",
    min = 0.0,
    max = 1.0
  )

  # check if columns are available
  cols <- unique(c(unlist(args[["aggregate"]]), args[["group"]]))
  if (!all(cols %in% colnames(args[["object"]]))) {
    missing_cols <- setdiff(cols, colnames(args[["object"]]))
    stop(
      sprintf(
        paste0(
          "The following columns specified in `aggregate` or ",
          "`group` were not found: %s"
        ),
        paste(sQuote(missing_cols, "'"), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

.ACAT_update_metadata <- function(
  object,
  metadata,
  agg,
  group,
  fixpval,
  fixpval_method,
  fixpval_maxP,
  fixpval_minP
) {
  # determine aggregation level
  if ("ACAT" %in% names(metadata)) {
    acat_level <- max(vapply(metadata$ACAT, `[[`, integer(1L), "level")) + 1L
  } else {
    metadata$ACAT <- list()
    acat_level <- 1L
  }

  metadata[["ACAT"]][[paste0("ACAT", acat_level)]] <- list(
    level = acat_level,
    aggregate = agg,
    group = group,
    args = list(
      fixpval = fixpval,
      fixpval_method = fixpval_method,
      fixpval_maxP = fixpval_maxP,
      fixpval_minP = fixpval_minP
    ),
    summary = summary(object, asList = TRUE)
  )

  metadata
}

.ACAT_fix_pvals <- function(
  object,
  group,
  fixpval = TRUE,
  fixpval_method = "minmax",
  fixpval_maxP = 0.99,
  fixpval_minP = 1e-32,
  warning = TRUE
) {
  ## Fix P-values that are exactly 1 or 0.
  # see: FAQ in https://github.com/yaowuliu/ACAT for other solutions
  if (fixpval && fixpval_method %in% c("minmax", "manual")) {
    if (fixpval_method == "minmax") {
      P_max <- max(object[["P"]][object[["P"]] != 1], na.rm = TRUE)
      P_min <- min(object[["P"]][object[["P"]] != 0], na.rm = TRUE)
    } else if (fixpval_method == "manual") {
      P_max <- fixpval_maxP
      P_min <- fixpval_minP
    }
    nP1 <- sum(object[["P"]] == 1, na.rm = TRUE)
    nP0 <- sum(object[["P"]] == 0, na.rm = TRUE)
    if (nP1 > 0L) {
      if (warning) {
        warning(
          sprintf(
            "%s P-values are exactly 1, these are set to %s",
            sum(nP1, na.rm = TRUE),
            P_max
          ),
          call. = FALSE
        )
      }
      object[["P"]] <- ifelse(object[["P"]] == 1, P_max, object[["P"]])
    }

    if (nP0 > 0L) {
      if (warning) {
        warning(
          sprintf(
            "%s P-values are exactly 0, these are set to %s",
            sum(nP0, na.rm = TRUE),
            P_min
          ),
          call. = FALSE
        )
      }
      object[["P"]] <- ifelse(object[["P"]] == 0, P_min, object[["P"]])
    }
    dat <- as.data.frame(object)[, c(group, "P")]
  } else if (fixpval_method == "Liu") {
    dat <- as.data.frame(object)[, c(group, "P")] %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>%
      dplyr::mutate(
        P = ifelse(P == 1, fixP_Liu(P, d = dplyr::n()), P)
      ) %>%
      dplyr::ungroup()

    nP0 <- sum(dat[["P"]] == 0, na.rm = TRUE)
    if (nP0 > 0L) {
      if (warning) {
        warning(
          sprintf(
            "%s P-values are exactly 0, these are set to %s",
            sum(nP0, na.rm = TRUE),
            fixpval_minP
          ),
          call. = FALSE
        )
      }
      dat[["P"]] <- ifelse(dat[["P"]] == 0, fixpval_minP, dat[["P"]])
    }
  }

  dat
}

.ACAT_parse_results <- function(object, acat, group, agg) {
  # fields and their types to summarize
  col_character <- c(
    "unit",
    "cohort",
    "varSetName",
    "name",
    "pheno",
    "covar",
    "geneticModel",
    "MAFweight",
    "test"
  )
  col_integer <- c(
    "nvar",
    "caseN",
    "ctrlN"
  )
  col_numeric <- c(
    "caseCarriers",
    "ctrlCarriers",
    "meanCaseScore",
    "meanCtrlScore",
    "caseCallRate",
    "ctrlCallRate"
  )
  col_character <- col_character[!col_character %in% group]
  col_integer <- col_integer[!col_integer %in% group]
  col_numeric <- col_numeric[!col_numeric %in% group]

  # summarize
  object <- as.data.frame(object) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>%
    dplyr::summarize(
      dplyr::across(
        dplyr::all_of(col_character),
        .ACAT_summarise_character_col,
        .names = "{.col}"
      ),
      dplyr::across(
        dplyr::all_of(col_integer),
        .ACAT_summarise_integer_col,
        .names = "{.col}"
      ),
      dplyr::across(
        dplyr::all_of(col_numeric),
        .ACAT_summarise_numeric_col,
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
    dplyr::left_join(acat, by = group)

  # set aggregate columns to "ACAT"
  object[, agg] <- "ACAT"

  # convert to rvbResult
  object <- rvbResult(object)

  # return
  object
}

.ACAT_summarise_character_col <- function(col) {
  if (length(unique(col)) == 1L) col[1L] else NA_character_
}

.ACAT_summarise_integer_col <- function(col) {
  if (length(unique(col)) == 1L) col[1L] else NA_integer_
}

.ACAT_summarise_numeric_col <- function(col) {
  if (
    !all(is.na(col)) &&
      sum(!is.na(col)) > 1L &&
      var(col, na.rm = TRUE) == 0
  ) {
    col[1L]
  } else {
    NA_real_
  }
}

fixP_Liu <- function(x, d) {
  if (d > 1) (1 - 1 / d) else (x)
}

# singlevarResults ------------------------------------------------------------

## Class def. and constructor --------------------------------------------------

# Class definition

proto_singlevarResult <- function(n = 0L) {
  object <- S4Vectors::DataFrame(
    VAR_id = rep(NA_character_, n),
    cohort = Rle(rep(NA_character_, n)),
    varSetName = Rle(rep(NA_character_, n)),
    name = Rle(rep(NA_character_, n)),
    pheno = Rle(rep(NA_character_, n)),
    covar = Rle(rep(NA_character_, n)),
    geneticModel = Rle(rep(NA_character_, n)),
    test = Rle(rep(NA_character_, n)),
    caseMAC = rep(NA_real_, n),
    ctrlMAC = rep(NA_real_, n),
    caseMAF = rep(NA_real_, n),
    ctrlMAF = rep(NA_real_, n),
    caseN = rep(NA_real_, n),
    ctrlN = rep(NA_real_, n),
    caseCallRate = rep(NA_real_, n),
    ctrlCallRate = rep(NA_real_, n),
    effectAllele = rep(NA_character_, n),
    otherAllele = rep(NA_character_, n),
    effect = rep(NA_real_, n),
    effectSE = rep(NA_real_, n),
    effectCIlower = rep(NA_real_, n),
    effectCIupper = rep(NA_real_, n),
    OR = rep(NA_real_, n),
    P = rep(NA_real_, n)
  )
  object
}


#' singlevarResult
#' @rdname rvatResult
#' @usage NULL
#' @export
singlevarResult <- function(object) {
  # return an empty singlevarResult if input is empty
  if (missing(object) || is.null(object)) {
    return(new("singlevarResult", proto_singlevarResult()))
  }

  # object should be a filepath or a data.frame/DFrame
  if (
    !is.character(object) && !is.data.frame(object) && !is(object, "DFrame")
  ) {
    stop(
      "`object` should be either a data.frame/DataFrame or a filepath.",
      call. = FALSE
    )
  }

  # assume filepath if object if of type character
  if (is.character(object)) {
    return(readResults(path = object, type = "singlevarResult", sep = "\t"))
  }

  # format and fill in missing fields if object is a data.frame/DFrame
  if (is.data.frame(object) || is(object, "DFrame")) {
    object <- .rvatResult_format(
      object,
      id_col = "VAR_id",
      columns_info = columns_singlevarResults,
      proto_func = proto_singlevarResult,
      rle_cols = columns_singlevarResults_rle,
      num_cols = columns_singlevarResults_numeric
    )

    # return singlevarResult object
    new("singlevarResult", object)
  }
}

# validity
setValidity2("singlevarResult", function(object) {
  .rvatResult_setvalidity(object, columns_info = columns_singlevarResults)
})

## Methods ---------------------------------------------------------------------

### summary --------------------------------------------------------------------

#' summarize singlevarResult
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("summary", "singlevarResult", function(object, asList = FALSE, ...) {
  info <- list()
  info[["Nvars"]] <- length(unique(object[["VAR_id"]]))

  info[["ctrlNmin"]] <- if ((all(is.na(object[["ctrlN"]])))) {
    NA_real_
  } else if (nrow(object) > 0L) {
    min(object[["ctrlN"]], na.rm = TRUE)
  } else {
    0L
  }
  info[["ctrlNmax"]] <- if ((all(is.na(object[["ctrlN"]])))) {
    NA_real_
  } else if (nrow(object) > 0L) {
    max(object[["ctrlN"]], na.rm = TRUE)
  } else {
    0L
  }
  info[["caseNmin"]] <- if ((all(is.na(object[["caseN"]])))) {
    NA_real_
  } else if (nrow(object) > 0L) {
    min(object[["caseN"]], na.rm = TRUE)
  } else {
    0L
  }
  info[["caseNmax"]] <- if ((all(is.na(object[["caseN"]])))) {
    NA_real_
  } else if (nrow(object) > 0L) {
    max(object[["caseN"]], na.rm = TRUE)
  } else {
    0L
  }
  info[["names"]] <- unique(object[["name"]])
  info[["varsets"]] <- unique(object[["varSetName"]])
  info[["covar"]] <- unique(object[["covar"]])
  info[["tests"]] <- unique(object[["test"]])
  info[["geneticModels"]] <- unique(object[["geneticModel"]])
  info[["pheno"]] <- unique(object[["pheno"]])
  info[["cohorts"]] <- unique(object[["cohort"]])

  # return as is if asList=TRUE
  if (asList) {
    return(info)
  }

  cat(sprintf("%s object\n-------------------\n", class(object)[1]))
  cat(sprintf(
    paste(
      "n vars = %s",
      "N cases = %s",
      "N controls = %s",
      "names = %s",
      "varSets = %s",
      "covars = %s",
      "tests = %s",
      "geneticModels = %s",
      "pheno = %s",
      "cohorts = %s",
      sep = "\n"
    ),
    info$Nvars,
    info$caseNmax,
    info$ctrlNmax,
    .summary_truncate(info$names),
    .summary_truncate(info$varsets),
    .summary_truncate(info$covar),
    .summary_truncate(info$tests),
    .summary_truncate(info$geneticModels),
    .summary_truncate(info$pheno),
    .summary_truncate(info$cohorts)
  ))
  cat("\n")
})
