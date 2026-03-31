#' @include allClasses.R
#' @include allGenerics.R
#' @include allInternalData.R

#' @rdname aggdbList
#'
#' @param filelist Character vector of file paths pointing to [`aggdb`] files.
#' @param checkDups Logical. If `TRUE` (the default), an error is raised when
#'   unit names are duplicated across the provided aggdbs. If `FALSE`, duplicated
#'   unit names are replaced with unique identifiers.
#'
#' @return An [`aggdbList-class`] object.
#' @export
aggdbList <- function(filelist, checkDups = TRUE) {
  # input validation
  check_character(filelist)
  check_length(filelist, min = 1L)
  check_bool(checkDups)

  # extract info from each aggdb one at a time to avoid opening too many
  # connections simultaneously
  units_list <- vector("list", length(filelist))
  samples_list <- vector("list", length(filelist))
  metadata_list <- vector("list", length(filelist))
  params_list <- vector("list", length(filelist))
  con <- NULL
  on.exit(
    if (!is.null(con) && DBI::dbIsValid(con)) try(close(con), silent = TRUE)
  )
  for (i in seq_along(filelist)) {
    con <- aggdb(filelist[[i]])
    units_list[[i]] <- listUnits(con)
    samples_list[[i]] <- listSamples(con)
    metadata_list[[i]] <- metadata(con)
    params_list[[i]] <- listParams(con)
    close(con)
  }
  con <- NULL

  # check whether duplicate units are included
  units_all <- unlist(units_list)
  if (anyDuplicated(units_all) != 0L) {
    if (checkDups) {
      stop(
        sprintf(
          "The following units are duplicated across databases: %s",
          paste(units_all[duplicated(units_all)], collapse = ",")
        ),
        call. = FALSE
      )
    } else {
      units_all <- paste0("X", seq_along(units_all))
    }
  }

  # check whether sample IDs are identical across dbs
  check_identical_samples <- unlist(lapply(samples_list, FUN = function(x) {
    identical(samples_list[[1L]], x)
  }))
  if (!all(check_identical_samples)) {
    stop(
      "aggdbs should include identical samples (in the same order).",
      call. = FALSE
    )
  }
  samples <- samples_list[[1L]]

  # check metadata
  metadata <- metadata_list
  rvatversion <- unique(unlist(lapply(
    metadata,
    function(x) x[["rvatVersion"]]
  )))
  gdbid <- unique(unlist(lapply(
    metadata,
    function(x) x[["gdbId"]]
  )))
  genomebuild <- unique(unlist(lapply(
    metadata,
    function(x) x[["genomeBuild"]]
  )))
  params <- params_list
  check_identical_params <- unlist(lapply(params, FUN = function(x) {
    identical(params[[1L]], x)
  }))
  if (!all(check_identical_params)) {
    stop(
      "Non-identical parameters were used to generate the input aggdbs.",
      call. = FALSE
    )
  }
  if (length(rvatversion) > 1L) {
    stop(
      "Aggdbs were generated using different rvat versions.",
      call. = FALSE
    )
  }
  if (length(gdbid) > 1L) {
    stop("Aggdbs were generated from different gdbs", call. = FALSE)
  }
  if (length(genomebuild) > 1L) {
    stop(
      sprintf(
        "Aggdbs were generated using different genome builds: %s.",
        paste(genomebuild, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  metadata <- list(
    rvatVersion = rvatversion,
    gdbId = gdbid,
    genomeBuild = if (length(genomebuild) > 0L) {
      genomebuild[1L]
    } else {
      NA_character_
    },
    creationDate = as.character(round(Sys.time(), units = "secs"))
  )

  # return aggdblist
  new(
    "aggdbList",
    paths = filelist,
    units = units_all,
    samples = samples,
    params = params[[1L]],
    metadata = metadata
  )
}

# aggdbList methods -------------------------------------------------------

#' @rdname aggdbList
#' @usage NULL
#' @export
setMethod("show", signature = "aggdbList", definition = function(object) {
  cat(sprintf(
    "aggdbList object\naggdbs: %s\nSamples: %s\nUnits: %s\n",
    length(object@paths),
    length(object@samples),
    length(object@units)
  ))
})

#' @rdname aggdbList
#' @usage NULL
#' @export
setMethod("listUnits", signature = "aggdbList", definition = function(object) {
  object@units
})

#' @rdname aggdbList
#' @usage NULL
#' @export
setMethod(
  "listSamples",
  signature = "aggdbList",
  definition = function(object) {
    object@samples
  }
)

#' @rdname aggdbList
#' @usage NULL
#' @export
setMethod("length", signature = "aggdbList", definition = function(x) {
  length(x@paths)
})

#' @rdname aggdbList
#' @usage NULL
#' @export
setMethod("metadata", signature = "aggdbList", definition = function(x) {
  x@metadata
})

#' @rdname aggdbList
#' @usage NULL
#' @export
setMethod("listParams", signature = "aggdbList", definition = function(object) {
  object@params
})
