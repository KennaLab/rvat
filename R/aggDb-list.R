#' @include allClasses.R
#' @include allGenerics.R
#' @include allInternalData.R

#' @rdname aggdbList
#' @export
aggdbList <- function(filelist, checkDups = TRUE) {
  # input validation
  check_character(filelist)
  check_length(filelist, min = 1L)
  check_bool(checkDups)

  # connect to aggdbs
  aggdb_list <- lapply(
    filelist,
    FUN = aggdb
  )
  ## close all gdbs on exit
  on.exit(lapply(
    aggdb_list,
    function(con) if (DBI::dbIsValid(con)) try(close(con), silent = TRUE)
  ))

  # check whether duplicate units are included
  units_all <- unlist(lapply(aggdb_list, listUnits))
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
  samples <- lapply(aggdb_list, listSamples)
  check_identical_samples <- unlist(lapply(samples, FUN = function(x) {
    identical(samples[[1L]], x)
  }))
  if (!all(check_identical_samples)) {
    stop(
      "aggdbs should include identical samples (in the same order).",
      call. = FALSE
    )
  }
  samples <- samples[[1L]]

  # check metadata
  metadata <- lapply(aggdb_list, metadata)
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
  params <- lapply(
    aggdb_list,
    listParams
  )
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
    stop(sprintf("Aggdbs were generated using different genome builds: %s.", paste(genomebuild, collapse = ", ")), call. = FALSE)
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

setMethod("listParams", signature = "aggdbList", definition = function(object) {
  object@params
})
