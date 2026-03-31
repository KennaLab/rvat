#' resamplingFile
#' @rdname resamplingFile
#' @usage NULL
#' @export
resamplingFile <- function(path) {
  con = gzfile(path, "r")
  on.exit(close(con), add = TRUE)
  methodResampling <- scan(con, nlines = 1, what = "character", quiet = TRUE)
  nSamples <- scan(con, nlines = 1, what = "character", quiet = TRUE)
  nResampling <- scan(con, nlines = 1, what = "character", quiet = TRUE)
  new(
    "resamplingFile",
    path = path,
    methodResampling = methodResampling,
    nSamples = as.numeric(nSamples),
    nResampling = as.numeric(nResampling)
  )
}

#' @rdname resamplingFile
#' @usage NULL
#' @export
setMethod("show", signature = "resamplingFile", definition = function(object) {
  cat(sprintf(
    "resamplingFile object\nPath:%s\nResampling method: %s\nN samples: %s\nN resamplings: %s",
    object@path,
    object@methodResampling,
    object@nSamples,
    object@nResampling
  ))
})

#' Build a resampling matrix
#'
#' Currently the 'permutation' approach is implemented. If `output` is specified,
#' the resampling matrix will be written to disk, and can be connected to with [`resamplingFile`],
#' otherwise a matrix is returned. Can be used in combination with [`assocTest`] to perform
#' resampling tests.
#'
#' @param nSamples Number of samples.
#' @param nResampling Number of resamplings. Defaults to 1000.
#' @param memlimit Maximum number of resamplings to generate at a time (chunk size).
#' Relevant when writing to `output` to manage memory usage. Defaults to 1000.
#' @param methodResampling Resampling method, currently 'permutation' is implemented.
#' @param output File path (.gz extension) to write output to.
#' If not specified, a matrix with resamplings is returned.
#' @return If `output` is `NULL`, returns a matrix of dimensions `nSamples` x `nResampling`.
#' If `output` is specified, the function writes to the file and returns `NULL`.
#' @example inst/examples/example-resamplingFile.R
#' @export
buildResamplingFile <- function(
  nSamples,
  nResampling = 1000,
  memlimit = 1000,
  methodResampling = "permutation",
  output = NULL
) {
  if (methodResampling == "permutation") {
    if (!is.null(output)) {
      output <- gzcon(file(output, open = "wb"))
      on.exit(close(output), add = TRUE)
      write(methodResampling, file = output, append = FALSE)
      write(nSamples, file = output, append = TRUE)
      write(nResampling, file = output, append = TRUE)

      chunks <- rep(memlimit, times = nResampling %/% memlimit)
      if (sum(chunks) - nResampling != 0) {
        chunks <- c(chunks, nResampling - sum(chunks))
      }
      for (i in seq_along(chunks)) {
        perms <- matrix(
          rep(1:nSamples, chunks[i]),
          ncol = chunks[i],
          byrow = FALSE
        )
        perms <- apply(perms, 2, sample)
        write.table(
          as.data.frame(t(perms)),
          col.names = FALSE,
          row.names = FALSE,
          sep = "\t",
          file = output,
          append = TRUE
        )
      }
    } else {
      perms <- matrix(
        rep(1:nSamples, nResampling),
        ncol = nResampling,
        byrow = FALSE
      )
      perms <- apply(perms, 2, sample)
      colnames(perms) <- paste0("perm", 1:nResampling)
      perms
    }
  } else {
    stop(
      sprintf("Unknown methodResampling: %s", methodResampling),
      call. = FALSE
    )
  }
}
