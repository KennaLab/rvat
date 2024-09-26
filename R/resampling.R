#==============================================================================
# reSampling objects and methods
#==============================================================================

#' resamplingFile
#' @rdname resamplingFile
#' @usage NULL
#' @export
resamplingFile=function(path)
{
  con=gzfile(path,"r")
  methodResampling <- scan(con, nlines = 1, what = "character", quiet = TRUE)
  nSamples <- scan(con, nlines = 1, what = "character", quiet = TRUE)
  nResampling <- scan(con, nlines = 1, what = "character", quiet = TRUE)
  close(con)
  new("resamplingFile", path=path, methodResampling=methodResampling, nSamples=as.numeric(nSamples), nResampling=as.numeric(nResampling))
}

setMethod("show", signature = "resamplingFile",
          definition = function(object) {
            cat(sprintf("resamplingFile object\nPath:%s\nResampling method: %s\nN samples: %s\nN resamplings: %s", 
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
#' @param nSamples Number of samples
#' @param nResampling Number of resamplings
#' @param memlimit Chunk sizes (when writing to output)
#' @param methodResampling Resampling method, currently 'permutation' is implemented.
#' @param output File path (.gz extension) to write output to. If not specified, a matrix with resamplings is returned.
#' @export
buildResamplingFile <- function(nSamples,
                                nResampling, 
                                memlimit=1000,
                                methodResampling = "permutation", 
                                output = NULL)  {
  
  if(methodResampling == "permutation") {
    if(!is.null(output)) {
      output <- gzcon(file(output,open='wb'))
      write(methodResampling, 
            file = output, append = FALSE) 
      write(nSamples, 
            file = output, append = TRUE)
      write(nResampling, 
            file = output, append = TRUE)
      
      chunks <-rep(memlimit, times = nResampling %/% memlimit)
      if(sum(chunks)-nResampling != 0) chunks <- c(chunks, nResampling-sum(chunks))
      for(i in 1:length(chunks))  {
        perms <- matrix(rep(1:nSamples, chunks[i]), ncol=chunks[i], byrow=FALSE)
        perms <- apply(perms, 2, sample)
        write.table( as.data.frame(t(perms)), col.names=FALSE, row.names=FALSE,
                     sep="\t",file = output, append = TRUE) 
      }
      close(output)
    }  else {
        perms <- matrix(rep(1:nSamples, nResampling), ncol=nResampling, byrow=FALSE)
        perms <- apply(perms, 2, sample)
        colnames(perms) <- paste0("perm", 1:nResampling)
        perms
    }
  }
}
