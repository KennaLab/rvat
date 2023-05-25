#' @include allClasses.R  
#' @include allGenerics.R
#' @include allInternalData.R


#' @rdname aggregateFile
#' @usage NULL
#' @export
aggregateFile=function(path)
{
  con=gzfile(path,"r")
  samples <- unlist(strsplit(scan(con, nlines = 1, what = "character", quiet = TRUE), split = ","))
  units <- unlist(strsplit(scan(con, nlines = 1, what = "character", quiet = TRUE), split = ","))
  close(con)
  new("aggregateFile", path=path, units=units, samples=samples)
}

setMethod("show", signature = "aggregateFile",
          definition = function(object) {
            message("rvat aggregateFile object")
            message(sprintf("Path:%s",object@path))
            message(sprintf("Samples:%s",length(object@samples)))
            message(sprintf("Units:%s",length(object@units)))
          })

setMethod("listUnits", signature = "aggregateFile",
          definition = function(object) {
            object@units
          })

setMethod("listSamples", signature = "aggregateFile",
          definition = function(object) {
            object@samples
          })


setMethod("getUnit", signature = "aggregateFile",
          definition = function(object, unit) {
            if(!all(unit %in% listUnits(object))) {
              stop(sprintf("The following units are not present in the aggregateFile: %s",
                           paste(unit[!unit %in% listUnits(object)], collapse=",")
              ))
            }
            indices <- sort(which(listUnits(object) %in% unit))
            unit <- listUnits(object)[indices]
            if(length(indices) > 1) indices[2:length(indices)] <- (dplyr::lead(indices)-indices)[1:(length(indices)-1)]
            indices[1] <- indices[1] + 2
            
            con <- gzfile(object@path,"r")
            dat <- lapply(indices,
                          FUN = function(i) {
                            scan(con, skip = i-1, nlines = 1, what = "character", quiet = TRUE)
                          })
            dat <- lapply(dat, FUN = function(x)  as.numeric(unlist(strsplit(memDecompress(as.raw(as.hexmode(unlist(strsplit(x, split = ",")))), type = "gzip", asChar = TRUE), split = "\n"))) )
            dat <- do.call(rbind, dat)
            rownames(dat) <- unit
            colnames(dat) <- listSamples(object)
            close(con)
            dat
          })


#' @rdname aggregateFileList
#' @usage NULL
#' @export
aggregateFileList <- function(filelist, checkDups = TRUE) {
  lst <- lapply(
    filelist, 
    FUN = aggregateFile
  )
  
  paths <- unlist(lapply(lst, FUN = function(x) x@path))
  
  # Check whether duplicate units are included
  units <- unlist(lapply(lst, listUnits))
  if(sum(duplicated(units)) > 0) {
    if (checkDups) {
      stop(sprintf("The following units are duplicated across aggregateFiles: %s",
                   paste(units[duplicated(units)], collapse = ",")))
    } else {
      units <- paste0("X", 1:length(units))
    }
  }
  
  # Samples -> check whether sample vectors are identical
  samples <- lapply(lst, listSamples)
  test <- unlist(lapply(samples, FUN = function(x) {identical(samples[[1]], x)}))
  if(!all(test)) {
    stop("Aggregate files should include identical samples (in the same order)")
  } 
  samples <- samples[[1]]
  
  new("aggregateFileList", 
      paths=paths, 
      units=units, 
      samples=samples)
}

setMethod("show", signature = "aggregateFileList",
          definition = function(object) {
            message("rvat aggregateFileList object")
            message(sprintf("aggregateFiles:%s",length(object@paths)))
            message(sprintf("Samples:%s",length(object@samples)))
            message(sprintf("Units:%s",length(object@units)))
          })

setMethod("listUnits", signature = "aggregateFileList",
          definition = function(object) {
            object@units
          })

#' @export
setMethod("listSamples", signature = "aggregateFileList",
          definition = function(object) {
            object@samples
          })

setMethod("length", signature = "aggregateFileList",
          definition = function(x) {
            length(x@paths)
          })

# mergeAggregateFiles -------------------------------------------------------------------------

#' @rdname mergeAggregateFiles
#' @usage NULL
#' @export
setMethod("mergeAggregateFiles", 
          signature = signature(object="aggregateFileList"),
          definition=function(
            object,
            collapse = TRUE,
            output = NULL,
            verbose = TRUE
          )
          {
            if(collapse) {
              agg <- vector(mode = "numeric", length = length(listSamples(object)))
              
              for(i in 1:length(object)) {
                if(verbose) message(sprintf("%s/%s", i, length(object)))
                dat <- aggregateFile(object@paths[i])
                dat <- getUnit(dat, unit = listUnits(dat))
                dat <- colSums(dat) 
                agg <- agg + dat
                }
              agg <- data.frame(
                IID = listSamples(object),
                aggregate = agg,
                stringsAsFactors = FALSE
              )
              if(is.null(output)) {
                return(agg)
              } else {
                write.table(
                  agg,
                  file = gzfile(output),
                  quote = FALSE,
                  sep = "\t",
                  row.names = FALSE
                )
              }
               } else {
                 if(!is.null(output)) {
                   output <- gzcon(file(output,open='wb'))
                   
                   write(paste(listSamples(object), collapse = ","), 
                         file = output, append = FALSE) 
                   write(paste(listUnits(object), collapse=","), 
                         file = output, append = TRUE)
                   dat <- lapply(
                     1:length(object),
                     FUN = function(i, object, verbose) {
                       if(verbose) message(sprintf("%s/%s", i, length(object)))
                       dat <- read.table(object@paths[i], skip = 2, header = FALSE)
                       write.table(
                         dat,
                         file = output,
                         quote = FALSE,
                         sep = "\t",
                         row.names = FALSE,
                         col.names = FALSE
                       )
                     },
                     object = object,
                     verbose = verbose
                   )
                   close(output)
                 } else {
                   dat <- lapply(
                     1:length(object),
                     FUN = function(i, object, verbose) {
                       if(verbose) message(sprintf("%s/%s", i, length(object)))
                       dat <- aggregateFile(object@paths[i])
                       dat <- getUnit(dat, unit = listUnits(dat))
                       dat
                     },
                     object = object,
                     verbose = verbose
                   )
                   
                   dat <- do.call(rbind, dat)
                   return(dat)
                 }
              }
            }
          )
