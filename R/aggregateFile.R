#' @include allClasses.R  
#' @include allGenerics.R
#' @include allInternalData.R


#' @rdname aggregateFile
#' @usage NULL
#' @export
aggregateFile=function(path)
{
  
  # read in metadata
  metadata <- .parse_rvat_header(path, 
                                 expected_metadata = metadata_aggregate,
                                 expected_filetype = "aggregateFile",
                                 n = length(metadata_aggregate) + 1 # file description + metadata
  )
  
  # get sample and unit IDs
  header <- readLines(path, n = length(metadata_aggregate) + 1) # filetype + metadata
  skip <- sum(startsWith(header, "#"))
  con <- gzfile(path,"r")
  if(skip > 0) skip <- readLines(con, n = skip)
  samples <- unlist(strsplit(scan(con, nlines = 1, what = "character", quiet = TRUE), split = ","))
  units <- unlist(strsplit(scan(con, nlines = 1, what = "character", quiet = TRUE), split = ","))
  close(con)
  new("aggregateFile", path=path, units=units, samples=samples, metadata=metadata)
}

#' @rdname aggregateFile
#' @usage NULL
#' @export
setMethod("show", signature = "aggregateFile",
          definition = function(object) {
            cat(sprintf("aggregateFile object\nPath: %s\nSamples: %s\nUnits: %s\n",
                        object@path,
                        length(object@samples),
                        length(object@units)))
          })

#' @rdname aggregateFile
#' @usage NULL
#' @export
setMethod("listUnits", signature = "aggregateFile",
          definition = function(object) {
            object@units
          })

#' @rdname aggregateFile
#' @usage NULL
#' @export
setMethod("listSamples", signature = "aggregateFile",
          definition = function(object) {
            object@samples
          })

#' @rdname aggregateFile
#' @usage NULL
#' @export
setMethod("metadata", signature="aggregateFile",
          definition=function(x)
          {
            x@metadata
          })

#' @rdname aggregateFile
#' @usage NULL
#' @export
setMethod("getGdbId", signature="aggregateFile",
          definition=function(object)
          {
            metadata(object)$gdbId
          })

#' @rdname aggregateFile
#' @usage NULL
#' @export
setMethod("getRvatVersion", signature="aggregateFile",
          definition=function(object)
          {
            metadata(object)$rvatVersion
          })

#' @rdname aggregateFile
#' @usage NULL
#' @export
setMethod("getUnit", signature = "aggregateFile",
          definition = function(object, unit) {
            
            # check if unit is present in aggregatefile
            if(!all(unit %in% listUnits(object))) {
              stop(sprintf("The following units are not present in the aggregateFile: %s",
                           paste(unit[!unit %in% listUnits(object)], collapse=",")
              ))
            }
            
            # parse metadata
            header <- readLines(object@path, n = length(metadata_aggregate) + 1) # metadata + filetype
            skip <- sum(startsWith(header, "#"))
            
            # get indices
            indices <- sort(which(listUnits(object) %in% unit))
            unit <- listUnits(object)[indices]
            if(length(indices) > 1) indices[2:length(indices)] <- (dplyr::lead(indices) - indices)[1:(length(indices)-1)]
            indices[1] <- indices[1] + skip + 2
            
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
  
  # Check and add metadata
  metadata <- lapply(lst, metadata)
  rvatversion <- unique(unlist(lapply(metadata, function(x) x[["rvatVersion"]])))
  gdbid <- unique(unlist(lapply(metadata, function(x) x[["gdbId"]])))
  genomebuild <- unique(unlist(lapply(metadata, function(x) x[["genomeBuild"]])))
  if (length(rvatversion) > 1) {stop("AggregateFiles were generated using different rvat versions.")}
  if (length(gdbid) > 1) {stop("AggregateFiles were generated from different gdbs")}
  
  metadata <- list(
    rvatVersion = rvatversion,
    gdbId = gdbid,
    genomeBuild = genomebuild[1],
    creationDate = as.character(round(Sys.time(), units = "secs"))
  )
  
  new("aggregateFileList", 
      paths=paths, 
      units=units, 
      samples=samples,
      metadata=metadata
      )
}

#' @rdname aggregateFile
#' @usage NULL
#' @export
setMethod("show", signature = "aggregateFile",
          definition = function(object) {
            cat(sprintf("aggregateFile object\nPath: %s\nSamples: %s\nUnits: %s\n",
                        object@path,
                        length(object@samples),
                        length(object@units)))
          })


#' @rdname aggregateFileList
#' @usage NULL
#' @export
setMethod("show", signature = "aggregateFileList",
          definition = function(object) {
            cat(sprintf("aggregateFileList object\naggregateFiles: %s\nSamples: %s\nUnits: %s\n",
                        length(object@paths),
                        length(object@samples),
                        length(object@units)))
            
          })

#' @rdname aggregateFileList
#' @usage NULL
#' @export
setMethod("listUnits", signature = "aggregateFileList",
          definition = function(object) {
            object@units
          })

#' @rdname aggregateFileList
#' @usage NULL
#' @export
setMethod("listSamples", signature = "aggregateFileList",
          definition = function(object) {
            object@samples
          })

#' @rdname aggregateFileList
#' @usage NULL
#' @export
setMethod("length", signature = "aggregateFileList",
          definition = function(x) {
            length(x@paths)
          })

#' @rdname aggregateFileList
#' @usage NULL
#' @export
setMethod("metadata", signature="aggregateFileList",
          definition=function(x)
          {
            x@metadata
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
                   
                   # write metadata
                   metadata <- metadata(object)
                   metadata$creationDate <- as.character(round(Sys.time(), units = "secs"))
                   .write_rvat_header(filetype = "aggregateFile", 
                                      metadata = metadata(object), 
                                      con = output)
                   
                   # write sample and unit IDs
                   write(paste(listSamples(object), collapse = ","), 
                         file = output, append = FALSE) 
                   write(paste(listUnits(object), collapse=","), 
                         file = output, append = TRUE)
                   dat <- lapply(
                     1:length(object),
                     FUN = function(i, object, verbose) {
                       if(verbose) message(sprintf("%s/%s", i, length(object)))
                       # skip header
                       header <- readLines(object@paths[i], n = length(metadata_aggregate) + 1) # metadata + filetype
                       skip <- sum(startsWith(header, "#"))
                       dat <- read.table(object@paths[i], skip = 2 + skip, header = FALSE)
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
