#' @include rvatResult.R
#' @include gsaResult.R
#==============================================================================
# geneSet objects
#==============================================================================

# geneSet ----------------------------------------------------------------------

#' @rdname geneSet
#' @usage NULL
#' @export
geneSet=function(geneSetName,units,w=NULL,metadata="")
{
  if(is.null(w) || is.na(w)) w <- paste(rep(1, length(unlist(strsplit(units, split = ",")))), collapse=",")
  if(is.null(metadata) || is.na(metadata)) metadata <- ""
  new("geneSet",geneSetName = geneSetName, units = units, w=w, metadata=metadata)
}

#' @rdname geneSet
#' @usage NULL
#' @export
setMethod("metadata", signature = "geneSet",
          definition = function(x) {
            x@metadata
          })

setMethod("length", signature = "geneSet",
          definition = function(x) {
            length(listUnits(x))
          })

setMethod("as.data.frame", signature = "geneSet",
          definition = function(x) {
            data.frame(geneSetName = x@geneSetName, 
                       units = x@units,
                       w = x@w,
                       metadata = x@metadata,
                       stringsAsFactors = FALSE
            )
          })

setMethod("listUnits", signature = "geneSet",
          definition = function(object) {
            unlist(strsplit(object@units,split=","))
          })


setMethod("listWeights", signature = "geneSet",
          definition = function(object) {
            x <- unlist(strsplit(object@w,split=","))
            names(x) <- listUnits(object)
            x
          })


# geneSetList ------------------------------------------------------------------

#' geneSetList
#' @rdname geneSetList
#' @usage NULL
#' @export
geneSetList <- function(geneSets) {
  geneSetNames <- unlist(lapply(geneSets, function(x) {x@geneSetName}))
  new("geneSetList", geneSets = geneSets, geneSetNames = geneSetNames)
}


setMethod("names", signature = "geneSetList",
          definition = function(x) {
            x@geneSetNames
          })

setMethod("length", signature = "geneSetList",
          definition = function(x) {
            length(x@geneSets)
          })

setMethod("lengths", signature = "geneSetList",
          definition = function(x) {
            lengths(as.list(x))
          })

setMethod("show", signature = "geneSetList",
          definition = function(object) {
            message("geneSetList")
            message(sprintf("Contains %s sets",length(object)))
          })

setMethod("metadata", signature = "geneSetList",
          definition = function(x) {
            unlist(lapply(x@geneSets, FUN = metadata))
          })

setMethod("listGeneSets", signature = "geneSetList",
          definition = function(object) {
            names(object)
          })

setMethod("listUnits", signature = "geneSetList",
          definition = function(object) {
            unique(unlist(as.list(object)))
          })

setMethod("[[", c("geneSetList", "ANY", "missing"),
          function(x, i, j, ...)
          {
            x@geneSets[[i, ...]]
          })

setMethod("[", c("geneSetList", "ANY", "ANY"),
          function(x, i, j, ..., drop=TRUE)
          {
            if (1L != length(drop) || (!missing(drop) && drop))
              warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
            
            if (missing(i) && missing(j))
              return(x)
            
            initialize(x, geneSets = x@geneSets[i,...], geneSetNames = x@geneSetNames[i,...])
          })


setMethod("sort", c("geneSetList"),
          function(x)
          {
            sorted <- sort(names(x))
            new("geneSetList", geneSets = x@geneSets[match(sorted, names(x))], geneSetNames = sorted)
          })


setMethod("as.data.frame", signature = "geneSetList",
          definition = function(x) {
            do.call(rbind, lapply(x@geneSets, FUN = as.data.frame))
          })


setMethod("as.list", signature = "geneSetList",
          definition = function(x) {
            names <- names(x)
            x <- strsplit(unlist(lapply(1:length(x), function(i) {x[[i]]@units})), split=",")
            names(x) <- names
            x
          })

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("mapToMatrix", signature = c("geneSetList", "rvbResult"),
          definition = function(object, results, ID = "unit", sparse = TRUE) {
            dat <- lapply(as.list(object),
                          FUN = function(gene_set, IDs) {
                            x <- IDs %in% gene_set
                          },
                          IDs =  unique(as.character(results[[ID]])))
            dat <- do.call(cbind, dat)
            rownames(dat) <- unique(as.character(results[[ID]]))
            colnames(dat) <- names(object)
            if(sparse) return(as(dat, "sparseMatrix")) else return(dat)
          })


setMethod("getGeneSet", signature = "geneSetList",
          definition = function(object, geneSet = NULL, unit = NULL) {
            if(is.null(geneSet) && is.null(unit)) message("At least one of `geneSet` or `unit` should be specified.")
            
            if(!is.null(unit)) {
              object <- object[unlist(lapply(object@geneSets, FUN = function(x) {
                any(unit %in% listUnits(x))
              }))]
            }
            if(!is.null(geneSet)) {
              object <- object[names(object) %in% geneSet]
            }
            return(object)
          })

setMethod("dropUnits", signature = "geneSetList",
          definition = function(object, unit = NULL) {
              object@geneSets <- lapply(object@geneSets,
                               function(x) {
                                 x@units <- paste(listUnits(x)[!listUnits(x) %in% unit], collapse = ",")
                                 x@w <-  paste(listWeights(x)[!listUnits(x) %in% unit], collapse = ",")
                                 return(x)
                               }
                              )
              object <- object[lengths(object) > 0]
              return(object)
          })

setMethod("remapIDs", signature = "geneSetList",
          definition = function(object, dict, targets = NULL, duplicate_ids = c("keep_all", "keep_first")) {
            
            # Check validity of dictionary
            if(ncol(dict) != 2) stop("`dict` should be a data.frame with two columns")
            colnames(dict) <- c("original_id", "new_id")
            dict$original_id <- as.character(dict$original_id)
            dict$new_id <- as.character(dict$new_id)
            
            # Number of IDs in geneSetList that are present in dictionary 
            original_ids <- listUnits(object)
            message(sprintf("%s/%s IDs in the geneSetList are present in the linker file.", 
                            sum(original_ids %in% dict$original_id), length(original_ids)
                            ))
            dict <- dict[dict$original_id %in% original_ids,,drop=FALSE]
            
            # Handle IDs that map to multiple IDs
            if(!is.null(targets)) {
              dict$target <- ifelse(dict$new_id %in% targets, TRUE, FALSE)
            } else {
             dict$target <- TRUE
            }
            duplicates <- dict %>% dplyr::count(original_id) %>% dplyr::filter(n > 1) %>% dplyr::ungroup()
            
            # handle duplicates
            ## if targets are specified, only need to handle duplicates which
            ## are both in the target list.
            if(!is.null(targets)) {
              check <- dict %>% 
                dplyr::filter(original_id %in% duplicates$original_id, target) %>%
                dplyr::count(original_id) 
              check_1 <- check %>% dplyr::filter(n == 1)
              duplicates1 <- dict %>%
                dplyr::filter(original_id %in% check_1$original_id, target)
              check_multiple <- check %>% dplyr::filter(n > 1)
            } else {
              check_multiple <- duplicates
              duplicates1 <- data.frame(original_id = character(), new_id = character(), target = logical(), stringsAsFactors = FALSE)
            }
      
            # if duplicate_ids == "keep_all", simply keep all duplicates
            # if duplicate_ids == "keep_first", keep first
            if(duplicate_ids == "keep_all") {
              duplicates2 <- dict %>% 
                dplyr::filter(original_id %in% check_multiple$original_id, target)
            } else if(duplicate_ids == "keep_first") {
              duplicates2 <- dict %>% 
                dplyr::filter(original_id %in% check_multiple$original_id, target) %>% 
                dplyr::group_by(original_id) %>%
                dplyr::mutate(keep = new_id[1]) %>%
                dplyr::ungroup() %>%
                dplyr::filter(new_id == keep) %>%
                dplyr::select(-keep)
            }
            duplicates_keep <- dplyr::bind_rows(duplicates1, duplicates2)
            
            # subset dictionary based on which duplicates to keep
            dict <- dplyr::bind_rows(
              dict %>% dplyr::filter(!original_id %in% duplicates$original_id),
              duplicates_keep
            )
            
            lengths_premapping <- lengths(object)
            
            # remap the IDs in each geneset based on the dictionary
            remapped <- lapply(
              object@geneSets,
              FUN = function(x, dct) {
                mapping <- data.frame(original_id = listUnits(x), w = listWeights(x), stringsAsFactors = FALSE) %>%
                          dplyr::left_join(dct, by = "original_id") %>%
                          dplyr::filter(!is.na(new_id))

                initialize(x,
                           units = paste(mapping$new_id, collapse=","),
                           w = paste(mapping$w, collapse=","))
                },
              dct = dict
            )
            object@geneSets <- remapped
            lengths_postmapping <- lengths(object)
            
            # return remapped object
            object
          })

#' @rdname geneSetList
#' @usage NULL
#' @export
setMethod("write", "geneSetList",
          function(x, file = "data", append = FALSE)
          {
            out <- gzfile(file,"w")
            write.table(as.data.frame(x), out, sep="|", append = append, row.names = FALSE, col.names = FALSE, quote = FALSE)
            close(out)
          })


# geneSetFile ------------------------------------------------------------------


#' geneSetFile
#'
#' Initialize \code{\link{geneSetFile-class}} object
#' @param path Path to geneSet file.
#' @param memlimit Maximum number of records to load at a time. 
#' @export
geneSetFile=function(path,memlimit=5000)
{
  con=gzfile(path,"r")
  sets=c()
  while (length(i <- readLines(con,n=memlimit)) > 0)
  {
    sets=c(sets,
           sapply(strsplit(i,split="\\|"),"[[",1))
  }
  close(con)
  new("geneSetFile", path=path, sets=sets)
}

setMethod("show", signature="geneSetFile",
          definition=function(object){
            message("rvat geneSetFile object")
            message(sprintf("Path:%s",object@path))
            message(sprintf("Sets:%s",length(object@sets)))
          })

setMethod("names", signature = "geneSetFile",
          definition = function(x) {
            x@sets
          })

setMethod("listGeneSets", signature = "geneSetFile",
          definition = function(object) {
            names(object)
          })

setMethod("length", signature = "geneSetFile",
          definition = function(x) {
            length(names(x))
          })


#' @describeIn geneSetFile-class getGeneSet
#' 
#' Extract a geneSet from a \code{\link{geneSetFile-class}} object.
#' 
#' @param object a \code{\link{geneSetFile-class}} object
#' @param geneSet a vector of geneSets to subset.
#' @export
setMethod("getGeneSet", signature="geneSetFile",
          definition=function(object, geneSet)
          {
            if (!all(geneSet %in% object@sets)) {warning("Not all specified geneSets are present in the geneSetFile, use `listGeneSets()` to check which geneSets are available.")}
            indices <- sort(which(object@sets %in% geneSet))
            set <- object@sets[indices]
            if(length(indices) > 1) indices[2:length(indices)] <- (dplyr::lead(indices)-indices)[1:(length(indices)-1)]
            con <- gzfile(object@path,"r")
            x <- lapply(indices,
                        FUN = function(i) {
                          i <- scan(con, skip = i-1, nlines = 1, what = "character", quiet = TRUE)
                          i <- unlist(strsplit(i,split="\\|"))
                          geneSet(
                            geneSetName = i[1],
                            units = i[2],
                            w = i[3],
                            metadata = i[4]
                          )
                        })
            x <- geneSetList(x)
            close(con)
            x
          })


#' @export
setMethod("as.geneSetList", signature="geneSetFile",
          definition=function(object)
          {
            con=gzfile(object@path,"r")
            setStart=0
            setEnd=length(object)
            x=vector(mode="list", length=setEnd)
            counter=1
            for (i in scan(con,skip=setStart, nlines=setEnd, what="character"))
            {
              i=unlist(strsplit(i,split="\\|"))
              x[[counter]]=geneSet(
                geneSetName = i[1],
                units = i[2],
                w = i[3],
                metadata = i[4]
              )
              counter=counter+1
            }
            close(con)
            return(geneSetList(geneSets = x))
          })


#' Build a geneSetList or geneSetFile
#'
#' Build a [`geneSetList`] or [`geneSetFile`] for use in gene set analyses ([`geneSetAssoc`] or [`assocTest-aggregateFile`])
#' Currently these can be build directly from GMT-files, data.frames and lists.
#'
#' @param data Can be 1) data.frame where the first column includes the names of geneSets,
#' the second column contains the genes included in the geneSet (comma-delimited).
#' The third and fourth column are optional: the third column contains weights (comma-delimited),
#' and the fourth column contains metadata. 
#' or 2) a list, where the names of the list represent the geneSet names and each element in the list contains
#' a vector with the genes included in the respective geneset.
#' @param gmtpath Path to a gmt-file
#' @param output Optional output file path (output will be gz compressed text).
#' Defaults to `NULL`, in which case a [`geneSetList`] is returned.
#' @param sep Separator used in input file. 
#' @export
buildGeneSet <- function(data=NULL, gmtpath=NULL, output=NULL, sep="\t") {
  if(!is.null(data)) {
    if(is.data.frame(data)) {
      ncol <- ncol(data)
      if(ncol > 4) {stop("data.frame should have at most 4 columns (geneSet names, units, weights, metadata)")}
      if(ncol == 3) data[,4] <- NA
      if(ncol == 2) data[,c(3,4)] <- NA
      
        sets <- lapply(1:nrow(data), 
          FUN = function(i, data) {
            geneSet(
              geneSetName = data[i, 1],
              units = data[i, 2],
              w = data[i,3],
              metadata = data[i,4]
            )
          },
          data=data)
        return(geneSetList(sets))
    } else if(is.list(data)) {
      classes <- unlist(lapply(data, FUN = class))
      if(all(classes == "character")) {
        sets <- lapply(
          1:length(data),
          FUN = function(i) {geneSet(
            geneSetName = names(data)[i],
            units = paste(data[[i]], collapse=","),
            w = NA,
            metadata = NA
          )
          }
        )
        return(geneSetList(sets))
      } else {
        stop("Each element in the list should be a character vector")
      }
    }
  }
    
  if(!is.null(gmtpath)) {
    genesets <- sort(readGMT(gmtpath, sep = sep))
  }
  if(!is.null(output)) {
    out <- gzfile(output,"w")
    write.table(as.data.frame(genesets), out, sep="|", row.names = FALSE,col.names = FALSE, quote = FALSE)
    close(out)
    message(sprintf("Generated geneSetFile: %s",output))
    return(geneSetFile(output))
  } else {
    return(genesets)
  }
}

# readers ----------------------------------------------------------------------

#' Read a gmt-file
#'
#' Read a gmt-file and return an object of class [`geneSetList`].
#' @param path File path
#' @param sep Delimiter used
#' @export
readGMT <- function(path, sep = "\t") {
  gmt <- readLines(path)
  gmt <- strsplit(gmt, sep)
  metadata <- unlist(lapply(gmt, FUN = function(x) x[2]))
  
  sets <- lapply(gmt, FUN = function(x) {
    geneSet(
      geneSetName = x[1],
      units = paste(x[3:length(x)], collapse = ","),
      w = paste(rep(1, times = (length(x)-2)), collapse = ","),
      metadata = x[2]
    )
  }
  )
  geneSetList(sets)
}

# geneSetAssoc -----------------------------------------------------------------

setMethod("checkDuplicates", signature = c("rvbResult"),
          definition = function(object, stop = TRUE) {
            if (length(unique(object$unit)) != nrow(object)) {
              if(stop)  {
                stop("Duplicated units are present in the rvbResult object")
              } else {
                warning("Duplicated units are present in the rvbResult object")
              }
            }
          })

# geneSetAssoc -----------------------------------------------------------------
# 
#' @export
setMethod("geneSetAssoc", signature=c("rvbResult"),
          definition=function(object,
                              geneSet = NULL,
                              scoreMatrix = NULL,
                              cormatrix = NULL,
                              condition = NULL,
                              covar = NULL,
                              test = c("lm", "mlm", "fisher", "ttest", "ztest", "ACAT"),
                              threshold = NULL,
                              Zcutoffs = NULL,
                              INT = FALSE,
                              scoreCutoffs = NULL,
                              minSetSize = 1,
                              maxSetSize = Inf,
                              oneSided = TRUE,
                              memlimit = 1000, 
                              ID = "unit",
                              output = NULL
          ) {
            # Check methods and available tests
            ## only 'lm' and 'mlm' are implemented for `scoreMatrix` object
            if(!is.null(scoreMatrix)) test <- test[test %in% geneSetAssoc_tests_score]
            
            # Check mlm / cormatrix
            # if ("mlm" %in% test) {
            #   warning("Unfortunately, `mlm` is currently not implemented yet, it will be soon!")
            #   test <- test[test!="mlm"]
            # }
            
            # Prepare data 
            if("mlm" %in% test) {
              
              if (oneSided) {
                warning("Note that currently mlm will return two-sided P-values, regardless of the `oneSided` argument.")
              }
              if("mlm" %in% test && is.null(cormatrix)) {
                stop("The mlm test requires specifying the 'cormatrix' argument, see the `buildCorMatrix` method.")
              }
              
              nullmodel <- fitNullModelGSA(object = object, 
                                           cormatrix = cormatrix, 
                                           covar = covar, 
                                           Zcutoffs = Zcutoffs)
              object <- getResults(nullmodel)
  
            } else {
              object <- .prepare_stats_GSA(
                object, covar, Zcutoffs, INT
              )
            }
            
            # Check `condition` parameter
            if (!is.null(condition)) {
              if ( is(condition, "geneSetList") || is(condition, "geneSetFile") ) {
                condition.type <- "geneSet"
              } else if ( is(condition, "matrix") ) {
                condition.type <- "matrix"
                
                ## checks
                nonoverlap <- sum(!object$unit %in% rownames(matrix))
                if(nonoverlap > 0) {
                  message(sprintf("%s/%s units in the results are not present in the condition matrix, these are excluded.",
                                  nonoverlap,
                                  nrow(object)
                  ))
                  object <- object[object$unit %in% rownames(matrix),]
                }
                nonoverlap <- sum(!rownames(matrix) %in% object$unit)
                if(nonoverlap > 0) {
                  message(sprintf("%s/%s units in the condition matrix are not present in the results, these are excluded",
                                  nonoverlap,
                                  nrow(object)))
                }
                
                # Make sure the order is correct
                matrix <- matrix[as.character(object$unit),,drop = FALSE]
                
                if(nrow(matrix) == 0 || nrow(object) == 0) {
                  stop("No overlapping units left to test!")
                }
                
              } else if ( is(condition, "vector")) {
                condition.type <- "vector"
                check <- sum(condition %in% names)
                if (check < length(condition)) {
                  message(sprintf("%s/%s specified in `condition` are present in the %s", check, length(condition), input.type))
                }
                condition <- condition[condition %in% names]
                
              } else {
                condition.type <- "vector"
                condition <- if(!is.null(geneSet)) {
                  names(geneSet) 
                } else if ( !is.null(scoreMatrix) ) {
                  colnames(scoreMatrix) 
                }
              }
            }
            
            ## scoreMatrix -----------------------------------------------------
            
            if(!is.null(scoreMatrix)) {
              # Check overlap between scoreMatrix and results file 
              # Check between results object and cormatrix
              nonoverlap <- sum(!object$unit %in% rownames(scoreMatrix))
              if(nonoverlap > 0) {
                message(sprintf("%s/%s units in the results are not present in the scoreMatrix, these are excluded.",
                                nonoverlap,
                                nrow(object)
                ))
                object <- object[object$unit %in% rownames(scoreMatrix),]
              }
              nonoverlap <- sum(!rownames(scoreMatrix) %in% object$unit)
              if(nonoverlap > 0) {
                message(sprintf("%s/%s units in the scoreMatrix are not present in the results, these are excluded",
                                nonoverlap,
                                nrow(object)))
              }
              
              # Make sure the order is correct
              scoreMatrix <- scoreMatrix[as.character(object$unit),,drop = FALSE]
              
              if(nrow(scoreMatrix) == 0 || nrow(object) == 0) {
                stop("No overlapping units left to test!")
              }
              
              if( !is.null(scoreCutoffs) ) {
                if ( length(scoreCutoffs) != 2 ) {
                  stop("`scoreCutoffs` should be a vector of length 2 (minimum and maximum sd)")
                }
                
                if ( any(scoreCutoffs < 0)) {
                  stop("scoreCutoffs should be a positive vector (the number of standard deviations below and above the mean)")
                }
                scoreMatrix <- apply(scoreMatrix, 
                                     2,
                                     function(x) {
                                       sdev <- sd(x)
                                       mn <- mean(x)
                                       x[x < (mn - (scoreCutoffs[1] * sdev))] <- (mn - (scoreCutoffs[1] * sdev))
                                       x[x > (mn + (scoreCutoffs[2] * sdev))] <- (mn + (scoreCutoffs[2] * sdev))
                                       x
                                     }
                )
              }
              
              chunks <- split(1:ncol(scoreMatrix), ceiling(seq_along(1:ncol(scoreMatrix))/memlimit))
              result_list <- list()
              i <- 1
              for(chunk in chunks) {
                scoreMatrix_chunk <- scoreMatrix[,chunk,drop=FALSE]
                
                result_list[[i]] <- enrich_test(
                  as.data.frame(object),
                  scorematrix=scoreMatrix_chunk,
                  nullmodel=if("mlm" %in% test) nullmodel else NULL,
                  covar=covar,
                  test=test,
                  oneSided=oneSided,
                  ID = ID
                )
                i <- i + 1
              }
              
              result_list <- do.call(rbind, result_list)
              
              if(!is.null(output)) {
                write.table(result_list, sep="\t", quote = FALSE, file = gzfile(output), row.names = FALSE)
              } else {
                return(result_list)
              }
            }
            
            ## GSA -------------------------------------------------------------
            if(sum(geneSetAssoc_tests_competitive_threshold %in% test) > 0) {
              if(is.null(threshold)) {
                threshold <- 0.05/nrow(object)
                message(sprintf("`threshold` not specified for defining significant genes, using a bonferroni threshold: %s", signif(threshold, 4)))
              }
            } else {
              threshold <- NULL
            }
            
            ## define chunks based on memlimit
            chunks <- split(1:length(geneSet), 
                            ceiling(seq_along(1:length(geneSet))/memlimit))
            result_list <- list()
            i <- 1
            
            for(chunk in chunks) {
              geneSetList_chunk <- getGeneSet(geneSet, geneSet = listGeneSets(geneSet)[chunk])
              mappedMatrix <- mapToMatrix(geneSetList_chunk, object, ID = ID, sparse = TRUE)
              
              if(minSetSize > 0 || maxSetSize < Inf) {
                keep <- Matrix::colSums(mappedMatrix)
                keep <- names(keep[keep >= minSetSize & keep <= maxSetSize])
              } else {
                keep <- names(geneSetList_chunk)
              }
              
              if(is.null(condition)) {
                if(length(keep) > 0) {
                  result_list[[i]] <- gsa(
                    as.data.frame(object),
                    genesetlist = getGeneSet(geneSetList_chunk, geneSet = keep),
                    mappedMatrix = as.matrix(mappedMatrix[,keep,drop=FALSE]),
                    nullmodel = if( "mlm" %in% test) nullmodel else NULL,
                    covar=covar,
                    test=test,
                    threshold=threshold,
                    oneSided=oneSided,
                    ID = ID
                  )
                } else {
                  result_list[[i]] <- setNames(data.frame(matrix(ncol = 12, nrow = 0)), 
                                               c("geneSetName", "test", "covar", "threshold", "geneSetSize", "genesObs",
                                                 "P", "effect", "effectSE", "effectCIlower", "effectCIupper"
                                               ))
                }
              } else {
                if( length(keep) > 0 ) {
                  result_list[[i]] <- gsa_conditional(
                    object,
                    condition = condition,
                    condition.type = condition.type,
                    genesetlist=getGeneSet(geneSetList_chunk, geneSet = keep),
                    mappedMatrix=as.matrix(mappedMatrix[,keep,drop=FALSE]),
                    geneSetFile = geneSetFile,
                    geneSetList = geneSetList,
                    nullmodel=if("mlm" %in% test) nullmodel else NULL,
                    covar=covar,
                    test=test,
                    threshold=threshold,
                    oneSided=oneSided,
                    ID = ID,
                    memlimit = memlimit
                  )
                } else {
                  result_list[[i]] <- setNames(data.frame(matrix(ncol = 12, nrow = 0)), 
                                               c("geneSetName", "test", "covar", "threshold", "geneSetSize", "genesObs",
                                                 "P", "effect", "effectSE", "effectCIlower", "effectCIupper", "condition"
                                               ))
                }
              }
              
              i <- i + 1
            }
            
            result_list <- do.call(rbind, result_list)
            result_gsa <- gsaResult(result_list)

            message(sprintf("%s out of %s sets are kept.", 
                            length(unique(result_list$geneSetName)), 
                            length(geneSet)
                            ))
            
            if(!is.null(output)) {
              write.table(result_gsa, sep="\t", quote = FALSE, file = gzfile(output), row.names = FALSE)
            } else {
              return(result_gsa)
            }
          }
)


gsa_conditional <- function(
  object,
  condition,
  condition.type,
  genesetlist,
  mappedMatrix,
  geneSetFile = NULL,
  geneSetList = NULL,
  nullmodel = NULL,
  covar = NULL,
  test = c("lm"),
  threshold = NULL,
  oneSided = TRUE,
  ID = "unit",
  memlimit = 1000
) {
  Pl=effectl=effectSEl=effectCIlowerl=effectCIupperl=list()
  Ngenes_available <- Matrix::colSums(mappedMatrix)
  
  if("lm" %in% test) {
    Plt=effectlt=effectSElt=effectCIlowerlt=effectCIupperlt=list()
    if (condition.type %in% c("geneSet", "vector")) {
      chunksCond <- split(1:length(condition), 
                      ceiling(seq_along(1:length(condition))/memlimit))
      if(condition.type == "geneSet") condNames <- names(condition)
      if(condition.type == "vector") condNames <- condition
    } else {
      chunksCond <- list(1:ncol(matrix))
      condNames <- colnames(matrix)
    }
    
    for ( chunk in chunksCond ) {
      
      if (condition.type == "geneSet") {
        
        geneSetList_chunk <- getGeneSet(condition, geneSet = names(condition)[chunk])
        mappedMatrixCond <- mapToMatrix(geneSetList_chunk, object, ID = ID, sparse = TRUE)
        
      } else if (condition.type == "vector") {
        names <- condition[chunk]
        if (!is.null(geneSet)) geneSetList_chunk <- getGeneSet(geneSet, geneSet = names)
        mappedMatrixCond <- mapToMatrix(geneSetList_chunk, object, ID = ID, sparse = TRUE)
        
      } else {
        mappedMatrixCond <- condition
      }
      
      for(geneset in colnames(mappedMatrixCond)) {
        stats <- cbind(as.data.frame(object), as.matrix(mappedMatrixCond)[,geneset,drop=FALSE])
        
        cov <- if(length(covar) == 0) geneset else c(covar, geneset)
        if (var(stats[,geneset]) == 0) cov <- cov[cov != geneset]
        
        if(length(cov) > 0) {
          X <- cbind(1,as.matrix(stats[,c(cov),drop=FALSE]))
        } else {
          X <- cbind(rep(1, nrow(stats)))
        }
        U1 <- crossprod(X, stats$Z)
        U2 <- solve(crossprod(X), U1)
        ytr <- stats$Z - X %*% U2
        U3 <- crossprod(X, as.matrix(mappedMatrix))
        U4 <- solve(crossprod(X), U3)
        Str <- mappedMatrix - X %*% U4
        effect <- as.vector(crossprod(ytr, Str)) / colSums(Str^2)
        Str2 <- colSums(Str^2)
        sig <- (sum(ytr^2) - effect^2 * Str2) / (nrow(stats)-ncol(X)-2)
        effectSE <- sqrt(sig * (1/Str2))
        
        if(oneSided) {
          P <- pt((effect/effectSE), nrow(stats) - ncol(X) - 1, lower.tail = FALSE)
          fac <- qt(0.95, df = nrow(stats) - ncol(X) - 1)
          effectCIlower <-  effect - (fac * effectSE)
          effectCIupper <-  rep(Inf, length(P))
        } else {
          P <- 2 * pt(abs(effect/effectSE), nrow(stats) - ncol(X) - 1, lower.tail = FALSE)
          fac <- qt(0.975, df = nrow(stats) - ncol(X) - 1)
          effectCIlower <-  effect - (fac * effectSE)
          effectCIupper <-  effect + (fac * effectSE)
        }
        Plt[[geneset]]=P;effectlt[[geneset]]=effect;effectSElt[[geneset]]=effectSE;effectCIlowerlt[[geneset]]=effectCIlower;effectCIupperlt[[geneset]]=effectCIupper
      }
    }
    Plt <- Plt[condNames]
    effectlt <- effectlt[condNames]
    effectSElt <- effectSElt[condNames]
    effectCIlowerlt <- effectCIlowerlt[condNames]
    effectCIupperlt <- effectCIupperlt[condNames]
    Pl[["lm"]]=unlist(Plt);effectl[["lm"]]=unlist(effectlt);effectSEl[["lm"]]=unlist(effectSElt);effectCIlowerl[["lm"]]=unlist(effectCIlowerlt);effectCIupperl[["lm"]]=unlist(effectCIupperlt)
  }
  
  Ngenes_available <- Matrix::colSums(mappedMatrix)
  results <- data.frame(
    geneSetName = c(rep(names(genesetlist), times = (length(test[test %in% geneSetAssoc_tests_competitive_condition]) * length(condNames)))),
    test = c(rep(geneSetAssoc_tests_competitive_condition[geneSetAssoc_tests_competitive_condition %in% test], each = (length(genesetlist) * length(condNames)))),
    covar = c(rep(paste(covar, collapse=","), times = (length(genesetlist) * length(condNames) * length(test[test %in% geneSetAssoc_tests_competitive_condition])))),
    condition = c(rep(condNames, each = (length(test[test %in% geneSetAssoc_tests_competitive_condition]) * length(genesetlist)))),
    threshold = c(rep(rep(NA_real_, times = (length(genesetlist) * length(condNames))))),
    geneSetSize = rep(lengths(genesetlist), times = (length(test[test %in% geneSetAssoc_tests_competitive_condition]) * length(condNames))),
    genesObs = c(rep(Ngenes_available, times = (length(test[test %in% geneSetAssoc_tests_competitive_condition]) * length(condNames)))),
    stringsAsFactors = FALSE, row.names = NULL
  )
  
  results <- cbind(results,
                   P = unname(c(Pl[["lm"]])),
                   effect = unname(c(effectl[["lm"]])),
                   effectSE = unname(c(effectSEl[["lm"]])),
                   effectCIlower = unname(c(effectCIlowerl[["lm"]])),
                   effectCIupper= unname(c(effectCIupperl[["lm"]]))
  )
  results
  
}


gsa <- function(stats,
                genesetlist,
                mappedMatrix,
                nullmodel = NULL,
                covar = NULL,
                test = c("lm", "mlm", "fisher", "ttest", "ztest", "ACAT"),
                threshold = NULL,
                oneSided = TRUE,
                ID = "unit"
) {
  Pl=effectl=effectSEl=effectCIlowerl=effectCIupperl=list()
  Ngenes_available <- Matrix::colSums(mappedMatrix)
  
  if(any(test %in% geneSetAssoc_tests_competitive)) {
    
    if("lm" %in% test) {
      ## based on: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-166
      if(length(covar) > 0) {
        X <- cbind(1,as.matrix(stats[,c(covar),drop=FALSE]))
      } else {
        X <- cbind(rep(1, nrow(stats)))
      }
      U1 <- crossprod(X, stats$Z)
      U2 <- solve(crossprod(X), U1)
      ytr <- stats$Z - X %*% U2
      U3 <- crossprod(X, as.matrix(mappedMatrix))
      U4 <- solve(crossprod(X), U3)
      Str <- mappedMatrix - X %*% U4
      effect <- as.vector(crossprod(ytr, Str)) / colSums(Str^2)
      Str2 <- colSums(Str^2)
      sig <- (sum(ytr^2) - effect^2 * Str2) / (nrow(stats)-ncol(X)-2)
      effectSE <- sqrt(sig * (1/Str2))
      
      if(oneSided) {
        P <- pt((effect/effectSE), nrow(stats) - ncol(X) - 1, lower.tail = FALSE)
        fac <- qt(0.95, df = nrow(stats) - ncol(X) - 1)
        effectCIlower <-  effect - (fac * effectSE)
        effectCIupper <-  rep(Inf, length(P))
      } else {
        P <- 2 * pt(abs(effect/effectSE), nrow(stats) - ncol(X) - 1, lower.tail = FALSE)
        fac <- qt(0.975, df = nrow(stats) - ncol(X) - 1)
        effectCIlower <-  effect - (fac * effectSE)
        effectCIupper <-  effect + (fac * effectSE)
      }
      
      Pl[["lm"]]=P;effectl[["lm"]]=effect;effectSEl[["lm"]]=effectSE;effectCIlowerl[["lm"]]=effectCIlower;effectCIupperl[["lm"]]=effectCIupper
    }
    
    if("mlm" %in% test) {
      P=effect=effectSE=effectCIupper=effectCIlower <- rep(NA_real_, length(genesetlist))
      
      if (is(getNullModel(nullmodel), "GENESIS.nullMixedModel") || is(getNullModel(nullmodel), "GENESIS.nullModel")) {
        
        mat <- GWASTools::MatrixGenotypeReader(
          genotype=t(as.matrix(mappedMatrix)), 
          snpID=as.integer(1:ncol(mappedMatrix)),
          chromosome=rep(1L, ncol(mappedMatrix)), 
          position=1:ncol(mappedMatrix), 
          scanID=1:nrow(mappedMatrix)
        )
        mat <- GWASTools::GenotypeBlockIterator(
          GWASTools::GenotypeData(mat,
                                  scanAnnot = GWASTools::ScanAnnotationDataFrame(
                                  data.frame(scanID = stats$scanID, stringsAsFactors = FALSE))), 
          snpBlock = ncol(mappedMatrix))
        res <- GENESIS::assocTestSingle(gdsobj = mat, null.model = getNullModel(nullmodel), test = c("Score"))
        ## double check
        if(!identical(res$variant.id, as.integer(1:ncol(mappedMatrix))))  {
          stop ("")
        }
        P <- res$Score.pval
        effect <- res$Est
        effectSE <- res$Est.SE
        effectCIlower <- NA_real_
        effectCIlupper <- NA_real_
      }
      
      Pl[["mlm"]]=P;effectl[["mlm"]]=effect;effectSEl[["mlm"]]=effectSE;effectCIlowerl[["mlm"]]=effectCIlower;effectCIupperl[["mlm"]]=effectCIupper
    }
    
    if("fisher" %in% test) {
      Plt=effectlt=effectSElt=effectCIlowerlt=effectCIupperlt=list()
      
      for(thresh in threshold) {
        
        ## multiple thresholds 
        stats$sig <- ifelse(stats$P < thresh, TRUE, FALSE)
        sig <- stats[[ID]][stats$sig]
        not_sig <- stats[[ID]][!stats$sig]
        P=effect=effectSE=effectCIupper=effectCIlower <- rep(NA_real_, length(genesetlist))
        
        for(i in 1:length(genesetlist)) {
          geneset <- names(genesetlist)[i]
          tryCatch(
            {
              pathway <- rownames(mappedMatrix)[mappedMatrix[,geneset]]
              not_pathway <- rownames(mappedMatrix)[!mappedMatrix[,geneset]]
              mat <- matrix(
                c(sum(sig %in% pathway),
                  sum(sig %in% not_pathway),
                  sum(not_sig %in% pathway),
                  sum(not_sig %in% not_pathway)),
                nrow = 2, byrow = TRUE
              )
              
              tst <- fisher.test(mat, alternative = if(oneSided) "greater" else "two.sided")
              
              effect[i] <- tst$estimate
              #effectSE[i] <- coef(sum)[geneset, 2]
              effectCIlower[i] <- tst$conf.int[1]
              effectCIupper[i] <-tst$conf.int[2]
              P[i] <-tst$p.value
            }
          )
        }
        Plt[[as.character(thresh)]]=P;effectlt[[as.character(thresh)]]=effect;effectSElt[[as.character(thresh)]]=effectSE;effectCIlowerlt[[as.character(thresh)]]=effectCIlower;effectCIupperlt[[as.character(thresh)]]=effectCIupper
      }
      
      Pl[["fisher"]]=unlist(Plt);effectl[["fisher"]]=unlist(effectlt);effectSEl[["fisher"]]=unlist(effectSElt);effectCIlowerl[["fisher"]]=unlist(effectCIlowerlt);effectCIupperl[["fisher"]]=unlist(effectCIupperlt)
    }
    
    results <- data.frame(
      geneSetName = c(rep(names(genesetlist), times = length(test[test %in% geneSetAssoc_tests_competitive_nothreshold])),
                      rep(names(genesetlist), times = (length(test[test %in% geneSetAssoc_tests_competitive_threshold]) * length(threshold)))),
      test = c(rep(geneSetAssoc_tests_competitive_nothreshold[geneSetAssoc_tests_competitive_nothreshold %in% test], each = length(genesetlist)),
               rep(geneSetAssoc_tests_competitive_threshold[geneSetAssoc_tests_competitive_threshold %in% test], each = (length(genesetlist) * length(threshold)))),
      covar = c(rep(paste(covar, collapse=","), times = (length(genesetlist) * length(test[test %in% geneSetAssoc_tests_competitive_nothreshold]))),
                rep(NA_character_, times = (length(genesetlist) * length(threshold) * length(test[test %in% c("fisher")])))),
      threshold = c(rep(rep(NA_real_, times = length(genesetlist)), times = length(test[test %in% geneSetAssoc_tests_competitive_nothreshold])),
                    rep(rep(as.character(threshold), each = length(genesetlist)), times = length(test[test %in% geneSetAssoc_tests_competitive_threshold]))),
      geneSetSize = c(rep(lengths(genesetlist), times = length(test[test %in% geneSetAssoc_tests_competitive_nothreshold])),
                 rep(lengths(genesetlist), times = (length(test[test %in% geneSetAssoc_tests_competitive_threshold]) * length(threshold)))),
      genesObs = c(rep(Ngenes_available, times = length(test[test %in% geneSetAssoc_tests_competitive_nothreshold])),
                           rep(Ngenes_available, times = (length(test[test %in% geneSetAssoc_tests_competitive_threshold]) * length(threshold)))),
      stringsAsFactors = FALSE, row.names = NULL
      
    )
    
    tst <- c(geneSetAssoc_tests_competitive_nothreshold, geneSetAssoc_tests_competitive_threshold)[c(geneSetAssoc_tests_competitive_nothreshold, geneSetAssoc_tests_competitive_threshold) %in% test]
    results_competitive <- cbind(results,
                                 P = unlist(Pl[tst]),
                                 effect = unlist(effectl[tst]),
                                 effectSE = unlist(effectSEl[tst]),
                                 effectCIlower = unlist(effectCIlowerl[tst]),
                                 effectCIupper= unlist(effectCIupperl[tst]))
    rownames(results_competitive) <- NULL
  }
  
  if ( any(test %in% geneSetAssoc_tests_selfcontained)) {
    
    if("ttest" %in% test) {
      P=effect=effectSE=effectCIupper=effectCIlower <- rep(NA_real_, length(genesetlist))
      Z <- stats$Z
      for(i in 1:length(genesetlist)) {
        geneset <- names(genesetlist)[i]
        Z.tmp <- Z[as.matrix(mappedMatrix[,geneset,drop=FALSE])[,1]]
        tryCatch(
          {
            if(oneSided) {
              fit <- t.test(Z.tmp, alternative = "greater")
            } else {
              fit <- t.test(Z.tmp, alternative = "two.sided")
            }
  
            effect[i] <- fit$estimate
            effectSE[i] <-  fit$estimate / fit$statistic
            effectCIlower[i] <- fit$conf.int[1]
            effectCIupper[i] <- fit$conf.int[2]
            P[i] <- fit$p.value
          }
        )
      }
      Pl[["ttest"]]=P;effectl[["ttest"]]=effect;effectSEl[["ttest"]]=effectSE;effectCIlowerl[["ttest"]]=effectCIlower;effectCIupperl[["ttest"]]=effectCIupper
    }
    
    if("ztest" %in% test) {
      P=effect=effectSE=effectCIupper=effectCIlower <- rep(NA_real_, length(genesetlist))
      Z <- stats$Z
      for(i in 1:length(genesetlist)) {
        geneset <- names(genesetlist)[i]
        Z.tmp <- Z[as.matrix(mappedMatrix[,geneset,drop=FALSE])[,1]]
        tryCatch(
          {
            if(oneSided) {
              fit <- z_test(Z.tmp, alternative = "greater")
            } else {
              fit <- z_test(Z.tmp, alternative = "two.sided")
            }
        
            effect[i] <- fit$mean
            effectSE[i] <- fit$SE
            effectCIlower[i] <- fit$CIlower
            effectCIupper[i] <- fit$CIupper
            P[i] <- fit$P
          }
        )
      }
      Pl[["ztest"]]=P;effectl[["ztest"]]=effect;effectSEl[["ztest"]]=effectSE;effectCIlowerl[["ztest"]]=effectCIlower;effectCIupperl[["ztest"]]=effectCIupper
    }
    
    
    if("ACAT" %in% test) {
      P=effect=effectSE=effectCIupper=effectCIlower <- rep(NA_real_, length(genesetlist))
      Pvals <- stats$P
      for(i in 1:length(genesetlist)) {
        geneset <- names(genesetlist)[i]
        P.tmp <- Pvals[as.matrix(mappedMatrix[,geneset,drop=FALSE])[,1]]
        tryCatch(
          {
            Pval = .rvat_ACAT(P.tmp)
            P[i] <- Pval
          }
        )
      }
      Pl[["ACAT"]]=P;effectl[["ACAT"]]=effect;effectSEl[["ACAT"]]=effectSE;effectCIlowerl[["ACAT"]]=effectCIlower;effectCIupperl[["ACAT"]]=effectCIupper
    }
    
    results <- data.frame(
      geneSetName = c(rep(names(genesetlist), times = length(test[test %in% geneSetAssoc_tests_selfcontained]))),
      test = c(rep(geneSetAssoc_tests_selfcontained[geneSetAssoc_tests_selfcontained %in% test], each = length(genesetlist))),
      covar = c(rep(NA_character_, times = (length(genesetlist) * length(test[test %in% geneSetAssoc_tests_selfcontained])))),
      threshold = c(rep(rep(NA_real_, times = length(genesetlist)), times = length(test[test %in% geneSetAssoc_tests_selfcontained]))),
      geneSetSize = c(rep(lengths(genesetlist), times = length(test[test %in% geneSetAssoc_tests_selfcontained]))),
      genesObs = c(rep(Ngenes_available, times = length(test[test %in% geneSetAssoc_tests_selfcontained]))),
      stringsAsFactors = FALSE, row.names = NULL
    )
    
    tst <- geneSetAssoc_tests_selfcontained[geneSetAssoc_tests_selfcontained %in% test]
    results_selfcontained <- cbind(results,
                     P = unlist(Pl[tst]),
                     effect = unlist(effectl[tst]),
                     effectSE = unlist(effectSEl[tst]),
                     effectCIlower = unlist(effectCIlowerl[tst]),
                     effectCIupper= unlist(effectCIupperl[tst]))
    rownames(results_selfcontained) <- NULL
  }
  
  if(any(test %in% geneSetAssoc_tests_selfcontained) && any(test %in% geneSetAssoc_tests_competitive)) {
    rbind(results_competitive, results_selfcontained)
  } else if(any(test %in% geneSetAssoc_tests_competitive)) {
    results_competitive
  } else if(any(test %in% geneSetAssoc_tests_selfcontained)) {
    results_selfcontained
  }
}

enrich_test <- function(stats,
                        scorematrix,
                        nullmodel = NULL,
                        covar = NULL,
                        test = c("lm", "mlm"),
                        oneSided = TRUE,
                        ID = "unit"
) {
  Pl=effectl=effectSEl=effectCIlowerl=effectCIupperl=list()
  
  if("lm" %in% test) {
    if(length(covar) > 0) {
      X <- cbind(1,as.matrix(stats[,c(covar),drop=FALSE]))
    } else {
      X <- cbind(rep(1, nrow(stats)))
    }
    U1 <- crossprod(X, stats$Z)
    U2 <- solve(crossprod(X), U1)
    ytr <- stats$Z - X %*% U2
    U3 <- crossprod(X, scorematrix)
    U4 <- solve(crossprod(X), U3)
    Str <- scorematrix - X %*% U4
    effect <- as.vector(crossprod(ytr, Str)) / colSums(Str^2)
    Str2 <- colSums(Str^2)
    sig <- (sum(ytr^2) - effect^2 * Str2) / (nrow(stats)-ncol(X)-2)
    effectSE <- sqrt(sig * (1/Str2))
    
    if(oneSided) {
      P <- pt((effect/effectSE), nrow(stats) - ncol(X) - 1, lower.tail = FALSE)
      fac <- qt(0.95, df = nrow(stats) - ncol(X) - 1)
      effectCIlower <-  effect - (fac * effectSE)
      effectCIupper <-  rep(Inf, length(P))
    } else {
      P <- 2 * pt(abs(effect/effectSE), nrow(stats) - ncol(X) - 1, lower.tail = FALSE)
      fac <- qt(0.975, df = nrow(stats) - ncol(X) - 1)
      effectCIlower <-  effect - (fac * effectSE)
      effectCIupper <-  effect + (fac * effectSE)
    }
    
    Pl[["lm"]]=P;effectl[["lm"]]=effect;effectSEl[["lm"]]=effectSE;effectCIlowerl[["lm"]]=effectCIlower;effectCIupperl[["lm"]]=effectCIupper
  }
  
  
  if("mlm" %in% test) {
    # ..
  }
  
  results <- data.frame(
    name = rep(colnames(scorematrix), times = length(test[test %in% geneSetAssoc_tests_score])),
    test = test[test %in% c("lm", "mlm")],
    covar = rep(paste(covar, collapse=","), times = (ncol(scorematrix) * length(test[test %in% geneSetAssoc_tests_score]))),
    genesObs = nrow(scorematrix),
    stringsAsFactors = FALSE, 
    row.names = NULL
  )
  
  results <- cbind(results,
                   P = c(Pl[["lm"]], Pl[["mlm"]]),
                   effect = c(effectl[["lm"]],  effectl[["mlm"]]),
                   effectSE = c(effectSEl[["lm"]], effectSEl[["mlm"]]),
                   effectCIlower = c(effectCIlowerl[["lm"]], effectCIlowerl[["mlm"]]),
                   effectCIupper= c(effectCIupperl[["lm"]], effectCIupperl[["mlm"]]))
  results
}


.prepare_stats_GSA <- function(object, covar, Zcutoffs, INT) {
  
  if(!is.null(Zcutoffs) && length(Zcutoffs) != 2) {
    stop("`Zcutoffs should be a vector of length 2 (minimum and maximum)")
  }
  
  if(!is.null(Zcutoffs) && INT) {
    warning("Z-score cutoffs are specified while INT = TRUE. Z-score cutoffs won't be applied.")
  }
  
  if(!all(covar %in% colnames(object))) {
    stop(sprintf("The following covariates are not present in the rvbResult: %s",
                 paste(covar[!covar %in% colnames(object)], collapse=",")))
  }
  
  # Check if there are duplicate units 
  if(sum(duplicated(object[["unit"]])) > 0) {
    stop("Units are duplicated; first filter the input so each row has a unique unit")
  }
  # Exclude rows with missing covariate values
  check <- complete.cases(as.data.frame(object)[,covar,drop=FALSE])
  if(sum(!check) > 0) {
    message(sprintf("%s row(s) are excluded because of missing covariate values.", sum(!check)))
    object <- object[check,]
  }
  
  # Exclude missing P-values
  if(sum(is.na(object$P)) > 0) {
    message(sprintf("%s P-values are missing, these are excluded.", sum(is.na(object$P))))
    object <- object[!is.na(object$P),]
  } 
  
  # Add Z-score
  object$Z <- qnorm(1-object$P)
  
  # Apply Z score cutoffs
  if(INT) {
    object$Z <- qnorm((rank(object$Z,na.last = "keep")-0.5)/sum(!is.na(object$Z)))
  } else if(!is.null(Zcutoffs)) {
    if(length(Zcutoffs) != 2) {stop("The length of `Zcutoffs` should be 2 (minimum and maximum).")}
    message(sprintf("%s Z-scores <%s are set to %s", sum(object$Z < Zcutoffs[1]), Zcutoffs[1], Zcutoffs[1]))
    message(sprintf("%s Z-scores >%s are set to %s", sum(object$Z > Zcutoffs[2]), Zcutoffs[2], Zcutoffs[2]))
    object$Z <- ifelse(object$Z < Zcutoffs[1], Zcutoffs[1], object$Z)
    object$Z <- ifelse(object$Z > Zcutoffs[2], Zcutoffs[2], object$Z)
  } else if (sum(is.infinite(object$Z)) > 0) {
    # If Z-score cutoffs are not specified, check if any Z-scores are ±infinite
    if(sum(is.infinite(object$Z) & object$Z < 0) > 0) {
      minZ <- min(object$Z[!is.infinite(object$Z)])
      message(sprintf("%s Z-scores are -Inf, these are set to the minimum observed Z-score: %s.", 
                      sum(is.infinite(object$Z)), 
                      signif(minZ, 4)))
      object$Z[is.infinite(object$Z) & object$Z < 0] <- minZ
    }
    if(sum(is.infinite(object$Z) & object$Z > 0) > 0) {
      maxZ <- max(object$Z[!is.infinite(object$Z)])
      message(sprintf("%s Z-scores are +Inf, these are set to the maximum observed Z-score: %s.", 
                      sum(is.infinite(object$Z)), 
                      signif(maxZ, 4)))
      object$Z[is.infinite(object$Z) & object$Z > 0] <- maxZ
    }
  }
  object
}


z_test <- function(x,mu=0, var=1, alternative = "two.sided") {
  se <- var/sqrt(length(x))
  b <- mean(x)
  list(
    mean = b - mu,
    SE = se, 
    P = if(alternative == "two.sided") 2*pnorm(abs((mean(x)-mu)/se), lower.tail = FALSE) else pnorm(((mean(x)-mu)/se), lower.tail = FALSE),
    CIlower = if(alternative == "two.sided") b - qnorm(0.975)*se else b - qnorm(0.95)*se,
    CIupper = if(alternative == "two.sided") b + qnorm(0.975)*se else Inf
  )
}