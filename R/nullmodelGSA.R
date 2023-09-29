## partition rvatResult for block-wise correlation matrix
## 
setMethod("addBlocks", signature = "rvatResult",
          definition = function(object, maxDist = 2.5e6) {
            # if already a genomicRanges objectif(..)
            
            # check if any are missing 
            check_missing <- !complete.cases(as.data.frame(object)[,c("CHROM", "start", "end")])
            if(sum(check_missing) > 0) {
              message(sprintf("%s/%s records have missing data for `CHROM`, `start` and/or `end`. These are excluded.", 
                              sum(check_missing), nrow(object)))
              object <- object[!check_missing,,drop=FALSE]
            }
            if(!is.numeric(object$start)) object$start <- as.numeric(as.character(object$start))
            if(!is.numeric(object$end)) object$end <- as.numeric(as.character(object$end))
            ranges <- GenomicRanges::GRanges(
              seqnames = paste0("chr", object$CHROM),
              ranges = IRanges::IRanges(
                start = object$start,
                end = object$end
              )
            )
            names(ranges) <- mcols(ranges)$unit <- object$unit
            ranges <- ranges[!duplicated(ranges$unit)]
            ranges <- GenomeInfoDb::sortSeqlevels(ranges)
            ranges <- GenomicRanges::sort(ranges)
            
            chroms <- unique(as.character(GenomicRanges::seqnames(ranges)))
            blocks <- vector(mode = "integer", length = length(ranges))
            names(blocks) <- names(ranges)
            
            for(chrom in chroms) {
              if(chrom == chroms[1]) current_block <- 1
              ranges_ <- ranges[GenomicRanges::seqnames(ranges) == chrom]
              if(length(ranges_) == 1) {
                blocks[names(ranges_)] <- current_block
                current_block <- current_block + 1L
              } else {
                dist <- IRanges::distance(ranges_[1:(length(ranges_)-1)], ranges_[2:length(ranges_)])
                gaps <- which(dist > maxDist)
                gaps <- c(gaps, length(ranges_))
                if(length(gaps) > 0) {
                  for(i in 1:length(gaps)) {
                    if(i == 1) {
                      index <- 1
                    } else {
                      index <- gaps[(i-1)]+1
                    }
                    blocks[names(ranges_)[index:gaps[i]]] <- current_block
                    current_block <- current_block + 1L
                  }
                } else {
                  blocks[names(ranges_)] <- current_block
                  current_block <- current_block + 1L
                }
              }
            }
            
            # Add to stats
            object <- merge(object, data.frame(
              unit = names(blocks),
              block = blocks
            ),by="unit")
            object
          })


## build block correlation matrix
## partly based on: https://github.com/opain/TWAS-GSEA/blob/master/TWAS-GSEA.V1.2.R
#' @rdname buildCorMatrix
#' @export
setMethod("buildCorMatrix", signature = c("rvatResult", "aggregateFile"),
          definition = function(object, 
                                aggregateFile, 
                                memlimit = 5000, 
                                minR2 = 1e-04, 
                                makePD = TRUE, 
                                absolute = TRUE, 
                                maxDist = 2.5e6,
                                verbose = TRUE
          ) {
            ## Checks
            if(minR2 < 0 || minR2 > 1) stop("`minR2` should be >=0 and <=1")
            
            if(mean(object$unit %in% listUnits(aggregateFile)) < 1) {
              message(sprintf("%s/%s units in the rvbResults object are not present in the aggregateFile, these will be excluded.", 
                              sum(!object$unit %in% listUnits(aggregateFile)),
                              nrow(object)
              ))
              object <- object[object$unit %in% listUnits(aggregateFile),]
            }
            
            ## block results based on specified maximum distance between blocks
            object <- addBlocks(object, maxDist = 2.5e6)
            
            corblocks <- list()
            units_per_block <- table(object$block)
            
            # generate chunks that will be loaded jointly based on memlimit
            chunks <- list()
            end <- TRUE
            i <- j <- 1
            while(end) {
              cmsum <- cumsum(units_per_block[i:length(units_per_block)])
              if(length(cmsum[cmsum <= memlimit]) == 0) {
                chunks[[j]] <- names(cmsum)[1]
              } else {
                chunks[[j]] <- names(cmsum[cmsum <= memlimit])
              }
              i <- as.numeric(tail(chunks[[j]], 1)) + 1
              j <- j + 1
              if(i >= length(units_per_block)) end <- FALSE
            }
            
            # calculate correlation per block
            chunk <- 1
            load_chunk <- "1"
            for(i in sort(unique(object$block))) {
              if(verbose) message(sprintf("Analyzing block: %s/%s", i, length(unique(object$block))))
              if(i == load_chunk) {
                if(verbose) message(sprintf("i = %s, loading new chunk", i))
                blocks <- chunks[[chunk]]
                load_chunk <- as.character(as.numeric(tail(blocks, 1))+1)
                units <- as.character(object$unit[object$block %in% blocks])
                mat <- getUnit(aggregateFile, unit = units)
                chunk <- chunk + 1
              }
              
              ## correlation if block consists of one gene (only diagonal)
              if(sum(object$block == i,na.rm=TRUE) == 1) {
                cor_block <- Matrix::Matrix(1, nrow = 1, ncol = 1, sparse = TRUE)
                colnames(cor_block)<- as.character(object$unit[object$block == i])
                rownames(cor_block)<- as.character(object$unit[object$block == i])
                corblocks[[i]] <- Matrix::Matrix(cor_block, sparse=TRUE)
              } else {
                
                ## calculate correlation using WGCNA::cor function
                if(absolute)  {
                  cor_block <- abs(WGCNA::cor(t(mat[as.character(object$unit[object$block == i]),,drop=FALSE]), method='pearson'))
                } else {
                  cor_block <- WGCNA::cor(t(mat[as.character(object$unit[object$block == i]),,drop=FALSE]), method='pearson')
                }
                
                ## set missing correlations to zero
                cor_block[is.na(cor_block)] <- 0
                
                ## Set correlations < minR2 to 0
                cor_block[abs(cor_block) < sqrt(minR2)] <- 0
                
                ## make positive definitive (if makePD = TRUE)
                if(makePD) {
                  if(!matrixcalc::is.positive.definite(as.matrix(cor_block))){
                    cor_block <- Matrix::nearPD(cor_block, corr = TRUE, maxit = 100)$mat
                  }
                }
                corblocks[[i]] <- Matrix::Matrix(cor_block, sparse=TRUE)
              }
            }
            names <- unlist(lapply(corblocks, FUN = function(x) rownames(x)))
            corMat <- Matrix::bdiag(corblocks)
            colnames(corMat) <- names
            rownames(corMat) <- names
            corMat
          })

## nullModelGSA onstructor
nullModelGSA <- function(nullmodel, ID, results, covar, method) {
  new("nullModelGSA",
      nullmodel = nullmodel,
      ID = ID,
      results = results,
      units = as.character(results$unit),
      covar = covar,
      method = method
  )
}

## methods

#' @rdname nullModelGSA-class
#' @param x \code{\link{nullModelGSA-class}} object
#' @export
setMethod("listUnits", signature = "nullModelGSA",
          definition = function(object) {
            object@units
          })


#' @rdname nullModelGSA-class
#' @param x \code{\link{nullModelGSA-class}} object
#' @export
setMethod("length", signature = "nullModelGSA",
          definition = function(x) {
            length(listUnits(x))
          })


#' @rdname nullModelGSA-class
#' @param object \code{\link{nullModelGSA-class}} object
#' @export
setMethod("show", 
          signature = "nullModelGSA",
          definition = function(object) {
            message(sprintf("nullModelGSA, generated with %s", object@method))
            message(sprintf("Contains a null model including %s units",length(object)))
            message(sprintf("Covar: %s", paste(getCovar(object), collapse=",")))
          })

#' @export
setMethod("getNullModel", signature = "nullModelGSA",
          definition = function(object) {
            object@nullmodel
          })


#' @export
setMethod("getResults", signature = "nullModelGSA",
          definition = function(object) {
            object@results
          })


#' @export
setMethod("getCovar", signature = "nullModelGSA",
          definition = function(object) {
            object@covar
          })

#' @rdname fitNullModelGSA
#' @export
setMethod("fitNullModelGSA", signature = "rvatResult",
          definition = function(object, 
                                cormatrix = NULL, 
                                covar = NULL, 
                                Zcutoffs = NULL, 
                                INT = FALSE, 
                                method = c("GENESIS"), 
                                ...) {
            
            method <- match.arg(method)
            
            # Prepare stats
            object <- .prepare_stats_GSA(
              object, covar, Zcutoffs, INT
            )
            
            # check if cormatrix is a matrix
            if (!is.null(cormatrix)) {
              # Check class of correlation matrix 
              if(!is(cormatrix, "Matrix") && !is(cormatrix, "matrix")) {
                stop("`cormatrix` should be a (sparse) matrix.")
              }
              
              ## Check overlap between results object and cormatrix
              nonoverlap <- sum(!object$unit %in% rownames(cormatrix))
              if(nonoverlap > 0) {
                message(sprintf("%s/%s units are not present in the correlation matrix, these are excluded.",
                                nonoverlap,
                                nrow(object)
                ))
                object <- object[object$unit %in% rownames(cormatrix),]
              }
              nonoverlap <- sum(!rownames(cormatrix) %in% object$unit)
              if(nonoverlap > 0) {
                message(sprintf("%s/%s units in the correlation matrix, are not present in the results, these are excluded",
                                nonoverlap,
                                nrow(object)
                ))
              }
              cormatrix <- cormatrix[as.character(object$unit), as.character(object$unit)]
            }
            
            # generate null model
            rownames(object)<- as.character(1:nrow(object))
            object$scanID <- as.integer(1:nrow(object))
            rownames(cormatrix) <- as.integer(1:nrow(object))
            colnames(cormatrix) <- as.integer(1:nrow(object))
            nullmodel <- GENESIS::fitNullModel(as.data.frame(object), 
                                               outcome = "Z", 
                                               covars = covar, 
                                               cov.mat = cormatrix, 
                                               family = "gaussian",
                                               ...)
            rownames(object) <- NULL
            
            # return null model 
            ID <- object$scanID
            names(ID) <- as.character(object$unit)
            nullModelGSA(nullmodel, 
                         ID = ID, 
                         results = object, 
                         covar = if(is.null(covar)) character(0) else covar, 
                         method = method)
          })

