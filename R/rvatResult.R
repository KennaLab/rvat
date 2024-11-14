# Methods ---------------------------------------------------------------------

## reading & writing ----------------------------------------------------------

#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("writeResult", "rvatResult",
          function(object, 
                   file = "",
                   append = FALSE,
                   quote = FALSE,
                   sep = "\t",
                   eol = "\n",
                   na = "NA",
                   dec = ".",
                   col.names = !append,
                   qmethod = c("escape", "double"),
                   fileEncoding = "")
          {
            if (!append) {
              file <- file(file, "w") 
              
              # write metadata
              metadata <- metadata(object)
              metadata <- metadata[names(metadata) %in% metadata_rvatresult]
              .write_rvat_header(filetype = as.character(class(object)[1]),
                                 metadata = metadata, 
                                 con = file)
              
              # write results
              write.table(object, file = file, append = FALSE, quote = quote,
                          sep = sep, eol = eol, na = na, dec = dec,
                          row.names = FALSE, col.names=col.names, 
                          qmethod=qmethod,
                          fileEncoding=fileEncoding)
              
              close(file) 
            } else {
              write.table(object, file = file, append = append, quote = quote,
                          sep = sep, eol = eol, na = na, dec = dec,
                          row.names = FALSE, col.names=col.names, 
                          qmethod=qmethod,
                          fileEncoding=fileEncoding)
            }
          })


#' Read association results
#' 
#' Read results generated using the \code{\link{assocTest}} method.

#' @param path File path to data
#' @param header A logical value indicating whether the data contains a header. Defaults to `TRUE`.
#' @param type Result type ('singlevarResult', 'rvbResult', 'gsaResult').
#' Defaults to `NULL` in which case the result type is inferred.
#' @param sep The field separator. Defaults to `\\t`, which is the default separator using in \code{\link{assocTest}}.
#' @return An object of type \code{\link{rvbResult}} or \code{\link{singlevarResult}}.
#' @export
readResults <- function(path, header = TRUE, type = NULL, sep = "\t") {
  
  metadata <- .parse_rvat_header(path, 
                                 expected_metadata = metadata_rvatresult,
                                 expected_filetype = if (!is.null(type)) type else names(columns_rvatResult),
                                 n = length(metadata_rvatresult) + 1 # file description + metadata
                                 )
  
  if (header) {
    type <- checkClassrvatResult(object = path, type = type, sep = sep)
    dat <- read.table(file = path,
                      header = TRUE,
                      sep = sep,
                      stringsAsFactors = FALSE)
    
  } else {
    dat <- read.table(file = path,
                      header = FALSE,
                      sep = sep,
                      stringsAsFactors = FALSE)
    
    if (is.null(type)) {
      if (ncol(dat) == length(columns_singlevarResults)){
        type="singlevarResult"
        message("Result type inferred from number of columns: singlevarResult")
        colnames(dat) <- names(columns_singlevarResults)
      } else if(ncol(dat) == length(columns_rvbResults)) {
        type="rvbResult"
        message("Result type inferred from number of columns: rvbResult")
        colnames(dat) <- names(columns_rvbResults)
      } else if (ncol(dat) == length(columns_gsaResults)) {
        type="gsaResult"
        message("Result type inferred from number of columns: gsaResult")
        colnames(dat) <- names(columns_gsaResults)
      } else {
        stop ("Couldn't infer the result type. The result type can be specified by setting the `type` parameter.")
      }
    } else {
      if (type == "singlevarResult") {
        if (ncol(dat) < length(columns_singlevarResults)) {
          stop("Number of columns is smaller than expected for an singlevarResult")
        }
        
        if (ncol(dat) > length(columns_singlevarResults)) {
          message(sprintf("Number of columns is larger than the default, assuming that the first %s columns 
are the default singlevarResult columns.",length(columns_singlevarResults)))
          colnames(dat)[1:length(columns_singlevarResults)] <- names(columns_singlevarResults)
        } else {
          colnames(dat) <- names(columns_singlevarResults)
        }
      } else if (type == "rvbResult") {
        if (ncol(dat) < length(columns_rvbResults)) {
          stop("Number of columns is smaller than expected for an rvbResult")
        }
        if (ncol(dat) > length(columns_rvbResults)) {
          message(sprintf("Number of columns is larger than the default, assuming that the first %s columns 
are the default rvbResult columns.",length(columns_rvbResults)))
          colnames(dat)[1:length(columns_rvbResults)] <- names(columns_rvbResults)
        } else {
          colnames(dat) <- names(columns_rvbResults)
        }
      } else if (type == "gsaResult") {
        if (ncol(dat) < length(columns_gsaResults)) {
          stop("Number of columns is smaller than expected for an gsaResult")
        }
        
        if (ncol(dat) > length(columns_gsaResults)) {
          message(sprintf("Number of columns is larger than the default, assuming that the first %s columns 
are the default gsaResult columns.",length(columns_gsaResults)))
          colnames(dat)[1:length(columns_gsaResults)] <- names(columns_gsaResults)
        } else {
          colnames(dat) <- names(columns_gsaResults)
        }
      }
    }
  }
  
  result <-  rvatResult(dat, class = type)
  metadata(result) <- metadata
  result
}


## getters ----------------------------------------------------------
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("getGdbId", signature="rvatResult",
          definition=function(object){
            value <- metadata(object)$gdbId
            if (length(value) == 0) value <- NA_character_
            value
          }
)

#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("getGenomeBuild", signature="rvatResult",
          definition=function(object){
            value <- metadata(object)$genomeBuild
            if (length(value) == 0) value <- NA_character_
            value
          }
)

## combine -------------------------------------------------------------------
#' merge.rvatResult.data.frame
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("merge", c("rvatResult", "data.frame"), function(x, y, by, ...) {
  class <- checkClassrvatResult(x)
  nrows <- nrow(x)
  metadata <- metadata(x)
  x <- dplyr::left_join(as.data.frame(x), y, by = by)
  
  if (nrow(x) != nrows) {
    if (nrow(x) != nrows) {
      message("Rows have been added to the rvatResult object")
    }
  }
  
  x <- rvatResult(x, class = class)
  metadata(x) <- metadata
  x
})

#' merge.rvatResult.DataFrame
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("merge", c("rvatResult", "DataFrame"), function(x, y, by, ...) {
  class <- checkClassrvatResult(x)
  nrows <- nrow(x)
  metadata <- metadata(x)
  x <- dplyr::left_join(as.data.frame(x), y, by = by)
  
  if (nrow(x) != nrows) {
    if (nrow(x) != nrows) {
      message("Rows have been added to the rvatResult object")
    }
  }
  
  x <- rvatResult(x, class = class)
  metadata(x) <- metadata
  x
})

## misc -----------------------------------------------------------------------

# topResult
### topResults --------------------------------------------------------------

#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("topResult", "rvatResult",
          function(object, n = 10) {
            head(object[order(object$P),],n)
          })

# Check child class, for internal use
checkClassrvatResult <- function(object, type = NULL, sep = "\t") {
  
  if(is(object, "character")) {
    if(!is.null(type) && ! type %in% names(columns_rvatResult)) {
      stop(sprintf("Type should be one of the following: %s", 
                   paste(paste0("'", names(columns_rvatResult), "'"), collapse=","))) 
    }
    
    cols <- read.table(object, header = TRUE, nrows = 1, sep = sep)
    cols <- colnames(cols)
    
    overlap <- vector(mode = "numeric", length = length(columns_rvatResult))
    names(overlap) <- names(columns_rvatResult)
    for(i in names(columns_rvatResult)) {
      overlap[i] <- mean(names(columns_rvatResult[[i]]) %in% cols)
    }
    
    if(!is.null(type)) {
      if(overlap[type] == 1) {
        class <- type
      } else {
        stop(sprintf("The following columns are missing for results of type %s:
                   %s
             ", 
             type,
             paste(names(columns_rvatResult[[type]])[!names(columns_rvatResult[[type]]) %in% cols], collapse=",")
        ))
      }
    } else if ( all(overlap < 1) ) {
      msg <- "The following columns are missing, depending on the result type:"
      for(i in names(columns_rvatResult)) {
        msg <- paste(msg, "\n", sprintf("%s: %s", 
                                        i, 
                                        paste(names(columns_rvatResult[[i]])[!names(columns_rvatResult[[i]]) %in% cols], collapse = ",")
        ))
      }
      stop(msg)
    } else if (sum(overlap == 1) > 1) {
      stop(sprintf("Can't distinguish whether the results are of type: %s.
          This can be specified by setting the `type` parameter.", 
          paste(paste0("'", names(overlap)[overlap==1], "'"), collapse=",")))
    } else {
      class <- names(overlap[overlap == 1])
    }
  } else if(is(object, "rvatResult")) {
    if(is(object, "rvbResult")) {
      class <- "rvbResult"
    } else if (is(object, "singlevarResult")) {
      class <- "singlevarResult"
    } else if (is(object, "gsaResult")) {
      class <- "gsaResult"
    } else {
      class <- NULL
    }
  }
  class
}

# rvatResult constructor, for internal use
rvatResult <- function(object, class) {
  if (class == "rvbResult") {
    rvbResult(object)
  } else if (class == "singlevarResult") {
    singlevarResult(object)
  } else if (class == "gsaResult") {
    gsaResult(object)
  } else {
    stop("Class not known")
  }
}

## plots -----------------------------------------------------------------------
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("qqplot", c("rvatResult"),
          function(object, 
                   title = "", 
                   label = "label", 
                   threshold = NULL, 
                   showThreshold = TRUE, 
                   labelThreshold = NULL, 
                   cex = 16, 
                   lambda1000 = FALSE, 
                   case = NULL, 
                   control = NULL)
          {
            if (lambda1000 & (checkClassrvatResult(object) == "gsaResult")) {
              message("lambda1000 does not apply to a gsaResult object and is set to FALSE")
              lambda1000 <- FALSE
            }
            
            if(lambda1000 & (is.null(case) | is.null(control))) {
              case = max(object$caseN, na.rm = TRUE)
              control = max(object$ctrlN, na.rm = TRUE)
            }
            
            if (label != "label") {
              P <- as.data.frame(object)[c("P", label)]
            } else {
              P <- as.data.frame(object)["P"]
            }
            
            .qqplot(
              P = P,
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
          })


.qqplot <- function(P, 
                    title = "", 
                    label = "label", 
                    threshold = NULL, 
                    showThreshold = TRUE, 
                    labelThreshold = NULL, 
                    cex = 16, 
                    lambda1000 = FALSE, 
                    case = NULL, 
                    control = NULL) {
  if (is.null(threshold)){threshold=-log10(0.05/nrow(P))} else {threshold <- -log10(threshold)}
  if (is.null(labelThreshold)){labelThreshold=threshold} else {labelThreshold=-log10(labelThreshold)}
  
  P=P[!is.na(P$P),,drop=FALSE]
  X=data.frame(logp=-log10(P$P))
  chisq=qchisq(P$P[!(is.na(P$P))],1,lower.tail=FALSE)
  lambda <- median(chisq) / qchisq(0.5,1)
  if(lambda1000) lambda1000_ <- 1 + (lambda - 1) * (1 / case + 1 / control) * 500
  
  l=data.frame(xpos=Inf,
               ypos=-Inf,
               annotateText=if(lambda1000) {c(sprintf("\u03BB = %s",
                                                           signif(lambda,5)),
                                                   sprintf("\u03BB1000 = %s",signif(lambda1000_,5)))} else  {
                                                     sprintf("\u03BB = %s",signif(lambda,5))},
               hjust=1,
               vjust=if(lambda1000) c(-3,-1) else -1)

  if (ncol(P) > 1) {
    X['labels'] = unlist(P[label])
  }

  X=X %>% dplyr::arrange(logp)
  X$null=sort(-log10(runif(nrow(X))))
  mx <- max(X$null)


  qqplot=ggplot2::ggplot(X,ggplot2::aes(x=null,y=logp)) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(expression("-log"[10]*"(Null)")) + ggplot2::ylab(expression("-log"[10]*"(Obs)")) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(color="black"),
                   axis.line.y = ggplot2::element_line(color="black"),
                   axis.text.x=ggplot2::element_text(size=cex),
                   axis.text.y=ggplot2::element_text(size=cex),
                   axis.title.x=ggplot2::element_text(size=cex+2),
                   axis.title.y=ggplot2::element_text(size=cex+2)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept=0,slope=1,linetype=2,col="slateblue") +
    ggplot2::geom_text(data=l,ggplot2::aes(x=xpos,y=ypos,hjust=hjust,vjust=vjust,label=annotateText),parse=FALSE,size=cex/2.5)

  if(showThreshold) {
    qqplot <- qqplot + ggplot2::geom_segment(ggplot2::aes(x=x,xend=xend,y=y,yend=yend),
                                             col = "red", 
                                             linetype = 2,
                                             data = data.frame(x = 0, xend = mx, y = threshold, yend = threshold)
                                             )
  }

  if (ncol(X) == 3){qqplot = qqplot +
    ggrepel::geom_text_repel(data= (X[!is.na(X$logp) & X$logp > labelThreshold,,drop=FALSE]),
                             ggplot2::aes_string(label = "labels"), min.segment.length = 0, box.padding = 0.5)}
  return(qqplot)
}

#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("manhattan", c("rvatResult"),
          function(object, 
                   highlight = NULL, 
                   label = "label", 
                   threshold = NULL, 
                   labelThreshold = NULL, 
                   labelRepel = FALSE, 
                   labelSize = 3.88,
                   contigs = c(),
                   title = ""
                   )
          {
            # Check if 'CHROM' and 'POS' are present in the rvatResult object
            if(mean(c("CHROM", "POS") %in% colnames(object)) < 1) {
              if(mean(c("CHROM", "start", "end") %in% colnames(object)) == 1) {
                object$POS <- (object$start + object$end)/2
              } else if("CHROM" %in% colnames(object)) {
                stop("'POS' column is missing!")
              } else {
                stop("'CHROM' and 'POS' columns are missing!")
              }
            }
            class <- checkClassrvatResult(object)
            
            # configuration
            nvar=nrow(object)
            if (is.null(threshold)){threshold=0.05/nvar}
            if (is.null(labelThreshold)){labelThreshold=threshold}
            
            if (length(contigs) == 0) {
              message("Contigs not specified, defaulting to GRCh37")
              contigs <- "GRCh37"
            }
            if (is.character(contigs) && length(contigs) == 1) {
              if(contigs == "GRCh37") {
                contigs=data.frame(matrix(c(1,249250621,2,243199373,3,198022430,4,191154276,5,180915260,6,171115067,7,159138663,8,146364022,9,141213431,10,135534747,11,135006516,12,133851895,13,115169878,14,107349540,15,102531392,16,90354753,17,81195210,18,78077248,19,59128983,20,63025520,21,48129895,22,51304566, 23, 155270560, 24, 59373566),ncol=2,byrow=TRUE))
                names(contigs)=c("CHROM","Length")
              } else if (contigs == "GRCh38") {
                contigs=data.frame(matrix(c(1,248956422,2,242193529,3,198295559,4,190214555,5,181538259,6,170805979,7,159345973,8,145138636,9,138394717,10,133797422,11,135086622,12,133275309,13,114364328,14,107043718,15,101991189,16,90338345,17,83257441,18,80373285,19,58617616,20,64444167,21,46709983,22,50818468, 23, 156040895, 24, 57227415),ncol=2,byrow=TRUE))
                names(contigs)=c("CHROM","Length")
              } else {
                message("Currently only contigs for build GRCh37 and GRCh38 are implemented, please provide custom contigs for your build.")
              }
            }
            increment=c(0)
            for (i in 2:nrow(contigs)){increment[i]=increment[i-1]+contigs$Length[i-1]}
            contigs$increment=increment
            
            # recode sex chromosomes
            object$CHROM=ifelse(object$CHROM %in% c("X", "chrX"), 23, object$CHROM) 
            object$CHROM=ifelse(object$CHROM %in% c("Y", "chrY"), 24, object$CHROM) 
                                
            # Remap positions for contiguous manhattan
            object$CHROM=as.numeric(gsub("chr","",object$CHROM))
            
            object=object[object$CHROM %in% contigs$CHROM,]
            object$code=as.numeric(as.character(object$CHROM)) %% 2
            object$POS=object$POS+contigs$increment[match(object$CHROM,contigs$CHROM)]
            
            # Convert p-values to logp
            object$logp=-log10(object$P)
            if(label %in% colnames(object)) labelvec <- object[[label]] else labelvec <- NULL
            
            # Make + filter dataframe with points that need to be highlighted
            if (!is.null(highlight)) {
              if (class == "rvbResult") {
                highlight <- data.frame(object[object$unit %in% highlight,])
              } else if (class == "singlevarResult") {
                highlight <- data.frame(object[object$VAR_id %in% highlight,])
              }
              highlight <- highlight[c("POS", "logp")]
            }
            
            object <- data.frame(POS = object$POS, P=object$P, logp=object$logp, code=object$code)
            if(!is.null(labelvec)) object[[label]] <- labelvec
            
            mplot=ggplot2::ggplot(object,ggplot2::aes(x=POS,y=logp,color=code)) +
              ggplot2::xlab("") + ggplot2::ylab(expression("-Log"[10]*"(p)")) +
              ggplot2::theme(panel.background = ggplot2::element_blank(),
                             axis.line.x = ggplot2::element_line(color="black"),
                             axis.line.y = ggplot2::element_line(color="black"),
                             axis.text.x= ggplot2::element_text(size=14, angle=90),
                             axis.text.y= ggplot2::element_text(size=16),
                             axis.title.y= ggplot2::element_text(size=18)) +
              ggplot2::scale_x_continuous(limits = c(min(object$POS),max(object$POS)), 
                                          expand = c(0, 0),breaks=(contigs$increment + contigs$Length/2),
                                          labels=c(1:22, "X", "Y")) +
              ggplot2::scale_y_continuous(limits = c(0,max(object$logp)+0.75), 
                                          expand = c(0, 0)) +
              ggplot2::guides(colour = "none") +
              ggplot2::geom_segment(ggplot2::aes(x=x, 
                                                 xend=xend, 
                                                 y=y, 
                                                 yend=yend),
                                    col="grey", 
                                    linetype=2,
                                    data = data.frame(x = min(object$POS), xend = max(object$POS), y = -log10(threshold), yend = -log10(threshold))
                                    ) +
               ggplot2::geom_point() 
              
            if (label %in% names(object)){
              if(labelRepel) {
                mplot=mplot + 
                  ggrepel::geom_text_repel(ggplot2::aes(x=POS,y=logp,label=.data[[label]]),
                                     data=(object[!is.na(object$logp) & object$P < labelThreshold,]),
                                     size = labelSize 
                                     )
              } else {
                mplot=mplot + 
                  ggplot2::geom_text(ggplot2::aes(x=POS,y=logp,label=.data[[label]]),
                                     data=(object[!is.na(object$logp) & object$P < labelThreshold,]), 
                                     nudge_y = 0.7, 
                                     angle=0, 
                                     vjust = "inward",
                                     hjust = "inward",
                                     size = labelSize
                                     )
              }
              }
            
            if (title!=""){mplot=mplot+ggplot2::ggtitle(title)}
            
            if (!is.null(highlight)) {mplot = mplot +
              ggplot2::geom_point(data = highlight,ggplot2::aes(x=POS,y=logp), color = "red")
            }
            
            return(mplot)
          })

#' @rdname densityPlot
#' @usage NULL
#' @export
setMethod("densityPlot", 
          signature = signature(object="rvatResult"),
          function(object, geneSet, geneSetList, showMeans = FALSE, INT = FALSE, Zcutoffs = NULL, title = "") {
            if (checkClassrvatResult(object) == "gsaResult") {
              stop("gsaResult objects do not have per-test P-values, only per gene set. For 'object', select an rvbResult or singlevarResult object.")
            }
            
            object <- .prepare_stats_GSA(object, INT = INT, Zcutoffs = Zcutoffs, covar=NULL)
            units <- listUnits(getGeneSet(geneSetList, geneSet = geneSet))
            object$in_geneset <- object$unit %in% units
            
            if(showMeans) {
              means <- as.data.frame(object) %>%
                dplyr::group_by(in_geneset) %>%
                dplyr::summarize(mean = mean(Z,na.rm=TRUE))
            }
            
            #https://github.com/clauswilke/colorblindr/blob/master/R/palettes.R
            colorBlindFriendly <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")
            
            ggplot2::ggplot(as.data.frame(object), ggplot2::aes(x=Z, color=in_geneset, fill=in_geneset)) +
              ggplot2::ggtitle(title) +
              ggplot2::geom_vline(xintercept=0, linetype="dashed") +
              ggplot2::geom_density(alpha=0.6) + 
              {if(showMeans) ggplot2::geom_vline(data=means, mapping=ggplot2::aes(xintercept=mean))} + 
              ggplot2::theme(text=ggplot2::element_text(size=12)) +
              ggplot2::scale_colour_manual(values = colorBlindFriendly) +
              ggplot2::scale_fill_manual(values = colorBlindFriendly) +
              ggplot2::theme_classic()
          })

# rvbResult --------------------------------------------------------------------

# Constructor
proto_rvbResult <- function(n=0) {
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
}


#' rvbResult
#' @rdname rvatResult
#' @usage NULL
#' @export
rvbResult <- function(object, header = TRUE) {
  if(missing(object) || is.null(object)) {
    object <- proto_rvbResult()
  } 
  
  if (is.character(object)) {
    readResults(path = object, type = "rvbResult", header = header, sep = "\t")
  } else if ( is.data.frame(object) || is(object, "DFrame")) {
    if (!all(c("unit", "P") %in% colnames(object))) stop("Object must have at least 'unit' and 'P' columns.")
    
    if (!all(names(columns_rvbResults) %in% colnames(object))) {
      proto <- proto_rvbResult(n = nrow(object))
      proto[,intersect(colnames(proto), colnames(object))] <- object[,intersect(colnames(proto), colnames(object))]
      object <- cbind(proto, object[,!colnames(object) %in% colnames(proto), drop = FALSE])
    } else {
      object <- cbind(object[,names(columns_rvbResults)], object[,!colnames(object) %in% names(columns_rvbResults), drop = FALSE])
    }
    if(is.data.frame(object)) object <- S4Vectors::DataFrame(object)
    object[,columns_rvbResults_rle] <- lapply(object[,columns_rvbResults_rle], S4Vectors::Rle)
    object[,columns_rvbResults_numeric] <- lapply(object[,columns_rvbResults_numeric], as.numeric)
    object[["unit"]] <- as.character(object[["unit"]])
    new("rvbResult", object)
  } else {
    stop("`object` should be either a data.frame/DataFrame or a filepath.")
  }
}

# validity
setValidity2("rvbResult", function(object) {
  msg <- NULL
  
  if (!all(names(columns_rvbResults) %in% colnames(object)) && ncol(object) != 1) {
    msg <- c(msg, sprintf("The following columns are missing: %s",
                          paste(names(columns_rvbResults)[!names(columns_rvbResults) %in% colnames(object)], collapse=",")))
  }
  
  if(ncol(object) != 1) {
    columns <- names(columns_rvbResults)[names(columns_rvbResults) %in% colnames(object)]
    types <- unlist(lapply(object@listData, class))
    type_check <- unlist(lapply(columns, FUN = function(x) types[x] %in% columns_rvbResults[[x]]))
    names(type_check) <- columns
    if (!all(type_check) && ncol(object) != 1) {
      type_msg <- lapply(
        names(type_check)[!type_check], FUN = function(x) sprintf("%s: %s", x, paste(columns_rvbResults[[x]], collapse=" or "))
      ) %>% paste(collapse="\n")
      msg <- c(msg,sprintf("The following columns should be of type:\n%s", type_msg))
    }
  }
  
  if (is.null(msg)) {
    TRUE
  } else msg
})


## Methods ---------------------------------------------------------------------

### summary --------------------------------------------------------------------

#' summarize rvbResult
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("summary", "rvbResult",
          function(object, asList=FALSE, ...)
          {
            info <- list()
            info[["Nunits"]] <- length(unique(object[["unit"]]))
            
            info[["ctrlNmin"]] <- if((all(is.na(object[["ctrlN"]])))) NA_real_ else if (nrow(object) > 0) min(object[["ctrlN"]], na.rm=TRUE) else 0
            info[["ctrlNmax"]] <- if((all(is.na(object[["ctrlN"]])))) NA_real_ else if (nrow(object) > 0) max(object[["ctrlN"]], na.rm=TRUE) else 0
            info[["caseNmin"]] <- if((all(is.na(object[["caseN"]])))) NA_real_ else if (nrow(object) > 0) min(object[["caseN"]], na.rm=TRUE) else 0
            info[["caseNmax"]] <- if((all(is.na(object[["caseN"]])))) NA_real_ else if (nrow(object) > 0) max(object[["caseN"]], na.rm=TRUE) else 0
            info[["name"]] <- unique(object[["name"]])
            info[["varSetName"]] <- unique(object[["varSetName"]])
            info[["test"]] <- unique(object[["test"]])
            info[["covar"]] <- unique(object[["covar"]])
            info[["geneticModel"]] <- unique(object[["geneticModel"]])
            info[["MAFweight"]] <- unique(object[["MAFweight"]])
            info[["pheno"]] <- unique(object[["pheno"]])
            info[["cohort"]] <- unique(object[["cohort"]])
            
            if(asList) {
              return(info)
            } else {
              cat(sprintf("%s object\n-------------------\n", class(object)[1]))
              cat(sprintf("n units = %s\nN cases = %s\nN controls = %s\nnames = %s\nvarSets = %s\ntests = %s\ncovars = %s\ngeneticModels = %s\nMAFweights = %s\npheno = %s\ncohorts = %s", 
                          info$Nunits,
                          info$caseNmax, 
                          info$ctrlNmax, 
                          if(length(info[["name"]]) <= 5)  paste(info[["name"]], collapse=";") else sprintf("%s..+ %s", paste(info[["name"]][1:5], collapse=";"),length(info[["name"]])-5),
                          if(length(info[["varSetName"]]) <= 5)  paste(info[["varSetName"]], collapse=";") else sprintf("%s..+ %s", paste(info[["varSetName"]][1:5], collapse=";"),length(info[["varSetName"]])-5),
                          paste(info[["test"]], collapse=";"),
                          if(length(info[["covar"]]) <= 5)  paste(info[["covar"]], collapse=";") else sprintf("%s..+ %s", paste(info[["covar"]][1:5], collapse=";"),length(info[["covar"]])-5),
                          paste(info[["geneticModel"]], collapse=";"),
                          paste(info[["MAFweight"]], collapse=";"),
                          if(length(info[["pheno"]]) <= 5)  paste(info[["pheno"]], collapse=";") else sprintf("%s..+ %s", paste(info[["pheno"]][1:5], collapse=";"),length(info[["pheno"]])-5),
                          if(length(info[["cohort"]]) <= 5)  paste(info[["cohort"]], collapse=";") else sprintf("%s..+ %s", paste(info[["cohort"]][1:5], collapse=";"),length(info[["cohort"]])-5)
              ))
            }
          })

### ACAT -----------------------------------------------------------------------

#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("ACAT", c("rvatResult"),
          function(object,
                   aggregate = "test",
                   group = c("unit", "cohort", "varSetName","name", "pheno", "covar", "geneticModel", "MAFweight", "test"),
                   fixpval = TRUE,
                   fixpval_method = c("minmax", "manual", "Liu"),
                   fixpval_maxP = 0.99,
                   fixpval_minP = 1e-32,
                   warning = TRUE
          )
          {
            fixpval_method <- match.arg(fixpval_method)
            
            ## Input check
            if(is.character(aggregate)) {
              aggregate <- list(aggregate)
            } else if(!is.list(aggregate)) {
              stop("`aggregate` should be either a list or a character vector.")
            }
            aggregate_all <- unlist(aggregate)
            if(!all(aggregate_all %in% group)) {
              group <- c(group, aggregate_all[!aggregate_all %in% group])
            }
            
            ## Check if results already include ACAT metadata
            metadata <- metadata(object)
            if("ACAT" %in% metadata) {
              levels <- unlist(lapply(metadata[["ACAT"]], function(x) x[["level"]] ))
              acat_level <- max(levels) + 1
            } else {
              acat_level <- 1
              metadata[["ACAT"]] <- list()
            }
            
            for(agg in aggregate) {
              
              group <- group[!group %in% agg]
              
              # Add summary for metadata
              metadata[["ACAT"]][[paste0("ACAT", acat_level)]] <- list(
                level = acat_level,
                aggregate = agg,
                group=group,
                args = list(
                  fixpval = fixpval,
                  fixpval_method = fixpval_method,
                  fixpval_maxP = fixpval_maxP,
                  fixpval_minP = fixpval_minP
                ),
                summary = summary(object, asList=TRUE)
              )
              
              # Remove missing P-values
              if(sum(is.na(object[["P"]])) > 0) {
                if(warning) message(sprintf("Removing %s missing P-values.", sum(is.na(object[["P"]]))))
                object <- object[!is.na(object[["P"]]),]
              }
              
              ## Fix P-values that are exactly 1 or 0.
              # see: FAQ in https://github.com/yaowuliu/ACAT for other solutions
              if( fixpval && fixpval_method %in% c("minmax", "manual")) {
                if( fixpval_method == "minmax" ) {
                  P_max <- max(object[["P"]][object[["P"]] != 1], na.rm = TRUE)
                  P_min <- min(object[["P"]][object[["P"]] != 0], na.rm = TRUE)
                } else if ( fixpval_method == "manual" ) {
                  P_max <- fixpval_maxP
                  P_min <- fixpval_minP
                }
                nP1 <- sum(object[["P"]] == 1, na.rm = TRUE)
                nP0 <- sum(object[["P"]] == 0, na.rm = TRUE)
                if( nP1 > 0 ) {
                  if( warning ) warning(sprintf("%s P-values are exactly 1, these are set to %s", 
                                                sum(nP1, na.rm = TRUE), P_max))
                  object[["P"]] <- ifelse(object[["P"]] == 1, P_max, object[["P"]])
                }
                
                if( nP0 > 0 ) {
                  if(warning) warning(sprintf("%s P-values are exactly 0, these are set to %s", 
                                              sum(nP0, na.rm = TRUE), P_min))
                  object[["P"]] <- ifelse(object[["P"]] == 0, P_min, object[["P"]])
                }
              }
              
              dat <- cbind(as.data.frame(object)[,group,drop=FALSE], 
                           P = object[["P"]])
              count <- dat %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>%
                dplyr::count()
              
              if(all(count$n == 1)) {
                warning("All groupings have length 1")
              }
              
              if( fixpval && fixpval_method == "Liu" ) {
                dat <- dat %>%
                  dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>%
                  dplyr::mutate(
                    P = ifelse(P == 1, fixP_Liu(P, d = dplyr::n()), P)) %>%
                  dplyr::ungroup()
                
                nP0 <- sum(dat[["P"]] == 0, na.rm = TRUE)
                if( nP0 > 0 ) {
                  if(warning) warning(sprintf("%s P-values are exactly 0, these are set to %s", 
                                              sum(nP0, na.rm = TRUE), fixpval_minP))
                  dat[["P"]] <- ifelse(dat[["P"]] == 0, fixpval_minP, dat[["P"]])
                }
                
                acat <- dat %>%
                  dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>%
                  dplyr::summarize(P = .rvat_ACAT(P),
                                   .groups = "drop"
                  )
              } else {
                acat <- dat %>%
                  dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>%
                  dplyr::summarize(P = .rvat_ACAT(P),
                                   .groups = "drop"
                  )
              }
              
              ## Parse info
              object <- as.data.frame(object) %>% 
                #  dplyr::select(dplyr::)
                dplyr::group_by(dplyr::across(dplyr::all_of(group))) %>% 
                dplyr::summarize(
                  unit = if(length(unique(unit)) == 1) unit[1] else NA_character_,
                  cohort = if(length(unique(cohort)) == 1) cohort[1] else NA_character_,
                  varSetName = if(length(unique(varSetName)) == 1) varSetName[1] else NA_character_,
                  name = if(length(unique(name)) == 1) name[1] else NA_character_,
                  pheno = if(length(unique(pheno)) == 1) pheno[1] else NA_character_,
                  covar = if(length(unique(covar)) == 1) covar[1] else NA_character_,
                  geneticModel = if(length(unique(geneticModel)) == 1) geneticModel[1] else NA_character_,
                  MAFweight = if(length(unique(MAFweight)) == 1) MAFweight[1] else NA_character_,
                  test = if(length(unique(test)) == 1) test[1] else NA_character_,
                  nvar = if(length(unique(nvar)) == 1) nvar[1] else NA_real_,
                  caseCarriers = if(!all(is.na(caseCarriers)) & sum(!is.na(caseCarriers)) > 1 & var(caseCarriers, na.rm=TRUE) == 0) caseCarriers[1] else NA_real_,
                  ctrlCarriers = if(!all(is.na(ctrlCarriers)) & sum(!is.na(ctrlCarriers)) > 1 & var(ctrlCarriers, na.rm=TRUE) == 0) ctrlCarriers[1] else NA_real_,
                  meanCaseScore = if(!all(is.na(meanCaseScore)) & sum(!is.na(meanCaseScore)) > 1 & var(meanCaseScore, na.rm=TRUE) == 0) meanCaseScore[1] else NA_real_,
                  meanCtrlScore = if(!all(is.na(meanCtrlScore)) & sum(!is.na(meanCtrlScore)) > 1 & var(meanCtrlScore, na.rm=TRUE) == 0) meanCtrlScore[1] else NA_real_,
                  caseN = if(length(unique(caseN)) == 1) caseN[1] else NA_real_,
                  ctrlN = if(length(unique(ctrlN)) == 1) ctrlN[1] else NA_real_,
                  caseCallRate = if(!all(is.na(caseCallRate)) &  sum(!is.na(caseCallRate)) > 1 & var(caseCallRate, na.rm=TRUE) == 0) caseCallRate[1] else NA_real_,
                  ctrlCallRate = if(!all(is.na(ctrlCallRate)) & sum(!is.na(ctrlCallRate)) > 1 & var(ctrlCallRate, na.rm=TRUE) == 0) ctrlCallRate[1] else NA_real_,
                  .groups = "drop"
                ) %>%
                dplyr::left_join(acat, by = group)
              
              ## Add remaining columns
              object$effect <- object$effectSE <- object$effectCIlower <- object$effectCIupper <- NA_real_
              for (y in agg) {
                object[[y]] <- "ACAT"
              }
              
              ## Convert to rvbResult
              object <- rvbResult(object)
              acat_level <- acat_level + 1
            }
            metadata(object) <- metadata
            object
          }
)

fixP_Liu <- function(x, d) {
  if(d>1)(1 - 1/d) else(x)
}

# singlevarResults ------------------------------------------------------------

## Class def. and constructor --------------------------------------------------

# Class definition

proto_singlevarResult <- function(n=0) {
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
}


#' singlevarResult
#' @rdname rvatResult
#' @usage NULL
#' @export
singlevarResult <- function(object, header = TRUE) {
  if(missing(object) || is.null(object)) {
    object <- proto_singlevarResult()
  } 
  
  if (is.character(object)) {
    readResults(path = object, type = "singlevarResult", header = header, sep = "\t")
  } else if ( is.data.frame(object) || is(object, "DFrame")) {
    if (!all(c("VAR_id", "P") %in% colnames(object))) stop("Object must have at least 'VAR_id' and 'P' columns.")
    
    if (!all(names(columns_singlevarResults) %in% colnames(object))) {
      proto <- proto_singlevarResult(n = nrow(object))
      proto[,intersect(colnames(proto), colnames(object))] <- object[,intersect(colnames(proto), colnames(object))]
      object <- cbind(proto, object[,!colnames(object) %in% colnames(proto), drop = FALSE])
    } else {
      object <- cbind(object[,names(columns_singlevarResults)], object[,!colnames(object) %in% names(columns_singlevarResults), drop = FALSE])
    }
    if(is.data.frame(object)) object <- S4Vectors::DataFrame(object)
    object[,columns_singlevarResults_rle] <- lapply(object[,columns_singlevarResults_rle], S4Vectors::Rle)
    object[,columns_singlevarResults_numeric] <- lapply(object[,columns_singlevarResults_numeric], as.numeric)
    object[["VAR_id"]] <- as.character(object[["VAR_id"]])
    new("singlevarResult", object)
  } else {
    stop("`object` should be either a data.frame/DataFrame or a filepath.")
  }
}

# validity
setValidity2("singlevarResult", function(object) {
  msg <- NULL
  
  if (!all(names(columns_singlevarResults) %in% colnames(object)) && ncol(object) != 1) {
    msg <- c(msg, sprintf("The following columns are missing: %s",
                          paste(names(columns_singlevarResults)[!names(columns_singlevarResults) %in% colnames(object)], collapse=",")))
  }
  
  if(ncol(object) != 1) {
    columns <- names(columns_singlevarResults)[names(columns_singlevarResults) %in% colnames(object)]
    types <- unlist(lapply(object@listData, class))
    type_check <- unlist(lapply(columns, FUN = function(x) types[x] %in% columns_singlevarResults[[x]]))
    names(type_check) <- columns
    if (!all(type_check) && ncol(object) != 1) {
      type_msg <- lapply(
        names(type_check)[!type_check], FUN = function(x) sprintf("%s: %s", x, paste(columns_singlevarResults[[x]], collapse=" or "))
      ) %>% paste(collapse="\n")
      msg <- c(msg,sprintf("The following columns should be of type:\n%s", type_msg))
    }
  }
  
  if (is.null(msg)) {
    TRUE
  } else msg
})

## Methods ---------------------------------------------------------------------

### summary --------------------------------------------------------------------

#' summarize singlevarResult
#' @rdname rvatResult
#' @usage NULL
#' @export
setMethod("summary", "singlevarResult",
          function(object, asList=FALSE, ...)
          {
            info <- list()
            info[["Nvars"]] <- length(unique(object[["VAR_id"]]))
            
            info[["ctrlNmin"]] <- if((all(is.na(object[["ctrlN"]])))) NA_real_ else if (nrow(object) > 0) min(object[["ctrlN"]], na.rm=TRUE) else 0
            info[["ctrlNmax"]] <- if((all(is.na(object[["ctrlN"]])))) NA_real_ else if (nrow(object) > 0) max(object[["ctrlN"]], na.rm=TRUE) else 0
            info[["caseNmin"]] <- if((all(is.na(object[["caseN"]])))) NA_real_ else if (nrow(object) > 0) min(object[["caseN"]], na.rm=TRUE) else 0
            info[["caseNmax"]] <- if((all(is.na(object[["caseN"]])))) NA_real_ else if (nrow(object) > 0) max(object[["caseN"]], na.rm=TRUE) else 0
            info[["names"]] <- unique(object[["name"]])
            info[["varsets"]] <- unique(object[["varSetName"]])
            info[["covar"]] <- unique(object[["covar"]])
            info[["tests"]] <- unique(object[["test"]])
            info[["geneticModels"]] <- unique(object[["geneticModel"]])
            info[["pheno"]] <- unique(object[["pheno"]])
            info[["cohorts"]] <- unique(object[["cohort"]])
            
            if(asList) {
              return(info)
            } else {
              cat(sprintf("%s object\n-------------------\n", class(object)[1]))
              cat(sprintf("n vars = %s\nN cases = %s\nN controls = %s\nnames = %s\nvarSets = %s\ncovars = %s\ntests = %s\ngeneticModels = %s\npheno = %s\ncohorts = %s", 
                          info$Nvars,
                          info$caseNmax, 
                          info$ctrlNmax, 
                          if(length(info[["names"]]) <= 5)  paste(info[["names"]], collapse=";") else sprintf("%s..+ %s", paste(info[["names"]][1:5], collapse=";"),length(info[["names"]])-5),
                          if(length(info[["varsets"]]) <= 5)  paste(info[["varsets"]], collapse=";") else sprintf("%s..+ %s", paste(info[["varsets"]][1:5], collapse=";"),length(info[["varsets"]])-5),
                          if(length(info[["covar"]]) <= 5)  paste(info[["covar"]], collapse=";") else sprintf("%s..+ %s", paste(info[["covar"]][1:5], collapse=";"),length(info[["covar"]])-5),
                          paste(info[["tests"]], collapse=";"),
                          paste(info[["geneticModels"]], collapse=";"),
                          if(length(info[["pheno"]]) <= 5)  paste(info[["pheno"]], collapse=";") else sprintf("%s..+ %s", paste(info[["pheno"]][1:5], collapse=";"),length(info[["pheno"]])-5),
                          if(length(info[["cohorts"]]) <= 5)  paste(info[["cohorts"]], collapse=";") else sprintf("%s..+ %s", paste(info[["cohorts"]][1:5], collapse=";"),length(info[["cohorts"]])-5)
              ))
            }
          })
