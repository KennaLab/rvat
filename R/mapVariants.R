#' @rdname mapVariants 
#' @usage NULL
#' @export
setMethod("mapVariants", 
          signature="gdb",
          definition=function(
                        object, 
                        ranges = NULL, 
                        gff = NULL,
                        bed = NULL,
                        bedCols = character(),
                        fields = NULL, 
                        uploadName = NULL,
                        output = NULL, 
                        sep = "\t",
                        skipIndexes = FALSE,
                        overWrite=FALSE,
                        verbose = TRUE
) {
  
  ## Check if at least one input is given
  if ( sum(c(is.null(ranges), is.null(gff), is.null(bed))) == 0 ) {
    stop("At least one of `ranges`, `gff` or `bed` should be specified.")
  }
  
  ## Check if not more than 1 input is give
  if ( sum(c(!is.null(ranges), !is.null(gff), !is.null(bed))) != 1 ) {
    stop("Multiple inputs specified, specifiy either `ranges`, `gff`, or `bed`.")
  }
  
  ## ranges: can be a data.frame, GRanges object, or a file path
  if ( !is.null(ranges) ) {
    if ( is.character(ranges) ) {
      ranges <- read.table(ranges, sep = sep, stringsAsFactors = FALSE, header = TRUE)
    } 
    if ( is.data.frame(ranges) ) {
      if (!all(c("CHROM", "start", "end") %in% colnames(ranges))) {
        stop("'CHROM', 'start' and 'end' should be present in ranges-file.")
      }
      ranges <- GenomicRanges::makeGRangesFromDataFrame(ranges, keep.extra.columns = TRUE)
    } else if (!is(ranges, "GRanges")) {
        stop("ranges should be either a data.frame or a GRanges object")
    }
  }
  
  ## import bed-file using rtracklayer
  if ( !is.null(bed) ) {
    if ( length(bedCols) > 0 && is.null(names(bedCols)) ) {
      bedCols <- setNames(object = rep("character", length(bedCols)), nm = bedCols)
      }
    ranges <- rtracklayer::import(bed, extraCols = bedCols, format = "bed")
  }
  
  ## import gff/gtf-file using rtracklayer
  if ( !is.null(gff) ) {
    ranges <- rtracklayer::import(gff)
  }
  
  ## set seqlevelsStyle
  GenomeInfoDb::seqlevelsStyle(ranges) <- "NCBI"
  
  ## if uploadName != NULL, check if name is valid and if already exists in gdb
  if (!is.null(uploadName)) {
    if (!is.character(uploadName)) {stop("`uploadName` should be a character string (indicating how the uploaded table in the gdb should be named)")}
    if (grepl("\\.", uploadName)) {stop("Table name shouldn't contain '.'")}
    if (uploadName %in% gdb_protected_tables) {stop(sprintf("'%s' already exists as a protected table in gdb and cannot be replaced", uploadName))}
    
    cohort <- listCohort(object)
    anno <- listAnno(object)
    if (uploadName %in% c(cohort$name, anno$name)){
      if (!overWrite) {
        stop(sprintf("Table '%s' already exists. Please drop this table before continuing, or set `overWrite=TRUE`.", uploadName))
      } else {
        message(sprintf("Table '%s' already exists, it will be overwritten (as `overWrite=TRUE`)", uploadName))
        dropTable(object, uploadName)
      }
    }
  } 
  
  ## if output!=NULL connect 
  if(!is.null(output)) {
    output <- gzcon(file(output,open='wb'))
  } 
  
  ## if both output and uploadName are not specified, output will be returned
  if(is.null(uploadName) && is.null(output)) {
    container <- list()
  }
  
  ## fields to keep from ranges object
  if(!is.null(fields)) ranges <- ranges[,colnames(mcols(ranges)) %in% fields]
  
  chroms <- getAnno(object, "var_ranges", fields="CHROM")$CHROM
  
  ## map per chromosome
  for ( chrom in chroms ) {
    
    if(verbose) message(sprintf("Mapping chromosome %s..", chrom))
    
    ## get genomicRanges for current chromosome
    gr <- unserialize(getAnno(object, "var_ranges", where = sprintf("CHROM = '%s'", chrom))$ranges[[1]])
    GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"
    
    ## generate overlaps 
    overlaps <- GenomicRanges::findOverlaps(gr, ranges)
    dat <- cbind(VAR_id=gr[S4Vectors::queryHits(overlaps)]$VAR_id, 
                 as.data.frame(GenomicRanges::mcols(ranges[S4Vectors::subjectHits(overlaps),])))
    dat <- dplyr::arrange(dat, VAR_id)
    
    ## write to output (if specified)
    if (!is.null(output)) {
      if(chrom == chroms[1]) {
        write.table(dat, file = output, row.names = FALSE, quote = FALSE, sep = "\t")
      } else {
        write.table(dat, file = output, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
      }
    } 
    
    ## write to gdb (if specified)
    if (!is.null(uploadName)) {
      if(chrom == chroms[1]) {
        DBI::dbCreateTable(con = object, name = uploadName, fields = dat)
        DBI::dbAppendTable(con = object, name = uploadName, value = dat)
      } else {
        DBI::dbAppendTable(con = object, name = uploadName, value = dat)
      }
    }
    
    ## store in container (if output and uploadName are both NULL)
    if (is.null(output) && is.null(uploadName)) {
      container[[chrom]] <- dat
    }
  }

  ## close output
  if (!is.null(output)) {
    close(output)
  } 
  
  ## if output is written to gdb, update metadata, indexes
  if (!is.null(uploadName)) {
    fields <- DBI::dbListFields(object,uploadName)
    message(sprintf('%s fields detected (%s)\n',length(fields),paste(fields,collapse=",")))
    
    # Update annotation meta-data table
    anno <- listAnno(object)
    if (uploadName %in% anno$name){DBI::dbExecute(object,"delete from anno where name=:name",params=list(name=uploadName))}
    DBI::dbExecute(object,"insert into anno values (:name, :value, :date)",params=list(name=uploadName,value="mapVariants",date=date()))
    
    # Index
    if (!skipIndexes)
    {
      if ("VAR_id" %in% fields){DBI::dbExecute(object,sprintf("create index %s_idx1 on %s (VAR_id)",uploadName,uploadName))}
      if (sum(c("CHROM","POS","REF","ALT") %in% fields)==4){DBI::dbExecute(object,sprintf("create index %s_idx2 on %s (CHROM, POS, REF,ALT)",uploadName,uploadName))}
    }}
  
  if ( is.null(output) && is.null(uploadName) ) {
    container <- do.call(rbind,container)
    rownames(container) <- NULL
    return(container)
   }
  }
)