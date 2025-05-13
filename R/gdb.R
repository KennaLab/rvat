# ===================================================================================
# gdb objects
# ===================================================================================

#' 
#' @rdname gdb
#' @usage NULL
#' @export
setMethod("show", signature="gdb",
          definition=function(object){
            cat("rvat gdb object\n",
                "Path:",object@dbname,"\n")})

# -----------------------------------------------------------------------------------
# Constructor

#' @rdname gdb
#' @usage NULL
#' @export
gdb=function(path)
{
  if(!file.exists(path)) stop(sprintf("'%s' doesn't exist.", path))
  tryCatch({con=DBI::dbConnect(DBI::dbDriver("SQLite"),path)}, error=function(e){stop(sprintf("Invalid gdb path '%s'",path))})
  new("gdb",con)
}

gdb_init=function(path)
{
  tryCatch({con=DBI::dbConnect(DBI::dbDriver("SQLite"),path)}, error=function(e){stop(sprintf("Invalid gdb path '%s'",path))})
  new("gdb",con)
}

# -----------------------------------------------------------------------------------
# Close gdb

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("close", signature="gdb",
          definition=function(con){
            DBI::dbDisconnect(con)
          }
)


# -----------------------------------------------------------------------------------
# Getters

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("listAnno", signature="gdb",
          definition=function(object){return(DBI::dbGetQuery(object,"select * from anno"))}
          )

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("listCohort", signature="gdb",
          definition=function(object){return(DBI::dbGetQuery(object,"select * from cohort"))})

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("getRvatVersion", signature="gdb",
          definition=function(object){
            value <- DBI::dbGetQuery(object,"select * from meta where name = 'rvatVersion'")$value
            if (length(value) == 0) value <- NA_character_
            value
            })

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("getGdbId", signature="gdb",
          definition=function(object){
            value <- DBI::dbGetQuery(object,"select * from meta where name = 'id'")$value
            if (length(value) == 0) value <- NA_character_
            value
  }
)

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("getGdbPath", signature="gdb",
          definition=function(object){
            object@dbname
          }
)

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("getCreationDate", signature="gdb",
          definition=function(object){
            value <- DBI::dbGetQuery(object,"select * from meta where name = 'creationDate'")$value
            if (length(value) == 0) value <- NA_character_
            value
          }
)

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("getGenomeBuild", signature="gdb",
          definition=function(object){
            value <- DBI::dbGetQuery(object,"select * from meta where name = 'genomeBuild'")$value
            if (length(value) == 0) value <- NA_character_
            value
            })



#' @rdname extractRanges 
#' @usage NULL
#' @export
setMethod("extractRanges", 
          signature="gdb",
          definition=function(object, ranges, padding = 250) {
            
            # add padding 
            if (!is.numeric(padding) || padding < 0) {
              stop("`padding` parameter should be a positive value")
            }
            
            # if ranges is a data.frame, it should contain 'CHROM', 'start' and 'end' columns
            if (is.data.frame(ranges) ) {
              if (!all(c("CHROM", "start", "end") %in% colnames(ranges))) {
                stop("'CHROM', 'start' and 'end' should be present in ranges-file.")
              }
            } else if (is(ranges, "GRanges"))  {
              
              # if ranges is a GRanges object, convert to a data.frame
              ranges <- data.frame(
                CHROM = as.character(seqnames(ranges)),
                start = start(ranges),
                end = end(ranges),
                stringsAsFactors = FALSE
              )
            } else {
              stop("ranges should be either a data.frame or a GRanges object")
            }
            
            # check whether supplied ranges uses chr.. notation in CHROM field
            # update based on format used in gdb if needed
            chroms <- unique(ranges$CHROM)
            chroms_var <-  unique(getAnno(object, "var_ranges", fields="CHROM")$CHROM)
            chrom_format_with_chr <- grepl("^chr", chroms_var[1]) 
            if (grepl("^chr", chroms[1])) {
              if (!chrom_format_with_chr) {
                ranges$CHROM <- gsub("chr", "", ranges$CHROM)
              }
            } else {
              if (chrom_format_with_chr) {
                ranges$CHROM <- paste0("chr", ranges$CHROM)
              }
            }
            chroms <- unique(ranges$CHROM)
            
            # if there is no chromosome overlap between ranges and gdb, raise a warning
            if(!any(chroms %in% chroms_var)) {
              warning ("None of the specified chromosomes are present in the gdb..")
              return(vector(mode = "character", length = 0))
            }
            
            # perform queries in chunks of 500, padding is added for variants where length(REF) > 1
            chunks <- split(1:nrow(ranges), ceiling(seq_along(1:nrow(ranges)) / 500))
            run_query <- function(x, ranges, gdb, padding) {
              query <- sprintf("(CHROM = '%s' and POS >= %s and POS <= %s)", 
                               ranges[x, ][["CHROM"]], 
                               ranges[x, ][["start"]] - padding, 
                               ranges[x, ][["end"]] + padding)
              query <- paste("select * from var where", paste(query, collapse = " or "))
              query <- RSQLite::dbGetQuery(gdb, query)
              query
            }
            queries <- do.call(rbind, lapply(chunks, FUN = run_query, ranges = ranges, gdb = object, padding = padding))
            rownames(queries) <- NULL
            
            # use GRanges to identify overlaps taking into account variants where length(REF) > 1
            gr <- GenomicRanges::GRanges(
              seqnames = queries$CHROM,
              ranges = IRanges::IRanges(
                start = queries$POS,
                end = queries$POS + stringr::str_length(queries$REF) - 1
              ),
              VAR_id = queries$VAR_id
            )
            names(gr) <- queries$VAR_id
            
            ranges <- GenomicRanges::GRanges(
              seqnames = ranges$CHROM,
              ranges = IRanges::IRanges(
                start = ranges$start,
                end = ranges$end
              )
            )
            overlaps <- GenomicRanges::findOverlaps(gr, ranges)
            
            # return overlapping VAR ids
            VAR_id <- as.character(gr[S4Vectors::queryHits(overlaps)]$VAR_id)
            VAR_id
          }
)

#' getAnno
#' @rdname getAnno
#' @usage NULL
#' @export
setMethod("getAnno", signature = "gdb",
          definition = function(object,table,fields="*",left=c(),inner=c(),VAR_id=NULL,ranges=NULL,padding = 250,where=c())
          {
            # Base query
            fields=paste(fields,collapse=",")
            query=sprintf("select %s from %s", fields, table)
            
            # ranges 
            if (!is.null(ranges) ) {
              VAR_id <- extractRanges(object, ranges = ranges, padding = padding)
              if (length(VAR_id) == 0) {
                message("No variants overlap with provided range")
              }
            }
            # Add left join operation
            for (i in left)
            {
              query=sprintf("%s left join %s using (VAR_id)",query,i)
            }
            # Add inner join operations
            for (i in inner)
            {
              query=sprintf("%s inner join %s using (VAR_id)",query,i)
            }
            # Build where statement
            if (!is.null(VAR_id))
            {
              VAR_id=sprintf("VAR_id in ('%s')",paste(VAR_id,collapse="','"))
              if (length(where)>0)
              {
                where=sprintf("(%s) AND (%s)",VAR_id,where)
              } else {where=VAR_id}
            } 
            if (length(where)>0){query=sprintf("%s where %s",query, where)}
            
            return(DBI::dbGetQuery(object,query))
          })


#' getCohort
#' @rdname getCohort
#'
#' @export
setMethod("getCohort", signature="gdb",
          definition=function(object, cohort, fields = "*", keepAll = FALSE)
            {
            # Base query
            fields=paste(fields,collapse=",")
            query=sprintf("select %s from %s", fields, cohort)
            cohort <- DBI::dbGetQuery(object, query)
            if (!keepAll) cohort <- cohort[!is.na(cohort[["IID"]]),]
            cohort
            })

#' getGT
#' @rdname getGT
#' @import GenomicRanges
#' @usage NULL
#' @export
setMethod("getGT", signature="gdb",
          definition=function(object, varSet = NULL, VAR_id = NULL, ranges = NULL, cohort = NULL, anno = NULL, annoFields = NULL, includeVarInfo = FALSE, checkPloidy = NULL, varSetName = "unnamed", unit = "unnamed", padding = 250, verbose = TRUE, strict = TRUE)
          {
            
            # Process varSet (if provided)/ranges/or VAR_id
            if( !is.null(varSet) ) {
              
              ## check if varSet was generated from the current gdb
              if (!is.null(varSet) && strict) {
                .check_gdb_ids(object, varSet, minVersion = "0.3.0")
              }
              
              if (is(varSet, "varSetList") && length(varSet) == 1) {
                ## check gdb id
                varSet <- varSet[[1]]
              }
              if(is(varSet, "varSet")) {
                VAR_id <- unlist(strsplit(varSet@VAR_id,split=","))
                w <- as.numeric(unlist(strsplit(varSet@w,split=",")))
                names(w) <- VAR_id
                unit <- varSet@unit
                varSetName <- varSet@varSetName
                
                ## check if duplicate VAR ids are present:
                if (length(unique(VAR_id)) < length(VAR_id)) {
                  stop("Duplicated VAR_ids are specified")
                }
                
              } else if (is(varSet, "varSetList")) {
                stop("A varSetList with >1 varSets is supplied. 
                     Only 1 varSet can be supplied to `getGT()`, use `collapseVarSetList` to merge the records of a varSetList")
              } else {
                stop("varSet parameter should be of class `varSet` or a `varSetList` of length 1")
              }
              
              ## if ranges are supplied extract those, set weights to 1
            } else if (!is.null(ranges)) {
              VAR_id <- extractRanges(object, ranges = ranges, padding = padding)
              w <- rep(1, length(VAR_id))
              names(w) <- VAR_id
              
              ## if VAR_ids are supplied, select those, set weights to 1
            } else if (!is.null(VAR_id)) {
              VAR_id <- unique(sort(as.integer(VAR_id)))
              w <- rep(1, length(VAR_id))
              names(w) <- VAR_id
            } else {
              stop("Specify either `varSet`, `VAR_id` or `ranges`")
            }

            
            # Process sample info
            
            ## if no cohort is supplied, load 'SM' table
            if ( is.null(cohort) ) {
              SM <- DBI::dbGetQuery(object,'select * from SM')
              cohort_name <- "SM"
              
            ## if cohort is supplied as either a data.frame or a DFrame, 
            ## infer cohort name from name of supplied object or retrieve from 
            ## metadata if class == DFrame
            } else if (is.data.frame(cohort) || is(cohort, "DFrame")) {
                if(is(cohort, "DFrame")) {
                  if("name" %in% names(metadata(cohort))) {
                    cohort_name <- metadata(cohort)$name 
                  } else {
                    cohort_name <- as.character(as.list(match.call())[-1][["cohort"]])
                  }
                  cohort <- as.data.frame(cohort)
                } else {
                    cohort_name <- as.character(as.list(match.call())[-1][["cohort"]])
                }
               
              ## check if cohort IIDs are present in gdb
              SM <- DBI::dbGetQuery(object,'select * from SM')
              SM$IID <- as.character(SM$IID)
              if (!all(cohort$IID %in% SM$IID)) {
                warning(
                  sprintf("%s/%s samples in the cohort are not present in the gdb.", 
                          sum(!unique(cohort$IID) %in% SM$IID), 
                          length(unique(cohort$IID))
                  )
                )
              }
              
              ## match cohort with SM IIDs
              cohort <- cohort[match(SM$IID,cohort$IID),,drop = FALSE]
              SM <- cohort
            } else {
              SM <- getCohort(object, cohort, keepAll = TRUE)
              cohort_name <- cohort
            }
            
            # Verify ploidy settings
            if (is.null(checkPloidy)) {
              build_gdb <- getGenomeBuild(object)
              if (build_gdb %in% names(nonPAR)) checkPloidy <- build_gdb 
            }
            
            if (is.null(checkPloidy))
            {
              ploidy <- "diploid"
            } else {
              if (!checkPloidy %in% names(nonPAR))
              {stop(sprintf("Unrecognized value for checkPloidy. Must be one of %s. 
                            Alternatively you can provide no value for checkPloidy and manually reset values from the default of 'diploid' using the resetPloidy function",paste(names(nonPAR),collapse=",")))}
              
              # get varinfo and overlap with nonPAR regions
              varInfo=DBI::dbGetQuery(object,
                                      sprintf("select VAR_id, 'diploid' as ploidy, CHROM, POS, REF from var where VAR_id in ('%s');",
                                              paste(VAR_id,collapse="','")))
              varInfo$end=varInfo$POS+nchar(varInfo$REF)
              varInfo=makeGRangesFromDataFrame(varInfo, seqnames.field = "CHROM", start.field = "POS", end.field="end", keep.extra.columns = TRUE)
              GenomeInfoDb::seqlevelsStyle(varInfo)="NCBI"
              withCallingHandlers(
                overlaps <- findOverlaps(varInfo, nonPAR[[checkPloidy]]), 
                warning = suppressSeqinfowarning)
              
              if (length(overlaps) > 0) {
                if (verbose) message(sprintf("Ploidy of non-pseudoautosomal regions of the sex chromosomes are being set based on build %s", checkPloidy))
              }
              ## update ploidy for variants that overlap nonPAR regions
              for(i in 1:length(overlaps))
              {
                varInfo$ploidy[queryHits(overlaps)[i]]=nonPAR[[checkPloidy]]$ploidy[subjectHits(overlaps)[i]]
              }
              ploidy=data.frame(VAR_id=varInfo$VAR_id, ploidy=varInfo$ploidy)
            }
            
            # query genotypes
            
            ## extract genotypes
            query <- sprintf("select VAR_id, GT from dosage where VAR_id in ('%s');",
                             paste(VAR_id,collapse="','"))
            GT <- RSQLite::dbGetQuery(object, query)
            success <- sum(VAR_id %in% GT[,1])
            if ( sum(success) < length(VAR_id) ) {warning(sprintf("Retrieved genotypes for only %s of %s variants",success,length(VAR_id)))}
            nvar <- nrow(GT)
            VAR_id <- GT[,1]
            w <- w[as.character(VAR_id)]
            if (is(ploidy, "data.frame")){ploidy <- ploidy$ploidy[match(VAR_id, ploidy$VAR_id)]}
            
            # if specified, extract annotations
            if (includeVarInfo) {
              anno <- getAnno(object, table = "var")
              if (sum(duplicated(anno$VAR_id)) > 0) stop("Annotations should contain 1 row per VAR_id to be included in a genoMatrix object.")
            } else if(!is.null(anno)) {
              anno <- getAnno(object, table = anno, fields = if (!is.null(annoFields)) annoFields else "*", VAR_id = VAR_id)
              if (sum(duplicated(anno$VAR_id)) > 0) stop("Annotations should contain 1 row per VAR_id to be included in a genoMatrix object.")
            }
            
            ## unpack genotypes and store in matrix
            GT <- as.integer(sapply(GT[,2],function(x){memDecompress(unlist(x),type="gzip",asChar=FALSE)}))-48
            GT[GT > 2] <- NA_real_
            GT <- matrix(GT,
                         nrow = nvar,
                         ncol = length(GT)/nvar,
                         byrow=TRUE)
            if(verbose) message(sprintf("Retrieved genotypes for %s variants",nvar))
            
            ## generate genoMatrix
            return(genoMatrix(GT=GT,SM=SM,anno=anno,VAR_id=VAR_id,w=w,ploidy=ploidy,varSetName=varSetName,unit=unit,cohortname=cohort_name,genomeBuild=getGenomeBuild(object),gdbpath=object@dbname,gdbid=getGdbId(object),verbose=verbose))
          })

#' @rdname subsetGdb
#' @usage NULL
#' @export
setMethod("subsetGdb", signature="gdb",
          definition=function(object, output, intersection, where, VAR_id, tables, skipIndexes, overWrite, verbose)
          {
            # Check for existence of output gdb
            if (file.exists(output))
              {
               if (overWrite){file.remove(output)} else {stop(sprintf("output gdb already exists '%s'. Must set overWrite=TRUE to replace existing files.",output))}
              }

            # table directory
            master=DBI::dbGetQuery(object,"select * from sqlite_master")
            tables.base=gdb_protected_tables[!gdb_protected_tables %in% c("var_ranges", "tmp")]
            tables.anno=DBI::dbGetQuery(object,"select name from anno")$name
            tables.cohort=DBI::dbGetQuery(object,"select name from cohort")$name
           
            for (i in DBI::dbListTables(object))
            {
              if (i %in% gdb_protected_tables){next}
              if (i %in% tables.anno){next}
              if (i %in% tables.cohort){next}
              warning(sprintf("Table '%s' is not a base table or listed in anno, cohort tables. This table will not be included in output gdb.",i))
            }
            
            if( !is.null(tables) ) {
              tables.anno <- tables.anno[tables.anno %in% tables]
              tables.cohort <- tables.cohort[tables.cohort %in% tables]
            }

            # create random output gdb handle
            tmp=paste(c("extract",sample(letters,28,replace=TRUE)),collapse="")
            while (TRUE)
            {
              if (!tmp %in% c(tables.base,tables.anno,tables.cohort)){break}
              tmp=paste(c("extract",sample(letters,28,replace=TRUE)),collapse="")
            }

            # Export new var table to output gdb
            DBI::dbExecute(object,sprintf("ATTACH '%s' as %s",output,tmp))
            query=sprintf("create table %s.var as select distinct var.* from var", tmp)
            if (!is.null(intersection)) {
              intersection <- unlist(strsplit(intersection,split=","))
            } else {
                intersection <- c()
              }
            for (i in intersection)
            {
              query=sprintf("%s inner join %s using (VAR_id)",query,i)
            }
            if ( !is.null(where) && !is.null(VAR_id) ){
              query <- sprintf("%s where %s and VAR_id in (%s)", query, where, sprintf("'%s'", paste(VAR_id, collapse = "','")))
            } else if ( !is.null(where) ) {
              query <- sprintf("%s where %s", query, where)
            } else if ( !is.null(VAR_id) ) {
              query <- sprintf("%s where VAR_id in (%s)", query, sprintf("'%s'", paste(VAR_id, collapse = "','")))
              }
            DBI::dbExecute(object, query)

            # Copy remaining base tables
            DBI::dbExecute(object,sprintf("create table %s.SM as select * from SM", tmp))
            DBI::dbExecute(object,sprintf("create table %s.dosage as select dosage.VAR_id,dosage.GT from dosage inner join %s.var using (VAR_id)", tmp, tmp))
            DBI::dbExecute(object,sprintf("create table %s.cohort as select * from cohort where name in (%s)", tmp, paste(paste0("'", tables.cohort,"'"),collapse=",")))
            DBI::dbExecute(object,sprintf("create table %s.anno as select * from anno where name in (%s)", tmp, paste(paste0("'", tables.anno,"'"),collapse=",")))

            # Copy cohort tables
            for (i in tables.cohort)
            {
              DBI::dbExecute(object,sprintf("create table %s.%s as select * from %s", tmp,i,i))
            }

            # Copy anno tables
            for (i in tables.anno)
            {
              DBI::dbExecute(object,sprintf("create table %s.%s as select %s.* from %s inner join %s.var using (VAR_id)", tmp, i, i, i, tmp))
            }

            # Copy indexes
            if (!skipIndexes) {
              for (copied in c(tables.base,tables.cohort,tables.anno))
              {
                for (index in (master %>% dplyr::filter(type=="index") %>% dplyr::filter(tbl_name==copied))$sql)
                {
                  DBI::dbExecute(object,gsub("CREATE INDEX ",sprintf("CREATE INDEX %s.",tmp),index))
                }
              }
             }
            
            # varRanges
            gdb <- gdb(output)
            addRangedVarinfo(gdb, overwrite=TRUE, verbose = verbose)
            
            # create meta table
            DBI::dbExecute(gdb,"create table meta (name text,value text)")
            
            # Add rvat version to meta table
            DBI::dbExecute(gdb, "insert into meta values (:name, :value)",
                           params = list(name = "rvatVersion", 
                                         value = as.character(packageVersion("rvat"))))
            
            # Add random identifier
            DBI::dbExecute(gdb,"insert into meta values (:name, :value)",
                           params = list(name = "id", 
                                         value = paste(sample(c(letters, 0:9), 28, replace = TRUE), collapse = "")))
            
            # Add genome build
            DBI::dbExecute(gdb,"insert into meta values (:name, :value)",
                           params = list(name = "genomeBuild", 
                                         value =  getGenomeBuild(object)))
            
            # add creation date
            DBI::dbExecute(gdb,"insert into meta values (:name, :value)",
                           params = list(name = "creationDate", 
                                         value = as.character(round(Sys.time(), units = "secs"))))
            

            message(sprintf("%s\tComplete", as.character(round(Sys.time(), units = "secs"))))
          })

#' writeVcf
#' @rdname writeVcf
#' @aliases writeVcf,gdb-method
#' @usage NULL
#' @export
setMethod("writeVcf", signature = "gdb",
          definition=function(object, 
                              output, 
                              VAR_id = NULL, 
                              IID = NULL, 
                              includeGeno = TRUE,
                              includeVarId = FALSE,
                              verbose = TRUE
                              )
            {

            # Open output connection
            tryCatch({output=gzfile(output,"w")}, error=function(e){stop(sprintf("Could not write to output path '%s'",output))})
            
            # what ID field to include
            id_field <- if (includeVarId) "VAR_id" else "ID"

            # Establish sample filtering rule
            SM=DBI::dbGetQuery(object, "select * from SM")
            if (is.null(IID))
              {
              if (verbose) message("No IID provided, all samples will be included in output vcf")
              IID=rep(TRUE,nrow(SM))
              } else
              {
              IID=SM$IID %in% IID
              }
            if (verbose) message(sprintf("%s/%s samples to be retained in output", sum(IID), nrow(SM)))

            # Variant retrieval queries with / without genotype data
            if (includeGeno)
              {
              recordQuery=sprintf("select CHROM, POS, %s, REF, ALT, QUAL, FILTER, '.', FORMAT, GT from var inner join dosage using (VAR_id)", id_field)
              } else
              {
              recordQuery=sprintf("select CHROM, POS, %s, REF, ALT, QUAL, FILTER, '.' from var", id_field)
              }

            # Send variant retrieval queries with / without VAR_id filtering
            if (is.null(VAR_id))
              {
              if (verbose) message("No VAR_id provided, all variants will be included in output vcf")
              records=DBI::dbSendQuery(object, recordQuery)
              } else
              {
                if (verbose) message(sprintf("%s variants to be retained in output",length(VAR_id)))
                records=DBI::dbSendQuery(object, paste(recordQuery,
                                                       sprintf("where VAR_id in ('%s')",paste(VAR_id,collapse="','"))))
              }

            # Write VCF
            write("##fileformat=VCFv4.2",output)

            if (includeGeno)
            {
              write(sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s",
                            paste(SM$IID[IID],collapse="\t")),output)
              while (!DBI::dbHasCompleted(records))
              {
                record=DBI::dbFetch(records,n=1)
                gt=unlist(
                    strsplit(
                      memDecompress(unlist(record$GT),type="gzip",asChar=TRUE)
                      ,split=""))
                gt[gt=="N"]="./."
                gt[gt=="0"]="0/0"
                gt[gt=="1"]="0/1"
                gt[gt=="2"]="1/1"
                record$GT=paste(gt[IID],collapse="\t")
                write.table(record,output,append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
              }
            } else
            {
              write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",output)
              while (!DBI::dbHasCompleted(records))
              {
                record=DBI::dbFetch(records,n=1)
                write.table(record,output,append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
              }
            }

            # Close connections
            DBI::dbClearResult(records)
            close(output)

          })

# -----------------------------------------------------------------------------------
# Setters

#' @rdname uploadAnno
#' @usage NULL
#' @export
setMethod("uploadAnno", signature="gdb",
          definition=function(object,name,value,sep,skipRemap,skipIndexes,ignoreAlleles,keepUnmapped,mapRef,verbose)
          {
            # check validity of table name:
            ## shouldn't be a protected table / no unsupported characters / shouldn't exist as cohort table
            if (name %in% gdb_protected_tables){stop(sprintf("'%s' already exists as a protected table in gdb and cannot be replaced", name))}
            if (grepl("\\+|\\-|\\.|\\,| ", name)) {stop("Table name shouldn't contain '.','+','-' or ','")}
            cohort=listCohort(object)
            if (name %in% cohort$name){stop(sprintf("%s already exists as a cohort table. Please drop this table before continuing", name))}
            
            # write table
            if (!is(value)[[1]] %in% c("character","data.frame")){stop("value must be a valid filename or data frame object")}
            if (is(value)[[1]]=="character")
            {
              if (verbose) message(sprintf("Loading table '%s' from '%s'\n",name, value))
              DBI::dbWriteTable(con=object,name=name,value=value,sep=sep,overwrite=TRUE)
            } else
            {
              if (verbose) message(sprintf("Loading table '%s' from interactive R session'\n",name))
              DBI::dbWriteTable(con=object,name=name,value=value,overwrite=TRUE)
              value="interactive_session"
            }
            fields=DBI::dbListFields(object,name)
            if (verbose) message(sprintf('%s fields detected (%s)\n',length(fields),paste(fields,collapse=",")))
            
            # map positions to VAR_id (if skipRemap=TRUE)
            if (!skipRemap)
            {
              if (!ignoreAlleles)
              {
                if (sum(c("CHROM","POS","REF","ALT") %in% fields)!=4){warning("Mapping not possible: CHROM, POS, REF, ALT not provided.")}
                DBI::dbExecute(object,sprintf("create index %s_idx2 on %s (CHROM, POS, REF,ALT)",name,name))
                
                if(keepUnmapped) {
                  DBI::dbExecute(object,sprintf("create table tmp as select %s.VAR_id, %s.* from %s left join %s using (CHROM,POS,REF,ALT)",mapRef,name,name,mapRef))
                  nullCount=DBI::dbGetQuery(object,sprintf("select count(1) as n from tmp where VAR_id is null"))$n[1]
                } else {
                  DBI::dbExecute(object,sprintf("create table tmp as select %s.VAR_id, %s.* from %s inner join %s using (CHROM,POS,REF,ALT)",mapRef,name,name,mapRef))
                  nullCount=DBI::dbGetQuery(object, sprintf("select count(1) as n from %s", name))$n[1] - 
                    DBI::dbGetQuery(object, "select count(1) as n from tmp")$n[1] 
                  
                }
                
              } else
              {
                if (sum(c("CHROM","POS") %in% fields)!=2){stop("Mapping not possible: CHROM, POS not provided.")}
                DBI::dbExecute(object,sprintf("create index %s_idx2 on %s (CHROM, POS)",name,name))
                
                if(keepUnmapped) {
                  DBI::dbExecute(object,sprintf("create table tmp as select %s.VAR_id, %s.* from %s left join %s using (CHROM,POS)",mapRef,name,name,mapRef))
                  nullCount=DBI::dbGetQuery(object,sprintf("select count(1) as n from tmp where VAR_id is null"))$n[1]
                } else {
                  DBI::dbExecute(object,sprintf("create table tmp as select %s.VAR_id, %s.* from %s inner join %s using (CHROM,POS)",mapRef,name,name,mapRef))
                  nullCount=DBI::dbGetQuery(object, sprintf("select count(1) as n from %s", name))$n[1] - 
                    DBI::dbGetQuery(object, "select count(1) as n from tmp")$n[1] # not fully correct because rows may map to multiple variants
                }
              }
              
              DBI::dbExecute(object,sprintf("drop table %s",name))
              DBI::dbExecute(object,sprintf("alter table tmp rename to %s",name))
              if (nullCount > 0){warning(sprintf("Warning: %s rows could not be mapped to variants in the gdb\n", nullCount))}
            }
            
            # index table
            if (!skipIndexes)
            {
              fields=DBI::dbListFields(object,name)
              if ("VAR_id" %in% fields){DBI::dbExecute(object,sprintf("create index %s_idx1 on %s (VAR_id)",name,name))}
              if (all(c("CHROM","POS","REF","ALT") %in% fields)){DBI::dbExecute(object,sprintf("create index %s_idx2 on %s (CHROM, POS, REF,ALT)",name,name))} else if(all(c("CHROM","POS") %in% fields)) {
                DBI::dbExecute(object,sprintf("create index %s_idx2 on %s (CHROM, POS)",name,name))
              }
            }
            
            # Update annotation meta-data table
            anno=listAnno(object)
            if (name %in% anno$name){DBI::dbExecute(object,"delete from anno where name=:name",params=list(name=name))}
            DBI::dbExecute(object,"insert into anno values (:name, :value, :date)",params=list(name=name,value=value,date=date()))
          })

#' @rdname uploadCohort
#' @usage NULL
#' @export
setMethod("uploadCohort", signature="gdb",
          definition = function(object,name,value,sep="\t",verbose=TRUE)
            {
            # check validity of table name:
            ## shouldn't be a protected table / no unsupported characters / shouldn't exist as anno table
            if (name %in% gdb_protected_tables){stop(sprintf("'%s' already exists as a protected table in gdb and cannot be replaced", name))}
            if (grepl("\\+|\\-|\\.|\\,| ", name)) {stop("Table name shouldn't contain '.','+','-' or ','")}
            anno=listAnno(object)
            if (name %in% anno$name){stop(sprintf("%s already exists as a variant annotation table. Please use drop this table before continuing", name))}

            # upload, mapping and sanity checks
            if (!is(value)[[1]] %in% c("character","data.frame")){stop("value must be a valid filename or data frame object")}
            if (is(value)[[1]]=="character")
              {
              if (verbose) message(sprintf("Loading cohort '%s' from '%s'\n",name, value))
              upload=read.table(value,sep=sep,as.is=TRUE,h=TRUE)
              } else
                {
                  if (verbose) message(sprintf("Loading cohort '%s' from interactive R session\n",name))
                  upload=value
                  value="interactive_session"
                }
            
            ## check and format IID,sex columns
            if (verbose) message(sprintf('%s fields detected (%s)\n', ncol(upload), paste(colnames(upload),collapse=",")))
            if (!("IID" %in% colnames(upload))){stop("No 'IID' field detected, aborting upload.")}
            if (!("sex" %in% colnames(upload))){stop("No 'sex' field detected, aborting upload.")}
            if (sum(duplicated(upload$IID)) > 0) {stop("The cohort contains duplicated IIDs.")}
            upload$sex[!(upload$sex %in% c(1,2))]=0
            if(verbose) {
              message(sprintf("%s males, %s females and %s unknown gender",
                            sum(upload$sex==1,na.rm=TRUE),
                            sum(upload$sex==2,na.rm=TRUE),
                            sum(upload$sex==0,na.rm=TRUE)))
            }
            
            ## upload
            SM=DBI::dbGetQuery(object,'select * from SM')
            SM$IID=as.character(SM$IID)
            m=nrow(upload)
            upload=upload[match(SM$IID,upload$IID),]
            if (is.null(nrow(upload)) | is.null(ncol(upload))){stop("Table filtering error")}
            if (verbose) message(sprintf("Retained %s of %s uploaded samples that could be mapped to dosage matrix",m,nrow(upload)))
            DBI::dbWriteTable(con=object,name=name,value=upload,overwrite=TRUE)

            # Update cohort meta-data table
            cohort=DBI::dbGetQuery(object,"select * from cohort")
            if (name %in% cohort$name){DBI::dbExecute(object,"delete from cohort where name=:name",params=list(name=name))}
            DBI::dbExecute(object,"insert into cohort values (:name, :value, :date)",params=list(name=name,value=value,date=date()))
            })

#' @rdname gdb
#' @usage NULL
#' @export
setMethod("dropTable", signature="gdb",
          definition=function(object, name, verbose = TRUE)
          {
            anno=listAnno(object)
            cohort=listCohort(object)
            anno=anno[anno$name!=name,]
            cohort=cohort[cohort$name!=name,]
            DBI::dbWriteTable(con=object,name="anno",value=anno,overwrite=TRUE)
            DBI::dbWriteTable(con=object,name="cohort",value=cohort,overwrite=TRUE)
            DBI::dbExecute(object,sprintf("drop table %s",name))
            if(verbose) message(sprintf("Table '%s' removed from gdb",name))
          }
        )

## misc ------------------------------------------------------------------------
suppressSeqinfowarning <- function(w) if( any( grepl( "sequence levels", w) ) ) invokeRestart( "muffleWarning" )
