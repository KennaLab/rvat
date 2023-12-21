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
  tryCatch({con=DBI::dbConnect(DBI::dbDriver("SQLite"),path)}, error=function(e){stop(sprintf("Invalid gdb path '%s'",path))})
  new("gdb",con)
}


# -----------------------------------------------------------------------------------
# Getters

#' @export 
setMethod("listAnno", signature="gdb",
          definition=function(object){return(DBI::dbGetQuery(object,"select * from anno"))}
          )

#' @export
setMethod("listCohort", signature="gdb",
          definition=function(object){return(DBI::dbGetQuery(object,"select * from cohort"))})

#' @export
setMethod("getRvatVersion", signature="gdb",
          definition=function(object){return(DBI::dbGetQuery(object,"select * from meta where name = 'rvatVersion'")$value)})

#' @rdname extractRanges 
#' @usage NULL
#' @export
setMethod("extractRanges", 
          signature="gdb",
          definition=function(object, ranges, return = "VAR_id") {
            if (is.data.frame(ranges) ) {
              if (!all(c("CHROM", "start", "end") %in% colnames(ranges))) {
                stop("'CHROM', 'start' and 'end' should be present in ranges-file.")
              }
              ranges <- GenomicRanges::makeGRangesFromDataFrame(ranges, keep.extra.columns = TRUE)
            } else if (!is(ranges, "GRanges")) {
              stop("ranges should be either a data.frame or a GRanges object")
            }
            chroms <- unique(as.character(seqnames(ranges)))
            chroms_var_ranges <-  getAnno(object, "var_ranges", fields="CHROM")$CHROM
            if(!any(chroms %in% chroms_var_ranges)) {
              warning ("None of the specified chromosomes are present in the gdb, you may have used the wrong chromosome specification, e.g. 'chr1' instead of '1', or vice versa.")
            }
            ## set seqlevelsStyle
            GenomeInfoDb::seqlevelsStyle(ranges) <- "NCBI"
            
            ## loop through chromosomes
            VAR_id <- list()
            for (chrom in chroms) {
              gr <- unserialize(getAnno(object, "var_ranges", where = sprintf("CHROM = '%s'", chrom))$ranges[[1]])
              GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"
              
              ## generate overlaps 
              overlaps <- GenomicRanges::findOverlaps(gr, ranges)
              VAR_id[[chrom]] <- gr[S4Vectors::queryHits(overlaps)]$VAR_id
            }
            VAR_id <- unique(unlist(VAR_id))
            VAR_id
          }
)

#' getAnno
#' @rdname getAnno
#' @usage NULL
#' @export
setMethod("getAnno", signature = "gdb",
          definition = function(object,table,fields="*",left=c(),inner=c(),VAR_id=c(),ranges=NULL,where=c())
          {
            # Base query
            fields=paste(fields,collapse=",")
            query=sprintf("select %s from %s", fields, table)
            
            # ranges 
            if (!is.null(ranges) ) {
              VAR_id <- extractRanges(object, ranges = ranges)
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
            if (length(VAR_id)>0)
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
          definition=function(object,cohort,fields="*")
            {
            # Base query
            fields=paste(fields,collapse=",")
            query=sprintf("select %s from %s", fields, cohort)
            return(DBI::dbGetQuery(object,query))
            })

#' getGT
#' @rdname getGT
#' @import GenomicRanges
#' @export
setMethod("getGT", signature="gdb",
          definition=function(object, varSet = NULL, VAR_id = NULL, ranges=NULL,cohort = NULL, checkPloidy = NULL, varSetName = "unnamed", unit = "unnamed", verbose=TRUE)
          {
            
            # Hardcoded ploidy configurations
            # https://www.ncbi.nlm.nih.gov/grc/human
            genomePloidy=list(
              hg19=GRanges(seqnames=c("chrX","chrX","chrX","chrY","chrY","chrY"),
                           ranges=IRanges::IRanges(start=c(1,2699521,155260561,1,2649521,59363567),
                                                   end=c(60000,154931043,155270560,10000,59034049,59373566)),
                           ploidy=c(rep("XnonPAR",3),rep("YnonPAR",3))),
              hg38=GRanges(seqnames=c("chrX","chrX","chrX","chrY","chrY","chrY"),
                           ranges=IRanges::IRanges(start=c(1,2781480,156030896,1,2781480,57217416),
                                                   end=c(10000,155701383,156040895,10000,56887902,57227415)),
                           ploidy=c(rep("XnonPAR",3),rep("YnonPAR",3))))
            genomePloidy[["GRCh37"]]=genomePloidy[["hg19"]]
            GenomeInfoDb::seqlevelsStyle(genomePloidy[["GRCh37"]])="NCBI"
            genomePloidy[["GRCh38"]]=genomePloidy[["hg38"]]
            GenomeInfoDb::seqlevelsStyle(genomePloidy[["GRCh38"]])="NCBI"
            
            # Process varSet (if provided)/ranges/or VAR_id
            if( !is.null(varSet) ) {
              if (is(varSet, "varSetList") && length(varSet) == 1) {
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
              VAR_id <- extractRanges(object, ranges = ranges)
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
              SM <- getCohort(object, cohort)
              cohort_name <- cohort
            }
            
            # Verify ploidy settings
            if (is.null(checkPloidy))
            {
              ploidy <- "diploid"
            } else {
              if (!checkPloidy %in% names(genomePloidy))
              {stop(sprintf("Unrecognized value for checkPloidy. Must be one of %s. 
                            Alternatively you can provide no value for checkPloidy and manually reset values from the default of 'diploid' using the resetPloidy function",paste(names(genomePloidy),collapse=",")))}
              
              # get varinfo and overlap with nonPAR regions
              varInfo=DBI::dbGetQuery(object,
                                      sprintf("select VAR_id, 'diploid' as ploidy, CHROM, POS, REF from var where VAR_id in ('%s');",
                                              paste(VAR_id,collapse="','")))
              varInfo$end=varInfo$POS+nchar(varInfo$REF)
              varInfo=makeGRangesFromDataFrame(varInfo, seqnames.field = "CHROM", start.field = "POS", end.field="end", keep.extra.columns = TRUE)
              GenomeInfoDb::seqlevelsStyle(varInfo)="NCBI"
              withCallingHandlers(
                overlaps <- findOverlaps(varInfo,genomePloidy[[checkPloidy]]), 
                warning = suppressSeqinfowarning)
              
              ## update ploidy for variants that overlap nonPAR regions
              for(i in 1:length(overlaps))
              {
                varInfo$ploidy[queryHits(overlaps)[i]]=genomePloidy[[checkPloidy]]$ploidy[subjectHits(overlaps)[i]]
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
            
            ## unpack genotypes and store in matrix
            GT <- as.integer(sapply(GT[,2],function(x){memDecompress(unlist(x),type="gzip",asChar=FALSE)}))-48
            GT[GT > 2] <- NA_real_
            GT <- matrix(GT,
                         nrow = nvar,
                         ncol = length(GT)/nvar,
                         byrow=TRUE)
            if(verbose) message(sprintf("Retrieved genotypes for %s variants",nvar))
            
            ## generate genoMatrix
            return(genoMatrix(GT=GT,SM=SM,VAR_id=VAR_id,w=w,ploidy=ploidy,varSetName=varSetName,unit=unit,cohortname=cohort_name,verbose=verbose))
          })

#' @rdname subsetGdb
#' @usage NULL
#' @export
setMethod("subsetGdb", signature="gdb",
          definition=function(object, output, intersection, where, tables, skipIndexes, overWrite)
          {
            # Check for existence of output gdb
            if (file.exists(output))
              {
               if (overWrite){file.remove(output)} else {stop(sprintf("output gdb already exits '%s'. Must set overWrite=TRUE to replace existing files.",output))}
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
            if (length(intersection) > 0){intersection=unlist(strsplit(intersection,split=","))}
            for (i in intersection)
            {
              query=sprintf("%s inner join %s using (VAR_id)",query,i)
            }
            if (length(where) > 0){query=sprintf("%s where %s",query, where)}
            DBI::dbExecute(object,query)

            # Copy remaining base tables
            DBI::dbExecute(object,sprintf("create table %s.SM as select * from SM", tmp))
            DBI::dbExecute(object,sprintf("create table %s.dosage as select dosage.VAR_id,dosage.GT from dosage inner join %s.var using (VAR_id)", tmp, tmp))
            DBI::dbExecute(object,sprintf("create table %s.cohort as select * from cohort where name in (%s)", tmp, paste(paste0("'", tables.cohort,"'"),collapse=",")))
            DBI::dbExecute(object,sprintf("create table %s.anno as select * from anno where name in (%s)", tmp, paste(paste0("'", tables.anno,"'"),collapse=",")))
            DBI::dbExecute(object,sprintf("create table %s.meta as select * from meta", tmp))

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
            addRangedVarinfo(gdb(output), overwrite=TRUE)

            message("Output gdb creation complete")
          })

#' writeVcf
#' @export
setMethod("writeVcf", signature="gdb",
          definition=function(object, output, VAR_id, IID, includeGeno)
            {

            # Open output connection
            tryCatch({output=gzfile(output,"w")}, error=function(e){stop(sprintf("Could not write to output path '%s'",output))})

            # Establish sample filtering rule
            SM=DBI::dbGetQuery(object, "select * from SM")
            if (missing(IID))
              {
              warning("WARNING: No IID provided, all samples will be included in output vcf")
              IID=rep(TRUE,nrow(SM))
              } else
              {
              IID=SM$IID %in% IID
              }
            message(sprintf("%s/%s samples to be retained in output", sum(IID), nrow(SM)))

            # Variant retrieval queries with / without genotype data
            if (includeGeno)
              {
              recordQuery="select CHROM, POS, VAR_id, REF, ALT, QUAL, FILTER, '.', '.', GT from var inner join dosage using (VAR_id)"
              } else
              {
              recordQuery="select CHROM, POS, VAR_id, REF, ALT, QUAL, FILTER, '.', '.' from var"
              }

            # Send variant retrieval queries with / without VAR_id filtering
            if (missing(VAR_id))
              {
              warning("WARNING: No VAR_id provided, all variants will be included in output vcf")
              records=DBI::dbSendQuery(object, recordQuery)
              } else
              {
                message(sprtinf("%s variants to be retained in output",length(VAR_id)))
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

setMethod("uploadAnno", signature="gdb",
          definition=function(object,name,value,sep,skipRemap,skipIndexes,ignoreAlleles,mapRef)
            {
            # check validity of table name:
            ## shouldn't be a protected table / no unsupported characters / shouldn't exist as cohort table
            if (name %in% gdb_protected_tables){stop(sprintf("'%s' already exists as a protected table in gdb and cannot be replaced", name))}
            if (grepl("[\\.\\+\\-\\,]", name)) {stop("Table name shouldn't contain '.','+','-' or ','")}
            cohort=listCohort(object)
            if (name %in% cohort$name){stop(sprintf("%s already exists as a cohort table. Please drop this table before continuing", name))}
            
            # write table
            if (!is(value)[[1]] %in% c("character","data.frame")){stop("value must be a valid filename or data frame object")}
            if (is(value)[[1]]=="character")
            {
              message(sprintf("Loading table '%s' from '%s'\n",name, value))
              DBI::dbWriteTable(con=object,name=name,value=value,sep=sep,overwrite=TRUE)
            } else
              {
                message(sprintf("Loading table '%s' from interactive R session'\n",name))
                DBI::dbWriteTable(con=object,name=name,value=value,overwrite=TRUE); value="interactive_session"}
            fields=DBI::dbListFields(object,name)
            message(sprintf('%s fields detected (%s)\n',length(fields),paste(fields,collapse=",")))

            # Update annotation meta-data table
            anno=listAnno(object)
            if (name %in% anno$name){DBI::dbExecute(object,"delete from anno where name=:name",params=list(name=name))}
            DBI::dbExecute(object,"insert into anno values (:name, :value, :date)",params=list(name=name,value=value,date=date()))

            # map positions to VAR_id (if skipRemap=TRUE)
            if (!skipRemap)
              {
              if (!ignoreAlleles)
                {
                if (sum(c("CHROM","POS","REF","ALT") %in% fields)!=4){warning("Mapping not possible: CHROM, POS, REF, ALT not provided.")}
                DBI::dbExecute(object,sprintf("create index %s_idx2 on %s (CHROM, POS, REF,ALT)",name,name))
                DBI::dbExecute(object,sprintf("create table tmp as select %s.VAR_id, %s.* from %s left join %s using (CHROM,POS,REF,ALT)",mapRef,name,name,mapRef))
                } else
                  {
                  if (sum(c("CHROM","POS") %in% fields)!=2){warning("Mapping not possible: CHROM, POS not provided.")}
                  DBI::dbExecute(object,sprintf("create index %s_idx2 on %s (CHROM, POS)",name,name))
                  DBI::dbExecute(object,sprintf("create table tmp as select %s.VAR_id, %s.* from %s left join %s using (CHROM,POS)",mapRef,name,name,mapRef))
                  }
              DBI::dbExecute(object,sprintf("drop table %s",name))
              DBI::dbExecute(object,sprintf("alter table tmp rename to %s",name))
              fields=DBI::dbListFields(object,name)
              nullCount=DBI::dbExecute(object,sprintf("select count(1) as n from %s where VAR_id is null",name))[1]
              if (as.numeric(nullCount)>0){message(sprintf("WARNING: %s rows where VAR_id assignment failed\n",nullCount))}
            }
            
            # index table
            if (!skipIndexes)
            {
              if ("VAR_id" %in% fields){DBI::dbExecute(object,sprintf("create index %s_idx1 on %s (VAR_id)",name,name))}
              if (sum(c("CHROM","POS","REF","ALT") %in% fields)==4){DBI::dbExecute(object,sprintf("create index %s_idx2 on %s (CHROM, POS, REF,ALT)",name,name))}
              if (ignoreAlleles & sum(c("CHROM","POS") %in% fields)==2){DBI::dbExecute(object,sprintf("create index %s_idx2 on %s (CHROM, POS)",name,name))}
            }
            })


setMethod("uploadCohort", signature="gdb",
          definition = function(object,name,value,sep="\t")
            {
            # check validity of table name:
            ## shouldn't be a protected table / no unsupported characters / shouldn't exist as anno table
            if (name %in% gdb_protected_tables){stop(sprintf("'%s' already exists as a protected table in gdb and cannot be replaced", name))}
            if (grepl("[\\.\\+\\-\\,]", name)) {stop("Table name shouldn't contain '.','+','-' or ','")}
            anno=listAnno(object)
            if (name %in% anno$name){stop(sprintf("%s already exists as a variant annotation table. Please use drop this table before continuing", name))}

            # upload, mapping and sanity checks
            if (!is(value)[[1]] %in% c("character","data.frame")){stop("value must be a valid filename or data frame object")}
            if (is(value)[[1]]=="character")
              {
              message(sprintf("Loading cohort '%s' from '%s'\n",name, value))
              upload=read.table(value,sep=sep,as.is=TRUE,h=TRUE)
              } else
                {
                  message(sprintf("Loading cohort '%s' from interactive R session\n",name))
                  upload=value
                  value="interactive_session"
                }
            
            ## check and format IID,sex columns
            message(sprintf('%s fields detected (%s)\n', ncol(upload), paste(colnames(upload),collapse=",")))
            if (!("IID" %in% colnames(upload))){stop("No 'IID' field detected, aborting upload.")}
            if (!("sex" %in% colnames(upload))){stop("No 'sex' field detected, aborting upload.")}
            upload$sex[!(upload$sex %in% c(1,2))]=0
            message(sprintf("%s males, %s females and %s unknown gender",
                            sum(upload$sex==1,na.rm=TRUE),
                            sum(upload$sex==2,na.rm=TRUE),
                            sum(upload$sex==0,na.rm=TRUE)))
            
            ## upload
            SM=DBI::dbGetQuery(object,'select * from SM')
            SM$IID=as.character(SM$IID)
            m=nrow(upload)
            upload=upload[match(SM$IID,upload$IID),]
            if (is.null(nrow(upload)) | is.null(ncol(upload))){stop("Table filtering error")}
            message(sprintf("Retained %s of %s uploaded samples that could be mapped to dosage matrix",m,nrow(upload)))
            DBI::dbWriteTable(con=object,name=name,value=upload,overwrite=TRUE)

            # Update cohort meta-data table
            cohort=DBI::dbGetQuery(object,"select * from cohort")
            if (name %in% cohort$name){DBI::dbExecute(object,"delete from cohort where name=:name",params=list(name=name))}
            DBI::dbExecute(object,"insert into cohort values (:name, :value, :date)",params=list(name=name,value=value,date=date()))
            })

setMethod("dropTable", signature="gdb",
          definition=function(object, name)
          {
            anno=listAnno(object)
            cohort=listCohort(object)
            anno=anno[anno$name!=name,]
            cohort=cohort[cohort$name!=name,]
            DBI::dbWriteTable(con=object,name="anno",value=anno,overwrite=TRUE)
            DBI::dbWriteTable(con=object,name="cohort",value=cohort,overwrite=TRUE)
            DBI::dbExecute(object,sprintf("drop table %s",name))
            message(sprintf("Table '%s' removed from gdb",name))
          }
        )

## misc ------------------------------------------------------------------------
suppressSeqinfowarning <- function(w) if( any( grepl( "sequence levels", w) ) ) invokeRestart( "muffleWarning" )
