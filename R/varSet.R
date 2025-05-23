#==============================================================================
# varSet objects
#==============================================================================

# general methods -------------------------------------------------------------

## show 

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("show", signature="varSet",
          definition=function(object)
            {
            cat(sprintf("unit=%s\nvarSetName=%s\nVAR_id=%s\nw=%s\n",
                            object@unit,
                            object@varSetName,
                            object@VAR_id,
                            object@w))
            })

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("metadata", signature="varSet",
          definition=function(x)
          {
            if (!is.null(x@metadata)) x@metadata else NULL
          })

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("getGdbId", signature="varSet",
          definition=function(object){
            metadata(object)$gdbId
          }
)

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("getRvatVersion", signature="varSet",
          definition=function(object){
            metadata(object)$rvatVersion
          }
)

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("show", signature="varSetList",
          definition=function(object)
            {
            cat(sprintf("varSetList\nContains %s records\n",
                        length(object)))
            print(head(object@varSets))
          })

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("show", signature="varSetFile",
          definition=function(object){
            cat(sprintf("varSetFile object\nPath:%s\nUnits:%s\n", object@path, length(object@units)))
          })

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("metadata", signature="varSetFile",
          definition=function(x)
          {
            x@metadata
          })

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("getGdbId", signature="varSetFile",
          definition=function(object){
            metadata(object)$gdbId
          }
)

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("getRvatVersion", signature="varSetFile",
          definition=function(object){
            metadata(object)$rvatVersion
          }
)

## length
### varSetList

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("length", signature = "varSetList",
          definition = function(x) {
            length(x@varSets)
          })

### varSetFile
#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("length", signature = "varSetFile",
          definition = function(x) {
            length(listUnits(x))
          })

## subsetting

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("[[", c("varSetList", "ANY", "missing"),
          function(x, i, j, ...)
          {
            y <- x@varSets[[i, ...]]
            y@metadata <- metadata(x)
            y
          })

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("[", c("varSetList", "ANY", "ANY"),
          function(x, i, j, ..., drop=TRUE)
          {
            if (1L != length(drop) || (!missing(drop) && drop))
              warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
            
            if (missing(i) && missing(j))
              return(x)
            
            #x@data <- x@data[i,j,...,drop=FALSE]
            initialize(x, varSets = x@varSets[i,...], units = x@units[i,...])
          })

## listing items

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("listUnits", signature="varSetFile",
          definition=function(object)
          {
            object@units
          })

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("listUnits", signature="varSetList",
          definition=function(object)
          {
            object@units
          })

#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("listVarSets", signature="varSetFile",
          definition=function(object, memlimit=5000)
          {
            # metadata
            header <- readLines(object@path, n = length(metadata_varsets) + 10)
            skip <- sum(startsWith(header, "#"))
            if ( skip > length(metadata_varsets) + 1) { # filetype + metadata
              stop ("File contains more header lines than expected.")
            }
            
            con <- gzfile(object@path,"r")
            if(skip > 0) skip <- readLines(con, n = skip)
            varsets <- c()
            while (length(i <- readLines(con,n=memlimit)) > 0)
            {
              varsets=c(varsets,
                        sapply(strsplit(i,split="\\|"),"[[",4))
            }
            close(con)
            if (length(varsets) != length(object)) {stop("Something's wrong..")}
            varsets
          })

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("listVarSets", signature="varSetList",
          definition=function(object)
          {
            unlist(lapply(object@varSets, FUN  = function(x) x@varSetName))
          })

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("listVars", signature="varSet",
          definition=function(object)
          {
            unlist(strsplit(object@VAR_id,split=","))
          })

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("listWeights", signature="varSet",
          definition=function(object)
          {
           as.numeric(unlist(strsplit(object@w,split=",")))
          })

## extract metadata
#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("getGdbId", signature="varSetList",
          definition=function(object){
            metadata(object)$gdbId
          }
)

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("getRvatVersion", signature="varSetList",
          definition=function(object){
            metadata(object)$rvatVersion
          }
)


# Constructors -------------------------------------------------------------------

# Initialize a varSet object
varSet=function(unit,varSetName,VAR_id,w)
{
  new("varSet",unit=unit,varSetName=varSetName,VAR_id=VAR_id,w=w)
}


#' varSetList
#' @rdname varSetList
#' @usage NULL
#' @export
varSetList <- function(x, metadata = list())
{
  
  if(missing(x)) {
    return(new("varSetList", varSets = list(), units = character(), metadata = list()))
  }
  
  # If list of varSet Records
  if (is(x, "list")) {
    units <- unlist(lapply(x, function(x) {x@unit}))
    return(new("varSetList", varSets = x, units = units, metadata = metadata))
  }
  
  # If data frame
  if (is(x, "data.frame"))
  {
    if (mean(colnames(x) %in% c("unit", "varSetName", "VAR_id", "w"))<1)
    {stop("Invalid input. Input must be either an object of class list that contains one or more varSet records, or else a data frame with fields unit, varSetName, VAR_id, w")}
    out=vector(mode="list",length=nrow(x))
    for (i in 1:nrow(x))
    {
      out[[i]]=varSet(unit=x$unit[i], varSet=x$varSetName[i], VAR_id=x$VAR_id[i], w=x$w[i])
    }
    return(new("varSetList",varSets = out, units = x$unit, metadata = metadata))
  }
  
  # Invalid input
  {stop("Invalid input. Input must be either an object of class list that contains one or more varSet records, or else a data frame with fields unit, varSetName, VAR_id, w")}
}


#' @importFrom S4Vectors setValidity2
S4Vectors::setValidity2("varSetList", function(object)
{
  msg=NULL

  # Ensure all records are of class varSet
  if(length(object) > 0) {
    for (i in 1:length(object@varSets))
    {
      if (!is(object@varSets[[i]], "varSet"))
      {
        msg=c(msg,
              "Invalid varSetList. Please supply only a list of varSet records or a valid data frame when running varSetList construction")
      }
    }
  }
  
  # Return
  if (length(msg)){msg} else{TRUE}
})


#' @rdname varSetFile
#' @usage NULL
#' @export
varSetFile=function(path, memlimit = 5000)
{
  # read in metadata
  metadata <- .parse_rvat_header(path, 
                                 expected_metadata = metadata_varsets,
                                 expected_filetype = "varSetFile",
                                 n = length(metadata_varsets) + 1 # file description + metadata
  )
  
  header <- readLines(path, n = length(metadata_varsets) + 10)
  skip <- sum(startsWith(header, "#"))
  
  # read units
  con=gzfile(path,"r")
  units=c()
  counter <- 1
  while (length(i <- readLines(con,n=memlimit)) > 0)
  {
    if (counter == 1) {
      i <- i[(skip + 1):length(i)]
    }
    units=c(units,
            sapply(strsplit(i,split="\\|"),"[[",1))
    counter <- counter + 1
  }
  close(con)
  new("varSetFile", path=path, units=units, metadata=metadata)
}



# Getters ----------------------------------------------------------------------
#' @rdname varSetFile
#' @usage NULL
#' @export
setMethod("getVarSet", signature="varSetFile",
          definition=function(object, unit, varSetName = NULL)
          {
            if (!all(unit %in% object@units)) {
              warning("Not all specified units are present in varSetFile, use `listUnits()` to check which units are available.")
            }
            
            # metadata
            header <- readLines(object@path, n = length(metadata_varsets) + 10)
            skip <- sum(startsWith(header, "#"))
            if ( skip > length(metadata_varsets) + 1 ) { # filetype + metadata
              stop ("File contains more header lines than expected.")
            }
            
            indices <- sort(which(object@units %in% unit))
            if( length(indices) > 1 ) indices[2:length(indices)] <- (dplyr::lead(indices) - indices)[1:(length(indices) - 1)]
            
            con <- gzfile(object@path,"r")
            if(skip > 0) skip <- readLines(con, n = skip)
            x <- lapply(indices,
                        FUN = function(i) {
                          i <- scan(con, skip = i - 1, nlines = 1, what = "character", quiet = TRUE)
                          i <- unlist(strsplit(i, split = "\\|"))
                          varSet(unit = i[1], varSet = i[4], VAR_id = i[2], w = i[3])
                        })
            x <- varSetList(x)
            close(con)
            
            # check 
            if ( !all(listUnits(x) %in% unit))  {
              stop ("Something's wrong..")
            }
            
            # additionally subset varSets, if specified 
            if( !is.null(varSetName) ) {
              x <- getVarSet(x, varSetName = varSetName)
            }
            x@metadata <- metadata(object)
            return(x)
          })

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("getVarSet", signature="varSetList",
          definition=function(object, unit = NULL, varSetName = NULL)
          {
            if(is.null(unit) && is.null(varSetName)) message("At least one of `unit` or `varSetName` should be specified.")

            if(!is.null(unit)) {
              if (!all(unit %in% listUnits(object))) {
                warning("Not all specified units are present in the varSetList")
              }
              object <- object[listUnits(object) %in% unit]
            }
            if(!is.null(varSetName)) {
              varsets <- listVarSets(object)
              if (!all(varSetName %in% varsets)) {
                warning("Not all specified varSets are present in the varSetList")
              }
              object <- object[varsets %in% varSetName]
            }
            return(object)
          })


# write ------------------------------------------------------------------------
#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("write", "varSetList",
          function(x, file = "data", append = FALSE)
          {
            out <- gzfile(file,"w")
            metadata <- metadata(x)
            metadata$rvatVersion <- as.character(packageVersion("rvat"))
            metadata$creationDate <- as.character(round(Sys.time(), units = "secs"))
            .write_rvat_header(filetype = "varSetFile", 
                               metadata = metadata, 
                               con = out)
            write.table(as.data.frame(x), out, sep="|", append = append, row.names = FALSE, col.names = FALSE, quote = FALSE)
            close(out)
          })


# varset utils -----------------------------------------------------------------

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("as.data.frame", signature = "varSetList",
          definition = function(x) {
            do.call(rbind, lapply(x@varSets, FUN = as.data.frame))
          })

#' @rdname varSet
#' @usage NULL
#' @export
setMethod("as.data.frame", signature = "varSet",
          definition = function(x) {
            data.frame(unit = x@unit,
                       VAR_id = x@VAR_id,
                       w = x@w,
                       varSetName = x@varSetName,
                       stringsAsFactors = FALSE
            )
          })

#' buildVarSet-gdb
#' 
#' Generate optionally weighted variant sets using annotation table(s) uploaded to the gdb.
#' See the tutorials for examples.
#' 
#' @rdname buildVarSet-gdb
#' @name buildVarSet-gdb
#' @aliases buildVarSet,gdb-method
#' @param object a [`gdb`] object.
#' @param varSetName Name to assign varSet grouping. This identifier column is used to allow for subsequent mergeing of multiple varSet files for coordinated analysis of multiple variant filtering/ weighting strategies)
#' @param unitTable Table containing aggregation unit mappings.
#' @param unitName Field to utilize for aggregation unit names.
#' @param output Output file name (output will be gz compressed text).
#' @param intersection Additional tables to filter through intersection (i.e. variants absent from intersection tables will not appear in output). Multiple tables should be ',' delimited.
#' @param where An SQL compliant where clause to filter output; eg: "CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')".
#' @param weightName Field name for desired variant weighting, must be a column within unitTable or other intersection table. Default value of 1 is equivalent to no weighting.
#' @param verbose Should the function be verbose? Defaults to `TRUE`.
#' 
#' @examples
#'
#' library(rvatData)
#'
#' # Build a varset including variants with a moderate predicted impact
#' gdb <- create_example_gdb()
#' varsetfile_moderate <- tempfile()
#' buildVarSet(object = gdb, 
#'             output = varsetfile_moderate,
#'             varSetName = "Moderate", 
#'             unitTable = "varInfo", 
#'             unitName = "gene_name",
#'             where = "ModerateImpact = 1")
#' 
#' # Build a varset that contains CADD scores
#' varsetfile_cadd <- tempfile()
#' buildVarSet(object = gdb, 
#'             output = varsetfile_cadd,
#'             varSetName = "CADD", 
#'             unitTable = "varInfo", 
#'             unitName = "gene_name",
#'             weightName = "CADDphred")
#'             
#' # connect to varsetfile and retrieve variant sets
#' varsetfile <- varSetFile(varsetfile_moderate)
#' varsets <- getVarSet(varsetfile, unit = c("SOD1", "FUS"))
#' 
#' # see ?getVarSet, ?varSetFile and ?varSetList for more details on connecting and handling varsetfiles.
#' # see e.g., ?assocTest and ?aggregate for downstream methods that can loop through varsetfiles and varsetlists.
#'
#' @export
setMethod("buildVarSet", 
          signature="gdb",
          definition = function(object,
                                varSetName,
                                unitTable,
                                unitName,
                                output = NULL,
                                intersection = NULL,
                                where = NULL,
                                weightName = "1",
                                verbose = TRUE
                                )
          {

            # Build query
            query=sprintf("select distinct %s as unit, VAR_id, %s as w from %s", unitName, weightName, unitTable)
            if (!is.null(intersection)){intersection=unlist(strsplit(intersection,split=","))}
            for (i in intersection)
            {
              query=sprintf("%s inner join %s using (VAR_id)",query,i)
            }
            if (!is.null(where)){query=sprintf("%s where %s",query, where)}
            query=sprintf("select unit, group_concat(VAR_id) as VAR_id, group_concat(w) as w, '%s' as varSetName from (%s) x group by unit",varSetName, query)

            # Output to varSetFile
            metadata <- list(
              rvatVersion = as.character(packageVersion("rvat")),
              gdbId = getGdbId(object),
              genomeBuild = getGenomeBuild(object),
              creationDate = as.character(round(Sys.time(), units = "secs"))
            )
            
            if (!is.null(output))
            {
              handle <- DBI::dbSendQuery(object,query)
              outhandle=gzfile(output,"w")
              
              # write metadata
              .write_rvat_header(filetype = "varSetFile", 
                                 metadata = metadata, 
                                 con = outhandle)
              
              while (!DBI::dbHasCompleted(handle))
              {
                i=DBI::dbFetch(handle,n=1)
                write.table(i, outhandle, append=TRUE, sep="|", row.names=FALSE, col.names=FALSE, quote=FALSE)
              }
              DBI::dbClearResult(handle)
              close(outhandle)
              if (verbose) message(sprintf("Generated varSetFile: %s",output))
              return(varSetFile(output))
            } else {
              # Output to varSetList
              varsetlist <- varSetList(DBI::dbGetQuery(object, query))
              varsetlist@metadata <- metadata
              return(varsetlist)
            }
          })


#' buildVarSet-data.frame
#' 
#' Generate optionally weighted variant sets using annotation table(s).
#' See the tutorials for examples.
#' 
#' @rdname buildVarSet-data.frame
#' @name buildVarSet-data.frame
#' @aliases buildVarSet,data.frame-method
#' @param object A data.frame containing variant annotations
#' @param unitName Field to utilize for aggregation unit names.
#' @param fields Which fields in the data.frame to use as variant annotations.
#' These fields can be 1) an indicator (0,1 or FALSE/TRUE) that flags variants with the annotation, 
#' e.g. a column named 'LOF' that indicates whether the variant is predicted to lead to loss-of-function.;
#' or 2) variant weights, e.g. a column named 'CADD' that contains CADD scores.
#' Multiple fields can be specified, which will result in multiple rows per aggregation unit in the resulting varSetFile.
#' @param output Optional output file name (output will be gz compressed text).
#' @param memlimit Defaults to `NULL`, currently not implemented.
#' @param mask If the annotation field is a 0,1 or FALSE/TRUE indicator, should variants with 0/FALSE be dropped? Defaults to `TRUE`.
#' If `FALSE`, all variants will be included, and will be assigned weights 0 and 1.
#'
#' @examples
#' 
#' library(rvatData)
#' 
#' gdb <- create_example_gdb()
#' anno <- getAnno(gdb, "varinfo", where = "gene_name in ('SOD1', 'FUS')")
#' varsetfile_from_df <- tempfile()
#' buildVarSet(
#'   anno,
#'   unitName = "gene_name",
#'   fields = c("HighImpact"),
#'   output = varsetfile_from_df
#' )
#' 
#' # connect to varsetfile and retrieve variant sets
#' varsetfile <- varSetFile(varsetfile_from_df)
#' varsets <- getVarSet(varsetfile, unit = c("SOD1", "FUS"))
#' 
#' # see ?getVarSet, ?varSetFile and ?varSetList for more details on connecting and handling varsetfiles.
#' # see e.g., ?assocTest and ?aggregate for downstream methods that can loop through varsetfiles and varsetlists.
#'
#' @export
setMethod("buildVarSet", 
          signature="data.frame",
          definition = function(object, 
                                unitName, 
                                fields, 
                                output = NULL, 
                                memlimit = NULL, 
                                mask = TRUE)
          {
            ## input checks
            if (!"VAR_id" %in% colnames(object)) {
              stop("'VAR_id' field should be present")
            }
            
            if (length(unitName) > 1) {
              stop("Specify 1 unitName")
            }
            
            if(!unitName %in% colnames(object)) {
              stop(sprintf("%s field should be present", unitName))
            }
            
            if( !all(fields %in% colnames(object))) {
              stop(sprintf("The following fields are not present: %s", paste(fields[!fields %in% colnames(object)], collapse = ",")))
            }
            
            ## subset relevant columns
            object <- object[,c("VAR_id", unitName, fields)]
            colnames(object)[2] <- "unit"
            
            ## test for each field whether it's 0,1 or weights 
            if(is.logical(mask) && mask) {
              is_mask <- vector("logical", length(fields))
              names(is_mask) <- fields
              for (field in fields) {
                if (all (object[[field]][!is.na(object[[field]])] %in% c(0,1))) {
                  is_mask[field] <- TRUE
                } else {
                  is_mask[field] <- FALSE
                }
              }
              is_mask <- data.frame(varSetName = fields, is_mask = unname(is_mask), stringsAsFactors = FALSE)
            } else if (is.logical(mask) && !mask) {
              is_mask <- data.frame(varSetName = fields, is_mask = FALSE, stringsAsFactors = FALSE)
            } else if (is.character(mask)) {
              is_mask <- data.frame(varSetName = fields, is_mask = fields %in% mask, stringsAsFactors = FALSE)
            } else {
              stop("`mask` should be either `TRUE`, `FALSE` or a character vector of fields that should be treated as mask.")
            }
            
            ## generate varset
            varSet <- object %>%
              tidyr::pivot_longer(
                cols = dplyr::all_of(fields),
                names_to = "varSetName",
                values_to = "value"
              ) %>%
              dplyr::arrange(unit) %>%
              dplyr::left_join(
                is_mask, by="varSetName"
              ) %>%
              dplyr::filter(
                (is_mask & !is.na(value) & value == 1) |
                  (!is_mask & !is.na(value))
              ) %>%
              dplyr::group_by(dplyr::across(c("unit", "varSetName"))) %>%
              dplyr::summarize(
                VAR_id = paste(VAR_id, collapse=","),
                w = paste(value, collapse = ","), .groups = "drop"
              ) %>%
              dplyr::select(unit, VAR_id, w, varSetName)
            
            ## save 
            metadata <- list(
              rvatVersion = as.character(packageVersion("rvat")),
              gdbId = NA_character_,
              genomeBuild = NA_character_,
              creationDate = as.character(round(Sys.time(), units = "secs"))
            )
            
            ## save 
            if (is.null(output)) {
              return(varSetList(varSet, metadata = metadata))
            } else {
              write(varSetList(varSet, metadata = metadata), file = output) 
            }
          }
)


## build varSetlist from VAR_ids with 'memlimit' option
.varsTovarSetList <- function(VAR_id, chunkSize = Inf) {
  
  chunks <- split(VAR_id, ceiling(seq_along(1:length(VAR_id))/chunkSize))
  varset <- varSetList(lapply(1:length(chunks),
         function(x) {
           varSet =  varSet(unit=paste0("chunk",x), 
                            varSet="none", 
                            VAR_id=paste(chunks[[x]], collapse=","), 
                            w=paste(rep("1",length(chunks[[x]])),collapse=","))
         }))
  varset
}

#' collapseVarSetList
#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("collapseVarSetList", signature="varSetList",
          definition=function(object, unit = "unnamed", varSetName = "unnamed", drop = TRUE)
          {
            VAR_id <- unique(unlist(lapply(object@varSets, listVars)))
            varset <- varSet(unit=unit, varSet=varSetName, VAR_id=paste(VAR_id, collapse=","), w=paste(rep("1",length(VAR_id)),collapse=","))
            if(drop) {
              varset
            } else {
              varSetList(list(varset))
            }
          })

#' @rdname varSetList
#' @usage NULL
#' @export
setMethod("metadata", signature="varSetList",
          definition=function(x)
          {
            x@metadata
          })

#' @rdname getRanges
#' @usage NULL
#' @export
setMethod(
  "getRanges",
  signature = "varSetFile",
  definition = function(object,
                        gdb,
                        output = NULL,
                        table = "var",
                        CHROM = "CHROM",
                        POS = "POS",
                        where = c()) {
    # read in varset data
    data <- readr::read_delim(
      object@path,
      delim = "|",
      col_names = c("unit", "VAR_id", "w", "varSetName"),
      col_types = "cccc",
      comment = "#"
    )
    
    # get ranges 
    .get_ranges(
      data = data,
      gdb = gdb,
      output = output,
      table = table,
      CHROM = CHROM,
      POS = POS,
      where = where
    )
    
  }
)

#' @rdname getRanges
#' @usage NULL
#' @export
setMethod(
  "getRanges",
  signature = "varSetList",
  definition = function(object,
                        gdb,
                        output = NULL,
                        table = "var",
                        CHROM = "CHROM",
                        POS = "POS",
                        where = c()) {
    # convert varsetlist to data.frame
    data <- as.data.frame(object)
    
    # get ranges 
    .get_ranges(
      data = data,
      gdb = gdb,
      output = output,
      table = table,
      CHROM = CHROM,
      POS = POS,
      where = where
    )
  }
)

.get_ranges <- function(data,
                        gdb,
                        output = NULL,
                        table = "var",
                        CHROM = "CHROM",
                        POS = "POS",
                        where = c()) {
  
  # check where statement, can be either of length 1 or same length as varSetList/varSetFile
  if (length(where) > 0) {
    if (length(where) > 1 & length(where) != nrow(data)) {
      stop(
        "The `where` parameter should be either of length 1, or of the same length as the provided varSetList/varSetFile"
      )
    }
  }
  
  # extract ranges for each record
  start <- vector(mode = "integer", length = nrow(data))
  end <- vector(mode = "integer", length = nrow(data))
  chrom <- vector(mode = "character", length = nrow(data))
  not_found <- 0L
  
  # check if required fields are available in specified table
  check_table <- DBI::dbGetQuery(gdb, sprintf("select * from %s limit 1", table))
  required_fields <- c("VAR_id", CHROM, POS)
  if (!all(required_fields %in% colnames(check_table)) ) {
    stop(sprintf("The following required fields are not present in table '%s': %s", 
                 table,
                 paste(required_fields[!required_fields %in% colnames(check_table)], collapse = ",")
    ))
  }
  
  for (i in 1:nrow(data)) {
    
    # extract var ids
    dat_ <- data[i, ]
    varids <- unlist(stringr::str_split(dat_$VAR_id, pattern = ","))
    
    # retrieve ranges
    pos <- getAnno(
      gdb,
      table = table,
      VAR_id = varids,
      where = if (length(where) > 1)
        where[i]
      else
        where
    )
    
    # raise error if multiple records are found for a single VAR_id
    if (sum(duplicated(pos$VAR_id)) > 0) {
      stop(
        sprintf(
          "More than one record per variant found in table '%s', first encountered at unit '%s'",
          table,
          dat_$unit
        )
      )
      
      # else retrieve ranges
    } else if (nrow(pos) > 0) {
      start[i] <- min(pos[[POS]])
      end[i] <- max(pos[[POS]])
      chrom[i] <- unique(pos[[CHROM]])
    } else {
      start[i] <- NA_integer_
      end[i] <- NA_integer_
      chrom[i] <- NA_character_
      not_found <- not_found + 1L
    }
  }
  
  # parse
  ranges <- data.frame(
    CHROM = chrom,
    start = start,
    end = end,
    unit = data$unit,
    varSetName = data$varSetName,
    stringsAsFactors = FALSE
  )
  
  # notify if no positions were found for records
  if (not_found > 0) {
    warning(
      sprintf(
        "Note: for %s/%s records in the varSetFile, no position mappings where found for any of the variants.",
        not_found,
        nrow(data)
      )
    )
  }
  
  # write output, if specified
  if (!is.null(output)) {
    readr::write_tsv(ranges, file = gzfile(output), col_names = TRUE)
    return(invisible(NULL))
  } else {
    ranges
  }
}