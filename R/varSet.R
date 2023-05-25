#==============================================================================
# varSet objects
#==============================================================================

# general methods -------------------------------------------------------------

## show 
setMethod("show", signature="varSet",
          definition=function(object)
            {
            message(sprintf("unit=%s\nvarSetName=%s\nVAR_id=%s\nw=%s",
                            object@unit,
                            object@varSetName,
                            object@VAR_id,
                            object@w))
            })

setMethod("show", signature="varSetList",
          definition=function(object)
            {
            message("varSetList")
            message(sprintf("Contains %s records",length(object)))
            print(head(object@varSets))
          })

setMethod("show", signature="varSetFile",
         definition=function(object){
           message("varSetFile object")
           message(sprintf("Path:%s",object@path))
           message(sprintf("Units:%s",length(object@units)))
         })


## length
### varSetList
setMethod("length", signature = "varSetList",
          definition = function(x) {
            length(x@varSets)
          })
### varSetFile
setMethod("length", signature = "varSetFile",
          definition = function(x) {
            length(listUnits(x))
          })

## subsetting

setMethod("[[", c("varSetList", "ANY", "missing"),
          function(x, i, j, ...)
          {
            x@varSets[[i, ...]]
          })


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

setMethod("listUnits", signature="varSetFile",
          definition=function(object)
          {
            object@units
          })

setMethod("listUnits", signature="varSetList",
          definition=function(object)
          {
            object@units
          })

setMethod("listVarSets", signature="varSetFile",
          definition=function(object, memlimit=5000)
          {
            con=gzfile(object@path,"r")
            varsets=c()
            while (length(i <- readLines(con,n=memlimit)) > 0)
            {
              varsets=c(varsets,
                        sapply(strsplit(i,split="\\|"),"[[",4))
            }
            close(con)
            varsets
          })

setMethod("listVarSets", signature="varSetList",
          definition=function(object)
          {
            unlist(lapply(object@varSets, FUN  = function(x) x@varSetName))
          })


setMethod("listVars", signature="varSet",
          definition=function(object)
          {
            unlist(strsplit(object@VAR_id,split=","))
          })


setMethod("listWeights", signature="varSet",
          definition=function(object)
          {
           as.numeric(unlist(strsplit(object@w,split=",")))
          })


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
varSetList=function(x)
  {

  # If list of varSet Records
  if (is(x, "list")) {
    units <- unlist(lapply(x, function(x) {x@unit}))
    return(new("varSetList", varSets = x, units = units))
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
    return(new("varSetList",varSets = out, units = x$unit))
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
varSetFile=function(path,memlimit=5000)
  {
  con=gzfile(path,"r")
  units=c()
  while (length(i <- readLines(con,n=memlimit)) > 0)
    {
    units=c(units,
            sapply(strsplit(i,split="\\|"),"[[",1))
  }
  close(con)
  new("varSetFile", path=path, units=units)
  }


# Getters ----------------------------------------------------------------------
setMethod("getVarSet", signature="varSetFile",
          definition=function(object, unit, varSetName = NULL)
          {
            if (!all(unit %in% object@units)) {warning("Not all specified units are present in varSetFile, use `listUnits()` to check which units are available.")}
            indices <- sort(which(object@units %in% unit))
            unit <- object@units[indices]
            if(length(indices) > 1) indices[2:length(indices)] <- (dplyr::lead(indices)-indices)[1:(length(indices)-1)]
            con <- gzfile(object@path,"r")
            x <- lapply(indices,
                          FUN = function(i) {
                            i <- scan(con, skip = i-1, nlines = 1, what = "character", quiet = TRUE)
                            i <- unlist(strsplit(i,split="\\|"))
                            varSet(unit= i[1], varSet=i[4], VAR_id=i[2], w=i[3])
                          })
            x <- varSetList(x)
            close(con)
            
            if( !is.null(varSetName) ) {
              x <- getVarSet(x, varSetName = varSetName)
            }
            return(x)
          })


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
            write.table(as.data.frame(x), out, sep="|", append = append, row.names = FALSE, col.names = FALSE, quote = FALSE)
            close(out)
          })


# varset utils -----------------------------------------------------------------

#' @describeIn varSetList as.data.frame
#' 
#' Return the the geneSets contained in a \code{\link{varSetList}} object as a data.frame.
#' 
#' @param x a [`varSetList`] object
#' @export
setMethod("as.data.frame", signature = "varSetList",
          definition = function(x) {
            do.call(rbind, lapply(x@varSets, FUN = as.data.frame))
          })


setMethod("as.data.frame", signature = "varSet",
          definition = function(x) {
            data.frame(unit = x@unit,
                       VAR_id = x@VAR_id,
                       w = x@w,
                       varSetName = x@varSetName,
                       stringsAsFactors = FALSE
            )
          })

setMethod("buildVarSet", 
          signature="gdb",
          definition = function(object,varSetName,unitTable,unitName,output=NULL,intersection=NULL,where=NULL,weightName="1")
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
            handle=DBI::dbSendQuery(object,query)
            if (!is.null(output))
            {
              outhandle=gzfile(output,"w")
              while (!DBI::dbHasCompleted(handle))
              {
                i=DBI::dbFetch(handle,n=1)
                write.table(i,outhandle,append=TRUE,sep="|",row.names=FALSE,col.names=FALSE,quote=FALSE)
              }
              DBI::dbClearResult(handle)
              close(outhandle)
              message(sprintf("Generated varSetFile: %s",output))
              return(varSetFile(output))
            }

            # Output to varSetList
            return(varSetList(DBI::dbGetQuery(object,query)))
          })


setMethod("buildVarSet", 
          signature="data.frame",
          definition = function(object, unitName, fields, output = NULL, memlimit=NULL, mask=TRUE)
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
            object <- object[,c("VAR_id", unitName,fields)]
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
              dplyr::select(unit,VAR_id, w, varSetName)
            
            ## save 
            if (is.null(output)) {
              return(varSetList(varSet))
            } else {
              write(varSetList(varSet), file = output) 
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
