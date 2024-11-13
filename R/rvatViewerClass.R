# ===================================================================================
#  rvatViewerClass class
# ===================================================================================

setClass("rvatViewerClass", 
         representation(
           rvatResults = "list",
           gdb = "character", # placeholder, implement later
           varSet = "character", # placeholder, implement later
           metadata = "list"
         ),
         prototype = list(
           rvatResults = list(),
           gdb = character(), # placeholder, implement later
           varSet = character(), # placeholder, implement later
           metadata = list()
         )
)

setClass("rvatViewerClassSQL", # think of a better class name
         representation(
           rvatResults = "SQLiteConnection",
           gdb = "character", # placeholder, implement later
           varSet = "character", # placeholder, implement later
           metadata = "list"
         ),
         prototype = list(
           gdb = character(), # placeholder, implement later
           varSet = character(), # placeholder, implement later
           metadata = list()
         )
)

# -----------------------------------------------------------------------------------
# constructor

## GT, varSet will be implemented later
prepareDataShiny <- function(
  rvatResults,  
  GT = NULL, 
  varSet = NULL
) {
  
  ## Check & parse input
  rvatResultsisSQL <- FALSE
  if(is.character(rvatResults)) {
    tryCatch({rvatResults=DBI::dbConnect(DBI::dbDriver("SQLite"),rvatResults)}, error=function(e){stop(sprintf("Invalid db path '%s'",rvatResults))})
    rvatResultsisSQL <- TRUE
  } else if (is.data.frame(rvatResults)) {
    if(all(names(columns_rvbResults) %in% colnames(rvatResults))) {
      rvatResults <- list(analysis1 = rvbResult(rvatResults))
    } else if(all(names(columns_singlevarResults) %in% colnames(rvatResults))) {
      rvatResults <- list(analysis1 = singlevarResult(rvatResults %>% dplyr::mutate(unit=VAR_id)))
    } else {
      stop("data.frame is not convertible to object of class 'singlevarResult' or 'rvbResult'")
    }
  } else if(is.list(rvatResults)) {
    if(is.null(names(rvatResults)) || sum(names(rvatResults) == "") > 0) {
      stop("the list should be named")
    }
    check_data.frame <- unlist(lapply(rvatResults, FUN = is.data.frame))
    check_rvatResults <- unlist(lapply(rvatResults, FUN = function(x) is(x, "rvatResult")))
    if(all(check_data.frame)) {
      names <- names(rvatResults)
      rvatResults <- lapply(rvatResults, FUN = function(x) {
          if(all(names(columns_rvbResults) %in% colnames(x))) {
            return(rvbResult(x))
          } else if(all(names(columns_singlevarResults) %in% colnames(x))) {
            return(singlevarResult(x %>% dplyr::mutate(unit=VAR_id)))
          } else {
            stop("Not all data.frames are convertible to objects of class 'singlevarResult' or 'rvbResult'")
          }
        })
      names(rvatResults) <- names
    } else if (all(check_rvatResults)) {
      names <- names(rvatResults)
      rvatResults <- lapply(rvatResults, FUN = function(x) {
        if(is(x, "singlevarResult")) {
          x$unit <- x$VAR_id
        }
        x
      })
      names(rvatResults) <- names
    } else {
      stop("the list should contain objects of either class `data.frame`,`rvbResult` or `singlevarResult`")
    }
  } else if(is(rvatResults, "rvatResult")) {
    if(is(rvatResults, "singlevarResult")) rvatResults$unit <- rvatResults$VAR_id
    rvatResults <- list(analysis1 = rvatResults)
  } else {
    stop("The `rvatResults` parameter should be one of the following: `data.frame`, `list`, `rvbResult`, `singlevarResult` or `character`")
  }
  
  ## Extract metadata
  if(rvatResultsisSQL) {
    metadata <- extract_metadata_rvatResults_sql(rvatResults)
  } else {
    metadata <- lapply(rvatResults, FUN = function(x) {
      meta <- list(
        columns = colnames(x),
        column_types = unlist(lapply(x, FUN = class)))
      for(var in c("cohort", "pheno", "covar", "test", "name", "varSetName", "geneticModel", "MAFweight")) {
        meta[[var]] <- unique(x[[var]])
      }
      meta
    })
    metadata$analyses <- names(rvatResults)
  }
  
  if(!rvatResultsisSQL) {
    new("rvatViewerClass",
        rvatResults = rvatResults,
        metadata = metadata
    )
  } else {
    new("rvatViewerClassSQL",
        rvatResults = rvatResults,
        metadata = metadata
    )
  }
}

extract_metadata_rvatResults_sql <- function(conn) {
  tables <- DBI::dbListTables(conn)
  metadata <- list() 
  
  for(table in tables) {
    tab <- conn %>% dplyr::tbl(table)
    fields <- tab %>% colnames()
    
    if(!all(names(columns_rvbResults) %in% fields) && !all(names(columns_singlevarResults) %in% fields)) {
      stop(sprintf("Table %s does not contain all required fields", table))
    }
    
    ## Extract column names and types 
    columns <- names(columns_rvbResults)[names(columns_rvbResults) %in% fields]
    types <- unlist(lapply(tab %>% head(5) %>% data.frame(), class))
    type_check <- unlist(lapply(columns, FUN = function(x) types[x] %in% columns_rvbResults[[x]]))
    names(type_check) <- columns
    if(!all(type_check)) {
      type_msg <- lapply(
        names(type_check)[!type_check], FUN = function(x) sprintf("%s: %s", x, paste(columns_rvbResults[[x]], collapse=" or "))
      ) %>% paste(collapse="\n")
    stop(sprintf("The following columns should be of type:\n%s", type_msg))
    }
    ## Check integers, booleans are written as 0,1 in database
    ## If type=="integer" and all values are either zero or one, we can treat the 
    ## column as a boolean. otherwise treat as numeric. 
    if(any(types == "integer")) {
      for(col in names(types)[types=="integer"]){
        if(all(tab %>% dplyr::pull(col) %in% c(0,1))) {
          types[col] <- "logical"
        } else {
          types[col] <- "numeric"
        }
      }
    }
    
    tmp <- list()
    for(var in c("cohort", "pheno", "covar", "test", "name", "varSetName", "MAFweight")) {
      tmp[[var]] <- tab %>% dplyr::pull(var) %>% unique() 
    }
    tmp[["columns"]] <- fields
    tmp[["column_types"]] <- types
    names(tmp[["column_types"]]) <- fields
    tmp[["Nunits"]] <- tab %>% dplyr::pull("unit") %>% dplyr::n_distinct() 
    tmp[["ctrlNmin"]] <- tab %>% dplyr::pull("ctrlN") %>% min(na.rm=TRUE) 
    tmp[["ctrlNmax"]] <- tab %>% dplyr::pull("ctrlN") %>% max(na.rm=TRUE)
    tmp[["caseNmin"]] <- tab %>% dplyr::pull("caseN") %>% min(na.rm=TRUE) 
    tmp[["caseNmax"]] <- tab %>% dplyr::pull("caseN") %>% max(na.rm=TRUE) 
    metadata[[table]] <- tmp
  }
  metadata$analyses <- tables
  metadata
}

# -----------------------------------------------------------------------------------
# convenience methods

#' @rdname rvatViewer
#' @usage NULL
setMethod("metadata", c("rvatViewerClass"),
          function(x)
          {
            x@metadata
          })

#' @rdname rvatViewer
#' @usage NULL
setMethod("metadata", c("rvatViewerClassSQL"),
          function(x)
          {
            x@metadata
          })


# -----------------------------------------------------------------------------------
# subset

setGeneric("subSet", function(x, analysis, filters) standardGeneric("subSet"))

setMethod("subSet", c("rvatViewerClass"),
          function(x, analysis, filters)
          {
            result <- x@rvatResults[[analysis]]
            as.data.frame(result[eval(parse(text = filter_parser(filters, data = "result"))),])
          })

setMethod("subSet", c("rvatViewerClassSQL"),
          function(x, analysis, filters)
          {
            RSQLite::dbGetQuery(x@rvatResults, statement = filter_parser_SQL(params = filters, analysis = analysis)) 
          })

## Parsers ------------------

## data.frame
filter_parse_logical <- function(param, data) {
  if(param$keepNA) sprintf("(is.na(%1$s[['%2$s']]) | %3$s%1$s[['%2$s']])", data, param$variable, if(param$bool) "" else "!")
  if(!param$keepNA) sprintf("(!is.na(%1$s[['%2$s']]) & %3$s%1$s[['%2$s']])", data, param$variable, if(param$bool) "" else "!")
}

filter_parse_character <- function(param, data) {
  if(param$keepNA) sprintf("(is.na(%1$s[['%2$s']]) | %3$s%1$s[['%2$s']] %%in%% c(%4$s))", data, param$variable, if(param$negate) "!" else "", paste(paste0("'",param$string, "'"), collapse = ", "))
  if(!param$keepNA) sprintf("(!is.na(%1$s[['%2$s']]) & %3$s%1$s[['%2$s']] %%in%% c(%4$s))", data, param$variable, if(param$negate) "!" else "", paste(paste0("'",param$string, "'"), collapse = ", "))
}

filter_parse_numeric <- function(param, data) {
  if(param$keepNA) return(sprintf("(is.na(%1$s[['%2$s']]) | (%1$s[['%2$s']] >= %3$s & (%1$s[['%2$s']] <= %4$s)))", data, param$variable,param$min, param$max))
  if(!param$keepNA) return(sprintf("(!is.na(%1$s[['%2$s']]) & (%1$s[['%2$s']] >= %3$s & %1$s[['%2$s']] <= %4$s))", data, param$variable,param$min, param$max))
}

filter_parser <- function(params, data) {
  types <- lapply(params, function(x) x[["type"]])
  parsed_filters <- paste(c(unlist(lapply(params[types == "logical"], FUN = filter_parse_logical, data)), 
                            unlist(lapply(params[types %in% c("character", "factor", "Rle")], FUN = filter_parse_character, data)), 
                            unlist(lapply(params[types %in% c("numeric", "integer")], FUN = filter_parse_numeric, data))),
                          collapse = " & ")
  parsed_filters
}

## SQL
filter_parse_logical_SQL <- function(param) {
  if(param$keepNA) return(sprintf("(%1$s is null or %1$s = %2$s)", param$variable, if(param$bool) "1" else "0"))
  if(!param$keepNA) return(sprintf("(%1$s is not null and %1$s = %2$s)", param$variable, if(param$bool) "1" else "0"))
}

filter_parse_character_SQL <- function(param) {
  if(param$keepNA) return(sprintf("(%1$s is null or %1$s %2$s (%3$s))", param$variable, if(param$negate) "not in" else "in", paste(paste0("'",param$string, "'"), collapse = ", ")))
  if(!param$keepNA) return(sprintf("(%1$s is not null and %1$s %2$s (%3$s))", param$variable, if(param$negate) "not in" else "in", paste(paste0("'",param$string, "'"), collapse = ", ")))
}

filter_parse_numeric_SQL <- function(param) {
  if(param$keepNA) {
    sprintf("((%1$s is null) or (%1$s >= %2$s and %1$s <= %3$s))", 
     param$variable, 
     if(is.infinite(param$min) && param$min < 0) "-'infinite'" else if(is.infinite(param$min) && param$min > 0) "'infinite'" else param$min, 
     if(is.infinite(param$max) && param$max < 0) "-'infinite'" else if(is.infinite(param$max) && param$max > 0) "'infinite'" else param$max
    )} else if(!param$keepNA) {sprintf("((%1$s is not null) and (%1$s >= %2$s and %1$s <= %3$s))", 
    param$variable,
    if(is.infinite(param$min) && param$min < 0) "-'infinite'" else if(is.infinite(param$min) && param$min > 0) "'infinite'" else param$min, 
    if(is.infinite(param$max) && param$max < 0) "-'infinite'" else if(is.infinite(param$max) && param$max > 0) "'infinite'" else param$max
  )}
}

filter_parser_SQL <- function(params, analysis) {
  types <- lapply(params, function(x) x[["type"]])
  where <- paste(c(unlist(lapply(params[types %in% c("logical", "integer")], FUN = filter_parse_logical_SQL)), 
                   unlist(lapply(params[types == "character"], FUN = filter_parse_character_SQL)), 
                   unlist(lapply(params[types == "numeric"], FUN = filter_parse_numeric_SQL))),
                          collapse = " and ")
 
  sprintf("select * from %s where %s;", analysis, where)
}