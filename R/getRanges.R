#' getRanges
#' 
#' Retrieve genomic ranges for variant sets.
#'
#' This function method retrieves genomic ranges (chromosome, start, and end positions) for variant sets defined in a `varSetFile` or `varSetList`.  
#' It will map positions based on a variant annotations table in provided `gdb` (default = "var").
#'
#' @param object Input [`varSetFile`] or [`varSetList`].
#' @param gdb A [`gdb`] object.
#' @param output Output file path. Output will be gz-compressed. Defaults to `NULL`,
#' in which case ranges are returned as a data.frame.
#' @param table Name of lookup table for extracting variant positions. Defaults to "var".
#' @param CHROM If the name of the column in `table` that specifies the chromosome is not named `CHROM`, 
#' the column name can be specified here.
#' @param POS If the name of the column in `table` that specifies variant position is not named `POS`, 
#' the column name can be specified here.
#' @param where An SQL compliant where clause to filter output; eg: "CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')". Can be either of length 1, in which case the same
#' where clause is applied for each varSet, or of length equal to the number of varSets in the provided varSetList or varSetFile,
#' in which case each where clause is applied to the corresponding varSet.
#' @examples
#' 
#' library(rvatData)
#' varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
#' gdb <- create_example_gdb()
#' ranges <- getRanges(varsetfile, gdb = gdb)
#' 
#' @export
setGeneric("getRanges", function(object,
                                 gdb,
                                 output = NULL,
                                 table = "var",
                                 CHROM = "CHROM",
                                 POS = "POS",
                                 where = c()
) standardGeneric("getRanges"))

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