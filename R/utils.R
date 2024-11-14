# metadata ---------------------------------------------------------------------
.parse_rvat_header <- function(path, expected_metadata = NULL, expected_filetype = NULL, n = 10) {
  lines <- readLines(path, n = n)
  
  # select commented lines
  header <- gsub("#\\s*", "", lines[grepl("^#", lines)])
  
  # if header is present, check filetype and parse metadata lines
  if (length(header) > 0) {
    
     # check if filetype is correctly specified
    if (!is.null(expected_filetype)) {
      filetype <- header[1]
      if (!filetype %in% sprintf("RVAT-%s", expected_filetype)) {
        stop(sprintf("Unexpected filetype: %s", filetype))
      }
    }
    
    # parse metadata
    metadata <- header[2:length(header)]
    metadata <- strsplit(metadata, ": ")
    if (!all(lengths(metadata) == 2)) {stop("Unexpected metadata")}
    metadata <- as.list(setNames(unlist(lapply(metadata, "[[", 2)), unlist(lapply(metadata, "[[", 1))))
    
  } else {
    metadata <- list()
  }
  
  # if expected metadata is specified, check if all metadata is found
  if (!is.null(expected_metadata)) {
    
    # check if other than expected metadata is included
    if (!all(names(metadata) %in% expected_metadata)) {
      stop(sprintf("The following unexpected metadata fields were found: %s", 
                   paste(names(metadata)[!names(metadata) %in% expected_metadata], collapse = ",")
      ))
    }
    
    # check if metadata is missing
    if (!all(expected_metadata %in% names(metadata))) {
      warning(sprintf("The following metadata fields are missing: %s", 
                      paste(expected_metadata[!expected_metadata %in% names(metadata)], collapse = ",")))
    }
    
    # parse
    metadata_ <- list()
    for (item in expected_metadata) {
      metadata_[[item]] <- if (item %in% names(metadata)) metadata[[item]] else NA_character_
    }
    metadata <- metadata_
  }
  
  # return metadata
  metadata
}


# write metadata
.write_rvat_filetype <- function(filetype, con) {
  writeLines(sprintf("# RVAT-%s", filetype), con = con)
}

.write_metadata <- function(metadata, con) {
  for (item in names(metadata)) {
    writeLines(sprintf("# %s: %s", item, metadata[[item]]), con = con)
  }
}

.write_rvat_header <- function(filetype, metadata, con) {
  .write_rvat_filetype(filetype, con)
  .write_metadata(metadata, con)
}

# input checks
.check_gdb_ids <- function(gdb, object) {
  # check if gdb id in object and gdb match
  if ( !is.null(getGdbId(object)) && !is.na(getGdbId(object)) &
       !is.null(getGdbId(gdb)) && !is.na(getGdbId(gdb))) {
    if (getGdbId(gdb) != getGdbId(object)) {
      stop (sprintf("The %s seems to be generated from a different gdb than supplied. Please check using `getGdbId`. Set `strict` = FALSE to ignore.", as.character(class(object))))
    }
  }
  
  if ( !is.null(getRvatVersion(object)) && !is.na(getRvatVersion(object)) &
       !is.null(getRvatVersion(gdb)) && !is.na(getRvatVersion(gdb))) {
    if (getRvatVersion(gdb) != getRvatVersion(object)) {
      warning (sprintf("The gdb and %s were generated using different RVAT versions", as.character(class(object))))
    }
  }
}