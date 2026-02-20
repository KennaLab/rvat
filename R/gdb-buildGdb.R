#' Create a gdb file.
#'
#' Creates a new [`gdb`] file.
#' The gdb can be structured and populated using a provided vcf file.
#'
#' @param vcf Input vcf file used to structure and populate gdb.
#' Warning, this function makes the following assumptions:
#'   1) strict adherence to vcf format (GT subfield first element in genotype fields),
#'   2) desired genotype QC has already been applied (DP,GQ filters),
#'   3) GT values conform to the set {0/0,0/1,1/0,1/1,./.,0|0,0|1,1|0,1|1,.|.}.
#'   Multiallelic parsing and genotype QC can be performed using vcftools and/or
#'   accompanying parser scripts included in the rvat repository.
#' @param output Path for output [`gdb`] file.
#' @param skipIndexes Flag to skip generation of indexes for var and dosage table (VAR_id;CHROM,POS,REF,ALT).
#' Typically only required if you plan to use [`concatGdb`] to concatenate a series of separately generated gdb files.
#' @param skipVarRanges Flag to skip generation of ranged var table.
#' Typically only useful (i.e., faster) if you plan to use [`concatGdb`] to concatenate a series of separately generated gdb files.
#' @param overWrite Overwrite if `output` already exists? Defaults to `FALSE`, in which case an error is raised.
#' @param genomeBuild Optional genome build to include in the gdb metadata.
#' If specified, it will be used to set ploidies (diploid, XnonPAR, YnonPAR) if the genome build is implemented in RVAT (currently: GRCh37, hg19, GRCh38, hg38).
#' @param memlimit Maximum number of vcf records to parse at a time, defaults to 1000.
#' @param verbose Should the function be verbose? (TRUE/FALSE), defaults to `TRUE`.
#' @example inst/examples/example-buildGdb.R
#' @export
buildGdb <- function(
  vcf,
  output,
  skipIndexes = TRUE, ## DuckDB does indexing internally, so no need to do it manually as well
  skipVarRanges = FALSE,
  overWrite = FALSE,
  genomeBuild = NULL,
  memlimit = 1000L,
  verbose = TRUE
) {
  # validate input
  arg <- as.list(environment())
  .buildgdb_validate_input(arg)

  # check if output exists and overwrite if `overWrite = TRUE`
  .check_output(
    output = output,
    overWrite = overWrite,
    verbose = verbose
  )

  # initialize gdb
  mygdb <- gdb_init(output)
  on.exit(close(mygdb), add = TRUE)

  # create schema
  if (verbose) {
    message(sprintf(
      "%s\tCreating gdb tables",
      as.character(round(Sys.time(), units = "secs"))
    ))
  }
  .buildgdb_create_schema(mygdb)

  # import variant records
  populateGdbFromVcf(mygdb, vcf, memlimit = memlimit, verbose = verbose)

  # generate indexes
  if (!skipIndexes) {
    .gdb_create_indexes(
      gdb = mygdb,
      verbose = verbose
    )
  }

  # generate var_ranges table
  if (!skipVarRanges) {
    if (verbose) {
      message(sprintf(
        "%s\tCreating ranged var table",
        as.character(round(Sys.time(), units = "secs"))
      ))
    }
    addRangedVarinfo(mygdb, overwrite = TRUE, verbose = verbose)
  }

  # populate meta table
  .gdb_populate_meta_table(
    gdb = mygdb,
    genomeBuild = genomeBuild,
    verbose = verbose
  )

  # finished
  if (verbose) {
    message(sprintf(
      "%s\tComplete",
      as.character(round(Sys.time(), units = "secs"))
    ))
  }

  invisible(NULL)
}


setMethod(
  "populateGdbFromVcf",
  signature = "gdb",
  definition = function(object, vcf, memlimit = 1000L, verbose = TRUE) {
    # open vcf connection
    if (vcf == "-") {
      vcf <- "stdin"
    } else if (!file.exists(vcf)) {
      stop(sprintf("Input vcf %s does not exist", vcf), call. = FALSE)
    }
    if (endsWith(vcf, ".gz")) {
      con <- gzcon(file(vcf, open = "rb"))
    } else {
      con <- file(vcf, open = "r")
    }
    ## close connection on exit
    if (!is.null(con) && vcf != "stdin") {
      on.exit(close(con), add = TRUE)
    }

    # skip over vcf meta-data
    header <- NULL
    while (length(i <- readLines(con, n = 1L)) > 0L) {
      if (substr(i, 1, 2) == "##") {
        next
      }
      header <- i
      break
    }

    # parse vcf header line
    ## last line should start with #CHROM
    if (is.null(header) || !grepl("^#CHROM", header)) {
      stop(
        "Invalid or missing vcf header: ",
        "expected line starting with #CHROM.",
        call. = FALSE
      )
    }
    header <- unlist(strsplit(header, split = "\t", fixed = TRUE))
    width <- length(header)
    m <- width - 9L
    if (m <= 0L) {
      stop(
        "Invalid parsing of vcf header line. No samples detected",
        call. = FALSE
      )
    }

    # write SM table
    if (verbose) {
      message(sprintf("%s sample IDs detected", m))
    }
    DBI::dbWriteTable(
      object,
      name = "SM",
      value = data.frame(
        IID = header[10:(m + 9)],
        sex = 0,
        stringsAsFactors = FALSE
      ),
      overwrite = TRUE
    )

    # parse vcf records
    if (verbose) {
      message(sprintf(
        "%s\tParsing vcf records",
        as.character(round(Sys.time(), units = "secs"))
      ))
    }
    counter <- 0L
    while (length(records <- readLines(con, n = memlimit)) > 0L) {
      counter <- counter + length(records)

      # parse records
      records <- stringi::stri_split_fixed(
        records,
        "\t",
        n = width,
        simplify = TRUE
      )

      for (i in seq_len(nrow(records))) {
        # parse multi-allelic site
        if (grepl(",", records[i, 5], fixed = TRUE)) {
          alleles <- unlist(strsplit(records[i, 5], split = ",", fixed = TRUE))
          gt <- substr(records[i, -(1:9), drop = FALSE], 1, 3)
          carrierIndex <- which(!(gt %in% c(".:0", "./.", "0/0", ".|.", "0|0")))
          carrierGt <- strsplit(gt[carrierIndex], split = "\\/|\\|")
          carrierN <- length(carrierIndex)
          for (ai in seq_along(alleles)) {
            # write variant info
            insertVarRecord(
              object,
              record = c(records[i, 1:4], alleles[ai], records[i, 6:9])
            )

            # write genotype data
            gtia <- gt
            for (carrier in seq_len(carrierN)) {
              gtia[carrierIndex[carrier]] <- paste(
                ifelse(carrierGt[[carrier]] == as.character(ai), "1", "0"),
                collapse = "/"
              )
            }
            insertDosageRecord(object, record = gtia)
          }
          # bi-allelic site
        } else {
          insertVarRecord(object, record = records[i, 1:9])
          insertDosageRecord(
            object,
            record = substr(records[i, -(1:9), drop = FALSE], 1, 3)
          )
        }
      }

      # log number of records processed
      if (verbose) {
        message(sprintf(
          "%s\tProcessing completed for %s records. Committing to db.",
          as.character(round(Sys.time(), units = "secs")),
          counter
        ))
      }
    }

    invisible(NULL)
  }
)

setMethod(
  "insertVarRecord",
  signature = "gdb",
  definition = function(object, record) {
    DBI::dbExecute(
      object,
      paste0(
        "insert into var(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) ",
        "values (:chrom, :pos, :id, :ref, :alt, :qual, :flt, :info, :form)"
      ),
      params = list(
        chrom = record[1],
        pos = record[2],
        id = record[3],
        ref = record[4],
        alt = record[5],
        qual = record[6],
        flt = record[7],
        info = record[8],
        form = record[9]
      )
    )
  }
)

setMethod(
  "insertDosageRecord",
  signature = "gdb",
  definition = function(object, record) {
    record[record == "0/0"] <- "0"
    record[record == "./0"] <- "0"
    record[record == "0/1"] <- "1"
    record[record == "1/0"] <- "1"
    record[record == "1/1"] <- "2"
    record[record == "./."] <- "N"
    record[record == "0|0"] <- "0"
    record[record == "0|1"] <- "1"
    record[record == "1|0"] <- "1"
    record[record == "./1"] <- "1"
    record[record == ".|1"] <- "1"
    record[record == "1|1"] <- "2"
    record[record == ".|."] <- "N"
    record[record == ".:0"] <- "N"
    obs <- sort(unique(c(record)))
    if (sum(!obs %in% c("0", "1", "2", "N")) > 0L) {
      stop(
        sprintf(
          paste0(
            "Invalid genotype codes after parsing. ",
            "Expected (0,1,2,N), observed (%s) "
          ),
          paste(obs, collapse = ",")
        ),
        call. = FALSE
      )
    }

    record <- I(list(memCompress(paste(record, collapse = ""), type = "gzip")))
    DBI::dbExecute(
      object,
      "insert into dosage(GT) values (:record)",
      params = list(record = record)
    )
  }
)

## add GenomicRanges var
setMethod(
  "addRangedVarinfo",
  signature = "gdb",
  definition = function(object, overwrite = FALSE, verbose = TRUE) {
    # check if 'var_ranges' is already present in the gdb
    if ("var_ranges" %in% DBI::dbListTables(object)) {
      if (!overwrite) {
        stop(
          "'var_ranges' is already present in the gdb. ",
          "Set `overwrite = TRUE` to overwrite.",
          call. = FALSE
        )
      }
      if (overwrite) {
        if (verbose) {
          message(
            "'var_ranges' is already present in the gdb, it will be overwritten."
          )
        }
        dropTable(object, "var_ranges")
      }
    }

    # generate GRanges for each chromosome and insert into var_ranges
    chroms <- unique(getAnno(object, "var", fields = "CHROM")$CHROM)
    DBI::dbExecute(object, "create table var_ranges (CHROM text, ranges blob)")
    for (chrom in chroms) {
      # get varinfo for current chromosome
      vars <- getAnno(object, "var", where = sprintf("CHROM = '%s'", chrom))

      # convert into genomicRanges
      gr <- GenomicRanges::GRanges(
        seqnames = vars$CHROM,
        ranges = IRanges::IRanges(
          start = vars$POS,
          end = vars$POS + stringr::str_length(vars$REF) - 1L
        ),
        VAR_id = vars$VAR_id
      )
      names(gr) <- vars$VAR_id

      # insert GRanges into gdb as blobs per chromosome
      DBI::dbExecute(
        object,
        statement = "insert into var_ranges(CHROM,ranges) values (:CHROM, :ranges)",
        params = list(
          CHROM = chrom,
          ranges = list(serialize(gr, connection = NULL))
        )
      )
    }
    invisible(NULL)
  }
)

.buildgdb_validate_input <- function(
  args
) {
  check_wrapper(check_character, args, "vcf", length_equal = 1L)
  check_wrapper(check_character, args, "output", length_equal = 1L)
  check_wrapper(check_bool, args, "skipIndexes", length_equal = 1L)
  check_wrapper(check_bool, args, "skipVarRanges", length_equal = 1L)
  check_wrapper(check_bool, args, "overWrite", length_equal = 1L)
  check_wrapper(
    check_character,
    args,
    "genomeBuild",
    allow_null = TRUE,
    length_equal = 1L
  )
  check_wrapper(check_number_whole, args, "memlimit", length_equal = 1L)
  check_positive(args[["memlimit"]], arg = "memlimit", length_equal = 1L)
  check_wrapper(check_bool, args, "verbose", length_equal = 1L)

  invisible(NULL)
}


.buildgdb_create_schema <- function(gdb) {
  DBI::dbExecute(
    gdb,
    paste0(
      "create table var (VAR_id integer primary key, CHROM text, POS int, ",
      "ID text, REF text, ALT text, QUAL text, FILTER text, INFO text, FORMAT text);"
    )
  )
  DBI::dbExecute(gdb, "create table SM (IID text, sex int)")
  DBI::dbExecute(
    gdb,
    "create table dosage (VAR_id integer primary key, GT BLOB);"
  )
  DBI::dbExecute(gdb, "create table anno (name text,value text,date text)")
  DBI::dbExecute(gdb, "create table cohort (name text,value text,date text)")
  DBI::dbExecute(gdb, "create table meta (name text,value text)")

  invisible(NULL)
}
