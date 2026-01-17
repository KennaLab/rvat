#' @rdname mapToCDS
#' @usage NULL
#' @export
setMethod(
  "mapToCDS",
  signature = signature(object = "gdb"),
  definition = function(
    object,
    gff,
    exonPadding = 12L,
    output = NULL,
    gene_id = NULL,
    transcript_id = NULL,
    biotype = NULL,
    verbose = TRUE
  ) {
    # input validation
    .mapToCDS_validate_input(as.list(environment()))

    # load and filter gtf/gtf -> can be a path or a GRanges object
    exons <- .mapToCDS_load_and_filter_gff(
      gff,
      gene_id,
      transcript_id,
      biotype,
      verbose
    )

    # check chromosome overlap
    chroms <- .mapToCDS_get_overlapping_chroms(object, exons)

    # initialize output
    if (!is.null(output)) {
      output_con <- file(output, open = "wb")
      on.exit(close(output_con), add = TRUE)
      # write header
      write.table(
        data.frame(
          VAR_id = character(),
          CHROM = character(),
          POS = integer(),
          gene_id = character(),
          transcript_id = character(),
          cdsPOS = integer(),
          padded = logical(),
          stringsAsFactors = FALSE
        ),
        file = output_con,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
      )
    }

    # loop through chromosomes and transcripts
    results_list <- list()

    has_chr <- any(grepl(
      "^chr",
      DBI::dbGetQuery(object, "SELECT CHROM FROM var_ranges LIMIT 1")$CHROM
    ))
    for (chr in chroms) {
      if (verbose) {
        message(sprintf("Mapping chromosome %s", chr))
      }

      # get variants for the current chromosome
      vars_chr <- .mapToCDS_get_variants_on_chrom(
        object,
        chr,
        has_chr = has_chr
      )

      # get exons for this chromosome
      exons_chr <- exons[seqnames(exons) == chr]

      # process all transcripts on this chromosome
      res_chr <- .mapToCDS_process_chromosome(
        vars_chr,
        exons_chr,
        exonPadding
      )

      # write to ouput if specified
      if (!is.null(res_chr)) {
        if (!is.null(output)) {
          write.table(
            res_chr,
            file = output_con,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            append = TRUE
          )
        } else {
          results_list[[chr]] <- res_chr
        }
      }
    }

    # return results
    if (is.null(output)) {
      results <- do.call(rbind, results_list)
      rownames(results) <- NULL
      results
    } else {
      invisible(NULL)
    }
  }
)

.mapToCDS_validate_input <- function(args) {
  check_wrapper(check_number_whole, args, "exonPadding")
  check_wrapper(check_character, args, "output", allow_null = TRUE)
  check_wrapper(check_character, args, "gene_id", allow_null = TRUE)
  check_wrapper(check_character, args, "transcript_id", allow_null = TRUE)
  check_wrapper(check_character, args, "biotype", allow_null = TRUE)
  check_wrapper(check_bool, args, "verbose")

  invisible(NULL)
}

.mapToCDS_load_and_filter_gff <- function(
  gff,
  gene_id,
  transcript_id,
  biotype,
  verbose
) {
  # load and filter gtf/gtf -> can be a path or a GRanges object
  if (is.character(gff)) {
    gtf <- rtracklayer::import(gff)
  } else if (is(gff, "GRanges")) {
    gtf <- gff
  } else {
    stop(
      "`gff` parameter should be either a path to a gff/gtf-file or a GRanges object",
      call. = FALSE
    )
  }

  # filter based on gene_id/transcript_id/biotype if specified
  if (!is.null(gene_id)) {
    gtf <- gtf[!is.na(gtf$gene_id) & gtf$gene_id %in% gene_id]
  }

  if (!is.null(transcript_id)) {
    gtf <- gtf[
      !is.na(gtf$transcript_id) & gtf$transcript_id %in% transcript_id
    ]
  }

  if (!is.null(biotype)) {
    gtf <- gtf[
      (!is.na(gtf$transcript_biotype) & gtf$transcript_biotype %in% biotype) &
        (!is.na(gtf$gene_biotype) & gtf$gene_biotype %in% biotype)
    ]
  }

  # subset to CDS only
  exons <- gtf[gtf$type == "CDS"]
  GenomeInfoDb::seqlevelsStyle(exons) <- "NCBI"

  if (length(exons) == 0L) {
    warning("No CDS features found after filtering.", call. = FALSE)
  }

  # return exons
  exons
}

.mapToCDS_get_overlapping_chroms <- function(object, exons) {
  chrom <- unique(as.character(seqnames(exons)))
  chrom_gdb <- RSQLite::dbGetQuery(object, "select CHROM from var_ranges")
  if (any(grepl("chr", chrom_gdb$CHROM, fixed = TRUE))) {
    chrom_gdb <- gsub("chr", "", chrom_gdb$CHROM, fixed = TRUE)
  }
  chrom <- intersect(chrom, chrom_gdb)
  if (length(chrom) == 0L) {
    stop(
      "No overlap between chromosomes included in gtf and chromosomes included in gdb.",
      call. = FALSE
    )
  }
  chrom
}

.mapToCDS_get_variants_on_chrom <- function(object, chrom, has_chr) {
  # handle potential 'chr' prefix for query
  query_chrom <- gsub("^chr", "", chrom)
  if (has_chr) {
    query_chrom <- paste0("chr", query_chrom)
  }

  # extract GRanges for variants on chromosome
  vars <- getAnno(
    object,
    "var_ranges",
    where = sprintf("CHROM = '%s'", query_chrom)
  )

  # return early with empty Granges if no variants found
  if (nrow(vars) == 0L) {
    return(GRanges())
  }

  # unserialze and normalize seqlevelsStyle
  gr <- unserialize(vars$ranges[[1]])
  GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"

  # return
  gr
}


.mapToCDS_process_chromosome <- function(vars, exons, exonPadding) {
  # Return early if no variants or exons to process
  if (length(vars) == 0L || length(exons) == 0L) {
    return(NULL)
  }

  # get unique transcripts
  transcripts <- unique(exons$transcript_id)
  results <- list()

  # process each transcript
  for (transcript in transcripts) {
    # subset and sort exons for current transcript
    exons_tx <- sort(exons[exons$transcript_id %in% transcript])
    strand_tx <- unique(as.character(strand(exons_tx)))

    # find variants overlapping exons (with optional padding)
    add <- if (!is.null(exonPadding)) exonPadding else 0L
    vars_tx <- subsetByOverlaps(vars, exons_tx, maxgap = add)
    if (length(vars_tx) == 0L) {
      next
    }

    # calculate gaps between exons for CDS pos adjustment
    gaps_tx <- GenomicRanges::gaps(exons_tx)
    gaps_tx <- gaps_tx[
      as.character(GenomicRanges::strand(gaps_tx)) == strand_tx &
      as.character(GenomicRanges::seqnames(gaps_tx)) == as.character(GenomicRanges::seqnames(exons_tx)[1])
    ]
    
    # calculate cumulative gap widths for position adjustment
    # deductionsPlus: for positive strand
    # deductionsMinus: for negative strand
    deductionsPlus <- cumsum(width(gaps_tx))
    
    if (length(deductionsPlus) == 1L) {
      deductionsMinus <- 0L
    } else if (length(deductionsPlus) > 1L) {
      deductionsMinus <- rev(c(
        0L,
        cumsum(rev(width(GenomicRanges::gaps(exons_tx)[2:length(exons_tx)])))
      ))
    }
    exons_tx$deductionsPlus <- deductionsPlus
    exons_tx$deductionsMinus <- deductionsMinus

    # find variants directly overlapping exons
    overlaps <- GenomicRanges::findOverlaps(vars_tx, exons_tx)

    # find variants in padded regions (if padding is specified)
    overlaps_padded <- NULL
    if (!is.null(exonPadding)) {
      overlaps_padded <- GenomicRanges::findOverlaps(
        vars_tx,
        exons_tx,
        maxgap = exonPadding
      )
      # exclude variants that directly overlap (already in 'overlaps')
      overlaps_padded <- overlaps_padded[
        !S4Vectors::queryHits(overlaps_padded) %in%
          S4Vectors::queryHits(overlaps)
      ]
    }

    # process padded variants: adjust their positions to nearest exon boundary
    if (length(overlaps_padded) > 0L) {
      vars_tx_padded <- vars_tx[S4Vectors::queryHits(overlaps_padded)]
      
      # calculate distance to start and end of nearest exon
      vars_tx_padded$distStart <- start(exons_tx[S4Vectors::subjectHits(
        overlaps_padded
      )]) -
        end(vars_tx_padded) -
        1L
      vars_tx_padded$distEnd <- start(vars_tx_padded) -
        end(exons_tx[S4Vectors::subjectHits(overlaps_padded)]) -
        1L
      
      # update position to nearest exon boundary
      vars_tx_padded$POS_updated <- ifelse(
        vars_tx_padded$distStart >= 0L &
          vars_tx_padded$distStart <= exonPadding,
        start(vars_tx_padded) + vars_tx_padded$distStart + 1L,
        start(vars_tx_padded) - vars_tx_padded$distEnd - 1L
      )
      
      # update ranges with adjusted positions
      ranges(vars_tx_padded) <- IRanges::IRanges(
        start = vars_tx_padded$POS_updated,
        end = vars_tx_padded$POS_updated
      )
      
      # clean up temporary columns
      vars_tx_padded$distStart <- vars_tx_padded$distEnd <- vars_tx_padded$POS_updated <- NULL
      
      # add exon metadata (gene_id, transcript_id, deductions)
      mcols(vars_tx_padded) <- cbind(
        mcols(vars_tx_padded),
        mcols(exons_tx[S4Vectors::subjectHits(overlaps_padded)])[, c(
          "gene_id",
          "transcript_id",
          "deductionsPlus",
          "deductionsMinus"
        )]
      )
      vars_tx_padded$padded <- TRUE
    }

    # process variants with direct exon overlaps
    if (length(overlaps) > 0L) {
      vars_tx <- vars_tx[S4Vectors::queryHits(overlaps)]
      
      # add exon metadata
      mcols(vars_tx) <- cbind(
        mcols(vars_tx),
        mcols(exons_tx[S4Vectors::subjectHits(overlaps)])[, c(
          "gene_id",
          "transcript_id",
          "deductionsPlus",
          "deductionsMinus"
        )]
      )
      vars_tx$padded <- FALSE
    } else {
      vars_tx <- vars_tx[0, ]
    }

    # combine direct and padded overlaps
    if (length(overlaps) > 0L && length(overlaps_padded) > 0L) {
      vars_tx <- sort(c(vars_tx, vars_tx_padded))
    } else if (length(overlaps_padded) > 0L && length(overlaps) == 0L) {
      vars_tx <- vars_tx_padded
    }

    # calculate CDS positions accounting for introns (gaps)
    if (length(vars_tx) > 0L) {
      if (strand_tx == "+") {
        # for positive strand: position minus cumulative gap widths
        mcols(vars_tx)[["cdsPOS"]] <- start(vars_tx) - vars_tx$deductionsPlus
      } else if (strand_tx == "-") {
        # for negative strand: reverse calculation from transcript end
        mcols(vars_tx)[["cdsPOS"]] <- rep(
          max(end(exons_tx)) + 1L,
          length(vars_tx)
        ) -
          start(vars_tx) -
          vars_tx$deductionsMinus
      }
      
      # convert to data frame with selected columns
      names(vars_tx) <- NULL
      vars_tx <- as.data.frame(vars_tx)[, c(
        "VAR_id",
        "seqnames",
        "start",
        "gene_id",
        "transcript_id",
        "cdsPOS",
        "padded"
      )]
      colnames(vars_tx) <- c(
        "VAR_id",
        "CHROM",
        "POS",
        "gene_id",
        "transcript_id",
        "cdsPOS",
        "padded"
      )
      vars_tx[["CHROM"]] <- as.character(vars_tx[["CHROM"]])

      results[[transcript]] <- vars_tx
    }
  }

  # combine results
  do.call(rbind, results)
}


