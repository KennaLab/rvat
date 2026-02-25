create_test_ranges <- function() {
  list(
    basic = data.frame(
      CHROM = c("chr21", "chr21", "chr4", "chrX"),
      start = c(31659820L, 31668483L, 169433614L, 1356275L),
      end = c(31666490L, 31668554L, 169438160L, 1365198L)
    ),
    indel_test = data.frame(
      CHROM = c("chr1", "chr4", "chr16", "chrX"),
      start = c(13112464L + 1L, 169433612L + 1L, 31184283L + 1L, 1382393L + 1L),
      end = c(13112464L + 1L, 169433612L + 1L, 31184283L + 1L, 1382393L + 1L),
      VAR_id = c(644L, 1362L, 1206L, 1765L)
    ),
    no_overlap = data.frame(
      CHROM = "chr21",
      start = 5000000L,
      end = 6000000L
    )
  )
}

test_that("extractRanges basic functionality works", {
  gdb <- create_example_gdb()
  ranges <- create_test_ranges()

  # check basic extraction
  check1 <- extractRanges(gdb, ranges = ranges$basic)
  expect_type(check1, "character")
  expect_gt(length(check1), 0L)
  expect_true(all(check1 %in% getAnno(gdb, "var", fields = "VAR_id")$VAR_id))

  # compare with manual extraction using GenomicRanges
  var <- getAnno(gdb, "var", fields = c("VAR_id", "CHROM", "POS", "REF", "ALT"))
  var_granges <- GRanges(
    seqnames = var$CHROM,
    ranges = IRanges::IRanges(
      start = var$POS,
      end = var$POS + stringr::str_length(var$REF) - 1L
    ),
    VAR_id = var$VAR_id
  )
  check2 <- subsetByOverlaps(
    var_granges,
    makeGRangesFromDataFrame(ranges$basic)
  )$VAR_id

  expect_identical(sort(as.numeric(check1)), sort(as.numeric(check2)))
})

test_that("extractRanges handles different input formats", {
  gdb <- create_example_gdb()

  # ensembl gene info
  gtf <- readr::read_tsv(
    test_path("data/Homo_sapiens.GRCh38.105.gene.txt.gz"),
    show_col_types = FALSE
  )
  colnames(gtf)[1] <- "CHROM"
  gtf <- gtf[gtf$gene_name %in% c("RIN3", "FUS", "ZNF483"), ]

  ## data.frame and GRanges should result in identical output
  check_dataframe <- extractRanges(gdb, ranges = gtf)
  check_granges <- extractRanges(
    gdb,
    ranges = makeGRangesFromDataFrame(gtf)
  )
  expect_identical(check_dataframe, check_granges)

  # 'chr' prefix shouldn't matter
  gtf_with_chr <- gtf %>% dplyr::mutate(CHROM = paste0("chr", CHROM))
  check_chr <- extractRanges(gdb, ranges = gtf_with_chr)
  expect_identical(check_dataframe, check_chr)

  # check if matches varinfo
  check_varinfo <- getAnno(
    gdb,
    "varInfo",
    where = "gene_name in ('RIN3', 'FUS', 'ZNF483')"
  )
  expect_identical(
    sort(check_dataframe),
    sort(as.character(check_varinfo$VAR_id))
  )
})

test_that("extractRanges padding functionality works", {
  gdb <- create_example_gdb()

  # check whether INDELs are extracted correctly
  # includes three deletions that should overlap the start coordinates of ranges
  ranges <- create_test_ranges()$indel_test

  ## without padding, none should be extracted
  check_no_padding <- extractRanges(
    gdb,
    ranges = ranges,
    padding = 0L
  )
  expect_identical(check_no_padding, vector(mode = "character", length = 0L))

  # with padding, should extract overlapping variants
  check_with_padding <- extractRanges(
    gdb,
    ranges = ranges,
    padding = 10L
  )
  expect_identical(sort(check_with_padding), sort(c("1362", "1206", "1765")))
})

test_that("extractRanges handles edge cases correctly", {
  gdb <- create_example_gdb()
  ranges <- create_test_ranges()

  # check no overlap message and empty result
  expect_message(
    var_no_overlap <- getAnno(
      gdb,
      "var",
      fields = c("VAR_id", "CHROM", "POS", "REF", "ALT"),
      ranges = ranges$no_overlap
    ),
    regexp = "No variants overlap with provided range"
  )

  var_no_overlap <- suppressMessages(getAnno(
    gdb,
    "var",
    fields = c("VAR_id", "CHROM", "POS", "REF", "ALT"),
    ranges = ranges$no_overlap
  ))
  expect_identical(nrow(var_no_overlap), 0L)

  # expect warning if no chromosomes overlap between ranges and gdb
  chroms <- unique(getAnno(gdb, "var", fields = "CHROM")$CHROM)
  chroms <- setdiff(paste0("chr", c(1:22, "X", "Y")), chroms)[1:2]

  expect_warning(
    extractRanges(
      gdb,
      ranges = data.frame(
        CHROM = chroms,
        start = c(1L, 2L),
        end = c(2L, 3L)
      )
    ),
    regexp = "None of the specified chromosomes"
  )
})

test_that("extractRanges input validation works correctly", {
  gdb <- create_example_gdb()
  ranges <- create_test_ranges()$basic

  # expect an error when invalid padding values are provided
  expect_error(
    extractRanges(gdb, ranges, padding = "250"),
    regexp = "should be a positive value"
  )

  expect_error(
    extractRanges(gdb, valid_ranges, padding = -5),
    regexp = "should be a positive value"
  )

  # expect error if required columns are missing
  expect_error(
    extractRanges(gdb, ranges = data.frame(CHROM = "chr21", start = 1000)),
    regexp = "should be present"
  )

  expect_error(
    extractRanges(gdb, ranges = data.frame(start = 1000, end = 1001)),
    regexp = "should be present"
  )

  expect_error(
    extractRanges(gdb, ranges = data.frame(CHROM = "chr21", end = 1001)),
    regexp = "should be present"
  )

  # expect error if input format is not supported
  expect_error(
    extractRanges(gdb, ranges = matrix()),
    regexp = "should be either"
  )

  expect_error(
    extractRanges(gdb, ranges = "hello"),
    regexp = "should be either"
  )
})
