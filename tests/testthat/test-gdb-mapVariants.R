test_that("mapVariants handles different input formats", {
  gdb <- create_example_gdb()

  # outputs for different formats
  output_gtf <- withr::local_tempfile()
  output_bed <- withr::local_tempfile()
  output_ranges <- withr::local_tempfile()
  output_granges <- withr::local_tempfile()

  # test whether different input formats produce identical

  ## gtf
  mapVariants(
    gdb,
    gff = "../data/protein_coding_genes.gtf",
    output = output_gtf,
    verbose = FALSE
  )

  ## bed file
  mapVariants(
    gdb,
    bed = "../data/protein_coding_genes.bed",
    bedCols = c("gene_id", "gene_name"),
    output = output_bed,
    verbose = FALSE
  )

  ## ranges file
  mapVariants(
    gdb,
    ranges = "../data/protein_coding_genes.ranges",
    output = output_ranges,
    verbose = FALSE
  )

  ## GRanges
  ranges <- GenomicRanges::makeGRangesFromDataFrame(
    readr::read_tsv(
      "../data/protein_coding_genes.ranges",
      show_col_types = FALSE
    ),
    keep.extra.columns = TRUE
  )
  mapVariants(
    gdb,
    ranges = ranges,
    output = output_granges,
    verbose = FALSE
  )

  # compare outputs
  output_gtf_data <- readr::read_tsv(output_gtf, show_col_types = FALSE)
  output_bed_data <- readr::read_tsv(output_bed, show_col_types = FALSE)
  output_ranges_data <- readr::read_tsv(output_ranges, show_col_types = FALSE)
  output_granges_data <- readr::read_tsv(output_granges, show_col_types = FALSE)

  # core fields should be identical
  expect_equal(
    output_gtf_data[, c("VAR_id", "gene_id", "gene_name")],
    output_bed_data[, c("VAR_id", "gene_id", "gene_name")]
  )

  expect_equal(
    output_gtf_data[, c("VAR_id", "gene_id", "gene_name")],
    output_ranges_data[, c("VAR_id", "gene_id", "gene_name")]
  )

  expect_equal(
    output_gtf_data[, c("VAR_id", "gene_id", "gene_name")],
    output_granges_data[, c("VAR_id", "gene_id", "gene_name")]
  )

  # also test w/o specifying output
  genes <- mapVariants(
    gdb,
    gff = "../data/protein_coding_genes.gtf",
    verbose = FALSE
  )
  output_gtf_data$score <- as.numeric(output_gtf_data$score)
  output_gtf_data$phase <- as.numeric(output_gtf_data$phase)
  expect_equal(as.data.frame(output_gtf_data), genes)
})

test_that("mapVariants matches vcfanno output", {
  gdb <- create_example_gdb()

  # run mapVariants
  output_file <- withr::local_tempfile()
  mapVariants(
    gdb,
    gff = "../data/protein_coding_genes.gtf",
    output = output_file,
    verbose = FALSE
  )

  output_data <- readr::read_tsv(output_file, show_col_types = FALSE)

  # compare with vcfanno
  vcfanno <- readr::read_tsv(
    "../data/rvatData.parsed.vcfAnno.gene.vcfInfo2Table",
    show_col_types = FALSE
  )

  expect_equal(
    output_data %>%
      dplyr::select(VAR_id, gene_id, gene_name) |>
      dplyr::arrange(VAR_id, gene_id),
    vcfanno %>%
      dplyr::select(VAR_id = ID, gene_id, gene_name) |>
      dplyr::arrange(VAR_id, gene_id)
  )
})


test_that("mapVariants field selection works", {
  gdb <- create_example_gdb()

  ranges <- data.frame(
    CHROM = "chr21",
    start = 31659666,
    end = 31668931,
    gene_name = "SOD1",
    gene_name2 = "SOD1"
  )
  output <- mapVariants(
    gdb,
    ranges = ranges,
    fields = c("gene_name"),
    verbose = FALSE,
    overWrite = TRUE
  )
  expect_equal(colnames(output), c("VAR_id", "gene_name"))

  # include fields that are standard GRanges fields
  output <- mapVariants(
    gdb,
    ranges = ranges,
    fields = c("seqnames", "start", "end", "gene_name"),
    verbose = FALSE,
    overWrite = TRUE
  )
  expect_equal(colnames(output), c("VAR_id", "seqnames", "start", "end", "gene_name"))
})

test_that("mapVariants produces identical results to uploadAnno for single positions", {
  gdb <- create_example_gdb()

  # write variant info as ranges file
  varInfo <- getAnno(gdb, "varInfo")
  varinfo_path <- withr::local_tempfile()
  readr::write_tsv(
    varInfo %>%
      dplyr::mutate(start = POS, end = POS) |>
      dplyr::select(CHROM, POS, start, end, dplyr::everything()) |>
      dplyr::select(-VAR_id),
    file = varinfo_path
  )

  # map using uploadAnno
  uploadAnno(
    gdb,
    name = "varInfo_uploadanno",
    value = varinfo_path,
    ignoreAlleles = TRUE,
    verbose = FALSE
  )

  # map using mapVariants
  mapVariants(
    gdb,
    ranges = varinfo_path,
    uploadName = "varInfo_mapvariants",
    verbose = FALSE
  )

  # compare results
  varInfo_uploadanno <- getAnno(gdb, "varInfo_uploadanno")
  varInfo_uploadanno$start <- NULL
  varInfo_uploadanno$end <- NULL
  varInfo_uploadanno$CHROM <- NULL

  varInfo_mapVariants <- getAnno(gdb, "varInfo_mapvariants")

  ## IDs for matching
  varInfo_mapVariants$ID_tmp <- sprintf(
    "%s_%s_%s",
    varInfo_mapVariants$VAR_id,
    varInfo_mapVariants$REF,
    varInfo_mapVariants$ALT
  )
  varInfo_uploadanno$ID_tmp <- sprintf(
    "%s_%s_%s",
    varInfo_uploadanno$VAR_id,
    varInfo_uploadanno$REF,
    varInfo_uploadanno$ALT
  )

  ## align
  varInfo_mapVariants <- varInfo_mapVariants[
    match(varInfo_uploadanno$ID_tmp, varInfo_mapVariants$ID_tmp),
  ]

  ## cleanup
  rownames(varInfo_mapVariants) <- NULL
  rownames(varInfo_uploadanno) <- NULL
  varInfo_uploadanno$ID <- ifelse(
    varInfo_uploadanno$ID == "NA",
    NA_character_,
    varInfo_uploadanno$ID
  )
  varInfo_uploadanno$FILTER <- ifelse(
    varInfo_uploadanno$FILTER == "NA",
    NA_character_,
    varInfo_uploadanno$FILTER
  )
  varInfo_mapVariants$FILTER <- as.character(varInfo_mapVariants$FILTER)

  # compare
  expect_equal(varInfo_uploadanno, varInfo_mapVariants)
})


test_that("mapVariants input validation works correctly", {
  gdb <- create_example_gdb()
  output_file <- withr::local_tempfile()

  # expect error when no input is specified
  expect_error(
    {
      mapVariants(gdb, output = output_file, verbose = FALSE)
    },
    regexp = "At least one of"
  )

  # expect error when multiple inputs are specified
  expect_error(
    mapVariants(
      gdb,
      gff = "../data/protein_coding_genes.gtf",
      bed = "../data/protein_coding_genes.bed",
      output = output_file,
      verbose = FALSE
    ),
    "Multiple region inputs specified"
  )

  # expect error when required fields are missing
  expect_error(
    {
      mapVariants(
        gdb,
        ranges = data.frame(chr = "chr1", start = 1000, end = 1001),
        verbose = FALSE
      )
    },
    "should be present in"
  )

  # expect error when ranges is of wrong type
  expect_error(
    {
      mapVariants(
        gdb,
        ranges = FALSE,
        verbose = FALSE
      )
    },
    "`ranges` must be a"
  )

  # expect error when ranges/gff/bed file does not exist
  tmpfile <- withr::local_tempfile()
  expect_error(
    {
      mapVariants(
        gdb,
        ranges = tmpfile,
        verbose = FALSE
      )
    },
    "does not exist"
  )
  expect_error(
    {
      mapVariants(
        gdb,
        bed = tmpfile,
        verbose = FALSE
      )
    },
    "does not exist"
  )
  expect_error(
    {
      mapVariants(
        gdb,
        gff = tmpfile,
        verbose = FALSE
      )
    },
    "does not exist"
  )

  # expect error when upload name includes punctuation marks
  for (mark in c(".", ",", "+", "-", " ")) {
    expect_error(
      {
        mapVariants(
          gdb,
          gff = "../data/protein_coding_genes.gtf",
          uploadName = sprintf("upload%sname", mark),
          verbose = FALSE
        )
      },
      "Table name may only contain"
    )
  }

  # expect error when a protected uploadName is specified
  for (name in c(
    "dosage",
    "meta",
    "anno",
    "var",
    "var_ranges",
    "SM",
    "cohort",
    "tmp"
  )) {
    expect_error(
      {
        mapVariants(
          gdb,
          gff = "../data/protein_coding_genes.gtf",
          uploadName = name,
          verbose = FALSE
        )
      },
      "already exists as a protected gdb table"
    )
  }

  # expect error when uploadName is of wrong type
  expect_error(
    {
      mapVariants(
        gdb,
        gff = "../data/protein_coding_genes.gtf",
        uploadName = TRUE,
        verbose = FALSE
      )
    },
    "must be a character vector"
  )

  # expect error when output exists and overWrite = FALSE
  expect_error(
    {
      mapVariants(
        gdb,
        gff = "../data/protein_coding_genes.gtf",
        uploadName = "varInfo",
        verbose = FALSE
      )
    },
    "already in use"
  )

  # overwrite with overWrite = TRUE
  uploadAnno(
    gdb,
    value = data.frame(VAR_id = 1:10, gene_name = LETTERS[1:10]),
    name = "gene",
    verbose = FALSE,
    skipRemap = TRUE
  )
  meta <- listAnno(gdb)
  output <- capture_messages({
    mapVariants(
      gdb,
      gff = "../data/protein_coding_genes.gtf",
      uploadName = "gene",
      verbose = TRUE,
      overWrite = TRUE
    )
  })
  expect_true(any(stringr::str_detect(output, "already exists")))
  expect_true(any(stringr::str_detect(output, "removed from gdb")))
  meta_updated <- listAnno(gdb)
  expect_true(
    meta_updated[meta_updated$name == "gene", ]$value == "mapVariants"
  )
})
