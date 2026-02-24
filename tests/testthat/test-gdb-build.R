test_that("buildGdb runs without errors", {
  # build gdb
  expect_no_error(suppressMessages(buildGdb(
    vcf = rvat_example("rvatData.vcf.gz"),
    output = testgdb,
    genomeBuild = "GRCh38"
  )))
  expect_no_error({
    gdb <- gdb(testgdb)
  })

  # required tables are present
  tables <- DBI::dbListTables(gdb)
  expected_tables <- c("meta", "var", "SM", "dosage", "var_ranges")
  expect_true(all(expected_tables %in% tables))
})

test_that("buildGdb creates correct metadata", {
  gdb <- gdb(testgdb)
  meta <- RSQLite::dbGetQuery(gdb, "SELECT * FROM meta")

  # check ID length (should be 28 characters)
  expect_equal(
    stringr::str_length(meta[meta$name == "id", ]$value),
    28L
  )

  # check genome build
  expect_equal(
    meta[meta$name == "genomeBuild", ]$value,
    "GRCh38"
  )

  # check creation date format (should have 2 dashes)
  expect_equal(
    stringr::str_count(meta[meta$name == "creationDate", ]$value, "-"),
    2L
  )

  # check RVAT version format (should have 2 dots)
  expect_equal(
    stringr::str_count(meta[meta$name == "rvatVersion", ]$value, "\\."),
    2L
  )
})


test_that("buildGdb creates correct database indexes", {
  gdb <- gdb(testgdb)
  expect_gdb_indexes(gdb)
})


test_that("buildGdb creates consistent variant ranges", {
  gdb <- gdb(testgdb)

  # check if unpacked ranged varinfo is identical to varinfo
  chroms <- unique(getAnno(gdb, "var", fields = "CHROM")$CHROM)
  granges <- lapply(chroms, function(chrom) {
    gr <- unserialize(getAnno(
      gdb,
      "var_ranges",
      where = sprintf("CHROM = '%s'", chrom)
    )$ranges[[1]])
    GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"
    gr
  })

  granges <- suppressWarnings(do.call("c", unlist(granges)))
  granges$ID <- paste(
    GenomicRanges::seqnames(granges),
    BiocGenerics::start(granges),
    sep = ":"
  )

  ## compare with variant table
  var <- getAnno(gdb, "var")
  var$ID <- paste(var$CHROM, var$POS, sep = ":")
  var$ID <- gsub("^chr", "", var$ID)

  expect_equal(sort(granges$ID), sort(var$ID))
})


test_that("buildGdb handles multi-allelic variants correctly", {
  # parse multi-allelic vcf with `parseMultiAllelic.py`
  out_parsed <- withr::local_tempfile()
  vcf_path <- test_path("data/vcf_multiallelic.vcf")
  pyscript <- system.file("exec", "parseMultiAllelic.py", package = "rvat")
  system(
    sprintf(
      "python3 %s %s > %s",
      pyscript,
      vcf_path,
      out_parsed
    )
  )
  gdb_path1 <- withr::local_tempfile()
  gdb_path2 <- withr::local_tempfile()

  # build from multi-allelic vcf directly
  buildGdb(
    test_path("data/vcf_multiallelic.vcf"),
    output = gdb_path1,
    genomeBuild = "GRCh38",
    verbose = FALSE
  )

  # build from parsed vcf
  buildGdb(
    out_parsed,
    output = gdb_path2,
    genomeBuild = "GRCh38",
    verbose = FALSE
  )

  # compare gdbs
  gdb1 <- gdb(gdb_path1)
  gdb2 <- gdb(gdb_path2)
  compare_gdbs(gdb1, gdb2)
})

test_that("buildGdb VCF validation works", {
  # expect error when VCF is missing header
  vcf_noheader <- withr::local_tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##contig=<ID=chr1,length=249250621>",
      "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1"
    ),
    vcf_noheader
  )

  expect_error(
    buildGdb(
      vcf = vcf_noheader,
      output = withr::local_tempfile(fileext = ".gdb"),
      genomeBuild = "GRCh38",
      verbose = FALSE
    ),
    regexp = "Invalid or missing vcf header"
  )

  # expect error when VCF is missing required FORMAT field
  vcf_nosamples <- withr::local_tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##contig=<ID=chr1,length=249250621>",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
      "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT"
    ),
    vcf_nosamples
  )

  expect_error(
    buildGdb(
      vcf = vcf_nosamples,
      output = withr::local_tempfile(fileext = ".gdb"),
      genomeBuild = "GRCh38",
      verbose = FALSE
    ),
    regexp = "No samples detected"
  )

  # expect error when VCF has invalid genotype format
  vcf_invalid_gt <- withr::local_tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##contig=<ID=chr1,length=249250621>",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1",
      "chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\tINVALID",
      "chr1\t200\t.\tC\tT\t30\tPASS\t.\tGT\t5/6",
      "chr1\t300\t.\tG\tA\t30\tPASS\t.\tGT\tabc"
    ),
    vcf_invalid_gt
  )

  expect_error(
    buildGdb(
      vcf = vcf_invalid_gt,
      output = withr::local_tempfile(fileext = ".gdb"),
      genomeBuild = "GRCh38",
      verbose = FALSE
    ),
    regexp = "Invalid genotype codes after parsing."
  )
})


test_that("buildGdb input validation works correctly", {
  # expect error when output gdb already exists
  expect_error(
    {
      buildGdb(
        vcf = rvat_example("rvatData.vcf.gz"),
        output = testgdb,
        genomeBuild = "GRCh38"
      )
    },
    "already exists"
  )

  # expect error when output is incorrectly specified
  expect_error(
    {
      buildGdb(
        vcf = rvat_example("rvatData.vcf.gz"),
        output = c(tempfile(), tempfile()),
        genomeBuild = "GRCh38"
      )
    },
    "should be of length 1"
  )

  # expect message when output is overwritten
  tmp_output <- withr::local_tempfile(fileext = ".gdb")
  readr::write_lines(1:5, tmp_output)
  suppressMessages(expect_message(
    {
      buildGdb(
        vcf = test_path("data/vcf_multiallelic.vcf"),
        output = tmp_output,
        genomeBuild = "GRCh38",
        overWrite = TRUE
      )
    },
    regexp = "already exists and is overwritten"
  ))

  # expect error when file exists and cannot be overwritten
  tmp_dir <- withr::local_tempdir()
  readr::write_lines(1:5, file = paste0(tmp_dir, "/test.txt"))
  expect_warning(expect_error(
    {
      buildGdb(
        vcf = rvat_example("rvatData.vcf.gz"),
        output = tmp_dir,
        overWrite = TRUE,
        genomeBuild = "GRCh38"
      )
    },
    "Failed to remove"
  ))

  # expect error when multiple vcfs are supplied
  expect_error(
    {
      buildGdb(
        vcf = c(
          rvat_example("rvatData.vcf.gz"),
          rvat_example("rvatData.vcf.gz")
        ),
        output = tempfile(),
        genomeBuild = "GRCh38",
        verbose = FALSE
      )
    },
    "should be of length 1"
  )

  # expect error when adding ranged varinfo when already present
  gdb <- gdb(testgdb)
  expect_error(
    {
      addRangedVarinfo(gdb)
    },
    regexp = "is already present"
  )

  # addRangedVarinfo works with overwrite = TRUE
  gdb <- create_example_gdb()
  suppressMessages(expect_message(
    addRangedVarinfo(
      gdb,
      overwrite = TRUE,
      verbose = TRUE
    ),
    regexp = "it will be overwritten"
  ))

  # expect warning when unsupported genome build is supplied
  ## small gdb for testing
  tmp <- withr::local_tempfile(fileext = ".vcf")
  writeVcf(gdb, VAR_id = 1:5, output = tmp, verbose = FALSE)
  expect_warning(
    {
      buildGdb(
        vcf = tmp,
        output = tempfile(),
        genomeBuild = "GRCh39",
        verbose = FALSE
      )
    },
    "not supported by RVAT"
  )
})
