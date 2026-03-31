get_anno_test_data <- function() {
  list(
    a = data.frame(
      VAR_id = 1:10,
      gene = LETTERS[10:1],
      width = c(5, 9, 59, 65, 39, 10, 12, 94, 52, 11)
    ),
    b = data.frame(
      VAR_id = 5:10,
      gene = LETTERS[6:1],
      transcript = LETTERS[16:21]
    )
  )
}

test_that("uploadAnno basic functionality works", {
  gdb <- create_example_gdb()

  # upload file-based annotations, expect no error
  expect_no_error(suppressMessages(uploadAnno(
    gdb,
    name = "varinfo2",
    value = rvat_example("rvatData.varinfo")
  )))

  # check if uploaded annotation appears in metadata
  annotations <- listAnno(gdb)
  expect_true("varinfo2" %in% annotations$name)

  # upload data without remapping (i.e. VAR_id is provided)
  test_anno <- get_anno_test_data()

  ## upload, expect no error
  expect_no_error(suppressMessages(uploadAnno(
    gdb,
    name = "test_anno",
    value = test_anno$a,
    skipRemap = TRUE
  )))

  ## verify uploaded data
  test_anno_retrieved <- getAnno(gdb, "test_anno")
  expect_equal(test_anno_retrieved, test_anno$a)
})

test_that("getAnno joins work", {
  gdb <- create_example_gdb()
  test_anno <- get_anno_test_data()
  suppressMessages(uploadAnno(
    gdb,
    name = "a",
    value = test_anno$a,
    skipRemap = TRUE
  ))
  suppressMessages(uploadAnno(
    gdb,
    name = "b",
    value = test_anno$b,
    skipRemap = TRUE
  ))

  # test inner join
  inner <- getAnno(gdb, "a", inner = "b")
  expect_equal(inner$VAR_id, 5:10)
  expect_equal(inner[, 2], inner[, 4])

  # test left join
  left <- getAnno(gdb, "var", left = c("a", "b"))
  expect_true(all(c("width", "gene", "transcript") %in% colnames(left)))
  var_count <- nrow(getAnno(gdb, "var"))
  expect_equal(nrow(left), var_count)
})

test_that("uploadAnno variant mapping works", {
  gdb <- create_example_gdb()

  # test variants for mapping
  var <- getAnno(gdb, "var", VAR_id = c(5, 6, 100, 193, 900))

  # annotation with CHROM/POS/REF/ALT (some unmappable)
  mapping_data <- data.frame(
    CHROM = c(var$CHROM, c("chrX", "chrY")),
    POS = c(var$POS, 1000, 10000),
    REF = c(var$REF, "A", "C"),
    ALT = c(var$ALT, "G", "T"),
    gene = LETTERS[7:1],
    width = c(5, 9, 59, 65, 39, 10, 12)
  )

  # upload with default behavior (ignore unmapped, should yield a warning)
  expect_warning(uploadAnno(
    gdb,
    name = "mapped_anno",
    value = mapping_data,
    verbose = FALSE
  ))

  # upload with keepUnmapped = TRUE (should warn about unmapped variants)
  expect_warning(
    uploadAnno(
      gdb,
      name = "mapped_with_unmapped",
      value = mapping_data,
      keepUnmapped = TRUE,
      verbose = FALSE
    ),
    regexp = "could not be mapped"
  )

  # compare results - should be identical after removing unmapped rows
  result1 <- getAnno(gdb, "mapped_anno")
  result2 <- getAnno(gdb, "mapped_with_unmapped")
  result2_filtered <- result2[!is.na(result2$VAR_id), ]
  expect_equal(result1, result2_filtered)

  # check ignoring alleles
  var <- getAnno(gdb, "var", VAR_id = c(5, 6, 100))[, c(
    "CHROM",
    "POS",
    "REF",
    "ALT"
  )]
  ## switch alleles
  ref <- var$REF
  var$REF <- var$ALT
  var$ALT <- ref

  ## when `ignoreAlleles = TRUE`, none should map
  expect_warning(suppressMessages(uploadAnno(
    gdb,
    name = "ignores_alleles_false",
    value = var,
    ignoreAlleles = FALSE,
    overWrite = TRUE
  )))
  var_ignore_alleles_false <- getAnno(gdb, "ignores_alleles_false")
  expect_equal(nrow(var_ignore_alleles_false), 0L)

  ## when `ignoreAlleles = FALSE`, all should map
  suppressMessages(uploadAnno(
    gdb,
    name = "ignores_alleles_true",
    value = var,
    ignoreAlleles = TRUE,
    overWrite = TRUE
  ))
  var_ignore_alleles_true <- getAnno(gdb, "ignores_alleles_true")
  var_ignore_alleles_true$VAR_id <- NULL
  expect_equal(var, var_ignore_alleles_true)
})

test_that("uploadAnno metadata tracking works", {
  gdb <- create_example_gdb()

  # upload data.frame annotations, and check if metadata is updated correctly
  original_annotations <- listAnno(gdb)$name
  df <- data.frame(VAR_id = 1:2, annotation = LETTERS[1:2])
  uploadAnno(gdb, value = df, name = "A", verbose = FALSE, skipRemap = TRUE)
  annotations_updated <- listAnno(gdb)
  expect_equal(setdiff(annotations_updated$name, original_annotations), "A")
  expect_equal(
    annotations_updated[annotations_updated$name == "A", ]$value,
    "interactive_session"
  )

  ## after dropping, should be back to original
  dropTable(gdb, "A", verbose = FALSE)
  expect_equal(sort(original_annotations), sort(listAnno(gdb)$name))

  # upload annotations from file, and check if metadata is updated correctly
  uploadAnno(
    gdb,
    value = rvat_example("rvatData.varinfo"),
    name = "A",
    verbose = FALSE
  )
  annotations_updated <- listAnno(gdb)
  expect_equal(setdiff(annotations_updated$name, original_annotations), "A")
  expect_true(stringr::str_detect(
    annotations_updated[annotations_updated$name == "A", ]$value,
    "rvatData.varinfo"
  ))
})

test_that("uploadAnno input validation works correctly", {
  gdb <- create_example_gdb()
  test_df <- data.frame(VAR_id = c("1", "2"), anno_field = c("A", "B"))

  # expect message when verbose = TRUE, and a table is overwritten
  varinfo <- getAnno(gdb, "varinfo")
  suppressMessages(expect_message(
    uploadAnno(
      gdb,
      value = varinfo,
      name = "varinfo",
      skipRemap = TRUE,
      verbose = TRUE,
      overWrite = TRUE
    ),
    regexp = "already exists and"
  ))

  # expect error when upload name is protected
  for (name in gdb_protected_tables) {
    expect_error(
      uploadAnno(gdb, value = test_df, name = name, skipRemap = TRUE),
      "already exists as a protected gdb table"
    )
  }

  # expect error when upload name already exists, and overWrite = FALSE
  expect_error(
    {
      uploadAnno(
        gdb,
        value = test_df,
        name = "varinfo",
        skipRemap = TRUE
      )
    },
    regexp = "already in use"
  )

  # expect error when upload name already exists as a cohort table
  expect_error(
    uploadAnno(gdb, value = test_df, name = "pheno", skipRemap = TRUE),
    regexp = "already in use"
  )

  # expect error when upload name includes punctuation marks
  for (mark in c(".", ",", "+", "-", " ")) {
    expect_error(
      {
        uploadAnno(
          gdb,
          value = test_df,
          name = sprintf("anno%stest", mark),
          skipRemap = TRUE
        )
      },
      regexp = "Table name may only contain"
    )
  }

  # expect error when value is neither a data.frame nor a file path
  expect_error(
    {
      uploadAnno(
        gdb,
        value = list(a = 1, b = 2),
        name = "test",
        skipRemap = TRUE
      )
    },
    regexp = "must be a data.frame or a character string"
  )

  expect_error(
    {
      uploadAnno(
        gdb,
        value = data.frame(VAR_id = 1:5),
        name = "test",
        skipRemap = FALSE,
        verbose = FALSE,
        overWrite = TRUE
      )
    },
    regexp = "Annotation table contains only"
  )

  # skipRemap = FALSE, but doesn't include required fields
  expect_error(
    {
      uploadAnno(
        gdb,
        value = data.frame(CHROM = c("chr1"), POS = c(1L)),
        name = "test",
        skipRemap = FALSE,
        verbose = FALSE,
        overWrite = TRUE
      )
    },
    regexp = "Mapping not possible: REF, ALT not provided."
  )

  ## the above should work when ignoreAlleles = TRUE
  ## but should return a warning given that the variant cannot be mapped
  expect_warning(
    {
      uploadAnno(
        gdb,
        value = data.frame(CHROM = c("chr1"), POS = c(1L)),
        name = "test",
        skipRemap = FALSE,
        ignoreAlleles = TRUE,
        verbose = FALSE,
        overWrite = TRUE
      )
    },
    regexp = "could not be mapped"
  )

  # the above should work when ignoreAlleles = TRUE
  expect_error(
    {
      uploadAnno(
        gdb,
        value = data.frame(CHROM = c("chr1")),
        name = "test",
        skipRemap = FALSE,
        ignoreAlleles = TRUE,
        verbose = FALSE,
        overWrite = TRUE
      )
    },
    regexp = "Mapping not possible: POS not provided."
  )
})
