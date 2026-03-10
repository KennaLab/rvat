test_that("buildVarSet basic functionality works", {
  gdb <- create_example_gdb()
  genes <- c("CYP19A1", "FUS", "OPTN", "SOD1")
  where_clause_genes <- sprintf(
    "gene_name in (%s)",
    paste(shQuote(genes), collapse = ", ")
  )

  # build moderate impact varset
  moderate_path <- withr::local_tempfile()
  expect_no_error(suppressMessages(buildVarSet(
    object = gdb,
    output = moderate_path,
    varSetName = "Moderate",
    unitTable = "varInfo",
    unitName = "gene_name",
    where = paste0(
      "(ModerateImpact = 1 or HighImpact = 1) and ",
      where_clause_genes
    )
  )))
  varsetfile_moderate <- varSetFile(moderate_path)
  expect_identical(sort(listUnits(varsetfile_moderate)), sort(genes))

  # test LOF varset
  LOF_path <- withr::local_tempfile()
  expect_no_error(buildVarSet(
    object = gdb,
    output = LOF_path,
    varSetName = "HighImpact",
    unitTable = "varInfo",
    unitName = "gene_name",
    where = paste("HighImpact = 1 and", where_clause_genes),
    verbose = FALSE
  ))
  varsetfile_LOF <- varSetFile(LOF_path)
  expect_identical(sort(listUnits(varsetfile_LOF)), sort(genes))

  # test CADD weights varset
  CADD_path <- withr::local_tempfile()
  expect_no_error(suppressMessages(buildVarSet(
    object = gdb,
    output = CADD_path,
    varSetName = "CADD",
    unitTable = "varInfo",
    unitName = "gene_name",
    weightName = "CADDphred",
    where = where_clause_genes
  )))
  varsetfile_CADD <- varSetFile(CADD_path)
  expect_identical(sort(listUnits(varsetfile_CADD)), sort(genes))

  # get annotation data
  anno <- getAnno(gdb, "varInfo")
  anno$Moderate <- ifelse(anno$ModerateImpact == 1 | anno$HighImpact == 1, 1, 0)
  anno$CADD <- anno$CADDphred

  # build varsets from data.frame
  moderate_df <- buildVarSet(
    anno,
    unitName = "gene_name",
    fields = c("Moderate")
  )

  LOF_df <- buildVarSet(
    anno,
    unitName = "gene_name",
    fields = c("HighImpact")
  )

  CADD_df <- buildVarSet(
    anno,
    unitName = "gene_name",
    fields = c("CADD")
  )

  # compare varsetfiles from gdb vs data.frame
  compare_varsetfile(varsetfile_moderate, moderate_df, units = genes)
  compare_varsetfile(varsetfile_LOF, LOF_df, units = genes)
  compare_varsetfile(varsetfile_CADD, CADD_df, units = genes)

  # write to file
  output <- withr::local_tempfile()
  buildVarSet(
    anno,
    unitName = "gene_name",
    fields = c("Moderate"),
    output = output
  )
  moderate_df_from_file <- varSetFile(output)
  compare_varsetfile(moderate_df, moderate_df_from_file, units = genes)

  close(gdb)
})

test_that("buildVarSet intersect functionality works", {
  gdb <- create_example_gdb()
  genes <- c("CYP19A1", "FUS", "OPTN", "SOD1")
  where_clause_genes <- sprintf(
    "gene_name in (%s)",
    paste(shQuote(genes), collapse = ", ")
  )

  # create random variant subset
  vars <- getAnno(gdb, "var", fields = "VAR_id")$VAR_id
  vars <- sort(sample(vars, size = 200))
  uploadAnno(
    gdb,
    value = data.frame(VAR_id = vars),
    name = "random_vars",
    skipRemap = TRUE,
    verbose = FALSE
  )

  # build varset with intersect
  moderate_intersect <- buildVarSet(
    object = gdb,
    varSetName = "Moderate",
    unitTable = "varInfo",
    unitName = "gene_name",
    intersect = "random_vars",
    where = paste(
      "(ModerateImpact = 1 or HighImpact = 1) and",
      where_clause_genes
    )
  )

  # verify intersection worked
  vars_check <- listVars(collapseVarSetList(moderate_intersect))
  expect_true(all(vars_check %in% vars))

  close(gdb)
})

test_that("buildVarSet-data.frame mask handling works", {
  df <- data.frame(
    VAR_id = c("1", "2", "3"),
    gene = c("A", "A", "B"),
    impact = c(1, 0, 1),
    impact2 = c(0, 1, 1),
    score = c(0.5, 0.8, 0.3)
  )

  # compare masked vs unmasked varset
  varset_masked <- buildVarSet(
    df,
    unitName = "gene",
    fields = "impact",
    mask = TRUE
  )
  varset_unmasked <- buildVarSet(
    df,
    unitName = "gene",
    fields = "impact",
    mask = FALSE
  )

  ## first varset should be different
  expect_identical(listVars(varset_masked[[1]]), "1")
  expect_identical(listVars(varset_unmasked[[1]]), c("1", "2"))

  ## second varSet should be identical
  compare_varsets(varset_masked[[2]], varset_unmasked[[2]])

  # specify multiple fields, one masked one unmasked
  varset_combined <- buildVarSet(
    df,
    unitName = "gene",
    fields = c("impact", "impact2"),
    mask = "impact"
  )
  varset_combined_masked <- getVarSet(varset_combined, varSetName = "impact")
  varset_combined_unmasked <- getVarSet(varset_combined, varSetName = "impact2")
  expect_equal(listWeights(varset_combined_masked[[1]]), 1)
  expect_equal(listWeights(varset_combined_unmasked[[1]]), c(0, 1))
  expect_equal(listWeights(varset_combined_masked[[2]]), 1)
  expect_equal(listWeights(varset_combined_unmasked[[2]]), 1)
})


test_that("buildVarSet-gdb input validation works", {
  gdb <- create_example_gdb()

  # expect error when unitTable does not exist in gdb
  expect_error(
    {
      buildVarSet(
        gdb,
        varSetName = "test",
        unitTable = "hello_world",
        unitName = "gene_name"
      )
    },
    regexp = "does not exist in the gdb"
  )

  # expect error when unitName does not exist in unitTable
  expect_error(
    {
      buildVarSet(
        gdb,
        varSetName = "test",
        unitTable = "varInfo",
        unitName = "hello_world"
      )
    },
    regexp = "is not present in"
  )

  # expect error when intersection table does not exist
  expect_error(
    {
      buildVarSet(
        gdb,
        varSetName = "test",
        unitTable = "varInfo",
        unitName = "gene_name",
        intersection = "hello, world"
      )
    },
    regexp = "Not all intersection tables present in gdb"
  )

  # expect error when output cannot be written
  tmpdir <- withr::local_tempdir()
  expect_error(
    {
      suppressWarnings(buildVarSet(
        gdb,
        varSetName = "test",
        unitTable = "varInfo",
        unitName = "gene_name",
        output = tmpdir
      ))
    },
    regexp = "Failed to open output file for writing"
  )
})

test_that("buildVarSet-data.frame input validation works", {
  df <- data.frame(
    VAR_id = c("1", "2"),
    gene = c("A", "B"),
    impact = c(1, 0)
  )

  # expect error when `VAR_id` column is missing
  expect_error(
    {
      buildVarSet(
        df %>% dplyr::select(-VAR_id),
        unitName = "gene",
        fields = "impact"
      )
    },
    regexp = "'VAR_id' column must be present"
  )

  # expect error when `unitName` does not exist in data.frame
  expect_error(
    {
      buildVarSet(df, unitName = "hello_world", fields = "impact")
    },
    regexp = "was not found"
  )

  # expect error when one of the `fields` does not exist in data.frame
  expect_error(
    {
      buildVarSet(df, unitName = "gene", fields = c("impact", "hello_world"))
    },
    regexp = "were not found"
  )

  # expect error when `mask` is of incorrect type
  expect_error(
    {
      buildVarSet(df, unitName = "gene", fields = "impact", mask = 123)
    },
    regexp = "`mask` must be"
  )
})
