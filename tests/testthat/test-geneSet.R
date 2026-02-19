genesetfile_c5 <- withr::local_tempfile()
genesetfile_c5_list <- withr::local_tempfile()
data(rvbresults)

test_that("buildGeneSet from list and gmt results in identical output", {
  # build GO gene set from gmt
  expect_no_error(
    suppressMessages(buildGeneSet(
      gmtpath = test_path("data/c5.go.mf.v2023.2.Hs.symbols.gmt"),
      output = genesetfile_c5
    ))
  )

  # build GO gene set from list
  genesetfile <- geneSetFile(genesetfile_c5)
  genesetlist <- as.list(as.geneSetList(genesetfile))
  expect_no_error(
    buildGeneSet(genesetlist, output = genesetfile_c5_list, verbose = FALSE)
  )

  # compare geneSetFile build from gmt and build from list
  genesetfile_from_list <- geneSetFile(genesetfile_c5_list)
  genesetlist1 <- as.geneSetList(genesetfile)
  genesetlist2 <- as.geneSetList(genesetfile_from_list)
  expect_equal(as.list(genesetlist1), as.list(genesetlist2))
})

buildGeneSet(
  gmtpath = test_path("data/c5.go.mf.v2023.2.Hs.symbols.gmt"),
  output = genesetfile_c5,
  verbose = FALSE
)
genesetfile <- geneSetFile(genesetfile_c5)
genesetlist <- as.geneSetList(genesetfile)
geneset <- getGeneSet(genesetlist, "GOMF_ACETYLCHOLINE_BINDING")[[1]]

test_that("core geneSet methods work", {
  # check if metadata includes source info
  expect_true(stringr::str_detect(
    metadata(geneset),
    "https.*GOMF_ACETYLCHOLINE_BINDING"
  ))

  # check if length is as expected
  expect_equal(length(geneset), 11)

  # convert to data.frame and back
  expect_equal(geneset, buildGeneSet(data = as.data.frame(geneset))[[1]])

  # check if units are as excpected
  expect_equal(
    listUnits(geneset),
    c(
      "ACHE",
      "CHRFAM7A",
      "CHRM3",
      "CHRNA1",
      "CHRNA3",
      "CHRNA4",
      "CHRNA7",
      "CHRNB1",
      "CHRNB2",
      "CHRNB3",
      "CHRND"
    )
  )

  # weights should be 1
  expect_equal(unname(listWeights(geneset)), rep("1", 11))
})

test_that("core geneSetList methods work", {
  # names should be identical to listGeneSets on geneSetFile
  expect_equal(names(genesetlist), listGeneSets(genesetfile))

  # length should equal that of geneSetFile
  expect_equal(length(genesetlist), length(genesetfile))

  # show
  expect_true(stringr::str_detect(
    capture_output({
      show(genesetlist)
    }),
    "geneSetList"
  ))

  # epect metadata to include source info
  expect_true(all(stringr::str_detect(
    listMetadata(genesetlist),
    "gsea-msigdb"
  )))

  # check if genesetlist includes all genes in gmt
  path <- test_path("data/c5.go.mf.v2023.2.Hs.symbols.gmt")
  genes <- system(
    sprintf("cut -f3- %s | tr '\t' '\n' | sort | uniq", path),
    intern = TRUE
  )
  expect_equal(sort(listUnits(genesetlist)), sort(genes))

  # check `[[` subsetting: should return a geneSet object
  # with correct name
  expect_true(is(genesetlist[[500]], "geneSet"))
  expect_equal(genesetlist[[500]]@geneSetName, names(genesetlist)[500])

  # check `[` subsetting:
  ## should return a geneSet object
  expect_true(is(genesetlist[c(5, 900, 1400)], "geneSetList"))
  ## check if length is as expected
  expect_equal(length(genesetlist[c(5, 900, 1400)]), 3)
  # check if names are as expected
  expect_equal(
    names(genesetlist[c(5, 900, 1400)]),
    names(genesetlist)[c(5, 900, 1400)]
  )

  # check if sorting shuffled genesetlist results in identical objects
  expect_equal(
    sort(genesetlist[sample(
      1:length(genesetlist),
      size = length(genesetlist)
    )]),
    sort(genesetlist)
  )

  # as.data.frame -> buildGeneSet -> genesetlist should be identical to original genesetlist
  genesetlist_tmp <- buildGeneSet(data = as.data.frame(genesetlist))
  genesetlist_tmp@metadata$source = genesetlist@metadata$source
  genesetlist_tmp@metadata$creationDate = genesetlist@metadata$creationDate
  expect_equal(
    genesetlist_tmp,
    genesetlist
  )

  # as.list -> buildGeneSet -> genesetlist should be identical to original genesetlist
  expect_equal(
    as.list(buildGeneSet(data = as.list(genesetlist))),
    as.list(genesetlist)
  )

  # mapToMatrix: should return a sparse matrix with correct dimensions
  mapped_matrix <- rvat:::mapToMatrix(genesetlist, results = rvbresults)
  expect_equal(
    dim(mapped_matrix),
    c(length(unique(rvbresults$unit)), length(genesetlist))
  )
  expect_true(is(mapped_matrix, "lgCMatrix"))

  # dropUnits: should return a geneSetList without specified units
  genesetlist_dropped <- dropUnits(
    genesetlist,
    c("ZFP28", "ZNF320", "CERS4", "TADA1", "MIR154")
  )
  expect_failure(expect_equal(genesetlist, genesetlist_dropped))
  expect_true(all(
    !c("ZFP28", "ZNF320", "CERS4", "TADA1", "MIR154") %in%
      listUnits(genesetlist_dropped)
  ))
  expect_equal(
    listUnits(genesetlist_dropped),
    listUnits(genesetlist)[
      !listUnits(genesetlist) %in%
        c("ZFP28", "ZNF320", "CERS4", "TADA1", "MIR154")
    ]
  )

  # remapIDs: check if dictionary mapping works correctly
  linker <- readr::read_tsv(
    test_path("data/Homo_sapiens.GRCh38.105.gene.txt.gz"),
    show_col_types = FALSE,
    progress = FALSE
  )
  linker <- linker[, c("gene_id", "gene_name")]
  genesetlist_remapped <- suppressMessages(remapIDs(
    genesetlist,
    dict = linker[, c(2, 1)],
    targets = unique(rvbresults$unit),
    duplicate_ids = "keep_first"
  ))
  check_geneset <- getGeneSet(
    genesetlist_remapped,
    "GOMF_ACETYLCHOLINE_BINDING"
  )[[1]]
  units1 <- listUnits(geneset)
  units2 <- listUnits(check_geneset)
  expect_equal(sort(linker[linker$gene_id %in% units2, ]$gene_name), units1)

  # remapIDs: compare keeping all duplicates vs keeping first
  genesetlist_remapped_keep_first <- remapIDs(
    genesetlist,
    dict = linker[, c(2, 1)],
    duplicate_ids = "keep_first",
    verbose = FALSE
  )
  genesetlist_remapped_keep_all <- remapIDs(
    genesetlist,
    dict = linker[, c(2, 1)],
    duplicate_ids = "keep_all",
    verbose = FALSE
  )
  ## expect more units when keeping all duplicates
  expect_gt(
    length(listUnits(genesetlist_remapped_keep_all)),
    length(listUnits(genesetlist_remapped_keep_first))
  )

  ## check for specific gene with duplicates
  dups <- linker[!is.na(linker$gene_name) & linker$gene_name == "RGS5", ]
  check <- getGeneSet(genesetlist, unit = "RGS5")[[1]]
  check_geneset_keep_first <- getGeneSet(
    genesetlist_remapped_keep_first,
    geneSet = check@geneSetName
  )
  check_geneset_keep_all <- getGeneSet(
    genesetlist_remapped_keep_all,
    geneSet = check@geneSetName
  )
  expect_equal(sum(dups$gene_id %in% listUnits(check_geneset_keep_first)), 1)
  expect_equal(sum(dups$gene_id %in% listUnits(check_geneset_keep_all)), 2)

  # write and read geneSetList: resulting object should be identical to original
  genesetlist_path <- withr::local_tempfile()
  write(genesetlist, genesetlist_path)
  genesetlist_path <- geneSetFile(genesetlist_path)
  genesetlist_path@metadata$creationDate <- genesetfile@metadata$creationDate
  expect_equal(as.geneSetList(genesetlist_path), as.geneSetList(genesetfile))
})


test_that("core geneSetFile methods work", {
  ## show
  expect_true(stringr::str_detect(
    capture_output({
      show(genesetfile)
    }),
    "Path"
  ))

  ## names
  expect_equal(names(genesetlist), listGeneSets(genesetfile))

  ## listGeneSets
  pathways <- readr::read_lines(
    test_path("data/c5.go.mf.v2023.2.Hs.symbols.pathways.txt")
  )
  expect_equal(listGeneSets(genesetlist), listGeneSets(genesetfile))
  expect_equal(sort(listGeneSets(genesetfile)), sort(pathways))

  ## length
  expect_equal(length(genesetfile), 1799)

  ## as.geneSetList
  expect_equal(as.geneSetList(genesetfile), genesetlist)
})

test_that("getGeneSet works", {
  # test whether extracting from genesetlist and genesetfile results in identical output
  genesets <- sample(listGeneSets(genesetlist), size = 10)
  genesets_fromlist <- getGeneSet(
    genesetlist,
    geneSet = genesets
  )
  genesets_fromfile <- getGeneSet(
    genesetfile,
    geneSet = genesets
  )
  expect_equal(genesets_fromlist, genesets_fromfile)

  # order of provided genesets should not matter
  genesets_fromlist_shuffled <- getGeneSet(
    genesetlist,
    geneSet = sample(genesets, size = length(genesets))
  )
  genesets_fromfile_shuffled <- getGeneSet(
    genesetfile,
    geneSet = sample(genesets, size = length(genesets))
  )
  expect_equal(genesets_fromlist_shuffled, genesets_fromlist)
  expect_equal(genesets_fromfile_shuffled, genesets_fromfile)

  # extract genesets that contain certain genes (currently only implemented for class geneSetList)
  units <- sample(listUnits(genesetlist), size = 10)
  genesets_fromlist <- getGeneSet(
    genesetlist,
    unit = units
  )
  expect_equal(sum(units %in% listUnits(genesets_fromlist)), 10)
})

test_that("geneSet/geneSetList/geneSetFile input validaton works", {
  genesets <- sample(listGeneSets(genesetlist), size = 10)
  genesets_fromlist <- getGeneSet(
    genesetlist,
    geneSet = genesets
  )
  genesets_fromfile <- getGeneSet(
    genesetfile,
    geneSet = genesets
  )
  # expect error when geneSet units and weights don't match 
  expect_error({
    geneSet(
      geneSetName = "test",
      units = "A,B,C",
      w = "1,2"
    )},
    regexp = "Number of weights"
  )
  # expect error when neither `geneSet` or `unit` is specified
  expect_error({
    getGeneSet(genesetlist)
  })

  expect_error({
    getGeneSet(genesetfile)
  })

  # expect warning when subsetting non-existing genesets
  expect_warning({
    genesets_fromlist2 <- getGeneSet(
      genesetlist,
      geneSet = c("A", genesets, "B")
    )
  })
  expect_warning({
    genesets_fromfile2 <- getGeneSet(
      genesetfile,
      geneSet = c("A", genesets, "B")
    )
  })
  ## should result in equal genestlists
  expect_equal(genesets_fromlist2, genesets_fromlist)
  expect_equal(genesets_fromfile2, genesets_fromfile)
  
  # expect error when buildGeneSet input is of wrong type
  expect_error(
    {
      buildGeneSet(data = "test")
    }, regexp = "`data` should be"
  )
  # buildgeneSet: either one of data or gmtpath must be provided
  expect_error(
    {
      buildGeneSet(data = NULL, gmtpath = NULL)
    }, regexp = "Specify one of"
  )

  # buildgeneSet: either one of data or gmtpath must be provided
  expect_error(
    {
      buildGeneSet(data = list("A" = c("gene1", "gene2"), "B" = TRUE))
    }, regexp = "Each element in the list"
  )

  # buildgeneSet: data.frame should contain at least two columns
  expect_error(
    {
      buildGeneSet(data = data.frame(geneSetName = c("geneSet1", "geneSet2")))
    }, regexp = "data.frame should have at least"
  )
})
