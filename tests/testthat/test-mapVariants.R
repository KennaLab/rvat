# gdb -----------------------------------
file.copy(rvatData::rvat_example("rvatData.gdb"), testgdb, overwrite = TRUE)

# check mapVariants
test_that("mapVariants works",{
  
  gdb <- gdb(testgdb)
  
  # test whether all formats (gff, bed, ranges, GenomicRanges) produce same output
  output_gtf <- withr::local_tempfile()
  output_bed <- withr::local_tempfile()
  output_ranges <- withr::local_tempfile()
  output_GRanges <- withr::local_tempfile()
  mapVariants(
    gdb, 
    gff = "../data/protein_coding_genes.gtf",
    output = output_gtf,
    verbose = FALSE
  )
  mapVariants(
    gdb, 
    bed = "../data/protein_coding_genes.bed",
    bedCols = c("gene_id", "gene_name"),
    output = output_bed,
    verbose = FALSE
  )
  mapVariants(
    gdb, 
    ranges = "../data/protein_coding_genes.ranges",
    output = output_ranges,
    verbose = FALSE
  )
  ranges <- makeGRangesFromDataFrame(readr::read_tsv("../data/protein_coding_genes.ranges", show_col_types = FALSE),keep.extra.columns = TRUE)
  mapVariants(
    gdb, 
    ranges = ranges,
    output = output_GRanges,
    verbose = FALSE
  )
  
  # compare
  output_gtf <- readr::read_tsv(output_gtf, show_col_types = FALSE)
  output_bed <- readr::read_tsv(output_bed, show_col_types = FALSE)
  expect_equal(output_gtf[,c("VAR_id", "gene_id", "gene_name")],
               output_bed[,c("VAR_id", "gene_id", "gene_name")]
  )
  rm(output_bed)
  output_ranges <- readr::read_tsv(output_ranges, show_col_types = FALSE)
  expect_equal(output_gtf[,c("VAR_id", "gene_id", "gene_name")],
               output_ranges[,c("VAR_id", "gene_id", "gene_name")]
  )
  rm(output_ranges)
  output_GRanges <- readr::read_tsv(output_GRanges, show_col_types = FALSE)
  
  # compare with vcfanno output
  expect_equal(output_gtf[,c("VAR_id", "gene_id", "gene_name")],
               output_GRanges[,c("VAR_id", "gene_id", "gene_name")])
  rm(output_GRanges)
  
  # compare with vcfanno output
  vcfanno <- readr::read_tsv("../data/rvatData.parsed.vcfAnno.gene.vcfInfo2Table", show_col_types = FALSE)
  expect_equal(output_gtf  %>% dplyr::select(VAR_id, gene_id, gene_name) %>% dplyr::arrange(VAR_id, gene_id),
               vcfanno %>% dplyr::select(VAR_id = ID, gene_id, gene_name) %>% dplyr::arrange(VAR_id, gene_id)
  )
  
  # also test w/o specifying output
  genes <- mapVariants(
    gdb, 
    gff = "../data/protein_coding_genes.gtf",
    verbose = FALSE
  )
  output_gtf$score <- as.numeric(output_gtf$score)
  output_gtf$phase <- as.numeric(output_gtf$phase)
  expect_equal(as.data.frame(output_gtf), genes)
  
  
  # expect error when:
  
  ## no gff/bed/ranges specified
  expect_error({mapVariants(
    gdb, 
    output = output_gtf,
    verbose = FALSE
  )}, "At least one of")
  
  ## multiple inputs specified
  expect_error({mapVariants(
    gdb, 
    gff = "../data/protein_coding_genes.gtf",
    bed = "../data/protein_coding_genes.bed",
    output = output_gtf,
    verbose = FALSE
  )}, "Multiple inputs specified")
  
  ## punctuation marks in uploadName
  for (mark in c(".", ",", "+", "-", " ")) {
    expect_error({mapVariants(
      gdb, 
      gff = "../data/protein_coding_genes.gtf",
      uploadName = sprintf("upload%sname", mark),
      verbose = FALSE
    )}, "Table name shouldn't contain")
  }
  
  ## protected tables
  for (name in c("dosage", "meta", "anno", "var", "var_ranges", "SM", "cohort", "tmp")) {
    expect_error({mapVariants(
      gdb, 
      gff = "../data/protein_coding_genes.gtf",
      uploadName = name,
      verbose = FALSE
    )}, "already exists as a protected table")
    }
  }
)

test_that("mapVariants and uploadAnno yield identical results for single positions",{
  gdb_tmp <- withr::local_tempfile()
  file.copy(rvat_example("rvatData.gdb"), gdb_tmp)
  gdb <- gdb(gdb_tmp)
  varInfo <- getAnno(gdb, "varInfo")
  varinfo_path <- withr::local_tempfile()
  readr::write_tsv(varInfo %>% dplyr::mutate(start=POS,end=POS) %>% dplyr::select(CHROM,POS,start,end,dplyr::everything()) %>% dplyr::select(-VAR_id),
                   file = varinfo_path)
  uploadAnno(gdb,
             name = "varInfo_uploadanno",
             value = varinfo_path,
             ignoreAlleles = TRUE,
             verbose = FALSE
  )
  suppressMessages(mapVariants(
    gdb,
    ranges = varinfo_path,
    uploadName = "varInfo_mapvariants",
    verbose = TRUE
  ))
  varInfo_uploadanno <- getAnno(gdb, "varInfo_uploadanno")
  varInfo_uploadanno$start=NULL
  varInfo_uploadanno$end=NULL
  varInfo_uploadanno$CHROM=NULL
  varInfo_mapVariants <- getAnno(gdb, "varInfo_mapvariants")
  varInfo_mapVariants$ID_tmp = sprintf("%s_%s_%s", varInfo_mapVariants$VAR_id, varInfo_mapVariants$REF, varInfo_mapVariants$ALT)
  varInfo_uploadanno$ID_tmp = sprintf("%s_%s_%s", varInfo_uploadanno$VAR_id, varInfo_uploadanno$REF, varInfo_uploadanno$ALT)
  varInfo_mapVariants <- varInfo_mapVariants[match(varInfo_uploadanno$ID_tmp, varInfo_mapVariants$ID_tmp),]
  
  rownames(varInfo_mapVariants) <- NULL
  rownames(varInfo_uploadanno) <- NULL
  varInfo_uploadanno$ID <- ifelse(varInfo_uploadanno$ID == "NA", NA_character_, varInfo_uploadanno$ID)
  varInfo_uploadanno$FILTER <- ifelse(varInfo_uploadanno$FILTER == "NA", NA_character_, varInfo_uploadanno$FILTER)
  varInfo_mapVariants$FILTER <- as.character(varInfo_mapVariants$FILTER)
  expect_equal(varInfo_uploadanno, varInfo_mapVariants)
}
)

