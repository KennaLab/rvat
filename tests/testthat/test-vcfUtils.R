vcf <- rvat_example("rvatData.vcf.gz")

test_that("vcfInfo2Table works",{
  # 
  output <- withr::local_tempfile()
  suppressMessages(vcfInfo2Table(
    vcf = vcf,
    output = output
  ))
  varinfo <- readr::read_tsv(rvat_example("rvatData.varinfo"), show_col_types = FALSE)
  varinfo$ID <- ifelse(is.na(varinfo$ID), ".", varinfo$ID)
  varinfo$FILTER <- ifelse(is.na(varinfo$FILTER), ".", varinfo$FILTER)
  varinfo2 <- readr::read_tsv(output, show_col_types = FALSE)
  varinfo2$PolyPhen <- ifelse(is.na(varinfo2$PolyPhen), ".", varinfo2$PolyPhen)
  varinfo2$SIFT <- ifelse(is.na(varinfo2$SIFT), ".", varinfo2$SIFT)
  varinfo2$CADDphred <- ifelse(is.na(varinfo2$CADDphred), ".", varinfo2$CADDphred)
  expect_equal(varinfo, varinfo2)
  
  check <- capture.output((suppressMessages(vcfInfo2Table(vcf = vcf, output = "-"))))
  }
) 

