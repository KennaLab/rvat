data(GT)

# show a genoMatrix in console
test_that("show works", {
  expect_true(stringr::str_detect(capture_output({show(GT)}), "rvat genoMatrix"))
  }
)

# test the recode method using snapshots
test_that("recoding and aggregate work" ,{
  
  ## simple aggregate
  expect_snapshot_value(aggregate(recode(GT, imputeMethod = "meanImpute"), returnGT=FALSE), style = "serialize")
  expect_snapshot_value(aggregate(recode(GT, imputeMethod = "missingToRef"), returnGT=FALSE), style = "serialize")
  expect_snapshot_value(aggregate(recode(GT, imputeMethod = "meanImpute",geneticModel="dominant"),returnGT=FALSE), style = "serialize")
  expect_snapshot_value(aggregate(recode(GT, imputeMethod = "meanImpute",geneticModel="recessive"),returnGT=FALSE), style = "serialize")
  expect_snapshot_value(aggregate(recode(GT, imputeMethod = "meanImpute",MAFweights="mb"),returnGT=FALSE), style = "serialize")
  expect_snapshot_value(aggregate(recode(GT, imputeMethod = "meanImpute",MAFweights="mb", geneticModel="dominant"),returnGT=FALSE), style = "serialize")
  expect_snapshot_value(aggregate(recode(GT, imputeMethod = "meanImpute",MAFweights="mb", geneticModel="recessive"),returnGT=FALSE), style = "serialize")
})


# check getCarriers by comparing with results of aggregate
test_that("getCarriers works" ,{

  ## carriers
  carriers <- getCarriers(GT)
  carriers_agg <- carriers %>%
    dplyr::group_by(IID) %>%
    dplyr::summarize(aggregate = sum(genotype))
  agg <- aggregate(recode(GT, imputeMethod = "meanImpute"), returnGT=FALSE)
  agg <- agg[agg >= 1]
  agg <- data.frame(
    IID = names(agg),
    aggregate = round(agg),
    stringsAsFactors=FALSE
  )  %>%
    dplyr::left_join(carriers_agg,by="IID")
  expect_identical(agg$aggregate.x, agg$aggregate.y)

  ## check adding rowDataFields/colDataFields
  carriers <- getCarriers(GT, colDataFields = c("sex", "age"))
  expect_true(all(c("sex", "age") %in% colnames(carriers)))
  carriers <- getCarriers(GT, rowDataFields = c("REF", "ALT"))
  expect_true(all(c("REF", "ALT") %in% colnames(carriers)))

  ### expect error when adding non-existing fields
  expect_error({carriers <- getCarriers(GT, colDataFields = c("A"))}, regexp = "Not all specified colDataFields are present in")
  expect_error({carriers <- getCarriers(GT, rowDataFields = c("A"))}, regexp = "Not all specified rowDataFields are present in")

  ## check aggregate/groupBy parameter
  carriers <- getCarriers(GT, aggregate = TRUE, groupBy = "superPop")
  GT_ <- recode(GT, imputeMethod = "meanImpute")
  check <- unlist(lapply(unique(GT$superPop), FUN = function(x) mean(aggregate(GT_[,GT_$superPop == x])$aggregate)))
  names(check) <- unique(GT$superPop)
  carriers <- carriers %>%
    dplyr::left_join(tibble(superPop = names(check), aggregate = unname(check)), by = "superPop")
  expect_equal(carriers$meanBurdenScore, carriers$aggregate)

})


test_that("updateGT works" ,{
  # updating annotations ---
  rowdata <- data.frame(
    VAR_id = rownames(GT),
    A = sample(LETTERS, size = nrow(GT), replace = TRUE),
    B = sample(LETTERS, size = nrow(GT), replace = TRUE),
    stringsAsFactors = FALSE
  )

  ## check if the input order doesn't matter
  GT1 <- GT
  GT1 <- updateGT(GT1, anno = rowdata)
  GT2 <- GT
  GT2 <- updateGT(GT2, anno = rowdata[sample(1:nrow(GT)),])
  expect_identical(GT1, GT2)

  ## check if error is thrown if anno includes duplicate VAR_ids
  expect_error({GT1 <- updateGT(GT1, anno = rbind(rowdata, rowdata))})

  # updating phenotypes ---

  ## add one field
  GT1 <- GT
  sm <- colData(GT1)
  sm$random_column <- sample(LETTERS, size = ncol(GT1), replace = TRUE)
  GT1 <- updateGT(GT, SM = sm)
  expect_identical(sm, colData(GT1))

  ## check if order matters
  GT1 <- GT
  sm <- colData(GT1)
  sm$random_column <- sample(LETTERS, size = ncol(GT1), replace = TRUE)
  GT1 <- updateGT(GT, SM = sm[sample(1:nrow(sm), size = nrow(sm)),])
  expect_identical(sm, colData(GT1))

  ## check uploading non-existent IIDs
  GT1 <- GT
  sm <- colData(GT1)
  sm <- rbind(sm[,c("IID", "sex")], DataFrame(IID = c("sample1", "sample2"), sex = c(0, 0), row.names = c("sample1", "sample2")))
  expect_error({GT1 <- updateGT(GT, SM = sm)})

  ## check uploading a subset
  GT1 <- GT
  sm <- colData(GT1)[2000:3000,]
  expect_error({GT1 <- updateGT(GT, SM = sm)})

})

# compare summariseGeno/getMAF etc. with plink
test_that("summariseGeno/getMAF/getNCarriers/getAC work" ,{
  gdb <- gdb(rvat_example("rvatData.gdb"))
  varids <- getAnno(gdb, "var", fields = "VAR_id")$VAR_id
  GT <- getGT(
    gdb,
    VAR_id = varids,
    cohort = "pheno",
    verbose = FALSE
  )

  # full sample size

  ## load plink count summary and format
  plink_freq <- readr::read_table("../data/rvatData.frq.gz", show_col_types = FALSE)
  plink_counts <- readr::read_table("../data/rvatData.frq.counts.gz", show_col_types = FALSE)
  plink_hwe <- suppressWarnings(readr::read_table("../data/rvatData.hwe.gz", show_col_types = FALSE))
  plink_cr <- readr::read_table("../data/rvatData.lmiss", show_col_types = FALSE)
  plink <- plink_freq %>%
    dplyr::select(SNP, MAF) %>%
    dplyr::left_join(plink_counts %>% dplyr::select(SNP, C1, C2), by = "SNP") %>%
    dplyr::left_join(plink_hwe %>% dplyr::select(SNP, P), by = "SNP") %>%
    dplyr::left_join(plink_cr %>% dplyr::select(SNP,F_MISS), by = "SNP") %>%
    dplyr::mutate(
      SNP = stringr::str_replace(SNP, "^25:", "X:"),
      SNP = stringr::str_replace(SNP, "^23:", "X:"),
    )

  ## RVAT summary
  var <- getAnno(gdb, "var", fields = c("VAR_id", "CHROM", "POS", "REF", "ALT"))
  var$VAR_id <- as.character(var$VAR_id)
  var$ID = sprintf("%s:%s:%s:%s", stringr::str_replace(var$CHROM, "chr", ""), var$POS, var$REF, var$ALT)
  sumgeno <- summariseGeno(GT)
  var <- var[match(rownames(GT), as.character(var$VAR_id)),]
  var$MAF <- getMAF(GT)
  var$AF <- getAF(GT)
  var$AC <- getAC(GT)
  var$CR <- getCR(GT)
  var$ncarriers <- getNCarriers(GT)
  var <- var %>%
    dplyr::left_join(sumgeno[,c("VAR_id", "geno0", "geno1", "geno2", "hweP")], by = "VAR_id") %>%
    dplyr::left_join(plink, by = c("ID" = "SNP"))

  ## compare
  expect_equal(var$MAF.x, var$MAF.y, tolerance = 1e-3)
  expect_equal(var$AC, var$C1)
  expect_equal(var$CR, (1-var$F_MISS), tolerance = 1e-3)
  expect_equal(sumgeno$AF, var$AF)
  var <- var %>%
    dplyr::mutate(check_carriers = AC - geno2)
  expect_equal(var$ncarriers, var$check_carriers)

  ## also check sample callrates
  plink_cr <- readr::read_table("../data/rvatData.imiss", show_col_types = FALSE)
  rvat_cr <- getCR(GT, var = FALSE)
  expect_equal(unname(rvat_cr), (1-plink_cr$F_MISS), tolerance = 1e-3)

  # compare for a random subset
  plink_freq <- readr::read_table("../data/rvatData.randomsubset.frq.gz", show_col_types = FALSE)
  plink_counts <- readr::read_table("../data/rvatData.randomsubset.frq.counts.gz", show_col_types = FALSE)
  plink_hwe <- suppressWarnings(readr::read_table("../data/rvatData.randomsubset.hwe.gz", show_col_types = FALSE))
  plink_cr <- readr::read_table("../data/rvatData.randomsubset.lmiss", show_col_types = FALSE)
  plink <- plink_freq %>%
    dplyr::select(SNP, MAF) %>%
    dplyr::left_join(plink_counts %>% dplyr::select(SNP, C1, C2), by = "SNP") %>%
    dplyr::left_join(plink_hwe %>% dplyr::select(SNP, P), by = "SNP") %>%
    dplyr::left_join(plink_cr %>% dplyr::select(SNP,F_MISS), by = "SNP") %>%
    dplyr::mutate(
      SNP = stringr::str_replace(SNP, "^25:", "X:"),
      SNP = stringr::str_replace(SNP, "^23:", "X:"),
    )

  var <- getAnno(gdb, "var", fields = c("VAR_id", "CHROM", "POS", "REF", "ALT"))
  var$VAR_id <- as.character(var$VAR_id)
  var$ID = sprintf("%s:%s:%s:%s", stringr::str_replace(var$CHROM, "chr", ""), var$POS, var$REF, var$ALT)
  keep <- readr::read_table("../data/rvatData.random.keep", col_names = FALSE, show_col_types = FALSE)
  GT_subset <- GT[,keep$X1]
  sumgeno <- summariseGeno(GT_subset)
  var <- var[match(rownames(GT_subset), as.character(var$VAR_id)),]
  var$MAF <- getMAF(GT_subset)
  var$AF <- getAF(GT_subset)
  var$AC <- getAC(GT_subset)
  var$CR <- getCR(GT_subset)
  var$ncarriers <- getNCarriers(GT_subset)
  var <- var %>%
    dplyr::left_join(sumgeno[,c("VAR_id", "geno0", "geno1", "geno2", "hweP")], by = "VAR_id") %>%
    dplyr::left_join(plink, by = c("ID" = "SNP"))

  ## compare
  expect_equal(var$MAF.x, var$MAF.y, tolerance = 1e-3)
  expect_equal(var$AC, var$C1)
  expect_equal(var$CR, (1-var$F_MISS), tolerance = 1e-3)
  expect_equal(sumgeno$AF, var$AF)
  var <- var %>%
    dplyr::mutate(check_carriers = AC - geno2)
  expect_equal(var$ncarriers, var$check_carriers)

  ## also check sample callrates
  plink_cr <- readr::read_table("../data/rvatData.randomsubset.imiss", show_col_types = FALSE)
  plink_cr <- plink_cr[match(colnames(GT_subset), plink_cr$IID),]
  rvat_cr <- getCR(GT_subset, var = FALSE)
  expect_equal(unname(rvat_cr), (1-plink_cr$F_MISS), tolerance = 1e-3)
})

# check flipToMinor
test_that("flipToMinor works" ,{
  gdb <- gdb(rvat_example("rvatData.gdb"))
  var <- getAnno(gdb, "varInfo", where = "gene_name = 'ABCA4'")
  GT1 <- getGT(
    gdb,
    VAR_id = var$VAR_id,
    anno = "var",
    verbose = FALSE
  )
  expect_equal(getMAC(GT1), getAC(GT1))
  GT1 <- flipToMinor(GT1)

  # manually flip some alleles
  flip <- rep(FALSE, nrow(GT1))
  flip[c(47, 52, 71)] <- rep(TRUE, 3)
  assays(GT1)$GT <- abs(assays(GT1)$GT - 2 * matrix(rep(flip, each = metadata(GT1)$m), nrow = metadata(GT1)$nvar, byrow = TRUE))
  rowData(GT1)$AF <- getAF(GT1)
  expect_failure(expect_equal(getMAC(GT1), getAC(GT1)))
  GT1_flipped <- flipToMinor(GT1)
  flip <- getAF(GT1) > 0.5
  expect_equal((1-getAF(GT1_flipped)[flip]), getAF(GT1)[flip])
  expect_equal((rowData(GT1_flipped)$effectAllele != rowData(GT1)$ALT),
               flip)
  expect_equal((rowData(GT1_flipped)$otherAllele != rowData(GT1)$REF),
               flip)

  # also test for some variants on the X-chromosome
  var <- getAnno(gdb, "varInfo", where = "gene_name = 'UBQLN2'")
  GT1 <- getGT(
    gdb,
    VAR_id = var$VAR_id,
    anno = "var",
    verbose = FALSE
  )
  expect_equal(getMAC(GT1), getAC(GT1))
  GT1 <- flipToMinor(GT1)

  # manually flip some alleles
  carriers <- getCarriers(GT1, colDataFields = "sex")
  ## subset a couple that include both male and female carriers
  vars <- carriers %>%
    dplyr::group_by(VAR_id) %>%
    dplyr::filter(all(c(1, 2) %in% sex)) %$%
    unique(VAR_id)[1:5]
  vars <- carriers %>%
    dplyr::filter(VAR_id %in% vars)
  for (varid in unique(vars$VAR_id)) {
    vars_male <- vars[vars$VAR_id == varid & vars$sex == 1,]
    vars_female <- vars[vars$VAR_id == varid & vars$sex == 2,]
    assays(GT1)$GT[varid,vars_male$IID] <- 0
    assays(GT1)$GT[varid,!colnames(GT1) %in% vars_male$IID] <- 1
    assays(GT1)$GT[varid,vars_female$IID] <- abs(assays(GT1)$GT[varid,vars_female$IID]) - 2
  }
  rowData(GT1)$AF <- getAF(GT1)
  expect_failure(expect_equal(getMAC(GT1), getAC(GT1)))
  GT1_flipped <- flipToMinor(GT1)
  flip <- getAF(GT1) > 0.5
  expect_equal((1 - getAF(GT1_flipped)[flip]), getAF(GT1)[flip])
  expect_equal((rowData(GT1_flipped)$effectAllele != rowData(GT1)$ALT),
               flip)
  expect_equal((rowData(GT1_flipped)$otherAllele != rowData(GT1)$REF),
               flip)

  }
)
