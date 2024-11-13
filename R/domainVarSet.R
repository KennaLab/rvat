#' domainVarSet
#'
#' Generate weighted variant sets for use in association testing, with partitioning by functional domains
#'
#' @param gdb gdb file path.
#' @param output Output file name (output will be gz compressed text).
#' @param varSetName Name to assign varSet grouping. This identifier column is used to allow for subsequent mergeing of multiple varSet files for coordinated analysis of multiple variant filtering/ weighting strategies)
#' @param unit_domainID Vector or file path to file with 1 column containing unit_domainIDs. It is assumed that the unit and domainID are separated with a '_' and that the first part is the unitID without other '_' in it.
#' @param unitName Column name in `unitTable` with the units of interest in it
#' @param unitTable Name of the table in the gdb with variant information of the units in `unit`. Must contain a `CHROM`, `VAR_id`, `unitName` and `positions` column (default = 'POS') 
#' @param positions Column name in `unitTable` with the coordinates of the variants
#' @param domains File path to file with domain coordinates, default column names are 'domain_ID', 'domain_start', and 'domain_end'.
#' @param domainsID Column name in `domains` with domain IDs, default = 'domain_ID'
#' @param domainsStart Column name in `domains` with domain start coordinates, default = 'domain_start' (Keep in mind that the coordinates of the domains and the variants are either both in bases or both in amino acids)
#' @param domainsEnd Column name in `domains` with domain end coordinates, default = 'domain_end' (Keep in mind that the coordinates of the domains and the variants are either both in bases or both in amino acids)
#' @export
domainVarSet=function(gdb,output,varSetName, unit_domainID, unitName, unitTable, positions = "POS", domains, domainsID = "domain_ID", domainsStart = "domain_start", domainsEnd = "domain_end") {
  gdb <- gdb(gdb)
  
  domains <- readr::read_delim(domains, col_names = TRUE)
  if (sum(c(domainsID, domainsStart, domainsEnd) %in% colnames(domains)) != 3) {
    stop("The `domains` file must contain columns with domain IDs, start positions, and end positions. Check whether the default columns `domain_ID`, `domain_start`, and `domain_end` are present, or define similar columns that are present in the `domains` file!")
  }
  
  if (length(unit_domainID) == 1) {
    if (substr(unit_domainID, nchar(unit_domainID)-2, nchar(unit_domainID)) %in% c("txt", "csv", "tsv", ".gz")) {
      unit_domainID <- as.character(read.table(unit_domainID)[,1])
    }
  }
  
  snpEff <- getAnno(gdb, unitTable)
  if(sum(c("CHROM", positions, "VAR_id", unitName) %in% colnames(snpEff)) != 4) {
    stop("The `unitTable` must contain the columns 'CHROM', `positions`, 'VAR_id', and `unitName`. Check whether the correct `positions` and `unitTable` names were given!")
  }
  
  snpEff_gr <- GenomicRanges::GRanges(seqnames = snpEff$CHROM,
                                      ranges = IRanges::IRanges(start = snpEff[[positions]], width = rep(1, nrow(snpEff))),
                                      VAR_id = snpEff$VAR_id, unit = snpEff[[unitName]])
  
  snpEff_chroms <- rbind(snpEff[,c("CHROM", unitName)]) %>% dplyr::distinct()
  colnames(snpEff_chroms) <- c("CHROM", "unit") 
    #assumption: the first part of unit_domainID is the same type of unit as unitName (e.g. both are Transcript_stable_id_version from ensembl)
  
  VAR_ids_W <- unname(sapply(unit_domainID, getVARidWeights, snpEff_chroms = snpEff_chroms, snpEff_gr = snpEff_gr, domains = domains,
                           domainsID = domainsID, domainsStart = domainsStart, domainsEnd = domainsEnd))
  
  varSetName <- rep(paste0("domains_", varSetName), length(unit_domainID))
  
  varSetsDF <- data.frame(unit = unit_domainID, VAR_id_W = VAR_ids_W, varSetName = varSetName)
  varSetsDF <- varSetsDF[varSetsDF$VAR_id_W != "",]
  varSetsDF <- varSetsDF %>% tidyr::unite(V1, unit, VAR_id_W, varSetName, sep = "|")
  
  readr::write_delim(varSetsDF, output, col_names = FALSE)
}

getVARidWeights <- function(unit_domainID, snpEff_chroms, snpEff_gr, domains, domainsID, domainsStart, domainsEnd) {
  transcript <- unlist(strsplit(unit_domainID, "_", fixed = TRUE))[1]
  
  if (transcript %in% snpEff_chroms$unit) {
    ranges <- GenomicRanges::GRanges(seqnames = unique(snpEff_chroms[snpEff_chroms$unit == transcript,"CHROM"]),
                                     ranges = IRanges::IRanges(start = domains[[domainsStart]][domains[[domainsID]] == unit_domainID],
                                                              end = domains[[domainsEnd]][domains[[domainsID]] == unit_domainID]))
    
    snpEff_gr_small <- snpEff_gr[grepl(transcript,snpEff_gr$unit, fixed = TRUE),]
    
    sub_by_overlap <- IRanges::subsetByOverlaps(snpEff_gr_small, ranges)
    if (length(sub_by_overlap) == 0) {
      return("")
    } else  {
      return(paste0(paste0(unique(sub_by_overlap$VAR_id), collapse = ","), "|", 
                    paste0(rep(1,length(unique(sub_by_overlap$VAR_id))), collapse = ",")))
    }
  } else {
    return("")
  }
}
