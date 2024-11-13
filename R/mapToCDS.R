#' @rdname varSet
#' @usage NULL
#' @export
setMethod("mapToCDS", 
          signature = signature(object="gdb"),
          definition=function(object,
                              gff,
                              exonPadding = 12,
                              output = NULL,
                              gene_id = NULL,
                              transcript_id = NULL,
                              biotype = NULL,
                              verbose = TRUE
                              
          ) {
            # load and filter gtf/gtf
            if(is.character(gff)) {
              gtf <- rtracklayer::import(gff)
            } else if(is(gff, "GRanges")) {
              gtf <- gff
            } else {
              stop("`gff` parameter should be either a path to a gff/gtf-file or a GRanges object")
            }
            
            ### filter based on geneid
            if (!is.null(gene_id)) {
              gtf <- gtf[!is.na(gtf$gene_id) & gtf$gene_id %in% gene_id]
            }
            
            if (!is.null(transcript_id)) {
              gtf <- gtf[!is.na(gtf$transcript_id) & gtf$transcript_id %in% transcript_id]
            }
            
            if (!is.null(biotype)) {
              gtf <- gtf[(!is.na(gtf$transcript_biotype) & gtf$transcript_biotype %in% biotype) & (!is.na(gtf$gene_biotype) & gtf$gene_biotype %in% biotype)]
            }
            
            exons <- gtf[gtf$type == "CDS"]
            transcripts <- unique(exons$transcript_id)
            GenomeInfoDb::seqlevelsStyle(exons)="NCBI"
            if (is.null(output)) {
              results <-  vector("list", length(transcripts))
              names(results) <- transcripts
            } else {
              output <- gzcon(file(output,open='wb'))
            }
            
            chrom <- unique(as.character(seqnames(exons)))
            chrom_gdb <- RSQLite::dbGetQuery(object, "select CHROM from var_ranges")
            if(any(grepl("chr", chrom_gdb$CHROM))) {
              chrom_gdb <- gsub("chr", "", chrom_gdb$CHROM)
              chrom_check <- TRUE
            } else {
              chrom_check <- FALSE
            }
            chrom <- intersect(chrom, chrom_gdb)
            if(length(chrom)==0) {
              stop("No overlap between chromosomes included in gtf and chromosomes included in gdb.")
            }
          
            i <- 1
            for(chr in chrom) {
              if(verbose) message(sprintf("Mapping chromosome %s", chr))
              transcripts_chr <- unique(exons[seqnames(exons)==chr]$transcript_id)
              vars_chr <- unserialize(getAnno(object, "var_ranges", where = sprintf("CHROM = '%s%s'", if(chrom_check) "chr" else "", chr))$ranges[[1]])
              GenomeInfoDb::seqlevelsStyle(vars_chr) <- "NCBI"
              
              for(transcript in transcripts_chr) {
                ## subset exons 
                exons_ <- sort(exons[exons$transcript_id %in% transcript])
                strand <- unique(as.character(strand(exons_)))
                
                ## load variants based on chr,POS
                ###  check if 'chr' needs to be added
                add <- if(!is.null(exonPadding)) exonPadding else 0
                vars <- subsetByOverlaps(vars_chr, exons_, maxgap = add)
                if (length(vars) > 0) {
                  ## gaps 
                  deductionsPlus <- cumsum(width(GenomicRanges::gaps(exons_)))
                  if (length(deductionsPlus) == 1) {
                    deductionsMinus <- c(0)
                  } else if (length(deductionsPlus) > 1) {
                    deductionsMinus <- rev(c(0,cumsum(rev(width(GenomicRanges::gaps(exons_)[2:length(exons_)])))))
                  }
                  exons_$deductionsPlus <- deductionsPlus
                  exons_$deductionsMinus <- deductionsMinus
                  
                  ## map variants to exons and add gaps
                  overlaps <- GenomicRanges::findOverlaps(vars,exons_)
                  
                  if (!is.null(exonPadding)) {
                    overlaps_padded <- GenomicRanges::findOverlaps(vars,exons_,maxgap = exonPadding)
                    overlaps_padded <- overlaps_padded[!S4Vectors::queryHits(overlaps_padded) %in% S4Vectors::queryHits(overlaps)]
                  }
                  
                  if ( length(overlaps_padded) > 0 ) {
                    
                    vars_padded <- vars[S4Vectors::queryHits(overlaps_padded)]
                    vars_padded$distStart <- start(exons_[S4Vectors::subjectHits(overlaps_padded)]) - end(vars_padded) - 1
                    vars_padded$distEnd <- start(vars_padded) - end(exons_[S4Vectors::subjectHits(overlaps_padded)]) - 1 
                    vars_padded$POS_updated = ifelse(vars_padded$distStart >= 0 & vars_padded$distStart <= exonPadding,
                                                     start(vars_padded) + vars_padded$distStart + 1,
                                                     start(vars_padded) - vars_padded$distEnd - 1
                    )
                    ranges(vars_padded) <- IRanges::IRanges(
                      start = vars_padded$POS_updated,
                      end = vars_padded$POS_updated
                    )
                    vars_padded$distStart <- vars_padded$distEnd <- vars_padded$POS_updated <- NULL
                    mcols(vars_padded) <- cbind(mcols(vars_padded), mcols(exons_[S4Vectors::subjectHits(overlaps_padded)])[,c("gene_id", "transcript_id", "deductionsPlus", "deductionsMinus")])
                    vars_padded$padded = TRUE
                  }
                  
                  if ( length(overlaps) > 0 ) {
                    # generate new variant positions 
                    vars <- vars[S4Vectors::queryHits(overlaps)]
                    mcols(vars) <- cbind(mcols(vars), mcols(exons_[S4Vectors::subjectHits(overlaps)])[,c("gene_id", "transcript_id", "deductionsPlus", "deductionsMinus")])
                    vars$padded <- FALSE
                  } else {
                    vars <- vars[0,]
                  }
                  
                  if (length(overlaps) > 0 && length(overlaps_padded) > 0) {
                    vars <- sort(c(vars,vars_padded))
                  }  else if (length(overlaps_padded) > 0 && length(overlaps) == 0 ) {
                    vars <- vars_padded
                  }
                  
                  if (length(vars) > 0) {
                    if (strand == "+") {
                      mcols(vars)[["cdsPOS"]] <- start(vars) - vars$deductionsPlus
                    } else if (strand == "-") {
                      mcols(vars)[["cdsPOS"]] <- rep(max(end(exons_)) + 1, length(vars)) - start(vars) - vars$deductionsMinus
                    }
                    names(vars) <- NULL
                    vars <- as.data.frame(vars)[,c("VAR_id", "seqnames", "start","gene_id", "transcript_id", "cdsPOS", "padded")]
                    colnames(vars) <- c("VAR_id", "CHROM", "POS", "gene_id", "transcript_id", "cdsPOS", "padded")
                    if(!is.null(output)) {
                      if(i == 1) {
                        write.table(vars, file = output, append = FALSE, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
                      } else {
                        write.table(vars, file = output, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
                      }
                    } else {
                      results[[transcript]] <- vars
                    }
                  }
                }
                i <- i+1
              }
            }
            if(!is.null(output)) {
              close(output)
            } else {
              results <- do.call(rbind,results)
              rownames(results) <- NULL
              results
            }
          }
)