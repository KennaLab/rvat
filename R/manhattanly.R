# Code adapted from https://github.com/sahirbhatnagar/manhattanly
# also see: https://cran.r-project.org/web/packages/manhattanly/index.html

manhattanly <- function(x,
                        # col = colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Set1"))(nchr),
                        # col = RColorBrewer::brewer.pal(n = 9, name = "Greys"),
                        ...,
                        col = c("#969696", "#252525"),
                        point_size = 5,
                        labelChr = NULL,
                        suggestiveline = -log10(1e-5),
                        suggestiveline_color = "blue",
                        suggestiveline_width = 1,
                        genomewideline = -log10(5e-8),
                        genomewideline_color = "red",
                        genomewideline_width = 1,
                        highlight = NULL,
                        highlight_color = "#00FF00",
                        showlegend = FALSE,
                        showgrid = FALSE,
                        xlab = NULL,
                        ylab = "-log10(p)",
                        title = "Manhattan Plot") {
  
  UseMethod("manhattanly")
  
}

manhattanly.default <- function(x,
                                ...,
                                col = c("#969696", "#252525"),
                                point_size = 5,
                                labelChr = NULL,
                                suggestiveline = -log10(1e-5),
                                suggestiveline_color = "blue",
                                suggestiveline_width = 1,
                                genomewideline = -log10(5e-8),
                                genomewideline_color = "red",
                                genomewideline_width = 1,
                                highlight = NULL,
                                highlight_color = "#00FF00",
                                showlegend = FALSE,
                                showgrid = FALSE,
                                xlab = NULL,
                                ylab = "-log10(p)",
                                title = "Manhattan Plot") {
  
  mh <- manhattanr(x, ...)
  nchr <- mh$nchr
  manhattanly.manhattanr(mh,
                         col = col,
                         labelChr = labelChr,
                         point_size = point_size,
                         suggestiveline = suggestiveline,
                         suggestiveline_color = suggestiveline_color,
                         suggestiveline_width = suggestiveline_width,
                         genomewideline = genomewideline,
                         genomewideline_color = genomewideline_color,
                         genomewideline_width = genomewideline_width,
                         highlight = highlight,
                         highlight_color = highlight_color,
                         showlegend = showlegend,
                         showgrid = showgrid,
                         xlab = xlab,
                         ylab = ylab,
                         title = title)
}


manhattanly.manhattanr <- function(x,
                                   ...,
                                   col = c("#969696", "#252525"),
                                   point_size = 5,
                                   labelChr = NULL,
                                   suggestiveline = -log10(1e-5),
                                   suggestiveline_color = "blue",
                                   suggestiveline_width = 1,
                                   genomewideline = -log10(5e-8),
                                   genomewideline_color = "red",
                                   genomewideline_width = 1,
                                   highlight = NULL,
                                   highlight_color = "#00FF00",
                                   showlegend = FALSE,
                                   showgrid = FALSE,
                                   xlab = NULL,
                                   ylab = "-log10(P)",
                                   title = "Manhattan Plot") {
  
  # x <- manhattanr(gwasResults)
  # x <- manhattanr(kk, annotation1 = "ZSCORE", annotation2 = "EFFECTSIZE")
  # x <- manhattanr(kk, annotation1 = "ZSCORE")
  # x <- manhattanr(kk, annotation1 = "ZSCORE", annotation2 = "EFFECTSIZE")
  # x <- manhattanr(HapMap, snp = "SNP", gene = "GENE")
  # 
  # x$data %>% head
  # str(x$data)
  # labelChr <- NULL
  # col <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="Set1")))(22)
  # showgrid <- TRUE
  # labelChr = NULL
  # point_size = 5
  # suggestiveline = -log10(1e-5)
  # genomewideline = -log10(5e-8)
  # suggestiveline_color = "blue"
  # genomewideline_color = "red"
  # suggestiveline_width = genomewideline_width = 1;
  # highlight_color = "#00FF00"
  # highlight = c(significantSNP, x$data$SNP[1:20])
  # showlegend = TRUE
  # showgrid = TRUE
  # ylab = "-log10(p)"
  # xlab = NULL
  # title = "Manhattan Plot"
  # col = c("#969696", "#252525")
  
  #########
  
  d <- x$data
  pName <- x$pName
  snpName <- x$snpName
  geneName <- x$geneName
  annotation1Name <- x$annotation1Name
  annotation2Name <- x$annotation2Name
  labs <- x$labs
  xlabel <- x$xlabel
  ticks <- x$ticks
  nchr <- x$nchr
  
  if (!is.null(highlight) & is.na(snpName)) stop("You're trying to highlight snps, but havent provided a snp column")
  
  # Initialize plot
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  
  # If manually specifying chromosome labels, ensure a character vector
  # and number of labels matches number chrs.
  if (!is.null(labelChr)) {
    if (is.character(labelChr)) {
      if (length(labelChr)==length(labs)) {
        labs <- labelChr
      } else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    } else {
      warning("If you're trying to specify chromosome labels, labelChr must be a character vector")
    }
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Initalize plotly
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  p <- plotly::plot_ly()
  
  # Add an axis.
  if (nchr == 1) {
    #If single chromosome, ticks and labels automatic.
    p <- p %>% plotly::layout(p,
                          title = title,
                          xaxis = list(
                            title = if(!is.null(xlab)) xlab else xlabel,
                            # title = "ll",
                            showgrid = showgrid,
                            range = c(xmin, xmax)
                          ),
                          yaxis = list(
                            title = ylab)#,
                          #range = c(0,ceiling(max(d$logp)))
                          #)
    )
  } else {
    # if multiple chrs, use the ticks and labels you created above.
    p <- p %>% plotly::layout(p,
                          title = title,
                          xaxis = list(
                            title = if(!is.null(xlab)) xlab else "Chromosome",
                            # title = "ll",
                            showgrid = showgrid,
                            range = c(xmin, xmax),
                            autotick = FALSE,
                            tickmode = "array",
                            tickvals = ticks,
                            ticktext = labs,
                            ticks = "outside"
                          ),
                          yaxis = list(
                            title = ylab)#,
                          #range = c(0,ceiling(max(d$logp)))
                          #)
    )
  }
  
  # Create a vector of alternatiting colors
  col <- rep(col, max(d$CHR))
  
  # Add points to the plot
  if (nchr==1) {
    
    # paste(if (!is.na(snpName)) paste0(snpName,": ",d[[snpName]],"<br>"),
    # if (!is.na(geneName)) paste0(geneName,": ",d[[geneName]],"<br>"),
    # if (!is.na(annotation1Name)) paste0(annotation1Name,": ",d[[annotation1Name]],"<br>")
    # if (!is.na(annotation2Name)) paste0(annotation2Name,": ",d[[annotation2Name]],"<br>")
    
    TEXT <- paste(if (!is.na(snpName)) paste0(snpName,": ",d[[snpName]]),
                  if (!is.na(geneName)) paste0(geneName,": ",d[[geneName]]),
                  if (!is.na(annotation1Name)) paste0(annotation1Name,": ",d[[annotation1Name]]),
                  if (!is.na(annotation2Name)) paste0(annotation2Name,": ",d[[annotation2Name]]), sep = "<br>")
    
    if (is.na(snpName) && is.na(geneName) && is.na(annotation1Name) && is.na(annotation2Name)) {
      p <- p %>% plotly::add_trace(x = d$pos, y = d$logp,
                               customdata = d$unit,
                               type = "scatter",
                               mode = "markers",
                               # text = TEXT,
                               showlegend = showlegend,
                               marker = list(color = col[1],
                                             size = point_size),
                               name = paste0("chr", unique(d$CHR))) 
    } else {
      
      p <- p %>% plotly::add_trace(x = d$pos, y = d$logp,
                               customdata = d$unit,
                               type = "scatter",
                               mode = "markers",
                               text = TEXT,
                               showlegend = showlegend,
                               marker = list(color = col[1],
                                             size = point_size),
                               name = paste0("chr", unique(d$CHR)))         
    }
    
  } else {
    
    icol <- 1
    
    for(i in unique(d$index)) {
      
      tmp <- d[d$index == unique(d$index)[i], ]
      
      TEXT <- paste(if (!is.na(snpName)) paste0(snpName,": ", tmp[[snpName]]),
                    if (!is.na(geneName)) paste0(geneName,": ", tmp[[geneName]]),
                    if (!is.na(annotation1Name)) paste0(annotation1Name,": ", tmp[[annotation1Name]]),
                    if (!is.na(annotation2Name)) paste0(annotation2Name,": ", tmp[[annotation2Name]]),
                    sep = "<br>")
      
      # get chromosome name for labeling
      chromo <- unique(tmp[which(tmp$index==i),"CHR"])
      
      if (is.na(snpName) && is.na(geneName) && is.na(annotation1Name) && is.na(annotation2Name)) {
        p <- p %>% plotly::add_trace(x = tmp$pos, y = tmp$logp, 
                                 customdata = tmp$unit,
                                 type = "scatter",
                                 mode = "markers", 
                                 showlegend = showlegend,
                                 marker = list(color = col[icol],
                                               size = point_size),
                                 name = paste0("chr",chromo)) 
      } else {
        
        p <- p %>% plotly::add_trace(x = tmp$pos, y = tmp$logp, 
                                 customdata = tmp$unit,
                                 type = "scatter",
                                 mode = "markers", 
                                 showlegend = showlegend,
                                 text = TEXT,
                                 marker = list(color = col[icol],
                                               size = point_size),
                                 name = paste0("chr",chromo))        
      }
      
      icol = icol + 1
    }
    
  }
  
  if (suggestiveline & genomewideline) {p <- p %>% plotly::layout(p,
                                                              shapes = list(
                                                                list(type = "line",
                                                                     fillcolor = suggestiveline_color,
                                                                     line = list(color = suggestiveline_color,
                                                                                 width = suggestiveline_width),
                                                                     x0 = xmin, x1 = xmax, xref = "x",
                                                                     y0 = suggestiveline, y1 = suggestiveline, yref = "y"),
                                                                list(type = "line",
                                                                     fillcolor = genomewideline_color,
                                                                     line = list(color = genomewideline_color,
                                                                                 width = genomewideline_width),
                                                                     x0 = xmin, x1 = xmax, xref = "x",
                                                                     y0 = genomewideline, y1 = genomewideline, yref = "y")
                                                              ))}
  
  if (suggestiveline & !(genomewideline)) {p <- p %>% plotly::layout(p,
                                                                 shapes = list(
                                                                   list(type = "line",
                                                                        fillcolor = suggestiveline_color,
                                                                        line = list(color = suggestiveline_color,
                                                                                    width = suggestiveline_width),
                                                                        x0 = xmin, x1 = xmax, xref = "x",
                                                                        y0 = suggestiveline, y1 = suggestiveline, yref = "y")
                                                                 ))}
  
  if (!(suggestiveline) & genomewideline) {p <- p %>% plotly::layout(p,
                                                                 shapes = list(
                                                                   list(type = "line",
                                                                        fillcolor = genomewideline_color,
                                                                        line = list(color = genomewideline_color,
                                                                                    width = genomewideline_width),
                                                                        x0 = xmin, x1 = xmax, xref = "x",
                                                                        y0 = genomewideline, y1 = genomewideline, yref = "y")
                                                                 ))}
  
  # Highlight snps from a character vector
  if (!is.na(snpName)) {
    if (!is.null(highlight)) {
      if (any(!(highlight %in% d[[snpName]]))) warning("You're trying to highlight SNPs that don't exist in your results.")
      
      d.highlight <- d[which(d[[snpName]] %in% highlight), ]
      
      
      # Add points to the plot
      if (nchr==1) {
        
        TEXT <- paste(if (!is.na(snpName)) paste0(snpName,": ",d.highlight[[snpName]]),
                      if (!is.na(geneName)) paste0(geneName,": ",d.highlight[[geneName]]),
                      if (!is.na(annotation1Name)) paste0(annotation1Name,": ",d.highlight[[annotation1Name]]),
                      if (!is.na(annotation2Name)) paste0(annotation2Name,": ",d.highlight[[annotation2Name]]), sep = "<br>")
        
        p <- p %>% plotly::add_trace(x = d$pos, y = d$logp,
                                 customdata = d$unit,
                                 type = "scatter",
                                 mode = "markers",
                                 text = TEXT,
                                 showlegend = showlegend,
                                 marker = list(color = highlight_color,
                                               size = point_size),
                                 name = "of interest")
        
      } else {
        
        # icol <- 1
        
        for(i in unique(d.highlight$index)) {
          
          tmp <- d.highlight[d.highlight$index == i, ]
          
          TEXT <- paste(if (!is.na(snpName)) paste0(snpName,": ", tmp[[snpName]]),
                        if (!is.na(geneName)) paste0(geneName,": ", tmp[[geneName]]),
                        if (!is.na(annotation1Name)) paste0(annotation1Name,": ", tmp[[annotation1Name]]),
                        if (!is.na(annotation2Name)) paste0(annotation2Name,": ", tmp[[annotation2Name]]),
                        sep = "<br>")
          
          # get chromosome name for labeling
          chromo <- unique(tmp[which(tmp$index==i),"CHR"])
          p <- p %>% plotly::add_trace(x = tmp$pos, 
                                   y = tmp$logp, 
                                   customdata = tmp$unit,
                                   type = "scatter",
                                   mode = "markers", 
                                   text = TEXT,
                                   showlegend = showlegend,
                                   marker = list(color = highlight_color,
                                                 size = point_size),
                                   name = "of interest")
          # icol = icol + 1
        }
        
      }
      
      
      # p %<>% plotly::add_trace(x = d.highlight$pos,
      #                  y = d.highlight$logp,
      #                  type = "scatter",
      #                  mode = "markers",
      #                  #evaluate = TRUE,
      #                  text = d.highlight[[snpName]],
      #                  showlegend = showlegend,
      #                  marker = list(color = highlight_color,
      #                                size = point_size),
      #                  name = "of interest")
    }
  }
  p
}

# jj <- manhattan_plotly(gwasResults, genomewideline = FALSE)
#
# jj
# str(jj)

# topHits = subset(d, P <= annotatePval)
# p %>% layout(annotations = list(x = topHits$pos[10],
#                                 y = -log10(topHits$P[10]),
#                                 text = topHits$SNP[10],
#                                 showarrow = T))

#"%ni%" <- Negate("%in%")


manhattanr <- function(x,
                       chr = "CHR",
                       bp = "BP",
                       p = "P",
                       snp,
                       gene,
                       annotation1,
                       annotation2,
                       logp = TRUE) {
  
  # NULLing out strategy
  # http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  CHR = BP = P = index = NULL
  # dotargs <- list(...)
  # print(dotargs)
  # message(paste(chr, bp, p,snp, gene, dotargs))
  
  # x = HapMap
  # chr = "CHR";
  # bp = "BP";
  # p = "P";
  # snp = "SNP"
  # browser()
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found in 'x' data.frame"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found in 'x' data.frame"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found 'x' data.frame"))
  
  ## warn if you don't have a snp column
  if (!missing(snp)) {
    if (!(snp %in% names(x))) stop(sprintf("snp argument specified as %s but this column not found in 'x' data.frame", snp))
  }
  
  if (!missing(gene)) {
    if(!(gene %in% names(x))) stop(sprintf("gene argument specified as %s but this column not found in 'x' data.frame", gene))
  }
  
  if (!missing(annotation1)) {
    if (!(annotation1 %in% names(x))) stop(sprintf("annotation1 argument specified as %s but this column not found in 'x' data.frame", annotation1))
  }
  
  if (!missing(annotation2)) {
    if (!(annotation2 %in% names(x))) stop(sprintf("annotation2 argument specified as %s but this column not found in 'x' data.frame", annotation2))
  }
  
  # if (!(gene %in% names(x))) warning(paste("No GENE column found. OK unless you're trying to annotate."))
  
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Create a new data.frame with columns called CHR, BP, and P.
  d <- data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  # str(d)
  # If the input data frame has a SNP column, add it to the new data frame
  # you're creating. Rename columns according to input
  if (!missing(snp)) {
    # using transform converts the snps to a factor variable!
    # d <- transform(d, SNP = x[[snp]])
    d[["SNP"]] <- x[[snp]]
    # str(d)
    colnames(d)[which(colnames(d) == "SNP")] <- snp
  }
  
  if (!missing(gene)) {
    # d <- transform(d, GENE = x[[gene]])
    d[["GENE"]] <- x[[gene]]
    colnames(d)[which(colnames(d) == "GENE")] <- gene
  }
  
  if (!missing(annotation1)) {
    # d <- transform(d, ANNOTATION1 = x[[annotation1]])
    d[["ANNOTATION1"]] <- x[[annotation1]]
    colnames(d)[which(colnames(d) == "ANNOTATION1")] <- annotation1
  }
  
  if (!missing(annotation2)) {
    # d <- transform(d, ANNOTATION2 = x[[annotation2]])
    d[["ANNOTATION2"]] <- x[[annotation2]]
    colnames(d)[which(colnames(d) == "ANNOTATION2")] <- annotation2
  }
  
  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  
  d$pos <- NA
  
  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  d$index <- NA
  ind <- 0
  for (i in unique(d$CHR)) {
    ind <- ind + 1
    d[d$CHR==i,]$index <- ind
  }
  
  # This section sets up positions and ticks. Ticks should be placed in the
  # middle of a chromosome. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr <- length(unique(d$CHR))
  if (nchr==1) {
    ## For a single chromosome
    d$pos <- d$BP
    ticks <- floor(length(d$pos))/2+1
    xlabel <- paste('Chromosome',unique(d$CHR),'position')
    labs <- ticks
  } else {
    ## For multiple chromosomes
    lastbase <- 0
    ticks <- NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos <- d[d$index==i, ]$BP
      } else {
        lastbase <- lastbase + utils::tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, ]$pos <- d[d$index==i, ]$BP + lastbase
      }
      # Old way: assumes SNPs evenly distributed
      # ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
      # New way: doesn't make that assumption
      ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
    }
    xlabel = 'Chromosome'
    labs <- unique(d$CHR)
  }
  
  manhattanr <- list(data = d, xlabel = xlabel, ticks = ticks, labs = labs,
                     nchr = nchr, pName = p,
                     snpName = if (missing(snp)) NA else snp,
                     geneName = if (missing(gene)) NA else gene,
                     annotation1Name = if (missing(annotation1)) NA else annotation1,
                     annotation2Name = if (missing(annotation2)) NA else annotation2)
  
  class(manhattanr) <- "manhattanr"
  
  manhattanr
  
}
