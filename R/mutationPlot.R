### mutationPlot ---------------------------------------------------------------

#' mutationPlot
#' 
#' Plot tracks and P-values of sub-gene units based on an \code{\link{singlevarResult}} object and `tracks` data.frame.
#' Currently implemented for ensembl gene or transcript IDs
#' @param varInfo data.frame object (or filepath) with variant info. Must include a `VAR_id` and `unit` column and a column with the positions of variants
#' @param unit the unit ID that needs to be plotted. This unit must be present in the `unit` column of `tracks` if tracks need to be plotted.
#' If this is a gene ID, the `varInfo` object must have a 'geneID' column
#' @param unitCol column in `varInfo` in which units are given, e.g.'unit', 'geneID', 'feature'. Default = 'unit'
#' @param mode Indicates which type of P-values are plotted: 'none', 'singlevar', or 'tracks'. Default = 'none'
#' If 'none', all `P` columns in `singlevar`/`tracks` will be ignored.
#' If 'singlevar', P-values from the `singlevar` object will be plotted.
#' If 'tracks', P-values from the `tracks` object will be plotted.
#' @param singlevar \code{\link{singlevarResult}} object or filepath.
#' @param tracks data.frame object (or filepath) with `start`, `end`, `track`, and `unit` column. 
#' If a `P` column is included, P-values of sub-gene units can be plotted (set `mode` to 'tracks')
#' @param Pvalue whether the P-values or -log10(P)-values are plotted if mode != 'none'. Default = 'P'; other option = 'logP'
#' @param POS name of the column with variant positions in `varInfo`. Default = 'POS'
#' @param title name of column in `varInfo` from which the title of the plot must be taken OR a string with the complete title.  
#' @param significanceGenomeWide a horizontal red line will be drawn at the genome-wide significance level to facilitate interpretation. A specific number can be given as the significance level.
#' @param significanceWholeGene a horizontal grey dotted line will be drawn at the significance level of the whole gene to facilitate interpretation. 
#' Make sure this matches with the `Pvalue` that is given (P or logP). The default is a bonferroni corrected P-value.
#' If more than one unique value is present in this column, the first is selected.

mutationPlot <- function(varInfo, unit, unitCol = "unit", mode = "none", singlevar = NULL, tracks = NULL, Pvalue = "P", 
                         POS = "POS", title = NULL, significanceGenomeWide = NULL, significanceWholeGene = NULL) {
  if (!(mode %in% c("none", "singlevar","tracks"))) {
    stop("`mode` must be 'none', 'singlevar', or 'tracks'")
  }
  
  if (mode == "singlevar" & is.null(singlevar)) {
    stop("`mode` is set to 'singlevar', but no singlevarResult object has been given")
  }
  
  if (mode == "tracks" & is.null(tracks)) {
    stop("`mode` is set to 'tracks', but no tracks object has been given")
  }
  
  #Prepare colours
  colours <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")
  
  #Prepare varInfo data.frame
  if (is.character(varInfo)) {
    varInfo <- read.table(varInfo, header = TRUE)
  } else if (is.data.frame(varInfo)) {
    varInfo <- varInfo
  } else if (is.list(varInfo)) {
    new_varInfo <- list()
    for (i in varInfo) {
      if (is.character(i) & length(i) == 1) {
        item <- read.table(i, header = TRUE)
      } else if (is.data.frame(i)) {
        item <- i
      } else {
        stop("The `varInfo` parameter should be one of the following: (a list of) `data.frame` or `character` (file names)") 
      }
      new_varInfo <- c(new_varInfo, item)
    }
    varInfo_full <- new_varInfo[[1]]
    for (i in 2:length(new_varInfo)) {
      varInfo_full <- rbind(varInfo_full, new_varInfo[[i]])
    }
    varInfo <- varInfo_full
  } else {
    stop("The `varInfo` parameter should be one of the following: (a list of) `data.frame` or `character` (file names)")
  }
  
  if (sum(c("VAR_id",unitCol, POS) %in% colnames(varInfo)) != 3) {
    stop("The `varInfo` parameter much have a `VAR_id`, unitCol and POS column Check whether the correct POS column in specified" )
  }
  
  if (!(POS %in% colnames(varInfo))) {
    stop("Choose a `POS` column that is in the given varInfo object")
  }
  
  varInfo$unitShort <- unname(sapply(varInfo[[unitCol]], function(x) {unlist(strsplit(x, "_", fixed = TRUE))[1]}))
  varInfo$unitShort <- unname(sapply(varInfo$unitShort, function(x) {unlist(strsplit(x, ".", fixed = TRUE))[1]}))
  
  varInfo <- varInfo[varInfo$unitShort == unit,, drop = FALSE]
  varInfo$VAR_id <- as.character(varInfo$VAR_id)

  #Prepare unit name
  unitShort <- unlist(strsplit(unit, "_", fixed = TRUE))[1]
  unitShort <- unlist(strsplit(unitShort, ".", fixed = TRUE))[1]
  
  #Prepare singlevar parameter
  if (is.null(singlevar)) {
    singlevar <- NULL
  } else if (mode == "singlevar"){
    if (is.character(singlevar)) {
      singlevar <- singlevarResult(singlevar)
    } else if (is(singlevar, "singlevarResult")) {
      singlevar <- singlevar
    } else if (is.data.frame(singlevar)) {
      singlevar <- singlevar
    } else if (is.list(singlevar)) {
      new_singlevar <- list()
      for (i in singlevar) {
        if (is.character(i) & length(i) == 1) {
          item <- read.table(i, header = TRUE)
        } else if (is(i, "singlevarResult")) {
          item <- i
        } else if (is.data.frame(i)) {
          item <- i
        } else {
          stop("The `singlevar` parameter should be one of the following: (a list of) `singlevarResult` or `character` (file names)") 
        }
        new_singlevar <- c(new_singlevar, item)
      }
      singlevar_full <- new_singlevar[[1]]
      for (i in 2:length(new_singlevar)) {
        singlevar_full <- rbind(singlevar_full, new_singlevar[[i]])
      }
      singlevar <- singlevar_full
    } else {
      stop("The `singlevar` parameter should be one of the following: (a list of) `singlevarResult` or `character` (file names)")
    } 
    
    if (sum(c("VAR_id", "P") %in% colnames(singlevar)) != 2) {
      stop("The `singlevar` parameter should have a `VAR_id` and `P` columns")
    }
    
    varInfo <- dplyr::left_join(varInfo, data.frame(singlevar)[,c("VAR_id", "P")], by = c("VAR_id" = "VAR_id"))
    varInfo$logP <- -log10(varInfo$P)
    varInfo <- varInfo[!is.na(varInfo$P),]
    
    if (is.null(significanceGenomeWide) & Pvalue == "P") {
      significanceGenomeWide <- 0.05/nrow(singlevar)
    } else if (is.null(significanceGenomeWide) & Pvalue == "logP") {
      significanceGenomeWide <- -log10(0.05/nrow(singlevar))
    } else {
      significanceGenomeWide <- significanceGenomeWide
    }
  }
  
  #Preparation tracks parameter
  if (is.null(tracks)) {
    tracks <- NULL
  } else if (mode == "tracks") {
    if (is.character(tracks)) {
      tracks <- read.table(tracks, header = TRUE)
    } else if (is.data.frame(tracks)) {
      tracks <- tracks
    } else if (is.list(tracks)) {
      new_tracks <- list()
      for (i in tracks) {
        if (is.character(i) & length(i) == 1) {
          item <- read.table(i, header = TRUE)
        } else if (is.data.frame(i)) {
          item <- i
        } else {
          stop("The `tracks` parameter should be one of the following: (a list of) `data.frame` or `character` (file names)") 
        }
        new_tracks <- c(new_tracks, item)
      }
      tracks_full <- new_tracks[[1]]
      for (i in 2:length(new_tracks)) {
        tracks_full <- rbind(tracks_full, new_tracks[[i]])
      }
      tracks <- tracks_full
    } else {
      stop("The `tracks` parameter should be one of the following: (a list of) `data.frame` or `character` (file names)")
    } 
  
    if (sum(c("CHROM", "start","end", "track", "unit") %in% colnames(tracks)) != 5) {
      stop("The `tracks` parameter much have a `CHROM`, `start`, `end`, `track`, and `unit` column" )
    }
    
    if (mode == "tracks" & !("P" %in% colnames(tracks))) {
      stop("The mode 'tracks' has been selected, but no P-value column in present in the given tracks object")
    }
    
    if (!is.null(tracks) & length(unique(tracks$track)) > 9) {
      stop("A maximum of 9 tracks can be plotted")
    }
  
    #Add columns necessary for plotting
    tracks$unitShort <- unname(sapply(tracks$unit, function(x) {unlist(strsplit(x, "_", fixed = TRUE))[1]}))
    tracks$unitShort <- unname(sapply(tracks$unitShort, function(x) {unlist(strsplit(x, ".", fixed = TRUE))[1]}))
  
    tracks$colours <- rep("", nrow(tracks))
    tracks$y0 <- rep(0, nrow(tracks))
    tracks$y1 <- rep(0, nrow(tracks))
    for (i in 1:length(unique(tracks$track))) {
      tracks[tracks$track == unique(tracks$track)[i],"colours"] <- colours[i]
      tracks[tracks$track == unique(tracks$track)[i],"y0"] <- (-i + 0.1)  
      tracks[tracks$track == unique(tracks$track)[i],"y1"] <- (-i + 0.9) 
    }
    tracks$posPlot <- (tracks$start + tracks$end)/2
    rownames(tracks) <- 1:nrow(tracks)
    
    if (unit != unitShort) {
      start = tracks[tracks$unit == unit, "start"]
      end = tracks[tracks$unit == unit, "end"]
      tracks <- tracks[tracks$unit == unit | ((tracks$track != track & tracks$unitShort == unitShort) &
                                                (data.table::inrange(tracks$start, start, end) | 
                                                   data.table::inrange(tracks$end, start, end))),]
    } else {
      tracks <- tracks[tracks$unitShort == unitShort,]
    }
    
    if ("P" %in% colnames(tracks)) {
      tracks$logP <- -log10(tracks$P)
    }
    
    if (nrow(tracks) == 0) {
      error("The `mode` 'tracks' was selected, but no tracks could be found. Do the unit and tracks have the same type of ID?")
    }
    
    if (is.null(significanceGenomeWide) & Pvalue == "P") {
      significanceGenomeWide <- 0.05/nrow(tracks)
    } else if (is.null(significanceGenomeWide) & Pvalue == "logP") {
      significanceGenomeWide <- -log10(0.05/nrow(tracks))
    } else {
      significanceGenomeWide <- significanceGenomeWide
    }
  }
  
  ##Make the elements of the plot
  plot <- plotly::plot_ly() %>%
    plotly::add_trace(data = varInfo, x = ~POS, y = 0, type = 'scatter', mode ='markers', name = "Variants", showlegend = FALSE)
  
  if (!is.null(title)) {
    if (title %in% colnames(varInfo)) {
      plot <- plot %>% plotly::layout(title = unique(varInfo[[title]]))
    } else {
      plot <- plot %>% plotly::layout(title = title)
    }
  }
  
  if (mode == "singlevar") {
    lines <- list()
    for (i in 1:nrow(varInfo)) {
      lines[[i]] <-
        list(x0 = varInfo[[i, "POS"]],
             x1 = varInfo[[i, "POS"]],
             y0 = 0,
             y1 = varInfo[[i, Pvalue]],
             line = list(color = "black", width = 1),
             type = "line",
             xref = "x",
             yref = "y")
    }
    
    if (!is.null(tracks)) {
      rectangles <- list()
      for (i in 1:nrow(tracks)) {
        rectangles[[i]] <- list(x0 = tracks[[i,"start"]], 
                                x1 = tracks[[i,"end"]], 
                                y0 = tracks[[i,"y0"]], 
                                y1 = tracks[[i,"y1"]],
                                type = "rect",
                                opacity = 0.1,
                                fillcolor = tracks[[i,"colours"]], 
                                line = list(color = tracks[[i,"colours"]]))
      }
    }
   
    if (is.null(significanceWholeGene)) {
      hline <- list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significanceGenomeWide, y1 = significanceGenomeWide, 
                         line = list(color = "red")))
    } else {
      hline <- list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significanceGenomeWide, y1 = significanceGenomeWide, 
                         line = list(color = "red")),
                    list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significanceWholeGene, y1 = significanceWholeGene,
                         line = list(color = "grey", dash = "dot")))
    }
    
    if (is.null(tracks)) {
      shapes <- c(hline, lines)
    } else {
      shapes <- c(hline, rectangles, lines)
    }
    plot <- plot %>% plotly::layout(shapes = shapes)
    
    if (Pvalue == "P") {
      plot <- plot %>%
        plotly::add_trace(data = varInfo, x = ~POS, y = ~P, type = "scatter", mode = "markers", name = "singlevar",
                          hoverinfo = 'text', text = ~paste('</br> Unit: ', unit, '</br> P: ', P, '</br> Effect: ', effect))
    } else if (Pvalue == "logP") {
      plot <- plot %>%
        plotly::add_trace(data = varInfo, x = ~POS, y = ~logP, type = "scatter", mode = "markers", name = "singlevar",
                          hoverinfo = 'text', text = ~paste('</br> Unit: ', unit, '</br> logP: ', logP, '</br> Effect: ', effect))
    } else {
      stop("`Pvalue` must be either 'P' or 'logP'")
    }
  } else if (mode == 'tracks') {
    lines <- list()
    for (i in 1:nrow(tracks)) {
      lines[[i]] <-
        list(x0 = tracks[[i, "posPlot"]],
             x1 = tracks[[i, "posPlot"]],
             y0 = 0,
             y1 = tracks[[i, Pvalue]],
             line = list(color = tracks[[i,"colours"]], width = 1),
             type = "line",
             xref = "x",
             yref = "y")
    }
    
    rectangles <- list()
    for (i in 1:nrow(tracks)) {
      rectangles[[i]] <- list(x0 = tracks[[i,"start"]], 
                              x1 = tracks[[i,"end"]], 
                              y0 = tracks[[i,"y0"]], 
                              y1 = tracks[[i,"y1"]],
                              type = "rect",
                              opacity = 0.1,
                              fillcolor = tracks[[i,"colours"]], 
                              line = list(color = tracks[[i,"colours"]]))
    }
    
    if (is.null(significanceWholeGene)) {
      hline <- list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significanceGenomeWide, y1 = significanceGenomeWide, 
                         line = list(color = "red")))
    } else {
      hline <- list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significanceGenomeWide, y1 = significanceGenomeWide, 
                         line = list(color = "red")),
                    list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significanceWholeGene, y1 = significanceWholeGene,
                         line = list(color = "grey", dash = "dot")))
    }
    
    shapes <- c(hline, rectangles, lines)
    plot <- plot %>% plotly::layout(shapes = shapes)
    
    if (Pvalue == 'P') {
      for (tr in unique(tracks$track)) {
        if ('effect' %in% colnames(tracks)) {
          plot <- plot %>% 
            plotly::add_trace(data = tracks[tracks$track == tr,], x = ~posPlot, y = ~P, type = 'scatter', mode = 'markers', name = tr,
                              marker = list(color = unique(tracks[tracks$track == tr, "colours"])), 
                              hoverinfo = 'text', text = ~paste('</br> Unit: ', unit, '</br> P: ', P, '</br> Effect: ', effect))
        } else {
          plot <- plot %>% 
            plotly::add_trace(data = tracks[tracks$track == tr,], x = ~posPlot, y = ~P, type = 'scatter', mode = 'markers', name = tr,
                              marker = list(color = unique(tracks[tracks$track == tr, "colours"])), 
                              hoverinfo = 'text', text = ~paste('</br> Unit: ', unit, '</br> P: ', P))
        }
      }
    } else if (Pvalue == 'logP') {
      for (tr in unique(tracks$track)) {
        if ('effect' %in% colnames(tracks)) {
          plot <- plot %>% 
            plotly::add_trace(data = tracks[tracks$track == tr,], x = ~posPlot, y = ~logP, type = 'scatter', mode = 'markers', name = tr,
                              marker = list(color = unique(tracks[tracks$track == tr, "colours"])), 
                              hoverinfo = 'text', text = ~paste('</br> Unit: ', unit,'</br> logP: ', logP, '</br> Effect: ', effect))
        } else {
          plot <- plot %>% 
            plotly::add_trace(data = tracks[tracks$track == tr,], x = ~posPlot, y = ~logP, type = 'scatter', mode = 'markers', name = tr,
                              marker = list(color = unique(tracks[tracks$track == tr, "colours"])), 
                              hoverinfo = 'text', text = ~paste('</br> Unit: ', unit,'</br> logP: ', logP))
        }
      }
    } else {
      stop("`Pvalue` must be either 'P' or 'logP'")
    }
  }
  
  plot
}