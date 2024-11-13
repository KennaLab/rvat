## mutationly
# Function to make a mutation plot in the viewer (all the filtering of the data is done beforehand)

mutationly <- function(varPos, tracks, plotMode, title, rvc, POS, significance, wholeGeneP, acatP) {
  #Set values
  if (!is.null(POS)) {
    varPos$POS <- varPos[[POS]]
  }
  
  if (!is.null(plotMode)) {
    if (plotMode == "Single-variant P-values" | plotMode == "Track P-values") {
      Pvalue <- "P"
    } else if (plotMode == "Single-variant logP-values" | plotMode == "Track logP-values") {
      Pvalue <- "logP"
    }
  }
  
  if (!is.null(wholeGeneP) & !is.null(plotMode)) {
    if (plotMode == "Single-variant P-values" | plotMode == "Track P-values") {
      wholeGeneP <- wholeGeneP
    } else if (plotMode == "Single-variant logP-values" | plotMode == "Track logP-values") {
      wholeGeneP <- -log10(wholeGeneP)
    }
  } 
  
  if (!is.null(acatP) & !is.null(plotMode)) {
    if (plotMode == "Single-variant P-values" | plotMode == "Track P-values") {
      acatP <- acatP
    } else if (plotMode == "Single-variant logP-values" | plotMode == "Track logP-values") {
      acatP <- -log10(acatP)
    }
  } 
  
  #Make basic plot
  plot <- plotly::plot_ly() %>%
    plotly::add_trace(data = varPos, x = ~POS, y = 0, type = 'scatter', mode ='markers', name = "Variants", showlegend = FALSE) %>%
    plotly::layout(title = unique(varPos[[title]])) 
  
  #Add a wholeGeneP text if appropriate
  if (is.null(plotMode)) {
    plot <- plot
  } else if (!is.null(wholeGeneP) & plotMode %in% c("Track P-values", "Single-variant P-values")) {
    plot <- plot %>% 
      plotly::add_annotations(x=1, y=1, xref = "paper", yref = "paper", text = paste0("<b>Gene P-value: </b>", 
                                                                                      signif(wholeGeneP,4)), showarrow = F)
  } else if (!is.null(wholeGeneP)) {
    plot <- plot %>% 
      plotly::add_annotations(x=1, y=1, xref = "paper", yref = "paper", text = paste0("<b>Gene logP-value: </b>", 
                                                                                      signif(wholeGeneP,4)), showarrow = F)
  }
  
  #Add an acatP text if appropriate
  if (is.null(plotMode)) {
    plot <- plot
  } else if (!is.null(acatP) & plotMode %in% c("Track P-values", "Single-variant P-values")) {
    plot <- plot %>% 
      plotly::add_annotations(x=1, y=0.95, xref = "paper", yref = "paper", text = paste0("<b>ACAT P-value: </b>", 
                                                                                         signif(acatP,4)), showarrow = F)
  } else if (!is.null(acatP)) {
    plot <- plot %>% 
      plotly::add_annotations(x=1, y=0.95, xref = "paper", yref = "paper", text = paste0("<b>ACAT logP-value: </b>", 
                                                                                      signif(acatP,4)), showarrow = F)
  } 
  
  varPos <- varPos[!is.na(varPos$P),]
  
  if (nrow(tracks) == 0) {
    if (length(plotMode) == 0) {
      plot <- plot
    } else if (plotMode == "Single-variant P-values" | plotMode == "Single-variant logP-values") {
      if (!("singlevarResult" %in% sapply(rvc@rvatResults, checkClassrvatResult))) {
        showNotification("Single-variant mode has been selected, but no singlevarResult is loaded!")
        plot <- plot
      } else {
        if(!is.null(significance) & plotMode == "Single-variant P-values") {
          significance <- significance
        } else if (!is.null(significance) & plotMode == "Single-variant logP-values") {
          significance <- -log10(significance)
        } else if (plotMode == "Single-variant P-values") {
          significance <-  0.05/nrow(varPos)
        } else if (plotMode == "Single-variant logP-values") {
          significance <- -log10(0.05/nrow(varPos))
        }
        
        lines <- list()
        for (i in 1:nrow(varPos)) {
          lines[[i]] <-
            list(x0 = varPos[[i, "POS"]],
                 x1 = varPos[[i, "POS"]],
                 y0 = 0,
                 y1 = varPos[[i, Pvalue]],
                 line = list(color = "black", width = 1),
                 type = "line",
                 xref = "x",
                 yref = "y")
        }
        
        if (!is.null(wholeGeneP)) {
          hline <- list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significance, y1 = significance, 
                             line = list(color = "red")),
                        list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = wholeGeneP, y1 = wholeGeneP, 
                             line = list(color = "grey", dash="dot")))
        } else {
          hline <- list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significance, y1 = significance, 
                             line = list(color = "red")))
        }
        
        shapes <- c(hline, lines)
        plot <- plot %>% plotly::layout(shapes = shapes)
        
        if (plotMode == "Single-variant P-values") {
          plot <- plot %>% 
            plotly::add_trace(data = varPos, x = ~POS, y = ~P, type = "scatter", mode = "markers", name = "singlevar",
                              hoverinfo = 'text', text = ~paste('</br> Unit: ', VAR_id, '</br> P: ', P, '</br> Effect: ', effect))
        } else if (plotMode == "Single-variant logP-values") {
          plot <- plot %>% 
            plotly::add_trace(data = varPos, x = ~POS, y = ~logP, type = "scatter", mode = "markers", name = "singlevar",
                              hoverinfo = 'text', text = ~paste('</br> Unit: ', VAR_id, '</br> logP: ', logP, '</br> Effect: ', effect))
        }
      }
    } else if (plotMode == "Track P-values" | plotMode == "Track logP-values") {
      if (!("P" %in% colnames(tracks))) {
        showNotification("Track mode has been selected, but no tracks could be found for this unit!")
        plot <- plot
      } 
    } 
  } else {
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
    
    if (length(plotMode) == 0) {
      plot <- plot %>% plotly::layout(shapes = rectangles, 
                                      yaxis = list(
                                        ticktext = as.list(unique(tracks$track)), 
                                        tickvals = as.list(0.5-1:length(unique(tracks$track))),
                                        tickmode = "array",
                                        title = ""))
    } else if (plotMode == "Single-variant P-values" | plotMode == "Single-variant logP-values") {
      if (!("singlevarResult" %in% sapply(rvc@rvatResults, checkClassrvatResult))) {
        showNotification("Single-variant P-values has been selected, but no singlevarResult is loaded!")
        plot <- plot %>% plotly::layout(shapes = rectangles, 
                                        yaxis = list(
                                          ticktext = as.list(unique(tracks$track)), 
                                          tickvals = as.list(0.5-1:length(unique(tracks$track))),
                                          tickmode = "array",
                                          title = ""))
      } else {
        if(!is.null(significance) & plotMode == "Single-variant P-values") {
          significance <- significance
        } else if (!is.null(significance) & plotMode == "Single-variant logP-values") {
          significance <- -log10(significance)
        } else if (plotMode == "Single-variant P-values" & nrow(tracks) > 0) {
          significance <- 0.05/nrow(tracks)
        } else if (plotMode == "Single-variant logP-values" & nrow(tracks) > 0) {
          significance <- -log10(0.05/nrow(tracks))
        } 
        
        lines <- list()
        for (i in 1:nrow(varPos)) {
          lines[[i]] <-
            list(x0 = varPos[[i, "POS"]],
                 x1 = varPos[[i, "POS"]],
                 y0 = 0,
                 y1 = varPos[[i, Pvalue]],
                 line = list(color = "black", width = 1),
                 type = "line",
                 xref = "x",
                 yref = "y")
        }
        
        if (!is.null(wholeGeneP)) {
          hline <- list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significance, y1 = significance, 
                             line = list(color = "red")),
                        list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = wholeGeneP, y1 = wholeGeneP, 
                             line = list(color = "grey", dash="dot")))
        } else {
          hline <- list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significance, y1 = significance, 
                             line = list(color = "red")))
        }
        
        shapes <- c(hline, rectangles, lines)
        plot <- plot %>% plotly::layout(shapes = shapes)
        
        if (plotMode == "Single-variant P-values") {
          plot <- plot %>% 
            plotly::add_trace(data = varPos, x = ~POS, y = ~P, type = "scatter", mode = "markers", name = "singlevar",
                              hoverinfo = 'text', text = ~paste('</br> Unit: ', VAR_id, '</br> P: ', P, '</br> Effect: ', effect))
        } else if (plotMode == "Single-variant logP-values") {
          plot <- plot %>% 
            plotly::add_trace(data = varPos, x = ~POS, y = ~logP, type = "scatter", mode = "markers", name = "singlevar",
                              hoverinfo = 'text', text = ~paste('</br> Unit: ', VAR_id, '</br> logP: ', logP, '</br> Effect: ', effect))
        }
      }
    } else if (plotMode == "Track P-values" | plotMode == "Track logP-values") {
      if (!("P" %in% colnames(tracks))) {
        showNotification("Track P-values has been selected, but no P-values have been loaded!")
        plot <- plot %>% plotly::layout(shapes = rectangles, 
                                        yaxis = list(
                                          ticktext = as.list(unique(tracks$track)), 
                                          tickvals = as.list(0.5-1:length(unique(tracks$track))),
                                          tickmode = "array",
                                          title = ""))
      } else {
        if(!is.null(significance) & plotMode == "Track P-values") {
          significance <- significance
        } else if (!is.null(significance) & plotMode == "Track logP-values") {
          significance <- -log10(significance)
        } else if (plotMode == "Track P-values" & nrow(tracks) > 0) {
          significance <- 0.05/nrow(tracks)
        } else if (plotMode == "Track logP-values" & nrow(tracks) > 0) {
          significance <- -log10(0.05/nrow(tracks))
        } 
        
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
        
        if (!is.null(wholeGeneP)) {
          hline <- list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significance, y1 = significance, 
                             line = list(color = "red")),
                        list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = wholeGeneP, y1 = wholeGeneP, 
                             line = list(color = "grey", dash="dot")))
        } else {
          hline <- list(list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = significance, y1 = significance, 
                             line = list(color = "red")))
        }
        
        shapes <- c(hline, rectangles, lines)
        plot <- plot %>% plotly::layout(shapes = shapes)
        
        if (plotMode == "Track P-values") {
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
        } else if (plotMode == "Track logP-values") {
          for (tr in unique(tracks$track)) {
            if ('effect' %in% colnames(tracks)) {
              plot <- plot %>% 
                plotly::add_trace(data = tracks[tracks$track == tr,], x = ~posPlot, y = ~logP, type = 'scatter', mode = 'markers', name = tr,
                                  marker = list(color = unique(tracks[tracks$track == tr, "colours"])),
                                  hoverinfo = 'text', text = ~paste('</br> Unit: ', unit, '</br> logP: ', logP, '</br> Effect: ', effect))
            } else {
              plot <- plot %>% 
                plotly::add_trace(data = tracks[tracks$track == tr,], x = ~posPlot, y = ~logP, type = 'scatter', mode = 'markers', name = tr,
                                  marker = list(color = unique(tracks[tracks$track == tr, "colours"])),
                                  hoverinfo = 'text', text = ~paste('</br> Unit: ', unit, '</br> logP: ', logP))
            }
          }
        }
      }
    } 
  }
  
  plot
}