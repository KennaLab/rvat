## Internal functions used in rvat_viewer()

shiny_compare_plot <- function(dat1, 
                               dat2, 
                               compare_variable,
                               transformation,
                               hover_fields, 
                               display_variable,
                               show_labels, 
                               use_custom_threshold, 
                               custom_threshold) {
  
  ## Calculate inflation for both datasets and combine them
  inflation1 <- signif(median(qchisq(1 - dat1$P,1), na.rm=TRUE)/qchisq(0.5,1),3)
  inflation2 <- signif(median(qchisq(1 - dat2$P,1), na.rm=TRUE)/qchisq(0.5,1),3)
  stats_combined <- dat1 %>% dplyr::left_join(dat2, by = "unit") %>%
    dplyr::filter(!is.na(P.x), !is.na(P.y))
  
  var1 <- paste0(compare_variable, ".x")
  var2 <- paste0(compare_variable, ".y")
  
  ## Add var
  if(transformation == "-log10") {
    stats_combined[["var.x"]] <- -log10(stats_combined[[var1]])
    stats_combined[["var.y"]] <- -log10(stats_combined[[var2]])
    labelx <- sprintf("-log10(%s)", var1)
    labely <- sprintf("-log10(%s)", var2)
  } else if (transformation == "log10") {
    stats_combined[["var.x"]] <- log10(stats_combined[[var1]])
    stats_combined[["var.y"]] <- log10(stats_combined[[var2]])
    labelx <- sprintf("log10(%s)", var1)
    labely <- sprintf("log10(%s)", var2)
  } else {
    stats_combined[["var.x"]] <- stats_combined[[var1]]
    stats_combined[["var.y"]] <- stats_combined[[var2]]
    labelx <- var1
    labely <- var2
  }
  
  min <- min(c(stats_combined[["var.x"]], 
               stats_combined[["var.y"]]),
             na.rm=TRUE)
  max <- max(c(stats_combined[["var.x"]], 
               stats_combined[["var.y"]]),
             na.rm=TRUE)
  
  ## Set threshold (defines which units are labelled in the plot)
  if(length(custom_threshold) > 0 && use_custom_threshold) {
    threshold <- custom_threshold
  } else {
    threshold <- 0.05/dplyr::n_distinct(stats_combined$unit)
  }
  
  ## Set hover fields
  if(length(hover_fields) >= 1) {
    for(i in 1:(length(hover_fields))){
      if(i == 1) col <- ""
      col <- sprintf("%s%s: %s<br>%s: %s<br>", 
                     col,
                     paste0(hover_fields[i], ".x"), 
                     stats_combined[[paste0(hover_fields[i], ".x")]],
                     paste0(hover_fields[i], ".y"), 
                     stats_combined[[paste0(hover_fields[i], ".y")]]
      )
    } 
  }
  
  
  ## Base plot
  plot <- plotly::plot_ly(data = stats_combined, hoverinfo = "text", 
                          hovertext = paste0("Unit: ", stats_combined$unit, "<br>",
                                             col)) %>%
    plotly::add_segments(x = 0, xend = max+1, y = 0, yend = max+1, line = list(dash = "dash", color = 'grey'), showlegend = FALSE) %>%
    plotly::add_markers(x = ~var.x, y = ~var.y, marker = list(color='#1f77b4'), showlegend = FALSE) 
  
  ## If one or more units pass the significance threshold, label them
  stats_combined_sig <- stats_combined %>% dplyr::filter(P.x < threshold | P.y < threshold)
  if(show_labels == 1 && nrow(stats_combined_sig) > 0) {
    
    plot <- plot %>% plotly::add_annotations(x = stats_combined_sig[["var.x"]],
                                             y = stats_combined_sig[["var.y"]],
                                             text = if(display_variable == "unit") stats_combined_sig[[display_variable]] else stats_combined_sig[[paste0(display_variable, ".x")]],
                                             xref = "x",
                                             yref = "y",
                                             showarrow = TRUE,
                                             arrowhead=2,
                                             ax = 10,
                                             ay = -10) 
  }     
  ## Add inflation to the title
  plot <- plot %>%
    plotly::layout(
      title = sprintf("inflation.x = %s; inflation.y = %s", inflation1, inflation2),
     xaxis = list(title = labelx, range = c(min-0.1*min,1.1*max), autotick = FALSE, tickmode = "array"), 
     yaxis = list(title = labely, range = c(min-0.1*min,1.1*max), autotick = FALSE, tickmode = "array")
    ) 
  plotly::toWebGL(plotly::style(plot, hoveron = NULL))
  
  
}