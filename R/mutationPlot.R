#' Mutation Plot
#'
#' Generates a mutation plot visualizing variant-level and gene-level association results along a transcript structure.  
#' Optionally, custom tracks such as protein domains or mutation clusters can be overlaid.
#'
#' @param singleVar A `singlevarResult` or data.frame containing single variant association results.
#' Required columns include `POS` (variant position, CDS coordinates), `P` (p-value), and `OR` (odds ratio). 
#' An optional `impact` column can be included to represent variant impact, which will be mapped to different point shapes.
#' @param cds A `GRanges` or `IRanges` object representing the coding sequence (CDS) regions of the transcript.  
#' @param customTracks An optional `rvbResult` or data.frame containing rare variant association statistics for custom tracks.
#' Should include `start` and `end` columns (CDS coordinates).
#' A `trackType` column is required if `splitByTrackType = TRUE`. 
#' Additional columns can be included for hover information (see `trackHoverFields`).
#' @param rvbGene An optional data.frame or `rvbResult` object containing gene-level association results.  
#' Should contain a `P` column for the gene-level p-value.
#' @param cdsGapSize The size of the gap to introduce between CDS exons, in base pairs. Defaults to 30.
#' @param cdsLimits Specifies the y-axis limits for the CDS track (a vector of length 2).  Defaults to `c(-0.05, 0.05)` (nominal significance).
#' @param pointRange Specifies the minimum and maximum size of the points representing variants (vector of length 2). 
#' Point size corresponds to the absolute log of the odds ratio (`OR`). Defaults to `c(0.25, 3)`.
#' @param impactScale Named character vector mapping impact levels to point shapes. 
#' For example: `c("HIGH" = 24, "MODERATE" = 21)`. Only used if `singleVar` contains an `impact` column.
#' @param trackSpacing Spacing between custom tracks. Defaults to 1.
#' @param trackHeight Height of each custom track. Defaults to 1.
#' @param trackOrder Optional vector of track names, specifying the order in which the track are plotted (if `splitByTrackType = TRUE`)
#' @param splitByTrackType Should tracks be grouped by type (TRUE/FALSE)?
#' If `TRUE`, the `customTracks` data frame must contain a `trackType` column. Defaults to `FALSE`.
#' @param panelsizes Relative sizes of the custom tracks vs. the other tracks. Defaults to `c(1,3)`.
#' @param interactive Should the plot be interactive (using plotly)? Defaults to `FALSE`.
#' @param svHoverFields If `interactive = TRUE`, an optional character vector that specifies which columns from `singleVar` to include in the hover information for single variant points.
#' @param trackHoverFields If `interactive = TRUE`, an optional character vector that specifies which columns from `customTracks` to include in the hover information for custom tracks.
#' @param cdsHoverFields If `interactive = TRUE`, an optional character vector that specifies which columns from the `cds` data frame to include in hover information for CDS regions.
#' 
#' @export
mutationPlot <- function(
    singleVar,
    cds,
    customTracks = NULL,
    rvbGene = NULL,
    cdsGapSize = 30,
    cdsLimits = c(log10(0.05), -log10(0.05)), 
    pointRange = c(0.25, 3),
    impactScale = NULL,
    trackSpacing = 1,
    trackHeight = 1,
    trackOrder = NULL,
    splitByTrackType = FALSE,
    panelsizes = c(1,3),
    interactive = FALSE,
    svHoverFields = NULL,
    trackHoverFields = NULL,
    cdsHoverFields = NULL
) {
  
  ## convert rvatResults to data.frame (if not already)
  if (!is(singleVar, "singlevarResult") && !is(singleVar, "data.frame")) {
    stop("`singleVar` should be of class singlevarResult or data.frame")
  }
  singleVar <- as.data.frame(singleVar)
  if (!is.null(customTracks)) {
    if (!is(customTracks, "rvbResult") && !is(customTracks, "data.frame")) {
      stop("`customTracks` should be of class rvbResult or data.frame")
    } 
  customTracks <- as.data.frame(customTracks)
  }
  if (!is.null(rvbGene)) {
    if (!is(rvbGene, "rvbResult") && !is(rvbGene, "data.frame")) {
      stop("`rvbGene` should be of class rvbResult or data.frame")
    }
    if (!is.null(rvbGene)) rvbGene <- as.data.frame(rvbGene)
  }
  
  # cds
  if (!is(cds, "GRanges") && !is(cds, "IRanges")) {
    stop("`cds` should be of class GRanges or IRanges")
  }
  cds <- sort(cds)
  mcols(cds)$gap <- cumsum(c(cdsGapSize, rep(cdsGapSize, (length(cds) - 1))))
  cds_remapped <- cds
  GenomicRanges::end(cds_remapped) <- GenomicRanges::end(cds_remapped) + mcols(cds_remapped)$gap
  GenomicRanges::start(cds_remapped) <- GenomicRanges::start(cds_remapped) + mcols(cds_remapped)$gap
  cds_remapped <- as.data.frame(cds_remapped)
  cds_remapped$type <- "single variants"
  if (!is.null(rvbGene)) {
    cds_remapped$P <- rvbGene$P
  }

  ### introduce gaps in custom tracks
  if (!is.null(customTracks)) {
    customTracks_remapped <- mapCDSGapsCustomTracks(customTracks, target_ranges = cds)
    customTracks_remapped$type <- "tracks"
  }
  
  ### introduce gaps in single variants
  singleVar <- mapCDSGapsSV(singleVar, target_ranges = cds)
  singleVar$type = "single variants"
  
  ## determine track order (if specified)
  if (!is.null(customTracks)) {
    if(!is.null(trackOrder)) {
      if(!"trackType" %in% colnames(customTracks_remapped)) {stop("In order to set `trackOrder`, a field named `trackType` should be included in customTracks")}
      factors <- list()
      for(track in trackOrder) {
        factors[[track]] <- customTracks_remapped[customTracks_remapped[["trackType"]] == track,]$unit
      }
      factors <- rev(unique(unlist(factors)))
      trackOrder <- unique(c("single variants", trackOrder))
    } else {
      factors <- rev(unique(customTracks_remapped$unit))
      if(splitByTrackType) {
        trackOrder <- c("single variants", unique(customTracks_remapped$trackType))
      } else {
        trackOrder <- c("single variants", "tracks")
      }
    }
  } else {
    trackOrder <- "single variants"
  }
  
  
  ## generate y-coordinates
  ### order this by trackOrder (if specified)
  if (!is.null(customTracks)) {
    coords <- dplyr::tibble(
      unit = factors,
      ymin = cumsum(c(trackSpacing, rep((trackSpacing + trackHeight), times = (dplyr::n_distinct(factors) - 1)))),
      ymax = cumsum(c((trackSpacing + trackHeight), rep((trackSpacing + trackHeight), times = (dplyr::n_distinct(factors) - 1))))
    )
    customTracks_remapped <- customTracks_remapped %>%
      dplyr::left_join(coords, by = "unit")
  }
  
  
  ## if interactive, add text 
  if(!is.null(svHoverFields)) {
    if(is.null(names(svHoverFields))) names(svHoverFields) <- svHoverFields
    singleVar <- addText(singleVar, svHoverFields)
  } else {
    singleVar$text <- ""
  }
  
  if (!is.null(customTracks)) {
    if(!is.null(trackHoverFields)) {
      if(is.null(names(trackHoverFields))) names(trackHoverFields) <- trackHoverFields
      customTracks_remapped <- addText(customTracks_remapped, trackHoverFields)
    } else {
      customTracks_remapped$text <- ""
    }
  }
  
  if(!is.null(cdsHoverFields)) {
    if(is.null(names(cdsHoverFields))) names(cdsHoverFields) <- trackHoverFields
    cds_remapped <- addText(cds_remapped, cdsHoverFields)
  } else {
    cds_remapped$text <- ""
  }
  
  ## splitByTrackType
  if(!is.null(customTracks) && splitByTrackType) {
    coords_list <- list()
    limits_list <- list()
    range_list <- list()
    for(track in unique(customTracks_remapped$trackType)) {
      units <- unique(customTracks_remapped[customTracks_remapped$trackType == track,]$unit)
      coords_list[[track]] <- coords[coords$unit %in% units,,drop=FALSE]
      range_list[[track]] <- dplyr::tibble(
        track = track, 
        min = min(coords[coords$unit %in% units,,drop = FALSE]$ymin),
        max = max(coords[coords$unit %in% units,,drop = FALSE]$ymax),
        range_total = (max - min))
    }
    range_list <- dplyr::bind_rows(range_list)
    range_total <- sum(range_list$range_total)
    
    ## determine padding based on largest range
    padding <- max(range_list$range_total)*0.05
    for(track in unique(customTracks_remapped$trackType)) {
      limits_list[[track]] <- c((range_list[range_list$track == track,]$min - padding),
                                (range_list[range_list$track == track,]$max + padding))
    }
    
    ## set panelsizes
    n <- padding*2*dplyr::n_distinct(customTracks_remapped$trackType) + range_total
    lengths <- range_list[match(trackOrder,range_list$track),]
    lengths <- lengths[!is.na(lengths$track),]$range_total
    panelsizes = c(panelsizes[1],
                   (((padding*2) + lengths)/n)*panelsizes[2]
    )
    customTracks_remapped$type = customTracks_remapped$trackType
  } else {
    panelsizes <- 1
  }
  
  ## build scales (formula style used in ggh4x)
  if (!is.null(customTracks)) {
    scales <- list()
    if(splitByTrackType) {
      for(track in unique(customTracks_remapped$trackType)) {
        scales[[track]] <- as.formula(sprintf('type == "%s" ~ ggplot2::scale_y_continuous(breaks=rowMeans(coords_list[["%s"]][,c("ymin", "ymax")]),labels=as.character(coords_list[["%s"]]$unit),limits=limits_list[["%s"]])',track,track,track,track ))
      }
    } else {
      scales[["tracks"]] <- type == "tracks" ~ ggplot2::scale_y_continuous(breaks = rowMeans(coords[,c("ymin","ymax")]),
                                                                  labels = as.character(coords$unit)
      )
    }
    customTracks_remapped$type <- factor(customTracks_remapped$type, levels = trackOrder)
  }
  
  
  ## set types
  cds_remapped$type <- factor(cds_remapped$type, levels = trackOrder)
  singleVar$type <- factor(singleVar$type, levels = trackOrder)
  
  sv_positive <- singleVar %>% dplyr::filter(OR > 1)
  sv_negative <- singleVar %>% dplyr::filter(OR < 1)
  sv_na <- singleVar %>% dplyr::filter(is.na(OR))
  
  ## flag if 'impact' field is included
  impactFlag <- if ("impact" %in% colnames(singleVar)) TRUE else FALSE
  
  ## plot ---------
  cdsplot <- ggplot2::ggplot(
    data = dplyr::tibble(type = factor(trackOrder, level = trackOrder)),
    mapping = ggplot2::aes(fill = -log10(P))
  ) +
    
    ## plot tracks 
    {if(!is.null(customTracks)) ggplot2::geom_rect(data = customTracks_remapped,
                                          mapping = ggplot2::aes(xmin = start, 
                                                        xmax = end,
                                                        ymin = ymin, 
                                                        ymax = ymax,
                                                        text = text))} +
    
    ## single variant segments
    
    ### OR > 1 (upper)
    {if(nrow(sv_positive) > 0) ggplot2::geom_segment(data = sv_positive, 
                                             mapping = ggplot2::aes(x = POS, xend = POS, y = 0, yend = -log10(P)), 
                                             color = "grey", 
                                             show.legend = FALSE)} +
    
    ### OR < 1 (lower)
    {if(nrow(sv_negative) > 0) ggplot2::geom_segment(data = sv_negative, 
                                            mapping = ggplot2::aes(x = POS, xend = POS, y = 0, yend = log10(P)), 
                                            color = "grey", 
                                            show.legend = FALSE)}+
    
    ## cds track
    ### if rvbGene is specified, color by -log10(P)
    {if(!is.null(rvbGene)) ggplot2::geom_rect(data = cds_remapped,
                                    ggplot2::aes(xmin = start,
                                         xmax = end,
                                         ymin = cdsLimits[1],
                                         ymax = cdsLimits[2],
                                         fill = -log10(P),
                                         text = text))} +
    
    {if(is.null(rvbGene)) ggplot2::geom_rect(data = cds_remapped,
                                    ggplot2::aes(xmin = start,
                                        xmax = end,
                                        ymin = cdsLimits[1],
                                        ymax = cdsLimits[2],
                                        text = text), 
                                    fill = "#FFD58A", 
                                    color = "darkgray")}+
    
    
    ## single variant dots
    
    {if(nrow(sv_positive) > 0) ggplot2::geom_point(data = sv_positive, 
                                          mapping = ggplot2::aes(x = POS,
                                                        y = -log10(P),
                                                        fill = -log10(P),
                                                        size = abs(log(OR)),
                                                        text = text,
                                                        shape = if(impactFlag) impact else NULL))} +
    
    ### OR < 1
    {if(nrow(sv_negative) > 0) ggplot2::geom_point(data = if(nrow(sv_negative) > 0) sv_negative else NULL, 
                                          mapping = ggplot2::aes(x = POS,
                                                        y = log10(P),
                                                        fill = -log10(P),
                                                        size = abs(log(OR)),
                                                        text = text,
                                                        shape = if(impactFlag) impact else NULL))}+
    
    ### is.na(OR) (singletons)
    {if(nrow(sv_na) > 0) ggplot2::geom_point(data = if(nrow(sv_na) > 0) sv_na else NULL, 
                                    mapping = ggplot2::aes(x = POS,
                                                  y = 0,
                                                  text = text,
                                                  shape = if(impactFlag) impact else NULL),
                                    fill = "black",
                                    size = 1)} +
    
    ## themes/scales/labels
    {if(!is.null(customTracks)) ggplot2::facet_grid(scales = "free_y", rows = ggplot2::vars(type), space = "free_y") } + 
    
    ### set point sizes
    ggplot2::scale_size_continuous(range = pointRange) +
    
    ### set variant impact scale
    {if(!is.null(impactScale)) ggplot2::scale_shape_manual(values = impactScale)} +
    
    ### facet
    {if(!is.null(customTracks)) ggh4x::facetted_pos_scales(y = scales)} +
    ggplot2::theme_classic() + 
    ggplot2::theme(
      text = ggplot2::element_text(size = 13),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    ) +
    ggplot2::xlab("CDS position") +
    ggplot2::ylab("") +
    ggplot2::labs(fill = bquote("-log"[10]~(italic(P))),
                  size = bquote("| Log"~(OR)~"|"),
                  shape = "impact"
    ) +
    ggh4x::force_panelsizes(rows = panelsizes)
  
  ## if interactive, use ggplotly to convert ggplot to plotly
  if(interactive) return(plotly::ggplotly(cdsplot, tooltip = "text")) else cdsplot
  }


mapCDSGapsSV <- function(input, target_ranges) {
  # convert sv input to GenomicRanges
  input_ranges <- IRanges(start = input$POS, end = input$POS)
  mcols(input_ranges)$VAR_id <- input$VAR_id
  overlaps <- GenomicRanges::findOverlaps(input_ranges, target_ranges)
  
  # add gaps
  mapped <- dplyr::tibble(
    VAR_id = mcols(input_ranges[queryHits(overlaps)])$VAR_id,
    POS = GenomicRanges::start(input_ranges[queryHits(overlaps)]) + mcols(target_ranges[subjectHits(overlaps)])$gap
  )
  input <- input %>%
    dplyr::select(-POS) %>%
    dplyr::left_join(mapped,by="VAR_id")
  input
}

mapCDSGapsCustomTracks <- function(input, target_ranges, warning = FALSE) {
  ## format
  input_ranges = IRanges(start = input$start, end = input$end)
  mcols(input_ranges)$unit <- input$unit
  mcols(input_ranges)$i <- 1:nrow(input)
  input$i <- 1:nrow(input)
  
  ## remap the start coordinates
  input_ <- IRanges(start = start(input_ranges), end = start(input_ranges))
  mcols(input_)$unit <- mcols(input_ranges)$unit
  mcols(input_)$i <- mcols(input_ranges)$i
  overlaps <- findOverlaps(input_, target_ranges)
  if(length(overlaps) < nrow(input)){
    if(warning) warning("Not all input ranges map to specified targets")
    if(!warning) stop("Not all input ranges map to specified targets")
  }
  
  ## remap start coordinates
  mapped_start <- dplyr::tibble(
    unit = mcols(input_[queryHits(overlaps)])$unit,
    i = mcols(input_[queryHits(overlaps)])$i,
    start = GenomicRanges::start(input_)[queryHits(overlaps)] + mcols(target_ranges[subjectHits(overlaps)])$gap
  )
  
  ## remap end coordinates
  input_ <- IRanges(start = GenomicRanges::end(input_ranges),end = GenomicRanges::end(input_ranges))
  mcols(input_)$unit <- mcols(input_ranges)$unit
  mcols(input_)$i <- mcols(input_ranges)$i
  overlaps <- GenomicRanges::findOverlaps(input_, target_ranges)
  mapped_end <- dplyr::tibble(
    unit = mcols(input_[queryHits(overlaps)])$unit,
    i = mcols(input_[queryHits(overlaps)])$i,
    end = GenomicRanges::start(input_) + mcols(target_ranges[subjectHits(overlaps)])$gap
  )
  
  ## combined
  mapped <- mapped_start %>%
    dplyr::left_join(mapped_end,by = c("unit", "i"))
  
  ## replace start and end 
  input <- input %>%
    dplyr::select(-c("start", "end")) %>% 
    dplyr::left_join(mapped, by = c("i", "unit"))
  input
}

addText <- function(df, fields) {
  base_string <- paste(rep("%s: %s\n", times = length(fields)), collapse = "")
  text <- list()
  z <- 1
  for(k in 1:length(fields)) {
    text[[z]] <- names(fields)[k]
    text[[(z+1)]] <- df[[fields[k]]]
    z <- z+2    
  }
  df$text <- do.call(sprintf, c(fmt = base_string, text))
  df
}