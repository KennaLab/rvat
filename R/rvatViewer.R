# ==============================================================================
#  results viewer
# ==============================================================================

#' rvatViewer
#'
#' Interactively explore results of variant/burden testing
#'
#' @param results Can be 1) an object of class \code{\link{rvbResult}} or \code{\link{singlevarResult}}, 
#' 2) a filepath to an SQL-database, 
#' 3) a `data.frame`, or 4) a list of objects of class \code{\link{rvbResult}}, \code{\link{singlevarResult}}, or `data.frame`
#' @param CHROM if the name of the column that specifies the chromosome is not named `CHROM`, 
#' the column name can be specified here.
#' @param POS if the name of the column that specifies the chromosomal position is not named `POS`, 
#' the column name can be specified here.
#' @import shiny
#' 
#' @export
#' 
rvatViewer <- function(results, CHROM = NULL, POS = NULL) {
  
  # Prepare data
  rvc <- prepareDataShiny(results)
  
  convert_vec_logical <- c(0, 1, 0, 1)
  names(convert_vec_logical) <- c("FALSE", "TRUE", "0", "1")
  core_filter_variables <-
    c("analysis",
      "cohort",
      "varSetName",
      "name",
      "pheno",
      "covar",
      "geneticModel",
      "MAFweight",
      "test")
  
  default_analysis <- metadata(rvc)$analyses[1]
  
  ## Check which core variables have more than 1 unique value
  ## Only when length(unique(var)) > 1, a filter option is included 
  core_filter_variables_default_analysis <-
    names(metadata(rvc)[[default_analysis]][core_filter_variables][lengths(metadata(rvc)[[default_analysis]][core_filter_variables]) > 1])
  if(length(core_filter_variables_default_analysis) == 0) {
    core_filter_variables_default_analysis = "test"
  }
  
  # ----------------------------------------------------------------------------
  # UI 
  ui <- fluidPage(navbarPage(HTML(
    paste0("RVAT viewer", '<sup>', 'v', packageVersion("rvat"))
  ),
  
  # ----------------------------------------------------------------------------
  # rvb panel
  
  tabPanel(
    "rvb",
    
    # sidebar rvb panel --------------------------------------------------------
    sidebarPanel(
      width = 4,
      
      # Info selected data
      h4("Selected data:"),
      htmlOutput("tab1_analysis_info"),
      hr(),
      h4("Options:"),
      checkboxInput("use_custom_threshold", label = "Use custom significance threshold?", 
                    value = FALSE),
      conditionalPanel(condition = "input.use_custom_threshold == 1",
                       fluidRow(
                         width = 6,
                         column(
                           4,
                           numericInput(
                             "custom_threshold",
                             label = NULL,
                             value = 5e-8,
                             min = 0,
                             max = 1
                           )
                         ),
                         column(2, actionButton("set_custom_threshold", "Set"))
                       )),
      
      # Edit hover variables 
      checkboxInput("tab1_edit_hover_variables", label = "Change hovertext?"),
      uiOutput("tab1_set_hover_fields_ui"),
      
      # Display variable
      selectInput("tab1_display_variable", 
                  label = "Displayed text for significant units:",
                  choices = names(metadata(rvc)[[default_analysis]]$column_types)[metadata(rvc)[[default_analysis]]$column_types %in% c("character", "factor", "Rle")],
                  selected = "unit"
                    ),
      
      hr(),
      h4("Filters:"),
      
      ## Core variables, dependent on selected analysis which inputs are included
      uiOutput("analysis_tab1_ui"),
      
      ## Core variables, dependent on selected analysis which inputs are included
      uiOutput("core_variables_tab1_main_ui"),
      
      # Additional user selected variables
      h5("Add variables:"),
      fluidRow(width = 12,
               column(
                 8,
                 selectInput(
                   "addvariables_select_tab1_main",
                   label = NULL,
                   choices = c(metadata(rvc)[[default_analysis]][["columns"]][!metadata(rvc)[[default_analysis]][["columns"]] %in% core_filter_variables]),
                   selected = "P"
                 )
               ),
               column(
                 4, actionButton("addvariables_add_tab1_main", label = "Add")
               )),
      uiOutput("addvariables_ui_tab1_main"),
      
      # Reset additional variables
      uiOutput("addvariables_ui_tab1_main_reset")
      
    ),
    
    # rvb main panel -----------------------------------------------------------
    
    mainPanel(width = 6,
              tabsetPanel(
                tabPanel("manhattan",
                         plotly::plotlyOutput("manhattan")
                         ),
                tabPanel("qqplot",
                         plotOutput("qqplot")
                ),
                tabPanel("topTable",
                         numericInput(
                           inputId = "tab1_toptable_nrows",
                           value = 10,
                           label = "Number of rows:",
                           min = 1,
                           max = Inf
                         ),
                         checkboxInput("tab1_edit_variables_toptable", label = "Select columns to show"),
                         uiOutput("tab1_edit_variables_toptable_ui"),
                         DT::dataTableOutput("topTable")
                ),
                tabPanel("compare",
                         fluidRow(width= 12,
                                  column(4,
                                         h4("Options:"),
                                         selectInput(
                                           "tab1_compare_variable",
                                           label = "Variable to compare:",
                                           selected = "P",
                                           choices = names(metadata(rvc)[[default_analysis]]$column_types)[metadata(rvc)[[default_analysis]]$column_types %in% c("numeric", "integer")]
                                         ),
                                         selectInput(
                                           "tab1_compare_transformation",
                                           label = "Transformation",
                                           selected = "-log10",
                                           choices = c("-log10", "log10", "none")
                                         ),
                                         checkboxInput("compare_show_labels", 
                                                       label = "Show labels?", value = 1),
                                         hr(),
                                         h4("Filters:"),
                                         uiOutput("core_variables_compare_tab_ui"),
                                         hr(),
                                         h4("Add variables:"),
                                         fluidRow(width = 12,
                                                  column(
                                                    8,
                                                    selectInput(
                                                      "addvariables_select_tab1_compare",
                                                      label = NULL,
                                                      choices = c(metadata(rvc)[[default_analysis]][["columns"]][!metadata(rvc)[[default_analysis]][["columns"]] %in% core_filter_variables]),
                                                      selected = "P"
                                                    )
                                                  ),
                                                  column(
                                                    4, actionButton("addvariables_add_tab1_compare", label = "Add")
                                                  )),
                                         uiOutput("addvariables_ui_tab1_compare"),
                                         # Reset additional variables
                                         uiOutput("addvariables_ui_tab1_compare_reset")
                                         ),
                                  column(8,
                                         br(),
                                         actionButton("tab1_compare_button", label = "Compare!"),
                                         br(),
                                         br(),
                                         plotly::plotlyOutput("compare_tab1")
                                        )
                         )
                )
              )
    ))))
  
  # ----------------------------------------------------------------------------
  # server 
  
  server <- function(input, output, session) {
    
    # ----------------------------------------------------------------------------
    # reactives 
    
    ## Reactive values 
    values <- reactiveValues(
      selected = c(),
      fields = c(),
      lookup = c(),
      unit_lookup = c(),
      alias_lookup = c(),
      analysis_tab1_main = default_analysis,
      core_variables_tab1_main = core_filter_variables_default_analysis,
      tab1_hover_fields = c("effect"),
      tab1_variables_toptable = c("unit", "P", "effect", "effectSE", "effectCIlower", "effectCIupper")
      # custom_threshold = c()
    )
    
    ## Selected data main tab1
    data_tab1_main <- reactive({
      req(length(filters_tab1_main[["filters"]]) > 0)
      dat <- subSet(rvc, analysis = values[["analysis_tab1_main"]], filters = filters_tab1_main[["filters"]])
      dat
    })
    
    data_tab1_compare <- eventReactive(input[["tab1_compare_button"]],{
      dat <- subSet(rvc, analysis = values[["analysis_tab1_main"]], filters = filters_tab1_compare[["filters"]])
      dat
    })
    
    ## Filters main tab 1
    ## Populated in observer 
    filters_tab1_main <- reactiveValues(
      filters = list(),
      add_filters = list()
    )
    
    # Filters compare tab
    ## Populated in observer 
    filters_tab1_compare <- reactiveValues(
      filters = list(),
      add_filters = list()
    )
    
    # ----------------------------------------------------------------------------
    # observers 
    
    ## rvb panel ---------------------
    
    ## Update filters 
    observe({
      
      ## Extract inputs core variables
      req(values[["core_variables_tab1_main"]])
      filters_core <- lapply(
        values[["core_variables_tab1_main"]],
        FUN = function(x) {
          type <-
            rvc@metadata[[values[["analysis_tab1_main"]]]]$column_types[x]
            req(input[[sprintf("%s_tab1_main", x)]])
            return(list(variable = x, type = type, string = input[[sprintf("%s_tab1_main", x)]], negate = FALSE, keepNA = FALSE))
        }
      )
      names(filters_core) <- values[["core_variables_tab1_main"]]
      
      ## Extract inputs additional variables
      filters_add <- lapply(
        filters_tab1_main[["add_filters"]],
        FUN = function(x) {
          type <-
            rvc@metadata[[values[["analysis_tab1_main"]]]]$column_types[x]
          if (type == "character") {
            req(input[[sprintf("%s_tab1_main", x)]])
            return(list(variable = x, type = type, string = input[[sprintf("%s_tab1_main", x)]], negate = FALSE, keepNA = FALSE))
            
          } else if (type == "logical") {
            req(input[[sprintf("%s_tab1_main", x)]])
            return(list(variable = x, type = type, bool = as.logical(convert_vec_logical[input[[sprintf("%s_tab1_main", x)]]]), 
                        keepNA = FALSE))
          } else if (type %in% c("numeric", "integer")) {
            req(input[[sprintf("%s_tab1_main_min", x)]])
            req(input[[sprintf("%s_tab1_main_max", x)]])
            return(list(
              variable = x,
              type = type,
              min = as.numeric(input[[sprintf("%s_tab1_main_min", x)]]),
              max = as.numeric(input[[sprintf("%s_tab1_main_max", x)]]),
              keepNA = FALSE
            ))
          }
        }
      )
      names(filters_add) <- filters_tab1_main[["add_filters"]]
      filters_tab1_main[["filters"]] <- c(filters_core, filters_add)
      
    })
    
    observe({
      
      ## Extract inputs core variables
      req(values[["core_variables_tab1_main"]])
      filters_core <- lapply(
        values[["core_variables_tab1_main"]],
        FUN = function(x) {
          type <-
            rvc@metadata[[values[["analysis_tab1_main"]]]]$column_types[x]
          req(input[[sprintf("%s_tab1_compare", x)]])
          return(list(variable = x, type = type, string = input[[sprintf("%s_tab1_compare", x)]], negate = FALSE, keepNA = FALSE))
        }
      )
      names(filters_core) <- values[["core_variables_tab1_main"]]
      
      ## Extract inputs additional variables
      filters_add <- lapply(
        filters_tab1_compare[["add_filters"]],
        FUN = function(x) {
          type <-
            rvc@metadata[[values[["analysis_tab1_main"]]]]$column_types[x]
          if (type == "character") {
            req(input[[sprintf("%s_tab1_compare", x)]])
            return(list(variable = x, type = type, string = input[[sprintf("%s_tab1_compare", x)]], negate = FALSE, keepNA = FALSE))
            
          } else if (type == "logical") {
            req(input[[sprintf("%s_tab1_compare", x)]])
            return(list(variable = x, type = type, bool = as.logical(convert_vec_logical[input[[sprintf("%s_tab1_compare", x)]]]), 
                        keepNA = FALSE))
          } else if (type %in% c("numeric", "integer")) {
            req(input[[sprintf("%s_tab1_compare_min", x)]])
            req(input[[sprintf("%s_tab1_compare_max", x)]])
            return(list(
              variable = x,
              type = type,
              min = as.numeric(input[[sprintf("%s_tab1_compare_min", x)]]),
              max = as.numeric(input[[sprintf("%s_tab1_compare_max", x)]]),
              keepNA = FALSE
            ))
          }
        }
      )
      names(filters_add) <- filters_tab1_compare[["add_filters"]]
      filters_tab1_compare[["filters"]] <- c(filters_core, filters_add)
    })
    
    ## Reset additional variables 
    observeEvent(input[["reset_addvariables_tab1_main"]],
                 {
                   filters_tab1_main[["add_filters"]] <- c()
                   updateSelectInput(session,
                                     inputId = "addvariables_select_tab1_main",
                                     choices = c(metadata(rvc)[[values[["analysis_tab1_main"]]]][["columns"]][!metadata(rvc)[[values[["analysis_tab1_main"]]]][["columns"]] %in% core_filter_variables]),
                                     selected = "P"
                   )
                 }
    )
    
    ## Reset additional variables - compare tab
    observeEvent(input[["reset_addvariables_tab1_compare"]],
                 {
                   filters_tab1_compare[["add_filters"]] <- c()
                   updateSelectInput(session,
                                     inputId = "addvariables_select_tab1_compare",
                                     choices = c(metadata(rvc)[[values[["analysis_tab1_main"]]]][["columns"]][!metadata(rvc)[[values[["analysis_tab1_main"]]]][["columns"]] %in% core_filter_variables]),
                                     selected = "P"
                   )
                 }
    )
    
    ## Set custom significance threshold 
    observeEvent(input[["set_custom_threshold"]],
                 {
                   values[["custom_threshold"]] <- input[["custom_threshold"]]
                 })
    
    ## Add additional filter variables-  sidebar
    observeEvent(input[["addvariables_add_tab1_main"]],
                 {
                   filters_tab1_main[["add_filters"]] <- c(filters_tab1_main[["add_filters"]], input[["addvariables_select_tab1_main"]])
                   updateSelectInput(session,
                                     inputId = "addvariables_select_tab1_main",
                                     choices = metadata(rvc)[[values[["analysis_tab1_main"]]]][["columns"]][!metadata(rvc)[[default_analysis]][["columns"]] %in% c(filters_tab1_main[["add_filters"]], core_filter_variables)])
                 })
    
    ## Add additional filter variables - compare tab
    observeEvent(input[["addvariables_add_tab1_compare"]],
                 {
                   filters_tab1_compare[["add_filters"]] <- c(filters_tab1_compare[["add_filters"]], input[["addvariables_select_tab1_compare"]])
                   updateSelectInput(session,
                                     inputId = "addvariables_select_tab1_compare",
                                     choices = metadata(rvc)[[values[["analysis_tab1_main"]]]][["columns"]][!metadata(rvc)[[default_analysis]][["columns"]] %in% c(filters_tab1_compare[["add_filters"]], core_filter_variables)])
                 })
    
    # Set selected samples on click
    observeEvent(plotly::event_data("plotly_click"), {
      d <- plotly::event_data("plotly_click")
      values$selected <- d
    })
    
    # Set values[["analysis_tab1_main]]
    observeEvent(input[["analysis_tab1"]],ignoreInit = TRUE,
                 {
                   values[["analysis_tab1_main"]] <- input[["analysis_tab1"]]
                   filters_tab1_main[["add_filters"]] <- c()
                   filters_tab1_compare[["add_filters"]] <- c()
                   values[["core_variables_tab1_main"]] <-
                     names(metadata(rvc)[[values[["analysis_tab1_main"]]]][core_filter_variables][lengths(metadata(rvc)[[values[["analysis_tab1_main"]]]][core_filter_variables]) > 1])
                   values[["tab1_hover_fields"]] <-
                     values[["tab1_hover_fields"]][values[["tab1_hover_fields"]] %in% metadata(rvc)[[values[["analysis_tab1_main"]]]][["columns"]]]
                   updateSelectInput(session,
                                     inputId = "tab1_compare_variable",
                                     choices = names(metadata(rvc)[[values[["analysis_tab1_main"]]]]$column_types)[metadata(rvc)[[values[["analysis_tab1_main"]]]]$column_types %in% c("numeric", "integer")],
                                     selected = "P"
                   )
                 })
    
    # Set hover fields
    observeEvent(input[["tab1_set_hover_fields"]],
                 {
                   values[["tab1_hover_fields"]] <- input[["tab1_set_hover_fields"]]
                 }
    )
    
    # Set toptable columns
    observeEvent(input[["tab1_set_variables_toptable"]],
                 {
                   values[["tab1_variables_toptable"]] <- input[["tab1_set_variables_toptable"]]
                 }
    )
    
    # ----------------------------------------------------------------------------
    # UI generation 
    
    ## rvb panel --------------------------------------------
    
    ## analysis
    output$analysis_tab1_ui <- renderUI({
      if(length(metadata(rvc)$analyses) > 1) {
        selectInput(
          inputId = "analysis_tab1",
          label = "Analysis",
          choices = metadata(rvc)$analyses,
          selected = metadata(rvc)$analyses[1]
        )
      }
    })
    
    ## inputs additional variables
    observeEvent(values[["core_variables_tab1_main"]],
                 {
                   output$core_variables_tab1_main_ui <- renderUI({
                     lapply(
                       values[["core_variables_tab1_main"]],
                       FUN = function(x) {
                         choices <-
                           metadata(rvc)[[values[["analysis_tab1_main"]]]][[x]]
                         return(selectInput(
                           inputId = sprintf("%s_tab1_main", x),
                           label = x,
                           choices = choices,
                           selected = choices[1]
                         ))
                       }
                     )
                   })
                 })
    
    observeEvent(values[["core_variables_tab1_main"]],
                 {
                   output$core_variables_compare_tab_ui <- renderUI({
                     lapply(
                       values[["core_variables_tab1_main"]],
                       FUN = function(x) {
                         choices <-
                           metadata(rvc)[[values[["analysis_tab1_main"]]]][[x]]
                         return(selectInput(
                           inputId = sprintf("%s_tab1_compare", x),
                           label = x,
                           choices = choices,
                           selected = choices[1]
                         ))
                       }
                     )
                   })
                 })
    
    
    ## Inputs additional filter variables
    output$addvariables_ui_tab1_main <- renderUI({
      if (length(filters_tab1_main[["add_filters"]]) > 0) {
        lapply(
          filters_tab1_main[["add_filters"]],
          FUN = function(x) {
            type <-
              rvc@metadata[[values[["analysis_tab1_main"]]]]$column_types[x]
            dat <- isolate(filters_tab1_main[["filters"]])
            if (type %in% c("numeric", "integer")) {
              return(fluidRow(column(
                6,
                numericInput(
                  inputId = sprintf("%s_tab1_main_min", x),
                  value = if (x %in% names(dat))
                    dat[[x]][["min"]]
                  else
                    0 ,
                  label = sprintf("%s_min", x),
                  min = -Inf,
                  max = Inf
                )
              ),
              column(
                6,
                numericInput(
                  inputId = sprintf("%s_tab1_main_max", x),
                  value = if (x %in% names(dat))
                    dat[[x]][["max"]]
                  else
                    1,
                  label = sprintf("%s_max", x),
                  min = -Inf,
                  max = Inf
                )
              )))
            } else if (type %in% c("character", "logical", "factor", "Rle")) {
              if("rvatViewerClassSQL" %in% class(rvc)) {
                choices <-
                  rvc@rvatResults %>% dplyr::tbl(values[["analysis_tab1_main"]]) %>% dplyr::pull(x) %>% unique()
              } else {
                choices <- unique(rvc@rvatResults[[values[["analysis_tab1_main"]]]][[x]])
              }
              
              return(
                selectInput(
                  inputId = sprintf("%s_tab1_main", x),
                  label = x,
                  choices = choices,
                  selected = choices[1]
                )
              )
            }
          }
        )
      }
    })
    
    output$addvariables_ui_tab1_compare <- renderUI({
      if (length(filters_tab1_compare[["add_filters"]]) > 0) {
        lapply(
          filters_tab1_compare[["add_filters"]],
          FUN = function(x) {
            type <-
              rvc@metadata[[values[["analysis_tab1_main"]]]]$column_types[x]
            dat <- isolate(filters_tab1_compare[["filters"]])
            if (type %in% c("numeric", "integer")) {
              return(fluidRow(column(
                6,
                numericInput(
                  inputId = sprintf("%s_tab1_compare_min", x),
                  value = if (x %in% names(dat))
                    dat[[x]][["min"]]
                  else
                    0 ,
                  label = sprintf("%s_min", x),
                  min = -Inf,
                  max = Inf
                )
              ),
              column(
                6,
                numericInput(
                  inputId = sprintf("%s_tab1_compare_max", x),
                  value = if (x %in% names(dat))
                    dat[[x]][["max"]]
                  else
                    1,
                  label = sprintf("%s_max", x),
                  min = -Inf,
                  max = Inf
                )
              )))
            } else if (type %in% c("character", "logical", "factor", "Rle")) {
              if("rvatViewerClassSQL" %in% class(rvc)) {
                choices <-
                  rvc@rvatResults %>% dplyr::tbl(values[["analysis_tab1_main"]]) %>% dplyr::pull(x) %>% unique()
              } else {
                choices <- unique(rvc@rvatResults[[values[["analysis_tab1_main"]]]][[x]])
              }
              
              return(
                selectInput(
                  inputId = sprintf("%s_tab1_compare", x),
                  label = x,
                  choices = choices,
                  selected = choices[1]
                )
              )
            }
          }
        )
      }
    })
    
    ## Reset additional variables button
    output$addvariables_ui_tab1_main_reset <- renderUI({
      if (length(filters_tab1_main[["add_filters"]]) > 0) {
        actionButton("reset_addvariables_tab1_main",
                     "Reset additional variables")
      }
    })
    
    ## Reset additional variables button - compare tab
    output$addvariables_ui_tab1_compare_reset <- renderUI({
      if (length(filters_tab1_compare[["add_filters"]]) > 0) {
        actionButton("reset_addvariables_tab1_compare",
                     "Reset additional variables")
      }
    })
    
    ## Set hover fields
    output$tab1_set_hover_fields_ui <- renderUI({
      if(input$tab1_edit_hover_variables) {
        checkboxGroupInput("tab1_set_hover_fields",
                           label = "Fields",
                           choices = setdiff(metadata(rvc)[[values[["analysis_tab1_main"]]]][["columns"]], c("unit", "alias")),
                           selected = values[["tab1_hover_fields"]],
                           inline = TRUE)
        
      }
    })
    
    ## Select columns topTable
    output$tab1_edit_variables_toptable_ui <- renderUI({
      if(input$tab1_edit_variables_toptable) {
        checkboxGroupInput("tab1_set_variables_toptable",
                           label = "Fields",
                           choices = metadata(rvc)[[values[["analysis_tab1_main"]]]][["columns"]],
                           selected = c("unit", "P", "effect", "effectSE", "effectCIlower", "effectCIupper"),
                           inline = TRUE)
        
      }
    })
    
    ## Output -----------------------------------------------------------------
    
    # rvb panel ---------------------------------------------------------------
    
    ## Manhattan plot
    
    output$manhattan <- plotly::renderPlotly({
      shiny::validate(need(nrow(data_tab1_main()) > 0, "No data available for this combination of parameters"))
      shiny::validate(need(all(c("CHROM", "POS") %in% rvc@metadata[[values[["analysis_tab1_main"]]]]$columns) ||
                           (!is.null(CHROM) && !is.null(POS)), 
                           "'CHROM' and 'POS' columns are required for manhattan plot.
If these columns are present, but have other column names, this can be specified
using the `CHROM` and `POS` parameters. 
                           "))
      if(is.null(CHROM)) CHROM <- "CHROM"
      if(is.null(POS)) POS <- "POS"
      dat <- data_tab1_main() 
      dat$effect <- round(dat$effect, 2)
      if(stringr::str_detect(dat[[CHROM]][1], "chr")) {
        dat[["CHR"]] <- as.numeric(stringr::str_replace(dat[[CHROM]], "chr", ""))
      } else {
        dat[["CHR"]] <- as.numeric(dat[[CHROM]])
      }
      dat[["BP"]] <- as.numeric(dat[[POS]])
      dat <- dat[,unique(c("unit", "CHR", "BP", input[["tab1_display_variable"]], "P", values[["tab1_hover_fields"]])),drop=FALSE]
      dat <- dat[complete.cases(dat[,c("CHR", "BP")]),]
      col <- paste0(dat[[values[["tab1_hover_fields"]][1]]], "\n")
      
      if(length(values[["tab1_hover_fields"]]) > 1) {
        for(i in 2:(length(values[["tab1_hover_fields"]]))){
          col <- paste0(col, values[["tab1_hover_fields"]][i], ": ", dat[[values[["tab1_hover_fields"]][i]]], "\n")
        } 
      }
      dat[[values[["tab1_hover_fields"]][[1]]]] <- col
      dat <- manhattanr(dat, 
                        gene = input[["tab1_display_variable"]], 
                        annotation1 =  values[["tab1_hover_fields"]][[1]])
      if(length(values[["custom_threshold"]]) > 0 && input[["use_custom_threshold"]]) {
        threshold <- values[["custom_threshold"]]
      } else {
        threshold <-  0.05/nrow(dat$data)
      }
      dat_sig = dat$data[dat$data[["P"]] < threshold,,drop=FALSE]
      
      if(nrow(dat_sig) > 0) {
        plot <- manhattanly(dat,
                           suggestiveline = FALSE, 
                           genomewideline = -log10(threshold),
                           title = "") %>%
          plotly::add_annotations(x = dat_sig$pos,
                                  y = dat_sig$logp,
                                  text = dat_sig[[input[["tab1_display_variable"]]]],
                                  xref = "x",
                                  yref = "y",
                                  showarrow = TRUE,
                                  arrowhead=2,
                                  ax = 10,
                                  ay = -10) 
      } else {
        plot <- manhattanly(dat,
                       suggestiveline = FALSE, 
                       genomewideline = -log10(threshold),
                       title = "") 
      }
      
      plotly::toWebGL(plot)
    })
    
    ## qqplot
    output$qqplot <- renderPlot({
      dat <- data_tab1_main()
      if(length(values[["custom_threshold"]]) > 0 && input[["use_custom_threshold"]]) {
        threshold <- values[["custom_threshold"]]
      } else {
        threshold <-  0.05/nrow(dat)
      }
      plot <- .qqplot(P = dat[,c("P"),drop=FALSE], threshold = threshold)
      plot
    })
    
    # 
    output$tab1_analysis_info <- renderText({
      dat <- data_tab1_main()
      sprintf(
        "<b>N</b> = %s<br><b>N cases</b> = %s<br><b>N controls</b> = %s<br><b>N units</b> = %s<br><b>&lambda;</b> = %s",
        max(dat$caseN + dat$ctrlN),
        max(dat$caseN),
        max(dat$ctrlN),
        dplyr::n_distinct(dat$unit),
        signif(median(qchisq(
          1 - dat$P, 1
        ), na.rm = TRUE) / qchisq(0.5, 1), 3)
      )
    })
    
    output$topTable <- DT::renderDataTable(
      {
        data_tab1_main()[,values[["tab1_variables_toptable"]],drop=FALSE] %>% dplyr::arrange(P) %>% head(input[["tab1_toptable_nrows"]])
      })
    
    output$compare_tab1 <- plotly::renderPlotly({
      
      dat1 <- data_tab1_main()
      dat2 <- data_tab1_compare()
      shiny::validate(need(nrow(dat1) > 0 && nrow(dat2) > 0, 
                           "No data available for this combination of parameters"))
      
      ## Check if rows in both data.frames are uniquely defined by the `unit` column
      check1 <- dplyr::n_distinct(dat1$unit) == nrow(dat1)
      check2 <- dplyr::n_distinct(dat2$unit) == nrow(dat2)
      shiny::validate(need(check1 && check2, 
                           "The join variable should uniquely define rows"))
      
      shiny_compare_plot(
        dat1 = dat1,
        dat2 = dat2,
        compare_variable = input[["tab1_compare_variable"]],
        transformation = input[["tab1_compare_transformation"]],
        custom_threshold = values[["custom_threshold"]],
        use_custom_threshold = input[["use_custom_threshold"]],
        hover_fields = values[["tab1_hover_fields"]],
        show_labels = input[["compare_show_labels"]],
        display_variable = input[["tab1_display_variable"]]
      )
    })
  }
  shinyApp(ui, server)
}
