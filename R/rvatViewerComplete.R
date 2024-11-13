# ==============================================================================
#  results viewer complete
# ==============================================================================

#' rvatViewerComplete
#'
#' Interactively explore variants or domains, the results of variant/burden testing, and the results of gene set analyses
#'
#' @param object Can be 1) an object of class \code{\link{rvbResult}}, \code{\link{singlevarResult}}, \code{\link{gsaResult}},
#' or \code{\link{gdb}}
#' 2) a filepath to an SQL-database, 
#' 3) a `data.frame`, or 
#' 4) a list of objects of class \code{\link{rvbResult}}, \code{\link{singlevarResult}}, or `data.frame`
#' @param gdb Can be 1) an object of class \code{\link{gdb}} or the name of the table in an SQL-database with the gdb data
#' @param varSet varSet as a varSetFile filepath, varSetList, data.frame/list with varSets, or varSet
#' @param varInfo variant info in a data.frame, list of data.frames, or filepath. If gdb != NULL, annotation tables will be loaded in as additional varInfo tables. 
#' The variant info tables should have a `VAR_id` column. If unit IDs are ensembl IDs, it is recommended to include a 'geneID' column.
#' @param mergeVarInfo whether the varInfo tables from `gdb` and `varInfo` should be merged by their overlapping columns, default = FALSE
#' @param geneSetList geneSetList object or filepath
#' @param tracks a (list of) file.name or data.frame in BED format (with a 'CHROM', 'start', 'end' column) with the addition of a 'track' and 'unit' column that is used to plot tracks in the unit viewer. 
#' @param CHROM if the name of the column that specifies the chromosome is not named `CHROM`, 
#' the column name can be specified here.
#' @param POS if the name of the column that specifies the chromosomal/variant position is not named `POS`, 
#' the column name can be specified here.
#' @param UNIT if the name of the column in variant info that specifies the unit is not named `unit`,
#' the column name can be specified here. 
#' @param ID if the name of the column in cohort tables in the gdb that specifies the sample IDs is not `IID`, 
#' the column name can be specified here. Default is `IID`
#' @import shiny
#' 
#' @export

rvatViewerComplete <- function(object, gdb = NULL, varSet = NULL, varInfo = NULL, mergeVarInfo = FALSE, geneSetList = NULL, 
                               tracks = NULL, CHROM = NULL, POS = NULL, UNIT = NULL, ID = 'IID') {
  # Data preparation -----------------------------------------------------------
  rvc <- prepareDataShiny(object, gdb = gdb, varSet = varSet, varInfo = varInfo, mergeVarInfo = mergeVarInfo, 
                          geneSetList = geneSetList, tracks = tracks) 
  
  convert_vec_logical <- c(0, 1, 0, 1)
  names(convert_vec_logical) <- c("FALSE", "TRUE", "0", "1")
  unit_core_filter_variables <-   
    c("start",
      "end",
      "track",
      "unit",
      "y0", "y1", "colours", "unitShort")
  
  assoc_core_filter_variables <-
    c("analysis",
      "cohort",
      "varSetName",
      "name",
      "pheno",
      "covar",
      "geneticModel",
      "MAFweight",
      "test")
  
  geneset_core_filter_variables <-
    c("analysis",
      "method",
      "test",
      "covar")
  
  analyses_types = c()
  for (analysis in metadata(rvc)$analyses) {
    analyses_types <- c(analyses_types, checkClassrvatResult(rvc@rvatResults[[analysis]]))
  }
    
  default_assoc <- metadata(rvc)$analyses[analyses_types == "rvbResult" | analyses_types == "singlevarResult"][1]
  if (length(analyses_types[analyses_types == "gsaResult"]) == 0) {
    default_geneset <- default_assoc
  } else {
    default_geneset <- metadata(rvc)$analyses[analyses_types == "gsaResult"][1]
  }
  
  ## Check which core variables have more than 1 unique value
  ## Only when length(unique(var)) > 1, a filter option is included 
  core_filter_variables_default_assoc <-
    names(metadata(rvc)[[default_assoc]][assoc_core_filter_variables][lengths(metadata(rvc)[[default_assoc]][assoc_core_filter_variables]) > 1])
  if(length(core_filter_variables_default_assoc) == 0) {
    core_filter_variables_default_assoc = "test"
  }
  
  core_filter_variables_default_geneset <-
    names(metadata(rvc)[[default_geneset]][geneset_core_filter_variables][lengths(metadata(rvc)[[default_geneset]][geneset_core_filter_variables]) > 1])
  if(length(core_filter_variables_default_geneset) == 0) {
    core_filter_variables_default_geneset = "test"
  }
  
  add_filter_variables_unit <- c()
  option <- colnames(tracks)[!(colnames(tracks) %in% unit_core_filter_variables)]
  for (opt in option) {
    if (length(unique(tracks[[opt]])) >1 & typeof(tracks[[opt]]) == "character") {
      add_filter_variables_unit <- c(add_filter_variables_unit, opt)
    }
  }
  
  # UI -------------------------------------------------------------------------
  ui <- fluidPage(
    navbarPage(
      HTML(
        paste0("RVAT viewer", '<sup>', 'v', packageVersion("rvat"))
      ),
      id = "main_menu",
      
      # Data overview ----------------------------------------------------------
      tabPanel(
        "Data overview",
        value = "data_overview",
        h4("Summary of the loaded data"),
        fluidRow(
          width = 12,
          column(
            width = 6,
            selectInput("rvat_analysis_input", label = "Analysis:", 
                        choices = rvc@metadata[["analyses"]], 
                        selected = rvc@metadata[["analyses"]][1]),
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                             tags$div("Info is loading...",id="loadmessage")),
            htmlOutput("data_summary_rvat")
          ),
          column(
            width = 6,
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                             tags$div("Info is loading...",id="loadmessage")),
            htmlOutput("data_summary_geneSetList")
          )),
        fluidRow(
          width = 12,
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Info is loading...",id="loadmessage")),
          htmlOutput("data_summary_gdb"),
          br()
        ),
        fluidRow(
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Info is loading...",id="loadmessage")),
          htmlOutput("data_summary_varSet"),
          br()
        ),
        fluidRow(
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Info is loading...",id="loadmessage")),
          htmlOutput("data_summary_varInfo")
        ),
        fluidRow(
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Info is loading...",id="loadmessage")),
          htmlOutput("data_summary_tracks")
        )
      ),
      
      # Unit viewer ------------------------------------------------------------
      tabPanel(
        "Unit viewer",
        value = "unit_viewer",
        
        sidebarPanel(
          width = 4,
          
          # Info selected data
          h4("Selected data:"),
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Info is loading...",id="loadmessage")),
          htmlOutput("unit_analysis_info"),
          hr(),
          h4("Options:"),
          checkboxInput("unit_use_custom_threshold", label = "Use custom significance threshold?", value = FALSE),
          conditionalPanel(condition = "input.unit_use_custom_threshold == 1",
                           fluidRow(
                             width = 12,
                             column(8, numericInput("unit_custom_threshold", label = NULL, value = 5e-8, min = 0, max = 1)),
                             column(4, actionButton("unit_set_custom_threshold", label = "Set")))),
          ## Add variables, all columns in tracks that aren't in the main columns and have a number of choices > 1 and < nrow(tracks)
          uiOutput("unit_additional_variables_main_ui")
        ),
        
        mainPanel(
          width = 6,
          tabsetPanel(
            tabPanel(
              "Mutation plot",
              htmlOutput("selected_unit_mut"),
              textInput("unit_selected_unit_mut", label = "Select unit to show", placeholder = "UNIT"),
              radioButtons("unit_selected_pvalues", label = "Select P-values to show", 
                           choices = c("Single-variant P-values", "Single-variant logP-values", 
                                       "Track P-values", "Track logP-values"), selected = character(0)),
              selectInput("unit_lollipop_title", label = "Select plot title:", choices = rvc@metadata$varInfo),
              h4("Mutation plot"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Plot is loading...",id="loadmessage")),
              plotly::plotlyOutput("unit_lollipop"),
              h4("Track table:"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::dataTableOutput("unit_track_table"),
              br(),
              h4("Variant table:"),
              checkboxInput("unit_edit_variables_variant_table", label = "Select columns to show"),
              uiOutput("unit_edit_variables_variant_table_ui"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::dataTableOutput("unit_variant_table")
            ),
            tabPanel(
              "Patient table",
              htmlOutput("selected_unit_pat"),
              em("To prevent long loading times, do not select units with very large numbers of variants!"),
              textInput("unit_selected_unit_GT", label = "Select unit to show", placeholder = "UNIT1 or UNIT1,UNIT2,UNIT3 etc"),
              uiOutput("unit_selected_cohort"),
              h4("Patient data"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::dataTableOutput("unit_patient_table"),
              br(),
              h4("GT assay"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::dataTableOutput("unit_GT_assay")
            )
          )
        )
      ),
      
      # Association results ----------------------------------------------------
      tabPanel(
        "Association results",
        value = "association_results",
        
        sidebarPanel(
          width = 4,
          
          # Info selected data
          h4("Selected data:"),
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Info is loading...",id="loadmessage")),
          htmlOutput("assoc_analysis_info"),
          hr(),
          h4("Options:"),
          
          #Set a custom threshold
          checkboxInput("assoc_use_custom_threshold", label = "Use custom significance threshold?", value = FALSE),
          conditionalPanel(condition = "input.assoc_use_custom_threshold == 1",
                           fluidRow(
                             width = 12,
                             column(8, numericInput("assoc_custom_threshold", label = NULL, value = 5e-8, min = 0, max = 1)),
                             column(4, actionButton("assoc_set_custom_threshold", label = "Set")))),

          # Edit hover variables 
          checkboxInput("assoc_edit_hover_variables", label = "Change hovertext?"),
          uiOutput("assoc_set_hover_fields_ui"),
          
          # Display variable for significant values
          selectInput("assoc_display_variable", label = "Displayed text for significant units:",
                      choices = names(metadata(rvc)[[default_assoc]]$column_types)[metadata(rvc)[[default_assoc]]$column_types %in% c("character", "factor", "Rle")],
                      selected = "unit"),
          
          hr(),
          h4("Filters:"),
          
          ## Core variables, dependent on selected analysis which inputs are included
          uiOutput("assoc_analysis_ui"),
          
          ## Core variables, dependent on selected analysis which inputs are included
          uiOutput("assoc_core_variables_main_ui"),
          
          # Additional user selected variables
          h5("Add variables:"),
          fluidRow(width = 12,
                   column(
                     8,
                     selectInput(
                       "assoc_addvariables_select_main",
                       label = NULL,
                       choices = c(metadata(rvc)[[default_assoc]][["columns"]][!metadata(rvc)[[default_assoc]][["columns"]] %in% assoc_core_filter_variables]),
                       selected = "P"
                     )
                   ),
                   column(
                     4, actionButton("assoc_addvariables_add_main", label = "Add")
                   )),
          uiOutput("assoc_addvariables_main_ui"),
          
          # Reset additional variables
          uiOutput("assoc_addvariables_main_ui_reset")
        ),
        
        mainPanel(
          width = 6,
          tabsetPanel(
            tabPanel(
              "Manhattan plot",
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Plot is loading...",id="loadmessage")),
              plotly::plotlyOutput("assoc_manhattan"),
              checkboxInput("assoc_edit_variables_toptable_manh", label = "Select columns to show"),
              uiOutput("assoc_edit_variables_toptable_manh_ui"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::dataTableOutput("assoc_toptable_manh")
            ),
            tabPanel(
              "qqplot",
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Plot is loading...",id="loadmessage")),
              plotly::plotlyOutput("assoc_qqplot"),
              checkboxInput("assoc_edit_variables_toptable_qq", label = "Select columns to show"),
              uiOutput("assoc_edit_variables_toptable_qq_ui"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::dataTableOutput("assoc_toptable_qq")
            ),
            tabPanel(
              "Forestplot",
              textInput("assoc_selected_units_forest", label = "Select units to plot", placeholder = "UNIT1 or UNIT1,UNIT2,UNIT3 etc"),
              plotOutput("assoc_forest_plot"),
              br(),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::dataTableOutput("assoc_forest_table"),
              h4("Columns with only 1 unique value:"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Info is loading...",id="loadmessage")),
              htmlOutput("assoc_forest_unique_columns")
            ),
            tabPanel(
              "Compare",
              fluidRow(width= 12,
                       column(4,
                              h4("Options:"),
                              selectInput(
                                "assoc_compare_variable",
                                label = "Variable to compare:",
                                selected = "P",
                                choices = names(metadata(rvc)[[default_assoc]]$column_types)[metadata(rvc)[[default_assoc]]$column_types %in% c("numeric", "integer")]
                              ),
                              selectInput(
                                "assoc_compare_transformation",
                                label = "Transformation",
                                selected = "-log10",
                                choices = c("-log10", "log10", "none")
                              ),
                              checkboxInput("assoc_compare_show_labels", 
                                            label = "Show labels?", value = 1),
                              hr(),
                              h4("Filters:"),
                              uiOutput("assoc_core_variables_compare_ui"),
                              hr(),
                              h4("Add variables:"),
                              fluidRow(width = 12,
                                       column(
                                         8,
                                         selectInput(
                                           "assoc_addvariables_select_compare",
                                           label = NULL,
                                           choices = c(metadata(rvc)[[default_assoc]][["columns"]][!metadata(rvc)[[default_assoc]][["columns"]] %in% assoc_core_filter_variables]),
                                           selected = "P"
                                         )
                                       ),
                                       column(
                                         4, actionButton("assoc_addvariables_add_compare", label = "Add")
                                       )),
                              uiOutput("assoc_addvariables_compare_ui"),
                              # Reset additional variables
                              uiOutput("assoc_addvariables_compare_ui_reset")
                       ),
                       column(8,
                              br(),
                              actionButton("assoc_compare_button", label = "Compare!"),
                              br(),
                              br(),
                              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                               tags$div("Plot is loading...",id="loadmessage")),
                              plotly::plotlyOutput("assoc_compare")
                       )
              )
            ),
            tabPanel(
              "tableViewer",
              checkboxInput("assoc_edit_variables_tableviewer", label = "Select columns to show"),
              uiOutput("assoc_edit_variables_tableviewer_ui"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::DTOutput("assoc_tableviewer")
            )
          )
        )
      ),
      
      # Gene set analysis ------------------------------------------------------
      tabPanel(
        "Gene set analysis",
        value = "gene_set_analysis",
        
        sidebarPanel(
          width = 4,
          
          # Info selected data
          h4("Selected data:"),
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div("Info is loading...",id="loadmessage")),
          htmlOutput("geneset_analysis_info"),
          hr(),
          h4("Options:"),
          
          #Set a custom threshold
          checkboxInput("geneset_use_custom_threshold", label = "Use custom significance threshold?", value = FALSE),
          conditionalPanel(condition = "input.geneset_use_custom_threshold == 1",
                           fluidRow(
                             width = 12,
                             column(8, numericInput("geneset_custom_threshold", label = NULL, value = 5e-8, min = 0, max = 1)),
                             column(4, actionButton("geneset_set_custom_threshold", label = "Set")))),
          
          # Edit hover variables 
          checkboxInput("geneset_edit_hover_variables", label = "Change hovertext?"),
          uiOutput("geneset_set_hover_fields_ui"),
          
          # Display variable for significant values
          selectInput("geneset_display_variable", label = "Displayed text for significant sets:",
                      choices = names(metadata(rvc)[[default_geneset]]$column_types)[metadata(rvc)[[default_geneset]]$column_types %in% c("character", "factor", "Rle")],
                      selected = "geneSetName"),
          
          hr(),
          h4("Filters:"),
          
          ## Core variables, dependent on selected analysis which inputs are included
          uiOutput("geneset_analysis_ui"),
          
          ## Core variables, dependent on selected analysis which inputs are included
          uiOutput("geneset_core_variables_main_ui"),
          
          # Additional user selected variables
          h5("Add variables:"),
          fluidRow(width = 12,
                   column(
                     8,
                     selectInput(
                       "geneset_addvariables_select_main",
                       label = NULL,
                       choices = c(metadata(rvc)[[default_geneset]][["columns"]][!metadata(rvc)[[default_geneset]][["columns"]] %in% geneset_core_filter_variables]),
                       selected = "P"
                     )
                   ),
                   column(
                     4, actionButton("geneset_addvariables_add_main", label = "Add")
                   )),
          uiOutput("geneset_addvariables_main_ui"),
          
          # Reset additional variables
          uiOutput("geneset_addvariables_main_ui_reset")
          
        ),
        
        mainPanel(
          width = 6,
          tabsetPanel(
            id = "geneset_tabs",
            tabPanel(
              "Manhattan plot",
              value = "geneset_manh",
              htmlOutput("selected_geneset"),
              textAreaInput("geneset_selected_units_manhattan", label = "Select gene set to highlight", placeholder = "Gene set"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Plot is loading...",id="loadmessage")),
              plotly::plotlyOutput("geneset_manhattan"),
              h4("Geneset info table:"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              tableOutput("geneset_info_table_manh"),
              h3("Geneset genes association test table:"),
              checkboxInput("geneset_edit_variables_genes_manh", label = "Select columns to show"),
              uiOutput("geneset_edit_variables_genes_manh_ui"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::dataTableOutput("geneset_genes_table_manh")
            ),
            tabPanel(
              "qqplot",
              value = "geneset_qq",
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Plot is loading...",id="loadmessage")),
              plotly::plotlyOutput("geneset_qqplot"),
              checkboxInput("geneset_edit_variables_toptable_qq", label = "Select columns to show"),
              uiOutput("geneset_edit_variables_toptable_qq_ui"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::dataTableOutput("geneset_toptable_qq")
            ),
            tabPanel(
              "Density plot",
              value = "geneset_density",
              textAreaInput("geneset_selected_units_density", label = "Select gene set to plot", placeholder = "Gene set"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Info is loading...",id="loadmessage")),
              plotly::plotlyOutput("geneset_density_plot"),
              br(),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::dataTableOutput("geneset_density_table"),
              h4("Columns with only 1 unique value:"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Info is loading...",id="loadmessage")),
              htmlOutput("geneset_density_unique_columns")
            ),
            tabPanel(
              "Forest plot",
              value = "geneset_forest",
              textAreaInput("geneset_selected_units_forest", label = "Select gene sets to plot", placeholder = "GS1 or GS1,GS2,GS3 etc"),
              plotOutput("geneset_forest_plot"),
              br(),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::dataTableOutput("geneset_forest_table"),
              h4("Columns with only 1 unique value:"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Info is loading...",id="loadmessage")),
              htmlOutput("geneset_forest_unique_columns")
            ),
            tabPanel(
              "tableViewer",
              value = "geneset_tableviewer",
              checkboxInput("geneset_show_index_tableviewer", label = "Show rownames?"),
              checkboxInput("geneset_edit_variables_tableviewer", label = "Select columns to show"),
              uiOutput("geneset_edit_variables_tableviewer_ui"),
              conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                               tags$div("Table is loading...",id="loadmessage")),
              DT::DTOutput("geneset_tableviewer")
            )
          )
        )
      )
    )
  ) 
  
  server <- function(input, output, session) {
    # General ------------------------------------------------------------------
    ## Selected data
    values <- reactiveValues(
      gdb_present = rvc@metadata$gdb,
      selected = c(),
      selected_cohort  = c(),
      selected_unit = c(),
      selected_geneset = c(),
      fields = c(),
      lookup = c(),
      unit_lookup = c(),
      alias_lookup = c(),
      unit_filter_variables_main = add_filter_variables_unit,
      unit_custom_threshold = c(),
      unit_variables_table = c(),
      assoc_analysis_main = default_assoc,
      assoc_core_variables_main = core_filter_variables_default_assoc,
      assoc_hover_fields = c("effect"),
      assoc_variables_table = c("unit", "P", "effect", "effectSE", "effectCIlower", "effectCIupper"),
      assoc_custom_threshold = c(),
      geneset_analysis_main = default_geneset,
      geneset_core_variables_main = core_filter_variables_default_geneset,
      geneset_hover_fields = c("effect"),
      geneset_variables_table = c("geneSetName","P", "effect", "effectSE", "effectCIlower", "effectCIupper"),
      geneset_custom_threshold = c()
    )
    
    # Data overview ------------------------------------------------------------
    output$data_summary_rvat <- renderText({
      summary_text <- ""
      dat_class = checkClassrvatResult(rvc@rvatResults[[input$rvat_analysis_input]])
      dat <- summary(rvc@rvatResults[[input$rvat_analysis_input]], asList = TRUE)
      if (identical(dat_class, "gsaResult")) {
        summary_text <- paste(summary_text, 
                              "<b>",input$rvat_analysis_input, ":</b><br>",
                              sprintf(
                                "<b>N gene sets</b> = %s<br><b>Tests</b> = %s<br><b>Methods</b> = %s<br><b>Covar</b> = %s<br><b>N significant positive effect (P<0.05; bonf. corrected)</b> = %s<br><b>N significant negative effect (P<0.05; bonf. corrected)</b> = %s",
                                dat[["Ngenesets"]],
                                paste(dat[["tests"]], collapse=","),
                                paste(dat[["methods"]],collapse=","),
                                paste(dat[["covar"]],collapse=","),
                                dat[["NsignPosEffect"]],
                                dat[["NsignNegEffect"]]
                              ),
                              "<br><br>")
      } else if (identical(dat_class, "rvbResult")) {
        summary_text <- paste(summary_text,
                              "<b>",input$rvat_analysis_input,":</b><br>",
                              sprintf(
                                "<b>N units</b> = %s<br><b>N cases</b> = %s<br><b>N controls</b> = %s<br><b>Names</b> = %s<br><b>varSets</b> = %s<br><b>Tests</b> = %s<br><b>geneticModels</b> = %s<br><b>MAFweights</b> = %s<br><b>Pheno</b> = %s<br><b>Cohorts</b> = %s",
                                dat[["Nunits"]],
                                dat[["caseNmax"]],
                                dat[["ctrlNmax"]],
                                paste(dat[["name"]], collapse = ","),
                                paste(dat[["varSetName"]], collapse = ","),
                                paste(dat[["test"]], collapse = ","),
                                paste(dat[["geneticModel"]], collapse = ","),
                                paste(dat[["MAFweight"]], collapse = ","),
                                paste(dat[["pheno"]], collapse = ","),
                                paste(dat[["cohort"]], collapse = ",")
                              ),
                              "<br><br>")
      } else {
        summary_text <- paste(summary_text,
                              "<b>",input$rvat_analysis_input,":</b><br>",
                              sprintf(
                                "<b>N vars</b> = %s<br><b>N cases</b> = %s<br><b>N controls</b> = %s<br><b>Names</b> = %s<br><b>varSets</b> = %s<br><b>Tests</b> = %s<br><b>geneticModels</b> = %s<br><b>Pheno</b> = %s<br><b>Cohorts</b> = %s",
                                dat[["Nvars"]],
                                dat[["caseNmax"]],
                                dat[["ctrlNmax"]],
                                paste(dat[["names"]], collapse = ","),
                                paste(dat[["varsets"]], collapse = ","),
                                paste(dat[["tests"]], collapse = ","),
                                paste(dat[["geneticModels"]], collapse = ","),
                                paste(dat[["pheno"]], collapse = ","),
                                paste(dat[["cohorts"]], collapse = ",")
                              ),
                              "<br><br>")
      }
    })
    
    output$data_summary_geneSetList <- renderText({
      summary_text = ""
      if (identical(rvc@geneSetList, geneSetList())) {
        summary_text <- paste(summary_text,
                              "<b>geneSetList</b><br>The loaded data does not contain a geneSetList<br><br>")
      } else {
        dat <- rvc@geneSetList
        summary_text <- paste(summary_text, "<b>geneSetList</b><br>",
                              sprintf("<b>N gene sets</b> = %s<br><b>N unique units</b> = %s<br><b>Min N units</b> = %s<br><b>Max N units</b> = %s<br><b>Unique weights</b> = %s",
                                      length(dat),
                                      length(unique(unlist(lapply(dat@geneSets, function(x) {strsplit(x@units, ",")})))),
                                      min(lengths(dat)),
                                      max(lengths(dat)),
                                      paste(unique(unlist(lapply(dat@geneSets, function(x){strsplit(x@w, ",")}))),collapse=",")))
      }
    })
    
    output$data_summary_gdb <- renderText({
      summary_text = ""
      if (substring(rvc@gdb@dbname, nchar(rvc@gdb@dbname)-3) == "NULL") {
        summary_text <- paste(summary_text,
                              "<b>GDB</b><br>The loaded data does not contain a gdb<br>")
      } else {
        dat <- rvc@gdb
        summary_text <- paste(summary_text,"<b>gdb</b><br>",
                              sprintf("<b>Filepath:</b> %s<br><b>Annotation tables</b> = %s<br><b>Cohorts</b> = %s",
                                      rvc@metadata$gdb,
                                      paste0(listAnno(dat)[,1], collapse = ","),
                                      paste0(listCohort(dat)[,1], collapse = ",")))
      }
      
      summary_text
    })
    
    output$data_summary_varSet <- renderText({
      summary_text = ""
      #if (is.null(rvc@varSet)) {
      if (identical(rvc@varSet, varSetList())) {
        summary_text <- paste(summary_text,
                              "<b>varSet</b><br>The loaded data does not contain varSets<br><br>")
      } else {
        if (length(rvc@metadata$varSet) == 2) {
          summary_text <- paste(summary_text, "<b>varSet</b><br>",
                                sprintf("<b>Filepath</b> = %s<br><b>Units</b> = %s",
                                        rvc@metadata$varSet[["filepath"]],
                                        as.character(length(rvc@metadata$varSet[["units"]]))))
        } else {
          summary_text <- paste(summary_text, "<b>varSet</b><br>",
                                sprintf("<b>Units: </b><br>%s",
                                        as.character(length(rvc@metadata$varSet[["units"]]))))
        }
      }
    })
    
    output$data_summary_varInfo <- renderText({
      summary_text = "<b>varInfo</b>"
      dat <- rvc@metadata$varInfo
      if (length(rvc@varInfo) == 0) {
        summary_text <- paste(summary_text,
                              "<br>The loaded data does not contain varInfo<br><br>")
      } else {
        columns <- rep("", ceiling(length(dat)/10))
        for (i in 1:ceiling(length(dat)/10)) {
          if (i < ceiling(length(dat)/10)) {
            columns[i] <- paste(dat[(1+(i-1)*10):(i*10)], collapse = ",")
          } else {
            columns[i] <- paste(dat[(1+(i-1)*10):length(dat)], collapse = ",")
          }
        }

        summary_text <- paste(summary_text,
                              sprintf("<b>Columns:</b><br>%s<br><br>",
                                      paste0(columns, collapse = "<br>")))
      }
      
      summary_text
    })
    
    output$data_summary_tracks <- renderText({
      summary_text = "<b>Tracks</b><br>"
      dat <- rvc@metadata$tracks
      if (rvc@metadata$tracks == 0) {
        summary_text <- paste(summary_text,
                              "The loaded data does not contain tracks<br><br>")
      } else {
        summary_text <- paste(summary_text,
                              sprintf("<b>Track names:</b><br> %s<br><br>",
                                      paste0(dat, collapse = ",")))
      }
      
      summary_text
    })
    
    # Unit viewer --------------------------------------------------------------
    ## Side panel --------------------------------------------------------------
    ### Reactives --------------------------------------------------------------
    data_unit_varInfo <- eventReactive(values$selected_unit, {
      unit <- unlist(strsplit(values$selected_unit, ","))
      unitShort <- unlist(strsplit(unit, "_", fixed = TRUE))[1]
      unitShort <- unlist(strsplit(unitShort, ".", fixed = TRUE))[1]
      
      varPos <- subVarPos(rvc = rvc, unit = unit, unitShort = unitShort, UNIT = UNIT, 
                          analysis = values[["assoc_analysis_main"]], filters = assoc_filters_main[["filters"]])
      varPos
    })
    
    data_unit_lollipop <- reactive({
      dat <- as.data.frame(rvc@tracks)
    
      dat
    })
    
    unit_filters_main <- reactiveValues(
      column = list(add_filter_variables_unit),
      value = list(rep("All", length(add_filter_variables_unit)))
    )
    
    data_unit_rvb <- reactive({
      if (length(assoc_filters_main[["filters"]]) > 0) {
        dat <- subSet(rvc, analysis = values[["assoc_analysis_main"]], filters = assoc_filters_main[["filters"]])
      } else {
        dat <- data.frame(rvc@rvatResults[[values[["assoc_analysis_main"]]]])
      }
      
      dat
    })
    
    ### Selected data info -----------------------------------------------------
    output$unit_analysis_info <- renderText({
      dat <- data_unit_rvb()
      dat <- dat[dat$unit == values[["selected_unit"]],]
      sprintf(
        "<b>Unit</b> = %s<br><b>Nvar</b> = %s<br><b>ctrlCarriers</b> = %s<br><b>caseCarriers</b> = %s<br><b>Tests</b> = %s",
        values[["selected_unit"]],
        max(dat$nvar),
        max(dat$ctrlCarriers),
        max(dat$caseCarriers),
        paste(unique(dat$test), collapse = ",")
      )
    })
    
    ### Set custom threshold ---------------------------------------------------
    observeEvent(input[["unit_set_custom_threshold"]],
                 {
                   values[["unit_custom_threshold"]] <- input[["unit_custom_threshold"]]
                 })
    
    ### Filter variables -------------------------------------------------------
    #Makes the buttons for each of the values$unit_filter_variables_main
    output$unit_additional_variables_main_ui <- renderUI({
        lapply(values[["unit_filter_variables_main"]],
               function(x){
                 return(
                    selectInput(inputId = sprintf("%s_unit_filter", x),
                                label = x,
                                choices = c("All", unique(tracks[[x]])),
                                selected = "All")
                 )
          
               })
    })
    
    #Have to add filters to unit_filters_main for each of the buttons
    observe({
      req(values[["unit_filter_variables_main"]])
      
      #Extract input for all additional columns
      column <- list(values[["unit_filter_variables_main"]])
        
      value <- lapply(
        values[["unit_filter_variables_main"]],
        function(x) {
          input[[sprintf("%s_unit_filter", x)]]
        }
      )
      
      #Update the unit_filter_main object
      unit_filters_main$column <- column
      unit_filters_main$value <- value
    })
    
    ## Tabs --------------------------------------------------------------------
    ### Mutation Plot ----------------------------------------------------------
    #### Reactives -------------------------------------------------------------
    observeEvent(input[["unit_selected_pvalues"]], {
      values[["selected_pvalues"]] <- input[["unit_selected_pvalues"]]
    })
    
    observeEvent(input[["unit_lollipop_title"]], {
      values[["unit_lollipop_title_column"]] <- input[["unit_lollipop_title"]]
    })
    
    observeEvent(input[["unit_selected_unit_mut"]], {
      values[["selected_unit"]] <- input[["unit_selected_unit_mut"]]
    })
    
    output$selected_unit_mut <- renderText({
      paste("<b>",values$selected_unit,"</b><br>")
    })
    
    #### Lollipop plot ----------------------------------------------------------
    output$unit_lollipop <- plotly::renderPlotly({
      unit <- unlist(strsplit(values$selected_unit, ","))
      shiny::validate(need(length(unit) == 1, "You need to select 1 unit!"))
      
      unitShort <- unlist(strsplit(unit, "_", fixed = TRUE))[1]
      unitShort <- unlist(strsplit(unitShort, ".", fixed = TRUE))[1]
    
      #Make the varPos table necessary for each type of plot
      varPos <- data_unit_varInfo()
      
      if ("singlevarResult" %in% sapply(rvc@rvatResults, checkClassrvatResult)) {
        varPos <- svJoin(rvc,varPos)
      }
      
      #Make the domain information necessary for each type of plot
      tracks <- data_unit_lollipop()
      
      shiny::validate(need(length(unique(tracks$track)) <= 8, "You can plot <= 8 tracks only!"))
      if (nrow(tracks) == 0 & length(values[["selected_pvalues"]]) > 0) {
        if (values[["selected_pvalues"]] == "Track P-values" | values[["selected_pvalues"]] == "Track logP-values") {
          showNotification("The `Tracks` mode has been selected, but no tracks were available for this unit")
        }
      }
      
      #For some reason this has to be here instead of in the reactive value data_unit_lollipop, otherwise the filtering doesn't work
      if (length(add_filter_variables_unit) > 0) {
        for (i in 1:length(unit_filters_main[["column"]])) {
          if (is.null(unit_filters_main$value[[i]])) {
            tracks <- tracks
          } else if (unit_filters_main$value[[i]] == "All") {
            tracks <- tracks
          } else {
            tracks <- tracks[tracks[[unit_filters_main$column[[i]]]] %in% unit_filters_main$value[[i]],]
          }
        }
      }
  
      #Edit tracks table so it has all the information for the different types of plots
      if (unit != unitShort) {
        track <- unique(tracks[tracks$unit == unit, "track"])
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
      
      ##Get significance scores for whole gene and ACATed
      rvc_data <- rvc@rvatResults[[values$assoc_analysis_main]]
      rvc_data <- rvc_data[rvc_data$unitShort == unitShort,]
      if (length(add_filter_variables_unit) > 0) {
        for (i in 1:length(unit_filters_main[["column"]])) {
          if (is.null(unit_filters_main$value[[i]])) {
            rvc_data <- rvc_data
          } else if (unit_filters_main$value[[i]] == "All") {
            rvc_data <- rvc_data
          } else {
            rvc_data <- rvc_data[rvc_data[[unit_filters_main$column[[i]]]] %in% unit_filters_main$value[[i]],]
          }
        }
      }
      
      rvc_dataWhole <- rvc_data[rvc_data$unit == rvc_data$unitShort,]
      if (!is.null(UNIT)) {
        wholeGene <- mean(rvc_dataWhole[rvc_dataWhole[[UNIT]] == unitShort,"P"])
      } else {
        wholeGene <- mean(rvc_dataWhole[rvc_dataWhole$unit == unitShort,"P"])
      }
      
      rvc_dataTrack <- rvc_data[rvc_data$unit != rvc_data$unitShort,]
      group = c("cohort","varSetName", "name", "unit", "pheno", "covar", "geneticModel", "MAFweight", "test")
      columns = unname(unlist(sapply(group, function(x) {ifelse(length(unique(rvc_dataTrack[[x]])) != 1, TRUE, FALSE)})))
      
      if (is.null(values[["selected_pvalues"]])) {
        acatP <- NULL
      } else if (values[["selected_pvalues"]] == "Track P-values" | values[["selected_pvalues"]] == "Track logP-values") {
        rvc_dataACAT <- ACAT(rvc_dataTrack, group = c(group, "unitShort"), aggregate = group[columns], fixpval_method = "Liu")
        acatP <- rvc_dataACAT$P
      } else if (values[["selected_pvalues"]] == "Single-variant P-values" | values[["selected_pvalues"]] == "Single-variant logP-values") {
        acatP <- NULL
      }
      
      plot <- mutationly(varPos = varPos, tracks = tracks, rvc = rvc, plotMode = values[["selected_pvalues"]], 
                         POS = POS, title = values[["unit_lollipop_title_column"]], significance = values[["unit_custom_threshold"]],
                         wholeGeneP = wholeGene, acatP = acatP)
      
      plotly::toWebGL(plot)
    })
    
    #### Table with track info -------------------------------------------------
    output$unit_track_table <- DT::renderDataTable({
      unit <- unlist(strsplit(values$selected_unit, ","))
      shiny::validate(need(length(unit) == 1, "You need to select 1 unit!"))
      
      unitShort <- unlist(strsplit(unit, "_", fixed = TRUE))[1]
      unitShort <- unlist(strsplit(unitShort, ".", fixed = TRUE))[1]
      
      varPos <- data_unit_varInfo()
      
      tracks <- data_unit_lollipop()
      
      if (length(add_filter_variables_unit) > 0) {
        for (i in 1:length(unit_filters_main[["column"]])) {
          if (is.null(unit_filters_main$value[[i]])) {
            tracks <- tracks
          } else if (unit_filters_main$value[[i]] == "All") {
            tracks <- tracks
          } else {
            tracks <- tracks[tracks[[unit_filters_main$column[[i]]]] %in% unit_filters_main$value[[i]],]
          }
        }
      }
      
      if (unit != unitShort) {
        #If the unit is not a whole gene, include all the information in tracks that falls between the start and 
        #end of the selected unit.
        start = tracks[tracks$unit == unit, "start"]
        end = tracks[tracks$unit == unit, "end"]
        tracks <- tracks[tracks$unit == unit | ((tracks$track != track & tracks$unitShort == unitShort) &
                                                  (data.table::inrange(tracks$start, start, end) | 
                                                     data.table::inrange(tracks$end, start, end))),]
      } else {
        tracks <- tracks[tracks$unitShort == unitShort,]
      }
      
      tracks <- tracks[,!names(tracks) %in% c("colours", "y0", "y1")]
      
      shiny::validate(need(nrow(tracks) > 0, "No tracks could be shown. Do your units all have the same type of ID?"))
      
      if (!is.null(values[["selected_pvalues"]])) {
        if (length(unique(tracks$unit)) < nrow(tracks) & (values[["selected_pvalues"]] == "Track P-values" | values[["selected_pvalues"]] == "Track logP-values")) {
          showNotification("Some of the plotted P-values belong to the same unit!")
        }
      }
      
      DT::datatable(tracks, filter = "top")
    })
    
    
    #### Table with variant info -----------------------------------------------
    observeEvent(input$unit_variant_table_selected_cols, {
      values$unit_variables_table <- input$unit_variant_table_selected_cols
    })
    
    output$unit_edit_variables_variant_table_ui <- renderUI({
      if (input$unit_edit_variables_variant_table) {
        unit <- unlist(strsplit(values$selected_unit, ","))
        shiny::validate(need(length(unit) == 1, "You need to select 1 unit!"))
        
        dat <- data_unit_varInfo()
        shiny::validate(need((!is.null(varSet) & unit %in% rvc@metadata$varSet$units) | "unit" %in% colnames(dat) | !is.null(UNIT), 
                             "No links could be made between the selected unit and the variant information. Either load a (more complete) varSet object, variant info with a 'unit' column, or specify UNIT"))
        
        if (!is.null(values[["selected_pvalues"]])) {
          if ((values[["selected_pvalues"]] == "Single-variant P-values" | 
               values[["selected_pvalues"]] == "Single-variant logP-values") & 
              "singlevarResult" %in% sapply(rvc@rvatResults, checkClassrvatResult)) {
            dat$VAR_id <- as.character(dat$VAR_id)
            dat <- svJoin(rvc, dat)
          }
        }
        
        choices <- colnames(dat)
        choices <- choices[!(choices %in% c("P", "logP"))]
        
        checkboxGroupInput("unit_variant_table_selected_cols", 
                           "Columns to show:", 
                           choices = choices, 
                           selected = choices[1:5], 
                           inline = TRUE)
      }
    })
    
    output$unit_variant_table <- DT::renderDataTable({
      unit <- unlist(strsplit(values$selected_unit, ","))
      shiny::validate(need(length(unit) == 1, "You need to select 1 unit!"))
      
      dat <- data_unit_varInfo()
      shiny::validate(need((!is.null(varSet) & unit %in% rvc@metadata$varSet$units) | "unit" %in% colnames(dat) | !is.null(UNIT), 
                           "No links could be made between the selected unit and the variant information. Either load a (more complete) varSet object, variant info with a 'unit' column, or specify UNIT"))
      
      
      if (!is.null(values[["selected_pvalues"]])) {
        if ((values[["selected_pvalues"]] == "Single-variant P-values" | 
             values[["selected_pvalues"]] == "Single-variant logP-values") & 
            "singlevarResult" %in% sapply(rvc@rvatResults, checkClassrvatResult)) {
          dat$VAR_id <- as.character(dat$VAR_id)
          dat <- svJoin(rvc, dat)
          
          shiny::validate(need(nrow(dat) > 0, "There are no matches. Did you use the correct unit name, e.g. unit vs gene_name?"))
          
          if (length(values$unit_variables_table) == 0) {
            DT::datatable(dat[,c(colnames(dat)[1:5],"P", "logP")], filter = "top")
          } else {
            DT::datatable(dat[,c(input$unit_variant_table_selected_cols, "P", "logP")], filter = "top")
          }
        } else {
          if (length(values$unit_variables_table) == 0) {
            DT::datatable(dat[,colnames(dat[1:5])], filter = "top")
          } else {
            DT::datatable(dat[,input$unit_variant_table_selected_cols], filter = "top")
          }
        }
      } else {
        if (length(values$unit_variables_table) == 0) {
          DT::datatable(dat[,colnames(dat[1:5])], filter = "top")
        } else {
          DT::datatable(dat[,input$unit_variant_table_selected_cols], filter = "top")
        }
      }
    })
    
    ### Patient table ----------------------------------------------------------
    #### Reactive values -------------------------------------------------------
    observeEvent(input[["unit_selected_unit_GT"]],
                 {
                   values[["selected_unit"]] <- input[["unit_selected_unit_GT"]]
                 })
    
    observeEvent(input[["unit_selected_cohort"]],
                 {
                   values[["selected_cohort"]] <- input[["unit_selected_cohort"]]
                 })
    
    output$selected_unit_pat <- renderText({
      paste("<b>",values$selected_unit,"</b><br>")
    })
    
    #### Tables ----------------------------------------------------------------
    output$unit_selected_cohort <- renderUI({
      shiny::validate(need(!is.null(values$gdb_present), "You need to load a gdb to use this tab!"))
      if (!is.null(values$gdb_present)) {
        selectInput("unit_selected_cohort", label = "Select cohort to use", choices = listCohort(rvc@gdb)[,1])
      }
    })
    
    output$unit_patient_table <- DT::renderDataTable({
      shiny::validate(need(!is.null(values$gdb_present), "You need to load a gdb to use this tab!"))
      unit <- unlist(strsplit(values$selected_unit, ","))
      shiny::validate(need(length(unit) == 1, "You need to select at least 1 unit!"))
      shiny::validate(need(!is.null(rvc@varSet[[unit]]),"No varSet could be found for the selected unit!"))
      
      dat <- rvc@gdb
      varPos <- data_unit_varInfo()
      
      if (!is.null(varSet) & unit %in% rvc@metadata$varSet$units) {
        VAR_ids <- unlist(strsplit(rvc@varSet[[unit]]@VAR_id, ","))
      } else if ("unit" %in% colnames(varPos)) {
        VAR_ids <- varPos[varPos$unit %in% unit, "VAR_id", drop = FALSE]
      } else if (!is.null(UNIT)) {
        rvb_data <- data_unit_rvb()[,c("unit", UNIT)]  
        VAR_ids <- varPos[varPos[[UNIT]] %in% rvb_data[rvb_data$unit == unit, UNIT],"VAR_id", drop = FALSE]
      }
      
      shiny::validate(need(!exists(VAR_ids), "No VAR_ids were found for this unit, do the VAR_ids in the annotation table match those in the assay?"))
      
      GT <- getGT(dat, VAR_id = VAR_ids, cohort = values$selected_cohort)
      GT <- flipToMinor(GT)
      
      GT <- GT[,colSums(SummarizedExperiment::assays(GT)$GT, na.rm = TRUE) != 0]
      samples <- colnames(GT)
      
      cohortTable <- getCohort(dat, values$selected_cohort)
      cohortTable <- cohortTable[cohortTable[[ID]] %in% samples,]
      
      cohortTable
    })
    
    output$unit_GT_assay <- DT::renderDataTable({
      shiny::validate(need(!is.null(values$gdb_present), "You need to load a gdb to use this tab!"))
      unit <- unlist(strsplit(values$selected_unit, ","))
      shiny::validate(need(length(unit) == 1, "You need to select at least 1 unit!"))
      
      dat <- rvc@gdb
      varPos <- data_unit_varInfo()
      
      if (!is.null(varSet) & unit %in% rvc@metadata$varSet$units) {
        VAR_ids <- unlist(strsplit(rvc@varSet[[unit]]@VAR_id, ","))
      } else if ("unit" %in% colnames(varPos)) {
        VAR_ids <- varPos[varPos$unit %in% unit, "VAR_id", drop = FALSE]
      } else if (!is.null(UNIT)) {
        rvb_data <- data_unit_rvb()[,c("unit", UNIT)]  
        VAR_ids <- varPos[varPos[[UNIT]] %in% rvb_data[rvb_data$unit == unit, UNIT],"VAR_id", drop = FALSE]
      }
      
      shiny::validate(need(!exists(VAR_ids), "No VAR_ids were found for this unit, do the VAR_ids in the annotation table match those in the assay?"))
      
      GT <- getGT(dat, VAR_id = VAR_ids, cohort = values$selected_cohort)
      GT <- flipToMinor(GT)
      
      GT <- GT[,colSums(SummarizedExperiment::assays(GT)$GT, na.rm = TRUE) != 0]
      assays <- data.frame(t(assays(GT)$GT))
      
      assays
    })
    
    # Association results ------------------------------------------------------
    ## Side panel --------------------------------------------------------------
    ### Reactives --------------------------------------------------------------
    assoc_data_main <- reactive({
      if (length(assoc_filters_main[["filters"]]) > 0) {
        dat <- subSet(rvc, analysis = values[["assoc_analysis_main"]], filters = assoc_filters_main[["filters"]])
      } else {
        dat <- data.frame(rvc@rvatResults[[values[["assoc_analysis_main"]]]])
      }
    })
    
    ### Filters main Assoc test (populated in observer)
    assoc_filters_main <- reactiveValues(
      filters = list(),
      add_filters = list()
    )
    
    ### Summary data -----------------------------------------------------------
    output$assoc_analysis_info <- renderText({
      dat <- assoc_data_main()
      
      if (identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "singlevarResult")) {
        dat$unit <- dat$VAR_id
      }
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

    ### Set custom significance threshold --------------------------------------
    observeEvent(input[["assoc_set_custom_threshold"]],
                 {
                   values[["assoc_custom_threshold"]] <- input[["assoc_custom_threshold"]]
                 })
    
    ### Set hover fields assoc -------------------------------------------------
    observeEvent(input[["assoc_set_hover_fields"]],
                 {
                   values[["assoc_hover_fields"]] <- input[["assoc_set_hover_fields"]]
                 }
    )
    
    output$assoc_set_hover_fields_ui <- renderUI({
      if(input$assoc_edit_hover_variables) {
        checkboxGroupInput("assoc_set_hover_fields",
                           label = "Fields",
                           choices = setdiff(metadata(rvc)[[values[["assoc_analysis_main"]]]][["columns"]], c("unit", "alias")),
                           selected = values[["assoc_hover_fields"]],
                           inline = TRUE)
        
      }
    })
    
    ### Addvariables -----------------------------------------------------------
    ### Selection dataset 
    output$assoc_analysis_ui <- renderUI({
      if(length(metadata(rvc)$analyses) > 1) {
        selectInput(
          inputId = "assoc_analysis",
          label = "Analysis",
          choices = metadata(rvc)$analyses,
          selected = metadata(rvc)$analyses[1]
        )
      }
    })
    
    # Set values[["assoc_analysis_main"]] after selection dataset
    observeEvent(input[["assoc_analysis"]],ignoreInit = TRUE,
                 {
                   values[["assoc_analysis_main"]] <- input[["assoc_analysis"]]
                   assoc_filters_main[["add_filters"]] <- c()
                   assoc_filters_compare[["add_filters"]] <- c()
                   
                   if (sum(lengths(metadata(rvc)[[values[["assoc_analysis_main"]]]][assoc_core_filter_variables]) > 1) == 0) {
                     values[["assoc_core_variables_main"]] <- core_filter_variables_default_assoc
                   } else {
                     values[["assoc_core_variables_main"]] <-
                       names(metadata(rvc)[[values[["assoc_analysis_main"]]]][assoc_core_filter_variables][lengths(metadata(rvc)[[values[["assoc_analysis_main"]]]][assoc_core_filter_variables]) > 1])
                   }
                   
                   values[["assoc_hover_fields"]] <-
                     values[["assoc_hover_fields"]][values[["assoc_hover_fields"]] %in% metadata(rvc)[[values[["assoc_analysis_main"]]]][["columns"]]]
                   updateSelectInput(session,
                                     inputId = "assoc_compare_variable",
                                     choices = names(metadata(rvc)[[values[["assoc_analysis_main"]]]]$column_types)[metadata(rvc)[[values[["assoc_analysis_main"]]]]$column_types %in% c("numeric", "integer")],
                                     selected = "P"
                   )
                 })
    
    ### Update filters 
    observe({
      
      ## Extract inputs core variables
      req(values[["assoc_core_variables_main"]])
      filters_core <- lapply(
        values[["assoc_core_variables_main"]],
        FUN = function(x) {
          type <-
            rvc@metadata[[values[["assoc_analysis_main"]]]]$column_types[x]
          req(input[[sprintf("%s_assoc_main", x)]])
          return(list(variable = x, type = type, string = input[[sprintf("%s_assoc_main", x)]], negate = FALSE, keepNA = FALSE))
        }
      )
      names(filters_core) <- values[["assoc_core_variables_main"]]
      
      ## Extract inputs additional variables
      filters_add <- lapply(
        assoc_filters_main[["add_filters"]],
        FUN = function(x) {
          type <-
            rvc@metadata[[values[["assoc_analysis_main"]]]]$column_types[x]
          if (type == "character") {
            req(input[[sprintf("%s_assoc_main", x)]])
            return(list(variable = x, type = type, string = input[[sprintf("%s_assoc_main", x)]], negate = FALSE, keepNA = FALSE))
            
          } else if (type == "logical") {
            req(input[[sprintf("%s_assoc_main", x)]])
            return(list(variable = x, type = type, bool = as.logical(convert_vec_logical[input[[sprintf("%s_assoc_main", x)]]]), 
                        keepNA = FALSE))
          } else if (type %in% c("numeric", "integer")) {
            req(input[[sprintf("%s_assoc_main_min", x)]])
            req(input[[sprintf("%s_assoc_main_max", x)]])
            return(list(
              variable = x,
              type = type,
              min = as.numeric(input[[sprintf("%s_assoc_main_min", x)]]),
              max = as.numeric(input[[sprintf("%s_assoc_main_max", x)]]),
              keepNA = FALSE
            ))
          }
        }
      )
      names(filters_add) <- assoc_filters_main[["add_filters"]]
      assoc_filters_main[["filters"]] <- c(filters_core, filters_add)
    })
    
    ### Filter variables 
    observeEvent(input[["assoc_addvariables_add_main"]],
                 {
                   assoc_filters_main[["add_filters"]] <- c(assoc_filters_main[["add_filters"]], input[["assoc_addvariables_select_main"]])
                   updateSelectInput(session,
                                     inputId = "assoc_addvariables_select_main",
                                     choices = metadata(rvc)[[values[["assoc_analysis_main"]]]][["columns"]][!metadata(rvc)[[default_assoc]][["columns"]] %in% c(assoc_filters_main[["add_filters"]], assoc_core_filter_variables)])
                 })
    
    ### Inputs additional variables - Observe
    observeEvent(values[["assoc_core_variables_main"]],
                 {
                   output$assoc_core_variables_main_ui <- renderUI({
                     lapply(
                       values[["assoc_core_variables_main"]],
                       FUN = function(x) {
                         choices <-
                           metadata(rvc)[[values[["assoc_analysis_main"]]]][[x]]
                         return(selectInput(
                           inputId = sprintf("%s_assoc_main", x),
                           label = x,
                           choices = choices,
                           selected = choices[1]
                         ))
                       }
                     )
                   })
                 })
    
    ### Inputs additional filter variables - Function
    output$assoc_addvariables_main_ui <- renderUI({
      if (length(assoc_filters_main[["add_filters"]]) > 0) {
        lapply(
          assoc_filters_main[["add_filters"]],
          FUN = function(x) {
            type <-
              rvc@metadata[[values[["assoc_analysis_main"]]]]$column_types[x]
            dat <- isolate(assoc_filters_main[["filters"]])
            if (type %in% c("numeric", "integer")) {
              return(fluidRow(column(
                6,
                numericInput(
                  inputId = sprintf("%s_assoc_main_min", x),
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
                  inputId = sprintf("%s_assoc_main_max", x),
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
                  rvc@rvatResults %>% dplyr::tbl(values[["assoc_analysis_main"]]) %>% dplyr::pull(x) %>% unique()
              } else {
                choices <- unique(rvc@rvatResults[[values[["assoc_analysis_main"]]]][[x]])
              }
              
              return(
                selectInput(
                  inputId = sprintf("%s_assoc_main", x),
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
    
    ### Reset additional variables ---------------------------------------------
    observeEvent(input[["assoc_addvariables_reset_main"]],
                 {
                   assoc_filters_main[["add_filters"]] <- c()
                   updateSelectInput(session,
                                     inputId = "assoc_addvariables_select_main",
                                     choices = c(metadata(rvc)[[values[["assoc_analysis_main"]]]][["columns"]][!metadata(rvc)[[values[["assoc_analysis_main"]]]][["columns"]] %in% assoc_core_filter_variables]),
                                     selected = "P"
                   )
                 }
    )
    
    ### Reset additional variables button
    output$assoc_addvariables_main_ui_reset <- renderUI({
      if (length(assoc_filters_main[["add_filters"]]) > 0) {
        actionButton("assoc_addvariables_reset_main",
                     "Reset additional variables")
      }
    })

    ## Tabs --------------------------------------------------------------------
    ### Manhattan --------------------------------------------------------------
    #### Manhattan plot using manhattanly::manhattanly -------------------------
    output$assoc_manhattan <- plotly::renderPlotly({
      dat <- assoc_data_main() 
      shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "rvbResult") | 
                             identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "singlevarResult"),
                           "You need an rvbResult or singlevarResult object to use this tab!"))
      shiny::validate(need(nrow(dat) > 0, "No data available for this combination of parameters"))
      shiny::validate(need(all(c("CHROM", "POS") %in% rvc@metadata[[values[["assoc_analysis_main"]]]]$columns) ||
                             (!is.null(CHROM) && !is.null(POS)), 
                           "'CHROM' and 'POS' columns are required for manhattan plot.
If these columns are present, but have other column names, this can be specified
using the `CHROM` and `POS` parameters. 
                           "))
      
      if(is.null(CHROM)) CHROM <- "CHROM"
      if(is.null(POS)) POS <- "POS"
      #dat <- assoc_data_main()
      dat$effect <- round(dat$effect, 2)
      dat[["CHR"]] <- as.numeric(dat[[CHROM]])
      dat[["BP"]] <- as.numeric(dat[[POS]])
      
      if (identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "singlevarResult")) {
        dat$unit <- dat$VAR_id
      }
      dat <- dat[,c("unit", "CHR", "BP", input[["assoc_display_variable"]], "P", values[["assoc_hover_fields"]]),drop=FALSE]
      dat <- dat[complete.cases(dat[,c("CHR", "BP")]),]
      col <- paste0(dat[[values[["assoc_hover_fields"]][1]]], "\n")
      
      if(length(values[["assoc_hover_fields"]]) > 1) {
        for(i in 2:(length(values[["assoc_hover_fields"]]))){
          col <- paste0(col, values[["assoc_hover_fields"]][i], ": ", dat[[values[["assoc_hover_fields"]][i]]], "\n")
        } 
      }
      dat[[values[["assoc_hover_fields"]][[1]]]] <- col
      dat <- manhattanr(dat, gene = input[["assoc_display_variable"]], 
                        annotation1 = values[["assoc_hover_fields"]][[1]])
      if(length(values[["assoc_custom_threshold"]]) > 0 && input[["assoc_use_custom_threshold"]]) {
        threshold <- values[["assoc_custom_threshold"]]
      } else {
        threshold <-  0.05/nrow(dat$data)
      }
      dat_sig = dat$data[dat$data[["P"]] < threshold,,drop=FALSE]
      
      if(nrow(dat_sig) > 0) {
        plot <- manhattanly(dat,
                            suggestiveline = FALSE, 
                            genomewideline = -log10(threshold),
                            title = "") %>% #paste0(input[["cohort"]], ";", input[["varSetName"]], "\n", input[["test"]], ";", input[["covar"]])
          plotly::add_annotations(x = dat_sig$pos,
                                  y = dat_sig$logp,
                                  text = dat_sig[[input[["assoc_display_variable"]]]],
                                  xref = "x",
                                  yref = "y",
                                  showarrow = TRUE,
                                  arrowhead=2,
                                  ax = 10,
                                  ay = -10) #%>%
        # plotly::layout(margin=list(l=10, r=10, b=10, t=60))
      } else {
        plot <- manhattanly(dat,
                            suggestiveline = FALSE, 
                            genomewideline = -log10(threshold),
                            title = "") #paste0(input[["cohort"]], ";", input[["varSetName"]], "\n", input[["test"]], ";", input[["covar"]], "\n\n\n")
      }
      
      #plotly::toWebGL(plotly::style(plot, hoveron = NULL))
      plotly::toWebGL(plot)
    })
    
    #### Toptable manhattan ----------------------------------------------------
    observeEvent(input[["assoc_set_variables_toptable_manh"]],
                 {
                   values[["assoc_variables_table"]] <- input[["assoc_set_variables_toptable_manh"]]
                 }
    )
    
    output$assoc_edit_variables_toptable_manh_ui <- renderUI({
      if(input$assoc_edit_variables_toptable_manh) {
        checkboxGroupInput("assoc_set_variables_toptable_manh",
                           label = "Fields",
                           choices = metadata(rvc)[[values[["assoc_analysis_main"]]]][["columns"]],
                           selected = c("unit", "P", "effect", "effectSE", "effectCIlower", "effectCIupper"),
                           inline = TRUE)
      }
    })
    
    output$assoc_toptable_manh <- DT::renderDataTable(
      {
        dat <- assoc_data_main()
        if(length(values[["assoc_custom_threshold"]]) > 0 && input[["assoc_use_custom_threshold"]]) {
          threshold <- values[["assoc_custom_threshold"]]
        } else {
          threshold <-  0.05/nrow(dat)
        }
        
        if (!("unit" %in% values[["assoc_variables_table"]])) {
          dat <- dat[dat$P < threshold, c("unit", values[["assoc_variables_table"]]), drop=FALSE] %>% dplyr::arrange(P)
        } else {
          dat <- dat[dat$P < threshold,values[["assoc_variables_table"]],drop=FALSE] %>% dplyr::arrange(P)
        }
        
        for (i in 1:nrow(dat)) {
          if(!is.null(UNIT)) {
            if ((!is.null(varSet) & dat[i,"unit"] %in% rvc@metadata$varSet$units) | 
                ("unit" %in% rvc@metadata$varInfo & dat[i,"unit"] %in% rvc@varInfo$unit) | 
                (!is.null(UNIT) & dat[i,"unit"] %in% rvc@varInfo[[UNIT]])) {
              dat[i,"unit"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_manh_assoc&quot;,this.id)">', dat[i,"unit"], '</a>')
            }
          } else {
            if ((!is.null(varSet) & dat[i,"unit"] %in% rvc@metadata$varSet$units) | 
                ("unit" %in% rvc@metadata$varInfo & dat[i,"unit"] %in% rvc@varInfo$unit)) {
              dat[i,"unit"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_manh_assoc&quot;,this.id)">', dat[i,"unit"], '</a>')
            }
          }
        }
        
        number_links <- nrow(dat[substr(dat[["unit"]],1,5) == "<a id",, drop= FALSE])
        if (number_links == 0) {
          showNotification("None of the units could be linked to its variants, so no links to the `Unit viewer` have been made")
        }

        table <- DT::datatable(dat, escape = FALSE, selection = "none")
      })
    
    ### Make buttons in the table to connect to other tabs
    observeEvent(input$select_button_manh_assoc, {
      dat <- assoc_data_main()
      
      if(length(values[["assoc_custom_threshold"]]) > 0 && input[["assoc_use_custom_threshold"]]) {
        threshold <- values[["assoc_custom_threshold"]]
      } else {
        threshold <-  0.05/nrow(dat)
      }
      
      if (!("unit" %in% values[["assoc_variables_table"]])) {
        dat <- dat[dat$P < threshold, c("unit", values[["assoc_variables_table"]]), drop=FALSE] %>% dplyr::arrange(P)
      } else {
        dat <- dat[dat$P < threshold,values[["assoc_variables_table"]],drop=FALSE] %>% dplyr::arrange(P)
      }
      selectedRow <- as.numeric(strsplit(input$select_button_manh_assoc, "_")[[1]][2])
      values$selected_unit <- dat[selectedRow, "unit"]
      updateNavbarPage(session, 
                       "main_menu",
                       "unit_viewer")
    })
    
    
    ### QQplot -----------------------------------------------------------------
    #### qqplot using manhattanly::qqly ----------------------------------------
    output$assoc_qqplot <- plotly::renderPlotly({
      shiny::validate(need(nrow(assoc_data_main()) > 0, "No data available for this combination of parameters"))
      shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "rvbResult") | 
                             identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "singlevarResult"),
                           "You need an rvbResult or singlevarResult object to use this tab!"))

      dat <- assoc_data_main()
      if(length(values[["assoc_custom_threshold"]]) > 0 && input[["assoc_use_custom_threshold"]]) {
        threshold <- values[["assoc_custom_threshold"]]
      } else {
        threshold <- 0.05/nrow(dat)
      }
      
      set.seed(1)
      dat_sig <- dplyr::select(dat, P, input[["assoc_display_variable"]])
      dat_sig$logp <- -log10(dat_sig$P)
      dat_sig <- dat_sig %>% dplyr::arrange(logp)
      dat_sig$null <- sort(-log10(runif(nrow(dat))))
      dat_sig = dat_sig[dat_sig[["logp"]] > -log10(threshold),,drop=FALSE]
      
      if (nrow(dat_sig) > 0) {
        plot <- qqly(dat, threshold = threshold, xlab = "-log10(NULL)", ylab = "-log10(Obs)", title = "") %>%
          plotly::add_annotations(x = dat_sig$null,
                                  y = dat_sig$logp,
                                  text = dat_sig[[input[["assoc_display_variable"]]]],
                                  xref = "x",
                                  yref = "y",
                                  showarrow = TRUE,
                                  arrowhead=2,
                                  ax = 10,
                                  ay = -10)
      } else {
        plot <- qqly(dat, threshold = threshold, xlab = "-log10(NULL)", ylab = "-log10(Obs)", title = "")
      }

      plotly::toWebGL(plot)
    })
    
    #### Toptable qqplot -------------------------------------------------------
    observeEvent(input[["assoc_set_variables_toptable_qq"]],
                 {
                   values[["assoc_variables_table"]] <- input[["assoc_set_variables_toptable_qq"]]
                 }
    )
    
    output$assoc_edit_variables_toptable_qq_ui <- renderUI({
      if(input$assoc_edit_variables_toptable_qq) {
        checkboxGroupInput("assoc_set_variables_toptable_qq",
                           label = "Fields",
                           choices = metadata(rvc)[[values[["assoc_analysis_main"]]]][["columns"]],
                           selected = c("unit", "P", "effect", "effectSE", "effectCIlower", "effectCIupper"),
                           inline = TRUE)
      }
    })
    
    output$assoc_toptable_qq <- DT::renderDataTable(
      {
        dat <- assoc_data_main()
        if(length(values[["assoc_custom_threshold"]]) > 0 && input[["assoc_use_custom_threshold"]]) {
          threshold <- values[["assoc_custom_threshold"]]
        } else {
          threshold <-  0.05/nrow(dat)
        }
        
        if (!("unit" %in% values[["assoc_variables_table"]])) {
          dat <- dat[dat$P < threshold, c("unit", values[["assoc_variables_table"]]), drop=FALSE] %>% dplyr::arrange(P)
        } else {
          dat <- dat[dat$P < threshold,values[["assoc_variables_table"]],drop=FALSE] %>% dplyr::arrange(P)
        }
        
        for (i in 1:nrow(dat)) {
          if(!is.null(UNIT)) {
            if ((!is.null(varSet) & dat[i,"unit"] %in% rvc@metadata$varSet$units) | 
                ("unit" %in% rvc@metadata$varInfo & dat[i,"unit"] %in% rvc@varInfo$unit) | 
                (!is.null(UNIT) & dat[i,"unit"] %in% rvc@varInfo[[UNIT]])) {
              dat[i,"unit"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_qq_assoc&quot;,this.id)">', dat[i,"unit"], '</a>')
            }
          } else {
            if ((!is.null(varSet) & dat[i,"unit"] %in% rvc@metadata$varSet$units) | 
                ("unit" %in% rvc@metadata$varInfo & dat[i,"unit"] %in% rvc@varInfo$unit)) {
              dat[i,"unit"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_qq_assoc&quot;,this.id)">', dat[i,"unit"], '</a>')
            }
          }
        }
        
        number_links <- nrow(dat[substr(dat[["unit"]],1,5) == "<a id",, drop= FALSE])
        if (number_links == 0) {
          showNotification("None of the units could be linked to its variants, so no links to the `Unit viewer` have been made")
        }
        
        table <- DT::datatable(dat, escape = FALSE, selection = "none")
      })
    
    ### Make buttons in the table to connect to other tabs
    observeEvent(input$select_button_qq_assoc, {
      dat <- assoc_data_main()
      
      if(length(values[["assoc_custom_threshold"]]) > 0 && input[["assoc_use_custom_threshold"]]) {
        threshold <- values[["assoc_custom_threshold"]]
      } else {
        threshold <-  0.05/nrow(dat)
      }
      
      if (!("unit" %in% values[["assoc_variables_table"]])) {
        dat <- dat[dat$P < threshold, c("unit", values[["assoc_variables_table"]]), drop=FALSE] %>% dplyr::arrange(P)
      } else {
        dat <- dat[dat$P < threshold,values[["assoc_variables_table"]],drop=FALSE] %>% dplyr::arrange(P)
      }
      selectedRow <- as.numeric(strsplit(input$select_button_qq_assoc, "_")[[1]][2])
      values$selected_unit <- dat[selectedRow, "unit"]
      updateNavbarPage(session, 
                       "main_menu",
                       "unit_viewer")
    })
    
    ### Forest plot ------------------------------------------------------------
    #### Forest plot -----------------------------------------------------------
    output$assoc_forest_plot <- renderPlot({
      shiny::validate(need(length(input[["assoc_selected_units_forest"]]) > 0, "You need to select at least 1 unit to plot first!"))
      shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "rvbResult") |
                             identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "singlevarResult"),
                           "You need an rvbResult or singlevarResult object to use this tab!"))
      
      dat <- assoc_data_main()
      if (identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "singlevarResult")) {
        dat$unit <- dat$VAR_id
      }
      unit <- unlist(strsplit(input[["assoc_selected_units_forest"]],","))
      unit_keep <- unit[unit %in% dat[[input[["assoc_display_variable"]]]] ]
      unit_drop <- unit[!(unit %in% dat[[input[["assoc_display_variable"]]]]) ]
      
      if (length(unit_drop) > 0) {
        showNotification(paste0("The units ", paste(unit_drop, collapse = ","), " are not in the dataset and have been dropped. Check if you selected the correct display variable!"))
      }
      
      forestData <- data.frame(unit = dat[[input[["assoc_display_variable"]]]],
                               effectSE = dat$effectSE,
                               effect = dat$effect,
                               effectCIlower = dat$effectCIlower,
                               effectCIupper = dat$effectCIupper)
      forestData <- forestData[forestData$unit %in% unit_keep,]
      forestData <- forestData[!is.na(forestData$effectSE),]
      forestData <- rbind(c(rep(NA,ncol(forestData))), forestData)
      
      labelData <- data.frame(unit = dat[[input[["assoc_display_variable"]]]], P = round(dat$P,5), effectSE = round(dat$effectSE,5))
      labelData <- labelData[labelData$unit %in% unit_keep,]
      labelData <- labelData[!is.na(labelData$effectSE),]
      
      labelData <- cbind(data.frame(Index = c(1:(nrow(labelData)))), labelData)
      labelData <- rbind(colnames(labelData), labelData)
      
      clip = c(floor(min(forestData$effect, na.rm = TRUE)), ceiling(max(forestData$effect, na.rm = TRUE)))
      
      forestplot::forestplot(forestData,
                             mean = effect,
                             lower = effectCIlower,
                             upper = effectCIupper,
                             clip = clip,
                             boxsize = 0.15,
                             xlab = "Effect size",
                             vertices = TRUE,
                             labeltext = labelData,
                             is.summary = c(TRUE, rep(FALSE, nrow(forestData))),
                             hrzl_lines = grid::gpar(lty = 2, col = "#444444"),
                             colgap = grid::unit(5,"mm"),
                             txt_gp = forestplot::fpTxtGp(label = grid::gpar(cex = 1),
                                                          ticks = grid::gpar(cex = 0.5)),
                             zero = 1,
                             cex = 2,
                             lineheight = "auto")
    })
    
    #### Forest plot table with columns with more than one unique value --------
    forest_table_data <- reactive({
      dat <- assoc_data_main()
      if (identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "singlevarResult")) {
        dat$unit <- dat$VAR_id
      }
      
      unit <- unlist(strsplit(input[["assoc_selected_units_forest"]],","))
      unit_keep <- unit[unit %in% dat[[input[["assoc_display_variable"]]]] ]
      
      tableData = data.frame(matrix(NA, ncol=1, nrow=nrow(dat)))[-1]
      for (column in colnames(dat)) {
        if (dplyr::n_distinct(data.frame(dat[column])) > 1) {
          tableData = cbind(tableData, dat[column])
        }
      } 
      
      if (!(input[["assoc_display_variable"]] %in% colnames(tableData))) {
        cbind(dat[[input[["assoc_display_variable"]]]], tableData)
      }
      
      if (!("unit" %in% colnames(tableData))) {
        cbind(dat[["unit"]], tableData)
      }
      
      tableData <- tableData[tableData[[input[["assoc_display_variable"]]]] %in% unit_keep,]
      tableData <- tableData[!is.na(tableData$effectSE),]
      tableData <- cbind(data.frame(Index = c(1:(nrow(tableData)))), tableData)
    })
    
    output$assoc_forest_table <- DT::renderDataTable({
      unit <- unlist(strsplit(input[["assoc_selected_units_forest"]],","))
      shiny::validate(need(length(unit) > 0, "You need to select at least 1 unit to plot first!"))
      
      tableData <- forest_table_data()
      for (i in 1:nrow(tableData)) {
        if (!is.null(UNIT)) {
          if ((!is.null(varSet) & tableData[i,"unit"] %in% rvc@metadata$varSet$units) | 
              ("unit" %in% rvc@metadata$varInfo & tableData[i,"unit"] %in% rvc@varInfo$unit) | 
              (!is.null(UNIT) & tableData[i,"unit"] %in% rvc@varInfo[[UNIT]])) {
            tableData[i,"unit"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_forest_assoc&quot;,this.id)">', tableData[i,"unit"], '</a>')
          }
        } else {
          if ((!is.null(varSet) & tableData[i,"unit"] %in% rvc@metadata$varSet$units) | 
              ("unit" %in% rvc@metadata$varInfo & tableData[i,"unit"] %in% rvc@varInfo$unit)) {
            tableData[i,"unit"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_forest_assoc&quot;,this.id)">', tableData[i,"unit"], '</a>')
          }
        }
      }
      
      number_links <- nrow(tableData[substr(tableData[["unit"]],1,5) == "<a id",, drop= FALSE])
      if (number_links == 0) {
        showNotification("None of the units could be linked to its variants, so no links to the `Unit viewer` have been made")
      }
      
      table <- DT::datatable(tableData, filter = "none", escape = FALSE, selection = "none")
    })
    
    ### Make buttons in the table to connect to other tabs
    observeEvent(input$select_button_forest_assoc, {
      unit <- unlist(strsplit(input[["assoc_selected_units_forest"]],","))
      shiny::validate(need(length(unit) > 0, "You need to select at least 1 unit to plot first!"))
      
      dat <- forest_table_data()
      
      if(length(values[["assoc_custom_threshold"]]) > 0 && input[["assoc_use_custom_threshold"]]) {
        threshold <- values[["assoc_custom_threshold"]]
      } else {
        threshold <-  0.05/nrow(dat)
      }
      
      selectedRow <- as.numeric(strsplit(input$select_button_forest_assoc, "_")[[1]][2])
      values$selected_unit <- dat[selectedRow, "unit"]
      updateNavbarPage(session, 
                       "main_menu",
                       "unit_viewer")
    })
    
    #### Forest plot unique column values --------------------------------------
    output$assoc_forest_unique_columns <- renderText({
      unit <- unlist(strsplit(input[["assoc_selected_units_forest"]],","))
      shiny::validate(need(length(unit) > 0, "You need to select at least 1 unit to plot first!"))
      
      dat <- assoc_data_main()
      textData <- list()
      for (column in colnames(dat)) {
        if (dplyr::n_distinct(data.frame(dat[column])) == 1) {
          textData[[column]] = dat[1,column]
        }
      } 
      
      text <- ""
      for (column in names(textData)) {
        text <- paste(text, "<b>", column, "</b>: ", textData[[column]], "<br>")
      }
      text
    })

    ### Compare ----------------------------------------------------------------
    ### A subset of the data fitting the filters the user has defined
    assoc_data_compare <- eventReactive(input[["assoc_compare_button"]],{
      dat <- subSet(rvc, analysis = values[["assoc_analysis_main"]], filters = assoc_filters_compare[["filters"]])
      dat
    })
    
    ### The filters the user has defined
    assoc_filters_compare <- reactiveValues(
      filters = list(),
      add_filters = list()
    )
    
    ### The user adds a filter and assoc_filters_compare is updated
    observeEvent(input[["assoc_addvariables_add_compare"]],
                 {
                   assoc_filters_compare[["add_filters"]] <- c(assoc_filters_compare[["add_filters"]], input[["assoc_addvariables_select_compare"]])
                   updateSelectInput(session,
                                     inputId = "assoc_addvariables_select_compare",
                                     choices = metadata(rvc)[[values[["assoc_analysis_main"]]]][["columns"]][!metadata(rvc)[[default_assoc]][["columns"]] %in% c(assoc_filters_compare[["add_filters"]], assoc_core_filter_variables)])
                 })
    
    ### The choices of variables you can compare between
    observeEvent(values[["assoc_core_variables_main"]],
                 {
                   output$assoc_core_variables_compare_ui <- renderUI({
                     lapply(
                       values[["assoc_core_variables_main"]],
                       FUN = function(x) {
                         choices <-
                           metadata(rvc)[[values[["assoc_analysis_main"]]]][[x]]
                         return(selectInput(
                           inputId = sprintf("%s_assoc_compare", x),
                           label = x,
                           choices = choices,
                           selected = choices[1]
                         ))
                       }
                     )
                   })
                 })
    
    ### Add filters
    observe({
      
      ## Extract inputs core variables
      req(values[["assoc_core_variables_main"]])
      filters_core <- lapply(
        values[["assoc_core_variables_main"]],
        FUN = function(x) {
          type <-
            rvc@metadata[[values[["assoc_analysis_main"]]]]$column_types[x]
          req(input[[sprintf("%s_assoc_compare", x)]])
          return(list(variable = x, type = type, string = input[[sprintf("%s_assoc_compare", x)]], negate = FALSE, keepNA = FALSE))
        }
      )
      names(filters_core) <- values[["assoc_core_variables_main"]]
      
      ## Extract inputs additional variables
      filters_add <- lapply(
        assoc_filters_compare[["add_filters"]],
        FUN = function(x) {
          type <-
            rvc@metadata[[values[["assoc_analysis_main"]]]]$column_types[x]
          if (type == "character") {
            req(input[[sprintf("%s_assoc_compare", x)]])
            return(list(variable = x, type = type, string = input[[sprintf("%s_assoc_compare", x)]], negate = FALSE, keepNA = FALSE))
            
          } else if (type == "logical") {
            req(input[[sprintf("%s_assoc_compare", x)]])
            return(list(variable = x, type = type, bool = as.logical(convert_vec_logical[input[[sprintf("%s_assoc_compare", x)]]]), 
                        keepNA = FALSE))
          } else if (type %in% c("numeric", "integer")) {
            req(input[[sprintf("%s_assoc_compare_min", x)]])
            req(input[[sprintf("%s_assoc_compare_max", x)]])
            return(list(
              variable = x,
              type = type,
              min = as.numeric(input[[sprintf("%s_assoc_compare_min", x)]]),
              max = as.numeric(input[[sprintf("%s_assoc_compare_max", x)]]),
              keepNA = FALSE
            ))
          }
        }
      )
      names(filters_add) <- assoc_filters_compare[["add_filters"]]
      assoc_filters_compare[["filters"]] <- c(filters_core, filters_add)
    })
    
    ### The user defines the precise filter using fields and buttons that appear
    output$assoc_addvariables_compare_ui <- renderUI({
      if (length(assoc_filters_compare[["add_filters"]]) > 0) {
        lapply(
          assoc_filters_compare[["add_filters"]],
          FUN = function(x) {
            type <-
              rvc@metadata[[values[["assoc_analysis_main"]]]]$column_types[x]
            dat <- isolate(assoc_filters_compare[["filters"]])
            if (type %in% c("numeric", "integer")) {
              return(fluidRow(column(
                6,
                numericInput(
                  inputId = sprintf("%s_assoc_compare_min", x),
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
                  inputId = sprintf("%s_assoc_compare_max", x),
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
                  rvc@rvatResults %>% dplyr::tbl(values[["assoc_analysis_main"]]) %>% dplyr::pull(x) %>% unique()
              } else {
                choices <- unique(rvc@rvatResults[[values[["assoc_analysis_main"]]]][[x]])
              }
              
              return(
                selectInput(
                  inputId = sprintf("%s_assoc_compare", x),
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
    
    ### Reset additional variables button - compare tab
    output$assoc_addvariables_compare_ui_reset <- renderUI({
      if (length(assoc_filters_compare[["add_filters"]]) > 0) {
        actionButton("assoc_addvariables_reset_compare",
                     "Reset additional variables")
      }
    })
    
    ### Reset additional variables based on what the user filled in - compare tab
    observeEvent(input[["assoc_addvariables_reset_compare"]],
                 {
                   assoc_filters_compare[["add_filters"]] <- c()
                   updateSelectInput(session,
                                     inputId = "assoc_addvariables_select_compare",
                                     choices = c(metadata(rvc)[[values[["assoc_analysis_main"]]]][["columns"]][!metadata(rvc)[[values[["assoc_analysis_main"]]]][["columns"]] %in% assoc_core_filter_variables]),
                                     selected = "P"
                   )
    })
    
    ### The actual comparison plot
    output$assoc_compare <- plotly::renderPlotly({
      shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "rvbResult") | 
                             identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "singlevarResult"),
                           "You need an rvbResult or singlevarResult object to use this tab!"))
      
      dat1 <- assoc_data_main()
      dat2 <- assoc_data_compare()
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
        compare_variable = input[["assoc_compare_variable"]],
        transformation = input[["assoc_compare_transformation"]],
        custom_threshold = values[["assoc_custom_threshold"]],
        use_custom_threshold = input[["assoc_use_custom_threshold"]],
        hover_fields = values[["assoc_hover_fields"]],
        show_labels = input[["assoc_compare_show_labels"]],
        display_variable = input[["assoc_display_variable"]]
      )
    })
    
    ### tableViewer ------------------------------------------------------------
    observeEvent(input[["assoc_set_variables_tableviewer"]],
                 {
                   values[["assoc_variables_table"]] <- input[["assoc_set_variables_tableviewer"]]
                 }
    )
    
    output$assoc_edit_variables_tableviewer_ui <- renderUI({
      if(input$assoc_edit_variables_tableviewer) {
        checkboxGroupInput("assoc_set_variables_tableviewer",
                           label = "Fields",
                           choices = metadata(rvc)[[values[["assoc_analysis_main"]]]][["columns"]],
                           selected = c("unit", "P", "effect", "effectSE", "effectCIlower", "effectCIupper"),
                           inline = TRUE)
        
      }
    })
    
    output$assoc_tableviewer <- DT::renderDT(
      {
        shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "rvbResult") | 
                               identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "singlevarResult"),
                             "You need an rvbResult or singlevarResult object to use this tab!"))
        
        dat <- as.data.frame(rvc@rvatResults[[values[["assoc_analysis_main"]]]])
        
        if (identical(checkClassrvatResult(rvc@rvatResults[[values[["assoc_analysis_main"]]]]), "singlevarResult")) {
          dat$unit <- dat$VAR_id
        }
        
        for (i in 1:nrow(dat)) {
          if (!is.null(UNIT)) {
            if ((!is.null(varSet) & dat[i,"unit"] %in% rvc@metadata$varSet$units) | 
                ("unit" %in% rvc@metadata$varInfo & dat[i,"unit"] %in% rvc@varInfo$unit) | 
                (!is.null(UNIT) & dat[i,"unit"] %in% rvc@varInfo[[UNIT]])) {
              dat[i,"unit"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_tableviewer_assoc&quot;,this.id)">', dat[i,"unit"], '</a>')
            }
          } else {
            if ((!is.null(varSet) & dat[i,"unit"] %in% rvc@metadata$varSet$units) | 
                ("unit" %in% rvc@metadata$varInfo & dat[i,"unit"] %in% rvc@varInfo$unit)) {
              dat[i,"unit"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_tableviewer_assoc&quot;,this.id)">', dat[i,"unit"], '</a>')
            }
          }
          
        }
        
        number_links <- nrow(dat[substr(dat[["unit"]],1,5) == "<a id",, drop= FALSE])
        if (number_links == 0) {
          showNotification("None of the units could be linked to its variants, so no links to the `Unit viewer` have been made")
        }
        
        if (!("unit" %in% values[["assoc_variables_table"]])) {
          DT::datatable(
            dat[,c("unit",values[["assoc_variables_table"]]),drop=FALSE], 
            filter = "top", escape = FALSE, selection = "none")
        } else {
          DT::datatable(
            dat[,values[["assoc_variables_table"]],drop=FALSE], 
            filter = "top", escape = FALSE, selection = "none")
        }
      })
    
    ### Make buttons in the table to connect to other tabs
    observeEvent(input$select_button_tableviewer_assoc, {
      dat <- as.data.frame(rvc@rvatResults[[values[["assoc_analysis_main"]]]])
      
      selectedRow <- as.numeric(strsplit(input$select_button_tableviewer_assoc, "_")[[1]][2])
      values$selected_unit <- dat[selectedRow, "unit"]
      updateNavbarPage(session, 
                       "main_menu",
                       "unit_viewer")
    })
    
    # Gene set analysis --------------------------------------------------------
    ## Side panel --------------------------------------------------------------
    ### Reactives --------------------------------------------------------------
    geneset_data_main <- reactive({
      if (length(geneset_filters_main[["filters"]]) > 0) {
        if (is.na(values[["geneset_analysis_main"]])) {
          values[["geneset_analysis_main"]] <- values[["assoc_analysis_main"]]
          
          if (length(assoc_filters_main[["filters"]]) > 0) {
            dat <- subSet(rvc, analysis = values[["geneset_analysis_main"]], filters = assoc_filters_main[["filters"]])
          } else {
            dat <- data.frame(rvc@rvatResults[[values[["geneset_analysis_main"]]]])
          }
        } else {
          if (checkClassrvatResult(rvc@rvatResults[[values[["geneset_analysis_main"]]]]) %in% c("rvbResult", "singlevarResult")) {
            if (length(assoc_filters_main[["filters"]]) > 0) {
              dat <- subSet(rvc, analysis = values[["geneset_analysis_main"]], filters = assoc_filters_main[["filters"]])
            } else {
              dat <- data.frame(rvc@rvatResults[[values[["geneset_analysis_main"]]]])
            }
          } else {
            dat <- subSet(rvc, analysis = values[["geneset_analysis_main"]], filters = geneset_filters_main[["filters"]])
          }

        }
      } else {
        if (is.na(values[["geneset_analysis_main"]])) {
          values[["geneset_analysis_main"]] <- values[["assoc_analysis_main"]]
        } 
        dat <- data.frame(rvc@rvatResults[[values[["geneset_analysis_main"]]]])
      }
      
      dat
    })
    
    ### Filters main Assoc test (populated in observer)
    geneset_filters_main <- reactiveValues(
      filters = list(),
      add_filters = list()
    )
    
    ### Summary data -----------------------------------------------------------
    output$geneset_analysis_info <- renderText({
      dat <- as.data.frame(rvc@rvatResults[[values[["geneset_analysis_main"]]]])
      sprintf(
        "<b>N genesets</b> = %s<br><b>N methods</b> = %s<br><b>N tests</b> = %s",
        dplyr::n_distinct(dat$geneSetName),
        dplyr::n_distinct(dat$method),
        dplyr::n_distinct(dat$test)
      )
    })
    
    ### Set custom significance threshold --------------------------------------
    observeEvent(input[["geneset_set_custom_threshold"]],
                 {
                   values[["geneset_custom_threshold"]] <- input[["geneset_custom_threshold"]]
                 })
    
    ### Set hover fields assoc -------------------------------------------------
    observeEvent(input[["geneset_set_hover_fields"]],
                 {
                   values[["geneset_hover_fields"]] <- input[["geneset_set_hover_fields"]]
                 }
    )
    
    output$geneset_set_hover_fields_ui <- renderUI({
      if(input$geneset_edit_hover_variables) {
        checkboxGroupInput("geneset_set_hover_fields",
                           label = "Fields",
                           choices = setdiff(metadata(rvc)[[values[["geneset_analysis_main"]]]][["columns"]], c("unit", "alias")),
                           selected = values[["geneset_hover_fields"]],
                           inline = TRUE)
        
      }
    })
    
    ### Addvariables -----------------------------------------------------------
    ### Selection dataset 
    output$geneset_analysis_ui <- renderUI({
      if(length(metadata(rvc)$analyses) > 1) {
        selectInput(
          inputId = "geneset_analysis",
          label = "Analysis",
          choices = metadata(rvc)$analyses,
          selected = metadata(rvc)$analyses[1]
        )
      }
    })
    
    # Set values[["geneset_analysis_main"]] after selection dataset
    observeEvent(input[["geneset_analysis"]],ignoreInit = TRUE,
                 {
                   values[["geneset_analysis_main"]] <- input[["geneset_analysis"]]
                   geneset_filters_main[["add_filters"]] <- c()
                   
                   if (checkClassrvatResult(rvc@rvatResults[[values[["geneset_analysis_main"]]]]) %in% c("rvbResult", "singlevarResult")) {
                     if (sum(lengths(metadata(rvc)[[values[["geneset_analysis_main"]]]][assoc_core_filter_variables]) > 1) == 0) {
                       values[["geneset_core_variables_main"]] <- core_filter_variables_default_assoc
                     } else {
                       values[["geneset_core_variables_main"]] <-
                         names(metadata(rvc)[[values[["geneset_analysis_main"]]]][geneset_core_filter_variables][lengths(metadata(rvc)[[values[["geneset_analysis_main"]]]][geneset_core_filter_variables]) > 1])
                     }
                   } else {
                     if (sum(lengths(metadata(rvc)[[values[["geneset_analysis_main"]]]][geneset_core_filter_variables]) > 1) == 0) {
                       values[["geneset_core_variables_main"]] <- core_filter_variables_default_geneset
                     } else {
                       values[["geneset_core_variables_main"]] <-
                         names(metadata(rvc)[[values[["geneset_analysis_main"]]]][geneset_core_filter_variables][lengths(metadata(rvc)[[values[["geneset_analysis_main"]]]][geneset_core_filter_variables]) > 1])
                     }
                   }
                   
                   values[["geneset_core_variables_main"]] <-
                     names(metadata(rvc)[[values[["geneset_analysis_main"]]]][geneset_core_filter_variables][lengths(metadata(rvc)[[values[["geneset_analysis_main"]]]][geneset_core_filter_variables]) > 1])
                   values[["geneset_hover_fields"]] <-
                     values[["geneset_hover_fields"]][values[["geneset_hover_fields"]] %in% metadata(rvc)[[values[["geneset_analysis_main"]]]][["columns"]]]
                 })
    
    ### Update filters 
    observe({
      
      ## Extract inputs core variables
      req(values[["geneset_core_variables_main"]])
      filters_core <- lapply(
        values[["geneset_core_variables_main"]],
        FUN = function(x) {
          type <-
            rvc@metadata[[values[["geneset_analysis_main"]]]]$column_types[x]
          req(input[[sprintf("%s_geneset_main", x)]])
          return(list(variable = x, type = type, string = input[[sprintf("%s_geneset_main", x)]], negate = FALSE, keepNA = FALSE))
        }
      )
      names(filters_core) <- values[["geneset_core_variables_main"]]
      
      ## Extract inputs additional variables
      filters_add <- lapply(
        geneset_filters_main[["add_filters"]],
        FUN = function(x) {
          type <-
            rvc@metadata[[values[["geneset_analysis_main"]]]]$column_types[x]
          if (type == "character") {
            req(input[[sprintf("%s_geneset_main", x)]])
            return(list(variable = x, type = type, string = input[[sprintf("%s_geneset_main", x)]], negate = FALSE, keepNA = FALSE))
            
          } else if (type == "logical") {
            req(input[[sprintf("%s_geneset_main", x)]])
            return(list(variable = x, type = type, bool = as.logical(convert_vec_logical[input[[sprintf("%s_geneset_main", x)]]]), 
                        keepNA = FALSE))
          } else if (type %in% c("numeric", "integer")) {
            req(input[[sprintf("%s_geneset_main_min", x)]])
            req(input[[sprintf("%s_geneset_main_max", x)]])
            return(list(
              variable = x,
              type = type,
              min = as.numeric(input[[sprintf("%s_geneset_main_min", x)]]),
              max = as.numeric(input[[sprintf("%s_geneset_main_max", x)]]),
              keepNA = FALSE
            ))
          }
        }
      )
      names(filters_add) <- geneset_filters_main[["add_filters"]]
      geneset_filters_main[["filters"]] <- c(filters_core, filters_add)
    })
    
    ### Filter variables 
    observeEvent(input[["geneset_addvariables_add_main"]],
                 {
                   geneset_filters_main[["add_filters"]] <- c(geneset_filters_main[["add_filters"]], input[["geneset_addvariables_select_main"]])
                   updateSelectInput(session,
                                     inputId = "geneset_addvariables_select_main",
                                     choices = metadata(rvc)[[values[["geneset_analysis_main"]]]][["columns"]][!metadata(rvc)[[default_geneset]][["columns"]] %in% c(geneset_filters_main[["add_filters"]], geneset_core_filter_variables)])
                 })
    
    ### Inputs additional variables - Observe
    observeEvent(values[["geneset_core_variables_main"]],
                 {
                   output$geneset_core_variables_main_ui <- renderUI({
                     lapply(
                       values[["geneset_core_variables_main"]],
                       FUN = function(x) {
                         choices <-
                           metadata(rvc)[[values[["geneset_analysis_main"]]]][[x]]
                         return(selectInput(
                           inputId = sprintf("%s_geneset_main", x),
                           label = x,
                           choices = choices,
                           selected = choices[1]
                         ))
                       }
                     )
                   })
                 })
    
    ### Inputs additional filter variables - Function
    output$geneset_addvariables_main_ui <- renderUI({
      if (length(geneset_filters_main[["add_filters"]]) > 0) {
        lapply(
          geneset_filters_main[["add_filters"]],
          FUN = function(x) {
            type <-
              rvc@metadata[[values[["geneset_analysis_main"]]]]$column_types[x]
            dat <- isolate(geneset_filters_main[["filters"]])
            if (type %in% c("numeric", "integer")) {
              return(fluidRow(column(
                6,
                numericInput(
                  inputId = sprintf("%s_geneset_main_min", x),
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
                  inputId = sprintf("%s_geneset_main_max", x),
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
                  rvc@rvatResults %>% dplyr::tbl(values[["geneset_analysis_main"]]) %>% dplyr::pull(x) %>% unique()
              } else {
                choices <- unique(rvc@rvatResults[[values[["geneset_analysis_main"]]]][[x]])
              }
              
              return(
                selectInput(
                  inputId = sprintf("%s_geneset_main", x),
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
    
    ### Reset additional variables
    observeEvent(input[["geneset_addvariables_reset_main"]],
                 {
                   geneset_filters_main[["add_filters"]] <- c()
                   updateSelectInput(session,
                                     inputId = "geneset_addvariables_select_main",
                                     choices = c(metadata(rvc)[[values[["geneset_analysis_main"]]]][["columns"]][!metadata(rvc)[[values[["geneset_analysis_main"]]]][["columns"]] %in% geneset_core_filter_variables]),
                                     selected = "P"
                   )
                 }
    )
    
    ### Reset additional variables button
    output$geneset_addvariables_main_ui_reset <- renderUI({
      if (length(geneset_filters_main[["add_filters"]]) > 0) {
        actionButton("geneset_addvariables_reset_main",
                     "Reset additional variables")
      }
    })
    
    
    ## Tabs --------------------------------------------------------------------
    ### Manhattan --------------------------------------------------------------
    output$selected_geneset <- renderText({
      paste("<b>",values[["selected_geneset"]],"</b><br>")
    })
    
    observeEvent(input[["geneset_selected_units_manhattan"]],
                 {
                   values[["selected_geneset"]] <- input[["geneset_selected_units_manhattan"]]
                 })
    
    #### Manhattan plot using manhattanly::manhattanly -------------------------
    output$geneset_manhattan <- plotly::renderPlotly({
      dat <- geneset_data_main()
      shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["geneset_analysis_main"]]]]), "rvbResult"),
                           "You need an rvbResult object to use this tab!"))
      shiny::validate(need(all(c("CHROM", "POS") %in% rvc@metadata[[values[["geneset_analysis_main"]]]]$columns) ||
                             (!is.null(CHROM) && !is.null(POS)),
                           "'CHROM' and 'POS' columns are required for manhattan plot.
If these columns are present, but have other column names, this can be specified
using the `CHROM` and `POS` parameters.
                           "))
      shiny::validate(need(!(identical(rvc@geneSetList, geneSetList())),
                           "You need a geneSetList object to use this tab!"))
      shiny::validate(need("gsaResult" %in% analyses_types, "You need to load a gsaResult object to use this tab!"))
      shiny::validate(need(nrow(dat) > 0, "No data available for this combination of parameters"))

      if(is.null(CHROM)) CHROM <- "CHROM"
      if(is.null(POS)) POS <- "POS"

      #Make rvbResult information
      dat$effect <- round(dat$effect, 2)
      dat[["CHR"]] <- as.numeric(dat[[CHROM]])
      dat[["BP"]] <- as.numeric(dat[[POS]])
      dat <- dat[,c("unit", "CHR", "BP", "P", values[["assoc_hover_fields"]]),drop=FALSE]
      dat <- dat[complete.cases(dat[,c("CHR", "BP")]),]
      col <- paste0(dat[[values[["assoc_hover_fields"]][1]]], "\n")

      if(length(values[["assoc_hover_fields"]]) > 1) {
        for(i in 2:(length(values[["assoc_hover_fields"]]))){
          col <- paste0(col, values[["assoc_hover_fields"]][i], ": ", dat[[values[["assoc_hover_fields"]][i]]], "\n")
        }
      }
      dat[[values[["assoc_hover_fields"]][[1]]]] <- col
      dat <- manhattanr(dat, gene = "unit",
                        annotation1 = values[["assoc_hover_fields"]][[1]])
      if(length(values[["geneset_custom_threshold"]]) > 0 && input[["geneset_use_custom_threshold"]]) {
        threshold <- values[["geneset_custom_threshold"]]
      } else {
        threshold <-  0.05/nrow(dat$data)
      }
      
      #Make the dataframe for labels of significant variants
      dat_sig = dat$data[dat$data[["P"]] < threshold,,drop=FALSE]

      #make a string with the P-values of the selected geneset
      geneset <- unlist(strsplit(values[["selected_geneset"]], ","))
      shiny::validate(need(length(geneset) == 1, "You need to select 1 gene set!"))

      dat_geneset <- as.data.frame(rvc@rvatResults[sapply(rvc@rvatResults, class) == "gsaResult"][[1]])
      
      shiny::validate(need(geneset %in% dat_geneset$geneSetName, "The selected geneset is not present in the loaded gsaResult"))
      dat_geneset <- dat_geneset[dat_geneset$geneSetName %in% geneset,"P",drop=FALSE]
      dat_geneset <- paste(round(-log10(dat_geneset$P),5), collapse = ",")
      
      #Make dataframe of units that must be highlighted
      units <- units(getGeneSet(rvc@geneSetList, geneset))
      dat_units <- dat$data[dat$data[["unit"]] %in% units,,drop=FALSE]
      
      #Make plot
      plot <- manhattanly(dat,
                          suggestiveline = FALSE,
                          genomewideline = -log10(threshold),
                          title = "") %>%
        plotly::add_trace(data = dat_units,
                          x = dat_units$pos,
                          y = dat_units$logp,
                          xref = "x", yref = "y",
                          color = I("red")) %>%
        plotly::layout(title = paste("-log10(P) values gene set:", dat_geneset))
      
      if(nrow(dat_sig) > 0) {
        plot <- plot %>%
          plotly::add_annotations(x = dat_sig$pos,
                                  y = dat_sig$logp,
                                  text = dat_sig[[input[["assoc_display_variable"]]]],
                                  xref = "x",
                                  yref = "y",
                                  showarrow = TRUE,
                                  arrowhead=2,
                                  ax = 10,
                                  ay = -10) 
      }

      plotly::toWebGL(plot)
    })

    #### Table geneset manhattan -----------------------------------------------
    #geneset table
    output$geneset_info_table_manh <- renderTable({
      geneset <- unlist(strsplit(values[["selected_geneset"]], ","))
      shiny::validate(need(length(geneset) == 1, "You need to select 1 gene set!"))
      
      dat <- as.data.frame(rvc@rvatResults[sapply(rvc@rvatResults, class) == "gsaResult"][[1]])
      shiny::validate(need(geneset %in% dat$geneSetName, "The selected geneset is not present in the loaded gsaResult"))
      dat <- dat[dat$geneSetName %in% geneset,,drop=FALSE]
      dat
    })
    
    #genes table
    
    ### Shown columns change after user input
    observeEvent(input[["geneset_set_variables_genes_manh"]],
                 {
                   values[["geneset_variables_table"]] <- input[["geneset_set_variables_genes_manh"]]
                 }
    )

    ### Checkboxes so the user can select which columns are shown
    output$geneset_edit_variables_genes_manh_ui <- renderUI({
      if(input$geneset_edit_variables_genes_manh) {
        checkboxGroupInput("geneset_set_variables_genes_manh",
                           label = "Fields",
                           choices = metadata(rvc)[[values[["geneset_analysis_main"]]]][["columns"]],
                           selected = c("unit", "P", "effect", "effectSE", "effectCIlower", "effectCIupper"),
                           inline = TRUE)
      }
    })
    
    #### DataTable that the user can see units in the selected gene set in -----
    output$geneset_genes_table_manh <- DT::renderDataTable(
      {
        geneset <- unlist(strsplit(values[["selected_geneset"]],","))
        shiny::validate(need(length(geneset) == 1, "You need to select 1 gene set to plot!"))
        
        dat <- geneset_data_main()  #in this case an rvbResult
        shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["geneset_analysis_main"]]]]), "rvbResult"),
                             "You need an rvbResult object to use this tab!"))
        
        genesetlist <- rvc@geneSetList
        shiny::validate(need(geneset %in% genesetlist@geneSetNames, "The selected geneset is not present in the loaded genesetlist"))
        
        genes_selected <- units(getGeneSet(genesetlist, set = geneset))
        
        if ("geneSetName" %in% values[["geneset_variables_table"]]) {
          values[["geneset_variables_table"]] <- values$geneset_variables_table[!(values$geneset_variables_table %in% "geneSetName")]
        }
        
        if (!("unit" %in% values[["geneset_variables_table"]])) {
          dat <- dat[dat$unit %in% genes_selected,c("unit",values[["geneset_variables_table"]]),drop=FALSE]
        } else {
          dat <- dat[dat$unit %in% genes_selected,values[["geneset_variables_table"]],drop=FALSE]
        }
        
        for (i in 1:nrow(dat)) {
          if (!is.null(UNIT)) {
            if ((!is.null(varSet) & dat[i,"unit"] %in% rvc@metadata$varSet$units) | 
                ("unit" %in% rvc@metadata$varInfo & dat[i,"unit"] %in% rvc@varInfo$unit) | 
                (!is.null(UNIT) & dat[i,"unit"] %in% rvc@varInfo[[UNIT]])) {
              dat[i,"unit"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_genes_manh_geneset&quot;,this.id)">', dat[i,"unit"], '</a>')
            }
          } else {
            if ((!is.null(varSet) & dat[i,"unit"] %in% rvc@metadata$varSet$units) | 
                ("unit" %in% rvc@metadata$varInfo & dat[i,"unit"] %in% rvc@varInfo$unit)) {
              dat[i,"unit"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_genes_manh_geneset&quot;,this.id)">', dat[i,"unit"], '</a>')
            }
          }
        }
        
        number_links <- nrow(dat[substr(dat[["unit"]],1,5) == "<a id",, drop= FALSE])
        if (number_links == 0) {
          showNotification("None of the units could be linked to its variants, so no links to the `Unit viewer` have been made")
        }

        table <- DT::datatable(dat, filter = "none",escape = FALSE, selection = "none")
      })

    ### Make buttons in the table to connect to other tabs
    observeEvent(input$select_button_genes_manh_geneset, {
      geneset <- unlist(strsplit(values[["selected_geneset"]], ","))
      shiny::validate(need(length(geneset) == 1, "You need to select 1 gene set!"))
      
      dat <- geneset_data_main()
      genesetlist <- rvc@geneSetList
      genes_selected <- units(getGeneSet(genesetlist, set = geneset))
      
      if (!("unit" %in% values[["geneset_variables_table"]])) {
        dat <- dat[dat$unit %in% genes_selected,c("unit",values[["assoc_variables_table"]]),drop=FALSE]
      } else {
        dat <- dat[dat$unit %in% genes_selected,values[["assoc_variables_table"]],drop=FALSE]
      }
      
      selectedRow <- as.numeric(strsplit(input$select_button_genes_manh_geneset, "_")[[1]][2])
      values$selected_unit <- dat[selectedRow, "unit"]
      updateNavbarPage(session,
                       "main_menu",
                       "unit_viewer")
    })
    
    ### QQplot -----------------------------------------------------------------
    #### qqplot using manhattanly::qqly ----------------------------------------
    output$geneset_qqplot <- plotly::renderPlotly({
      shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["geneset_analysis_main"]]]]), "gsaResult"),
                           "You need a gsaResult object to use this tab!"))
      shiny::validate(need(nrow(geneset_data_main()) > 0, "No data available for this combination of parameters"))
      
      dat <- geneset_data_main()
      if(length(values[["geneset_custom_threshold"]]) > 0 && input[["geneset_use_custom_threshold"]]) {
        threshold <- values[["geneset_custom_threshold"]]
      } else {
        threshold <- 0.05/nrow(dat)
      }
      
      set.seed(1)
      dat_sig <- dplyr::select(dat, P, input[["geneset_display_variable"]])
      dat_sig$logp <- -log10(dat_sig$P)
      dat_sig <- dat_sig %>% dplyr::arrange(logp)
      dat_sig$null <- sort(-log10(runif(nrow(dat))))
      dat_sig = dat_sig[dat_sig[["logp"]] > -log10(threshold),,drop=FALSE]
      
      if (nrow(dat_sig) > 0) {
        plot <- qqly(dat, threshold = threshold, xlab = "-log10(NULL)", ylab = "-log10(Obs)", title = "") %>%
          plotly::add_annotations(x = dat_sig$null,
                                  y = dat_sig$logp,
                                  text = dat_sig[[input[["geneset_display_variable"]]]],
                                  xref = "x",
                                  yref = "y",
                                  showarrow = TRUE,
                                  arrowhead=2,
                                  ax = 10,
                                  ay = -10)
      } else {
        plot <- qqly(dat, threshold = threshold, xlab = "-log10(NULL)", ylab = "-log10(Obs)", title = "")
      }
      
      plotly::toWebGL(plot)
    })
    
    #### Toptable qqplot - user-determined columns -----------------------------
    observeEvent(input[["geneset_set_variables_toptable_qq"]],
                 {
                   values[["geneset_variables_table"]] <- input[["geneset_set_variables_toptable_qq"]]
                 }
    )
    
    ### Checkboxes so user can select which columns are shown
    output$geneset_edit_variables_toptable_qq_ui <- renderUI({
      if(input$geneset_edit_variables_toptable_qq) {
        checkboxGroupInput("geneset_set_variables_toptable_qq",
                           label = "Fields",
                           choices = metadata(rvc)[[values[["geneset_analysis_main"]]]][["columns"]],
                           selected = c("unit", "P", "effect", "effectSE", "effectCIlower", "effectCIupper"),
                           inline = TRUE)
      }
    })
    
    ### DataTable that user can see the significant values in, refreshes after custom threshold
    output$geneset_toptable_qq <- DT::renderDataTable(
      {
        shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["geneset_analysis_main"]]]]), "gsaResult"),
                             "You need a gsaResult object to use this tab!"))
        dat <- geneset_data_main()   #gsaResult
        
        if(length(values[["geneset_custom_threshold"]]) > 0 && input[["geneset_use_custom_threshold"]]) {
          threshold <- values[["geneset_custom_threshold"]]
        } else {
          threshold <-  0.05/nrow(dat)
        }
        
        if (!("geneSetName" %in% values[["geneset_variables_table"]])) {
          dat <- dat[dat$P < threshold,c("geneSetName",values[["geneset_variables_table"]]),drop=FALSE] %>% dplyr::arrange(P)
        } else {
          dat <- dat[dat$P < threshold,values[["geneset_variables_table"]],drop=FALSE] %>% dplyr::arrange(P)
        }
        
        if (nrow(dat) > 0){
          for (i in 1:nrow(dat)) {
            if (!is.null(geneSetList) & dat[i,"geneSetName"] %in% rvc@metadata$geneSetList) {
              dat[i,"geneSetName"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_qq_geneset&quot;,this.id)">', dat[i,"geneSetName"], '</a>')
            }
          }
        }
        
        number_links <- nrow(dat[substr(dat[["geneSetName"]],1,5) == "<a id",, drop= FALSE])
        if (number_links == 0) {
          showNotification("None of the genesets could be linked to its genes, so no links to `Gene set analysis - manhattan plot` have been made")
        }
      
        table <- DT::datatable(dat, filter = "none",escape = FALSE, selection = "none")
      })
    
    ### Make buttons in the table to connect to other tabs
    observeEvent(input$select_button_qq_geneset, {
      dat <- geneset_data_main()
      
      if(length(values[["geneset_custom_threshold"]]) > 0 && input[["geneset_use_custom_threshold"]]) {
        threshold <- values[["geneset_custom_threshold"]]
      } else {
        threshold <-  0.05/nrow(dat)
      }
      
      if (!("geneSetName" %in% values[["geneset_variables_table"]])) {
        dat <- dat[dat$P < threshold,c("geneSetName",values[["geneset_variables_table"]]),drop=FALSE] %>% dplyr::arrange(P)
      } else {
        dat <- dat[dat$P < threshold,values[["geneset_variables_table"]],drop=FALSE] %>% dplyr::arrange(P)
      }
      
      selectedRow <- as.numeric(strsplit(input$select_button_qq_geneset, "_")[[1]][2])
      values$selected_geneset <- dat[selectedRow, "geneSetName"]
      updateTabsetPanel(session, 
                       "geneset_tabs",
                       "geneset_manh")
    })
    
    ### Density plot -----------------------------------------------------------
    #### Density plot ----------------------------------------------------------
    output$geneset_density_plot <- plotly::renderPlotly({
      shiny::validate(need(!(identical(rvc@geneSetList, geneSetList())), 
                           "You need a geneSetList object to use this tab!"))
      shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["geneset_analysis_main"]]]]), "rvbResult"),
                           "You need an rvbResult object to use this tab!"))
      
      geneset <- unlist(strsplit(input[["geneset_selected_units_density"]],","))
      shiny::validate(need(length(geneset) == 1, "You need to select 1 gene set to plot!"))
      
      dat <- rvc@rvatResults[[values[["geneset_analysis_main"]]]] #Not main, but 'selected' (similar to _compare)
      
      genesetlist <- rvc@geneSetList
      
      geneset_keep <- names(genesetlist)[geneset %in% names(genesetlist)]
      geneset_drop <- names(genesetlist)[!(geneset %in% names(genesetlist))]
      if (length(geneset_drop > 0)) {
        showNotification(paste0("The gene set ", paste(geneset_drop, collapse = ","), " is not in the geneSetList and has been dropped"))
      }
      
      plot <- densityplot(dat, geneset, genesetlist)
      plot <- plotly::ggplotly(plot)
      
      plotly::toWebGL(plot)
    })
    
    #### Density toptable ------------------------------------------------------
    output$geneset_density_table <- DT::renderDataTable({
      shiny::validate(need(!(identical(rvc@geneSetList, geneSetList())), 
                           "You need a geneSetList object to use this tab!"))
      shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["geneset_analysis_main"]]]]), "rvbResult"),
                           "You need an rvbResult object to use this tab!"))
      
      geneset <- unlist(strsplit(input[["geneset_selected_units_density"]],","))
      shiny::validate(need(length(geneset) == 1, "You need to select 1 gene set to plot!"))
      
      if ("gsaResult" %in% analyses_types) {
        dat <- as.data.frame(rvc@rvatResults[analyses_types == "gsaResult"][[1]])
        tableData = data.frame(matrix(NA, ncol=1, nrow=nrow(dat)))[-1]
        for (column in colnames(dat)) {
          if (dplyr::n_distinct(data.frame(dat[column])) > 1) {
            tableData = cbind(tableData, dat[column])
          } 
        } 
        if (!("geneSetName" %in% colnames(tableData))) {
          cbind(dat[["geneSetName"]], tableData)
        }
        tableData <- tableData[tableData[["geneSetName"]] %in% geneset,]
      } else {
        tableData <- data.frame(geneSetName = c(input[["geneset_selected_units_density"]]))
        
      }
      
      if (nrow(tableData) > 0){
        for (i in 1:nrow(tableData)) {
          if (!is.null(geneSetList) & tableData[i,"geneSetName"] %in% rvc@metadata$geneSetList) {
            tableData[i,"geneSetName"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_density_geneset&quot;,this.id)">', tableData[i,"geneSetName"], '</a>')
          }
        }
      }
      
      number_links <- nrow(tableData[substr(tableData[["geneSetName"]],1,5) == "<a id",, drop= FALSE])
      if (number_links == 0) {
        showNotification("None of the genesets could be linked to its genes, so no links to `Gene set analysis - manhattan plot` have been made")
      }
      
      tableData <- DT::datatable(tableData, filter = "none", escape = FALSE, selection = "none")
    })
    
    ### Make buttons in the table to connect to other tabs
    observeEvent(input$select_button_density_geneset, {
      geneset <- unlist(strsplit(input[["geneset_selected_units_density"]],","))
      if ("gsaResult" %in% analyses_types) {
        dat <- as.data.frame(rvc@rvatResults[analyses_types == "gsaResult"][[1]])
        tableData = data.frame(matrix(NA, ncol=1, nrow=nrow(dat)))[-1]
        for (column in colnames(dat)) {
          if (dplyr::n_distinct(data.frame(dat[column])) > 1) {
            tableData = cbind(tableData, dat[column])
          }
        } 
        if (!("geneSetName" %in% colnames(tableData))) {
          cbind(dat[["geneSetName"]], tableData)
        }
        tableData <- tableData[tableData[["geneSetName"]] %in% geneset,]
      } else {
        tableData <- data.frame(geneSetName = c(input[["geneset_selected_units_density"]]))
      }
      
      selectedRow <- as.numeric(strsplit(input$select_button_density_geneset, "_")[[1]][2])
      values$selected_geneset <- tableData[selectedRow, "geneSetName"]
      updateTabsetPanel(session, 
                        "geneset_tabs",
                        "geneset_manh")
    })
    
    output$geneset_density_unique_columns <- renderText({
      dat <- geneset_data_main()
      textData <- list()
      for (column in colnames(dat)) {
        if (dplyr::n_distinct(data.frame(dat[column])) == 1) {
          textData[[column]] = dat[1,column]
        }
      } 
      
      text <- ""
      for (column in names(textData)) {
        text <- paste(text, "<b>", column, "</b>: ", textData[[column]], "<br>")
      }
      text
    })
    
    ### Forest plot ------------------------------------------------------------
    #### Forest plot -----------------------------------------------------------
    output$geneset_forest_plot <- renderPlot({
      unit <- unlist(strsplit(input[["geneset_selected_units_forest"]],","))
      
      shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["geneset_analysis_main"]]]]), "gsaResult"),
                           "You need a gsaResult object to use this tab!"))
      shiny::validate(need(length(unit) > 0, "You need to select gene sets to plot first!"))
      
      dat <- geneset_data_main()
      #unit <- unlist(strsplit(input[["geneset_selected_units_forest"]],","))
      unit_keep <- unit[unit %in% dat[[input[["geneset_display_variable"]]]] ]
      unit_drop <- unit[!(unit %in% dat[[input[["geneset_display_variable"]]]]) ]
      
      if (length(unit_drop) > 0) {
        showNotification(paste0("The units ", paste(unit_drop, collapse = ","), " are not in the dataset and have been dropped. Check if you selected the correct display variable!"))
        #paste0("The units ", paste(unit_drop, collapse = ","), "are not in the dataset and have been dropped. Check if you selected the correct display variable!")
      }
      
      forestData <- data.frame(unit = dat[[input[["geneset_display_variable"]]]],
                               effectSE = dat$effectSE,
                               effect = dat$effect,
                               effectCIlower = dat$effectCIlower,
                               effectCIupper = dat$effectCIupper)
      forestData <- forestData[forestData$unit %in% unit_keep,]
      forestData <- forestData[!is.na(forestData$effectSE),]
      forestData <- rbind(c(rep(NA,ncol(forestData))), forestData)
      
      labelData <- data.frame(unit = dat[[input[["geneset_display_variable"]]]], P = round(dat$P,5), effectSE = round(dat$effectSE,5))
      labelData <- labelData[labelData$unit %in% unit_keep,]
      labelData <- labelData[!is.na(labelData$effectSE),]
      
      labelData <- rbind(colnames(labelData), labelData)
      
      clip = c(floor(min(forestData$effect, na.rm = TRUE)), ceiling(max(forestData$effect, na.rm = TRUE)))
      
      forestplot::forestplot(forestData,
                             mean = effect,
                             lower = effectCIlower,
                             upper = effectCIupper,
                             clip = clip,
                             boxsize = 0.15,
                             xlab = "Effect size",
                             vertices = TRUE,
                             labeltext = labelData,
                             is.summary = c(TRUE, rep(FALSE, nrow(forestData))),
                             hrzl_lines = grid::gpar(lty = 2, col = "#444444"),
                             colgap = grid::unit(5,"mm"),
                             txt_gp = forestplot::fpTxtGp(label = grid::gpar(cex = 1),
                                                          ticks = grid::gpar(cex = 0.5)),
                             zero = 1,
                             cex = 2,
                             lineheight = "auto")
    })
    
    #### Forest plot table -----------------------------------------------------
    forest_table_data_geneset <- reactive({
      shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["geneset_analysis_main"]]]]), "gsaResult"),
                           "You need a gsaResult object to use this tab!"))
      
      unit <- unlist(strsplit(input[["geneset_selected_units_forest"]],","))
      shiny::validate(need(length(unit) >= 1, "You need to select at least 1 gene set to plot!"))
      
      dat <- geneset_data_main()
      unit_keep <- unit[unit %in% dat[[input[["geneset_display_variable"]]]] ]
      
      tableData = data.frame(matrix(NA, ncol=1, nrow=nrow(dat)))[-1]
      for (column in colnames(dat)) {
        if (dplyr::n_distinct(data.frame(dat[column])) > 1) {
          tableData = cbind(tableData, dat[column])
        }
      } 
      
      if (!(input[["geneset_display_variable"]] %in% colnames(tableData))) {
        cbind(dat[[input[["geneset_display_variable"]]]], tableData)
      }
      
      if (!("geneSetName" %in% colnames(tableData))) {
        cbind(dat[["geneSetName"]], tableData)
      }
      
      tableData <- tableData[tableData[[input[["geneset_display_variable"]]]] %in% unit_keep,]
      tableData <- tableData[!is.na(tableData$effectSE),]
    })
    
    output$geneset_forest_table <- DT::renderDataTable({
      unit <- unlist(strsplit(input[["geneset_selected_units_forest"]],","))
      shiny::validate(need(length(unit) > 0, "You need to select at least 1 gene set to plot!"))
      
      tableData <- forest_table_data_geneset()
      
      if (nrow(tableData) > 0){
        for (i in 1:nrow(tableData)) {
          if (!is.null(geneSetList) & tableData[i,"geneSetName"] %in% rvc@metadata$geneSetList) {
            tableData[i,"geneSetName"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_forest_geneset&quot;,this.id)">', tableData[i,"geneSetName"], '</a>')
          }
        }
      }
      
      number_links <- nrow(tableData[substr(tableData[["geneSetName"]],1,5) == "<a id",, drop= FALSE])
      if (number_links == 0) {
        showNotification("None of the genesets could be linked to its genes, so no links to `Gene set analysis - manhattan plot` have been made")
      }
      
      table <- DT::datatable(tableData, filter = "none", rownames = FALSE,
                             escape = FALSE, selection = "none")
    })

    ### Make buttons in the table to connect to other tabs
    observeEvent(input$select_button_forest_geneset, {
      dat <- forest_table_data_geneset()
      selectedRow <- as.numeric(strsplit(input$select_button_forest_geneset, "_")[[1]][2])
      values$selected_geneset <- dat[selectedRow, "geneSetName"]
      updateTabsetPanel(session, 
                        "geneset_tabs",
                        "geneset_manh")
    })
    
    output$geneset_forest_unique_columns <- renderText({
      unit <- unlist(strsplit(input[["geneset_selected_units_forest"]],","))
      shiny::validate(need(length(unit) > 0, "You need to select at least 1 gene set to plot!"))
      
      dat <- geneset_data_main()
      textData <- list()
      for (column in colnames(dat)) {
        if (dplyr::n_distinct(data.frame(dat[column])) == 1) {
          textData[[column]] = dat[1,column]
        }
      } 
      
      text <- ""
      for (column in names(textData)) {
        text <- paste(text, "<b>", column, "</b>: ", textData[[column]], "<br>")
      }
      text
    })
    
    ### tableViewer ------------------------------------------------------------
    observeEvent(input[["geneset_set_variables_tableviewer"]],
                 {
                   values[["geneset_variables_table"]] <- input[["geneset_set_variables_tableviewer"]]
                 }
    )
    
    ### Checkboxes with variables the user can choose from
    output$geneset_edit_variables_tableviewer_ui <- renderUI({
      if(input$geneset_edit_variables_tableviewer) {
        checkboxGroupInput("geneset_set_variables_tableviewer",
                           label = "Fields",
                           choices = metadata(rvc)[[values[["geneset_analysis_main"]]]][["columns"]],
                           selected = c("geneSetName", "P", "effect", "effectSE", "effectCIlower", "effectCIupper"),
                           inline = TRUE)
        
      }
    })
    
    ### DataTable that the user can filter in
    output$geneset_tableviewer <- DT::renderDT(
      {
        shiny::validate(need(identical(checkClassrvatResult(rvc@rvatResults[[values[["geneset_analysis_main"]]]]), "gsaResult"),
                             "You need a gsaResult object to use this tab!"))
        dat <- as.data.frame(rvc@rvatResults[[values[["geneset_analysis_main"]]]])
        
        if (nrow(dat) > 0){
          for (i in 1:nrow(dat)) {
            if (!is.null(geneSetList) & dat[i,"geneSetName"] %in% rvc@metadata$geneSetList) {
              dat[i,"geneSetName"] <- paste0('<a id="button_', i, '" href="#" class="action-button" onclick="Shiny.onInputChange(&quot;select_button_tableviewer_geneset&quot;,this.id)">', dat[i,"geneSetName"], '</a>')
            }
          }
        }
        
        number_links <- nrow(dat[substr(dat[["geneSetName"]],1,5) == "<a id",, drop= FALSE])
        if (number_links == 0) {
          showNotification("None of the genesets could be linked to its genes, so no links to `Gene set analysis - manhattan plot` have been made")
        }
        
        if (!("geneSetName" %in% values[["geneset_variables_table"]])) {
          DT::datatable(
            dat[,c("geneSetName",values[["geneset_variables_table"]]),drop=FALSE], 
            filter = "top",
            rownames = input[["geneset_show_index_tableviewer"]],
            escape = FALSE, selection = "none")
        } else {
          DT::datatable(
            dat[,values[["geneset_variables_table"]],drop=FALSE], 
            filter = "top",
            rownames = input[["geneset_show_index_tableviewer"]],
            escape = FALSE, selection = "none")
        }
      })

    ### Make buttons in the table to connect to other tabs
    observeEvent(input$select_button_tableviewer_geneset, {
    dat <- as.data.frame(rvc@rvatResults[[values[["geneset_analysis_main"]]]])
    selectedRow <- as.numeric(strsplit(input$select_button_tableviewer_geneset, "_")[[1]][2])
    values$selected_geneset <- dat[selectedRow, "geneSetName"]
    updateTabsetPanel(session, 
                      "geneset_tabs",
                      "geneset_manh")
    })
    
  }

  shinyApp(ui, server)
}
