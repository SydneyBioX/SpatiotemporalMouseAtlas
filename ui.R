library(shinyjs)
library(shinythemes)

shinyUI(fluidPage(
    useShinyjs(),
    theme=shinytheme("cerulean"),

    tags$head(
      tags$style(HTML("
      .flexibleActionButton {
        flext: 1;
        word-wrap: break-word;
        white-space: normal;
      }
    "))
    ),

    titlePanel("Spatiotemporal Mouse Atlas"),

    fluidRow(

        column(1,

               checkboxGroupInput("embryo_subset",
                                  "Subset by embryo:",
                                  choiceNames = c(paste0("Embryo ", 1:7, c(rep(" (E6.5)", 2), rep(" (E7.5)", 2), rep(" (E8.5)", 3)))),
                                  choiceValues = c(paste0("embryo", 7:1)),
                                  selected = "embryo4"),
               
               checkboxInput("embryo_same_size",
                             "Make embryos same size",
                             value = FALSE),
               
               

               # HTML("z-slices are contiguous cell layers 12 um apart."),

               # commented 2024
               # checkboxGroupInput("z_subset",
               #                    "Subset by z-slice:",
               #                    choiceNames = c("Slice 1", "Slice 2"),
               #                    choiceValues = zvals,
               #                    selected = 2),

               # checkboxInput('show_segmentation',
               #               'Display cell segmentation in Spatial plots',
               #               value = FALSE)

        ),
        column(11,

               tabsetPanel(id = "tabs",
                           selected = "Virtual dissection", ## REMOVE
                           tabPanel("Landing page", fluid = TRUE,
                                    includeMarkdown("README.md")
                           ),

                           tabPanel("Spatial plot", fluid = TRUE,

                                    column(2,

                                           selectInput("colour_by",
                                                       "Colour by:",
                                                       choices = c("Expression logcounts",
                                                                   "Mapped cell type"
                                                       ),
                                                       selected = "Expression logcounts",
                                                       multiple = FALSE),

                                           selectizeInput("gene_name",
                                                          "Gene name:",
                                                          choices = genes,
                                                          selected = "T"),

                                           checkboxInput("celltype_subset_all_SP",
                                                         "Show all cell types",
                                                         value = TRUE),

                                           selectizeInput("celltype_subset",
                                                          "Subset by cell type:",
                                                          choices = celltypes,
                                                          selected = NULL,
                                                          multiple = TRUE),

                                           sliderInput('spatialPlotPointSize',
                                                       "Point size",
                                                       min = 0.05,
                                                       max = 1.5,
                                                       value = 0.3,
                                                       step = 0.05),

                                    ),

                                    # Show a plot of the generated distribution
                                    column(10,

                                               HTML("<p style=\"text-align: center;\"><em>Plot of spatial expression or cell types, and distribution of expression per cell type (below). Check option on left panel for plot of cells' segmentation. Use the icon in the top-right corner to download each image.</em></p>"),

                                               girafeOutput("spatialPlot",
                                                            height = "100%",
                                                            width = "100%"
                                               ) %>% withSpinner(color="#0dc5c1"),

                                               plotOutput("spatialPlot_leg",
                                                          height = "200px"),

                                               girafeOutput("violinPlot",
                                                            height = "100%",
                                                            width = "100%"
                                               ) %>% withSpinner(color="#0dc5c1")
                                    )
                           ),
                           tabPanel("Spatial plot (Imputed)", fluid = TRUE,

                                    column(2,

                                           # textInput("gene_name_imp",
                                           #           "Please type in a gene name:",
                                           #           value = "T"),
                                           
                                           selectizeInput("gene_name_imp",
                                                          "Gene name:",
                                                          choices = sort(genes_imp),
                                                          selected = c("T"),
                                                          options= list(maxOptions = length(genes_imp))),

                                           textOutput("gene_name_imp_parse_status",
                                                      inline = TRUE),

                                           checkboxInput("celltype_subset_all_SPImputed",
                                                         "Show all cell types",
                                                         value = TRUE),

                                           selectizeInput("celltype_subset_imp",
                                                          "Subset by cell type:",
                                                          choices = celltypes,
                                                          selected = NULL,
                                                          multiple = TRUE),
                                           
                                           sliderInput('spatialPlotPointSizeImp',
                                                       "Point size",
                                                       min = 0.05,
                                                       max = 1.5,
                                                       value = 0.3,
                                                       step = 0.05)
                                    ),

                                    column(10,

                                               HTML("<p style=\"text-align: center;\"><em>Failed QC cells are removed here. Plot of imputed spatial expression, and distribution of imputed expression per cell type (below). Check option on left panel for plot of cells' segmentation.</em></p>"),

                                               girafeOutput("spatialPlot_imp",
                                                            height = "100%",
                                                            width = "100%"
                                               ) %>% withSpinner(color="#0dc5c1"),

                                               plotOutput("spatialPlot_leg_imp",
                                                          height = "200px"),

                                               girafeOutput("violinPlot_imp",
                                                            height = "100%",
                                                            width = "100%"
                                               ) %>% withSpinner(color="#0dc5c1")
                                    )
                           ),
                           
                           tabPanel("Virtual dissection", fluid = TRUE,

                                    column(2,

                                           selectizeInput('virtualdissect_choice',
                                                          'Virtually dissect using...',
                                                          choices = c("Physical", "UMAP"),
                                                          selected = "Physical",
                                                          multiple = FALSE),

                                          selectInput("virtual_colour_by",
                                                      "Colour by:",
                                                      choices = c("Expression logcounts",
                                                                  "Imputed expression",
                                                                  "Mapped cell type",
                                                                  "AP axis",
                                                                  "DV axis"
                                                      ),
                                                      selected = "Expression logcounts",
                                                      multiple = FALSE),

                                          selectizeInput("virtual_gene_name",
                                                         "Gene name (Panel):",
                                                         choices = genes,
                                                         selected = "T"),
                                          
                                          selectizeInput("virtual_gene_name_imputed",
                                                         "Gene name (Imputed):",
                                                         choices = sort(genes_imp),
                                                         selected = "T",
                                                         options= list(maxOptions = length(genes_imp))),

                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                          

                                           actionButton("add", "Add selection to Group A",
                                                        width = "100%",
                                                        class = "flexibleActionButton",
                                                        style="color:#ff0000"),


                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),

                                           actionButton("add2", "Add selection to Group B",
                                                        width = "100%",
                                                        class = "flexibleActionButton",
                                                        style="color:#0000ff"),

                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),

                                           sliderInput('dissectionPointSize',
                                                       "Point size",
                                                       min = 0.1,
                                                       max = 1.5,
                                                       value = 0.5,
                                                       step = 0.1),

                                          sliderInput('dissectionStrokeSize',
                                                      "Selected outline width",
                                                      min = 0.1,
                                                      max = 0.5,
                                                      value = 0.3,
                                                      step = 0.05),

                                           HTML("<p style=\"margin-bottom:1cm;\"> </p>"),

                                           selectizeInput('virtualDissectionCelltypes',
                                                          'Mapped cell types',
                                                          choices = celltypes,
                                                          selected = NULL,
                                                          multiple = TRUE),

                                           actionButton("addCT", "Add listed cell types to Group A",
                                                        width = "100%",
                                                        class = "flexibleActionButton",
                                                        style="color:#ff0000"),

                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),

                                           actionButton("addCT2", "Add listed cell types to Group B",
                                                        width = "100%",
                                                        class = "flexibleActionButton",
                                                        style="color:#0000ff"),

                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),

                                           actionButton("removeCT", "Remove listed cell types from either Group",
                                                        width = "100%",
                                                        class = "flexibleActionButton",
                                                        style="font-size:12px"),

                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),


                                          actionButton("removeAllButCT", "Remove all except listed cell types from either Group",
                                                       width = "100%",
                                                       class = "flexibleActionButton",
                                                       style="font-size:12px"),

                                          HTML("<p style=\"margin-bottom:5mm;\"> </p>"),

                                           HTML("Only selected cells within the selected embryos will be downloaded."),

                                           downloadButton("virtualDissectionCellDownload",
                                                          "Download Group A cell names",
                                                          class = "flexibleActionButton",
                                                          style="width:100%"),

                                           downloadButton("virtualDissectionCellDownload2",
                                                          "Download Group B cell names",
                                                          class = "flexibleActionButton",
                                                          style="width:100%"),

                                           HTML("<p style=\"margin-bottom:1cm;\"> </p>"),

                                           HTML("You can upload multiple cell pre-selection files. Once uploaded press \"Add pre-selected\" to either Group."),

                                           fileInput("preselected", "Upload pre-selected cells",
                                                     multiple = TRUE,
                                                     accept = c(".Rds")),

                                           actionButton("add_pre", "Add pre-selected to Group A",
                                                        width = "100%",
                                                        class = "flexibleActionButton",
                                                        style = "color:#ff0000"),

                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),

                                           actionButton("add_pre2", "Add pre-selected to Group B",
                                                        width = "100%",
                                                        class = "flexibleActionButton",
                                                        style = "color:#0000ff"),

                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),

                                           actionButton("remove_pre", "Remove pre-selected cells from either Group",
                                                        width = "100%",
                                                        class = "flexibleActionButton",
                                                        style = "font-size:12px")
                                    ),

                                    # Show a plot of the generated distribution
                                    column(10,
                                               HTML("<p style=\"text-align: center;\"><em>Use the lasso selection icon (above the plot) to select cells, then click each button to assign them to a group. Once you have selected two groups, an MA-plot will appear below.</em></p>"),

                                           fluidRow(
                                             actionButton("removeAll",
                                                          HTML("<i>Reset all cells to Unselected</i>"),
                                                          class = "flexibleActionButton"),
                                             
                                             actionButton("remove",
                                                          "Remove selection from either Group",
                                                          class = "flexibleActionButton"),

                                             actionButton("virtualDissection_reset",
                                                          label = "Clear lasso selection",
                                                          class = "flexibleActionButton"),

                                           ),
                                           
                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),

                                           girafeOutput("virtualDissection",
                                                        height = "100%",
                                                        width = "100%"
                                             ) %>% withSpinner(color="#0dc5c1"),

                                           plotOutput("virtualDissection_leg",
                                                      height = "200px"),

                                           HTML("<p style=\"margin-bottom:5mm;\"> </p>"),
                                           
                                           plotOutput("VirtualDissectionBarPlot",
                                                      height = "700px") |> 
                                             withSpinner(color="#0dc5c1"),
                                           
                                           girafeOutput("VirtualDissectionMAPlot",
                                                        height = "700px",
                                                        width = "100%") %>% 
                                             withSpinner(color="#0dc5c1")
                                    )
                           ),
                           
                           
                           # AP DV visualisation panel
                           tabPanel("AP or DV axis", fluid = TRUE,
                                    

                                    
                                    column(2,
                                           
                                           selectizeInput("ap_or_dv",
                                                          "Plot axis:",
                                                          choices = c("AP", "DV"),
                                                          selected = "AP"),
                                           
                                           checkboxInput("celltype_subset_all_apdv",
                                                         "Show all cell types",
                                                         value = TRUE),
                                           
                                           selectizeInput("apdv_plot_celltype",
                                                       "Subset by cell type:",
                                                       choices = celltypes,
                                                       selected = NULL,
                                                       multiple = TRUE),
                                           
                                           selectizeInput("apdv_plot_genename",
                                                          "Gene name:",
                                                          choices = sort(genes_imp),
                                                          selected = c("T"),
                                                          options= list(maxOptions = length(genes_imp)),
                                                          multiple = TRUE),
                                           
                                           
                                           
                                           ),
                                    
                                    # Show a plot of the generated distribution
                                    column(10,
                                           
                                           HTML("<p style=\"text-align: center;\"><em>Plot of gene expression along AP or DV axis.</em></p>"),
                                           
                                           plotlyOutput(outputId= "APDV_plot", height = '1000px')
                                           
                                    )
                                    
                                    
                                    
                                    ) 
                           )
               )
        )
    )
)