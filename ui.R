


ui <- dashboardPage(
  dashboardHeader(title = "Pathway Analysis App"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload", tabName = "data_upload"),
      menuItem("KEGG", tabName = "KEGG"),
      menuItem("Reactome", tabName = "Reactome"),
      menuItem("GO BP", tabName = "GO_BP"),
      menuItem("GO MF", tabName = "GO_MF"),
      menuItem("GO CC", tabName = "GO_CC")
    )
  ),
  dashboardBody(
    tabItems(
      # Data Upload Panel --------------------------------------------------
      tabItem(
        tabName = "data_upload",
        titlePanel("Data Upload"),
        sidebarLayout(
          sidebarPanel(
            fileInput(
              "file", "Upload CSV file",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")
            ),
            textInput(
              "GeneID_column", "GeneID column",
              value = "Enter GeneID column name"
            )
          ),
          mainPanel(
            tabsetPanel(
              type = "tabs",
              tabPanel(
                "Table", 
                DT::dataTableOutput("table"))
            ),
          )
        )
      ),
      # KEGG Result Panel --------------------------------------------------
      tabItem(
        tabName = "KEGG",
        titlePanel("KEGG"),
        fluidRow(
          box(
            title = "bar plot",
            plotlyOutput("KEGG_barPlot"),
            sliderInput("KEGG_barPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12, step = 1,
                        value = 8)
          ),
          box(
            title = "dot plot",
            plotlyOutput("KEGG_dotPlot"),
            sliderInput("KEGG_dotPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 8)
          )
        ),
        fluidRow(
          box(
            title = "emap plot",
            plotlyOutput("KEGG_emapPlot"),
            sliderInput("KEGG_emapPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 8)
          ),
          box(
            title = "cnet plot",
            plotlyOutput("KEGG_cnetPlot"),
            sliderInput("KEGG_cnetPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 3)
          )
        ),
        fluidRow(
          tabPanel(
            "Table", 
            DT::dataTableOutput("kegg_table")
          )
        )
      ), # tabItem
      
      # Reactome Result Panel --------------------------------------------------
      tabItem(
        tabName = "Reactome",
        titlePanel("Reactome"),
        fluidRow(
          box(
            title = "bar plot",
            plotlyOutput("Reactome_barPlot"),
            sliderInput("Reactome_barPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12, step = 1,
                        value = 8)
          ),
          box(
            title = "dot plot",
            plotlyOutput("Reactome_dotPlot"),
            sliderInput("Reactome_dotPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 8)
          )
        ),
        fluidRow(
          box(
            title = "emap plot",
            plotlyOutput("Reactome_emapPlot"),
            sliderInput("Reactome_emapPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 8)
          ),
          box(
            title = "cnet plot",
            plotlyOutput("Reactome_cnetPlot"),
            sliderInput("Reactome_cnetPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 3)
          )
        ),
        fluidRow(
          tabPanel(
            "Table", 
            DT::dataTableOutput("reactome_table")
          )
        )
      ), # tabItem
      # GO_BP Result Panel --------------------------------------------------
      tabItem(
        tabName = "GO_BP",
        titlePanel("GO BP"),
        fluidRow(
          box(
            title = "bar plot",
            plotlyOutput("GO_BP_barPlot"),
            sliderInput("GO_BP_barPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12, step = 1,
                        value = 8)
          ),
          box(
            title = "dot plot",
            plotlyOutput("GO_BP_dotPlot"),
            sliderInput("GO_BP_dotPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 8)
          )
        ),
        fluidRow(
          box(
            title = "emap plot",
            plotlyOutput("GO_BP_emapPlot"),
            sliderInput("GO_BP_emapPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 8)
          ),
          box(
            title = "cnet plot",
            plotlyOutput("GO_BP_cnetPlot"),
            sliderInput("GO_BP_cnetPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 3)
          )
        ),
        fluidRow(
          tabPanel(
            "Table", 
            DT::dataTableOutput("GO_BP_table")
          )
        )
      ), # tabItem
      # GO_MF Result Panel --------------------------------------------------
      tabItem(
        tabName = "GO_MF",
        titlePanel("GO MF"),
        fluidRow(
          box(
            title = "bar plot",
            plotlyOutput("GO_MF_barPlot"),
            sliderInput("GO_MF_barPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12, step = 1,
                        value = 8)
          ),
          box(
            title = "dot plot",
            plotlyOutput("GO_MF_dotPlot"),
            sliderInput("GO_MF_dotPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 8)
          )
        ),
        fluidRow(
          box(
            title = "emap plot",
            plotlyOutput("GO_MF_emapPlot"),
            sliderInput("GO_MF_emapPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 8)
          ),
          box(
            title = "cnet plot",
            plotlyOutput("GO_MF_cnetPlot"),
            sliderInput("GO_MF_cnetPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 3)
          )
        ),
        fluidRow(
          tabPanel(
            "Table", 
            DT::dataTableOutput("GO_MF_table")
          )
        )
      ), # tabItem
      # GO_CC Result Panel --------------------------------------------------
      tabItem(
        tabName = "GO_CC",
        titlePanel("GO CC"),
        fluidRow(
          box(
            title = "bar plot",
            plotlyOutput("GO_CC_barPlot"),
            sliderInput("GO_CC_barPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12, step = 1,
                        value = 8)
          ),
          box(
            title = "dot plot",
            plotlyOutput("GO_CC_dotPlot"),
            sliderInput("GO_CC_dotPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 8)
          )
        ),
        fluidRow(
          box(
            title = "emap plot",
            plotlyOutput("GO_CC_emapPlot"),
            sliderInput("GO_CC_emapPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 8)
          ),
          box(
            title = "cnet plot",
            plotlyOutput("GO_CC_cnetPlot"),
            sliderInput("GO_CC_cnetPlot_showCategory", 
                        "Show Categories",
                        width = "50%",
                        min = 1, max = 12,
                        value = 3)
          )
        ),
        fluidRow(
          tabPanel(
            "Table", 
            DT::dataTableOutput("GO_CC_table")
          )
        )
      ) # tabItem
    ) # tabItems
  ) # dashboardBody
) # dashboardPage
  
  
  
  
  
  # Define UI for application that draws a histogram
  # shinyUI(fluidPage(
  #   
  #   # Application title
#   titlePanel("Pathway App"),
#   
#   sidebarLayout(
#     sidebarPanel(
#       
#       fileInput("file", "Upload CSV file",
#                 accept = c(
#                   "text/csv",
#                   "text/comma-separated-values,text/plain",
#                   ".csv")
#       ),
#       sliderInput("show_category",
#                   "Show Category",
#                   min = 1,
#                   max = 20,
#                   value = 8),
#     ),
#     
#     
#     # Show a plot of the generated distribution
#     mainPanel(
#       
#       plotlyOutput("barPlot"),
#       tabsetPanel(type = "tabs",
#                   tabPanel("Table", tableOutput("table"))
#       ),
#     )
#   )
# ))
