library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Pathway App"),

    sidebarLayout(
        sidebarPanel(
            
            fileInput("file", "Upload CSV file",
                      accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv")
            ),
            sliderInput("show_category",
                        "Show Category",
                        min = 1,
                        max = 20,
                        value = 8),
        ),
        

        # Show a plot of the generated distribution
        mainPanel(
            
            plotlyOutput("barPlot"),
            tabsetPanel(type = "tabs",
                        tabPanel("Table", tableOutput("table"))
            ),
        )
    )
))





# 
# ui <- dashboardPage(
#     dashboardHeader(title = "Pathway Analysis App"),
#     
#     dashboardSidebar(
#         sidebarMenu(
#             menuItem("Data Upload", tabName = "data_upload"),
#             menuItem("KEGG", tabName = "KEGG"),
#             menuItem("Reactom", tabName = "Reactom")
#         )
#     ),
#     dashboardBody(
#         tabItems(
#             # Data Upload Panel --------------------------------------------------
#             tabItem(
#                 tabName = "data_upload",
#                 titlePanel("Data Upload"),
#                 sidebarLayout(
#                     sidebarPanel(
#                         fileInput(
#                             "file", "Upload CSV file",
#                             accept = c(
#                                 "text/csv",
#                                 "text/comma-separated-values,text/plain",
#                                 ".csv")
#                         ),
#                         textInput(
#                             "GeneID_column", "GeneID column",
#                             value = "Enter GeneID column name"
#                         )
#                     ),
#                     mainPanel(
#                         tabsetPanel(
#                             type = "tabs",
#                             tabPanel(
#                                 "Table", 
#                                 tableOutput("table"))
#                         ),
#                     )
#                 )
#             ),
#             # KEGG Result Panel --------------------------------------------------
#             tabItem(
#                 tabName = "KEGG",
#                 titlePanel("KEGG barplot"),
#                 fluidRow(
#                     box(
#                         title = "barplot",
#                         plotlyOutput("KEGGbarPlot")
#                     ),
#                     box(
#                         title = "barplot",
#                     )
#                 ),
#                 fluidRow(
#                     box(
#                         title = "barplot",
#                         plotlyOutput("KEGGbarPlot")
#                     ),
#                     box(
#                         title = "barplot",
#                     )
#                 )
#             ),
#             
#             # Reactome Result Panel --------------------------------------------------
#             tabItem(tabName = "Reactom",
#                     titlePanel("Reactom barplot"))
#         )
#     )
# )
# 
