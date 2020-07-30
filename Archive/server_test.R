#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(clusterProfiler)
library(ReactomePA)
library(clusterProfiler)

library(tidyverse)
library(DOSE)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  observeEvent(input$file, {
    csv_file <- reactive(read_csv(input$file$datapath))
    output$table <- renderTable(csv_file())
  })
  
  output$KEGGbarPlot <- renderPlotly({
    
    sample_data <- read.csv(input$file$datapath)
    gene_IDs<-sample_data$Entrez_Gene_ID
    #print(gene_IDs)
    kegg_enrichment_result <- enrichKEGG(gene=gene_IDs,pvalueCutoff=0.05)
    print("dane")
    
    kegg_enrichment_bp　<-ggplotly( 
      barplot(kegg_enrichment_result, drop=TRUE, showCategory=8)#input$show_category)
    )
    
  })
  
  
  
})



