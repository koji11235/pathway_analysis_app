

shinyUI(
    fluidPage(
        # Application title
        titlePanel("a"),
        
        # Sidebar with a slider input for number of bins
        sidebarLayout(
            sidebarPanel(
            ),
            
            # Show a plot of the generated distribution
            mainPanel(
                DT::dataTableOutput("kegg_table"),
                imageOutput("KeggImage"),
                downloadButton('downloadKEGGImage', 'Download kegg pathway image')
            )
        )
    )
)