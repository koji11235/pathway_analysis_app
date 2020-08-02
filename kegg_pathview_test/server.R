
enrichment_result <- read.csv("sample_enrichment_resul.csv") %>% head(5)

shinyServer(function(input, output) {

    # enrichment_result <- reactive({
    #     return (enrichment_result)
    # })
    
    output$kegg_table <- DT::renderDataTable(
        datatable(
            {enrichment_result %>% 
                mutate(pvalue = pvalue %>% signif(3),
                       p.adjust = p.adjust %>% signif(3),
                       qvalue = qvalue %>% signif(3))},
            extensions="Buttons",
            options = list(
                buttons = c('csv', 'excel')
            ),
        ),
        selection = 'multiple'
    )
    
    output$KeggImage <- renderImage({
        temp_dir <- paste0(getwd(),"/temp")
        if (!dir.exists(temp_dir)){
            dir.create(temp_dir)
        }
        
        #rows_selected <- input$kegg_table_rows_selected
        row_selected <- input$kegg_table_row_last_clicked
        pathID <- enrichment_result$ID[[row_selected]]
        genes <- enrichment_result$geneID[[row_selected]] %>% str_split("/") %>% .[[1]]
        
        pv.out <- mypathview(gene.data = genes, pathway.id = pathID, kegg.dir = temp_dir)
        
        unnecessary.png.file <- paste0(temp_dir,"/",pathID,".png")
        unnecessary.xml.file <- paste0(temp_dir,"/",pathID,".xml")
        file.remove(unnecessary.png.file)
        file.remove(unnecessary.xml.file)
        
        #outfile <- paste( tempFolder,"/",pathID,".",randomString,".png",sep="")
        outfile <- paste(temp_dir,"/",pathID,".","pathview",".png",sep="")
        # outfile_tempdir <- paste(temp_dir,"/",pathID,".","pathview",".png",sep="")
        # file.rename(outfile, outfile_tempdir)
        # file.remove(outfile)
        
        img <- readPNG(outfile)
        width <- ncol(img)* 1/2
        height <- nrow(img)* 1/2
        
        list(src = outfile_tempdir,
             contentType = 'image/png',
             width = width,
             height = height,
             alt = "KEGG pathway image.")
    }, deleteFile = FALSE)
    
    output$downloadKEGGImage <- downloadHandler(
        filename = function(){
            row_selected <- input$kegg_table_row_last_clicked
            pathID <- enrichment_result$ID[[row_selected]]
            paste(pathID,".png",sep="")
        },
        content = function(file) {
            temp_dir <- paste0(getwd(),"/temp")
            row_selected <- input$kegg_table_row_last_clicked
            pathID <- enrichment_result$ID[[row_selected]]
            content_file <- paste(temp_dir,"/",pathID,".","pathview",".png",sep="")
            #print(content_file)
            file.copy(content_file,file, overwrite = TRUE)
        },
        contentType = 'image/png'
    )
    

})
