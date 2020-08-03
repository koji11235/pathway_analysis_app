shinyServer(function(input, output, session) {

    observeEvent(input$file, {
        input_data <- reactive(read_csv(input$file$datapath))
        output$table <- DT::renderDataTable(input_data())
    })
    
    input_data <- reactive({
        return (read_csv(input$file$datapath))
    })
    
    gene_IDs <- reactive({
        # columnを変更できるようにする
        print(input$GeneID_column)
        if (input$GeneID_column == ""){
            gene_IDs = input_data()[[1]]
        }else{
            gene_IDs = input_data()[[input$GeneID_column]]
        }
        return (gene_IDs)
    })

    # Plot KEGG results ------------------------------------------------------------------
    kegg_enrichment_result <- reactive({
        kegg_enrichment_result <- enrichKEGG(gene=gene_IDs(), pvalueCutoff=0.05)
        return (kegg_enrichment_result)
    })
    
    output$kegg_table <- DT::renderDataTable(
        as.data.frame(kegg_enrichment_result()) %>% 
            mutate(pvalue = pvalue %>% signif(3),
                   p.adjust = p.adjust %>% signif(3),
                   qvalue = qvalue %>% signif(3))
    )
    
    output$KEGG_barPlot <- renderPlotly({     
        kegg_enrichment_bp <- barplot_plotly(kegg_enrichment_result(), 
                                             showCategory = input$KEGG_barPlot_showCategory,
                                             selected_rows = input$kegg_table_rows_selected)
    })
    
    output$KEGG_dotPlot <- renderPlotly({
        kegg_enrichment_dp <- dotplot_plotly(kegg_enrichment_result(), 
                                             showCategory = input$KEGG_dotPlot_showCategory,
                                             selected_rows = input$kegg_table_rows_selected)
    })
    
    output$KEGG_emapPlot <- renderPlotly({
        kegg_enrichment_emap <- emap_plotly(kegg_enrichment_result(), 
                                            showCategory = input$KEGG_emapPlot_showCategory, 
                                            line_scale=3, color="p.adjust",
                                            selected_rows = input$kegg_table_rows_selected)
    })
    
    output$KEGG_cnetPlot <- renderPlotly({
        kegg_enrichment_cnet <- cnet_plotly(kegg_enrichment_result(),
                                            showCategory = input$KEGG_cnetPlot_showCategory,
                                            selected_rows = input$kegg_table_rows_selected,
                                            foldChange = NULL, 
                                            is_kegg = TRUE) 
    })
    
    output$KeggImage <- renderImage({
        temp_dir <- paste0(getwd(),"/temp")
        if (!dir.exists(temp_dir)){
            dir.create(temp_dir)
        }
        
        #rows_selected <- input$kegg_table_rows_selected
        row_selected <- input$kegg_table_row_last_clicked
        pathID <- kegg_enrichment_result()$ID[[row_selected]]
        genes <- kegg_enrichment_result()$geneID[[row_selected]] %>% str_split("/") %>% .[[1]]
        
        pv.out <- mypathview(gene.data = genes, pathway.id = pathID, kegg.dir = temp_dir)
        
        unnecessary.png.file <- paste0(temp_dir,"/",pathID,".png")
        unnecessary.xml.file <- paste0(temp_dir,"/",pathID,".xml")
        file.remove(unnecessary.png.file)
        file.remove(unnecessary.xml.file)
        
        #outfile <- paste( tempFolder,"/",pathID,".",randomString,".png",sep="")
        outfile <- paste(temp_dir,"/",pathID,".","pathview",".png",sep="")
        
        list(src = outfile,
             contentType = 'image/png',
             width = "100%",
             height = "100%",
             alt = "KEGG pathway image.")
    }, deleteFile = FALSE)
    
    output$downloadKEGGImage <- downloadHandler(
        filename = function(){
            row_selected <- input$kegg_table_row_last_clicked
            pathID <- kegg_enrichment_result()$ID[[row_selected]]
            paste(pathID,".png",sep="")
        },
        content = function(file) {
            temp_dir <- paste0(getwd(),"/temp")
            row_selected <- input$kegg_table_row_last_clicked
            pathID <- kegg_enrichment_result()$ID[[row_selected]]
            content_file <- paste(temp_dir,"/",pathID,".","pathview",".png",sep="")
            #print(content_file)
            file.copy(content_file,file, overwrite = TRUE)
        },
        contentType = 'image/png'
    )
    
    
    
    # Plot Reactome results ------------------------------------------------------------------
    reactome_enrichment_result <- reactive({
        return (enrichPathway(gene=gene_IDs(), pvalueCutoff=0.05, readable = TRUE))
    })
    
    output$reactome_table <- DT::renderDataTable(
        as.data.frame(reactome_enrichment_result()) %>% 
            mutate(pvalue = pvalue %>% signif(3),
                   p.adjust = p.adjust %>% signif(3),
                   qvalue = qvalue %>% signif(3))
    )
    
    
    output$Reactome_barPlot <- renderPlotly({     
        reactome_enrichment_bp <- barplot_plotly(reactome_enrichment_result(), 
                                             showCategory = input$Reactome_barPlot_showCategory,
                                             selected_rows = input$reactome_table_rows_selected)
    })
    
    output$Reactome_dotPlot <- renderPlotly({
        reactome_enrichment_dp <- dotplot_plotly(reactome_enrichment_result(), 
                                             showCategory = input$Reactome_dotPlot_showCategory,
                                             selected_rows = input$reactome_table_rows_selected)
    })
    
    output$Reactome_emapPlot <- renderPlotly({
        reactome_enrichment_emap <- emap_plotly(reactome_enrichment_result(), 
                                            showCategory = input$Reactome_emapPlot_showCategory,
                                            selected_rows = input$reactome_table_rows_selected,
                                            line_scale=1, color="p.adjust")
    })
    
    output$Reactome_cnetPlot <- renderPlotly({
        reactome_enrichment_cnet <- cnet_plotly(reactome_enrichment_result(),
                                            showCategory = input$Reactome_cnetPlot_showCategory,
                                            selected_rows = input$reactome_table_rows_selected,
                                            foldChange = NULL) 
    })
    
    
    # Plot GO_BP results ------------------------------------------------------------------
    GO_BP_enrichment_result <- reactive({
        return (enrichGO(gene=gene_IDs(), pvalueCutoff=0.05, readable = TRUE,
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP"))
    })
    
    output$GO_BP_table <- DT::renderDataTable(
        as.data.frame(GO_BP_enrichment_result()) %>% 
            mutate(pvalue = pvalue %>% signif(3),
                   p.adjust = p.adjust %>% signif(3),
                   qvalue = qvalue %>% signif(3))
    )
    
    
    output$GO_BP_barPlot <- renderPlotly({     
        GO_BP_enrichment_bp <- barplot_plotly(GO_BP_enrichment_result(), 
                                              showCategory = input$GO_BP_barPlot_showCategory,
                                              selected_rows = input$GO_BP_table_rows_selected)
    })
    
    output$GO_BP_dotPlot <- renderPlotly({
        GO_BP_enrichment_dp <- dotplot_plotly(GO_BP_enrichment_result(), 
                                              showCategory = input$GO_BP_dotPlot_showCategory,
                                              selected_rows = input$GO_BP_table_rows_selected)
    })
    
    output$GO_BP_emapPlot <- renderPlotly({
        GO_BP_enrichment_emap <- emap_plotly(GO_BP_enrichment_result(), 
                                             showCategory = input$GO_BP_emapPlot_showCategory, 
                                             selected_rows = input$GO_BP_table_rows_selected,
                                             line_scale=1, color="p.adjust")
    })
    
    output$GO_BP_cnetPlot <- renderPlotly({
        GO_BP_enrichment_cnet <- cnet_plotly(GO_BP_enrichment_result(),
                                             showCategory = input$GO_BP_cnetPlot_showCategory,
                                             selected_rows = input$GO_BP_table_rows_selected,
                                             foldChange = NULL) 
    })
    
    # Plot GO_MF results ------------------------------------------------------------------
    GO_MF_enrichment_result <- reactive({
        return (enrichGO(gene=gene_IDs(), pvalueCutoff=0.05, readable = TRUE,
                         OrgDb = org.Hs.eg.db, 
                         ont = "MF"))
    })
    
    output$GO_MF_table <- DT::renderDataTable(
        as.data.frame(GO_MF_enrichment_result()) %>% 
            mutate(pvalue = pvalue %>% signif(3),
                   p.adjust = p.adjust %>% signif(3),
                   qvalue = qvalue %>% signif(3))
    )
    
    output$GO_MF_barPlot <- renderPlotly({     
        GO_MF_enrichment_bp <- barplot_plotly(GO_MF_enrichment_result(), 
                                              showCategory = input$GO_MF_barPlot_showCategory,
                                              selected_rows = input$GO_MF_table_rows_selected)
    })
    
    output$GO_MF_dotPlot <- renderPlotly({
        GO_MF_enrichment_dp <- dotplot_plotly(GO_MF_enrichment_result(), 
                                              showCategory = input$GO_MF_dotPlot_showCategory,
                                              selected_rows = input$GO_MF_table_rows_selected)
    })
    
    output$GO_MF_emapPlot <- renderPlotly({
        GO_MF_enrichment_emap <- emap_plotly(GO_MF_enrichment_result(), 
                                             showCategory = input$GO_MF_emapPlot_showCategory, 
                                             selected_rows = input$GO_MF_table_rows_selected,
                                             line_scale=1, color="p.adjust")
    })
    
    output$GO_MF_cnetPlot <- renderPlotly({
        GO_MF_enrichment_cnet <- cnet_plotly(GO_MF_enrichment_result(),
                                             showCategory = input$GO_MF_cnetPlot_showCategory,
                                             selected_rows = input$GO_MF_table_rows_selected,
                                             foldChange = NULL) 
    })
    
    
    # Plot GO_CC results ------------------------------------------------------------------
    GO_CC_enrichment_result <- reactive({
        return (enrichGO(gene=gene_IDs(), pvalueCutoff=0.05, readable = TRUE,
                         OrgDb = org.Hs.eg.db, 
                         ont = "CC"))
    })
    
    output$GO_CC_table <- DT::renderDataTable(
        as.data.frame(GO_CC_enrichment_result()) %>% 
            mutate(pvalue = pvalue %>% signif(3),
                   p.adjust = p.adjust %>% signif(3),
                   qvalue = qvalue %>% signif(3))
    )
    
    output$GO_CC_barPlot <- renderPlotly({     
        GO_CC_enrichment_bp <- barplot_plotly(GO_CC_enrichment_result(), 
                                              showCategory = input$GO_CC_barPlot_showCategory,
                                              selected_rows = input$GO_CC_table_rows_selected)
    })
    
    output$GO_CC_dotPlot <- renderPlotly({
        GO_CC_enrichment_dp <- dotplot_plotly(GO_CC_enrichment_result(), 
                                              showCategory = input$GO_CC_dotPlot_showCategory,
                                              selected_rows = input$GO_CC_table_rows_selected)
    })
    
    output$GO_CC_emapPlot <- renderPlotly({
        GO_CC_enrichment_emap <- emap_plotly(GO_CC_enrichment_result(), 
                                             showCategory = input$GO_CC_emapPlot_showCategory, 
                                             selected_rows = input$GO_CC_table_rows_selected,
                                             line_scale=1, color="p.adjust")
    })
    
    output$GO_CC_cnetPlot <- renderPlotly({
        GO_CC_enrichment_cnet <- cnet_plotly(GO_CC_enrichment_result(),
                                             showCategory = input$GO_CC_cnetPlot_showCategory,
                                             selected_rows = input$GO_CC_table_rows_selected,
                                             foldChange = NULL) 
    })


})
