library(plotly)
library(shiny)
library(shinydashboard)
library(DT)

library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(DOSE)
library(reshape2)
library(igraph)
library(MetamapsDB)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(viridis)


setwd("~/Project/20200719_Pathway_app")

# previous version
barplot_plotly <- function(enrichment_result, showCategory){
  # create barplot object -----------------------------------------
  bp <-barplot(enrichment_result, showCategory=showCategory)
  
  # create df for plot -----------------------------------------
  plotly_df <- bp$data %>% 
    mutate(Pathway = Description,
           Description = str_wrap(Description, width = 30), 
    ) %>% 
    arrange(p.adjust)
  
  # set colorbar tick  -----------------------------------------
  colorbar_tickvals <- seq(floor(min(-log10(bp$data["p.adjust"]))), ceiling(max(-log10(bp$data["p.adjust"]))),  length=6)
  colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
  
  # plot  -----------------------------------------
  p <- plot_ly(
    x = plotly_df$Count, y = plotly_df$Description,
    type = "bar",
    color = ~-log10(plotly_df$p.adjust),
    colors = viridis_pal(option = "C", direction = 1)(3), 
    marker = list(
      size = plotly_df$node_size, 
      alpha = 1,
      line = list(width = 0)
    ),
    hovertemplate = paste(plotly_df$Description, '\n',
                          'p.adjust:', signif(plotly_df$p.adjust, 3), '\n',
                          'enriched genes:', plotly_df$Count, '\n',
                          "enriched gene ratio:", signif(plotly_df$GeneRatio, 2)),
    showlegend=FALSE) 
  p <- p %>% 
    layout(autosize = TRUE,
           margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4),
           xaxis = list(title = "Enriched Genes"),
           yaxis = list(categoryorder = "array",
                        categoryarray = rev(plotly_df$Descriptio)))%>% 
    colorbar(title = "p.adjust",
             len = 0.5,
             tickmode = "array",
             tickvals = colorbar_tickvals,
             ticktext = colorbar_ticktext,
             exponentformat="e")
  return(p)
}

# enrichment result preparation -------------------------------------------------------------
sample_data<-read.csv("sample_gene_data.csv")

#サンプルデータの中のEntrez_Gene_IDを抽出
gene_IDs<-sample_data$Entrez_Gene_ID %>% as.character()

#KEGGエンリッチメント解析を実行
kegg_enrichment_result <- enrichKEGG(gene=gene_IDs,pvalueCutoff=0.05)
reactome_enrichment_result <- enrichPathway(gene=gene_IDs,pvalueCutoff=0.05, readable=T)


#cn <- cnetplot(reactome_enrichment_result, showCategory = 2)
cn <- cnetplot(kegg_enrichment_result, showCategory = 2, circular=T)
cn <- cnetplot(kegg_enrichment_result, showCategory = 2, circular=F, colorEdge = T)
cn <- cnetplot(reactome_enrichment_result, showCategory = 2,  node_label = "gene")

#em <- emapplot(reactome_enrichment_result, showCategory = 5)
#emapplot(kegg_enrichment_result, showCategory = 5)
cnetplot.enrichResult(kegg_enrichment_result, showCategory = 2,  node_label = "none")

enrichment_result <- kegg_enrichment_result

# row selection test ----------------
showCategory = 5
selected_rows = c()

barplot_plotly <- function(enrichment_result, showCategory, selected_rows){
  # create barplot object -----------------------------------------
  bp <-barplot(enrichment_result, showCategory=showCategory)
  
  # retrieve selected rows -----------------------------------------
  if(length(selected_rows) == 0) {
    selected_discription = bp$data$Description
  }else{
    selected_discription <- bp$data$Description[selected_rows]
  }
  # print(selected_discription)
  
  # create df for plot -----------------------------------------
  plotly_df <- bp$data %>% 
    mutate(Pathway = Description,
           Description = str_wrap(Description, width = 30), 
           Opacity = if_else(Pathway %in% selected_discription, 1, 0.2)
    ) %>% 
    arrange(p.adjust)
  
  # set colorbar tick  -----------------------------------------
  colorbar_tickvals <- seq(floor(min(-log10(bp$data["p.adjust"]))), ceiling(max(-log10(bp$data["p.adjust"]))),  length=6)
  colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
  
  # plot  -----------------------------------------
  p <- plot_ly(
    x = plotly_df$Count, y = plotly_df$Description,
    type = "bar",
    color = ~-log10(plotly_df$p.adjust),
    colors = viridis_pal(option = "C", direction = 1)(3), 
    marker = list(
      size = plotly_df$node_size, 
      opacity = plotly_df$Opacity,
      line = list(width = 0)
    ),
    hovertemplate = paste(plotly_df$Description, '\n',
                          'p.adjust:', signif(plotly_df$p.adjust, 3), '\n',
                          'enriched genes:', plotly_df$Count, '\n',
                          "enriched gene ratio:", signif(plotly_df$GeneRatio, 2)),
    showlegend=FALSE) 
  p <- p %>% 
    layout(autosize = TRUE,
           margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4),
           xaxis = list(title = "Enriched Genes"),
           yaxis = list(categoryorder = "array",
                        categoryarray = rev(plotly_df$Descriptio)))%>% 
    colorbar(title = "p.adjust",
             len = 0.5,
             tickmode = "array",
             tickvals = colorbar_tickvals,
             ticktext = colorbar_ticktext,
             exponentformat="e")
  return(p)
}


barplot_plotly(enrichment_result, showCategory = 12, selected_rows = c(1,4))
