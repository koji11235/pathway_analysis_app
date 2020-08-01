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
dotplot_plotly <- function(enrichment_result, showCategory){
  # set gene node size arguments ----------------------------------------
  gene_node_size <- 15
  node_size_scale <- c(min=10, max=40)
  
  # create dotplot object -----------------------------------------
  dp <-dotplot(enrichment_result, showCategory=showCategory)
  
  # create df for plot -----------------------------------------
  plotly_df <- dp$data %>% 
    mutate(Pathway = Description,
           Description = str_wrap(Description, width = 30), 
           node_size = (Count - min(dp$data$Count)) / max(dp$data$Count) * node_size_scale["max"] + node_size_scale["min"]) %>% 
    arrange(desc(GeneRatio))
  
  # set colorbar tick  -----------------------------------------
  colorbar_tickvals <- seq(floor(min(-log10(dp$data["p.adjust"]))), ceiling(max(-log10(dp$data["p.adjust"]))),  length=6)
  colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
  
  p <- plot_ly(
    x = plotly_df$GeneRatio, y = plotly_df$Description,
    type = "scatter",
    mode = "markers",
    color = ~-log10(plotly_df$p.adjust),
    colors = viridis_pal(option = "C", direction = 1)(3), 
    marker = list( 
      size = plotly_df$node_size,
      alpha = 1,
      line = list(width = 0)#,
    ),
    hovertemplate = paste(plotly_df$Description, '\n',
                          'p.adjust:', signif(plotly_df$p.adjust, 3), '\n',
                          'enriched genes:', plotly_df$Count, '\n',
                          "enriched gene ratio:", signif(plotly_df$GeneRatio, 2)),
    showlegend=FALSE) %>% 
    add_markers(alpha=1)
  p <- p %>% 
    layout(autosize = TRUE,
           margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4),
           xaxis = list(title = "Enriched Gene Ratio"),
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
selected_rows = c(1,2)

dotplot_plotly <- function(enrichment_result, showCategory, selected_rows=c()){
  # set gene node size arguments ----------------------------------------
  gene_node_size <- 15
  node_size_scale <- c(min=10, max=40)
  
  # create dotplot object -----------------------------------------
  dp <-dotplot(enrichment_result, showCategory=showCategory)
  
  # retrieve selected rows -----------------------------------------
  if(length(selected_rows) == 0) {
    selected_discription = dp$data$Description
  }else{
    selected_discription <- dp$data$Description[selected_rows]
  }
  #print(selected_discription)
  
  # create df for plot -----------------------------------------
  plotly_df <- dp$data %>% 
    mutate(Pathway = Description,
           Description = str_wrap(Description, width = 30), 
           node_size = (Count - min(dp$data$Count)) / max(dp$data$Count) * node_size_scale["max"] + node_size_scale["min"],
           Opacity = if_else(Pathway %in% selected_discription, 1, 0.1)) %>% 
    arrange(desc(GeneRatio))
  
  # set colorbar tick  -----------------------------------------
  colorbar_tickvals <- seq(floor(min(-log10(dp$data["p.adjust"]))), ceiling(max(-log10(dp$data["p.adjust"]))),  length=6)
  colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
  
  p <- plot_ly(
    x = plotly_df$GeneRatio, y = plotly_df$Description,
    type = "scatter",
    mode = "markers",
    color = ~-log10(plotly_df$p.adjust),
    colors = viridis_pal(option = "C", direction = 1)(3), 
    marker = list( 
      size = plotly_df$node_size,
      opacity = plotly_df$Opacity,
      line = list(width = 0)#,
    ),
    hovertemplate = paste(plotly_df$Description, '\n',
                          'p.adjust:', signif(plotly_df$p.adjust, 3), '\n',
                          'enriched genes:', plotly_df$Count, '\n',
                          "enriched gene ratio:", signif(plotly_df$GeneRatio, 2)),
    showlegend=FALSE) %>% 
    add_markers(alpha=1)
  p <- p %>% 
    layout(autosize = TRUE,
           margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4),
           xaxis = list(title = "Enriched Gene Ratio"),
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


dotplot_plotly(enrichment_result, showCategory = 12, selected_rows = c(5,7,8))
