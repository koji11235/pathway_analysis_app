library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(viridis)

setwd("~/Project/20200719_Pathway_app")

sample_data<-read.csv("sample_gene_data.csv")

#サンプルデータの中のEntrez_Gene_IDを抽出
gene_IDs<-sample_data$Entrez_Gene_ID %>% as.character()

#KEGGエンリッチメント解析を実行
kegg_enrichment_result <- enrichKEGG(gene=gene_IDs,pvalueCutoff=0.05)

enrichment_result <- kegg_enrichment_result

#kegg_enrichment_result@result$Description <- str_wrap(kegg_enrichment_result@result$Description, width = 30)
bp <-barplot(kegg_enrichment_result, drop=TRUE, showCategory=8)#input$showCategory)
# 表示する内容 p.adjust, description, count

bp <-barplot(enrichment_result, drop=TRUE, showCategory=8)#input$showCategory)

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

barplot_plotly(enrichment_result, showCategory = 12)

#################################################################
# ggplotly ver
barplot_plotly <- function(enrichment_result, showCategory){
  bp <-barplot(enrichment_result, drop=TRUE, showCategory=showCategory)#input$showCategory)
  p <- ggplotly(
    bp$data %>% 
      mutate(Pathway = Description,
             Description = str_wrap(Description, width = 30)) %>% 
      ggplot(aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust,
                 text = paste("Pathway:", Pathway, "\n",
                              "p.adjust:", signif(p.adjust, 3), "\n",
                              "Count:", Count))) + 
      geom_bar(stat = "identity") + 
      theme_minimal() + 
      xlab("Gene Count") + 
      ylab("") + 
      viridis::scale_fill_viridis(option = "C") + #, guide = guide_colorbar(reverse = TRUE)) + # plotlyがcolor barの反転に対応してない
      theme(text = element_text(size=12),
            axis.text.x = element_text(size=12),
            axis.title.x = element_text(size=12),
            axis.text.y = element_text(size=12, lineheight=1.3)),
    tooltip = c("text")
  ) %>% 
    layout(autosize = TRUE, margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4))
  
  return(p)
}

barplot_plotly(enrichment_result, showCategory = 8)
  


g <- bp$data %>% 
  mutate(Pathway = Description,
         Description = str_wrap(Description, width = 30)) %>% 
  ggplot(aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust,
             text = paste("Pathway:", Pathway, "\n",
                          "p.adjust:", signif(p.adjust, 3), "\n",
                          "Count:", Count))) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  xlab("Gene Count") + 
  ylab("") + 
  #scale_fill_gradient(low = "red", high = "blue")+
  viridis::scale_fill_viridis(option = "C") + #, guide = guide_colorbar(reverse = TRUE)) + # plotlyがcolor barの反転に対応してない
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=12, lineheight=1.3))
  
ggplotly(g, tooltip = c("text"),dynamicTicks = )






###############################
# tests
g<- bp$data %>% 
  mutate(Pathway = Description,
         Description = str_wrap(Description, width = 30)) %>% 
  ggplot(aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust,
             text = paste("Pathway:", Pathway, "\n",
                          "p.adjust:", p.adjust, "\n",
                          "Count:", Count))) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  xlab("Gene Count") + 
  ylab("") + 
  scale_fill_gradient(low = "red", high = "blue")+
  scale_x_continuous(expand = c(0, 0)) + 
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=12, lineheight=1.3))
ggplotly(g, tooltip = c("text"))
  

ggplotly(
  barplot(kegg_enrichment_result, drop=TRUE, showCategory=8)#input$showCategory)
)

bp$data$Description %>% factor(levels = bp$data$Description)

g<- bp$data %>% 
  mutate(Pathway = Description,
         Description = str_wrap(Description, width = 30) %>% factor(levels = bp$data$Description)) %>% 
  # arrange(-p.adjust) %>% 
  ggplot(aes(x = Count, y = Description, fill = p.adjust,
             text=sprintf("Pathway: %s<br>", Pathway))) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  xlab("Gene Count") + 
  ylab("") + 
  scale_fill_gradient(low = "red", high = "blue")+
  scale_x_continuous(expand = c(0, 0)) + 
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=12, lineheight=1.3))
ggplotly(g, tooltip = c("text"))
ggplotly(g, tooltip = c("Pathway", "Count", "p.adjust"))

# https://stackoverflow.com/questions/34605919/formatting-mouse-over-labels-in-plotly-when-using-ggplotly

bp$data %>% 
  mutate(Description = str_wrap(Description, width = 30) %>% factor(levels = bp$data$Description)) %>% 
  # arrange(-p.adjust) %>% 
  ggplot(aes(x = Count, y = Description, fill = p.adjust,
             text=sprintf("letter: %s<br>Letter: %s", Count, Description))) + 
  geom_bar(stat = "identity") + 
  theme_minimal()


g<- bp$data %>% 
  mutate(Description = str_wrap(Description, width = 30) %>% factor(levels = bp$data$Description)) %>% 
  # arrange(-p.adjust) %>% 
  ggplot(aes(x = Count, y = Description, fill = p.adjust,
             text=sprintf("letter: %s<br>Letter: %s", Count, Description))) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  xlab("Gene Count") + 
  ylab("") + 
  scale_fill_gradient(low = "red", high = "blue")+
  scale_x_continuous(expand = c(0, 0)) + 
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=12, lineheight=1.3))
ggplotly(g, tooltip = c("Description", "Count", "p.adjust"))


g<- bp$data %>% 
  mutate(Pathway = Description,
         Description = str_wrap(Description, width = 30) %>% factor(levels = rev(bp$data$Description))) %>% 
  ggplot(aes(text=paste("Pathway:", Pathway, "\n",
                        "p.adjust:", p.adjust, "\n",
                        "Count:", Count))) + 
  geom_bar(aes(x = Count, y = Description, fill = p.adjust),
           stat = "identity") + 
  theme_minimal() + 
  xlab("Gene Count") + 
  ylab("") + 
  scale_fill_gradient(low = "red", high = "blue")+
  scale_x_continuous(expand = c(0, 0)) + 
  theme(text = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=12, lineheight=1.3))
ggplotly(g,)
ggplotly(g, tooltip = c("text"))



