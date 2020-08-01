library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(viridis)

dotplot_plotly <- function(enrichment_result, showCategory){
  dp <-dotplot(enrichment_result, showCategory=showCategory, orderBy = "GeneRatio")
  
  gene_node_size <- 15
  node_size_scale <- c(min=10, max=40)
  
  plotly_df <- dp$data %>% 
    mutate(Pathway = Description,
           Description = str_wrap(Description, width = 30), 
           node_size = (Count - min(dp$data$Count)) / max(dp$data$Count) * node_size_scale["max"] + node_size_scale["min"]) %>% 
    arrange(desc(GeneRatio))
  
  colorbar_tickvals <- seq(floor(min(-log10(dp$data["p.adjust"]))), ceiling(max(-log10(dp$data["p.adjust"]))),  length=6)
  colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
  
  p <- plot_ly(
    x = plotly_df$GeneRatio, y = plotly_df$Description,
    type = "scatter",
    mode = "markers",
    color = ~-log10(plotly_df$p.adjust),
    colors = viridis_pal(option = "C", direction = 1)(3), 
    alpha = 1,
    marker = list(
      size = plotly_df$node_size, 
      alpha = 1,
      line = list(width = 0)
    ),
    hovertemplate = paste(plotly_df$Description, '\n',
                          'p.adjust:', signif(plotly_df$p.adjust, 3), '\n',
                          'num genes:', plotly_df$Count, '\n',
                          "gene ratio:", signif(plotly_df$GeneRatio, 2)),
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

dotplot_plotly(enrichment_result, showCategory = 12)


#################################################################################
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

dotplot_plotly(enrichment_result, showCategory = 12)


#################################################################################
dotplot_plotly <- function(enrichment_result, showCategory){
  # create dotplot object -----------------------------------------
  dp <-dotplot(enrichment_result, showCategory=showCategory)
  
  # set colorbar tick  -----------------------------------------
  colorbar_tickvals <- seq(floor(min(-log10(dp$data["p.adjust"]))), ceiling(max(-log10(dp$data["p.adjust"]))),  length=6)
  colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
  
  p <- ggplotly(
    dp$data %>% 
      mutate(Pathway = Description,
             Description = str_wrap(Description, width = 30)) %>% 
      ggplot(aes(x = x, y = reorder(Description, x), 
                 color = -log10(p.adjust),
                 size = Count,
                 text = paste("Pathway:", Pathway, "\n",
                              "p.adjust:", signif(p.adjust, 3), "\n",
                              "Count:", Count, "\n",
                              "Gene Ratio:", signif(x, 3)))) + 
      geom_point() + 
      theme_minimal() + 
      xlab("Gene Ratio") + 
      ylab("") + 
      viridis::scale_color_viridis(option = "C", direction = 1,) + #, guide = guide_colorbar(reverse = TRUE)) + # plotlyがcolor barの反転に対応してない
      scale_size_continuous(range = c(3, 10)) + 
      #guides(colour = guide_colorbar(reverse = TRUE, title = "p.adjust", ) ,size = "none") +
      theme(text = element_text(size=12),
            axis.text.x = element_text(size=12),
            axis.title.x = element_text(size=12),
            axis.text.y = element_text(size=12, lineheight=1.3)),
    tooltip = c("text")) %>% 
    layout(autosize = TRUE, 
           margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4)) %>% 
    colorbar(title = "p.adjust",
             len = 0.5,
             tickmode = "array", 
             tickvals = colorbar_tickvals,
             ticktext = colorbar_ticktext, 
             exponentformat="e"#,reversescale=TRUE
    )
  return(p)
}

dotplot_plotly(enrichment_result, showCategory = 12)
  

dp <-dotplot(enrichment_result, showCategory=12)#input$showCategory)
ggplotly(dp)
p <- ggplotly(
  dp$data %>% 
    mutate(Pathway = Description,
           Description = str_wrap(Description, width = 30)) %>% 
    ggplot(aes(x = x, y = reorder(Description, x), 
               color = p.adjust,
               size = Count,
               text = paste("Pathway:", Pathway, "\n",
                            "p.adjust:", signif(p.adjust, 3), "\n",
                            "Count:", Count, "\n",
                            "Gene Ratio:", signif(x, 3)))) + 
    geom_point() + 
    theme_minimal() + 
    xlab("Gene Count") + 
    ylab("") + 
    viridis::scale_color_viridis(option = "C") + #, guide = guide_colorbar(reverse = TRUE)) + # plotlyがcolor barの反転に対応してない
    scale_size_continuous(range = c(3, 10)) + 
    theme(text = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.text.y = element_text(size=12, lineheight=1.3)),
  tooltip = c("text")) %>% 
  layout(autosize = TRUE, margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4))
p
