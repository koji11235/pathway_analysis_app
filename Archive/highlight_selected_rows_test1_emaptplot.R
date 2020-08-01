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
# functions required for emap_plotly -----------------------------------
update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  
  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (nrow(x) < n) {
    n <- nrow(x)
  }
  
  return(n)
}

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}

emap_graph_build <- function(y,geneSets,color,line_scale) {
  if (is.null(dim(y)) | nrow(y) == 1) {
    g <- graph.empty(0, directed=FALSE)
    g <- add_vertices(g, nv = 1)
    V(g)$name <- as.character(y$Description)
    V(g)$color <- "red"
    ##return(ggraph(g))
  } else {
    id <- y[,"ID"]
    geneSets <- geneSets[id]
    n <- nrow(y) #
    w <- matrix(NA, nrow=n, ncol=n)
    colnames(w) <- rownames(w) <- y$Description
    
    for (i in seq_len(n-1)) {
      for (j in (i+1):n) {
        w[i,j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
      }
    }
    
    wd <- melt(w)
    wd <- wd[wd[,1] != wd[,2],]
    wd <- wd[!is.na(wd[,3]),]
    g <- graph.data.frame(wd[,-3], directed=FALSE)
    E(g)$width=sqrt(wd[,3] * 5) * line_scale 
    g <- delete_edges(g, E(g)[wd[,3] < 0.2])
    ## g <- delete.edges(g, E(g)[wd[,3] < 0.05])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == y$Description)))
    
    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt
    
    colVar <- y[idx, color]
    V(g)$color <- colVar
  }
  
  return(g)
}

add_segments_variable_width <- function(p, x, xend, y, yend, width, opacity){
  for (i in 1:length(x)){
    p <- p %>% 
      add_segments(data=tibble(x_column = x[i],
                               xend_column = xend[i],
                               y_column = y[i],
                               yend_column = yend[i],
                               width_column = width[i]), 
                   x=~x_column,  xend=~xend_column,
                   y=~y_column, yend=~yend_column, 
                   showlegend=FALSE,
                   color = I("#AAAAAA"),
                   opacity = opacity[i],
                   inherit = FALSE,
                   line = list(width = width[i])) 
  }
  return(p)
}


emap_plotly <- function(enrichment_result, showCategory = 5, line_scale = 1, color="p.adjust"){
  
  n <- update_n(enrichment_result, showCategory)
  geneSets <- geneInCategory(enrichment_result)
  
  y <- as.data.frame(enrichment_result)
  if (is.numeric(n)) {
    y <- y[1:n,]
  } else {
    y <- y[match(n, y$Description),]
    n <- length(n)
  }
  
  if (n == 0) {
    stop("no enriched term found...")
  }
  
  g <- emap_graph_build(y=y,geneSets=geneSets,color=color, line_scale=line_scale)
  
  node_coordinates <- layout.circle(g)
  
  label <- names(V(g))
  edge_list <- as.data.frame(get.edgelist(g))
  
  node2index <- c(1:length(V(g)))
  names(node2index) <- names(V(g))
  
  X_node <- node_coordinates[,1]
  Y_node <- node_coordinates[,2]
  edge_df <- tibble(
    x = X_node[node2index[as.character(edge_list$V1)]],
    xend = X_node[node2index[as.character(edge_list$V2)]],
    y = Y_node[node2index[as.character(edge_list$V1)]],
    yend = Y_node[node2index[as.character(edge_list$V2)]],
    width = E(g)$width
  )
  
  axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
  
  node_size_scale = c(min = 15, max = 50)
  #marker_size = (V(g)$size - min(V(g)$size)) / (max(V(g)$size) - min(V(g)$size)) * (marker_size_scale["max"] - marker_size_scale["min"]) + marker_size_scale["min"]
  # max minに合わせるのでなくてratioは保ってmaxを50に合わせる → この手法だと小さいものがあまりに小さくなってしまう
  # ratioは保ってmaxを50に合わせる、下は下限以下はすべてminにクリップする
  marker_size = pmax((V(g)$size / max(V(g)$size)) * marker_size_scale["max"], marker_size_scale["min"])
  colorbar_tickvals <- seq(floor(min(-log10(V(g)$color))), ceiling(max(-log10(V(g)$color))),  length=6)
  colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
  
  p <- plot_ly(x = ~X_node, y = ~Y_node,
               type = "scatter",
               mode = "markers",
               color = ~-log10(V(g)$color),
               colors = viridis_pal(option = "C", direction = 1)(3), 
               marker = list(
                 size = marker_size, 
                 alpha = 1,
                 line = list(width = 0)#,
               ),
               hovertemplate = paste(label, '\n',
                                     'p.adjust:', signif(V(g)$color, 3), '\n',
                                     'enriched genes:', V(g)$size),
               showlegend=FALSE)%>% 
    add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                                y=edge_df$y, yend=edge_df$yend, 
                                width = edge_df$width) %>% 
    add_markers(showlegend=FALSE) %>% 
    add_text(text = str_wrap(label, 30), textposition = "bottom", 
             textfont = list(size = 12, family = 'Helvetica') ) %>% 
    layout(
      autosize = TRUE, margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4),
      xaxis = list(range = c(min(X_node) - (max(X_node) - min(X_node)) * 0.25,
                             max(X_node) + (max(X_node) - min(X_node)) * 0.25), 
                   title = "", showgrid = FALSE, 
                   showticklabels = FALSE, zeroline = FALSE),
      yaxis = list(range = c(min(Y_node) - (max(Y_node) - min(Y_node)) * 0.15,
                             max(Y_node) + (max(X_node) - min(X_node)) * 0.1), 
                   title = "", showgrid = FALSE, 
                   showticklabels = FALSE, zeroline = FALSE)
    ) %>% 
    colorbar(title = "p.adjust",
             len = 0.5,
             tickmode = "array", 
             tickvals = colorbar_tickvals,
             ticktext = colorbar_ticktext, 
             exponentformat="e"
    )
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

emap_plotly <- function(enrichment_result, showCategory = 5, line_scale = 1, color="p.adjust", selected_rows=c()){
  
  selected_node_opacity <- 1
  not_selected_node_opacity <- 0.1
  
  selected_edge_opacity <- 1
  not_selected_edge_opacity <- 0.1
  
  # set actual_showCategory -----------------------------------------
  actual_showCategory <- update_n(enrichment_result, showCategory)
  geneSets <- geneInCategory(enrichment_result)
  
  enrichment_result_df <- as.data.frame(enrichment_result)
  if (is.numeric(actual_showCategory)) {
    enrichment_result_df <- enrichment_result_df[1:actual_showCategory,]
  } else {
    enrichment_result_df <- enrichment_result_df[match(actual_showCategory, enrichment_result_df$Description),]
    actual_showCategory <- length(actual_showCategory)
  }
  
  if (actual_showCategory == 0) {
    stop("no enriched term found...")
    # 一応空のplotlyオブジェクトを返すようにしておく
    return(plot_ly())
  }
  
  # greate emap graph -----------------------------------------
  g <- emap_graph_build(y=enrichment_result_df,geneSets=geneSets,color=color, line_scale=line_scale)
  
  # set node coordinates -----------------------------------------
  node_coordinates <- layout.circle(g)
  X_node <- node_coordinates[,1]
  Y_node <- node_coordinates[,2]
  
  # set node label -----------------------------------------
  label <- names(V(g))
  
  # set node size -----------------------------------------
  node_size_scale = c(min = 15, max = 50)
  #node_size = (V(g)$size - min(V(g)$size)) / (max(V(g)$size) - min(V(g)$size)) * (node_size_scale["max"] - node_size_scale["min"]) + node_size_scale["min"]
  # max minに合わせるのでなくてratioは保ってmaxを50に合わせる → この手法だと小さいものがあまりに小さくなってしまう
  # ratioは保ってmaxを50に合わせる、下は下限以下はすべてminにクリップする
  node_size = pmax((V(g)$size / max(V(g)$size)) * node_size_scale["max"], node_size_scale["min"])
  
  

  # set node opacity -----------------------------------------
    # retrieve selected rows -----------------------------------------
  if(length(selected_rows) == 0) {
    selected_discription <- enrichment_result@result$Description
  }else{
    selected_discription <- enrichment_result@result$Description[selected_rows]
  }
  # print(selected_discription)
  
  node_opacity <- if_else(names(V(g)) %in% selected_discription, selected_node_opacity, not_selected_node_opacity)
  
  
  # set edge arguments -----------------------------------------
  node2index <- c(1:length(V(g)))
  names(node2index) <- names(V(g))
  edge_list <- as.data.frame(get.edgelist(g))
  edge_df <- tibble(
    x = X_node[node2index[as.character(edge_list$V1)]],
    xend = X_node[node2index[as.character(edge_list$V2)]],
    y = Y_node[node2index[as.character(edge_list$V1)]],
    yend = Y_node[node2index[as.character(edge_list$V2)]],
    width = E(g)$width,
    opacity = if_else(
      (edge_list$V1 %in% selected_discription) | (edge_list$V2 %in% selected_discription),
      true = selected_edge_opacity,
      false = not_selected_edge_opacity
      )
  )
  
  # set colorbar tick  -----------------------------------------
  colorbar_tickvals <- seq(floor(min(-log10(V(g)$color))), ceiling(max(-log10(V(g)$color))),  length=6)
  colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
  
  
  axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
  p <- plot_ly(x = ~X_node, y = ~Y_node,
               type = "scatter",
               mode = "markers",
               color = ~-log10(V(g)$color),
               colors = viridis_pal(option = "C", direction = 1)(3), 
               marker = list(
                 size = node_size, 
                 opacity = node_opacity,
                 line = list(width = 0)#,
               ),
               hovertemplate = paste(label, '\n',
                                     'p.adjust:', signif(V(g)$color, 3), '\n',
                                     'enriched genes:', V(g)$size),
               showlegend=FALSE)%>% 
    add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                                y=edge_df$y, yend=edge_df$yend, 
                                width = edge_df$width,
                                opacity = edge_df$opacity) %>% 
    add_markers(showlegend=FALSE) %>% 
    add_text(text = str_wrap(label, 30), textposition = "bottom", 
             textfont = list(size = 12, family = 'Helvetica') ) %>% 
    layout(
      autosize = TRUE, margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4),
      xaxis = list(range = c(min(X_node) - (max(X_node) - min(X_node)) * 0.25,
                             max(X_node) + (max(X_node) - min(X_node)) * 0.25), 
                   title = "", showgrid = FALSE, 
                   showticklabels = FALSE, zeroline = FALSE),
      yaxis = list(range = c(min(Y_node) - (max(Y_node) - min(Y_node)) * 0.15,
                             max(Y_node) + (max(X_node) - min(X_node)) * 0.1), 
                   title = "", showgrid = FALSE, 
                   showticklabels = FALSE, zeroline = FALSE)
    ) %>% 
    colorbar(title = "p.adjust",
             len = 0.5,
             tickmode = "array", 
             tickvals = colorbar_tickvals,
             ticktext = colorbar_ticktext, 
             exponentformat="e"
    )
  return(p)
}



emap_plotly(enrichment_result, showCategory = 12, selected_rows = c())
emap_plotly(enrichment_result, showCategory = 12, selected_rows = c(1,2), line_scale = 3)
