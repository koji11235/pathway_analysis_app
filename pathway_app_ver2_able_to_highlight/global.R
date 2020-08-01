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


# barplot function  -----------------------------------
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

# previous ver 
# barplot_plotly <- function(enrichment_result, showCategory){
#   # create barplot object -----------------------------------------
#   bp <-barplot(enrichment_result, showCategory=showCategory)
#   
#   # create df for plot -----------------------------------------
#   plotly_df <- bp$data %>% 
#     mutate(Pathway = Description,
#            Description = str_wrap(Description, width = 30), 
#     ) %>% 
#     arrange(p.adjust)
#   
#   # set colorbar tick  -----------------------------------------
#   colorbar_tickvals <- seq(floor(min(-log10(bp$data["p.adjust"]))), ceiling(max(-log10(bp$data["p.adjust"]))),  length=6)
#   colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
#   
#   # plot  -----------------------------------------
#   p <- plot_ly(
#     x = plotly_df$Count, y = plotly_df$Description,
#     type = "bar",
#     color = ~-log10(plotly_df$p.adjust),
#     colors = viridis_pal(option = "C", direction = 1)(3), 
#     marker = list(
#       size = plotly_df$node_size, 
#       alpha = 1,
#       line = list(width = 0)
#     ),
#     hovertemplate = paste(plotly_df$Description, '\n',
#                           'p.adjust:', signif(plotly_df$p.adjust, 3), '\n',
#                           'enriched genes:', plotly_df$Count, '\n',
#                           "enriched gene ratio:", signif(plotly_df$GeneRatio, 2)),
#     showlegend=FALSE) 
#   p <- p %>% 
#     layout(autosize = TRUE,
#            margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4),
#            xaxis = list(title = "Enriched Genes"),
#            yaxis = list(categoryorder = "array",
#                         categoryarray = rev(plotly_df$Descriptio)))%>% 
#     colorbar(title = "p.adjust",
#              len = 0.5,
#              tickmode = "array",
#              tickvals = colorbar_tickvals,
#              ticktext = colorbar_ticktext,
#              exponentformat="e")
#   return(p)
# }


# dotplot function  -----------------------------------
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
  # print(selected_discription)
  
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

## previous ver
# dotplot_plotly <- function(enrichment_result, showCategory){
#   # set gene node size arguments ----------------------------------------
#   gene_node_size <- 15
#   node_size_scale <- c(min=10, max=40)
#   
#   # create dotplot object -----------------------------------------
#   dp <-dotplot(enrichment_result, showCategory=showCategory)
#   
#   # create df for plot -----------------------------------------
#   plotly_df <- dp$data %>% 
#     mutate(Pathway = Description,
#            Description = str_wrap(Description, width = 30), 
#            node_size = (Count - min(dp$data$Count)) / max(dp$data$Count) * node_size_scale["max"] + node_size_scale["min"]) %>% 
#     arrange(desc(GeneRatio))
#   
#   # set colorbar tick  -----------------------------------------
#   colorbar_tickvals <- seq(floor(min(-log10(dp$data["p.adjust"]))), ceiling(max(-log10(dp$data["p.adjust"]))),  length=6)
#   colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
#   
#   p <- plot_ly(
#     x = plotly_df$GeneRatio, y = plotly_df$Description,
#     type = "scatter",
#     mode = "markers",
#     color = ~-log10(plotly_df$p.adjust),
#     colors = viridis_pal(option = "C", direction = 1)(3), 
#     marker = list( 
#       size = plotly_df$node_size,
#       alpha = 1,
#       line = list(width = 0)#,
#     ),
#     hovertemplate = paste(plotly_df$Description, '\n',
#                           'p.adjust:', signif(plotly_df$p.adjust, 3), '\n',
#                           'enriched genes:', plotly_df$Count, '\n',
#                           "enriched gene ratio:", signif(plotly_df$GeneRatio, 2)),
#     showlegend=FALSE) %>% 
#     add_markers(alpha=1)
#   p <- p %>% 
#     layout(autosize = TRUE,
#            margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4),
#            xaxis = list(title = "Enriched Gene Ratio"),
#            yaxis = list(categoryorder = "array",
#                         categoryarray = rev(plotly_df$Descriptio)))%>% 
#     colorbar(title = "p.adjust",
#              len = 0.5,
#              tickmode = "array",
#              tickvals = colorbar_tickvals,
#              ticktext = colorbar_ticktext,
#              exponentformat="e")
#   return(p)
# }




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

# previous ver
# emap_plotly <- function(enrichment_result, showCategory = 5, line_scale = 1, color="p.adjust"){
#   
#   n <- update_n(enrichment_result, showCategory)
#   geneSets <- geneInCategory(enrichment_result)
#   
#   y <- as.data.frame(enrichment_result)
#   if (is.numeric(n)) {
#     y <- y[1:n,]
#   } else {
#     y <- y[match(n, y$Description),]
#     n <- length(n)
#   }
#   
#   if (n == 0) {
#     stop("no enriched term found...")
#   }
#   
#   g <- emap_graph_build(y=y,geneSets=geneSets,color=color, line_scale=line_scale)
#   
#   node_coordinates <- layout.circle(g)
#   
#   label <- names(V(g))
#   edge_list <- as.data.frame(get.edgelist(g))
#   
#   node2index <- c(1:length(V(g)))
#   names(node2index) <- names(V(g))
#   
#   X_node <- node_coordinates[,1]
#   Y_node <- node_coordinates[,2]
#   edge_df <- tibble(
#     x = X_node[node2index[as.character(edge_list$V1)]],
#     xend = X_node[node2index[as.character(edge_list$V2)]],
#     y = Y_node[node2index[as.character(edge_list$V1)]],
#     yend = Y_node[node2index[as.character(edge_list$V2)]],
#     width = E(g)$width
#   )
#   
#   axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
#   
#   marker_size_scale = c(min = 15, max = 50)
#   #marker_size = (V(g)$size - min(V(g)$size)) / (max(V(g)$size) - min(V(g)$size)) * (marker_size_scale["max"] - marker_size_scale["min"]) + marker_size_scale["min"]
#   # max minに合わせるのでなくてratioは保ってmaxを50に合わせる → この手法だと小さいものがあまりに小さくなってしまう
#   # ratioは保ってmaxを50に合わせる、下は下限以下はすべてminにクリップする
#   marker_size = pmax((V(g)$size / max(V(g)$size)) * marker_size_scale["max"], marker_size_scale["min"])
#   colorbar_tickvals <- seq(floor(min(-log10(V(g)$color))), ceiling(max(-log10(V(g)$color))),  length=6)
#   colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
#   
#   p <- plot_ly(x = ~X_node, y = ~Y_node,
#                type = "scatter",
#                mode = "markers",
#                color = ~-log10(V(g)$color),
#                colors = viridis_pal(option = "C", direction = 1)(3), 
#                marker = list(
#                  size = marker_size, 
#                  alpha = 1,
#                  line = list(width = 0)#,
#                ),
#                hovertemplate = paste(label, '\n',
#                                      'p.adjust:', signif(V(g)$color, 3), '\n',
#                                      'enriched genes:', V(g)$size),
#                showlegend=FALSE)%>% 
#     add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
#                                 y=edge_df$y, yend=edge_df$yend, 
#                                 width = edge_df$width) %>% 
#     add_markers(showlegend=FALSE) %>% 
#     add_text(text = str_wrap(label, 30), textposition = "bottom", 
#              textfont = list(size = 12, family = 'Helvetica') ) %>% 
#     layout(
#       autosize = TRUE, margin = list(l = 0, r = 0, b = 0, t = 20, pad = 4),
#       xaxis = list(range = c(min(X_node) - (max(X_node) - min(X_node)) * 0.25,
#                              max(X_node) + (max(X_node) - min(X_node)) * 0.25), 
#                    title = "", showgrid = FALSE, 
#                    showticklabels = FALSE, zeroline = FALSE),
#       yaxis = list(range = c(min(Y_node) - (max(Y_node) - min(Y_node)) * 0.15,
#                              max(Y_node) + (max(X_node) - min(X_node)) * 0.1), 
#                    title = "", showgrid = FALSE, 
#                    showticklabels = FALSE, zeroline = FALSE)
#     ) %>% 
#     colorbar(title = "p.adjust",
#              len = 0.5,
#              tickmode = "array", 
#              tickvals = colorbar_tickvals,
#              ticktext = colorbar_ticktext, 
#              exponentformat="e"
#     )
#   return(p)
# }



# functions required for cnet_plotly-----------------------------------
extract_geneSets <- function(x, n) {
  n <- update_n(x, n)
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  geneSets <- geneSets[y$ID]
  names(geneSets) <- y$Description
  if (is.numeric(n)) {
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}

list2graph <- function(inputList) {
  x <- list2df(inputList)
  g <- graph.data.frame(x, directed=FALSE)
  return(g)
}


list2df <- function(inputList) {
  ldf <- lapply(1:length(inputList), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })
  
  do.call('rbind', ldf)
}

fc_readable <- function(x, foldChange = NULL) {
  if (is.null(foldChange))
    return(NULL)
  
  if(x@readable) {
    gid <- names(foldChange)
    if (is(x, 'gseaResult')) {
      ii <- gid %in% names(x@geneList)
    } else {
      ii <- gid %in% x@gene
    }
    gid[ii] <- x@gene2Symbol[gid[ii]]
    names(foldChange) <- gid
  }
  return(foldChange)
}

add_segments_different_color <- function(p, x, xend, y, yend, color, opacity, line_width=1){
  for (i in 1:length(x)){
    p <- p %>% 
      add_segments(data=tibble(x_column = x[i],
                               xend_column = xend[i],
                               y_column = y[i],
                               yend_column = yend[i]), 
                   x=~x_column,  xend=~xend_column,
                   y=~y_column, yend=~yend_column, 
                   showlegend=FALSE, 
                   opacity = opacity[i],
                   color = I(color[i]),
                   inherit = FALSE,
                   line = list(width = line_width)) 
  }
  return(p)
}
# previous version
# add_segments_different_color <- function(p, x, xend, y, yend, color, line_width=1){
#   for (i in 1:length(x)){
#     p <- p %>% 
#       add_segments(data=tibble(x_column = x[i],
#                                xend_column = xend[i],
#                                y_column = y[i],
#                                yend_column = yend[i]), 
#                    x=~x_column,  xend=~xend_column,
#                    y=~y_column, yend=~yend_column, 
#                    showlegend=FALSE, 
#                    alpha = 0.5,
#                    color = I(color[i]),
#                    inherit = FALSE,
#                    line = list(width = line_width)) 
#   }
#   return(p)
# }


# calculate_gene_color <- function(categories, categories2color){
#   result_color <- "#FFFFFF"
#   minimum_rgb_value = c(red=54, green=54, blue=54)
#   hsvs = lapply(categories2color, function(color){rgb2hsv(col2rgb(color))}) 
#   v_list =c()
#   for (hsv in hsvs){
#     v_list <- c(v_list, hsv["v",1])
#   }
#   mean_v = mean(v_list)
#   
#   for (category in categories){
#     result_color_rgb <- mapply(min, 
#                                col2rgb(categories2color[[category]])[,1], 
#                                col2rgb(result_color)[,1])
#     result_color <- rgb(red = result_color_rgb["red"], 
#                         green = result_color_rgb["green"], 
#                         blue = result_color_rgb["blue"], 
#                         maxColorValue = 255)
#   }
#   result_color_rgb <- mapply(max, 
#                              minimum_rgb_value, 
#                              col2rgb(result_color)[,1])
#   result_color <- rgb(red = result_color_rgb["red"], 
#                       green = result_color_rgb["green"], 
#                       blue = result_color_rgb["blue"], 
#                       maxColorValue = 255)
#   
#   if (length(categories) > 1){
#     result_color_hsv <- col2rgb(result_color) %>% rgb2hsv() 
#     result_color_hsv["v",1] <- mean_v
#     result_color <- hsv(result_color_hsv["h",1], result_color_hsv["s",1],result_color_hsv["v",1])
#   }
#   return(result_color)
# }

reorder_palette <- function(palette, showCategory){
  reordered <- c()
  separated_palettes <- lapply(seq_len(showCategory), function(i){NULL})
  
  block_size <- ceiling(length(palette) / showCategory)
  for (i in seq_len(showCategory)){
    separated_palettes[[i]] <- palette[((i - 1) * block_size + 1) : (i * block_size)]
  }
  
  for (i in seq_len(block_size)){
    for (separated_palette in separated_palettes){
      reordered <- c(reordered, separated_palette[[i]])
    }
  }
  reordered <- reordered[!is.na(reordered)]
  return(reordered)
}
cnet_plotly <- function(enrichment_result, showCategory, selected_rows = c(), foldChange=NULL, is_kegg=FALSE){
  
  node_label_wrap_width <- 40
  gene_font_size <- 8
  category_font_size <- 12
  gene_node_size <- 10
  node_size_scale <- c(gene_node_size=gene_node_size, category_size_min=20, category_size_max=40)
  
  gene_node_default_opacity <- 0.3
  geneset_node_default_opacity <- 1.0
  selected_gene_node_opacity <- 0.5
  selected_geneset_node_opacity <- 1.0
  not_selected_node_opacity <- 0.1
  
  node_line_width <- 2
  node_line_color <- "#FFFFFF"
  
  edge_width <- 2
  selected_edge_opacity <- 0.6
  not_selected_edge_opacity <- 0.1
  
  
  # set actual_showCategory -----------------------------------------
  actual_showCategory <- update_n(enrichment_result, showCategory)
  # category_color_palette <- brewer.pal(actual_showCategory, color_palette)
  # node_label <- match.arg(node_label, c("category", "gene", "all", "none"))
  foldChange <- fc_readable(enrichment_result, foldChange)
  
  # create igraph object----------------------------------------
  geneSets <- extract_geneSets(enrichment_result, actual_showCategory)
  g <- list2graph(geneSets)
  
  # set node label----------------------------------------
  node_names <- names(V(g))
  if (is_kegg){
    node_gene_symbols <- mapIds(org.Hs.eg.db, keys=node_names[(actual_showCategory+1):length(node_names)], column="SYMBOL", keytype="ENTREZID", multiVals="first")
    V(g)$label <- str_wrap(c(node_names[1:actual_showCategory], node_gene_symbols), width = 40)
  }else{
    V(g)$label <- str_wrap(node_names, width = node_label_wrap_width)
  }
  
  # set node size ----------------------------------------
  #(gene_set_sizes - min(gene_set_sizes)) / (max(gene_set_sizes) - min(gene_set_sizes)) * (node_size_scale[["category_size_max"]] - node_size_scale[["category_size_min"]]) + node_size_scale[["category_size_min"]]
  gene_set_sizes <- sapply(geneSets, length)
  gene_set_sizes_scaled <- (gene_set_sizes - min(gene_set_sizes)) / (max(gene_set_sizes) - min(gene_set_sizes)) * (node_size_scale[["category_size_max"]] - node_size_scale[["category_size_min"]]) + node_size_scale[["category_size_min"]]
  V(g)$size <- gene_node_size
  V(g)$size[1:actual_showCategory] <- gene_set_sizes_scaled
  
  # set node text font size ----------------------------------------
  V(g)$font_size <- gene_font_size
  V(g)$font_size[1:actual_showCategory] <- category_font_size
  
  
  # set node categorie ----------------------------------------
  gene2category_df <- as.data.frame(get.edgelist(g)) %>% 
    group_by(V2) %>%
    summarise(categories = str_c(V1, collapse  = ";"),
              n_categories = n())
  categories_with_combination <- gene2category_df %>% 
    pull(categories) %>% unique()
  categories_combination_only <- categories_with_combination[!(categories_with_combination %in% names(geneSets))]
  
  gene2category_df <- gene2category_df %>% 
    mutate(category_id = gene2category_df %>% 
             pull(categories) %>% 
             factor(levels = c(levels(factor(names(geneSets))), categories_combination_only)) %>% 
             as.numeric()
    )
  gene2category <- gene2category_df %>% pull(categories)
  names(gene2category) <- gene2category_df %>% pull(V2)
  V(g)$category <- c(names(geneSets), gene2category[names(V(g))[(actual_showCategory+1):length(V(g))]])
  
  #head(V(g)$category)
  #head(V(g)$label)
  
  
  # set node color ----------------------------------------
  if (!is.null(foldChange)) {
    fc <- foldChange[V(g)$name[(actual_showCategory+1):length(V(g))]]
    V(g)$color <- NA
    V(g)$color[(actual_showCategory+1):length(V(g))] <- fc
  } else {
    category_color_palette <- viridis_pal(option = "C", direction = -1)(max(unique(gene2category_df$category_id))) 
    category_color_palette <- reorder_palette(palette = category_color_palette, showCategory = actual_showCategory)
    
    V(g)$category_id <- V(g)$category %>% 
      factor(levels = c(levels(factor(names(geneSets))), categories_combination_only)) %>% 
      as.numeric()
    V(g)$color <- category_color_palette[V(g)$category_id]
    
    # gene2color_df <- as.data.frame(get.edgelist(g)) %>% 
    #   group_by(V2) %>%
    #   summarise(categories = str_c(V1, collapse  = ";"),
    #             n_categories = n())
    # categories_with_combination <- gene2color_df %>% 
    #   pull(categories) %>% unique()
    # categories_combination_only <- categories_with_combination[!(categories_with_combination %in% names(geneSets))]
    # 
    # gene2color_df <- gene2color_df %>% 
    #   mutate(category_id = gene2color_df %>% 
    #            pull(categories) %>% 
    #            factor(levels = c(levels(factor(names(geneSets))), categories_combination_only)) %>% 
    #            as.numeric()
    #   )
    # category_color_palette <- viridis_pal(option = "C", direction = -1)(max(unique(gene2color_df$category_id))) 
    # category_color_palette <- reorder_palette(palette = category_color_palette, showCategory = actual_showCategory)
    # 
    # gene2color <- category_color_palette[gene2color_df %>% pull(category_id)]
    # names(gene2color) <- gene2color_df %>% pull(V2)
    # V(g)$color <- c(category_color_palette[factor(names(geneSets))], gene2color[names(V(g))[(actual_showCategory+1):length(V(g))]])
  }
  
  # set node opacity ----------------------------------------
  if(length(selected_rows) == 0) {
    selected_discription <- enrichment_result@result$Description
    
    V(g)$opacity <- gene_node_default_opacity
    V(g)$opacity[1:actual_showCategory] <- geneset_node_default_opacity
  }else{
    selected_discription <- enrichment_result@result$Description[selected_rows]
    #V(g)$opacity <- if_else(V(g)$category %in% selected_discription, 1, 0.1)
    # V(g)$opacity <- if_else(
    #   sapply(V(g)$category, 
    #          function(categories){
    #            category_vector <- str_split(categories, ";")[[1]]
    #            return(any(category_vector %in% selected_discription))}),
    #   true=1, false=0.1
    # )
    V(g)$opacity <- if_else(
      sapply(V(g)$category,
             function(categories){
               category_vector <- str_split(categories, ";")[[1]]
               return(any(category_vector %in% selected_discription))}),
      true=if_else(V(g)$size == gene_node_size,
                   selected_gene_node_opacity,
                   selected_geneset_node_opacity),
      false=not_selected_node_opacity
    )
  }
  
  
  # set edge color ----------------------------------------
  node_coordinates <- layout_with_kk(g)
  
  node2index <- c(1:length(V(g)))
  names(node2index) <- names(V(g))
  
  X_node <- node_coordinates[,1]
  Y_node <- node_coordinates[,2]
  
  edge_list <- as.data.frame(get.edgelist(g))
  
  edge_df <- tibble(
    x = X_node[node2index[as.character(edge_list$V1)]],
    xend = X_node[node2index[as.character(edge_list$V2)]],
    y = Y_node[node2index[as.character(edge_list$V1)]],
    yend = Y_node[node2index[as.character(edge_list$V2)]],
    category = edge_list$V1,
    category_color = category_color_palette[as.numeric(factor(edge_list$V1))],
    opacity = if_else(edge_list$V1 %in% selected_discription, 
                      selected_edge_opacity, 
                      not_selected_edge_opacity)
  )
  
  # set hover info text ----------------------------------------
  hoverinfo_text <- tibble(node_name = names(V(g)),
                           node_label = V(g)$label) %>% 
    left_join(as.data.frame(get.edgelist(g)) %>% 
                group_by(V2) %>%
                summarise(categories = str_c(V1, collapse  = "\n"),
                          n_categories = n()),
              by = c("node_name" = "V2")) %>% 
    left_join(enrichment_result@result,
              by = c("node_name"="Description")) %>% 
    mutate(hoverinfo = if_else(is.na(categories), 
                               paste0(node_label, "\n",
                                      "p.adjust: ", signif(p.adjust, 3), "\n",
                                      "enriched genes: ", Count),
                               paste0(node_label, "\n",
                                      "belongs to ", n_categories, if_else(n_categories == 1, "category", "categories"), "\n",
                                      "category:", "\n", 
                                      categories, "\n"))) %>% 
    pull(hoverinfo)
  
  # plot cnet plot by plotly ----------------------------------------
  axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
  
  #print(V(g)$opacity)
  
  p <- plot_ly(x = ~X_node, y = ~Y_node,
               type = "scatter",
               mode = "markers",
               text = V(g)$label,
               marker = list(color = ~V(g)$color, 
                             opacity = ~V(g)$opacity,
                             size = ~V(g)$size,
                             line = list(width = node_line_width, color=node_line_color)),
               hovertemplate = hoverinfo_text,
               showlegend=F) %>% 
    add_segments_different_color(x=edge_df$x,  xend=edge_df$xend,
                                 y=edge_df$y, yend=edge_df$yend, 
                                 color = edge_df$category_color,
                                 opacity = edge_df$opacity,
                                 line_width = edge_width) %>% 
    add_markers() %>% 
    add_text(text = V(g)$label, 
             textposition = "center",
             #opacity = V(g)$opacity,
             textfont = list(
               size = V(g)$font_size, 
               family = 'Helvetica')) %>% 
    layout(xaxis = axis,yaxis = axis)
  return(p)
}


# previous version
# cnet_plotly <- function(enrichment_result, showCategory, foldChange=NULL, is_kegg=FALSE){
#   gene_font_size <- 8
#   category_font_size <- 12
#   gene_node_size <- 10
#   node_size_scale <- c(gene_node_size=gene_node_size, category_size_min=20, category_size_max=40)
#   
#   actual_showCategory <- update_n(enrichment_result, showCategory)
#   # category_color_palette <- brewer.pal(actual_showCategory, color_palette)
#   # node_label <- match.arg(node_label, c("category", "gene", "all", "none"))
#   foldChange <- fc_readable(enrichment_result, foldChange)
#   
#   # create igraph object----------------------------------------
#   geneSets <- extract_geneSets(enrichment_result, actual_showCategory)
#   g <- list2graph(geneSets)
#   
#   # set node label----------------------------------------
#   node_names <- names(V(g))
#   if (is_kegg){
#     node_gene_symbols <- mapIds(org.Hs.eg.db, keys=node_names[(actual_showCategory+1):length(node_names)], column="SYMBOL", keytype="ENTREZID", multiVals="first")
#     V(g)$label <- str_wrap(c(node_names[1:actual_showCategory], node_gene_symbols), width = 40)
#   }else{
#     V(g)$label <- str_wrap(node_names, width = 40)
#   }
#   # set node size ----------------------------------------
#   #(gene_set_sizes - min(gene_set_sizes)) / (max(gene_set_sizes) - min(gene_set_sizes)) * (node_size_scale[["category_size_max"]] - node_size_scale[["category_size_min"]]) + node_size_scale[["category_size_min"]]
#   
#   
#   gene_set_sizes <- sapply(geneSets, length)
#   gene_set_sizes_scaled <- (gene_set_sizes - min(gene_set_sizes)) / (max(gene_set_sizes) - min(gene_set_sizes)) * (node_size_scale[["category_size_max"]] - node_size_scale[["category_size_min"]]) + node_size_scale[["category_size_min"]]
#   V(g)$size <- gene_node_size
#   V(g)$size[1:actual_showCategory] <- gene_set_sizes_scaled
#   
#   # set node text font size ----------------------------------------
#   V(g)$font_size <- gene_font_size
#   V(g)$font_size[1:actual_showCategory] <- category_font_size
#   
#   
#   # set gene node color ----------------------------------------
#   if (!is.null(foldChange)) {
#     fc <- foldChange[V(g)$name[(actual_showCategory+1):length(V(g))]]
#     V(g)$color <- NA
#     V(g)$color[(actual_showCategory+1):length(V(g))] <- fc
#   } else {
#     gene2color_df <- as.data.frame(get.edgelist(g)) %>% 
#       group_by(V2) %>%
#       summarise(categories = str_c(V1, collapse  = ";"),
#                 n_categories = n())
#     categories_with_combination <- gene2color_df %>% 
#       pull(categories) %>% unique()
#     categories_combination_only <- categories_with_combination[!(categories_with_combination %in% names(geneSets))]
#     
#     gene2color_df <- gene2color_df %>% 
#       mutate(category_id = gene2color_df %>% 
#                pull(categories) %>% 
#                factor(levels = c(levels(factor(names(geneSets))), categories_combination_only)) %>% 
#                as.numeric()
#       )
#     category_color_palette <- viridis_pal(option = "C", direction = -1)(max(unique(gene2color_df$category_id))) 
#     category_color_palette <- reorder_palette(palette = category_color_palette, showCategory = actual_showCategory)
#     
#     gene2color <- category_color_palette[gene2color_df %>% pull(category_id)]
#     names(gene2color) <- gene2color_df %>% pull(V2)
#     V(g)$color <- c(category_color_palette[factor(names(geneSets))], gene2color[names(V(g))[(actual_showCategory+1):length(V(g))]])
#     
#   }
#   
#   # set edge color ----------------------------------------
#   node_coordinates <- layout_with_kk(g)
#   
#   node2index <- c(1:length(V(g)))
#   names(node2index) <- names(V(g))
#   
#   X_node <- node_coordinates[,1]
#   Y_node <- node_coordinates[,2]
#   
#   edge_list <- as.data.frame(get.edgelist(g))
#   
#   edge_df <- tibble(
#     x = X_node[node2index[as.character(edge_list$V1)]],
#     xend = X_node[node2index[as.character(edge_list$V2)]],
#     y = Y_node[node2index[as.character(edge_list$V1)]],
#     yend = Y_node[node2index[as.character(edge_list$V2)]],
#     category = edge_list$V1,
#     category_color = category_color_palette[as.numeric(factor(edge_list$V1))]
#   )
#   
#   # set hover info text ----------------------------------------
#   as.data.frame(get.edgelist(g)) %>% 
#     group_by(V2) %>%
#     summarise(categories = str_c(V1, collapse  = "; "),
#               n_categories = n())
#   
#   hoverinfo_text <- tibble(node_name = names(V(g)),
#                            node_label = V(g)$label) %>% 
#     left_join(as.data.frame(get.edgelist(g)) %>% 
#                 group_by(V2) %>%
#                 summarise(categories = str_c(V1, collapse  = "\n"),
#                           n_categories = n()),
#               by = c("node_name" = "V2")
#     ) %>% left_join(enrichment_result@result,
#                     by = c("node_name"="Description")) %>% 
#     mutate(hoverinfo = if_else(is.na(categories), 
#                                paste0(node_label, "\n",
#                                       "p.adjust: ", signif(p.adjust, 3), "\n",
#                                       "enriched genes: ", Count),
#                                paste0(node_label, "\n",
#                                       "belongs to ", n_categories, if_else(n_categories == 1, "category", "categories"), "\n",
#                                       "category:", "\n", 
#                                       categories, "\n")
#     )
#     ) %>% pull(hoverinfo)
#   
#   # plot cnet plot by plotly ----------------------------------------
#   axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
#   
#   p <- plot_ly(x = ~X_node, y = ~Y_node,
#                type = "scatter",
#                mode = "markers",
#                text = V(g)$label,
#                #color = ~-log10(V(g)$color),
#                #colors = viridis_pal(option = "C")(3), 
#                
#                marker = list(color = V(g)$color, 
#                              size = V(g)$size,
#                              line = list(width = 0)),
#                #hoverinfo = "text",
#                hovertemplate = hoverinfo_text,
#                showlegend=FALSE) %>% 
#     add_segments_different_color(x=edge_df$x,  xend=edge_df$xend,
#                                  y=edge_df$y, yend=edge_df$yend, 
#                                  color = edge_df$category_color,
#                                  line_width = 2) %>% 
#     add_markers(showlegend=T) %>% 
#     add_text(text = V(g)$label, textposition = "center",
#              textfont = list(size = V(g)$font_size, family = 'Helvetica')) %>% 
#     layout(xaxis = axis,yaxis = axis)
#   return(p)
# }
