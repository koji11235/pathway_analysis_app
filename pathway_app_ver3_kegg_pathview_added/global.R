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

library(pathview)
library(png)

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




# KEGG pathview functions ---------------------------------------

mypathview <- function (gene.data = NULL, cpd.data = NULL, pathway.id, species = "hsa", 
                        kegg.dir = ".", cpd.idtype = "kegg", gene.idtype = "entrez", 
                        gene.annotpkg = NULL, min.nnodes = 3, kegg.native = TRUE, 
                        map.null = TRUE, expand.node = FALSE, split.group = FALSE, 
                        map.symbol = TRUE, map.cpdname = TRUE, node.sum = "sum", 
                        discrete = list(gene = FALSE, cpd = FALSE), limit = list(gene = 1, 
                                                                                 cpd = 1), bins = list(gene = 10, cpd = 10), both.dirs = list(gene = T, 
                                                                                                                                              cpd = T), trans.fun = list(gene = NULL, cpd = NULL), 
                        low = list(gene = "green", cpd = "blue"), mid = list(gene = "gray", 
                                                                             cpd = "gray"), high = list(gene = "red", cpd = "yellow"), 
                        na.col = "transparent", ...) 
{
  dtypes = !is.null(gene.data) + (!is.null(cpd.data))
  cond0 = dtypes == 1 & is.numeric(limit) & length(limit) > 
    1
  if (cond0) {
    if (limit[1] != limit[2] & is.null(names(limit))) 
      limit = list(gene = limit[1:2], cpd = limit[1:2])
  }
  if (is.null(trans.fun)) 
    trans.fun = list(gene = NULL, cpd = NULL)
  arg.len2 = c("discrete", "limit", "bins", "both.dirs", "trans.fun", 
               "low", "mid", "high")
  for (arg in arg.len2) {
    obj1 = eval(as.name(arg))
    if (length(obj1) == 1) 
      obj1 = rep(obj1, 2)
    if (length(obj1) > 2) 
      obj1 = obj1[1:2]
    obj1 = as.list(obj1)
    ns = names(obj1)
    if (length(ns) == 0 | !all(c("gene", "cpd") %in% ns)) 
      names(obj1) = c("gene", "cpd")
    assign(arg, obj1)
  }
  if (is.character(gene.data)) {
    gd.names = gene.data
    gene.data = rep(1, length(gene.data))
    names(gene.data) = gd.names
    both.dirs$gene = FALSE
    ng = length(gene.data)
    nsamp.g = 1
  }
  else if (!is.null(gene.data)) {
    if (length(dim(gene.data)) == 2) {
      gd.names = rownames(gene.data)
      ng = nrow(gene.data)
      nsamp.g = 2
    }
    else if (is.numeric(gene.data) & is.null(dim(gene.data))) {
      gd.names = names(gene.data)
      ng = length(gene.data)
      nsamp.g = 1
    }
    else stop("wrong gene.data format!")
  }
  else if (is.null(cpd.data)) {
    stop("gene.data and cpd.data are both NULL!")
  }
  gene.idtype = toupper(gene.idtype)
  data(bods)
  if (species != "ko") {
    species.data = kegg.species.code(species, na.rm = T, 
                                     code.only = FALSE)
  }
  else {
    species.data = c(kegg.code = "ko", entrez.gnodes = "0", 
                     kegg.geneid = "K01488", ncbi.geneid = NA, ncbi.proteinid = NA, 
                     uniprot = NA)
    gene.idtype = "KEGG"
    msg.fmt = "Only KEGG ortholog gene ID is supported, make sure it looks like \"%s\"!"
    msg = sprintf(msg.fmt, species.data["kegg.geneid"])
    message("Note: ", msg)
  }
  if (length(dim(species.data)) == 2) {
    message("Note: ", "More than two valide species!")
    species.data = species.data[1, ]
  }
  species = species.data["kegg.code"]
  entrez.gnodes = species.data["entrez.gnodes"] == 1
  if (is.na(species.data["ncbi.geneid"])) {
    if (!is.na(species.data["kegg.geneid"])) {
      msg.fmt = "Mapping via KEGG gene ID (not Entrez) is supported for this species,\nit looks like \"%s\"!"
      msg = sprintf(msg.fmt, species.data["kegg.geneid"])
      message("Note: ", msg)
    }
    else {
      stop("This species is not annotated in KEGG!")
    }
  }
  if (is.null(gene.annotpkg)) 
    gene.annotpkg = bods[match(species, bods[, 3]), 1]
  if (length(grep("ENTREZ|KEGG|NCBIPROT|UNIPROT", gene.idtype)) < 
      1 & !is.null(gene.data)) {
    if (is.na(gene.annotpkg)) 
      stop("No proper gene annotation package available!")
    if (!gene.idtype %in% gene.idtype.bods[[species]]) 
      stop("Wrong input gene ID type!")
    gene.idmap = id2eg(gd.names, category = gene.idtype, 
                       pkg.name = gene.annotpkg, unique.map = F)
    gene.data = mol.sum(gene.data, gene.idmap)
    gene.idtype = "ENTREZ"
  }
  if (gene.idtype != "KEGG" & !entrez.gnodes & !is.null(gene.data)) {
    id.type = gene.idtype
    if (id.type == "ENTREZ") 
      id.type = "ENTREZID"
    kid.map = names(species.data)[-c(1:2)]
    kid.types = names(kid.map) = c("KEGG", "ENTREZID", "NCBIPROT", 
                                   "UNIPROT")
    kid.map2 = gsub("[.]", "-", kid.map)
    kid.map2["UNIPROT"] = "up"
    if (is.na(kid.map[id.type])) 
      stop("Wrong input gene ID type for the species!")
    message("Info: Getting gene ID data from KEGG...")
    gene.idmap = keggConv(kid.map2[id.type], species)
    message("Info: Done with data retrieval!")
    kegg.ids = gsub(paste(species, ":", sep = ""), "", names(gene.idmap))
    in.ids = gsub(paste0(kid.map2[id.type], ":"), "", gene.idmap)
    gene.idmap = cbind(in.ids, kegg.ids)
    gene.data = mol.sum(gene.data, gene.idmap)
    gene.idtype = "KEGG"
  }
  if (is.character(cpd.data)) {
    cpdd.names = cpd.data
    cpd.data = rep(1, length(cpd.data))
    names(cpd.data) = cpdd.names
    both.dirs$cpd = FALSE
    ncpd = length(cpd.data)
  }
  else if (!is.null(cpd.data)) {
    if (length(dim(cpd.data)) == 2) {
      cpdd.names = rownames(cpd.data)
      ncpd = nrow(cpd.data)
    }
    else if (is.numeric(cpd.data) & is.null(dim(cpd.data))) {
      cpdd.names = names(cpd.data)
      ncpd = length(cpd.data)
    }
    else stop("wrong cpd.data format!")
  }
  if (length(grep("kegg", cpd.idtype)) < 1 & !is.null(cpd.data)) {
    data(rn.list)
    cpd.types = c(names(rn.list), "name")
    cpd.types = tolower(cpd.types)
    cpd.types = cpd.types[-grep("kegg", cpd.types)]
    if (!tolower(cpd.idtype) %in% cpd.types) 
      stop("Wrong input cpd ID type!")
    cpd.idmap = cpd2kegg(cpdd.names, in.type = cpd.idtype)
    cpd.data = mol.sum(cpd.data, cpd.idmap)
  }
  warn.fmt = "Parsing %s file failed, please check the file!"
  if (length(grep(species, pathway.id)) > 0) {
    pathway.name = pathway.id
    pathway.id = gsub(species, "", pathway.id)
  }
  else pathway.name = paste(species, pathway.id, sep = "")
  kfiles = list.files(path = kegg.dir, pattern = "[.]xml|[.]png")
  npath = length(pathway.id)
  out.list = list()
  tfiles.xml = paste(pathway.name, "xml", sep = ".")
  tfiles.png = paste(pathway.name, "png", sep = ".")
  if (kegg.native) 
    ttype = c("xml", "png")
  else ttype = "xml"
  xml.file <- paste(kegg.dir, "/", tfiles.xml, sep = "")
  for (i in 1:npath) {
    if (kegg.native) 
      tfiles = c(tfiles.xml[i], tfiles.png[i])
    else tfiles = tfiles.xml[i]
    if (!all(tfiles %in% kfiles)) {
      dstatus = download.kegg(pathway.id = pathway.id[i], 
                              species = species, kegg.dir = kegg.dir, file.type = ttype)
      if (dstatus == "failed") {
        warn.fmt = "Failed to download KEGG xml/png files, %s skipped!"
        warn.msg = sprintf(warn.fmt, pathway.name[i])
        message("Warning: ", warn.msg)
        return(invisible(0))
      }
    }
    if (kegg.native) {
      node.data = try(node.info(xml.file[i]), silent = T)
      if (class(node.data) == "try-error") {
        warn.msg = sprintf(warn.fmt, xml.file[i])
        message("Warning: ", warn.msg)
        return(invisible(0))
      }
      node.type = c("gene", "enzyme", "compound", "ortholog")
      sel.idx = node.data$type %in% node.type
      nna.idx = !is.na(node.data$x + node.data$y + node.data$width + 
                         node.data$height)
      sel.idx = sel.idx & nna.idx
      if (sum(sel.idx) < min.nnodes) {
        warn.fmt = "Number of mappable nodes is below %d, %s skipped!"
        warn.msg = sprintf(warn.fmt, min.nnodes, pathway.name[i])
        message("Warning: ", warn.msg)
        return(invisible(0))
      }
      node.data = lapply(node.data, "[", sel.idx)
    }
    else {
      gR1 = try(parseKGML2Graph2(xml.file[i], genes = F, 
                                 expand = expand.node, split.group = split.group), 
                silent = T)
      node.data = try(node.info(gR1), silent = T)
      if (class(node.data) == "try-error") {
        warn.msg = sprintf(warn.fmt, xml.file[i])
        message("Warning: ", warn.msg)
        return(invisible(0))
      }
    }
    if (species == "ko") 
      gene.node.type = "ortholog"
    else gene.node.type = "gene"
    if ((!is.null(gene.data) | map.null) & sum(node.data$type == 
                                               gene.node.type) > 1) {
      plot.data.gene = node.map(gene.data, node.data, node.types = gene.node.type, 
                                node.sum = node.sum, entrez.gnodes = entrez.gnodes)
      kng = plot.data.gene$kegg.names
      kng.char = gsub("[0-9]", "", unlist(kng))
      if (any(kng.char > "")) 
        entrez.gnodes = FALSE
      if (map.symbol & species != "ko" & entrez.gnodes) {
        if (is.na(gene.annotpkg)) {
          warn.fmt = "No annotation package for the species %s, gene symbols not mapped!"
          warn.msg = sprintf(warn.fmt, species)
          message("Warning: ", warn.msg)
        }
        else {
          plot.data.gene$labels = NA # Try to fix this error: Error in $<-.data.frame: replacement has 97 rows, data has 103
          plot.data.gene$labels = eg2id(as.character(plot.data.gene$kegg.names), 
                                        category = "SYMBOL", pkg.name = gene.annotpkg)[, 
                                                                                       2]
          mapped.gnodes = rownames(plot.data.gene)
          node.data$labels[mapped.gnodes] = plot.data.gene$labels
        }
      }
      cols.ts.gene = node.color(plot.data.gene, limit$gene, 
                                bins$gene, both.dirs = both.dirs$gene, trans.fun = trans.fun$gene, 
                                discrete = discrete$gene, low = low$gene, mid = mid$gene, 
                                high = high$gene, na.col = na.col)
    }
    else plot.data.gene = cols.ts.gene = NULL
    if ((!is.null(cpd.data) | map.null) & sum(node.data$type == 
                                              "compound") > 1) {
      plot.data.cpd = node.map(cpd.data, node.data, node.types = "compound", 
                               node.sum = node.sum)
      if (map.cpdname & !kegg.native) {
        plot.data.cpd$labels = cpdkegg2name(plot.data.cpd$labels)[, 
                                                                  2]
        mapped.cnodes = rownames(plot.data.cpd)
        node.data$labels[mapped.cnodes] = plot.data.cpd$labels
      }
      cols.ts.cpd = node.color(plot.data.cpd, limit$cpd, 
                               bins$cpd, both.dirs = both.dirs$cpd, trans.fun = trans.fun$cpd, 
                               discrete = discrete$cpd, low = low$cpd, mid = mid$cpd, 
                               high = high$cpd, na.col = na.col)
    }
    else plot.data.cpd = cols.ts.cpd = NULL
    if (kegg.native) {
      pv.pars = my.keggview.native(plot.data.gene = plot.data.gene, 
                                   cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd, 
                                   cols.ts.cpd = cols.ts.cpd, node.data = node.data, 
                                   pathway.name = pathway.name[i], kegg.dir = kegg.dir, 
                                   limit = limit, bins = bins, both.dirs = both.dirs, 
                                   discrete = discrete, low = low, mid = mid, high = high, 
                                   na.col = na.col, ...)
    }
    else {
      pv.pars = keggview.graph(plot.data.gene = plot.data.gene, 
                               cols.ts.gene = cols.ts.gene, plot.data.cpd = plot.data.cpd, 
                               cols.ts.cpd = cols.ts.cpd, node.data = node.data, 
                               path.graph = gR1, pathway.name = pathway.name[i], 
                               map.cpdname = map.cpdname, split.group = split.group, 
                               limit = limit, bins = bins, both.dirs = both.dirs, 
                               discrete = discrete, low = low, mid = mid, high = high, 
                               na.col = na.col, ...)
    }
    plot.data.gene = cbind(plot.data.gene, cols.ts.gene)
    if (!is.null(plot.data.gene)) {
      cnames = colnames(plot.data.gene)[-(1:8)]
      nsamp = length(cnames)/2
      if (nsamp > 1) {
        cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp + 
                                                          1):(2 * nsamp)], "col", sep = ".")
      }
      else cnames[2] = "mol.col"
      colnames(plot.data.gene)[-(1:8)] = cnames
    }
    plot.data.cpd = cbind(plot.data.cpd, cols.ts.cpd)
    if (!is.null(plot.data.cpd)) {
      cnames = colnames(plot.data.cpd)[-(1:8)]
      nsamp = length(cnames)/2
      if (nsamp > 1) {
        cnames[(nsamp + 1):(2 * nsamp)] = paste(cnames[(nsamp + 
                                                          1):(2 * nsamp)], "col", sep = ".")
      }
      else cnames[2] = "mol.col"
      colnames(plot.data.cpd)[-(1:8)] = cnames
    }
    out.list[[i]] = list(plot.data.gene = plot.data.gene, 
                         plot.data.cpd = plot.data.cpd)
  }
  if (npath == 1) 
    out.list = out.list[[1]]
  else names(out.list) = pathway.name
  return(invisible(out.list))
}
my.keggview.native <- function (plot.data.gene = NULL, plot.data.cpd = NULL, cols.ts.gene = NULL, 
                                cols.ts.cpd = NULL, node.data, pathway.name, out.suffix = "pathview", 
                                kegg.dir = ".", multi.state = TRUE, match.data = TRUE, same.layer = TRUE, 
                                res = 400, cex = 0.25, discrete = list(gene = FALSE, cpd = FALSE), 
                                limit = list(gene = 1, cpd = 1), bins = list(gene = 10, cpd = 10), 
                                both.dirs = list(gene = T, cpd = T), low = list(gene = "green", 
                                                                                cpd = "blue"), mid = list(gene = "gray", cpd = "gray"), 
                                high = list(gene = "red", cpd = "yellow"), na.col = "transparent", 
                                new.signature = TRUE, plot.col.key = TRUE, key.align = "x", 
                                key.pos = "topright", ...) 
{
  img <- readPNG(paste(kegg.dir, "/", pathway.name, ".png", 
                       sep = ""))
  width <- ncol(img)
  height <- nrow(img)
  cols.ts.gene = cbind(cols.ts.gene)
  cols.ts.cpd = cbind(cols.ts.cpd)
  nc.gene = max(ncol(cols.ts.gene), 0)
  nc.cpd = max(ncol(cols.ts.cpd), 0)
  nplots = max(nc.gene, nc.cpd)
  pn.suffix = colnames(cols.ts.gene)
  if (length(pn.suffix) < nc.cpd) 
    pn.suffix = colnames(cols.ts.cpd)
  if (length(pn.suffix) < nplots) 
    pn.suffix = 1:nplots
  if (length(pn.suffix) == 1) {
    pn.suffix = out.suffix
  }
  else pn.suffix = paste(out.suffix, pn.suffix, sep = ".")
  na.col = colorpanel2(1, low = na.col, high = na.col)
  if ((match.data | !multi.state) & nc.gene != nc.cpd) {
    if (nc.gene > nc.cpd & !is.null(cols.ts.cpd)) {
      na.mat = matrix(na.col, ncol = nplots - nc.cpd, nrow = nrow(cols.ts.cpd))
      cols.ts.cpd = cbind(cols.ts.cpd, na.mat)
    }
    if (nc.gene < nc.cpd & !is.null(cols.ts.gene)) {
      na.mat = matrix(na.col, ncol = nplots - nc.gene, 
                      nrow = nrow(cols.ts.gene))
      cols.ts.gene = cbind(cols.ts.gene, na.mat)
    }
    nc.gene = nc.cpd = nplots
  }
  out.fmt = "Working in directory %s"
  wdir = getwd()
  out.msg = sprintf(out.fmt, wdir)
  message("Info: ", out.msg)
  out.fmt = "Writing image file %s"
  multi.state = multi.state & nplots > 1
  if (multi.state) {
    nplots = 1
    pn.suffix = paste(out.suffix, "multi", sep = ".")
    if (nc.gene > 0) 
      cols.gene.plot = cols.ts.gene
    if (nc.cpd > 0) 
      cols.cpd.plot = cols.ts.cpd
  }
  for (np in 1:nplots) {
    # img.file = paste(pathway.name, pn.suffix[np], "png", 
    #    sep = ".")
    img.file = paste(kegg.dir,"/",pathway.name, ".",pn.suffix[np], ".png", 
                     sep = "")
    out.msg = sprintf(out.fmt, img.file)
    message("Info: ", out.msg)
    png(img.file, width = width, height = height, res = res)
    op = par(mar = c(0, 0, 0, 0))
    plot(c(0, width), c(0, height), type = "n", xlab = "", 
         ylab = "", xaxs = "i", yaxs = "i")
    if (new.signature) 
      img[height - 4:25, 17:137, 1:3] = 1
    if (same.layer != T) 
      rasterImage(img, 0, 0, width, height, interpolate = F)
    if (!is.null(cols.ts.gene) & nc.gene >= np) {
      if (!multi.state) 
        cols.gene.plot = cols.ts.gene[, np]
      if (same.layer != T) {
        render.kegg.node(plot.data.gene, cols.gene.plot, 
                         img, same.layer = same.layer, type = "gene", 
                         cex = cex)
      }
      else {
        img = render.kegg.node(plot.data.gene, cols.gene.plot, 
                               img, same.layer = same.layer, type = "gene")
      }
    }
    if (!is.null(cols.ts.cpd) & nc.cpd >= np) {
      if (!multi.state) 
        cols.cpd.plot = cols.ts.cpd[, np]
      if (same.layer != T) {
        render.kegg.node(plot.data.cpd, cols.cpd.plot, 
                         img, same.layer = same.layer, type = "compound", 
                         cex = cex)
      }
      else {
        img = render.kegg.node(plot.data.cpd, cols.cpd.plot, 
                               img, same.layer = same.layer, type = "compound")
      }
    }
    if (same.layer == T) 
      rasterImage(img, 0, 0, width, height, interpolate = F)
    pv.pars = list()
    pv.pars$gsizes = c(width = width, height = height)
    pv.pars$nsizes = c(46, 17)
    pv.pars$op = op
    pv.pars$key.cex = 2 * 72/res
    pv.pars$key.lwd = 1.2 * 72/res
    pv.pars$sign.cex = cex
    off.sets = c(x = 0, y = 0)
    align = "n"
    ucol.gene = unique(as.vector(cols.ts.gene))
    na.col.gene = ucol.gene %in% c(na.col, NA)
    if (plot.col.key & !is.null(cols.ts.gene) & !all(na.col.gene)) {
      off.sets = col.key(limit = limit$gene, bins = bins$gene, 
                         both.dirs = both.dirs$gene, discrete = discrete$gene, 
                         graph.size = pv.pars$gsizes, node.size = pv.pars$nsizes, 
                         key.pos = key.pos, cex = pv.pars$key.cex, lwd = pv.pars$key.lwd, 
                         low = low$gene, mid = mid$gene, high = high$gene, 
                         align = "n")
      align = key.align
    }
    ucol.cpd = unique(as.vector(cols.ts.cpd))
    na.col.cpd = ucol.cpd %in% c(na.col, NA)
    if (plot.col.key & !is.null(cols.ts.cpd) & !all(na.col.cpd)) {
      off.sets = col.key(limit = limit$cpd, bins = bins$cpd, 
                         both.dirs = both.dirs$cpd, discrete = discrete$cpd, 
                         graph.size = pv.pars$gsizes, node.size = pv.pars$nsizes, 
                         key.pos = key.pos, off.sets = off.sets, cex = pv.pars$key.cex, 
                         lwd = pv.pars$key.lwd, low = low$cpd, mid = mid$cpd, 
                         high = high$cpd, align = align)
    }
    if (new.signature) 
      pathview.stamp(x = 17, y = 20, on.kegg = T, cex = pv.pars$sign.cex)
    par(pv.pars$op)
    dev.off()
  }
  return(invisible(pv.pars))
}

colorpanel2<-function (n, low, mid, high){
  if (missing(mid) || missing(high)) {
    low <- col2rgb(low)
    if (missing(high))
      high <- col2rgb(mid)
    else high <- col2rgb(high)
    red <- seq(low[1, 1], high[1, 1], length = n)/255
    green <- seq(low[3, 1], high[3, 1], length = n)/255
    blue <- seq(low[2, 1], high[2, 1], length = n)/255
  }
  else {
    isodd <- n%%2 == 1
    if (isodd) {
      n <- n + 1
    }
    low <- col2rgb(low)
    mid <- col2rgb(mid)
    high <- col2rgb(high)
    lower <- floor(n/2)
    upper <- n - lower
    red <- c(seq(low[1, 1], mid[1, 1], length = lower),
             seq(mid[1, 1], high[1, 1], length = upper))/255
    green <- c(seq(low[3, 1], mid[3, 1], length = lower),
               seq(mid[3, 1], high[3, 1], length = upper))/255
    blue <- c(seq(low[2, 1], mid[2, 1], length = lower),
              seq(mid[2, 1], high[2, 1], length = upper))/255
    if (isodd) {
      red <- red[-(lower + 1)]
      green <- green[-(lower + 1)]
      blue <- blue[-(lower + 1)]
    }
  }
  rgb(red, blue, green)
}

render.kegg.node <-function(plot.data, cols.ts, img, same.layer=TRUE, type=c("gene","compound")[1], text.col="black", cex=0.25){
  width=ncol(img)
  height=nrow(img)
  nn=nrow(plot.data)
  pwids=plot.data$width
  if(!all(pwids==max(pwids))){
    message("Info: ", "some node width is different from others, and hence adjusted!")
    wc=table(pwids)
    pwids=plot.data$width=as.numeric(names(wc)[which.max(wc)])
  }
  
  if(type=="gene"){
    if(same.layer!=T){
      rect.out=sliced.shapes(plot.data$x+0.5, height-plot.data$y, plot.data$width/2-0.5, plot.data$height/2-0.25,  cols=cols.ts, draw.border=F, shape="rectangle")
      text(plot.data$x+0.5, height-plot.data$y, labels = as.character(plot.data$labels),
           cex = cex, col = text.col)
      return(invisible(1))
    } else{
      img2=img
      pidx=cbind(ceiling(plot.data$x-plot.data$width/2)+1,
                 floor(plot.data$x+plot.data$width/2)+1,
                 ceiling(plot.data$y-plot.data$height/2)+1,
                 floor(plot.data$y+plot.data$height/2)+1)
      cols.ts=cbind(cols.ts)
      ns=ncol(cols.ts)
      brk.x= sapply(plot.data$width/2, function(wi) seq(-wi, wi, length = ns+1))
      for(k in 1:ns){
        col.rgb=col2rgb(cols.ts[,k])/255
        pxr=t(apply(pidx[,1:2], 1, function(x) x[1]:x[2]))-plot.data$x-1
        sel=pxr>=ceiling(brk.x[k,]) & pxr<=floor(brk.x[k+1,])
        for(i in 1:nn){
          sel.px=(pidx[i,1]:pidx[i,2])[sel[i,]]
          node.rgb=img[pidx[i,3]:pidx[i,4],sel.px, 1:3]
          node.rgb.sum=apply(node.rgb,c(1,2), sum)
          blk.ind=which(node.rgb.sum==0|node.rgb.sum==1,arr.ind=T)
          node.rgb=array(col.rgb[,i],dim(node.rgb)[3:1])
          node.rgb=aperm(node.rgb, 3:1)
          for(j in 1:3) node.rgb[cbind(blk.ind,j)]=0
          img2[pidx[i,3]:pidx[i,4],sel.px, 1:3]=node.rgb
        }
      }
      return(img2)
    }
  } else if(type=="compound"){
    if(same.layer!=T){
      nc.cols=ncol(cbind(cols.ts))
      if(nc.cols>2){#block the background circle
        na.cols=rep("#FFFFFF", nrow(plot.data))
        cir.out=sliced.shapes(plot.data$x, height-plot.data$y, plot.data$width[1], plot.data$width[1], cols=na.cols, draw.border=F, shape="ellipse", lwd=0.2)
      }
      cir.out=sliced.shapes(plot.data$x, height-plot.data$y, plot.data$width[1], plot.data$width[1], cols=cols.ts, shape="ellipse", blwd=0.2)
      return(invisible(1))
    } else{
      #    col.rgb=col2rgb(cols.ts)/255
      blk=c(0,0,0)
      img2=img
      w=ncol(img) #repeat
      h=nrow(img) #repeat
      cidx=rep(1:w, each=h)
      ridx=rep(1:h, w)
      pidx=lapply(1:nn, function(i){
        ii=which((cidx-plot.data$x[i])^2+(ridx-plot.data$y[i])^2<(plot.data$width[i])^2)
        imat=cbind(cbind(ridx, cidx)[rep(ii,each=3),],1:3)
        imat[,1:2]=imat[,1:2]+1
        ib=which(abs((cidx-plot.data$x[i])^2+(ridx-plot.data$y[i])^2-(plot.data$width[i])^2)<=8)
        ibmat=cbind(cbind(ridx, cidx)[rep(ib,each=3),],1:3)
        ibmat[,1:2]=ibmat[,1:2]+1
        return(list(fill=imat,border=ibmat))
      })
      
      cols.ts=cbind(cols.ts)
      ns=ncol(cols.ts)
      brk.x= sapply(plot.data$width, function(wi) seq(-wi, wi, length = ns+1))
      for(i in 1:nn){
        pxr=pidx[[i]]$fill[,2]-1-plot.data$x[i]
        col.rgb=col2rgb(cols.ts[i,])/255
        for(k in 1:ns){
          sel=pxr>=brk.x[k,i] & pxr<=brk.x[k+1,i]
          img2[pidx[[i]]$fill[sel,]]=col.rgb[,k]
        }
        img2[pidx[[i]]$border]=blk
      }
      return(img2)
    }
  } else stop("unrecognized node type!")
}

pathview.stamp <-function(x=NULL, y=NULL, position="bottomright", graph.sizes, on.kegg=TRUE, cex=1){
  if(on.kegg)    labels ="Data on KEGG graph\nRendered by Pathview"
  else labels="-Data with KEGG pathway-\n-Rendered  by  Pathview-"
  if(is.null(x)| is.null(y)){
    x=graph.sizes[1]*.80
    y=graph.sizes[2]/40
    if(length(grep('left',position))==1)  x=graph.sizes[1]/40
    if(length(grep('top', position))==1)  y=graph.sizes[2]-y
  }
  text(x=x, y=y, labels=labels, adj=0, cex = cex, font=2)
}
