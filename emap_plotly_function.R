library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(DOSE)
library(reshape2)
library(igraph)
library(MetamapsDB)


setwd("~/Project/20200719_Pathway_app")

# define functions -----------------------------------
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

add_segments_variable_width <- function(p, x, xend, y, yend, width){
  tmp <- p
  for (i in 1:length(x)){
    tmp <- tmp %>% 
      add_segments(data=tibble(x_column = x[i],
                               xend_column = xend[i],
                               y_column = y[i],
                               yend_column = yend[i],
                               width_column = width[i]), 
                   x=~x_column,  xend=~xend_column,
                   y=~y_column, yend=~yend_column, 
                   showlegend=FALSE,
                   color = I("#AAAAAA"),
                   inherit = FALSE,
                   line = list(width = width[i])) 
  }
  return(tmp)
}


# create reactome_enrichment_result -----------------------------------
sample_data<-read.csv("sample_gene_data.csv")

#サンプルデータの中のEntrez_Gene_IDを抽出
gene_IDs<-sample_data$Entrez_Gene_ID

#KEGGエンリッチメント解析を実行
kegg_enrichment_result <- enrichKEGG(gene=gene_IDs,pvalueCutoff=0.05)
reactome_enrichment_result <- enrichPathway(gene=gene_IDs,pvalueCutoff=0.05, readable=T)
#cn <- cnetplot(reactome_enrichment_result, showCategory = 5)
em <- emapplot(reactome_enrichment_result, showCategory = 5)
emapplot(kegg_enrichment_result, showCategory = 5)

##################################################
showCategory = 5
line_scale = 1
color="p.adjust"

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
  
  marker_size_scale = c(min = 15, max = 50)
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
emap_plotly(kegg_enrichment_result, showCategory = 12, line_scale = 1)
emap_plotly(reactome_enrichment_result, showCategory = 12, line_scale = 1)

##################################################
# previous version
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
  
  marker_size_scale = c(min = 10, max = 50)
  marker_size = (V(g)$size - min(V(g)$size)) / (max(V(g)$size) - min(V(g)$size)) * (marker_size_scale["max"] - marker_size_scale["min"]) + marker_size_scale["min"]
  # max minに合わせるのでなくてratioは保ってmaxを40に合わせる → この手法だと小さいものがあまりに小さくなってしまう
  # marker_size = (V(g)$size / max(V(g)$size)) * marker_size_scale["max"]
  
  #textpositions = c(rep("bottom", 2), rep("top", 3))
  colorbar_tickvals <- seq(floor(min(-log10(V(g)$color))), ceiling(max(-log10(V(g)$color))),  length=6)
  colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
  
  p <- plot_ly(x = ~X_node, y = ~Y_node,
               type = "scatter",
               mode = "markers",
               color = ~-log10(V(g)$color),
               colors = viridis_pal(option = "C", direction = -1)(3), #"RdYlBu",
               #size = ~V(g)$size,
               #sizes = c(100,800),
               #text = str_wrap(label, 30),
               #hoverinfo = "text",
               marker = list(#color = V(g)$color, 
                 size = marker_size, #V(g)$size, 
                 alpha = 1,
                 line = list(width = 0)#,
                 # color = -log10(V(g)$color),
                 # colorscale = viridis_pal(option = "C", direction = -1)(3),
                 # reversescale =T
                 ),
               # hovertemplate = paste('Pathway: %{text}<br>',
               #                       'p.adjust: %{color:.2f}<br>',
               #                       'Gene Set Size: %{~vs$size}'),
               hovertemplate = paste(label, '\n',
                                     'p.adjust:', signif(V(g)$color, 3), '\n',
                                     'num genes:', V(g)$size),
               showlegend=FALSE)%>% 
    add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                                y=edge_df$y, yend=edge_df$yend, 
                                width = edge_df$width) %>% 
    add_markers(showlegend=FALSE) %>% 
    add_text(text = str_wrap(label, 30), textposition = "bottom", #textpositions,#"auto",
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
             exponentformat="e"#,reversescale=TRUE
             )
  return(p)
}
emap_plotly(kegg_enrichment_result, showCategory = 12, line_scale = 1)

##################################################
# previous version
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
  
  #gray_axis <- '#dadada'
  #font_size <- list(size = 12, family = 'Lato')
  #width <- 0.5
  #legend_name <- Hmisc::capitalize( quo_name(zenc) ) # WATCH OUT! it works only with the package with Hmisc
  #decimal <- ',.2f'
  #sep <- ','
  
  
  L <- layout.circle(g)
  
  vs <- V(g)
  label <- names(vs)
  es <- as.data.frame(get.edgelist(g))
  
  Nv <- length(vs)
  Ne <- length(es[1]$V1)
  
  node2index <- c(1:Nv)
  names(node2index) <- names(vs)
  
  Xn <- L[,1]
  Yn <- L[,2]
  edge_df <- tibble(
    x = Xn[node2index[as.character(es$V1)]],
    xend = Xn[node2index[as.character(es$V2)]],
    y = Yn[node2index[as.character(es$V1)]],
    yend = Yn[node2index[as.character(es$V2)]],
    width = E(g)$width
  )
  
  axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
  
  
  p <- plot_ly(x = ~Xn, y = ~Yn,
          type = "scatter",
          mode = "markers",
          color = ~vs$color,
          colors = viridis_pal(option = "C")(3), #"RdYlBu",
          size = ~vs$size,
          sizes = c(100,800),
          text = label,
          #hoverinfo = "text",
          hovertemplate = paste('<i>Price</i>: %{y:.2f}',
                                '<br><b>X</b>: %{marker.size}<br>',
                                '<b>%{text}</b>'),
          showlegend=FALSE)%>% 
    add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                                y=edge_df$y, yend=edge_df$yend, 
                                width = edge_df$width) %>% 
    add_markers(showlegend=FALSE) %>% 
    add_text(text = label, textposition = "center",
             textfont = list(size = 12, family = 'Helvetica') ) %>% 
  layout(
    xaxis = axis,
    yaxis = axis
  )
  
  return(p)
}

kegg_emap <- emap_plotly(reactome_enrichment_result, showCategory, line_scale, color)
emap_plotly(kegg_enrichment_result, showCategory = 5, line_scale = 1)
emapplot(kegg_enrichment_result, showCategory = 5)

# colorbarの検討
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
  
  marker_size_scale = c(min = 10, max = 40)
  #marker_size = (V(g)$size - min(V(g)$size)) / (max(V(g)$size) - min(V(g)$size)) * (marker_size_scale["max"] - marker_size_scale["min"]) + marker_size_scale["min"]
  # max minに合わせるのでなくてratioは保ってmaxを40に合わせる
  marker_size = (V(g)$size / max(V(g)$size)) * marker_size_scale["max"]
  
  #textpositions = c(rep("bottom", 2), rep("top", 3))
  colorbar_tickvals <- seq(floor(min(-log10(V(g)$color))), ceiling(max(-log10(V(g)$color))),  length=6)
  colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
  
  p <- plot_ly(x = ~X_node, y = ~Y_node,
               type = "scatter",
               mode = "markers",
               color = ~-log10(V(g)$color),
               colors = viridis_pal(option = "C")(3), #"RdYlBu",
               #size = ~V(g)$size,
               #sizes = c(100,800),
               #text = str_wrap(label, 30),
               #hoverinfo = "text",
               marker = list(#color = V(g)$color, 
                             size = marker_size, #V(g)$size, 
                             alpha = 1,
                             line = list(width = 0)),
               # hovertemplate = paste('Pathway: %{text}<br>',
               #                       'p.adjust: %{color:.2f}<br>',
               #                       'Gene Set Size: %{~vs$size}'),
                hovertemplate = paste(label, '\n',
                                     'p.adjust:', signif(V(g)$color, 3), '\n',
                                     'num genes:', V(g)$size),
               showlegend=FALSE)%>% 
    add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                                y=edge_df$y, yend=edge_df$yend, 
                                width = edge_df$width) %>% 
    add_markers(showlegend=FALSE) %>% 
    add_text(text = str_wrap(label, 30), textposition = "bottom", #textpositions,#"auto",
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
    ) %>% #colorbar(title = "",exponentformat="e")#tickformat=".1e")
    colorbar(title = "", tickmode = "array", tickvals = colorbar_tickvals, ticktext = colorbar_ticktext, exponentformat="e")
  return(p)
}
emap_plotly(kegg_enrichment_result, showCategory = 12, line_scale = 1)

# https://plotly.com/r/reference/#scattergl-marker-colorbar-exponentformat
-log10(V(g)$color)
floor(min(-log10(V(g)$color)))
ceiling(max(-log10(V(g)$color)))
colorbar_tickvals <- seq(floor(min(-log10(V(g)$color))), ceiling(max(-log10(V(g)$color))),  length=6)
colorbar_ticktext <- 10^(-colorbar_tickvals)
colorbar_tickvals <- c(floor(min(-log10(V(g)$color))) - 0.5, seq(floor(min(-log10(V(g)$color))), ceiling(max(-log10(V(g)$color))),  length=5), ceiling(max(-log10(V(g)$color))) + 0.5)
colorbar_ticktext <- signif(10^(-colorbar_tickvals), 3)
##############
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
  
  marker_size_scale = c(min = 10, max = 40)
  #marker_size = (V(g)$size - min(V(g)$size)) / (max(V(g)$size) - min(V(g)$size)) * (marker_size_scale["max"] - marker_size_scale["min"]) + marker_size_scale["min"]
  # max minに合わせるのでなくてratioは保ってmaxを40に合わせる
  marker_size = (V(g)$size / max(V(g)$size)) * marker_size_scale["max"]
  
  #textpositions = c(rep("bottom", 2), rep("top", 3))
  colorbar_tickvals <- seq(floor(min(-log10(V(g)$color))), ceiling(max(-log10(V(g)$color))),  length=6)
  colorbar_ticktext <- signif(10^(-colorbar_tickvals), 2)
  
  p <- plot_ly(x = ~X_node, y = ~Y_node,
               type = "scatter",
               mode = "markers",
               color = ~-log10(V(g)$color),
               colors = viridis_pal(option = "C", direction = 1)(3), #"RdYlBu",
               #size = ~V(g)$size,
               #sizes = c(100,800),
               #text = str_wrap(label, 30),
               #hoverinfo = "text",
               marker = list(#color = V(g)$color, 
                 size = marker_size, #V(g)$size, 
                 alpha = 1,
                 line = list(width = 0)#,
                 #color = ~-log10(V(g)$color),
                 #colorscale = viridis::plasma(99),#viridis_pal(option = "C", direction = -1)(3),
                 #reversescale =T
               ),
               # hovertemplate = paste('Pathway: %{text}<br>',
               #                       'p.adjust: %{color:.2f}<br>',
               #                       'Gene Set Size: %{~vs$size}'),
               hovertemplate = paste(label, '\n',
                                     'p.adjust:', signif(V(g)$color, 3), '\n',
                                     'num genes:', V(g)$size),
               showlegend=FALSE)%>% 
    add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                                y=edge_df$y, yend=edge_df$yend, 
                                width = edge_df$width) %>% 
    add_markers(showlegend=FALSE) %>% 
    add_text(text = str_wrap(label, 30), textposition = "bottom", #textpositions,#"auto",
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
             exponentformat="e"#,reversescale=TRUE
    )
  return(p)
}

emap_plotly(kegg_enrichment_result, showCategory = 12, line_scale = 1)

