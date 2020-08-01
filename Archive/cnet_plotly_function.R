# different category color but same gene color


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



# funnction for cnetplot -------------------------
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

add_segments_different_color <- function(p, x, xend, y, yend, color){
  tmp <- p
  for (i in 1:length(x)){
    tmp <- tmp %>% 
      add_segments(data=tibble(x_column = x[i],
                               xend_column = xend[i],
                               y_column = y[i],
                               yend_column = yend[i]), 
                   x=~x_column,  xend=~xend_column,
                   y=~y_column, yend=~yend_column, 
                   showlegend=FALSE,
                   color = I(color[i]),
                   inherit = FALSE,
                   line = list(width = 1)) 
  }
  return(tmp)
}


calculate_gene_color <- function(categories, categories2color){
  result_color <- "#FFFFFF"
  minimum_rgb_value = c(red=54, green=54, blue=54)
  hsvs = lapply(categories2color, function(color){rgb2hsv(col2rgb(color))}) 
  v_list =c()
  for (hsv in hsvs){
    v_list <- c(v_list, hsv["v",1])
  }
  mean_v = mean(v_list)
  
  for (category in categories){
    result_color_rgb <- mapply(min, 
                               col2rgb(categories2color[[category]])[,1], 
                               col2rgb(result_color)[,1])
    result_color <- rgb(red = result_color_rgb["red"], 
                        green = result_color_rgb["green"], 
                        blue = result_color_rgb["blue"], 
                        maxColorValue = 255)
  }
  result_color_rgb <- mapply(max, 
                             minimum_rgb_value, 
                             col2rgb(result_color)[,1])
  result_color <- rgb(red = result_color_rgb["red"], 
                      green = result_color_rgb["green"], 
                      blue = result_color_rgb["blue"], 
                      maxColorValue = 255)
  
  if (length(categories) > 1){
    result_color_hsv <- col2rgb(result_color) %>% rgb2hsv() 
    result_color_hsv["v",1] <- mean_v
    result_color <- hsv(result_color_hsv["h",1], result_color_hsv["s",1],result_color_hsv["v",1])
  }
  return(result_color)
}

# reorder_palette <- function(palette){
#   reordered <- c()
#   for (i in seq_along(palette)){
#     #print(i)
#     if (i%%2 == 1){
#       # print(i)
#       reordered[i] <- palette[(i - 1) / 2 + 1]
#     }else{
#       reordered[i] <- palette[length(palette) - (i / 2 - 1)]
#     }
#   }
#   return(reordered)
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

# create reactome_enrichment_result -----------------------------------
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
##################################################
# arguments of cnetplot ----------------------------
showCategory = 5
foldChange   = NULL
# layout = "kk"
# colorEdge = FALSE
# circular = FALSE
# node_label = "all"


cnet_plotly <- function(enrichment_result, showCategory, foldChange=NULL, is_kegg=FALSE){
  gene_font_size = 8
  category_font_size = 16
  color_palette = "Set2"
  
  actual_showCategory <- update_n(enrichment_result, showCategory)
  category_color_palette <- brewer.pal(actual_showCategory, color_palette)
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
    V(g)$label <- str_wrap(node_names, width = 40)
  }
  # set node size ----------------------------------------
  #(gene_set_sizes - min(gene_set_sizes)) / (max(gene_set_sizes) - min(gene_set_sizes)) * (node_size_scale[["category_size_max"]] - node_size_scale[["category_size_min"]]) + node_size_scale[["category_size_min"]]
  
  gene_node_size <- 15
  node_size_scale <- c(gene_node_size=gene_node_size, category_size_min=45, category_size_max=90)
  gene_set_sizes <- sapply(geneSets, length)
  gene_set_sizes_scaled <- (gene_set_sizes - min(gene_set_sizes)) / (max(gene_set_sizes) - min(gene_set_sizes)) * (node_size_scale[["category_size_max"]] - node_size_scale[["category_size_min"]]) + node_size_scale[["category_size_min"]]
  V(g)$size <- gene_node_size
  V(g)$size[1:actual_showCategory] <- gene_set_sizes_scaled
  
  # set node text font size ----------------------------------------
  V(g)$font_size <- gene_font_size
  V(g)$font_size[1:actual_showCategory] <- category_font_size
  
  # set gene node color ----------------------------------------
  if (!is.null(foldChange)) {
    fc <- foldChange[V(g)$name[(actual_showCategory+1):length(V(g))]]
    V(g)$color <- NA
    V(g)$color[(actual_showCategory+1):length(V(g))] <- fc
  } else {
    #V(g)$color <- "#B3B3B3"
    #V(g)$color[1:n] <- category_color_palette[factor(names(geneSets))]
    categories2color <- category_color_palette[as.numeric(factor(names(geneSets)))]
    names(categories2color) <- names(geneSets)
    
    gene2color_df <- as.data.frame(get.edgelist(g)) %>% 
      group_by(V2) %>% 
      summarise(categories = str_c(V1, collapse  = ";"),
                n = n(),
                gene_color = calculate_gene_color(V1, categories2color)) %>% 
      arrange(desc(n)) 
    gene2color <- gene2color_df %>% pull(gene_color)
    names(gene2color) <- gene2color_df %>% pull(V2)
    
    V(g)$color <- c(category_color_palette[factor(names(geneSets))], gene2color[names(V(g))[(actual_showCategory+1):length(V(g))]])
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
    category_color = category_color_palette[as.numeric(factor(edge_list$V1))]
  )
  
  
  # plot cnet plot by plotly ----------------------------------------
  axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
  
  p <- plot_ly(x = ~X_node, y = ~Y_node,
          type = "scatter",
          mode = "markers",
          text = V(g)$label,
          marker = list(color = V(g)$color, 
                        size = V(g)$size,
                        line = list(width = 0)),
          hoverinfo = "text",
          showlegend=FALSE) %>% 
    add_segments_different_color(x=edge_df$x,  xend=edge_df$xend,
                                 y=edge_df$y, yend=edge_df$yend, 
                                 color = edge_df$category_color) %>% 
    add_markers(showlegend=T) %>% 
    add_text(text = V(g)$label, textposition = "center",
             textfont = list(size = V(g)$font_size, family = 'Helvetica')) %>% 
    layout(xaxis = axis,yaxis = axis)
  return(p)
} 
cnet_plotly(reactome_enrichment_result, showCategory = 5, foldChange   = NULL)

cnet_plotly(kegg_enrichment_result, showCategory = 5, foldChange   = NULL, is_kegg = TRUE)

cnet_plotly <- function(enrichment_result, showCategory, foldChange=NULL, is_kegg=FALSE){
  gene_font_size = 8
  category_font_size = 16
  color_palette = "Set2"
  
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
    V(g)$label <- str_wrap(node_names, width = 40)
  }
  # set node size ----------------------------------------
  #(gene_set_sizes - min(gene_set_sizes)) / (max(gene_set_sizes) - min(gene_set_sizes)) * (node_size_scale[["category_size_max"]] - node_size_scale[["category_size_min"]]) + node_size_scale[["category_size_min"]]
  
  gene_node_size <- 15
  node_size_scale <- c(gene_node_size=gene_node_size, category_size_min=45, category_size_max=90)
  gene_set_sizes <- sapply(geneSets, length)
  gene_set_sizes_scaled <- (gene_set_sizes - min(gene_set_sizes)) / (max(gene_set_sizes) - min(gene_set_sizes)) * (node_size_scale[["category_size_max"]] - node_size_scale[["category_size_min"]]) + node_size_scale[["category_size_min"]]
  V(g)$size <- gene_node_size
  V(g)$size[1:actual_showCategory] <- gene_set_sizes_scaled
  
  # set node text font size ----------------------------------------
  V(g)$font_size <- gene_font_size
  V(g)$font_size[1:actual_showCategory] <- category_font_size
  

  # set gene node color ----------------------------------------
  if (!is.null(foldChange)) {
    fc <- foldChange[V(g)$name[(actual_showCategory+1):length(V(g))]]
    V(g)$color <- NA
    V(g)$color[(actual_showCategory+1):length(V(g))] <- fc
  } else {
    #V(g)$color <- "#B3B3B3"
    #V(g)$color[1:n] <- category_color_palette[factor(names(geneSets))]
    
    
    # categories2color <- category_color_palette[as.numeric(factor(names(geneSets)))]
    # names(categories2color) <- names(geneSets)
    # 
    # gene2color_df <- as.data.frame(get.edgelist(g)) %>% 
    #   group_by(V2) %>% 
    #   summarise(categories = str_c(V1, collapse  = ";"),
    #             n = n(),
    #             gene_color = calculate_gene_color(V1, categories2color)) %>% 
    #   arrange(desc(n)) 
    # gene2color <- gene2color_df %>% pull(gene_color)
    # names(gene2color) <- gene2color_df %>% pull(V2)
    # 
    # V(g)$color <- c(category_color_palette[factor(names(geneSets))], gene2color[names(V(g))[(actual_showCategory+1):length(V(g))]])
    gene2color_df <- as.data.frame(get.edgelist(g)) %>% 
        group_by(V2) %>%
        summarise(categories = str_c(V1, collapse  = ";"),
                  n_categories = n())
    categories_with_combination <- gene2color_df %>% 
      pull(categories) %>% unique()
    categories_combination_only <- categories_with_combination[!(categories_with_combination %in% names(geneSets))]
    
    gene2color_df <- gene2color_df %>% 
      mutate(category_id = gene2color_df %>% 
               pull(categories) %>% 
               factor(levels = c(levels(factor(names(geneSets))), categories_combination_only)) %>% 
               as.numeric()
      )
    #gene2color_df %>% pull(category_id) %>% unique() %>% length()
    category_color_palette <- viridis_pal(option = "C", direction = -1)(length(categories_with_combination)) 
    category_color_palette <- reorder_palette(palette = category_color_palette, showCategory = actual_showCategory)
      
    gene2color <- category_color_palette[gene2color_df %>% pull(category_id)]
    names(gene2color) <- gene2color_df %>% pull(V2)
    V(g)$color <- c(category_color_palette[factor(names(geneSets))], gene2color[names(V(g))[(actual_showCategory+1):length(V(g))]])
    
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
    category_color = category_color_palette[as.numeric(factor(edge_list$V1))]
  )
  

  # plot cnet plot by plotly ----------------------------------------
  axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
  
  p <- plot_ly(x = ~X_node, y = ~Y_node,
               type = "scatter",
               mode = "markers",
               text = V(g)$label,
               #color = ~-log10(V(g)$color),
               #colors = viridis_pal(option = "C")(3), 
               marker = list(color = V(g)$color, 
                             size = V(g)$size,
                             line = list(width = 0)),
               hoverinfo = "text",
               showlegend=FALSE) %>% 
    add_segments_different_color(x=edge_df$x,  xend=edge_df$xend,
                                 y=edge_df$y, yend=edge_df$yend, 
                                 color = edge_df$category_color) %>% 
    add_markers(showlegend=T) %>% 
    add_text(text = V(g)$label, textposition = "center",
             textfont = list(size = V(g)$font_size, family = 'Helvetica')) %>% 
    layout(xaxis = axis,yaxis = axis)
  return(p)
}
cnet_plotly(kegg_enrichment_result, showCategory = 3, foldChange   = NULL, is_kegg = TRUE)
cnet_plotly(kegg_enrichment_result, showCategory = 4, foldChange   = NULL, is_kegg = TRUE)
cnet_plotly(kegg_enrichment_result, showCategory = 5, foldChange   = NULL, is_kegg = TRUE)

category_color_palette
palette <- viridis_pal(option = "C")(10)

reorder_palette <- function(palette){
  reordered <- c()
  palette <- rev(palette)
  for (i in seq_along(palette)){
    #print(i)
    if (i%%2 == 1){
      # print(i)
      reordered[i] <- palette[(i - 1) / 2 + 1]
    }else{
      reordered[i] <- palette[length(palette) - (i / 2 - 1)]
    }
  }
  return(reordered)
}

palette <- viridis_pal(option = "C")(14)
reorder_palette <- function(palette, showCategory){
  reordered <- c()
  separated_palettes <- lapply(seq_len(showCategory), function(i){NULL})
  
  block_size <- ceiling(length(palette) / showCategory)
  for (i in seq_len(showCategory)){
    separated_palettes[[i]] <- palette[((i - 1) * block_size + 1) : (i * block_size)]
  }
  
  for (i in seq_len(block_size)){
    for (separated_palette in separated_palettes){
      #print(paste(separated_palette, i))
      reordered <- c(reordered, separated_palette[[i]])
    }
  }
  reordered <- reordered[!is.na(reordered)]
  return(reordered)
}

reordered <- reorder_palette(palette, 4)
  



reorder_palette <- function(palette, showCategory){
  reordered <- c()
  #separated_palettes <- list(c()*showCategory)
  separated_palettes <- lapply(seq_len(showCategory), function(i){NULL})

  block_size <- ceiling(length(palette) / showCategory)
  for (i in seq_len(showCategory)){
    print(i)
    separated_palettes[[i]] <- palette[((i - 1) * block_size + 1) : (i * block_size)]
  }
  
  for (i in seq_len(block_size)){
    for (separated_palette in separated_palettes){
      reordered <- c(reordered, separated_palette[[i]])
    }
  }
  reordered <- reordered[seq_along(palette)]
  
  
  for (i in seq_along(palette)){
    #i %% showCategory 
    print(i)
    separated_palettes[[((i-1) %% showCategory + 1)]] <- c(separated_palettes[[((i-1) %% showCategory + 1)]] , palette[[i]])
  }
  
  
  
  for (i in seq_along(palette)){
    #print(i)
    if (i%%2 == 1){
      reordered[i] <- palette[(i - 1) / 2 + 1]
    }else{
      reordered[i] <- palette[length(palette) - (i / 2 - 1)]
    }
  }
  return(reordered)
}