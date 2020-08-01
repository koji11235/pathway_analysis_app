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

cnet_plotly <- function(enrichment_result, showCategory, foldChange=NULL, is_kegg=FALSE){
  gene_font_size <- 8
  category_font_size <- 12
  gene_node_size <- 10
  node_size_scale <- c(gene_node_size=gene_node_size, category_size_min=20, category_size_max=40)
  
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
    category_color_palette <- viridis_pal(option = "C", direction = -1)(max(unique(gene2color_df$category_id))) 
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
  
  # set hover info text ----------------------------------------
  as.data.frame(get.edgelist(g)) %>% 
    group_by(V2) %>%
    summarise(categories = str_c(V1, collapse  = "; "),
              n_categories = n())
  
  hoverinfo_text <- tibble(node_name = names(V(g)),
                           node_label = V(g)$label) %>% 
    left_join(as.data.frame(get.edgelist(g)) %>% 
                group_by(V2) %>%
                summarise(categories = str_c(V1, collapse  = "\n"),
                          n_categories = n()),
              by = c("node_name" = "V2")
    ) %>% left_join(enrichment_result@result,
                    by = c("node_name"="Description")) %>% 
    mutate(hoverinfo = if_else(is.na(categories), 
                               paste0(node_label, "\n",
                                      "p.adjust: ", signif(p.adjust, 3), "\n",
                                      "enriched genes: ", Count),
                               paste0(node_label, "\n",
                                      "belongs to ", n_categories, if_else(n_categories == 1, "category", "categories"), "\n",
                                      "category:", "\n", 
                                      categories, "\n")
    )
    ) %>% pull(hoverinfo)
  
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
               #hoverinfo = "text",
               hovertemplate = hoverinfo_text,
               showlegend=FALSE) %>% 
    add_segments_different_color(x=edge_df$x,  xend=edge_df$xend,
                                 y=edge_df$y, yend=edge_df$yend, 
                                 color = edge_df$category_color,
                                 line_width = 2) %>% 
    add_markers(showlegend=T) %>% 
    add_text(text = V(g)$label, textposition = "center",
             textfont = list(size = V(g)$font_size, family = 'Helvetica')) %>% 
    layout(xaxis = axis,yaxis = axis)
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



cnet_plotly(enrichment_result, showCategory = 3, selected_rows = c(2), is_kegg = TRUE)

cnet_plotly(enrichment_result, showCategory = 3, selected_rows = c(1,2), is_kegg = TRUE)
cnet_plotly(enrichment_result, showCategory = 3, selected_rows = c(), is_kegg = TRUE)
