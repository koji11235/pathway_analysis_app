
library(org.Hs.eg.db)

# create reactome_enrichment_result -----------------------------------
sample_data<-read.csv("sample_gene_data.csv")

#サンプルデータの中のEntrez_Gene_IDを抽出
gene_IDs<-sample_data$Entrez_Gene_ID

#KEGGエンリッチメント解析を実行
GO_enrichment_result <- enrichGO(gene=gene_IDs, pvalueCutoff=0.05, readable = TRUE,
                                 OrgDb = org.Hs.eg.db, 
                                 ont = "MF")


##################################################
showCategory = 5
enrichment_result <- GO_enrichment_result
barplot_plotly(enrichment_result, showCategory = 8)
dotplot_plotly(enrichment_result, showCategory = 8)
emap_plotly(enrichment_result, showCategory = 8)
cnet_plotly(enrichment_result, showCategory = 3)


##################################################
# 4 functions
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
