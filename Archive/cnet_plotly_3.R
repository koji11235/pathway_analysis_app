# different category color but same gene color


library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(DOSE)
library(reshape2)
library(igraph)
library(MetamapsDB)
library(RColorBrewer)


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

# create reactome_enrichment_result -----------------------------------
sample_data<-read.csv("sample_gene_data.csv")

#サンプルデータの中のEntrez_Gene_IDを抽出
gene_IDs<-sample_data$Entrez_Gene_ID

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



enrichment_result <- reactome_enrichment_result
##################################################
# arguments of cnetplot ----------------------------
showCategory = 2
foldChange   = NULL
layout = "kk"
colorEdge = FALSE
circular = FALSE
node_label = "all"
category_color_palette <- brewer.pal(showCategory, "Set2")
gene_font_size = 12
category_font_size = 20

node_label <- match.arg(node_label, c("category", "gene", "all", "none"))


# create igraph object----------------------------------------

# circularに対応できてない
# if (circular) {
#   layout <- "linear"
#   geom_edge <- geom_edge_arc
# } else {
#   geom_edge <- geom_edge_link
# }

geneSets <- extract_geneSets(enrichment_result, showCategory)

g <- list2graph(geneSets)

foldChange <- fc_readable(x, foldChange)

# set node size ----------------------------------------
size <- sapply(geneSets, length)
V(g)$size <- min(size)/2
n <- length(geneSets)
V(g)$size[1:n] <- size

# set node text font size ----------------------------------------
V(g)$font_size <- gene_font_size
V(g)$font_size[1:n] <- category_font_size


# if (colorEdge) {
#   E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
#   #edge_layer <- geom_edge(aes_(color = ~category), alpha=.8)
# } else {
#   #edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
# }


# set gene node color ----------------------------------------
if (!is.null(foldChange)) {
  fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
  V(g)$color <- NA
  V(g)$color[(n+1):length(V(g))] <- fc
  #palette <- fc_palette(fc)
  # p <- ggraph(g, layout=layout, circular = circular) +
  #   edge_layer +
  #   geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
  #   #scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494")
  #   scale_colour_gradient2(name = "fold change", low = "green", mid = "blue", high = "red")
} else {
  V(g)$color <- "#B3B3B3"
  #V(g)$color[1:n] <- "#E5C494"
  V(g)$color[1:n] <- category_color_palette[factor(names(geneSets))]
  # p <- ggraph(g, layout=layout, circular=circular) +
  #   edge_layer +
  #   geom_node_point(aes_(color=~I(color), size=~size))
}

g

# p <- p + scale_size(range=c(3, 10), breaks=unique(round(seq(min(size), max(size), length.out=4)))) +
#  theme_void()


# if (node_label == "category") {
#   p <- p + geom_node_text(aes_(label=~name), data = p$data[1:n,])
# } else if (node_label == "gene") {
#   p <- p + geom_node_text(aes_(label=~name), data = p$data[-c(1:n),], repel=TRUE)
# } else if (node_label == "all") {
#   p <- p + geom_node_text(aes_(label=~name), repel=TRUE)
# } 

L <- layout_with_kk(g)

vs <- V(g)
label <- names(vs)
es <- as.data.frame(get.edgelist(g))

Nv <- length(vs)
Ne <- length(es[1]$V1)

node2index <- c(1:Nv)
names(node2index) <- names(vs)

Xn <- L[,1]
Yn <- L[,2]

# geneSet2index <- c(1:n)
# names(geneSet2index) <- names(geneSets)

edge_df <- tibble(
  x = Xn[node2index[as.character(es$V1)]],
  xend = Xn[node2index[as.character(es$V2)]],
  y = Yn[node2index[as.character(es$V1)]],
  yend = Yn[node2index[as.character(es$V2)]],
  category = es$V1,
  category_color = category_color_palette[as.numeric(factor(es$V1))]
)


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

V(g)$color_test <- c(category_color_palette[factor(names(geneSets))], gene2color[names(V(g))[(n+1):length(V(g))]])


axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)

plot_ly(x = ~Xn, y = ~Yn,
        type = "scatter",
        mode = "markers",
        #color = ~vs$color,
        #size = ~vs$size,
        #sizes = c(100,800),
        text = label,
        marker = list(color = V(g)$color_test, 
                      size = vs$size,
                      line = list(width = 0)),
        hoverinfo = "text",
        showlegend=FALSE) %>% 
  add_segments_different_color(x=edge_df$x,  xend=edge_df$xend,
                               y=edge_df$y, yend=edge_df$yend, 
                               color = edge_df$category_color) %>% 
  add_markers(showlegend=T) %>% 
  add_text(text = label, textposition = "bottom center",
           textfont = list(size = V(g)$font_size, family = 'Helvetica')) %>% 
  layout(xaxis = axis,yaxis = axis)


