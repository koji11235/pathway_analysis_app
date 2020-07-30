library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(tidygraph)
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


# create reactome_enrichment_result -----------------------------------
sample_data<-read.csv("sample_gene_data.csv")

#サンプルデータの中のEntrez_Gene_IDを抽出
gene_IDs<-sample_data$Entrez_Gene_ID

#KEGGエンリッチメント解析を実行
kegg_enrichment_result <- enrichKEGG(gene=gene_IDs,pvalueCutoff=0.05, readable=T)
reactome_enrichment_result <- enrichPathway(gene=gene_IDs,pvalueCutoff=0.05, readable=T)
#cn <- cnetplot(reactome_enrichment_result, showCategory = 5)
em <- emapplot(reactome_enrichment_result, showCategory = 5)


# set arguments -----------------------------------
showCategory = 5
line_scale = 1
color="p.adjust"

# create graph -----------------------------------
x = reactome_enrichment_result

n <- update_n(x, showCategory)
geneSets <- geneInCategory(x)

y <- as.data.frame(x)
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
## emap_graph_buildの中身




font_size <- list(size = 12, family = 'Lato')
width <- 0.5
gray_axis <- '#dadada'
font_size <- list(size = 12, family = 'Lato')
width <- 0.5
#legend_name <- Hmisc::capitalize( quo_name(zenc) ) # WATCH OUT! it works only with the package with Hmisc
decimal <- ',.2f'
sep <- ','


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

es$V1
edge_df <- tibble(
  x = Xn[node2index[as.character(es$V1)]],
  xend = Xn[node2index[as.character(es$V2)]],
  y = Yn[node2index[as.character(es$V1)]],
  yend = Yn[node2index[as.character(es$V2)]],
  width = E(g)$width
)


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

axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)


plot_ly(x = ~Xn, y = ~Yn,
        type = "scatter",
        mode = "markers",
        color = ~vs$color,
        colors = "RdYlBu",
        size = ~vs$size,
        sizes = c(100,800),
        text = label,
        hoverinfo = "text",
        showlegend=FALSE)%>% 
  add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                              y=edge_df$y, yend=edge_df$yend, 
                              width = edge_df$width) %>% 
  add_markers(showlegend=FALSE) %>% 
  add_text(text = label, textposition = "bottom center",
           textfont = list(size = 12, family = 'Helvetica') )# %>% 
  layout(
    xaxis = axis,
    yaxis = axis
  )
