library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
#library(ggraph)
library(tidygraph)
library(DOSE)
library(reshape2)
library(igraph)
library(MetamapsDB)

#library(org.Hs.eg.db)

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

g


mygraph <- g #igraph::graph_from_data_frame(ELM) #create a graph
mylayout <- igraph::layout_as_tree(g, circular = T) #create a circular layout

mygraph$layout = mylayout #store layout as a graph attribute


#another gotcha is that ig2ggplot needs both vertex names and vertex labels. 
#as of now you just have vertex names. 
V(mygraph)$label = V(mygraph)$name #store label as a vertex attrbute

#https://plotly.com/r/network-graphs/

L <- layout.circle(g)

vs <- V(g)

es <- as.data.frame(get.edgelist(g))

Nv <- length(vs)
Ne <- length(es[1]$V1)

node2index <- c(1:Nv)
names(node2index) <- names(vs)

Xn <- L[,1]
Yn <- L[,2]

network <- plot_ly(x = ~Xn, y = ~Yn, 
                   mode = "text", 
                   text = vs$label,
                   hoverinfo = "text") %>% 
  add_markers()

edge_shapes <- list()
for(i in 1:Ne) {
  v0 <- es[i,]$V1
  v1 <- es[i,]$V2
  
  edge_shape = list(
    type = "line",
    line = list(color = "#030303", width = 0.3),
    x0 = Xn[node2index[[as.character(v0)]]],
    y0 = Yn[node2index[[as.character(v0)]]],
    x1 = Xn[node2index[[as.character(v1)]]],
    y1 = Yn[node2index[[as.character(v1)]]]
  )
  
  edge_shapes[[i]] <- edge_shape
}

axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)

fig <- layout(
  network,
  title = 'Karate Network',
  shapes = edge_shapes,
  xaxis = axis,
  yaxis = axis
)

fig
#####################################################
## 以下だとIsolated nodeが表示されない
plot(g)
plot(mygraph)
g2 <- MetamapsDB::ig2ggplot(mygraph, 
                           dfOnly = FALSE, 
                           labels = TRUE, 
                           metab = TRUE ) + 
  theme_graph(legend.position = 'none')

ggplotly(g2)

g$layout <- igraph::layout_with_drl(g)
V(g)$label = V(g)$name
MetamapsDB::ig2ggplot(g, 
                      dfOnly = FALSE, 
                      labels = TRUE, 
                      metab = TRUE )+ 
  theme(legend.position = 'none')





if(n == 1) {
  return(ggraph(g) + geom_node_point(color="red", size=5) + geom_node_text(aes_(label=~name)))
}
##} else {
p <- ggraph(g, layout=layout)
if (length(E(g)$width) > 0) {
  p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)), colour='darkgrey')
}
p + geom_node_point(aes_(color=~color, size=~size)) +
  geom_node_text(aes_(label=~name), repel=TRUE) + theme_void() +
  scale_color_continuous(low="red", high="blue", name = color, guide=guide_colorbar(reverse=TRUE)) +
  scale_size(range=c(3, 8) * pie_scale)











