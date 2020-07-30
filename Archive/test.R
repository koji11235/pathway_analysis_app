library(plotly)
library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(MetamapsDB)

#library(DOSE)
#library(org.Hs.eg.db)

setwd("~/Project/20200719_Pathway_app")

sample_data<-read.csv("sample_gene_data.csv")

#サンプルデータの中のEntrez_Gene_IDを抽出
gene_IDs<-sample_data$Entrez_Gene_ID

#KEGGエンリッチメント解析を実行
kegg_enrichment_result <- enrichKEGG(gene=gene_IDs,pvalueCutoff=0.05, readable=T)

cnetplot(kegg_enrichment_result, )
reactome_enrichment_result <- enrichPathway(gene=gene_IDs,pvalueCutoff=0.05, readable=T)
cn <- cnetplot(reactome_enrichment_result, showCategory = 5)
ggplotly(cn)
layout <- ggraph::create_layout(cn)
cn+ theme_classic()
tidygraph::as.igraph(cn$layers)

ggraph(cn)
em <- emapplot(reactome_enrichment_result, showCategory = 5)
em$data
###########################################################
#geneListデータセットから発現変動遺伝子のGene IDを抽出
data(geneList)
gene <- names(geneList)[abs(geneList) > 2]
gene_IDs<-sample_data$Entrez_Gene_ID

#KEGGエンリッチメント解析を実行
kegg_enrichment_result <- enrichKEGG(gene=gene_IDs,pvalueCutoff=0.05)

#結果の表示
head(kegg_enrichment_result, n=10)

kegg_enrichment_bp　<- barplot(kegg_enrichment_result, drop=TRUE, showCategory=12)
kegg_enrichment_bp_plotly <- ggplotly(kegg_enrichment_bp)


ggplotly(emapplot(kegg_enrichment_result, showCategory=12))



##################################################################
#サンプルデータ読み込み
sample_data<-read.csv("sample_gene_data.csv")

#サンプルデータの中のEntrez_Gene_IDを抽出
gene_IDs<-sample_data$Entrez_Gene_ID

#KEGGエンリッチメント解析を実行
kegg_enrichment_result <- enrichKEGG(gene=gene_IDs,pvalueCutoff=0.05)

#結果の表示
head(kegg_enrichment_result, n=10)

kegg_enrichment_bp　<- barplot(kegg_enrichment_result, drop=TRUE, showCategory=12)

kegg_enrichment_bp
kegg_enrichment_bp_plotly <- ggplotly(kegg_enrichment_bp)





##################################################################
library(plotly)

p <- ggplot(fortify(forecast::gold), aes(x, y)) + geom_line()

fig <- ggplotly(p)

fig <- style(fig, line = list(color = 'gold'), hoverinfo = "y", traces = 1)

fig

###############################
library(MetamapsDB)
library(ggplot2)
library(plotly)
library(ITNr)

data("ELEnet16")

mygraph <- ELEnet16#igraph::graph_from_data_frame(ELM) #create a graph
mylayout <- igraph::layout_as_tree(ELEnet16, circular = T) #create a circular layout

mygraph$layout = mylayout #store layout as a graph attribute


#another gotcha is that ig2ggplot needs both vertex names and vertex labels. 
#as of now you just have vertex names. 
V(mygraph)$label = V(mygraph)$name #store label as a vertex attrbute


g <- MetamapsDB::ig2ggplot(mygraph, 
                           dfOnly = FALSE, 
                           labels = FALSE, 
                           metab = TRUE ) + 
  theme(legend.position = 'none')
g
ggplotly(g)


###############################
library(plotly)
library(igraph)
library(ITNr)

data("ELEnet16")

vs <- V(ELEnet16) #Get node list
es <- as.data.frame(get.edgelist(ELEnet16)) # Get edgelist
node.data<-get.data.frame(ELEnet16,what="vertices") # Get node attribute dataframe

Nv <- length(vs) #number of nodes
Ne <- length(es[1]$V1) #number of edges

#Coordinates for nodes
L <- layout.fruchterman.reingold(ELEnet16)
Xn <- L[,1]
Yn <- L[,2]

#Creates the nodes (plots the points)
network <- plot_ly(x = ~Xn, y = ~Yn, #Node points
                   mode = "markers", 
                   text = vs$name, 
                   hoverinfo = "text",
                   color =as.factor(node.data$region) )

#Create edges
edge_shapes <- list()
for(i in 1:Ne) {
  v0 <- es[i,]$V1
  v1 <- es[i,]$V2
  
  edge_shape = list(
    type = "line",
    line = list(color = "gray", width = 0.3),
    x0 = Xn[v0],
    y0 = Yn[v0],
    x1 = Xn[v1],
    y1 = Yn[v1]
  )
  
  edge_shapes[[i]] <- edge_shape
}

axis <- list(title = "", showgrid = FALSE, 
             showticklabels = FALSE, zeroline = FALSE)

p <- layout(
  network,
  title = 'Networks & Plotly',
  shapes = edge_shapes,
  xaxis = axis,
  yaxis = axis,
  showlegend=FALSE
)
p





