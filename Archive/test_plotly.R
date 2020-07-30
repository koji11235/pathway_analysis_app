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
    g <- delete.edges(g, E(g)[wd[,3] < 0.2])
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

#network <- 
plot_ly(x = ~Xn, y = ~Yn,
        mode = "markers",
        #text = vs$label, textposition = "bottom center",
        color = vs$color, colors = "RdYlBu",
        size = vs$size, marker = list(size = vs$size),
        hoverinfo = "text") %>% 
add_text(text = vs$label, textposition = "bottom center",
         textfont = list(size = 12, family = 'Helvetica') )


font_size <- list(size = 12, family = 'Lato')
width <- 0.5
gray_axis <- '#dadada'
font_size <- list(size = 12, family = 'Lato')
width <- 0.5
#legend_name <- Hmisc::capitalize( quo_name(zenc) ) # WATCH OUT! it works only with the package with Hmisc
decimal <- ',.2f'
sep <- ','

p <- plot_ly(x = ~Xn, y = ~Yn, name = '', type = "scatter",  mode = 'markers', 
             hoverlabel = list(font = font_size), size = vs$size, color = vs$color, hoverinfo = "text+y", colors= "RdYlBu",
             marker = list(size = vs$size)) %>% 
  add_text(text = vs$label, textfont = font_size, textposition = "bottom center") %>% hide_colorbar()

  

  
network <- plot_ly()
for (i in seq_along(vs)){
  # print(i)
  fig <- fig %>% add_markers(
    x = Xn[i],
    y = Yn[i],
    showlegend = F
  )
  fig <- fig %>% add_annotations(
    x = Xn[i],
    y = Yn[i],
    xref = "x",
    yref = "y",
    text = vs$label[i],
    xanchor = 'right',
    showarrow = F
  )
}
  




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
  p,
  title = 'Karate Network',
  shapes = edge_shapes,
  xaxis = axis,
  yaxis = axis
)

fig
p %>%  layout(shapes = edge_shapes)


p %>%  layout(shapes = edge_shapes,
  xaxis = list(zeroline = F,
                              title = '',
                              linecolor = gray_axis,
                              titlefont = font_size,
                              tickfont  = font_size,
                              rangemode='tozero',
                              gridcolor = gray_axis,
                              gridwidth = width,
                              hoverformat = decimal,
                              mirror = "ticks",
                              tickmode = 'array',
                              tickcolor = gray_axis,
                              linewidth = width,
                              showgrid = F ),
                 yaxis = list(#title = y_name,
                              zerolinecolor = gray_axis,
                              linecolor = gray_axis,
                              mirror = "ticks",
                              hoverformat = '.2f',
                              linewidth = width,
                              tickcolor = gray_axis,
                              tickformat = '.2f',
                              titlefont = font_size,
                              tickfont  = font_size,
                              showgrid = FALSE) ) %>%
  config(displayModeBar = F)

#################
# https://stackoverflow.com/questions/58737891/set-text-size-within-marker-in-r-plotly-bubble-chart
library(plotly)
library(dplyr)

DF = data.frame( x = c( 3, 2, 2, 4, 6, 8 ), y = c( 3, 2, 3, 4, 6, 7 ), size = c( 20, 20, 20, 30, 40, 40 ), labels = letters[1:6] )

evo_bubble <- function(plot_data, x_var, y_var, z_var, t_var) {
  # browser()
  # Trasform data into dataframe and quos
  DF <- data.frame(plot_data)
  xenc <- enquo(x_var)
  yenc <- enquo(y_var)
  zenc <- enquo(z_var)
  tenc <- enquo(t_var)
  DF <- DF %>% mutate( bubble_size = !!zenc*50 ) # Modify the denominator if you want to change the dimension of the bubble
  
  
  # Set parameters for the plot 
  bubble_pal <- c("white", "#AECEE8" )
  gray_axis <- '#dadada'
  font_size <- list(size = 12, family = 'Lato')
  width <- 0.5
  legend_name <- Hmisc::capitalize( quo_name(zenc) ) # WATCH OUT! it works only with the package with Hmisc
  decimal <- ',.2f'
  sep <- ','
  #x_name <- capitalize(quo_name(xenc))
  y_name <- Hmisc::capitalize(quo_name(yenc))
  
  p <- plot_ly(DF, x = xenc, y = yenc, name = '', type = "scatter",  mode = 'markers', 
               hoverlabel = list(font = font_size), size = zenc, color = zenc, hoverinfo = "text+y", colors= bubble_pal,
               marker = list(size = DF$bubble_size, line = list(color = gray_axis))) %>% 
    add_text(text = tenc, textfont = font_size, textposition = "middle center") %>% hide_colorbar()
  print(p)
  p <- p %>%  layout(xaxis = list(zeroline = F,
                                  title = '',
                                  linecolor = gray_axis,
                                  titlefont = font_size,
                                  tickfont  = font_size,
                                  rangemode='tozero',
                                  gridcolor = gray_axis,
                                  gridwidth = width,
                                  hoverformat = decimal,
                                  mirror = "ticks",
                                  tickmode = 'array',
                                  tickcolor = gray_axis,
                                  linewidth = width,
                                  showgrid = F ),
                     yaxis = list(title = y_name,
                                  zerolinecolor = gray_axis,
                                  linecolor = gray_axis,
                                  mirror = "ticks",
                                  hoverformat = '.2f',
                                  linewidth = width,
                                  tickcolor = gray_axis,
                                  tickformat = '.2f',
                                  titlefont = font_size,
                                  tickfont  = font_size,
                                  showgrid = FALSE) ) %>%
    config(displayModeBar = F)
  
  #return(p)
}


evo_bubble( DF, x, y, size, labels )

