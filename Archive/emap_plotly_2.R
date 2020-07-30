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


add_segments_variable_width <- function(p, x, xend, y, yend, width, color){
  tmp <- p
  for (i in 1:length(x)){
    edge_df_tmp <- tibble(x = x[i],
                          xend = xend[i],
                          y = y[i],
                          yend = yend[i],
                          width = width[i])
    print(edge_df_tmp)
    tmp <- tmp %>% 
      add_segments(data=edge_df_tmp, 
                 x=~edge_df_tmp$x,  xend=~edge_df_tmp$xend,
                 y=~edge_df_tmp$y, yend=~edge_df_tmp$yend, 
                 line = list(width = edge_df_tmp$width)) 
    print(tmp)
  }
  return(tmp)
}




plot_ly() %>% 
  add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                              y=edge_df$y, yend=edge_df$yend, width = edge_df$width)

x=edge_df$x
xend=edge_df$xend
y=edge_df$y
yend=edge_df$yend
width = edge_df$width


p <- plot_ly()
for (i in 1:length(x)){
  p <- p %>% 
    add_segments(data=tibble(x_column = x[i],
                             xend_column = xend[i],
                             y_column = y[i],
                             yend_column = yend[i],
                             width_column = width[i]), 
                 x=~x_column,  xend=~xend_column,
                 y=~y_column, yend=~yend_column, 
                 color = I("#AAAAAA"),
                 line = list(width = width[i])) 
}

add_segments_variable_width <- function(p, x, xend, y, yend, width#, color="#AAAAAA"
                                        ){
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

plot_ly() %>% 
  add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                              y=edge_df$y, yend=edge_df$yend, width = edge_df$width)


plot_ly() %>% 
  add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                              y=edge_df$y, yend=edge_df$yend, 
                              width = edge_df$width) %>% 
  add_markers(x = ~Xn, y = ~Yn,
              #mode = "markers",
              #text = vs$label, textposition = "bottom center",
              color = vs$color, 
              colors = "RdYlBu",
              #size = vs$size,
              #marker = list(size = vs$size),
              #sizes= c(100, 800),
              hoverinfo = "text",
              inherit = FALSE) %>% 
  add_text(x = ~Xn, y = ~Yn, text = vs$label, textposition = "bottom center",
           textfont = list(size = 12, family = 'Helvetica') )

plot_ly() %>% 
  # add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
  #                             y=edge_df$y, yend=edge_df$yend, 
  #                           width = edge_df$width) %>% 
  add_markers(x = ~Xn, y = ~Yn,
              #mode = "markers",
              #text = vs$label, textposition = "bottom center",
              color = ~vs$color, 
              colors = "RdYlBu",
              size = vs$size,
              marker = list(size = vs$size),
              #sizes= c(100, 800),
              hoverinfo = "text",
              inherit = FALSE)


plot_ly(x = ~Xn, y = ~Yn,
        #mode = "markers",
        #stroke = vs$color,
        #text = vs$label, textposition = "bottom center",
        color = vs$color, colors = "RdYlBu",
        #size = vs$size, 
        marker = list(size = vs$size),
        hoverinfo = "text") %>% 
  add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                              y=edge_df$y, yend=edge_df$yend, 
                              width = edge_df$width) %>% 
  add_text(text = vs$label, textposition = "bottom center",
           textfont = list(size = 12, family = 'Helvetica') )
plot_ly(x = ~Xn, y = ~Yn,
        #mode = "markers",
        #stroke = vs$color,
        #text = vs$label, textposition = "bottom center",
        color = vs$color, colors = "RdYlBu",
        #size = vs$size, 
        marker = list(size = vs$size),
        hoverinfo = "text") %>% 
  # add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
  #                             y=edge_df$y, yend=edge_df$yend, 
  #                             width = edge_df$width) %>% 
  add_text(text = vs$label, textposition = "bottom center",
           textfont = list(size = 12, family = 'Helvetica') ) %>% 
  add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                              y=edge_df$y, yend=edge_df$yend, 
                              width = edge_df$width) 
plot_ly(x = ~Xn, y = ~Yn,
        mode = "markers",
        #stroke = vs$color,
        #text = vs$label, textposition = "bottom center",
        color = vs$color, colors = "RdYlBu",
        #size = vs$size, 
        marker = list(size = vs$size),
        hoverinfo = "text") %>% 
  # add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
  #                             y=edge_df$y, yend=edge_df$yend, 
  #                             width = edge_df$width) %>% 
  add_text(text = vs$label, textposition = "bottom center",
           textfont = list(size = 12, family = 'Helvetica') ) %>% 
  add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                              y=edge_df$y, yend=edge_df$yend, 
                              width = edge_df$width) %>% 
  add_markers()

plot_ly(x = ~Xn, y = ~Yn,type = "scatter",
        mode = "markers",
        #stroke = vs$color,
        #text = vs$label, textposition = "bottom center",
        color = ~vs$color, colors = "RdYlBu",
        size = ~vs$size * 100, 
        #sizes = c(1,100),
        marker = list(size = vs$size),
        hoverinfo = "text") %>% 
  add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                              y=edge_df$y, yend=edge_df$yend, 
                              width = edge_df$width) %>% 
  add_text(text = vs$label, textposition = "bottom center",
           textfont = list(size = 12, family = 'Helvetica') )




plot_ly(x = ~Xn, y = ~Yn,
        type = "scatter",
        mode = "markers",
        color = ~vs$color, 
        colors = "RdYlBu",
        size = ~vs$size * 100, 
        #sizes = c(1,100),
        marker = list(size = vs$size),
        hoverinfo = "text") 

plot_ly(x = ~Xn, y = ~Yn,
        type = "scatter",
        mode = "markers",
        color = ~vs$color,
        colors = "RdYlBu",
        size = ~vs$size,
        sizes = c(100,800),
        text = vs$label,
        hoverinfo = "text",
        showlegend=FALSE
        )%>% 
  add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                              y=edge_df$y, yend=edge_df$yend, 
                              width = edge_df$width) %>% 
  add_markers(showlegend=FALSE) %>% 
  add_text(text = vs$label, textposition = "bottom center",
           textfont = list(size = 12, family = 'Helvetica') ) 


###############################################

i=1
edge_df_tmp <- tibble(x = x[i],
                      xend = xend[i],
                      y = y[i],
                      yend = yend[i],
                      width = width[i])


tmp <- plot_ly() %>% 
  add_segments(data=edge_df_tmp, 
               x=~edge_df_tmp$x,  xend=~edge_df_tmp$xend,
               y=~edge_df_tmp$y, yend=~edge_df_tmp$yend, 
               line = list(width = edge_df_tmp$width)) 
i = 2
edge_df_tmp <- tibble(x_ = x[i],
                     xend_ = xend[i],
                     y_ = y[i],
                     yend_ = yend[i],
                     width_ = width[i])
i = 1
#tmp %>% 
tmp <- plot_ly() %>% 
  add_segments(data=tibble(x_ = x[i],
                           xend_ = xend[i],
                           y_ = y[i],
                           yend_ = yend[i],
                           width_ = width[i]), 
               x=~x_,  xend=~xend_,
               y=~y_, yend=~yend_, 
               line = list(width = width)) 
i = 2

tmp %>% 
  add_segments(data=tibble(x_ = x[i],
                           xend_ = xend[i],
                           y_ = y[i],
                           yend_ = yend[i],
                           width_ = width[i]), 
               x=~x_,  xend=~xend_,
               y=~y_, yend=~yend_, 
               line = list(width = width)) 





i = 2
#plot_ly() %>% 
tmp %>% 
  add_segments(data=tibble(x = x[i],
                           xend = xend[i],
                           y = y[i],
                           yend = yend[i],
                           width = width[i]), 
               x=~x,  xend=~xend,
               y=~y, yend=~yend, 
               line = list(width = width)) 



plot_ly() %>%
  add_segments()

plot_ly() %>%
  add_segments(data=edge_df, 
               x=~edge_df$x,  xend=~edge_df$xend,
               y=~edge_df$y, yend=~edge_df$yend, 
               #width = edge_df$width,
               #size=edge_df$width
               color=edge_df$width
               #line = list(width = edge_df$width)
               ) %>% 
  add_markers(x = ~Xn, y = ~Yn,
              #mode = "markers",
              #text = vs$label, textposition = "bottom center",
              color = vs$color, colors = "RdYlBu",
              size = vs$size, marker = list(size = vs$size),
              hoverinfo = "text") %>% 
  add_text(x = ~Xn, y = ~Yn, text = vs$label, textposition = "bottom center",
           textfont = list(size = 12, family = 'Helvetica') )

edge_df2 <- edge_df[1,]
temp <- plot_ly() %>%
  add_segments(data=edge_df2, 
               x=~edge_df2$x,  xend=~edge_df2$xend,
               y=~edge_df2$y, yend=~edge_df2$yend, 
               width = edge_df$width[1],
               #size=edge_df$width
               #color="gray",
               line = list(width = 10, color = "#909090")
  )
edge_df3 <- edge_df[2,]

temp%>%
  add_segments(data=edge_df3, 
               x=~edge_df3$x,  xend=~edge_df3$xend,
               y=~edge_df3$y, yend=~edge_df3$yend, 
               #width = edge_df$width,
               #size=edge_df$width
               #color="gray",
               line = list(width = 5, color = "#030303")
  )




edge_shapes <- list()
for(i in 1:Ne) {
  v0 <- es[i,]$V1
  v1 <- es[i,]$V2
  
  edge_shape = list(
    type = "line",
    line = list(color = "#030303", width = E(g)$width[[i]]),
    x0 = Xn[node2index[[as.character(v0)]]],
    y0 = Yn[node2index[[as.character(v0)]]],
    x1 = Xn[node2index[[as.character(v1)]]],
    y1 = Yn[node2index[[as.character(v1)]]]
  )
}


















###########################
#network <- 
plot_ly(x = ~Xn, y = ~Yn,
        mode = "markers",
        stroke = vs$color,
        #text = vs$label, textposition = "bottom center",
        color = vs$color, colors = "RdYlBu",
        size = vs$size, 
        marker = list(size = vs$size),
        hoverinfo = "text") %>% 
  add_segments_variable_width(x=edge_df$x,  xend=edge_df$xend,
                              y=edge_df$y, yend=edge_df$yend, 
                              width = edge_df$width) %>% 
  add_text(text = vs$label, textposition = "bottom center",
           textfont = list(size = 12, family = 'Helvetica') )













edge_shapes <- list()
for(i in 1:Ne) {
  v0 <- es[i,]$V1
  v1 <- es[i,]$V2
  
  edge_shape = list(
    type = "line",
    line = list(color = "#030303", width = E(g)$width[[i]]),
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
plot_ly() %>%
  layout(shapes = edge_shapes) %>% 
  add_markers(x = ~Xn, y = ~Yn,
              #mode = "markers",
              #text = vs$label, textposition = "bottom center",
              color = vs$color, colors = "RdYlBu",
              size = vs$size, marker = list(size = vs$size),
              hoverinfo = "text") %>% 
  add_text(x = ~Xn, y = ~Yn, text = vs$label, textposition = "bottom center",
           textfont = list(size = 12, family = 'Helvetica') )

test = tibble(x=c(1,-1), xend = c(0,-0.5), y = c(0, -1), yend = c(0, -1))
plot_ly() %>%
  add_segments(data=test, x=~test$x,  xend=~test$xend, y=~test$y, yend=~test$yend, color="#000000") %>% 
  add_markers(x = ~Xn, y = ~Yn,
              #mode = "markers",
              #text = vs$label, textposition = "bottom center",
              color = vs$color, colors = "RdYlBu",
              size = vs$size, marker = list(size = vs$size),
              hoverinfo = "text") %>% 
  add_text(x = ~Xn, y = ~Yn, text = vs$label, textposition = "bottom center",
           textfont = list(size = 12, family = 'Helvetica') )



plot_ly() %>%
  add_markers(x = ~Xn, y = ~Yn,
              mode = "markers",
              #text = vs$label, textposition = "bottom center",
              color = vs$color, colors = "RdYlBu",
              size = vs$size, marker = list(size = vs$size),
              hoverinfo = "text") 


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



s <- read.csv("https://raw.githubusercontent.com/plotly/datasets/master/school_earnings.csv")
# order factor levels by men's income (plot_ly() will pick up on this ordering)
s$School <- factor(s$School, levels = s$School[order(s$Men)])

library(plotly)
fig <- plot_ly(s, color = I("gray80"))
fig <- fig %>% add_segments(x = ~Women, xend = ~Men, y = ~School, yend = ~School, showlegend = FALSE)
plot_ly() %>% 
  add_segments(x = s$Women, xend = s$WomenMen,
               y = s$WomenSchool, yend = s$WomenSchool,
               showlegend = FALSE)
plot_ly() %>% add_segments(data=s, x = ~Women, xend = ~Men, y = ~School, yend = ~School, showlegend = FALSE)
