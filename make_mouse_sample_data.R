library(tidyverse)


counts <- read.table("http://bowtie-bio.sourceforge.net/recount/countTables/katzmouse_count_table.txt", header = TRUE, row.names = 1)
counts <- counts[rowSums(counts) > 10, ]    # 低発現量のデータを取り除く
group <- c("CUGBP1", "control", "CUGBP1", "control")

library(edgeR)
d <- DGEList(counts = counts, group = group)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
result <- exactTest(d)
table <- as.data.frame(topTags(result, n = nrow(counts)))

library(clusterProfiler)
library(org.Mm.eg.db)   # マウスのアノテーションデータ

all.genes <- rownames(table)
is.degs <- all.genes[table$FDR < 0.01]

all.genes.entrez <- bitr(all.genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
is.degs.entrez <- bitr(is.degs, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Mm.eg.db)

is.degs.entrez %>% 
  dplyr::select(ENTREZID, ENSEMBL) %>% 
  write_csv("./sample_data/mouse_sample_data.csv")
