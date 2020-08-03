library(clusterProfiler)
library(ReactomePA)
library(tidyverse)

library(org.Hs.eg.db)
library(org.Mm.eg.db) 

mouse_gene_set <- read_csv("./sample_data/mouse_sample_data.csv") %>% 
  pull(ENTREZID)


species2keggorganism <- c("hsa", "mmu")
names(species2keggorganism) <- c("Human", "Mouse")

species2OrgDB <- c(org.Hs.eg.db, org.Mm.eg.db)
names(species2OrgDB) <- c("Human", "Mouse")

input_species <- "Mouse"

species2keggorganism[[input_species]]
species2OrgDB[[input_species]]

mouse_enrich_kegg <- enrichKEGG(gene=mouse_gene_set, pvalueCutoff=0.05, organism = species2keggorganism[[input_species]])# "hsa"
mouse_enrich_reactome <- enrichPathway(gene=mouse_gene_set, pvalueCutoff=0.05, readable = TRUE, organism = input_species)
mouse_enrich_GO <- enrichGO(gene=mouse_gene_set, pvalueCutoff=0.05, readable = TRUE,OrgDb = species2OrgDB[[input_species]], ont = "BP")

