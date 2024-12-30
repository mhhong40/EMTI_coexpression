library(tidyverse)
library(here)
library(hdWGCNA)

setwd(here::here("Documents/tnbc_coexpression"))

# Read in TNBC datasets and obtain their common genes
EpiMeta <- readRDS("EpiMeta.rds")
CycEpiMeta <- readRDS("CycEpiMeta.rds")
tnbc_bulk <- readRDS("tnbc_final.rds")

e <- GetMetacellObject(EpiMeta)
e <- e[, e$cell_type == "Epithelial"]  # 5445 epithelial metacells

c <- GetMetacellObject(CycEpiMeta)
c <- c[, c$cell_type == "Cycling Epithelial"]  # 2898 cycling epithelial metacells

# Split single-cell into datasets based on patient-based cohort
e1 <- e[, e$group == "TN_B1_4031"]
e2 <- e[, e$group == "TN_B1_0131"]
e3 <- e[, !(e$group %in% c("TN_B1_4031", "TN_B1_0131"))]  # Cohort 3 comprises numerous patient IDs

c1 <- c[, c$group == "TN_B1_4031"]
c2 <- c[, c$group == "TN_B1_0131"]
c3 <- c[, !(c$group %in% c("TN_B1_4031", "TN_B1_0131"))]

# Clear up unused memory
rm(EpiMeta, CycEpiMeta, e, c); gc()

# For each single-cell dataset, keep genes whose variance > 0. A more liberal threshold than bulk.
# Don't need to do this for bulk because variance thresholding was performed beforehand
library(matrixStats)
sc <- list(e1, e2, e3, c1, c2, c3)
scv <- list()
vkeep <- list()
for(i in 1:length(sc)){
  scv[[i]] <- rowVars(as.matrix(sc[[i]]@assays$RNA@data))
  vkeep[[i]] <- names(scv[[i]][which(scv[[i]] > 0)]) 
} 

# Subset metacell expression matrices to only EMT genes
emts <- read_tsv("EMTome_gene_list.tsv")
emts <- append(emts$gene, "BACH1") # LIN28A and NFE2L2 not in the data, so ignore

# The genes whose expression variance > 0 in all metacell expression matrices.
scgk <- Reduce(intersect, vkeep) 
rm(scv, vkeep)

# We can just take the intersection between EMT gene names in bulk and single-cell, 
# as all single-cell datasets have expression information for the same genes.
emtfinal <- rownames(tnbc_bulk)[rownames(tnbc_bulk) %in% intersect(rownames(tnbc_bulk), scgk)]
tnbc_bulk <- tnbc_bulk[rownames(tnbc_bulk) %in% emtfinal, ] 
scgk <- scgk[which(scgk %in% emtfinal)]
gs <- emtfinal
saveRDS(gs, "gs_tnbc.rds")

# Repeat but for normal breast datasets
LPMeta <- readRDS("LPMeta.rds")
MLMeta <- readRDS("MLMeta.rds")
BasalMeta <- readRDS("BasalMeta.rds")
norm_bulk <- readRDS("norm_final.rds")

lp <- GetMetacellObject(LPMeta)
lp <- lp[, lp$cell_type == "LP"]

ml <- GetMetacellObject(MLMeta)
ml <- ml[, ml$cell_type == "ML"]

basal <- GetMetacellObject(BasalMeta)
basal <- basal[, basal$cell_type == "Basal"]

sc <- list(lp, ml, basal)
scv <- list()
vkeep <- list()
for(i in 1:length(sc)){
  scv[[i]] <- rowVars(as.matrix(sc[[i]]@assays$RNA@data))
  vkeep[[i]] <- names(scv[[i]][which(scv[[i]] > 0)]) 
} 

scgk <- Reduce(intersect, vkeep) 
rm(scv, vkeep)

emtfinal <- rownames(norm_bulk)[rownames(norm_bulk) %in% intersect(rownames(norm_bulk), scgk)]
norm_bulk <- norm_bulk[rownames(norm_bulk) %in% emtfinal, ] 
scgk <- scgk[which(scgk %in% emtfinal)]
gs <- emtfinal
saveRDS(gs, "gs_norm.rds")

# Take the intersection of both
gs_tnbc <- readRDS("gs_tnbc.rds")
gs_norm <- readRDS("gs_norm.rds")
gs_final <- intersect(gs_tnbc, gs_norm)

saveRDS(gs_final, "gs_final.rds")
