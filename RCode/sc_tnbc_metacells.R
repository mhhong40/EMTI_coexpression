library(tidyverse)
library(here)
library(R.utils)
library(Seurat)

## Load data
setwd(here::here("Documents/tnbc_coexpression"))
TNBCs <- readRDS("SeuratObject_TNBC.rds")
TNBCs <- RunUMAP(TNBCs, reduction = "pca", dims = 1:30)

DimPlot(TNBCs, reduction = "tsne", group.by = "seurat_clusters")

# Subset epithelial cells (high EPCAM expression)
FeaturePlot(TNBCs, "EPCAM", reduction = "umap")

# Note that this increments each seurat_clusters up 1 (0 --> 1, 1 --> 2, etc.)
Cluster <- as.integer(TNBCs@meta.data$seurat_clusters)
Epi <- Cluster %in% c(1, 3)

Samples <- c("TN-B1-0554", "TN-B1-0177", "TN-0135", "TN-B1-4031", 
             "TN-B1-0131", "TN-0126", "TN-0106", "TN-0114-T2")
SamplesComb <- gsub("-","_",Samples)

### Read in the data 
## This code is directly taken from Pal et al. (2021)'s TNBC script.
## Default Seurat parameters for all functions. 
DGE <- paste0("dge_", SamplesComb)
DD <- paste0("dd_", SamplesComb)

cellNamesEpi <- rownames(TNBCs@meta.data)[Epi]
DGEEpi <- paste0("dge_epi_", SamplesComb)
DDEpi <- paste0("dd_epi_", SamplesComb)
for(i in 1:length(SamplesComb)) {
  d <- TNBCs[, grep(SamplesComb[i], colnames(TNBCs))] 
  d <- d[, colnames(d) %in% cellNamesEpi] 
  keep1 <- rowSums(d@assays$RNA@counts > 0) >= ncol(d)*0.01
  eval( parse(text=paste0(DDEpi[i],"<- d")) )
  d <- d[keep1, ]
  eval( parse(text=paste0(DGEEpi[i],"<- d")) )
}

# HVGs
SamplesCombEpi <- paste0(SamplesComb, "_Epi")
CombSeuratEpi <- list()
for(i in 1:length(SamplesCombEpi)) {
  d <- get(DGEEpi[i])
  so <- CreateSeuratObject(d@assays$RNA@counts, project=SamplesCombEpi[i])
  so <- NormalizeData(so)
  so <- FindVariableFeatures(so, selection.method="vst", nfeatures=1500)
  so <- ScaleData(so)
  so@meta.data$group <- SamplesComb[i]
  CombSeuratEpi[[i]] <- so
  eval( parse(text=paste0(SamplesCombEpi[i],"<- so")) )
}
names(CombSeuratEpi) <- SamplesCombEpi

# Integration
dimUsed <- 30
CSEBig <- CombSeuratEpi[c(1:6)]
AnchorsEpi <- FindIntegrationAnchors(object.list = CSEBig, dims = 1:dimUsed,
                                     anchor.features = 1000, scale = TRUE, k.anchor = 5, k.filter = 30,
                                     k.score = 20, max.features = 100)
TNBCEpi <- IntegrateData(anchorset = AnchorsEpi, dims = 1:dimUsed, k.weight = 100)
rm(AnchorsEpi)

# Dimensional reduction
dimUsed <- 30
DefaultAssay(TNBCEpi) <- "integrated"
TNBCEpi <- ScaleData(TNBCEpi, verbose=FALSE)
TNBCEpi <- RunPCA(TNBCEpi, npcs=dimUsed, verbose=FALSE)
TNBCEpi <- RunUMAP(TNBCEpi, reduction = "pca", dims=1:dimUsed)

GroupEpi <- factor(TNBCEpi@meta.data$group, levels=SamplesComb)
plotOrd <- sample(ncol(TNBCEpi))

# Clustering
resolution <- 0.1
TNBCEpi <- FindNeighbors(TNBCEpi, dims=1:dimUsed, verbose=FALSE)
TNBCEpi <- FindClusters(TNBCEpi, resolution=0.1, verbose=FALSE)
ClusterEpi <- as.integer(TNBCEpi@meta.data$seurat_clusters)
ncls <- length(table(ClusterEpi)) # number of clusters

DimPlot(TNBCEpi, group.by = "integrated_snn_res.0.1")
p3 <- FeaturePlot(TNBCEpi, c("MKI67", "CDK1"))
tm <- FindAllMarkers(TNBCEpi, only.pos = TRUE)  # cluster 1 is a discrete cluster of cycling cells

NewCIs <- c("Epithelial", "Cycling Epithelial", "Epithelial")
names(NewCIs) <- levels(TNBCEpi)
TNBCEpi <- RenameIdents(TNBCEpi, NewCIs)
p4 <- DimPlot(TNBCEpi, reduction = "umap", label = FALSE, pt.size = 0.5)
saveRDS(TNBCEpi, file = "TNBCEpi.rds", compress = "xz")

p3 + p4 

### Form metacells to reduce sparsity
library(hdWGCNA)
set.seed(12345)
Cluster <- as.integer(TNBCEpi@meta.data$integrated_snn_res.0.1)

# Clusters 0 and 2 are non-cycling; 1 is the distinct cycling epithelial cluster
NewCIs <- c("Epithelial", "Cycling Epithelial", "Epithelial")
names(NewCIs) <- levels(TNBCEpi)
TNBCEpi <- RenameIdents(TNBCEpi, NewCIs)
TNBCEpi@meta.data$cell_type <- Idents(TNBCEpi)

# SetupForWGCNA needed in order to form metacells
# This will not affect our downstream analyses.
TNBCEpi1 <- SetupForWGCNA(TNBCEpi, gene_select = "variable", wgcna_name = "emt")

# Non-integrated values desired for analyses 
DefaultAssay(TNBCEpi1) <- "RNA"
TEMetacells <- MetacellsByGroups(seurat_obj = TNBCEpi1,
                             group.by = c("cell_type", "group"),
                             reduction = "pca",
                             k = 20,
                             max_shared = 10,
                             ident.group = "cell_type")

# Divide each expression value by total counts for the cell, multiply by 10^4, then ln(x + 1)
TEMetacells <- NormalizeMetacells(TEMetacells)

TEMetacells <- ScaleMetacells(TEMetacells, features = VariableFeatures(TEMetacells))
TEMetacells <- RunPCAMetacells(TEMetacells, features = VariableFeatures(TEMetacells@assays$integrated))
TEMetacells <- RunHarmonyMetacells(TEMetacells, group.by.vars = 'group')
TEMetacells <- RunUMAPMetacells(TEMetacells, reduction = 'harmony', dims=1:15)

# SF 1 (metacell segregation by patient/cycling status)
p1 <- DimPlotMetacells(TEMetacells, reduction = "umap", group.by= "group") + ggtitle("Patient of Origin")
p2 <- DimPlotMetacells(TEMetacells, reduction = "umap", group.by= "cell_type") + ggtitle("Cell Cycling Status")
p1 | p2

EpiMeta <- SetDatExpr(TEMetacells,
                          group_name = "Epithelial",
                          group.by = "cell_type",
                          assay = "RNA",
                          slot = "data")
saveRDS(EpiMeta, "EpiMeta.rds", compress = "xz")

CycEpiMeta <- SetDatExpr(TEMetacells,
                         group_name = "Cycling Epithelial",
                         group.by = "cell_type",
                         assay = "RNA",
                         slot = "data")
saveRDS(CycEpiMeta, "CycEpiMeta.rds", compress = "xz")