library(tidyverse)
library(here)
library(R.utils)
library(Seurat)

## Load data
setwd(here::here("Documents/tnbc_coexpression"))
Norms <- readRDS("SeuratObject_NormEpiSub.rds")

### Form metacells to reduce sparsity ----
library(hdWGCNA)
set.seed(12345)
Cluster <- as.integer(Norms@meta.data$integrated_snn_res.0.015)
DimPlot(Norms, group.by = "seurat_clusters")
Norms <- Norms[, -which(Norms@meta.data$seurat_clusters == 3)]

# Clusters 0 and 2 are non-cycling; 1 is the distinct cycling epithelial cluster
NewCIs <- c("LP", "ML", "Basal")
names(NewCIs) <- levels(Norms)
Norms <- RenameIdents(Norms, NewCIs)
Norms@meta.data$cell_type <- Idents(Norms)

# SetupForWGCNA needed in order to form metacells
# This will not affect our downstream analyses.
Norms <- SetupForWGCNA(Norms, gene_select = "variable",
                          fraction = 0.05, wgcna_name = "emt")

# Non-integrated values desired for analyses 
DefaultAssay(Norms) <- "RNA"
NMetacells <- MetacellsByGroups(seurat_obj = Norms,
                                 group.by = c("cell_type", "group"),
                                 reduction = "pca",
                                 k = 20,
                                 max_shared = 10,
                                 ident.group = "cell_type")

# Divide each expression value by total counts for the cell, multiply by 10^4, then ln(x + 1)
NMetacells <- NormalizeMetacells(NMetacells)

NMetacells <- ScaleMetacells(NMetacells, features = VariableFeatures(NMetacells))
NMetacells <- RunPCAMetacells(NMetacells, features = VariableFeatures(NMetacells@assays$integrated))
NMetacells <- RunHarmonyMetacells(NMetacells, group.by.vars = 'group')
NMetacells <- RunUMAPMetacells(NMetacells, reduction = 'harmony', dims=1:15)

rm(Norms); gc()

## Normal cells are integrated well
# p1 <- DimPlotMetacells(NMetacells, reduction = "umap", group.by= "group") + ggtitle("Patient of Origin")
# p2 <- DimPlotMetacells(NMetacells, reduction = "umap", group.by= "cell_type") + ggtitle("Cell Lineage")
# p1 | p2

LPMeta <- SetDatExpr(NMetacells,
                      group_name = "LP",
                      group.by = "cell_type",
                      assay = "RNA",
                      slot = "data")
saveRDS(LPMeta, "LPMeta.rds", compress = "xz")

MLMeta <- SetDatExpr(NMetacells,
                         group_name = "ML",
                         group.by = "cell_type",
                         assay = "RNA",
                         slot = "data")
saveRDS(MLMeta, "MLMeta.rds", compress = "xz")

BasalMeta <- SetDatExpr(NMetacells,
                     group_name = "Basal",
                     group.by = "cell_type",
                     assay = "RNA",
                     slot = "data")
saveRDS(BasalMeta, "BasalMeta.rds", compress = "xz")