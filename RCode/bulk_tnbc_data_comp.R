library(tidyverse)
library(here)

setwd(here::here("Documents/tnbc_coexpression"))

# Gather TCGA data that's harmonized with TARGET/GTEx
# Synced TPMs 
library(UCSCXenaTools)
library(matrixStats)
library(data.table)
library(R.utils)
library(rtracklayer)

matchedtpms <- XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGtex_rsem_gene_tpm")

XenaQuery(matchedtpms) %>%
  XenaDownload(destdir = getwd())

# Merge the two datasets
primtnbcs <- readRDS("primtnbcs.rds")
filt <- c(primtnbcs, "sample")
tpmbysamp <- fread("TcgaTargetGtex_rsem_gene_tpm.gz", select = filt)  # TCGA-AR-A0U1-01 not found in col name header

## Include only protein-coding and lncRNA genes
probemap <- fread("zz_gencode.v23.annotation.txt", select = c(1,2))
exprALL <- merge(probemap, tpmbysamp, by.x = "id", by.y = "sample")

# Ensure no duplicates
exprALL <- exprALL[!(duplicated(exprALL$gene) | duplicated(exprALL$gene, fromLast = TRUE)),]
write.csv(exprALL, file = "exprALL.csv")

# Protein-coding
genesPC <- fread("zz_gene.protein.coding.csv")
exprPC <- subset(exprALL, gene %in% genesPC$Gene_Symbol)

# LncRNA
genesLNR <- rtracklayer::import("gencode.v23.long_noncoding_RNAs.gtf")
genesLNR_df <- as.data.frame(genesLNR)
exprLNR <- subset(exprALL, gene %in% genesLNR_df$gene_name)

# Combined expression matrix
exprFinal <- as.data.frame(rbind(exprPC, exprLNR))
exprFinal <- exprFinal[!(duplicated(exprFinal$gene) | duplicated(exprFinal$gene, fromLast = TRUE)),]
rownames(exprFinal) <- exprFinal$gene
exprFinal <- exprFinal[, -c(1:2)]
saveRDS(exprFinal, "tnbc_exprFinal.rds", compress = "xz")

# Remove genes with variance < 25th percentile of variance to reduce noise
gv <- data.frame(genes = rownames(exprFinal),
                 variances = rowVars(as.matrix(exprFinal)))
gq25 <- quantile(gv$variances, probs = 0.25)

# Histogram of gene expression variances, if desired
# ggplot(gv, aes(x = variances)) + 
#  geom_histogram(binwidth = 0.4, color="black", fill="lightgrey") + 
#  geom_vline(aes(xintercept = gq25),
#             color = "red", linetype = "dashed") +
#  theme_classic()

gn <- gv[gv$variances > gq25, ]$genes # 25432 genes remaining
exprF2 <- as.matrix(exprFinal[rownames(exprFinal) %in% gn, ])

# Correct for covariate (age) using limma. Minimal cleaning here.
library(limma)
TCGAclinical <- read_tsv("TNBC_clinicalMatrix.tsv")
TCGAc1 <- TCGAclinical[which(!is.na(TCGAclinical$Age_at_Initial_Pathologic_Diagnosis_nature2012)), ]  # 100 samples left
TCGAc1 <- TCGAc1[TCGAc1$sampleID %in% colnames(exprF2), ]

# Scale expression matrix
tnbc_scaled <- scale(exprF2[, which(colnames(exprF2) %in% TCGAc1$sampleID)])
saveRDS(exprF2, file = "tnbc_scaled.rds", compress = "xz")

ages <- TCGAc1$Age_at_Initial_Pathologic_Diagnosis_nature2012
tnbc_scaled_corrected <- removeBatchEffect(tnbc_scaled, covariates = ages)

# Visualize removal of batch effect, if desired
oldpar <- par(mfrow = c(1, 2))
plotMDS(tnbc_scaled, main = "Original", gene.selection = "common", pch = 16)
plotMDS(tnbc_scaled_corrected, main = "Age corrected", gene.selection = "common", pch = 16)
par(oldpar) 

# Subset expression matrix to EMT genes (1153 initially)
emts <- read_tsv("EMTome_gene_list.tsv")
emts <- append(emts$gene, c("BACH1", "LIN28A", "NFE2L2"))  # Don't forget
tnbc_scaled_corrected <- tnbc_scaled_corrected[rownames(tnbc_scaled_corrected) %in% emts, ]  # 860 genes retained
saveRDS(tnbc_scaled_corrected, file = "tnbc_final.rds", compress = "xz")