library(tidyverse)
library(here)
library(data.table)

# Using the UCSCXenaTools package to download the gene expected counts matrix yielded an incomplete one (only 2456 genes)
# Used this link instead: https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz

## TNBC patients
primtnbcs <- readRDS("primtnbcs.rds")
filt <- c(primtnbcs, "sample")
tc <- fread("TcgaTargetGtex_gene_expected_count.gz", select = filt)  # TCGA-AR-A0U1-01 not found in col name header

# Include only protein-coding and lncRNA genes
probemap <- fread("zz_gencode.v23.annotation.txt", select = c(1,2))
texprALL <- merge(probemap, tc, by.x = "id", by.y = "sample")

# Ensure no duplicates
texprALL <- texprALL[!(duplicated(texprALL$gene) | duplicated(texprALL$gene, fromLast = TRUE)),]
# write.csv(texprALL, file = "texprALL.csv")

# Protein-coding
genesPC <- fread("zz_gene.protein.coding.csv")
exprPC <- subset(texprALL, gene %in% genesPC$Gene_Symbol)

# LncRNA
genesLNR <- rtracklayer::import("gencode.v23.long_noncoding_RNAs.gtf")
genesLNR_df <- as.data.frame(genesLNR)
exprLNR <- subset(texprALL, gene %in% genesLNR_df$gene_name)

# Combined expression matrix
texprFinal <- as.data.frame(rbind(exprPC, exprLNR))
texprFinal <- texprFinal[!(duplicated(texprFinal$gene) | duplicated(texprFinal$gene, fromLast = TRUE)),]
rownames(texprFinal) <- texprFinal$gene
texprFinal <- texprFinal[, -c(1:2)]

# Convert to integer counts
texprFinal <- round(((2^texprFinal) - 1), digits = 0)
saveRDS(texprFinal, "texprFinal.rds", compress = "xz")

## Normal patients
f2GTEx <- readRDS("f2GTEx.rds")
f2GTEx <- f2GTEx[which(f2GTEx$gender == "Female"),]
filt <- c(f2GTEx$sample, "sample")

nc <- fread("TcgaTargetGtex_gene_expected_count.gz", select = filt)  
nexprALL <- merge(probemap, nc, by.x = "id", by.y = "sample")

# Ensure no duplicates
nexprALL <- nexprALL[!(duplicated(nexprALL$gene) | duplicated(nexprALL$gene, fromLast = TRUE)),]
write.csv(nexprALL, file = "nexprALL.csv")

exprPC <- subset(nexprALL, gene %in% genesPC$Gene_Symbol)
exprLNR <- subset(nexprALL, gene %in% genesLNR_df$gene_name)

nexprFinal <- as.data.frame(rbind(exprPC, exprLNR))
nexprFinal <- nexprFinal[!(duplicated(nexprFinal$gene) | duplicated(nexprFinal$gene, fromLast = TRUE)),]
rownames(nexprFinal) <- nexprFinal$gene
nexprFinal <- nexprFinal[, -c(1:2)]
colnames(nexprFinal) <- gsub("\\-[[:digit:]]{4}-SM-[[:alnum:]]{5}", "", colnames(nexprFinal))

# Convert to integer counts
nexprFinal <- round(((2^nexprFinal) - 1), digits = 0)
saveRDS(nexprFinal, "nexprFinal.rds", compress = "xz")

## Perform differential gene expression analysis using DESeq2
library(DESeq2)

# Merged counts matrix
fc <- cbind(nexprFinal, texprFinal)
saveRDS(fc, "fc.rds", compress = "xz")

bt <- readRDS("tnbc_final.rds")
n_norm <- ncol(nexprFinal) # equals 79
n_tnbc <- ncol(bt) # equals 86

fc <- fc[, c(1:ncol(nexprFinal), which(colnames(fc) %in% colnames(bt)))]
cd <- data.frame(row.names = colnames(fc), condition = c(rep("normal", n_norm), rep("TNBC", n_tnbc)))
dds <- DESeqDataSetFromMatrix(countData = fc,
                              colData = cd,
                              design = ~ condition)

# Taken from the DESeq2 vignette - pre-filtering to save memory + possibly improve visualizations
smallestGroupSize <- n_norm # the number of normal samples, the smallest group size
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "normal")

dds <- DESeq(dds)
res05 <- results(dds, lfcThreshold = 1, altHypothesis = "greaterAbs", alpha = 0.05)
summary(res05)

# ylim <- c(-15,15)
# plotMA(res05, ylim=ylim); drawLines()

library(EnhancedVolcano)
gs <- readRDS("gs_final.rds")
EnhancedVolcano(res05,
                lab = rownames(res05),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = gs,
                title = "TNBC versus normal breast",
                pCutoff = 10e-32, 
                max.overlaps = 30)
resEMT <- res05[rownames(res05) %in% gs,]
resEMT <- resEMT[order(resEMT$pvalue),]
resEMT <- rownames_to_column(as.data.frame(resEMT))
write_csv(resEMT, file = "tnbc_norm_resEMT.csv")
