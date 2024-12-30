library(tidyverse)
library(here)
library(readxl)
library(data.table)

setwd(here::here("Documents/tnbc_coexpression"))

# Using the UCSCXenaTools package to download the gene expected counts matrix yielded an incomplete one (only 2456 genes)
# Used this link instead: https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz

## Get ER+ patients
clinical <- read_tsv("SuppTable1-Table 1.tsv", skip = 1, col_names = TRUE)
erpat <- filter(clinical, Gender == "FEMALE" & ER_Status == "Positive")
erpat <- erpat$tcga_id %>% 
  paste("-01", sep = "")
saveRDS(erpat, file = "erpat.rds")

probemap <- fread("zz_gencode.v23.annotation.txt", select = c(1,2))
erpc <- fread("TcgaTargetGtex_gene_expected_count.gz", select = c(erpat, "sample"))  # TCGA-BH-A0AY-01, E2-A108-01, and E2-A1IP-01 not found
erexprALL <- merge(probemap, erpc, by.x = "id", by.y = "sample")

## Obtain gene lists
# Protein-coding
genesPC <- fread("zz_gene.protein.coding.csv")
exprPC <- subset(erexprALL, gene %in% genesPC$Gene_Symbol)

# LncRNA
genesLNR <- rtracklayer::import("gencode.v23.long_noncoding_RNAs.gtf")
genesLNR_df <- as.data.frame(genesLNR)
exprLNR <- subset(erexprALL, gene %in% genesLNR_df$gene_name)

# Ensure no duplicates
erexprALL <- erexprALL[!(duplicated(erexprALL$gene) | duplicated(erexprALL$gene, fromLast = TRUE)),]
# write.csv(erexprALL, file = "erexprALL.csv")

exprPC <- subset(erexprALL, gene %in% genesPC$Gene_Symbol)
exprLNR <- subset(erexprALL, gene %in% genesLNR_df$gene_name)

erexprFinal <- as.data.frame(rbind(exprPC, exprLNR))
erexprFinal <- erexprFinal[!(duplicated(erexprFinal$gene) | duplicated(erexprFinal$gene, fromLast = TRUE)),]
rownames(erexprFinal) <- erexprFinal$gene
erexprFinal <- erexprFinal[, -c(1:2)]

# Convert to integer counts
erexprFinal <- round(((2^erexprFinal) - 1), digits = 0)
saveRDS(erexprFinal, "erexprFinal.rds", compress = "xz")

## Perform differential gene expression analysis using DESeq2
library(DESeq2)

# Merged counts matrix
texprFinal <- readRDS("texprFinal.rds")
fc <- cbind(erexprFinal, texprFinal)
saveRDS(fc, "er_t_fc.rds", compress = "xz")

bt <- readRDS("tnbc_final.rds")
n_er <- ncol(erexprFinal) # equals 590
n_tnbc <- ncol(bt) # equals 86

fc <- fc[, c(1:n_er, which(colnames(fc) %in% colnames(bt)))]
cd <- data.frame(row.names = colnames(fc), condition = c(rep("ER", n_er), rep("TNBC", n_tnbc)))
dds <- DESeqDataSetFromMatrix(countData = fc,
                              colData = cd,
                              design = ~ condition)

# Taken from the DESeq2 vignette - pre-filtering to save memory + possibly improve visualizations
smallestGroupSize <- n_tnbc # the number of TNBC samples, the smallest group size
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "ER") # desire for our LFCs to be log2(TNBC/ER) 

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
                title = "TNBC versus ER positive",
                pCutoff = 10e-32, 
                max.overlaps = 30)
resEMT <- res05[rownames(res05) %in% gs,]
resEMT <- resEMT[order(resEMT$pvalue),]
resEMT <- rownames_to_column(as.data.frame(resEMT))
write_csv(resEMT, file = "er_tnbc_resEMT.csv")
