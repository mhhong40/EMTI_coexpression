library(tidyverse)
library(here)

setwd(here::here("Documents/tnbc_coexpression"))

# Gather GTEx data that's harmonized with TARGET/TCGA
# Synced TPMs 
library(UCSCXenaTools)
library(matrixStats)
library(data.table)
library(R.utils)
library(rtracklayer)

# Merge the two datasets
# TO DO: the expected counts matrix has one fewer desired sample than TPM; re-run normal scripts :/
f2GTEx <- readRDS("f2GTEx.rds")
f2GTEx <- f2GTEx[which(f2GTEx$gender == "Female"),] # Only female patients are desired

# GTEX-12WSK-2226-SM-5GCO5 is missing from the expected count data used for DEA
# So it needs to be excluded from coexpression analysis, too
filt <- c(f2GTEx$sample, "sample")
filt <- filt[which(filt != "GTEX-12WSK-2226-SM-5GCO5")]
tpmbysamp <- fread("TcgaTargetGtex_rsem_gene_tpm.gz", select = filt)

## Include only protein-coding and lncRNA genes
probemap <- fread("zz_gencode.v23.annotation.txt", select = c(1,2))
exprALL <- merge(probemap, tpmbysamp, by.x = "id", by.y = "sample")

# Ensure no duplicates
exprALL <- exprALL[!(duplicated(exprALL$gene) | duplicated(exprALL$gene, fromLast = TRUE)),]
# write.csv(exprALL, file = "exprALL.csv")

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
saveRDS(exprFinal, "gtex_exprFinal.rds", compress = "xz")

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

# Creating a unified metadata matrix just in case.
# Sample ID
sampatt <- read.delim("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header = TRUE) 
sampatt <- sampatt[sampatt$SAMPID %in% f2GTEx$sample, ]
sampatt <- sampatt[match(colnames(exprF2), sampatt$SAMPID), ]
sampatt$SAMPID <- gsub("\\-[[:digit:]]{4}-SM-[[:alnum:]]{5}", "", sampatt$SAMPID)
colnames(exprF2) <- gsub("\\-[[:digit:]]{4}-SM-[[:alnum:]]{5}", "", colnames(exprF2))

# Sample sex
# Technically not necessary because males were filtered out...? Whatever
pheno <- read.delim("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", header = TRUE)
pheno <- pheno[pheno$SUBJID %in% sampatt$SAMPID,]
pheno <- pheno[match(sampatt$SAMPID, pheno$SUBJID), ]
pheno$SEX <- gsub("1", "Male", pheno$SEX)
pheno$SEX <- gsub("2", "Female", pheno$SEX)
rownames(pheno) <- pheno$SUBJID
rownames(sampatt) <- sampatt$SAMPID

# Complete metadata matrix
metadata <- cbind(pheno, sampatt)
saveRDS(metadata, file = "metadata.rds")

# Scale resulting expression matrix
norm_scaled <- scale(exprF2)
saveRDS(norm_scaled, file = "norm_scaled.rds", compress = "xz")

# Subset expression matrix to EMT genes (1153 initially)
emts <- read_tsv("EMTome_gene_list.tsv")
emts <- append(emts$gene, c("BACH1", "LIN28A", "NFE2L2"))  # Don't forget
norm_scaled_corrected <- norm_scaled[rownames(norm_scaled) %in% emts, ]  # 816 genes retained
saveRDS(norm_scaled_corrected, file = "norm_final.rds", compress = "xz")
