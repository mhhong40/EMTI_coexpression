library(tidyverse)
library(here)

setwd(here::here("~/Documents/tnbc_coexpression"))

bt <- readRDS("tnbc_final.rds")
bn <- readRDS("norm_final.rds")
gs <- readRDS("gs_final.rds")

n_tnbc <- ncol(bt)
bt <- bt[rownames(bt) %in% gs,]

n_norm <- ncol(bn)
bn <- bn[rownames(bn) %in% gs,]
bn <- bn[match(rownames(bt), rownames(bn)),]

tnem <- cbind(bt, bn)

TNBC <- c(rep(1, each = n_tnbc), rep(0, each = n_norm))
Norm <- c(rep(0, each = n_tnbc), rep(1, each = n_norm))

TNBC <- factor(TNBC)
Norm <- factor(Norm)

design <- model.matrix(~ TNBC + Norm) 
design <- design[, -1]
colnames(design) <- c("Bulk TNBC", "Bulk Normal")

library(DGCA)
# Obtain differential coexpression z-scores
tn_norm <- ddcorAll(inputMat = tnem,
                   design = design,
                   compare = c("Bulk TNBC", "Bulk Normal"),
                   adjust = "perm", nPerm = 10)
tn_norm[, c("Gene1", "Gene2")] <- lapply(tn_norm[, c("Gene1", "Gene2")], as.character)
tn_norm$Pair <- paste(tn_norm$Gene1, tn_norm$Gene2)
tn_norm <- tn_norm[order(tn_norm$Pair), ]


## KRT5 (non-differentially expressed between TNBC and normal)
KRT5 <- tn_norm[which(grepl("KRT5", tn_norm$Pair)), ]
resEMT <- read_csv("tnbc_norm_resEMT.csv")
changeInd <- !!(KRT5$Gene1 != "KRT5") # KRT5 must be Gene1
KRT5[changeInd, c("Gene1", "Gene2")] <- KRT5[changeInd, c("Gene2", "Gene1")]

resEMT <- resEMT[match(KRT5$Gene2, resEMT$rowname), ]
dv_KRT5 <- data.frame(KRT5$Gene2, KRT5$zScoreDiff, resEMT$stat)
colnames(dv_KRT5) <- c("Gene", "dc_zscore", "de_wald")

# Figure 6A (all points)
krt5 <- ggplot(dv_KRT5, aes(x = dc_zscore, y = de_wald)) + 
  ggtitle("KRT5 (non-differentially expressed), all genes") +
  labs(x = "Differential coexpression z-score",
       y = "Differential expression Wald statistic") +
  theme_linedraw() +
  geom_point() +
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_smooth(method = 'lm')
krt5_cor <- cor.test(dv_KRT5$dc_zscore, dv_KRT5$de_wald, method = "spearman", exact = FALSE)
krt5_cor

# Removing non-differentially expressed genes (Wald statistic != 0, padj < 0.05)
# Note that all genes with |LFC| <= 1 have Wald statistic = 0 
dv_KRT5_DE <- cbind(dv_KRT5, resEMT$padj)
colnames(dv_KRT5_DE)[4] <- "DE_padj"
dv_KRT5_DE <- dv_KRT5_DE[which(dv_KRT5_DE$de_wald != 0 & dv_KRT5_DE$DE_padj < 0.05), ]

krt5_DE <- ggplot(dv_KRT5_DE, aes(x = dc_zscore, y = de_wald)) + 
  ggtitle("KRT5, differentially expressed genes only") +
  labs(x = "Differential coexpression z-score",
       y = "Differential expression Wald statistic") +
  theme_linedraw() +
  geom_point() +
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_smooth(method = 'lm')
krt5_DE_cor <- cor.test(dv_KRT5_DE$dc_zscore, dv_KRT5_DE$de_wald, method = "spearman", exact = FALSE)
krt5_DE_cor


## EPCAM (upregulated in TNBC)
EPCAM <- tn_norm[which(grepl("EPCAM", tn_norm$Pair)), ]
changeInd <- !!(EPCAM$Gene1 != "EPCAM") # EPCAM must be Gene1
EPCAM[changeInd, c("Gene1", "Gene2")] <- EPCAM[changeInd, c("Gene2", "Gene1")]

resEMT <- resEMT[match(EPCAM$Gene2, resEMT$rowname), ]
dv_EPCAM <- data.frame(EPCAM$Gene2, EPCAM$zScoreDiff, resEMT$stat)
colnames(dv_EPCAM) <- c("Gene", "dc_zscore", "de_wald")

# Figure 6B
epcam <- ggplot(dv_EPCAM, aes(x = dc_zscore, y = de_wald)) + 
  ggtitle("EPCAM (upregulated in TNBC), all genes") +
  labs(x = "Differential coexpression z-score",
       y = "Differential expression Wald statistic") +
  theme_linedraw() +
  geom_point() +
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_smooth(method = 'lm')
epcam_cor <- cor.test(dv_EPCAM$dc_zscore, dv_EPCAM$de_wald, method = "spearman", exact = FALSE)
epcam_cor

# Removing NDE genes
dv_EPCAM_DE <- cbind(dv_EPCAM, resEMT$padj)
colnames(dv_EPCAM_DE)[4] <- "DE_padj"
dv_EPCAM_DE <- dv_EPCAM_DE[which(dv_EPCAM_DE$de_wald != 0 & dv_EPCAM_DE$DE_padj < 0.05), ]

epcam_DE <- ggplot(dv_EPCAM_DE, aes(x = dc_zscore, y = de_wald)) + 
  ggtitle("EPCAM, differentially expressed genes only") +
  labs(x = "Differential coexpression z-score",
       y = "Differential expression Wald statistic") +
  theme_linedraw() +
  geom_point() +
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_smooth(method = 'lm')
epcam_DE_cor <- cor.test(dv_EPCAM_DE$dc_zscore, dv_EPCAM_DE$de_wald, method = "spearman", exact = FALSE)
epcam_DE_cor


## TGFB1I1 (downregulated in TNBC)
TGFB1I1 <- tn_norm[which(grepl("TGFB1I1", tn_norm$Pair)), ]
changeInd <- !!(TGFB1I1$Gene1 != "TGFB1I1") # TGFB1I1 must be Gene1
TGFB1I1[changeInd, c("Gene1", "Gene2")] <- TGFB1I1[changeInd, c("Gene2", "Gene1")]

resEMT <- read_csv("tnbc_norm_resEMT.csv")
resEMT <- resEMT[match(TGFB1I1$Gene2, resEMT$rowname), ]
dv_TGFB1I1 <- data.frame(TGFB1I1$Gene2, TGFB1I1$zScoreDiff, resEMT$stat)
colnames(dv_TGFB1I1) <- c("Gene", "dc_zscore", "de_wald")

# Figure 6C
tgfb1i1 <- ggplot(dv_TGFB1I1, aes(x = dc_zscore, y = de_wald)) + 
  ggtitle("TGFB1I1 (downregulated in TNBC), all genes") +
  labs(x = "Differential coexpression z-score",
       y = "Differential expression Wald statistic") +
  theme_linedraw() +
  geom_point() +
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_smooth(method = 'lm')
tgfb1i1_cor <- cor.test(dv_TGFB1I1$dc_zscore, dv_TGFB1I1$de_wald, method = "spearman", exact = FALSE)
tgfb1i1_cor

# Removing NDE genes
dv_TGFB1I1_DE <- cbind(dv_TGFB1I1, resEMT$padj)
colnames(dv_TGFB1I1_DE)[4] <- "DE_padj"
dv_TGFB1I1_DE <- dv_TGFB1I1_DE[which(dv_TGFB1I1_DE$de_wald != 0 & dv_TGFB1I1_DE$DE_padj < 0.05), ]

tgfb1i1_DE <- ggplot(dv_TGFB1I1_DE, aes(x = dc_zscore, y = de_wald)) + 
  ggtitle("TGFB1I1, differentially expressed genes only") +
  labs(x = "Differential coexpression z-score",
       y = "Differential expression Wald statistic") +
  theme_linedraw() +
  geom_point() +
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_smooth(method = 'lm')
tgfb1i1_DE_cor <- cor.test(dv_TGFB1I1_DE$dc_zscore, dv_TGFB1I1_DE$de_wald, method = "spearman", exact = FALSE)
tgfb1i1_DE_cor

# Figure 6
library(ggpubr)
# TO DO: get plots with non-differentially expressed genes removed
figure6 <- ggarrange(krt5, epcam, tgfb1i1, krt5_DE, epcam_DE, tgfb1i1_DE,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 2)
pdf(file = "Fig6.pdf", height = 8, width = 14)
figure6
dev.off()

save(list = c("krt5_cor", "krt5_DE_cor", "epcam_cor", "epcam_DE_cor", "tgfb1i1_cor", "tgfb1i1_DE_cor"),
     file = "diff_coexp_data.rda")
