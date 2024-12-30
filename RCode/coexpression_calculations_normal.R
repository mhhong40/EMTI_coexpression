library(tidyverse)
library(here)
library(hdWGCNA)

# Load TNBC datasets
# Scale the single-cell datasets
LPMeta <- readRDS("LPMeta.rds")
MLMeta <- readRDS("MLMeta.rds")
BasalMeta <- readRDS("BasalMeta.rds")
norm_bulk <- readRDS("norm_final.rds")

lp <- GetMetacellObject(LPMeta)
lp <- lp[, lp$cell_type == "LP"]  # 9196 LP metacells

ml <- GetMetacellObject(MLMeta)
ml <- ml[, ml$cell_type == "ML"]  # 7185 ML metacells

ba <- GetMetacellObject(BasalMeta)
ba <- ba[, ba$cell_type == "Basal"] # 2932 basal metacells

# Clear up unused memory
rm(LPMeta, MLMeta, BasalMeta); gc()

# Scale datasets, subset data to final gene list, alphabetize gene names
nsc <- list(lp, ml, ba)
gs <- readRDS("gs_final.rds")
nsc <- lapply(nsc, function(x) {
  
  x <- x@assays$RNA@scale.data
  x <- x[which(rownames(x) %in% gs), ]
  x <- x[order(rownames(x)), ]
})
norm_bulk <- norm_bulk[which(rownames(norm_bulk) %in% gs), ]
norm_bulk <- norm_bulk[order(rownames(norm_bulk)), ]

# Correlation calculation 
library(rstatix)
all_norm <- append(list(norm_bulk), nsc)
names(all_norm) <- c("Bulk Normal", "Luminal Progenitor", "Mature Luminal", "Basal Epithelial")
ct <- cdist_norm <- list()
for(i in 1:length(all_norm)){
  d <- as.data.frame(t(all_norm[[i]]))
  ct[[i]] <- d %>% cor(method = "pearson")
  cdist_norm[[i]] <- melt(ct[[i]])
  colnames(cdist_norm[[i]]) <- c("Gene1", "Gene2", "Correlation")
  cdist_norm[[i]] <- cdist_norm[[i]][which(cdist_norm[[i]]$Gene1 != cdist_norm[[i]]$Gene2), ] # Remove self-correlations
}
names(ct) <- names(cdist_norm) <- names(all_norm)

# Retrieve only top 5% highly correlated gene pairs
sig <- list()
for(i in 1:length(all_norm)) {
  
  # Calculate cutoffs
  lower <- quantile(cdist_norm[[i]]$Correlation, probs = 0.025)
  upper <- quantile(cdist_norm[[i]]$Correlation, probs = 0.975)
  
  # Retrieve significant correlations
  sig[[i]] <- cdist_norm[[i]][cdist_norm[[i]]$Correlation < lower | cdist_norm[[i]]$Correlation > upper, ]
}
names(sig) <- names(all_norm)

# Alphabetize gene pair names
for(i in 1:length(all_norm)){
  
  s <- sig[[i]]
  
  s[, c("Gene1", "Gene2")] <- lapply(s[, c("Gene1", "Gene2")], as.character)
  changeInd <- !!(s$Gene1 > s$Gene2)
  s[changeInd, c("Gene1", "Gene2")] <- s[changeInd, c("Gene2", "Gene1")]
  
  s$Pair <- NULL
  s$Pair <- paste(s$Gene1, s$Gene2)
  s <- s[order(s$Pair), ]
  s <- distinct(s)
  
  sig[[i]] <- s
}

# Number concordances between datasets
n_pairs <- length(all_norm) * (length(all_norm) - 1) * 1/2
conc_only_norm <- list()
n_conc <- data.frame(matrix(nrow = n_pairs, ncol = 3))
colnames(n_conc) <- c("Dataset1", "Dataset2", "Num_concordant")
k <- 1 # Counter for the conc_only_norm list
for(i in 1:(length(all_norm) - 1)) {
  
  for(j in (i + 1):length(all_norm)) {
    
    # Retrieve significant correlations for each set
    d1 <- sig[[i]]
    d2 <- sig[[j]]
    
    # Gene pairs significantly correlated in both sets
    d1 <- d1[which(d1$Pair %in% intersect(d1$Pair, d2$Pair)), ]
    d2 <- d2[which(d2$Pair %in% intersect(d1$Pair, d2$Pair)), ]
    
    # Calculate number of concordant pairs between datasets
    conc_only_norm[[k]] <- cbind(d1[, 1:2], d1$Pair, d1$Correlation, d2$Correlation)
    
    co <- conc_only_norm[[k]]
    colnames(co)[3:5] <- c("Pair", names(sig)[i], names(sig)[j])
    rownames(co) <- co$Pair
    co <- co[which(rowProds(as.matrix(co[, 4:5])) > 0), ]
    conc_only_norm[[k]] <- co
    
    n_conc$Dataset1[k] <- names(sig)[i]
    n_conc$Dataset2[k] <- names(sig)[j]
    n_conc$Num_concordant[k] <- nrow(conc_only_norm[[k]])
    
    k <- k + 1
  }
}

n_conc2 <- data.frame(Dataset1 = n_conc$Dataset2, 
                      Dataset2 = n_conc$Dataset1,
                      Num_concordant = n_conc$Num_concordant)
self <- data.frame(Dataset1 = names(all_norm),
                   Dataset2 = names(all_norm),
                   Num_concordant = rep(nrow(sig[[1]]), length(all_norm)))
n_conc_m_norm <- rbind(n_conc, n_conc2, self)

n_conc_m_norm <- mutate(n_conc_m_norm, Dataset2 = factor(Dataset2, levels = names(cdist_norm)),
                   Dataset1 = factor(Dataset1, levels = names(cdist_norm)))

n_conc_m_norm <- t(dcast(n_conc_m_norm, Dataset1 ~ Dataset2)) # Wide format
colnames(n_conc_m_norm) <- n_conc_m_norm[1, ]
n_conc_m_norm <- n_conc_m_norm[-1, ]
n_conc_m_norm[upper.tri(n_conc_m_norm)] <- NA # To remove the upper triangle

n_conc_m_norm <- as.data.frame(n_conc_m_norm) # Needed for IDs when reverting to long format
n_conc_m_norm$Dataset1 <- colnames(n_conc_m_norm)

n_conc_m_norm <- na.omit(melt(n_conc_m_norm, "Dataset1", variable_name = "Dataset2")) # Back to long format
colnames(n_conc_m_norm)[2] <- "Dataset2"
n_conc_m_norm$value <- as.numeric(n_conc_m_norm$value)

n_conc_m_norm <- mutate(n_conc_m_norm, Dataset2 = factor(Dataset2, levels = names(cdist_norm)),
                   Dataset1 = factor(Dataset1, levels = rev(names(cdist_norm))))

# Figure 3
pdf("Fig3.pdf", height = 6, width = 9)
# postscript("Fig3.eps", height = 6, width = 9)
ggplot(data = n_conc_m_norm, aes(x = Dataset1, y = Dataset2, fill = value)) + 
  ggtitle("Number of concordant gene-gene correlations \n between normal breast datasets") + 
  geom_tile(color = "white") +
  guides(fill = guide_colorbar(title = "Number of \ngene pairs")) +
  scale_fill_gradient2(low = "white", high = "red", 
                       midpoint = 100, limit = c(0, 2000)) + 
  geom_text(aes(x = Dataset1, y = Dataset2, label = value), 
            color = "black", size = 4) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5)) 
dev.off()

## Find concordances across multiple datasets
# Between all normal
pairs_only <- lapply(conc_only_norm, function(x) x[, "Pair", drop = FALSE])
c_all_norm <- Reduce(intersect, pairs_only[1:3]) # All bulk/[single-cell] pairs
ref <- conc_only_norm[[1]]
View(ref[which(ref$Pair %in% c_all_norm$Pair), ])

# All single-cell, not bulk
c_sc_norm <- Reduce(intersect, pairs_only[4:6])
c_sc_norm <- c_sc_norm[which(!(c_sc_norm$Pair %in% c_all_norm$Pair)), ]
ref <- conc_only_norm[[6]]
length(which(sign(ref[ref$Pair %in% c_sc_norm, ]$`Mature Luminal`) > 0))

# All bulk only
c_bulk <- sig$`Bulk Normal`[which(!(sig$`Bulk Normal`$Pair %in% 
                                           unlist(pairs_only))), ]
length(which(sign(c_bulk$Correlation) > 0))

# All LP single-cell only
c_lp <- sig$`Luminal Progenitor`[which(!(sig$`Luminal Progenitor`$Pair %in% 
                                           unlist(pairs_only))), ]
length(which(sign(c_lp$Correlation) > 0))

# All ML single-cell only
c_ml <- sig$`Mature Luminal`[which(!(sig$`Mature Luminal`$Pair %in% 
                                           unlist(pairs_only))), ]
length(which(sign(c_ml$Correlation) > 0))

# All basal single-cell only
c_ba <- sig$`Basal Epithelial`[which(!(sig$`Basal Epithelial`$Pair %in% 
                                       unlist(pairs_only))), ]
length(which(sign(c_ba$Correlation) > 0))

save(list = c("cdist_norm", "conc_only_norm", "n_conc_m_norm",
              "c_all_norm", "c_sc_norm", "c_lp", "c_ml", "c_ba"), file = "norm_coexp_data.rda")
