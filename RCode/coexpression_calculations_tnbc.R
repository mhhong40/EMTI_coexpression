library(tidyverse)
library(here)
library(hdWGCNA)

setwd(here::here("~/Documents/tnbc_coexpression"))

# Load TNBC datasets
# Scale the single-cell datasets
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

# Scale datasets, subset data to final gene list, alphabetize gene names
tsc <- list(e1, e2, e3, c1, c2, c3)
gs <- readRDS("gs_final.rds")
tsc <- lapply(tsc, function(x) {
  
  x <- x@assays$RNA@scale.data
  x <- x[which(rownames(x) %in% gs), ]
  x <- x[order(rownames(x)), ]
})
tnbc_bulk <- tnbc_bulk[which(rownames(tnbc_bulk) %in% gs), ]
tnbc_bulk <- tnbc_bulk[order(rownames(tnbc_bulk)), ]

# Correlation calculation 
library(rstatix)
library(reshape2)
all_tnbc <- append(list(tnbc_bulk), tsc)
names(all_tnbc) <- c("Bulk TNBC", "Epithelial 1", "Epithelial 2", "Epithelial 3",
                     "Cycling Epithelial 1", "Cycling Epithelial 2", "Cycling Epithelial 3")
ct <- cdist <- list()
for(i in 1:length(all_tnbc)){
  d <- as.data.frame(t(all_tnbc[[i]]))
  ct[[i]] <- d %>% cor(method = "pearson")
  cdist[[i]] <- melt(ct[[i]])
  colnames(cdist[[i]]) <- c("Gene1", "Gene2", "Correlation")
  cdist[[i]] <- cdist[[i]][which(cdist[[i]]$Gene1 != cdist[[i]]$Gene2), ] # Remove self-correlations
}
names(ct) <- names(cdist) <- names(all_tnbc)

# Retrieve only top 5% highly correlated gene pairs
sig <- list()
for(i in 1:length(all_tnbc)) {
  
  # Calculate cutoffs
  lower <- quantile(cdist[[i]]$Correlation, probs = 0.025)
  upper <- quantile(cdist[[i]]$Correlation, probs = 0.975)
  
  # Retrieve significant correlations
  sig[[i]] <- cdist[[i]][cdist[[i]]$Correlation < lower | cdist[[i]]$Correlation > upper, ]
}
names(sig) <- names(all_tnbc)

# Alphabetize gene pair names
for(i in 1:length(all_tnbc)){
  
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
library(matrixStats)
n_pairs <- length(all_tnbc) * (length(all_tnbc) - 1) * 1/2
conc_only <- list()
n_conc <- data.frame(matrix(nrow = n_pairs, ncol = 3))
colnames(n_conc) <- c("Dataset1", "Dataset2", "Num_concordant")
k <- 1 # Counter for the conc_only list
for(i in 1:(length(all_tnbc) - 1)) {
  
  for(j in (i + 1):length(all_tnbc)) {
    
    # Retrieve significant correlations for each set
    d1 <- sig[[i]]
    d2 <- sig[[j]]
    
    # Gene pairs significantly correlated in both sets
    d1 <- d1[which(d1$Pair %in% intersect(d1$Pair, d2$Pair)), ]
    d2 <- d2[which(d2$Pair %in% intersect(d1$Pair, d2$Pair)), ]
    
    # Calculate number of concordant pairs between datasets
    conc_only[[k]] <- cbind(d1[, 1:2], d1$Pair, d1$Correlation, d2$Correlation)
    
    co <- conc_only[[k]]
    colnames(co)[3:5] <- c("Pair", names(sig)[i], names(sig)[j])
    rownames(co) <- co$Pair
    co <- co[which(rowProds(as.matrix(co[, 4:5])) > 0), ]
    conc_only[[k]] <- co
    
    n_conc$Dataset1[k] <- names(sig)[i]
    n_conc$Dataset2[k] <- names(sig)[j]
    n_conc$Num_concordant[k] <- nrow(conc_only[[k]])
    
    k <- k + 1
  }
}
  
n_conc2 <- data.frame(Dataset1 = n_conc$Dataset2, 
                      Dataset2 = n_conc$Dataset1,
                      Num_concordant = n_conc$Num_concordant)
self <- data.frame(Dataset1 = names(all_tnbc),
                   Dataset2 = names(all_tnbc),
                   Num_concordant = rep(nrow(sig[[1]]), length(all_tnbc)))
n_conc_m <- rbind(n_conc, n_conc2, self)

n_conc_m <- mutate(n_conc_m, Dataset2 = factor(Dataset2, levels = names(cdist)),
                   Dataset1 = factor(Dataset1, levels = names(cdist)))

n_conc_m <- t(dcast(n_conc_m, Dataset1 ~ Dataset2)) # Wide format
colnames(n_conc_m) <- n_conc_m[1, ]
n_conc_m <- n_conc_m[-1, ]
n_conc_m[upper.tri(n_conc_m)] <- NA # To remove the upper triangle

n_conc_m <- as.data.frame(n_conc_m) # Needed for IDs when reverting to long format
n_conc_m$Dataset1 <- colnames(n_conc_m)

n_conc_m <- na.omit(melt(n_conc_m, "Dataset1", variable_name = "Dataset2")) # Back to long format
colnames(n_conc_m)[2] <- "Dataset2"
n_conc_m$value <- as.numeric(n_conc_m$value)

n_conc_m <- mutate(n_conc_m, Dataset2 = factor(Dataset2, levels = names(cdist)),
                   Dataset1 = factor(Dataset1, levels = rev(names(cdist))))

# Figure 1
library(ggplot2)
pdf("Fig1.pdf", height = 6, width = 9)
# postscript("Fig1.eps", height = 6, width = 9)
ggplot(data = n_conc_m, aes(x = Dataset1, y = Dataset2, fill = value)) + 
  ggtitle("Number of concordant gene-gene co-expressions \n between TNBC datasets") + 
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
# Between all TNBC
pairs_only <- lapply(conc_only, function(x) x[, "Pair", drop = FALSE])
c_all <- Reduce(intersect, pairs_only[1:6]) # All bulk/[single-cell] pairs

# Bulk only
c_bulk <- sig$`Bulk TNBC`[which(!(sig$`Bulk TNBC`$Pair %in% 
                                      unlist(pairs_only))), ]
length(which(sign(c_bulk$Correlation) > 0))

# All single-cell, not bulk
c_sc <- Reduce(intersect, pairs_only[7:21])
c_sc <- c_sc[which(!(c_sc$Pair %in% c_all$Pair)), ]
ref <- conc_only[[7]]
View(ref[which(ref$Pair %in% c_sc), ])

# All non-cycling single-cell only
c_epi <- Reduce(intersect, pairs_only[c(7, 8, 12)])
c_epi <- c_epi[which(!(c_epi$Pair %in% c(c_all$Pair, c_sc, c_cyc))), ]
View(ref[which(ref$Pair %in% c_epi), ])

# All cycling single-cell only
c_cyc <- Reduce(intersect, pairs_only[19:21])
c_cyc <- c_cyc[which(!(c_cyc$Pair %in% c(c_all$Pair, c_sc))), ]
ref <- conc_only[[21]]
View(ref[which(ref$Pair %in% c_cyc), ])

save(list = c("cdist", "conc_only", "n_conc_m",
              "c_all", "c_sc", "c_cyc", "c_epi"), file = "tnbc_coexp_data.rda")
