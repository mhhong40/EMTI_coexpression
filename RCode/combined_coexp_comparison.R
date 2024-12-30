library(tidyverse)
library(here)

setwd(here::here("~/Documents/tnbc_coexpression"))

load("tnbc_coexp_data.rda")
load("norm_coexp_data.rda")

## Normal vs TNBC
# Merged data objects
cdist_all <- c(cdist, cdist_norm)

# Retrieve only top 5% highly correlated gene pairs
sig <- list()
for(i in 1:length(cdist_all)) {
  
  # Calculate cutoffs
  lower <- quantile(cdist_all[[i]]$Correlation, probs = 0.025)
  upper <- quantile(cdist_all[[i]]$Correlation, probs = 0.975)
  
  # Retrieve significant correlations
  sig[[i]] <- cdist_all[[i]][cdist_all[[i]]$Correlation < lower | cdist_all[[i]]$Correlation > upper, ]
}
names(sig) <- names(cdist_all)

# Alphabetize gene pair names
for(i in 1:length(cdist_all)){
  
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

# Concordances
n_pairs <- length(cdist_all) * (length(cdist_all) - 1) * 1/2
conc_only <- list()
n_conc <- data.frame(matrix(nrow = n_pairs, ncol = 3))
colnames(n_conc) <- c("Dataset1", "Dataset2", "Num_concordant")
k <- 1 # Counter for the conc_only list

library(matrixStats)
for(i in 1:(length(cdist_all) - 1)) {
  
  for(j in (i + 1):length(cdist_all)) {
    
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
self <- data.frame(Dataset1 = names(cdist_all),
                   Dataset2 = names(cdist_all),
                   Num_concordant = rep(nrow(sig[[1]]), length(cdist_all)))
n_conc_f <- rbind(n_conc, n_conc2, self)

n_conc_f <- distinct(n_conc_f) # Remove duplicate self-correlations
n_conc_f <- mutate(n_conc_f, Dataset1 = factor(Dataset1, levels = names(cdist_all)),
                 Dataset2 = factor(Dataset2, levels = names(cdist_all)))

library(reshape2)
n_conc_f <- t(dcast(n_conc_f, Dataset1 ~ Dataset2)) # Wide format
colnames(n_conc_f) <- n_conc_f[1, ]
n_conc_f <- n_conc_f[-1, ]
n_conc_f[upper.tri(n_conc_f)] <- NA # To remove the upper triangle

n_conc_f <- as.data.frame(n_conc_f) # Needed for IDs when reverting to long format
n_conc_f$Dataset1 <- colnames(n_conc_f)

n_conc_f <- na.omit(melt(n_conc_f, "Dataset1", variable_name = "Dataset2")) # Back to long format
colnames(n_conc_f)[2] <- "Dataset2"
n_conc_f$value <- as.numeric(n_conc_f$value)

n_conc_f <- mutate(n_conc_f, Dataset2 = factor(Dataset2, levels = names(cdist_all)),
                   Dataset1 = factor(Dataset1, levels = rev(names(cdist_all))))

# Figure 5
pdf("Fig5.pdf", height = 6, width = 9)
# postscript("Fig5.eps", height = 6, width = 9)
ggplot(data = n_conc_f, aes(x = Dataset1, y = Dataset2, fill = value)) + 
  ggtitle("Number of concordant gene-gene correlations \nbetween normal and TNBC breast datasets") + 
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

### Concordances between various mixed groupings of datasets
# All 11 datasets
pairs_only <- lapply(conc_only, function(x) x[, "Pair", drop = FALSE])
c_all_all <- intersect(c_all$Pair, c_all_norm$Pair) # CLDN4/ELF3 is positive in all

# All single-cell datasets
c_sc_all <- Reduce(intersect, pairs_only[c(11:15, 17:19)]) # All epithelial 1/[rest TNBC single-cell], all epithelial 1/[normal single-cell]
# Equal to 0 

## All single-cell TNBC with each normal epithelial lineage
# LP
c_sc_lp <- sig$`Luminal Progenitor`[which(sig$`Luminal Progenitor`$Pair %in% c_sc), ] 
ref <- sig$`Epithelial 1`[which(sig$`Epithelial 1`$Pair %in% c_sc_lp$Pair), ]
c_sc_lp <- c_sc_lp[which(sign(c_sc_lp$Correlation) == sign(ref$Correlation)), ] # Check that the signs match

# ML
c_sc_ml <- sig$`Mature Luminal`[which(sig$`Mature Luminal`$Pair %in% c_sc), ]
ref <- sig$`Epithelial 1`[which(sig$`Epithelial 1`$Pair %in% c_sc_ml$Pair), ]
c_sc_ml <- c_sc_ml[which(sign(c_sc_ml$Correlation) == sign(ref$Correlation)), ] 

# Basal
c_sc_ba <- sig$`Basal Epithelial`[which(sig$`Basal Epithelial`$Pair %in% c_sc), ]
ref <- sig$`Epithelial 1`[which(sig$`Epithelial 1`$Pair %in% c_sc_ba$Pair), ]
c_sc_ba <- c_sc_ba[which(sign(c_sc_ba$Correlation) == sign(ref$Correlation)), ]


## All non-cycling epithelial only with each normal epithelial lineage
# LP
c_epi_lp <- sig$`Luminal Progenitor`[which(sig$`Luminal Progenitor`$Pair %in% c_epi), ] 
ref <- sig$`Epithelial 1`[which(sig$`Epithelial 1`$Pair %in% c_epi_lp$Pair), ]
c_epi_lp <- c_epi_lp[which(sign(c_epi_lp$Correlation) == sign(ref$Correlation)), ] # Check that the signs match

# ML
c_epi_ml <- sig$`Mature Luminal`[which(sig$`Mature Luminal`$Pair %in% c_epi), ]
ref <- sig$`Epithelial 1`[which(sig$`Epithelial 1`$Pair %in% c_epi_ml$Pair), ]
c_epi_ml <- c_epi_ml[which(sign(c_epi_ml$Correlation) == sign(ref$Correlation)), ] 

# Basal
c_epi_ba <- sig$`Basal Epithelial`[which(sig$`Basal Epithelial`$Pair %in% c_epi), ]
ref <- sig$`Epithelial 1`[which(sig$`Epithelial 1`$Pair %in% c_epi_ba$Pair), ]
c_epi_ba <- c_epi_ba[which(sign(c_epi_ba$Correlation) == sign(ref$Correlation)), ]


## All cycling epithelial only with each normal epithelial lineage
c_cyc_lp <- sig$`Luminal Progenitor`[which(sig$`Luminal Progenitor`$Pair %in% c_cyc), ] 
ref <- sig$`Cycling Epithelial 1`[which(sig$`Cycling Epithelial 1`$Pair %in% c_cyc_lp$Pair), ]
c_cyc_lp <- c_cyc_lp[which(sign(c_cyc_lp$Correlation) == sign(ref$Correlation)), ] # Check that the signs match

# ML
c_cyc_ml <- sig$`Mature Luminal`[which(sig$`Mature Luminal`$Pair %in% c_cyc), ]
ref <- sig$`Cycling Epithelial 1`[which(sig$`Cycling Epithelial 1`$Pair %in% c_cyc_ml$Pair), ]
c_cyc_ml <- c_cyc_ml[which(sign(c_cyc_ml$Correlation) == sign(ref$Correlation)), ] 

# Basal
c_cyc_ba <- sig$`Basal Epithelial`[which(sig$`Basal Epithelial`$Pair %in% c_cyc), ]
ref <- sig$`Cycling Epithelial 1`[which(sig$`Cycling Epithelial 1`$Pair %in% c_cyc_ba$Pair), ]
c_cyc_ba <- c_cyc_ba[which(sign(c_cyc_ba$Correlation) == sign(ref$Correlation)), ] 
                     
## Both bulk datasets only
c_bulk_all <- pairs_only[[7]]
c_bulk_all <- c_bulk_all[which(!(c_bulk_all$Pair %in% unlist(pairs_only[-7]))), ] # Pairs that aren't concordant anywhere else
ref <- conc_only[[7]][which(conc_only[[7]]$Pair %in% c_bulk_all), ]
length(which(sign(ref$`Bulk TNBC`) > 0))

save(list = c("cdist_all", "conc_only", "pairs_only", 
              "c_all_all", "c_sc_all", "c_sc_lp", "c_sc_ml", "c_sc_ba",
              "c_epi_lp", "c_epi_ml", "c_epi_ba",
              "c_cyc_lp", "c_cyc_ml", "c_cyc_ba", "c_bulk_all"),
     file = "combined_coexp_data.rda")
