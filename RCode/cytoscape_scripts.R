library(RCy3)
library(igraph)

cytoscapePing()
cytoscapeVersionInfo()

# Load genes
network_gs <- read_tsv("EMT network genes - Sheet1.tsv")
network_gs <- network_gs$gene
network_gs <- append(network_gs, "MMP14")

## TNBC networks
load("tnbc_coexp_data.rda")
n_tnbc <- length(cdist)

tnbc_ngs <- list()
for(i in 1:n_tnbc){
  
  cors <- cdist[[i]]
  
  # Alphabetize gene pair names
  cors[, c("Gene1", "Gene2")] <- lapply(cors[, c("Gene1", "Gene2")], as.character)
  changeInd <- !!(cors$Gene1 > cors$Gene2)
  cors[changeInd, c("Gene1", "Gene2")] <- cors[changeInd, c("Gene2", "Gene1")]
  
  cors$Pair <- NULL
  cors$Pair <- paste(cors$Gene1, cors$Gene2)
  cors <- cors[order(cors$Pair), ]
  cors <- distinct(cors)
  
  # Calculate cutoffs
  lower <- quantile(cors$Correlation, probs = 0.025)
  upper <- quantile(cors$Correlation, probs = 0.975)
  
  # Mark significant correlations
  cors$Significant <- NULL
  cors$Significant <- cors$Correlation > upper | cors$Correlation < lower
  
  tnbc_ngs[[i]] <- cors[cors$Gene1 %in% network_gs & cors$Gene2 %in% network_gs, ]
  tnbc_ngs[[i]] <- tnbc_ngs[[i]][order(tnbc_ngs[[i]]$Pair),]

}
names(tnbc_ngs) <- names(cdist)
saveRDS(tnbc_ngs, "tnbc_ngs.rds")

tg <- list()
ref <- tnbc_ngs[[1]]
for(i in 1:n_tnbc){
  
  t <- tnbc_ngs[[i]]
  md <- network_gs[which(network_gs %in% ref$Gene1 | network_gs %in% ref$Gene2)]
  cors <- data.frame(from = t$Gene1,
                     to = t$Gene2,
                     correlation = t$Correlation,
                     significant = t$Significant)
  tg[[i]] <- graph_from_data_frame(cors, directed = FALSE, vertices = md)
}
names(tg) <- names(tnbc_ngs)
createNetworkFromIgraph(tg[[1]], "Bulk TNBC")
createNetworkFromIgraph(tg[[2]], "Epithelial 1")
createNetworkFromIgraph(tg[[3]], "Epithelial 2")
createNetworkFromIgraph(tg[[4]], "Epithelial 3")
createNetworkFromIgraph(tg[[5]], "Cycling Epithelial 1")
createNetworkFromIgraph(tg[[6]], "Cycling Epithelial 2")
createNetworkFromIgraph(tg[[7]], "Cycling Epithelial 3")


## Normal breast networks
load("norm_coexp_data.rda")
n_norm <- length(cdist_norm)
norm_ngs <- list()
for(i in 1:n_norm){
  
  cors <- cdist_norm[[i]]
  
  # Alphabetize gene pair names
  cors[, c("Gene1", "Gene2")] <- lapply(cors[, c("Gene1", "Gene2")], as.character)
  changeInd <- !!(cors$Gene1 > cors$Gene2)
  cors[changeInd, c("Gene1", "Gene2")] <- cors[changeInd, c("Gene2", "Gene1")]
  
  cors$Pair <- NULL
  cors$Pair <- paste(cors$Gene1, cors$Gene2)
  cors <- cors[order(cors$Pair), ]
  cors <- distinct(cors)
  
  # Calculate cutoffs
  lower <- quantile(cors$Correlation, probs = 0.025)
  upper <- quantile(cors$Correlation, probs = 0.975)
  
  # Mark significant correlations
  cors$Significant <- NULL
  cors$Significant <- cors$Correlation > upper | cors$Correlation < lower
  
  norm_ngs[[i]] <- cors[cors$Gene1 %in% network_gs & cors$Gene2 %in% network_gs, ]
  norm_ngs[[i]] <- norm_ngs[[i]][order(norm_ngs[[i]]$Pair),]
  
}
names(norm_ngs) <- names(cdist_norm)
saveRDS(norm_ngs, "norm_ngs.rds")

ng <- list()
ref <- norm_ngs[[1]]
for(i in 1:n_norm){
  
  n <- norm_ngs[[i]]
  md <- network_gs[which(network_gs %in% ref$Gene1 | network_gs %in% ref$Gene2)]
  cors <- data.frame(from = n$Gene1,
                     to = n$Gene2,
                     correlation = n$Correlation,
                     significant = n$Significant)
  ng[[i]] <- graph_from_data_frame(cors, directed = FALSE, vertices = md)
}
names(ng) <- names(norm_ngs)
createNetworkFromIgraph(ng[[1]], "Bulk Normal")
createNetworkFromIgraph(ng[[2]], "Luminal Progenitor")
createNetworkFromIgraph(ng[[3]], "Mature Luminal")
createNetworkFromIgraph(ng[[4]], "Basal")
