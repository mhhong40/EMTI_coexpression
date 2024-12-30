library(tidyverse)
library(here)
library(UCSCXenaTools)

# GTEx phenotypic data
gtexPheno <- XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TCGATargetGTEx_phenotype")

XenaQuery(gtexPheno) %>% 
  XenaDownload(destdir = getwd())

f1GTEx <- fread("TcgaTargetGTEX_phenotype.txt.gz")
names(f1GTEx) <- gsub("\\_", "", names(f1GTEx))

paraStudy <- "GTEX"
paraPrimarySiteGTEx <- "Breast"
paraPrimaryTissueGTEx <- "^Breast"

# Subset GTEx breast tissue RNA-seq data
f2GTEx <- subset(f1GTEx, study == paraStudy &
                  primarysite == paraPrimarySiteGTEx &
                  grepl(paraPrimaryTissueGTEx, f1GTEx$"primary disease or tissue"))
saveRDS(f2GTEx, file = "f2GTEx.rds")