library(tidyverse)
library(here)
library(readxl)
library(UCSCXenaTools)

# TCGA primary tumor IDs that are TNBC
clinicalTCGA <- XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterCohorts = "TCGA Breast Cancer") %>%
  XenaFilter(filterDatasets = "BRCA_clinicalMatrix")

XenaQuery(clinicalTCGA) %>%
  XenaDownload(destdir = getwd()) # Download file called BRCA_clinicalMatrix.txt

# Load in the manually filtered spreadsheet (TNBC_clinicalMatrix.tsv)
TCGAclinical <- read_tsv("TNBC_clinicalMatrix.tsv")
tnbcs <- TCGAclinical$sampleID

f2TCGA <- TCGAclinical[TCGAclinical$sampleID %in% tnbcs, ]
primtnbcs <- tnbcs[-which(grepl("-11", tnbcs) | grepl("-06", tnbcs))] # Only TNBC primary tumors
saveRDS(primtnbcs, "primtnbcs.rds")