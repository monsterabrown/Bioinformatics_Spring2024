

#Installing and Loading Packages ####

# BiocManager::install("UniprotR")
# BiocManager::install("protti")

pacman::p_load("Biostrings","msa","tidyr","seqinr","UniprotR","protti")
p_loaded("Biostrings","msa","tidyr","seqinr","UniprotR","protti")

#Read in Data ####

setwd("/Class020824/Homework6")
# check out the 'here' package. Helps with directory syncing across machines

UM1 <- read.table("Data/UniprotM1.csv")
UM2 <- read.table("Data/UniprotM2.csv")
UM3 <- read.table("Data/UniprotM3.csv")
UM4 <- read.table("Data/UniprotM4.csv")
UM5 <- read.table("Data/UniprotM5.csv")
UM <- c(UM1,UM2,UM3,UM4,UM5)
UM
print(UM)

#Get Gene Ontology ####
UMGO <- GetProteinGOInfo(UM)

# output plot to a new directory
PlotGoInfo(UMGO, directorypath = "PlotGoInfo")

#Get Gene Pathology ####
UMP <-  GetPathology_Biotech(UM)
print(UMP)

#Structural Bioinformatics using protti ####
UMProvided <- c("1ZMR","2HWG")
fetch_uniprot(UM)
fetch_pdb(UMProvided)
fetch_alphafold_prediction(UM)












