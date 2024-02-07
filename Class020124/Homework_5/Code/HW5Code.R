

# Loading Package and Checking ####
pacman::p_load("Biostrings","msa","tidyr","seqinr")
p_loaded("Biostrings","msa","tidyr","seqinr")

#set WD #####
setwd(~/Documents/GitHub/Bioinformatics_Spring2024/Class020124/Homework_5"")

#Reading and Combining Sequences ####
S_auriculata <- readDNAStringSet("Data/Smilaxauriculata.fasta")
S_bonanox <- readDNAStringSet("Data/Smilaxbonanox.fasta")
S_glauca <- readDNAStringSet("Data/Smilaxglauca.fasta")
S_laurifolia <- readDNAStringSet("Data/Smilaxlaurifolia.fasta")
S_rotundifolia <- readDNAStringSet("Data/Smilaxrotundifolia.fasta")


cSmilax <- c(S_auriculata,S_bonanox,S_glauca,S_laurifolia,S_rotundifolia)
cSmilax

#Doing a MSA w/ incomplete framing ####
cSmilax_msa <- msa(cSmilax,"ClustalW")

#Re-assinging row names and completely printing MSA ####
rownames(cSmilax_msa) <- c("Sglauca", "Slaurifolia", "Srotundifolia", "Sauriculata", "Sbona_nox")
print(cSmilax_msa, show="complete")

#Stats on MSA ####
alphabetFrequency(cSmilax_msa)
GC_Smilax <- letterFrequency(cSmilax, letters="CG")
head(GC_Smilax)

GCP_Smilax <- GC_Smilax/720
GCP_Smilax

#Converting to seqinr#

SmilaxAlnSIR <- msaConvert(cSmilax_msa, type="seqinr::alignment")
SmilaxAlnSIR

Smilax_DA <- dist.alignment(SmilaxAlnSIR, "identity")
Smilax_DA


