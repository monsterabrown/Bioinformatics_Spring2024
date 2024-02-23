#START OF MIDTERM ####

#PACKAGE MANAGMENT####
#First I enable the pacman package  using the library() function, as I use pacman to load and check library usage for packages I use. I then load in all the packages that we have used so far using the function p_load, and check that they are loaded using the function p_loaded.

library(pacman)
pacman::p_load("Biostrings","msa","tidyr","seqinr","UniprotR","protti")
p_loaded("Biostrings","msa","tidyr","seqinr","UniprotR","protti")


#READING IN DATA ####
#To read in the data I use the Biostrings function readDNAStringSet which allows the sequences.fasta file to be read. To allow the function to know that the sequence file is in fasta format, you use the format="fasta" argument within the function. I put this onto a variable called "PatientData", so I could use it in downstream processes. I type the variable out after to test that the function worked properly. 

PatientData <- readDNAStringSet(file="Data/sequences.fasta", format="fasta")
PatientData

#DOING THE MSA ####
#I will use the msa function to align the sequences into an alignment matrix. I will use the argument method = "Muscle" to use Muscle algorithms for my multiple sequence alignment. I will assign this function to PD_Aligned to use it in later processes. I then write the variable on the next line to double check that the function in the previous line worked correctly. Lastly, I printed the whole alignment using the print function print with the argument show = "complete" to make sure that alignment wasn't faulty or incorrect. 

PD_Aligned <- msa(PatientData,method="Muscle")
PD_Aligned
print(PD_Aligned,show="complete")

#DISTANCE MATRIX ####
#I first converted the Patient Data MSA into a seqinr format for compatibility within the other seqinr package functions. I then wrote out the variable to make sure the function worked. Next, I performed a pairwise distance alignment by writing the dist.alignment function to the PDA_Dist variable. Then I wrote the PDA_Dist variable on the next line so you could view the pairwise distance matrix when the line was run.

PD_Aligned_SIR <- msaConvert(PD_Aligned, type="seqinr::alignment")
PD_Aligned_SIR

PDA_Dist <- dist.alignment(PD_Aligned_SIR, "identity")
PDA_Dist

#EXPORTING MSA TO FASTA ####
# I then wrote the Patient Data MSA to a variable through the DNAStringSet function to get it back into a format I could export into a fasta. Then using the Biostrings package's function writeXStringSet, I exported the PD_MSA variable to a fasta file format and then moved it into my output folder.

PD_MSA <-  DNAStringSet(PD_Aligned)
Biostrings::writeXStringSet(PD_MSA, "PD_MSA.fasta")

#RESULTS FROM BLASTING ####
#I ran Homo_sapiens_1 as it was one of the sequences within the 20 that were in large similar to the rest. The result from the BLAST was that this gene was an hbb gene for beta globin. The ascension number for the top result was LC121775.1, with an 100% match to our Homo_sapiens_1 sequence.

#CONVERTING Homo_sapiens_6 TO PROTEIN AND THEN WRITING IT TO FASTA ####
#Using BBEdit, I extracted the Homo_sapiens_6 sequence into its own fasta file and read it into the variable Hs6 using the Biostrings package function readDNAStringSet. Then I converted it to an amino acid sequence using Biostrings translate function. The argument genetic.code = GENETIC_CODE makes it the default genetic code in the function arguments. The no.init.codon = TRUE argument specifies that there isn't a specified initiation codon (AUG) within the protein sequence. Then I wrote the amino acid sequence to a fasta file with the write.fasta function, and put the file into the output folder.

Hs6 <- readDNAStringSet(file="Data/Homo_sapiens_6.fasta", format="fasta")
Hs6

Hs6_AA <- Biostrings::translate(Hs6, genetic.code = GENETIC_CODE,no.init.codon = TRUE)
Hs6_AA

write.fasta(names = "Hs66_AA",sequences=Hs6_AA,file.out="Hs6_AA.fasta" )

#GETTING GENE PATHOLOGY ####
#I first assigned the accession number that was achieved through BLASTing the Homo_sapiens_6 amino acid sequence in UniProt to the varible Hs6_A. Then I assigned the GetPathology_Biotech function to the variable hbb_P. This function allows the scraping of pathology data from UniProt. Then I printed the varible hbb_P to see the results, in which there was none. I don't think this function is very good at what it does because this type of result has happened multiple times even though the accession clearly has pathology related to it on the UniProt page.

Hs6_A <- "A0A0J9YWK4"
hbb_P <-  GetPathology_Biotech(Hs6_A)
print(hbb_P)

#DOING AN MSA BETWEEN Homo_sapiens_6 AND THE HBB GENE ####
#So the reason I decided to do an alignment between the Homo_sapiens_6 entry and the hbb gene was to check for protein mutations so I could track pathology. This alignment was ultimately unhelpful though in solving what kind of mutation this individual had, and resulting diseases or variants that resulted in such mutation. The first part of this section is following the same process in "READING IN DATA" except it then uses the combine function to create a data frame of both the hbb gene and the Homo_sapiens_6 amino acid sequence. This combine function is assigned to the variable HG_Hs6_AA. Then I did the multiple sequence alignment for the two sequences which follows the exact same process as outlined in the section "DOING THE MSA".

HG <- readAAStringSet(file="Data/hbbgene.fasta", format="fasta")
HG
HG_Hs6_AA <-  c(HG,Hs6_AA)

HGHs6_Aligned <- msa(HG_Hs6_AA,method="Muscle")
HGHs6_Aligned
print(HGHs6_Aligned, show="complete")

#END OF MIDTERM ####