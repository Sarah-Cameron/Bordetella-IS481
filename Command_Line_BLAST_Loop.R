#First import the BLAST data file to the environment 
args = commandArgs(trailingOnly = TRUE)
#args[1]="GCF_004635725.1_ASM463572v1_genomic.fna.BLAST.Result.txt"
IS481_BLAST <- read.delim(args[1],stringsAsFactors=F,header=F)

#Assign the columns names to easily identify  
colnames(IS481_BLAST) <- c("IS481_ID", "Strain_ID", "%_Identical", "Align_Length", "No_Mismatches", "No_Gaps", "Start_Align_IS481", "End_Align_IS481", "Start_Align_Strain", "End_Align_Strain", "EValue", "BitScore")

#Now we are going to make an extra column that shows whether the IS is on the opposite strand
IS481_BLAST$Orientation <- IS481_BLAST$Start_Align_Strain > IS481_BLAST$End_Align_Strain

#We then want to separate the truncated insertions into another table and remove them from the IS481_BLAST table
Truncated <- which(IS481_BLAST$Align_Length < 1047)
Truncated_ISE <- IS481_BLAST[Truncated,]
IS481_BLAST <- IS481_BLAST[-Truncated,]

#Now we want to separate the IS481_BLAST results into Fwd and Rev strands 
Rev <- which(IS481_BLAST$Orientation == TRUE)
Rev_Strand <- IS481_BLAST[Rev,]
Fwd_Strand <- IS481_BLAST[-Rev,]

#Now we are going to work out the up and downstream target sequences in the Fwd dataframe for those alignments that aligned the full insertion seq
Fwd_Strand$End_Upstream_Target <- NA
Fwd_Strand$Start_Downstream_Target <- NA 
Fwd_Strand$Start_Upstream_Target <- NA  
Fwd_Strand$End_Downstream_Target <- NA

#Because the start and end positions are the start and end of the IS481 sequence we need to create two columns one -1 from the start and one +1 from the end
Full_1 <- which(Fwd_Strand$Start_Align_IS481 == 1)
Full_1053 <- which(Fwd_Strand$End_Align_IS481 >= 1053)
Fwd_Strand$End_Upstream_Target[Full_1] <- Fwd_Strand$Start_Align_Strain[Full_1] + 4
Fwd_Strand$Start_Downstream_Target[Full_1053] <- Fwd_Strand$End_Align_Strain[Full_1053] - 4
#Now to find the start location of upstream target sequence
Fwd_Strand$Start_Upstream_Target[Full_1] <- Fwd_Strand$"Start_Align_Strain"[Full_1] - 1 
#Now to find the end location of the downstream target
Fwd_Strand$End_Downstream_Target[Full_1053] <- Fwd_Strand$"End_Align_Strain"[Full_1053] + 1

#Now to work out the 6bp for those few that didn't align at the ends because the seq didn't match but were likely to be full length 
Odd_2 <- which(Fwd_Strand$Start_Align_IS481 == 2)
Odd_3 <- which(Fwd_Strand$Start_Align_IS481 == 3)
Odd_4 <- which(Fwd_Strand$Start_Align_IS481 == 4)
Odd_5 <- which(Fwd_Strand$Start_Align_IS481 == 5)
Odd_6 <- which(Fwd_Strand$Start_Align_IS481 == 6)
Odd_7 <- which(Fwd_Strand$Start_Align_IS481 == 7)
Odd_1048 <- which(Fwd_Strand$End_Align_IS481 == 1048)
Odd_1049 <- which(Fwd_Strand$End_Align_IS481 == 1049)
Odd_1050 <- which(Fwd_Strand$End_Align_IS481 == 1050)
Odd_1051 <- which(Fwd_Strand$End_Align_IS481 == 1051)
Odd_1052 <- which(Fwd_Strand$End_Align_IS481 == 1052)
Fwd_Strand$End_Upstream_Target[Odd_2] <- Fwd_Strand$Start_Align_Strain[Odd_2] + 3
Fwd_Strand$End_Upstream_Target[Odd_3] <- Fwd_Strand$Start_Align_Strain[Odd_3] + 2
Fwd_Strand$End_Upstream_Target[Odd_4] <- Fwd_Strand$Start_Align_Strain[Odd_4] + 1
Fwd_Strand$End_Upstream_Target[Odd_5] <- Fwd_Strand$Start_Align_Strain[Odd_5]
Fwd_Strand$End_Upstream_Target[Odd_6] <- Fwd_Strand$Start_Align_Strain[Odd_6] - 1
Fwd_Strand$End_Upstream_Target[Odd_7] <- Fwd_Strand$Start_Align_Strain[Odd_7] - 2
Fwd_Strand$Start_Downstream_Target[Odd_1048] <- Fwd_Strand$End_Align_Strain[Odd_1048] + 1
Fwd_Strand$Start_Downstream_Target[Odd_1049] <- Fwd_Strand$End_Align_Strain[Odd_1049]
Fwd_Strand$Start_Downstream_Target[Odd_1050] <- Fwd_Strand$End_Align_Strain[Odd_1050] - 1
Fwd_Strand$Start_Downstream_Target[Odd_1051] <- Fwd_Strand$End_Align_Strain[Odd_1051] - 2
Fwd_Strand$Start_Downstream_Target[Odd_1052] <- Fwd_Strand$End_Align_Strain[Odd_1052] - 3
#No to find the start location of the upstream targets
Fwd_Strand$Start_Upstream_Target[Odd_2] <- Fwd_Strand$"Start_Align_Strain"[Odd_2] - 2 
Fwd_Strand$Start_Upstream_Target[Odd_3] <- Fwd_Strand$"Start_Align_Strain"[Odd_3] - 3
Fwd_Strand$Start_Upstream_Target[Odd_4] <- Fwd_Strand$"Start_Align_Strain"[Odd_4] - 4
Fwd_Strand$Start_Upstream_Target[Odd_5] <- Fwd_Strand$"Start_Align_Strain"[Odd_5] - 5 
Fwd_Strand$Start_Upstream_Target[Odd_6] <- Fwd_Strand$"Start_Align_Strain"[Odd_6] - 6
Fwd_Strand$Start_Upstream_Target[Odd_7] <- Fwd_Strand$"Start_Align_Strain"[Odd_7] - 7
#Now to find the end location of the downstream target
Fwd_Strand$End_Downstream_Target[Odd_1048] <- Fwd_Strand$"End_Align_Strain"[Odd_1048] + 6
Fwd_Strand$End_Downstream_Target[Odd_1049] <- Fwd_Strand$"End_Align_Strain"[Odd_1049] + 5
Fwd_Strand$End_Downstream_Target[Odd_1050] <- Fwd_Strand$"End_Align_Strain"[Odd_1050] + 4
Fwd_Strand$End_Downstream_Target[Odd_1051] <- Fwd_Strand$"End_Align_Strain"[Odd_1051] + 3
Fwd_Strand$End_Downstream_Target[Odd_1052] <- Fwd_Strand$"End_Align_Strain"[Odd_1052] + 2
#Then we need to extract the sequence between the start and end positions each side 
#We are going to use the package Biostrings to do this 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#                      install.packages("BiocManager")
#BiocManager::install("Biostrings")

#Import Biostrings library and using a FASTA version of the genome sequence (change from .fna to .fasta on server and import to environment)
library(Biostrings)

#args[2]="GCF_004635725.1_ASM463572v1_genomic.fna"
Genome_Seq <- readDNAStringSet(args[2], format="fasta")

#Now create two new columns that will hold the target sequences and leave them black (NA)
Fwd_Strand$Upstream_Target <- NA
Fwd_Strand$Downstream_Target <- NA

#Now we run a loop of this code in order to generate the target sequences for each of the rows.
for(line in c(1:nrow(Fwd_Strand)))
{
  
  Fwd_Strand$Upstream_Target[line]=as.character(Genome_Seq[[1]][Fwd_Strand$Start_Upstream_Target[line]:Fwd_Strand$End_Upstream_Target[line]])
  Fwd_Strand$Downstream_Target[line]=as.character(Genome_Seq[[1]][Fwd_Strand$Start_Downstream_Target[line]:Fwd_Strand$End_Downstream_Target[line]])
}

#Now we need to give each IS element a 'tag' for it's location so this is 20bp up and downstream from the 6bp target sequence
Fwd_Strand$End_Up_20 <- Fwd_Strand$Start_Upstream_Target - 1 
Fwd_Strand$Start_Up_20 <- Fwd_Strand$End_Up_20 - 19
Fwd_Strand$Start_Dwn_20 <- Fwd_Strand$End_Downstream_Target + 1 
Fwd_Strand$End_Dwn_20 <- Fwd_Strand$Start_Dwn_20 + 19
Fwd_Strand$Upstream_Tag <- NA
Fwd_Strand$Downstream_Tag <- NA 

for (line in c(1:nrow(Fwd_Strand)))
{
Fwd_Strand$Upstream_Tag[line] <- as.character(Genome_Seq[[1]][Fwd_Strand$Start_Up_20[line]:Fwd_Strand$End_Up_20[line]])
Fwd_Strand$Downstream_Tag[line] <- as.character(Genome_Seq[[1]][Fwd_Strand$Start_Dwn_20[line]:Fwd_Strand$End_Dwn_20[line]])
}


Rev_Strand$End_Upstream_Target <- NA
Rev_Strand$Start_Downstream_Target <- NA 
Rev_Strand$Start_Upstream_Target <- NA  
Rev_Strand$End_Downstream_Target <- NA

#Now we are going to work out the same for the Rev strand but the opposite orientation
Full_1 <- which(Rev_Strand$Start_Align_IS481 == 1)
Full_1053 <- which(Rev_Strand$End_Align_IS481 >= 1053)
Rev_Strand$End_Upstream_Target[Full_1] <- Rev_Strand$Start_Align_Strain[Full_1] - 4
Rev_Strand$Start_Downstream_Target[Full_1053] <- Rev_Strand$End_Align_Strain[Full_1053] + 4
#Now to find the start location of upstream target sequence
Rev_Strand$Start_Upstream_Target[Full_1] <- Rev_Strand$"Start_Align_Strain"[Full_1] + 1 
#Now to find the end location of the downstream target
Rev_Strand$End_Downstream_Target[Full_1053] <- Rev_Strand$"End_Align_Strain"[Full_1053] - 1

#Now to work out the 6bp for those few that didn't align at the ends because the seq didn't match but were likely to be full length 
Odd_2 <- which(Rev_Strand$Start_Align_IS481 == 2)
Odd_3 <- which(Rev_Strand$Start_Align_IS481 == 3)
Odd_4 <- which(Rev_Strand$Start_Align_IS481 == 4)
Odd_5 <- which(Rev_Strand$Start_Align_IS481 == 5)
Odd_6 <- which(Rev_Strand$Start_Align_IS481 == 6)
Odd_1048 <- which(Rev_Strand$End_Align_IS481 == 1048)
Odd_1049 <- which(Rev_Strand$End_Align_IS481 == 1049)
Odd_1050 <- which(Rev_Strand$End_Align_IS481 == 1050)
Odd_1051 <- which(Rev_Strand$End_Align_IS481 == 1051)
Odd_1052 <- which(Rev_Strand$End_Align_IS481 == 1052)
Rev_Strand$End_Upstream_Target[Odd_2] <- Rev_Strand$Start_Align_Strain[Odd_2] - 3
Rev_Strand$End_Upstream_Target[Odd_3] <- Rev_Strand$Start_Align_Strain[Odd_3] - 2
Rev_Strand$End_Upstream_Target[Odd_4] <- Rev_Strand$Start_Align_Strain[Odd_4] - 1
Rev_Strand$End_Upstream_Target[Odd_5] <- Rev_Strand$Start_Align_Strain[Odd_5] 
Rev_Strand$End_Upstream_Target[Odd_6] <- Rev_Strand$Start_Align_Strain[Odd_6] + 1
Rev_Strand$Start_Downstream_Target[Odd_1048] <- Rev_Strand$End_Align_Strain[Odd_1048] - 1
Rev_Strand$Start_Downstream_Target[Odd_1049] <- Rev_Strand$End_Align_Strain[Odd_1049]
Rev_Strand$Start_Downstream_Target[Odd_1050] <- Rev_Strand$End_Align_Strain[Odd_1050] + 1 
Rev_Strand$Start_Downstream_Target[Odd_1051] <- Rev_Strand$End_Align_Strain[Odd_1051] + 2
Rev_Strand$Start_Downstream_Target[Odd_1052] <- Rev_Strand$End_Align_Strain[Odd_1052] + 3
#Now to find the start location of the upstream targets
Rev_Strand$Start_Upstream_Target[Odd_2] <- Rev_Strand$"Start_Align_Strain"[Odd_2] + 2
Rev_Strand$Start_Upstream_Target[Odd_3] <- Rev_Strand$"Start_Align_Strain"[Odd_3] + 3
Rev_Strand$Start_Upstream_Target[Odd_4] <- Rev_Strand$"Start_Align_Strain"[Odd_4] + 4
Rev_Strand$Start_Upstream_Target[Odd_5] <- Rev_Strand$"Start_Align_Strain"[Odd_5] + 5 
Rev_Strand$Start_Upstream_Target[Odd_6] <- Rev_Strand$"Start_Align_Strain"[Odd_6] + 6
#Now to find the end location of the downstream target
Rev_Strand$End_Downstream_Target[Odd_1048] <- Rev_Strand$"End_Align_Strain"[Odd_1048] - 6
Rev_Strand$End_Downstream_Target[Odd_1049] <- Rev_Strand$"End_Align_Strain"[Odd_1049] - 5
Rev_Strand$End_Downstream_Target[Odd_1050] <- Rev_Strand$"End_Align_Strain"[Odd_1050] - 4
Rev_Strand$End_Downstream_Target[Odd_1051] <- Rev_Strand$"End_Align_Strain"[Odd_1051] - 3
Rev_Strand$End_Downstream_Target[Odd_1052] <- Rev_Strand$"End_Align_Strain"[Odd_1052] - 2

Rev_Strand$Upstream_Target <- NA 
Rev_Strand$Downstream_Target <- NA

#Now we run a loop of this code in order to generate the target sequences for each of the rows.
for(line in c(1:nrow(Rev_Strand)))
{
  Upstream_Target <- DNAStringSet(Genome_Seq[[1]][Rev_Strand$End_Upstream_Target[line]:Rev_Strand$Start_Upstream_Target[line]])
  Downstream_Target <- DNAStringSet(Genome_Seq[[1]][Rev_Strand$End_Downstream_Target[line]:Rev_Strand$Start_Downstream_Target[line]])
  Rev_Strand$Upstream_Target[line] <- reverseComplement(Upstream_Target)
  Rev_Strand$Downstream_Target[line] <- reverseComplement(Downstream_Target)
}

#Make the 20bp tags either side
Rev_Strand$End_Up_20 <- Rev_Strand$Start_Upstream_Target + 1 
Rev_Strand$Start_Up_20 <- Rev_Strand$End_Up_20 + 19 
Rev_Strand$Start_Dwn_20 <- Rev_Strand$End_Downstream_Target - 1 
Rev_Strand$End_Dwn_20 <- Rev_Strand$Start_Dwn_20 - 19
Rev_Strand$Upstream_Tag <- NA 
Rev_Strand$Downstream_Tag <- NA

for(line in c(1:nrow(Rev_Strand)))
{
  Upstream_Tag <- DNAStringSet(Genome_Seq[[1]][Rev_Strand$End_Up_20[line]:Rev_Strand$Start_Up_20[line]])
  Downstream_Tag <- DNAStringSet(Genome_Seq[[1]][Rev_Strand$End_Dwn_20[line]:Rev_Strand$Start_Dwn_20[line]])
  Rev_Strand$Upstream_Tag[line] <- reverseComplement(Upstream_Tag)
  Rev_Strand$Downstream_Tag[line] <- reverseComplement(Downstream_Tag)
}

#Now we want to merge the two tables back together 
Target_Sequences <- merge(Fwd_Strand, Rev_Strand, all = TRUE)

#Save full insertions and truncated insertions into csvs 
write.csv(Truncated_ISE, paste(args[1],"_Truncated_Insertions.csv",sep=""))
write.csv(Target_Sequences, paste(args[1],"_Full_Insertions.csv",sep=""))
