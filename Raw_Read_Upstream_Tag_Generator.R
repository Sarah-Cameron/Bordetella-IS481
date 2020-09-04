#First import the BLAST data file to the environment 
args = commandArgs(trailingOnly = TRUE)
#args[1]="ERR212378.fasta.BLAST.Result.First.txt"
IS481_BLAST <- read.delim(args[1],stringsAsFactors=F, head=F)

#Assign the column names 
colnames(IS481_BLAST) <- c("IS481_ID", "Strain_ID", "%_Identical", "Align_Length", "No_Mismatches", "No_Gaps", "Start_Align_IS481", "End_Align_IS481", "Start_Align_Strain", "End_Align_Strain", "EValue", "BitScore", "Strand")

#Remove the matches that haven't aligned at the start 
Truncated <- which(IS481_BLAST$Start_Align_IS481 > 7)
Truncated_ISE <- IS481_BLAST[Truncated,]
IS481_BLAST <- IS481_BLAST[-Truncated,]

#Now we want to separate the IS481_BLAST results into Fwd and Rev strands 
Rev <- which(IS481_BLAST$Strand == "minus")
Rev_Strand <- IS481_BLAST[Rev,]
Fwd_Strand <- IS481_BLAST[-Rev,]

#Assign the ends 
Full_1 <- which(Fwd_Strand$Start_Align_IS481 == 1)
Odd_2 <- which(Fwd_Strand$Start_Align_IS481 == 2)
Odd_3 <- which(Fwd_Strand$Start_Align_IS481 == 3)
Odd_4 <- which(Fwd_Strand$Start_Align_IS481 == 4)
Odd_5 <- which(Fwd_Strand$Start_Align_IS481 == 5)
Odd_6 <- which(Fwd_Strand$Start_Align_IS481 == 6)
Odd_7 <- which(Fwd_Strand$Start_Align_IS481 == 7)

#Add some columns to work out the upstream tag 
Fwd_Strand$Start_Up_Tag <- NA
Fwd_Strand$End_Up_Tag <- NA 
Fwd_Strand$Upstream_Tag <- NA 

#Calculate the end of the upstream tag 
Fwd_Strand$End_Up_Tag[Full_1] <- Fwd_Strand$Start_Align_Strain[Full_1] - 2 
Fwd_Strand$End_Up_Tag[Odd_2] <- Fwd_Strand$Start_Align_Strain[Odd_2] - 3
Fwd_Strand$End_Up_Tag[Odd_3] <- Fwd_Strand$Start_Align_Strain[Odd_3] - 4 
Fwd_Strand$End_Up_Tag[Odd_4] <- Fwd_Strand$Start_Align_Strain[Odd_4] - 5 
Fwd_Strand$End_Up_Tag[Odd_5] <- Fwd_Strand$Start_Align_Strain[Odd_5] - 6 
Fwd_Strand$End_Up_Tag[Odd_6] <- Fwd_Strand$Start_Align_Strain[Odd_6] - 7
Fwd_Strand$End_Up_Tag[Odd_7] <- Fwd_Strand$Start_Align_Strain[Odd_7] - 8

#Calculate the start of the upstream tag
Fwd_Strand$Start_Up_Tag[Full_1] <- Fwd_Strand$End_Up_Tag[Full_1] - 19 
Fwd_Strand$Start_Up_Tag[Odd_2] <- Fwd_Strand$End_Up_Tag[Odd_2] - 19
Fwd_Strand$Start_Up_Tag[Odd_3] <- Fwd_Strand$End_Up_Tag[Odd_3] - 19 
Fwd_Strand$Start_Up_Tag[Odd_4] <- Fwd_Strand$End_Up_Tag[Odd_4] - 19
Fwd_Strand$Start_Up_Tag[Odd_5] <- Fwd_Strand$End_Up_Tag[Odd_5] - 19
Fwd_Strand$Start_Up_Tag[Odd_6] <- Fwd_Strand$End_Up_Tag[Odd_6] - 19
Fwd_Strand$Start_Up_Tag[Odd_7] <- Fwd_Strand$End_Up_Tag[Odd_7] - 19

#Now we need to add a control step where the script will not try and work out any tags for alignments that have ended at the end of the read (aka have minus numbers in the positon columns)
End_of_Reads <- which(Fwd_Strand$Start_Up_Tag < 1 | Fwd_Strand$End_Up_Tag < 20)
Fwd_Strand <- Fwd_Strand[-End_of_Reads,]

#Load package Biostrings to extract the sequences between the start and end of the tag
library(Biostrings)
#args[2]="SRR12105034.fasta"
Genome_Seq <- readDNAStringSet(args[2], format="fasta")

#Loop to extract the upstream tags 
for (line in c(1:nrow(Fwd_Strand)))
{
  Read_No <- Fwd_Strand$Strain_ID[line]
  Fwd_Strand$Upstream_Tag[line] <- as.character(Genome_Seq[[Read_No]][Fwd_Strand$Start_Up_Tag[line]:Fwd_Strand$End_Up_Tag[line]])
}

#Do the same for the reverse strand 
#Assign the ends 
Full_1 <- which(Rev_Strand$Start_Align_IS481 == 1)
Odd_2 <- which(Rev_Strand$Start_Align_IS481 == 2)
Odd_3 <- which(Rev_Strand$Start_Align_IS481 == 3)
Odd_4 <- which(Rev_Strand$Start_Align_IS481 == 4)
Odd_5 <- which(Rev_Strand$Start_Align_IS481 == 5)
Odd_6 <- which(Rev_Strand$Start_Align_IS481 == 6)
Odd_7 <- which(Rev_Strand$Start_Align_IS481 == 7)

#Add some columns to work out the upstream tag 
Rev_Strand$Start_Up_Tag <- NA
Rev_Strand$End_Up_Tag <- NA 
Rev_Strand$Upstream_Tag <- NA 

#Calculate the end of the upstream tag 
Rev_Strand$End_Up_Tag[Full_1] <- Rev_Strand$Start_Align_Strain[Full_1] + 2 
Rev_Strand$End_Up_Tag[Odd_2] <- Rev_Strand$Start_Align_Strain[Odd_2] + 3
Rev_Strand$End_Up_Tag[Odd_3] <- Rev_Strand$Start_Align_Strain[Odd_3] + 4 
Rev_Strand$End_Up_Tag[Odd_4] <- Rev_Strand$Start_Align_Strain[Odd_4] + 5 
Rev_Strand$End_Up_Tag[Odd_5] <- Rev_Strand$Start_Align_Strain[Odd_5] + 6 
Rev_Strand$End_Up_Tag[Odd_6] <- Rev_Strand$Start_Align_Strain[Odd_6] + 7
Rev_Strand$End_Up_Tag[Odd_7] <- Rev_Strand$Start_Align_Strain[Odd_7] + 8

#Calculate the start of the upstream tag
Rev_Strand$Start_Up_Tag[Full_1] <- Rev_Strand$End_Up_Tag[Full_1] + 19 
Rev_Strand$Start_Up_Tag[Odd_2] <- Rev_Strand$End_Up_Tag[Odd_2] + 19
Rev_Strand$Start_Up_Tag[Odd_3] <- Rev_Strand$End_Up_Tag[Odd_3] + 19 
Rev_Strand$Start_Up_Tag[Odd_4] <- Rev_Strand$End_Up_Tag[Odd_4] + 19
Rev_Strand$Start_Up_Tag[Odd_5] <- Rev_Strand$End_Up_Tag[Odd_5] + 19
Rev_Strand$Start_Up_Tag[Odd_6] <- Rev_Strand$End_Up_Tag[Odd_6] + 19
Rev_Strand$Start_Up_Tag[Odd_7] <- Rev_Strand$End_Up_Tag[Odd_7] + 19

#Now we need to add a control step where the script will not try and work out any tags for alignments that have ended at the end of the read
for (line in c(1:nrow(Rev_Strand)))
{
  Read_No <- Rev_Strand$Strain_ID[line]
  Rev_Strand$Max_Length[line] <- length(Genome_Seq[[Read_No]])
}

End_of_Reads <- which(Rev_Strand$Start_Up_Tag > Rev_Strand$Max_Length | Rev_Strand$End_Up_Tag > Rev_Strand$Max_Length)
Rev_Strand <- Rev_Strand[-End_of_Reads,]

#Loop to extract the upstream tags 
for (line in c(1:nrow(Rev_Strand)))
{
  Read_No <- Rev_Strand$Strain_ID[line]
  Up_Tag <- DNAStringSet(Genome_Seq[[Read_No]][Rev_Strand$End_Up_Tag[line]:Rev_Strand$Start_Up_Tag[line]])
  Rev_Strand$Upstream_Tag[line] <- reverseComplement(Up_Tag)
}

Fwd_Strand <- Fwd_Strand[,c(1, 2, 16)]
Rev_Strand <- Rev_Strand[,c(1, 2, 16)]

#Merge both upstream tables together 
Upstream_Both <- merge(Fwd_Strand, Rev_Strand, all = TRUE)

#Write to csv 
write.table(Upstream_Both, paste(args[1],"_Upstream_Tags.csv", sep=""),row.names = F,quote = F,sep = ",")

