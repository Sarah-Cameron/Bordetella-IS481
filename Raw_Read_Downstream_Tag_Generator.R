#First import the BLAST data file to the environment 
args = commandArgs(trailingOnly = TRUE)
#args[1]="SRR12105034.fasta.BLAST.Result.Last.txt"
IS481_BLAST <- read.delim(args[1],stringsAsFactors=F, head=F)

#Assign the column names 
colnames(IS481_BLAST) <- c("IS481_ID", "Strain_ID", "%_Identical", "Align_Length", "No_Mismatches", "No_Gaps", "Start_Align_IS481", "End_Align_IS481", "Start_Align_Strain", "End_Align_Strain", "EValue", "BitScore", "Strand")

#Remove the matches that haven't aligned at the start 
Truncated <- which(IS481_BLAST$End_Align_IS481 < 34)
Truncated_ISE <- IS481_BLAST[Truncated,]
IS481_BLAST <- IS481_BLAST[-Truncated,]

#Now we want to separate the IS481_BLAST results into Fwd and Rev strands 
Rev <- which(IS481_BLAST$Strand == "minus")
Rev_Strand <- IS481_BLAST[Rev,]
Fwd_Strand <- IS481_BLAST[-Rev,]

#Assign the ends 
Odd_1047 <- which(Fwd_Strand$End_Align_IS481 == 34)
Odd_1048 <- which(Fwd_Strand$End_Align_IS481 == 35)
Odd_1049 <- which(Fwd_Strand$End_Align_IS481 == 36)
Odd_1050 <- which(Fwd_Strand$End_Align_IS481 == 37)
Odd_1051 <- which(Fwd_Strand$End_Align_IS481 == 38)
Odd_1052 <- which(Fwd_Strand$End_Align_IS481 == 39)
Full_1053 <- which(Fwd_Strand$End_Align_IS481 == 40)

#Add some columns to work out the upstream tag 
Fwd_Strand$Start_Dwn_Tag <- NA
Fwd_Strand$End_Dwn_Tag <- NA 
Fwd_Strand$Downstream_Tag <- NA 

#Calculate the end of the downstream tag 
Fwd_Strand$Start_Dwn_Tag[Odd_1047] <- Fwd_Strand$End_Align_Strain[Odd_1047] + 8 
Fwd_Strand$Start_Dwn_Tag[Odd_1048] <- Fwd_Strand$End_Align_Strain[Odd_1048] + 7
Fwd_Strand$Start_Dwn_Tag[Odd_1049] <- Fwd_Strand$End_Align_Strain[Odd_1049] + 6
Fwd_Strand$Start_Dwn_Tag[Odd_1050] <- Fwd_Strand$End_Align_Strain[Odd_1050] + 5
Fwd_Strand$Start_Dwn_Tag[Odd_1051] <- Fwd_Strand$End_Align_Strain[Odd_1051] + 4 
Fwd_Strand$Start_Dwn_Tag[Odd_1052] <- Fwd_Strand$End_Align_Strain[Odd_1052] + 3
Fwd_Strand$Start_Dwn_Tag[Full_1053] <- Fwd_Strand$End_Align_Strain[Full_1053] + 2 

#Calculate the start of the upstream tag
Fwd_Strand$End_Dwn_Tag[Odd_1047] <- Fwd_Strand$Start_Dwn_Tag[Odd_1047] + 19 
Fwd_Strand$End_Dwn_Tag[Odd_1048] <- Fwd_Strand$Start_Dwn_Tag[Odd_1048] + 19
Fwd_Strand$End_Dwn_Tag[Odd_1049] <- Fwd_Strand$Start_Dwn_Tag[Odd_1049] + 19 
Fwd_Strand$End_Dwn_Tag[Odd_1050] <- Fwd_Strand$Start_Dwn_Tag[Odd_1050] + 19
Fwd_Strand$End_Dwn_Tag[Odd_1051] <- Fwd_Strand$Start_Dwn_Tag[Odd_1051] + 19
Fwd_Strand$End_Dwn_Tag[Odd_1052] <- Fwd_Strand$Start_Dwn_Tag[Odd_1052] + 19
Fwd_Strand$End_Dwn_Tag[Full_1053] <- Fwd_Strand$Start_Dwn_Tag[Full_1053] + 19

#Load package Biostrings to extract the sequences between the start and end of the tag
library(Biostrings)
#args[2]="SRR12105034.fasta"
Genome_Seq <- readDNAStringSet(args[2], format="fasta")

#Now we need to add a control step where the script will not try and work out any tags for alignments that have ended at the end of the read (aka have minus numbers in the positon columns)
for (line in c(1:nrow(Fwd_Strand)))
{
  Read_No <- Fwd_Strand$Strain_ID[line]
  Fwd_Strand$Max_Length[line] <- length(Genome_Seq[[Read_No]])
}

End_of_Reads <- which(Fwd_Strand$Start_Dwn_Tag > Fwd_Strand$Max_Length | Fwd_Strand$End_Dwn_Tag > Fwd_Strand$Max_Length)
Fwd_Strand <- Fwd_Strand[-End_of_Reads,]


#Loop to extract the upstream tags 
for (line in c(1:nrow(Fwd_Strand)))
{
  Read_No <- Fwd_Strand$Strain_ID[line]
  Fwd_Strand$Downstream_Tag[line] <- as.character(Genome_Seq[[Read_No]][Fwd_Strand$Start_Dwn_Tag[line]:Fwd_Strand$End_Dwn_Tag[line]])
}

#Do the same for the reverse strand 
#Assign the ends 
Odd_1047 <- which(Rev_Strand$End_Align_IS481 == 34)
Odd_1048 <- which(Rev_Strand$End_Align_IS481 == 35)
Odd_1049 <- which(Rev_Strand$End_Align_IS481 == 36)
Odd_1050 <- which(Rev_Strand$End_Align_IS481 == 37)
Odd_1051 <- which(Rev_Strand$End_Align_IS481 == 38)
Odd_1052 <- which(Rev_Strand$End_Align_IS481 == 39)
Full_1053 <- which(Rev_Strand$End_Align_IS481 == 40)

#Add some columns to work out the upstream tag 
Rev_Strand$Start_Dwn_Tag <- NA
Rev_Strand$End_Dwn_Tag <- NA 
Rev_Strand$Downstream_Tag <- NA 

#Calculate the end of the downstream tag 
Rev_Strand$Start_Dwn_Tag[Odd_1047] <- Rev_Strand$End_Align_Strain[Odd_1047] - 8 
Rev_Strand$Start_Dwn_Tag[Odd_1048] <- Rev_Strand$End_Align_Strain[Odd_1048] - 7
Rev_Strand$Start_Dwn_Tag[Odd_1049] <- Rev_Strand$End_Align_Strain[Odd_1049] - 6 
Rev_Strand$Start_Dwn_Tag[Odd_1050] <- Rev_Strand$End_Align_Strain[Odd_1050] - 5 
Rev_Strand$Start_Dwn_Tag[Odd_1051] <- Rev_Strand$End_Align_Strain[Odd_1051] - 4 
Rev_Strand$Start_Dwn_Tag[Odd_1052] <- Rev_Strand$End_Align_Strain[Odd_1052] - 3
Rev_Strand$Start_Dwn_Tag[Full_1053] <- Rev_Strand$End_Align_Strain[Full_1053] - 2

#Calculate the start of the upstream tag
Rev_Strand$End_Dwn_Tag[Odd_1047] <- Rev_Strand$Start_Dwn_Tag[Odd_1047] - 19 
Rev_Strand$End_Dwn_Tag[Odd_1048] <- Rev_Strand$Start_Dwn_Tag[Odd_1048] - 19
Rev_Strand$End_Dwn_Tag[Odd_1049] <- Rev_Strand$Start_Dwn_Tag[Odd_1049] - 19 
Rev_Strand$End_Dwn_Tag[Odd_1050] <- Rev_Strand$Start_Dwn_Tag[Odd_1050] - 19
Rev_Strand$End_Dwn_Tag[Odd_1051] <- Rev_Strand$Start_Dwn_Tag[Odd_1051] - 19
Rev_Strand$End_Dwn_Tag[Odd_1052] <- Rev_Strand$Start_Dwn_Tag[Odd_1052] - 19
Rev_Strand$End_Dwn_Tag[Full_1053] <- Rev_Strand$Start_Dwn_Tag[Full_1053] - 19

#Now we need to add a control step where the script will not try and work out any tags for alignments that have ended at the end of the read
End_of_Reads <- which(Rev_Strand$Start_Dwn_Tag < 20 | Rev_Strand$End_Dwn_Tag < 1)
Rev_Strand <- Rev_Strand[-End_of_Reads,]

#Loop to extract the upstream tags 
for (line in c(1:nrow(Rev_Strand)))
{
  Read_No <- Rev_Strand$Strain_ID[line]
  Dwn_Tag <- DNAStringSet(Genome_Seq[[Read_No]][Rev_Strand$End_Dwn_Tag[line]:Rev_Strand$Start_Dwn_Tag[line]])
  Rev_Strand$Downstream_Tag[line] <- reverseComplement(Dwn_Tag)
}

Fwd_Strand <- Fwd_Strand[,c(1, 2, 16)]
Rev_Strand <- Rev_Strand[,c(1, 2, 16)]

#Merge both upstream tables together 
Upstream_Both <- merge(Fwd_Strand, Rev_Strand, all = TRUE)

#Write to csv 
write.table(Upstream_Both, paste(args[1],"_Downstream_Tags.csv", sep=""),row.names = F,quote = F,sep = ",")

