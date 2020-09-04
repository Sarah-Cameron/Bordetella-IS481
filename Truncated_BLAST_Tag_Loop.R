#First import the BLAST data file to the environment 
args = commandArgs(trailingOnly = TRUE)
#args[1]="GCF_004636965.1_ASM463696v1_genomic.fna.BLAST.Result.txt_Truncated_Insertions.csv"
IS481_BLAST <- read.csv(args[1],stringsAsFactors=F)
IS481_BLAST$File <- "Truncated_Insertion"

#Now we want to separate the IS481_BLAST results into Fwd and Rev strands 
Rev <- which(IS481_BLAST$Orientation == TRUE)
Rev_Strand <- IS481_BLAST[Rev,]
Fwd_Strand <- IS481_BLAST[-Rev,]

#Create two tables - one where there is the upstream end intact and one where there is the downstream end intact 
Up_Tags <- which(Fwd_Strand$Start_Align_IS481 <= 7)
Fwd_Strand_Upstream <- Fwd_Strand[Up_Tags,]
Down_Tags <- which(Fwd_Strand$End_Align_IS481 >= 1049)
Fwd_Strand_Downstream <- Fwd_Strand[Down_Tags,]

#Assign the different starting ends 
Full_1 <- which(Fwd_Strand_Upstream$Start_Align_IS481 == 1)
Odd_2 <- which(Fwd_Strand_Upstream$Start_Align_IS481 == 2) 
Odd_6 <- which(Fwd_Strand_Upstream$Start_Align_IS481 == 6)

#Make the columns for filling the positions and the tags  
Fwd_Strand_Upstream$Start_Up_Tag <- NA 
Fwd_Strand_Upstream$End_Up_Tag <- NA
Fwd_Strand_Upstream$Upstream_Tag <- NA 

#Find the positions of the end and start of the upstream tag 
Fwd_Strand_Upstream$End_Up_Tag[Full_1] <- Fwd_Strand_Upstream$Start_Align_Strain[Full_1] - 2 
Fwd_Strand_Upstream$End_Up_Tag[Odd_2] <- Fwd_Strand_Upstream$Start_Align_Strain[Odd_2] - 3 
Fwd_Strand_Upstream$End_Up_Tag[Odd_6] <- Fwd_Strand_Upstream$Start_Align_Strain[Odd_6] - 7 

Fwd_Strand_Upstream$Start_Up_Tag[Full_1] <- Fwd_Strand_Upstream$End_Up_Tag[Full_1] - 19 
Fwd_Strand_Upstream$Start_Up_Tag[Odd_2] <- Fwd_Strand_Upstream$End_Up_Tag[Odd_2] - 19
Fwd_Strand_Upstream$Start_Up_Tag[Odd_6] <- Fwd_Strand_Upstream$End_Up_Tag[Odd_6] - 19

#Load package Biostrings to extract the sequences between these points 
library(Biostrings)
#args[2]="GCF_004636965.1_ASM463696v1_genomic.fna"
Genome_Seq <- readDNAStringSet(args[2], format="fasta")

#Loop to extract the upstream tags 
for (line in c(1:nrow(Fwd_Strand_Upstream)))
{
  Fwd_Strand_Upstream$Upstream_Tag[line] <- as.character(Genome_Seq[[1]][Fwd_Strand_Upstream$Start_Up_Tag[line]:Fwd_Strand_Upstream$End_Up_Tag[line]])
}

#Now to do the same for the downstream tags 
#Subset rows based on the end of the alignment to the IS481
Odd_1049 <- which(Fwd_Strand_Downstream$End_Align_IS481 == 1049)
Odd_1050 <- which(Fwd_Strand_Downstream$End_Align_IS481 == 1050)
Full_1053 <- which(Fwd_Strand_Downstream$End_Align_IS481 == 1053)

#Add the columns 
Fwd_Strand_Downstream$Start_Dwn_Tag <- NA 
Fwd_Strand_Downstream$End_Dwn_Tag <- NA 
Fwd_Strand_Downstream$Downstream_Tag <- NA

#Find the positions of the start and end of the downstream tag 
Fwd_Strand_Downstream$Start_Dwn_Tag[Odd_1049] <- Fwd_Strand_Downstream$End_Align_Strain[Odd_1049] + 6 
Fwd_Strand_Downstream$Start_Dwn_Tag[Odd_1050] <- Fwd_Strand_Downstream$End_Align_Strain[Odd_1050] + 5 
Fwd_Strand_Downstream$Start_Dwn_Tag[Full_1053] <- Fwd_Strand_Downstream$End_Align_Strain[Full_1053] + 2 

Fwd_Strand_Downstream$End_Dwn_Tag[Odd_1049] <- Fwd_Strand_Downstream$Start_Dwn_Tag[Odd_1049] + 19
Fwd_Strand_Downstream$End_Dwn_Tag[Odd_1050] <- Fwd_Strand_Downstream$Start_Dwn_Tag[Odd_1050] + 19
Fwd_Strand_Downstream$End_Dwn_Tag[Full_1053] <- Fwd_Strand_Downstream$Start_Dwn_Tag[Full_1053] + 19

#Loop to extract the downstream tags 
for (line in c(1:nrow(Fwd_Strand_Downstream)))
{
  Fwd_Strand_Downstream$Downstream_Tag[line] <- as.character(Genome_Seq[[1]][Fwd_Strand_Downstream$Start_Dwn_Tag[line]:Fwd_Strand_Downstream$End_Dwn_Tag[line]])
}

#Do the same as above for the reverse strand 
#Create two tables - one where there is the upstream end intact and one where there is the downstream end intact 
Up_Tags <- which(Rev_Strand$Start_Align_IS481 <= 7)
Rev_Strand_Upstream <- Rev_Strand[Up_Tags,]
Down_Tags <- which(Rev_Strand$End_Align_IS481 >= 1049)
Rev_Strand_Downstream <- Rev_Strand[Down_Tags,]

#Assign the different starting ends 
Full_1 <- which(Rev_Strand_Upstream$Start_Align_IS481 == 1)
Odd_2 <- which(Rev_Strand_Upstream$Start_Align_IS481 == 2)
Odd_3 <- which(Rev_Strand_Upstream$Start_Align_IS481 == 3)

#Make the columns for filling the positions and the tags  
Rev_Strand_Upstream$Start_Up_Tag <- NA 
Rev_Strand_Upstream$End_Up_Tag <- NA
Rev_Strand_Upstream$Upstream_Tag <- NA 

#Find the positions of the end and start of the upstream tag 
Rev_Strand_Upstream$End_Up_Tag[Full_1] <- Rev_Strand_Upstream$Start_Align_Strain[Full_1] + 2 
Rev_Strand_Upstream$End_Up_Tag[Odd_2] <- Rev_Strand_Upstream$Start_Align_Strain[Odd_2] + 3 
Rev_Strand_Upstream$End_Up_Tag[Odd_3] <- Rev_Strand_Upstream$Start_Align_Strain[Odd_3] + 4 

Rev_Strand_Upstream$Start_Up_Tag[Full_1] <- Rev_Strand_Upstream$End_Up_Tag[Full_1] + 19 
Rev_Strand_Upstream$Start_Up_Tag[Odd_2] <- Rev_Strand_Upstream$End_Up_Tag[Odd_2] + 19
Rev_Strand_Upstream$Start_Up_Tag[Odd_3] <- Rev_Strand_Upstream$End_Up_Tag[Odd_3] + 19

#Loop to extract the upstream tags 
for (line in c(1:nrow(Rev_Strand_Upstream)))
{
  Upstream_Tag <- DNAStringSet(Genome_Seq[[1]][Rev_Strand_Upstream$End_Up_Tag[line]:Rev_Strand_Upstream$Start_Up_Tag[line]])
  Rev_Strand_Upstream$Upstream_Tag[line] <- reverseComplement(Upstream_Tag)
}

#Now to do the same for the downstream tags 
#Subset rows based on the end of the alignment to the IS481
Full_1053 <- which(Rev_Strand_Downstream$End_Align_IS481 == 1053)

#Add the columns 
Rev_Strand_Downstream$Start_Dwn_Tag <- NA 
Rev_Strand_Downstream$End_Dwn_Tag <- NA 
Rev_Strand_Downstream$Downstream_Tag <- NA

#Find the positions of the start and end of the downstream tag 
Rev_Strand_Downstream$Start_Dwn_Tag[Full_1053] <- Rev_Strand_Downstream$End_Align_Strain[Full_1053] - 2 

Rev_Strand_Downstream$End_Dwn_Tag[Full_1053] <- Rev_Strand_Downstream$Start_Dwn_Tag[Full_1053] - 19

#Loop to extract the downstream tags 
for (line in c(1:nrow(Rev_Strand_Downstream)))
{
  Downstream_Tag <- DNAStringSet(Genome_Seq[[1]][Rev_Strand_Downstream$End_Dwn_Tag[line]:Rev_Strand_Downstream$Start_Dwn_Tag[line]])
  Rev_Strand_Downstream$Downstream_Tag[line] <- reverseComplement(Downstream_Tag)
}


#Merge tables back together to get upstream and downstream tags together
Upstream_Both <- merge(Fwd_Strand_Upstream, Rev_Strand_Upstream, all = TRUE)
Downstream_Both <- merge(Fwd_Strand_Downstream, Rev_Strand_Downstream, all = TRUE)

#To get the columns needed
Upstream_Both <- Upstream_Both[,c(3, 15, 1, 18)]
Downstream_Both <- Downstream_Both[,c(3, 15, 1, 18)]

#Write to csvs 
write.table(Upstream_Both, paste(args[1],"_Upstream_Tags.csv", sep=""),row.names = F,quote = F,sep = ",")
write.table(Downstream_Both, paste(args[1],"_Downstream_Tags.csv", sep=""),row.names = F,quote = F,sep = ",")
