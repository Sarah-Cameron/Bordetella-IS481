#First we need to load the csv files and use args 
args = commandArgs(trailingOnly = TRUE)
#args[1]='GCF_000306945.1_ASM30694v1_genomic.fna.BLAST.Result.txt_Full_Insertions.csv'
Target_Sequences <- read.csv(args[1], stringsAsFactors = FALSE)

#Make a column and find out whether the Up and Downstream targets match 
Target_Sequences$Matches <- Target_Sequences$Upstream_Target==Target_Sequences$Downstream_Target

#Save these column numbers as a variable
Matching_Target_Seqs <- which(Target_Sequences$Matches == TRUE)

#Then make a table of the targets and frequency of use from Upstream and Downstream Target Seqs
All_Up_Target_Seqs <- table(Target_Sequences$Upstream_Target)
All_Dwn_Target_Seqs <- table(Target_Sequences$Downstream_Target)

#Then save as data frames 
Up_Targets <- data.frame(All_Up_Target_Seqs)
Dwn_Targets <- data.frame(All_Dwn_Target_Seqs)

#To make a data frame of just the target sequences that match 
Matching_Targets <- table(Target_Sequences$Upstream_Target[Matching_Target_Seqs])
Matching_Targets <- data.frame(Matching_Targets)
colnames(Matching_Targets) <- c("Target_Sequences", args[1])
write.table(Matching_Targets, paste(args[1],"_Matching_Targets.csv", sep=""),row.names = F,quote = F,sep = ",")

#Merge Up and Down Targets together
Merge <- merge(Up_Targets, Dwn_Targets, by.x = "Var1", by.y = "Var1", all = TRUE)
colnames(Merge) <- c("Sequence", "Upstream_Freq", "Downstream_Freq")
Merge[is.na(Merge)] <- 0
Merge$Total_Freq <- NA

for(line in c(1:nrow(Merge)))
{
  Merge$Total_Freq[line] <- sum(Merge$Upstream_Freq[line]+Merge$Downstream_Freq[line])
  #print(Merge[line,])
}
Merge <- Merge[order(Merge$Total_Freq, decreasing = TRUE),]

#Give column names to new data frame 
Freq <- Merge[, c("Sequence", "Total_Freq")]
colnames(Freq) <- c("Target_Sequence", args[1])

#Assign strain name 
write.table(Freq, paste(args[1],"_Merged_Targets.csv", sep=""),row.names = F,quote = F,sep = ",")

