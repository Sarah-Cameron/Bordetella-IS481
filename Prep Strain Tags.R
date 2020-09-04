#First set to use args throughout the script 
args = commandArgs(trailingOnly = TRUE)
#args[1]='GCF_000193595.2_ASM19359v3_genomic.fna.BLAST.Result.txt_Full_Insertions.csv'
Insertions <- read.csv(args[1], stringsAsFactors = FALSE)

#Create a new table of just the strain id and the upstream and downstream location tags 
Tags_Only <- Insertions[, c(3, 25, 26)]

#Save as a new file 
write.table(Tags_Only, paste(args[1],"_All_Tags.csv", sep=""),row.names = F,quote = F,sep = ",")
