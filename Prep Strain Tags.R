#First set to use args throughout the script 
args = commandArgs(trailingOnly = TRUE)
#args[1]='GCF_000193595.2_ASM19359v3_genomic.fna.BLAST.Result.txt_Full_Insertions.csv'
Insertions <- read.csv(args[1], stringsAsFactors = FALSE)
Insertions$File <- "Full_Insertion"

#Create a new table of just the strain id and the upstream and downstream location tags 
Upstream_Tags <- Insertions[, c(3, 27, 1, 25)]
Downstream_Tags <- Insertions[, c(3, 27, 1, 26)]

#Save as a new file 
write.table(Upstream_Tags, paste(args[1],"_Upstream_Tags.csv", sep=""),row.names = F,quote = F,sep = ",")
write.table(Downstream_Tags, paste(args[1],"_Downstream_Tags.csv", sep=""),row.names = F,quote = F,sep = ",")
