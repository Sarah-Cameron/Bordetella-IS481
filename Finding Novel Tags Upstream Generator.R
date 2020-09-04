#Argument 1 is a database of all the known tags in the closed genomes of 552 strains
args = commandArgs(trailingOnly = TRUE)
args[1]='ALL_TAGS_UPSTREAM.csv'
Known_Tags <- read.csv(args[1], stringsAsFactors = FALSE)

#Argument 2 is a database of all the tags calculated from the raw reads 
args[2]='SRR12168675.fasta.BLAST.Result.First.txt_Upstream_Tags.csv'
Raw_Tags <- read.csv(args[2], stringsAsFactors = FALSE)

#Load Biostrings package
library(Biostrings)

#As there are multiple repeats of the same tag in each list we need to only have a dataframe of the unique ones
Library <- table(Known_Tags$Upstream_Tag)
Library <- data.frame(Library)

#This to get rid of the title in one of the rows
Library <- Library[c(-496),]

#Now need to make pattern dictionary which we can use to search the raw reads for matches, setting the max.mismatch to 3 allows a variation of up to 3 nucleotides in a sequence to be identified as the same pattern
#Pre-processing it helps with the search 
DNA_Library <- DNAStringSet(Library$Var1)
DNA_Library <- PDict(DNA_Library, max.mismatch = 5)

#We then need to get all the unique tags from the raw reads and make it into a DNAStringSet 
Uniq_Raw_Tags <- table(Raw_Tags$Upstream_Tag)
Uniq_Raw_Tags <- data.frame(Uniq_Raw_Tags)
DNA_URAW <- DNAStringSet(Uniq_Raw_Tags$Var1)

#We can use vwhichPDict to use the dictionary of patterns to search the DNAStringSet for which rows have that pattern or a inexact match of it 
what <- vwhichPDict(DNA_Library, DNA_URAW, max.mismatch = 5, min.mismatch = 0, with.indels = FALSE, fixed = TRUE, algorithm = "auto")

#Using the tibble package and function enframe we can create a dataframe of the list generated 
#install.packages("tibble")
library(tibble)
WHAT <- enframe(what)

#Then using which we can find out the rows that got no matches i.e. are seemingly unique to the raw reads
rows_unique <- which(WHAT$value == "integer(0)")

#Then create a dataframe of these unique tags 
Unique_Tags <- as.character(DNA_URAW[rows_unique,])
Unique_Tags <- data.frame(Unique_Tags)

#Write the dataframe to csv 
write.table(Unique_Tags, paste(args[2], "_Unique_Raw_Upstream_Tags.csv", sep=""),row.names = F,quote = F,sep = ",")
write.csv(Unique_Tags, "UK36_Upstream_Unique_Tags.csv")
