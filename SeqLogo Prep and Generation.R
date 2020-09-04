#Import the top 7 sequence targets on all the strains, asking it to not 'stringAsFactors' is important otherwise it creates problems later 
All_Targets <- read.csv("All_Targets_50_Strains.csv",stringsAsFactors = F)

#Now we need the total frequency of targets across the strains 
All_Targets$Total_Frequency <- rowSums(All_Targets[, -1])

#Now we just need the target sequence and the total frequency columns - adjust 
All_Targets <- All_Targets[c(1, 52)]

#Merge 2 halves together if necessary, make NAs = 0 and find rowSums again to leave 2 columns 
#Total_Targets <- merge(All_Targets, All_Targets2, by.x = "Target_Sequence", by.y = "Target_Sequence", all = TRUE)
#Total_Targets[is.na(Total_Targets)] <- 0
#Total_Targets$Total_Frequency <- rowSums(Total_Targets[, -1])
#Total_Targets <- Total_Targets[c(1, 4)]
#All_Targets <- All_Targets[order(All_Targets$Total_Frequency, decreasing = TRUE),]
#To get top 7 
#All_Targets <- All_Targets[c(1:7),]
#Remove all rows that have ambiguity letters M, K and S in 
#Total_Targets <- Total_Targets[-c(68, 69, 70, 71, 81, 82),]

#Now we are left with top 7 target sequences in each strain which we have combined and now we need to determine the frequency of each base at each of the 6 positions 
#We need to make a results table that we can fill in this information 
Results_Table <- data.frame(Pos1=rep(0,4),Pos2=rep(0,4),Pos3=rep(0,4),Pos4=rep(0,4),Pos5=rep(0,4),Pos6=rep(0,4))

#Label characters (needs to be in A, C, G, T format for the SeqLogo step)
rownames(Results_Table) <- c("A", "C", "G", "T")

#Now we want to create a loop that will go through each sequence and split the sequence into the individual letters
for(line in c(1:nrow(All_Targets)))
{
  #This will split it 
  Letters <- strsplit(All_Targets$Target_Sequence[line],"")
  #This will make them individual and not in a list 
  Letters <- unlist(Letters)
  #Now we need to create another loop (within this loop!) to loop through each of the letters we have just unlisted and add the totals to the results table
  for(L in c(1:length(Letters)))
  {
  #This line is asking (letter by letter along the length of the sequence) which are at each position and assigns it a variable. 
  Appropriate_row=which(Letters[L]==row.names(Results_Table))
  #This line of code then says take that position and add what is already there to what is in the current line of the data table (the frequency that sequence is used)
  Results_Table[Appropriate_row,L]=sum(Results_Table[Appropriate_row,L]+All_Targets$Total_Frequency[line])
  }
}

#Then to create a weighted table we then need each column to be a percentage - to see as a percentage how often each base was used at that position
Total_Freq <- sum(Results_Table$Pos1)
for (line in c(1:nrow(Results_Table)))
{ 
  Results_Table$Pos1[line] <- Results_Table$Pos1[line]/Total_Freq
  Results_Table$Pos2[line] <- Results_Table$Pos2[line]/Total_Freq
  Results_Table$Pos3[line] <- Results_Table$Pos3[line]/Total_Freq
  Results_Table$Pos4[line] <- Results_Table$Pos4[line]/Total_Freq
  Results_Table$Pos5[line] <- Results_Table$Pos5[line]/Total_Freq
  Results_Table$Pos6[line] <- Results_Table$Pos6[line]/Total_Freq
}

#It is worth checking that the columns add up to 1 before using SeqLogo otherwise will throw up errors 
sum(Results_Table$Pos1)

#Install SeqLogo Package if needed  
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("seqLogo")

#Making the consensus motif plot using SeqLogo
library(seqLogo)
PWM <- Results_Table
Weight_Table <- makePWM(PWM)
seqLogo(Weight_Table)