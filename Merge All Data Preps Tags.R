#Merge all data from all the strains into one dataframe 
#FOR R Studio: setwd("~/MSc Molecular Biosciences (Microbiology)/Research Project 2/Finding Insertions Project 1/Merged_csvs")
setwd("~/Sarah_project/Project/URLs/Merge")
temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.csv)

install.packages("reshape")
library(reshape)
All_Tags <- merge_recurse(myfiles)

#Make NAs = 0 
All_Tags[is.na(All_Tags)] <- 0

#Write as csv 
write.table(All_Tags, "All_Tags.csv",row.names = F,quote = F,sep = ",")
