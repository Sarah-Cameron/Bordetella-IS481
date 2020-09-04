#First import the All_Targets.csv which should have all the target sequences and their frequency in each strain  
All_Targets <- read.csv("All_Targets.csv", stringsAsFactors = FALSE)

#Create a column to collect the averages in 
All_Targets$Average <- rowMeans(All_Targets[, -1])

#Now we just need the target sequence and the average - adjust 
All_Targets <- All_Targets[c(1, 8)]

#Merge 2 halves together if necessary, make NAs = 0 and find rowSums again to leave 2 columns 
#Total_Targets <- merge(All_Targets, All_Targets2, by.x = "Target_Sequences", by.y = "Target_Sequences", all = TRUE)
#Total_Targets[is.na(Total_Targets)] <- 0
#Total_Targets$Total_Frequency <- rowSums(Total_Targets[, -1])
#Total_Targets <- Total_Targets[c(1, 4)]
#Total_Targets$Average <- Total_Targets$Total_Frequency/552
#Total_Targets <- Total_Targets[c(1, 3)]

#Put the average column in descending order 
#Total_Targets <- Total_Targets[order(Total_Targets$Average, decreasing = TRUE),]
#Exclude <- which(Total_Targets$Average < 1)
#Total_Targets <- Total_Targets[-Exclude,]

#Name the columns
#colnames(All_Targets) <- c("Target_Sequences", "B203", "E150", "E476", "E976", "I344", "J090", "Average")

#Create the graph using the average 
library(ggplot2)
ggplot(All_Targets, aes(x=reorder(Target_Sequence,-Average), Average)) +
  geom_bar(stat='identity') +
  scale_x_discrete(name = "Target Sequences") +
  scale_y_continuous(name = "Number per Genome", limits = c(0, 140)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 