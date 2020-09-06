#Generate 10000 fake 300bp reads.
#half normal reads


library(Biostrings)
#328684 reads generated
UK38_fasta=readDNAStringSet("UK38_assem.fasta")
IS481=readDNAStringSet("M22031.1_IS481_Sequence.txt")
#Take 1/4
Ten_k_reads=list()
for(q in c(1:10000))
{
randy_start=sample(1:length(UK38_fasta[[1]])-300,1)
Ten_k_reads[q]=UK38_fasta[[1]][c(randy_start:(randy_start+149))]

}
UK38_fasta[[1]][1:100]
#First quarter will be reads which are made up by the first half of the read being normal DNA and the second half being the beginning of the IS.
First_quarter=Ten_k_reads[1:(length(Ten_k_reads)/4)]
#Second quarter will be reads which are made up by the first half of the read being the end of the IS
Second_quarter=Ten_k_reads[(1+(length(Ten_k_reads)/4)):(length(Ten_k_reads)/2)]

Normal_reads=Ten_k_reads[((length(Ten_k_reads)/2)+1):length(Ten_k_reads)]

for(q in c(1:(length(Ten_k_reads)/4)))
{
First_quarter[[q]][76:150]=IS481[[1]][1:75]
Second_quarter[[q]][1:75]=IS481[[1]][979:1053]
#print(q)
Normal_reads[[q]]=reverseComplement(Normal_reads[[q]])
if(q>(length(Ten_k_reads)/8))
{
First_quarter[[q]]=reverseComplement(First_quarter[[q]])
Second_quarter[[q]]=reverseComplement(Second_quarter[[q]])

}
}

#Rename them
names(First_quarter)[c(1:(length(First_quarter)/2))]=paste("Read",c(1:(length(First_quarter)/2)),"Forward_orient","IS_begining",sep="_")
names(First_quarter)[c((length(First_quarter)/2)+1):length(First_quarter)]=paste("Read",c((length(First_quarter)/2)+1):length(First_quarter),"Reverse_orient","IS_begining",sep="_")
#
names(Second_quarter)[c(1:(length(First_quarter)/2))]=paste("Read",c(1:(length(First_quarter)/2)),"Forward_orient","IS_end",sep="_")
names(Second_quarter)[c((length(First_quarter)/2)+1):length(First_quarter)]=paste("Read",c((length(First_quarter)/2)+1):length(First_quarter),"Reverse_orient","IS_end",sep="_")
names(Normal_reads)[1:c((length(Normal_reads)/2))]=paste("Read",c(1:(length(Normal_reads)/2)),"Normal","Forward_orient",sep="_")

names(Normal_reads)[c((length(Normal_reads)/2)+1):length(Normal_reads)]=paste("Read",c(c((length(Normal_reads)/2)+1):length(Normal_reads)),"Normal","Reverse_orient",sep="_")
#reverseComplement(First_quarter[[1:5]])

all_final_reads=c(First_quarter,Second_quarter,Normal_reads)
names(all_final_reads)
all_final_reads=DNAStringSet(all_final_reads)
writeXStringSet(all_final_reads, "Fake_reads.fa",format="fasta")

Key_dataframe=data.frame(Read_names=names(all_final_reads),Orientation=names(all_final_reads),IS=names(all_final_reads),stringsAsFactors=F)

Key_dataframe$IS[grepl("Normal",names(all_final_reads))]="None"
Key_dataframe$IS[grepl("IS_begining",names(all_final_reads))]="Begining"
Key_dataframe$IS[grepl("IS_end",names(all_final_reads))]="End"
Key_dataframe$Orientation[grepl("Forward",names(all_final_reads))]="Forward"
Key_dataframe$Orientation[grepl("Reverse",names(all_final_reads))]="Reverse"
write.csv(Key_dataframe,"Fake_data_key_data_frame.csv")
