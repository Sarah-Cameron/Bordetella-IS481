#!/bin/bash
#Loop for making each genome into a database
files=./Raw_Reads/SRR12168673.fasta
for file in $files
do
makeblastdb -in $file -parse_seqids -dbtype nucl;
done

#Loop for conducting BLAST search in each database
for database in ./Raw_Reads/SRR12168673.fasta; do
blastn \
-query First_40_IS481.txt \
-db $database \
-max_target_seqs 50000 \
-outfmt "6 std sstrand" \
-out $database.BLAST.Result.First.txt;
done
