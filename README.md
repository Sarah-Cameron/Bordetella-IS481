# Bordetella-IS481

Target Sequence Analysis of IS481 in Bordetella pertussis. 

1. BLAST_Commands.sh
These commands are for a nucleotide BLAST+ written in loop to allow it to be conducted on multiple genomes. Adapt this script to the directory of your genomes/short reads, the IS481 sequence you wish to use and the name of the output file. This will give the BLAST results in tabular format 6.

For Linux use this command to generate a list of all the genomes. This will help later to make sure each file in the next script is run against the correct corresponding file for that strain. 
$ find ./ -iname "*.fna"|sed 's/FIND/REPLACE/g'|sed 's/FIND/REPLACE/g' >All_Genomes.txt

2. Command_Line_BLAST_Loop.R 
This script will go through each of the BLAST results from each file, separate the full insertions (>1047) from those that are truncated, it will then work out the positions for the target sequence within the IS and extract its sequence (the same applies for flanking 20bp tags). If the IS has aligned on the reverse strand it will work out the reverse complement. 
$ cat All_Genomes.txt| xargs -d '\n' -P 8 -n 1 -I file Rscript Command_Line_BLAST_Loop.R file.BLAST.Result.txt file

3. Data Prep to All Targets + Matching Targets.R
This script will take each genome's individual full insertions.csv and process to an output file with just all the target sequences and their frequency within that genome.

4. Merge All Data Preps.R 
This script will merge all of the previously created files together. This will give one table with all target sequences in individual rows and column a different genome.

5. SeqLogo Prep and Generation.R 
To generate a consensus motif of all the target sequences this script will first create a position weight matrix of the 6 bp target sequence and then using seqLogo will generate the consensus motif.

6. Frequency_Graph.R
This will take the All_Targets.csv file and find the average use of each target sequence per genome and plot it into a histogram using ggplot2.


Creating Library of Known Locations for IS481 from Closed Genomes 

1. BLAST_Commands.sh
These commands are for a nucleotide BLAST+ written in loop to allow it to be conducted on multiple genomes. Adapt this script to the directory of your genomes, the IS481 sequence you wish to use and the name of the output file. This will give the BLAST results in tabular format 6.

2. Command_Line_BLAST_Loop.R 
This script will go through each of the BLAST results from each file, separate the full insertions (>1047) from those that are truncated, it will then work out the positions for the 20 bp flanking tags and extract its sequence. If the IS has aligned on the reverse strand it will work out the reverse complement. 

3. Truncated_BLAST_Tag_Loop.R
This script was used to go back and extract the tags from all the files that were deemed truncated in the previous script. Tags are only created from ends of the IS481 that aligned. 

4. Prep Strain Tags.R 
Used to create 2 files for each genome one with the upstream tags and one with the downstream tags. 
Separate all upstream tags into one directory and all downstream tags into another.

5. Merge All Data Preps Tags.R 
Adapt directories to files merging and output names etc. Merges all the individual files to create one large file once with the upstream tags and again with the downstream files. 



Creating Tags From Short Read Sequence Data

Download the paired-end short read sequence data to your server. 

Linux commands:
If the files download as .gz files, use ‘gunzip *.gz’ to unzip them.

Convert each of the fastq files into fasta files with:
sed -n '1~4s/^@/>/p;2~4p' SRR12168682_1.fastq > SRR12168682_1.fasta

Rename each of the reads within these files to allow for traceability back to the individual read.
awk '/^>/{print ">File_1." ++i; next}{print}' SRR12168682_1.fasta >SRR12168682_F1.fasta

Merge both of these files together using
cat SRR12168682_F1.fasta SRR12168682_F2.fasta >SRR12168682.fasta 

1. BLAST_Commands.sh 
Adapt to BLAST the first 40 bp of the IS481 sequence against the raw short read sequence data. Adapt again and run for the last 40 bp changing the output file name.
$ sh BLAST_Commands.sh 

2. Raw_Read_Upstream_Tag_Generator.R
Takes all the BLAST results from using the first part of the IS481 and generates the 20 bp flanking tag.
$ Rscript Raw_Read_Upstream_Tag_Generator.R SRR12168682.fasta.BLAST.Result.First.txt SRR12168682.fasta

3. Raw_Read_Downstream_Tag_Generator.R
Takes all the BLAST results from using the last part of the IS481 and generates the 20bp flanking tag.
$ Rscript Raw_Read_Downstream_Tag_Generator.R SRR12168682.fasta.BLAST.Result.Last.txt SRR12168682.fasta

4. Finding Novel Tags Upstream Generator.R 
Creates a pattern dictionary of all the known upstream tags and compares them to the tags generated from raw reads. Giving a list of the ones that are unique to the raw reads.

5. Finding Novel Tags Downstream Generator.R 
Creates a pattern dictionary of all the known downstream tags and compares them to the tags generated from raw reads. Giving a list of the ones that are unique to the raw reads.

6. Finding Novel Tags Upstream Same Strain.R / Findind Novel Tags Downstream Same Strain.R 
Same as before but only creates a pattern dictionary of the known tags of that strain. Need to enter the strain ID.


Method Test With Fake Reads 

1. Generate random fake reads from a closed genome using snippy. 

2. Generating Fake Reads.R 
To create positives for the data set, add random portions of the genome next to either the front / back of the IS481 sequence. 

3. Finding Fake Dataset Files.R 
Use this script to locate which reads have been found from the list of unique tags previously generated.  
