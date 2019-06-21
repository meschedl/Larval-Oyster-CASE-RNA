##!/bin/bash

#Specify working directory
F=/home/mschedl/Working-CASE-RNA/histat/

#Indexing a reference genome and no annotation file (allowing for novel transcript discovery)
#Build HISAT index with Cvirginica genome file, this file is from Erin Roberts, it has the extra space in the headers of each chromosome removed

hisat2-build -f $F/cvir_edited.fa $F/cvirginica_hisat_edited
#-f indicates that the reference input files are FASTA files

#Aligning paired end reads
#Has the F in here because the sed in the for loop changes it to a R. SAM files are of both forward and reverse reads
array1=($(ls $F/*F.t.fq.gz))

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

for i in ${array1[@]}; do
        hisat2 --dta -x $F/cvirginica_hisat_edited -1 ${i} -2 $(echo ${i}|sed s/F.t/R.t/) -S ${i}.sam
        samtools sort ${i}.sam > ${i}.s.bam
    		echo "${i}_bam"
        rm ${i}.sam
        echo "HISAT2 PE ${i}" $(date)
done
