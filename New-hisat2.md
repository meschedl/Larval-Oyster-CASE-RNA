cd Working-CASE-RNA
mkdir 2020-New-Analysis
cd 2020-New-Analysis

the trimmed files from April 2019 are in /RAID_STORAGE2/mschedl/RNA-CASE/

`ln -s /RAID_STORAGE2/mschedl/RNA-CASE/*.trim.fq.gz .`

also want the reference genome

`ln -s /RAID_STORAGE2/mschedl/RNA-CASE/cvir_edited.fa .`

This should be what I need for running hisat2

building index file shouldn't be any different
-x is the index file just made
fast q is the default so I don't really need to specify it
-- rna-strandedness "Specify strand-specific information: the default is unstranded.
For single-end reads, use F or R. 'F' means a read corresponds to a transcript. 'R' means a read corresponds to the reverse complemented counterpart of a transcript. For paired-end reads, use either FR or RF.
With this option being used, every read alignment will have an XS attribute tag: '+' means a read belongs to a transcript on '+' strand of genome. '-' means a read belongs to a transcript on '-' strand of genome."

I think I want to use FR. I am assuming that the way the reads are labeled means the right stranded ness

want to print out a summary file

after it will sort the SAM file into a bam file and remove the SAM file

`nano 2020hisat.sh`

```
##!/bin/bash

#Specify working directory
F=/home/mschedl/Working-CASE-RNA/2020-New-Analysis/

#Indexing a reference genome and no annotation file (allowing for novel transcript discovery)
#Build HISAT index with Cvirginica genome file, this file is from Erin Roberts, it has the extra space in the headers of each chromosome removed

hisat2-build -f $F/cvir_edited.fa $F/cvirginica_hisat_edited
#-f indicates that the reference input files are FASTA files

#Aligning paired end reads
#Has the F in here because the sed in the for loop changes it to a R. SAM files are of both forward and reverse reads
array1=($(ls $F/*F.trim.fq.gz))

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

for i in ${array1[@]}; do
        hisat2 --dta -x $F/cvirginica_hisat_edited -1 ${i} -2 $(echo ${i}|sed s/F.trim/R.trim/) --rna-strandness FR --new-summary 2020hisat2summary -S ${i}.sam
        samtools sort ${i}.sam > ${i}.s.bam
    		echo "${i}_bam"
        rm ${i}.sam
        echo "HISAT2 PE ${i}" $(date)
done
```

nohup bash 2020hisat.sh

ok no immediate errors

ok so potentially no errors with building the index but errors when running Hisat2

can't figure out what though

Extra parameter(s) specified: "FR", "2020hisat2summary"
Note that if <mates> files are specified using -1/-2, a <singles> file cannot
also be specified.  Please run bowtie separately for mates and singles.
Error: Encountered internal HISAT2 exception (#1)
Command: /home/mschedl/miniconda3/bin/hisat2-align-s --wrapper basic-0 --new-summary -S /home/mschedl/Working-CASE-RNA/2020-New-Analysis//CA_J06
.F.trim.fq.gz.sam -1 /tmp/75852.inpipe1 -2 /tmp/75852.inpipe2 rna-strandedness FR 2020hisat2summary
(ERR): hisat2-align exited with value 1
[E::hts_open_format] Failed to open file "/home/mschedl/Working-CASE-RNA/2020-New-Analysis//CA_J06.F.trim.fq.gz.sam" : No such file or directory
samtools sort: can't open "/home/mschedl/Working-CASE-RNA/2020-New-Analysis//CA_J06.F.trim.fq.gz.sam": No such file or directory
/home/mschedl/Working-CASE-RNA/2020-New-Analysis//CA_J06.F.trim.fq.gz_bam
rm: cannot remove ‘/home/mschedl/Working-CASE-RNA/2020-New-Analysis//CA_J06.F.trim.fq.gz.sam’: No such file or directory

let me see if I can try with one sample

hisat2 --dta -x cvirginica_hisat_edited -1 CA_J06.F.trim.fq.gz -2 CA_J06.R.trim.fq.gz -- rna-strandedness FR --new-summary 2020hisat2summary -S CA_J06.sam

this has the same problem but I don't know what it is

what if I remove the summary

hisat2 --dta -x cvirginica_hisat_edited -1 CA_J06.F.trim.fq.gz -2 CA_J06.R.trim.fq.gz --rna-strandness FR -S CA_J06.sam






```
##!/bin/bash

#Specify working directory
F=/home/mschedl/Working-CASE-RNA/2020-New-Analysis/

array1=($(ls $F/*F.trim.fq.gz))

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore

for i in ${array1[@]}; do
        hisat2 --dta -x $F/cvirginica_hisat_edited -1 ${i} -2 $(echo ${i}|sed s/F.trim/R.trim/) --rna-strandness FR -S ${i}.sam
        samtools sort ${i}.sam > ${i}.s.bam
    		echo "${i}_bam"
        rm ${i}.sam
        echo "HISAT2 PE ${i}" $(date)
done
```

problem was trying to name the summary file? might not be able to name it
