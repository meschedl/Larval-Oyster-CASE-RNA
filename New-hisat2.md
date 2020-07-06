### Re-Running Histat2 with strandedness parameters

Create Directory
```
cd Working-CASE-RNA
mkdir 2020-New-Analysis
cd 2020-New-Analysis
```

The trimmed files from April 2019 are in /RAID_STORAGE2/mschedl/RNA-CASE/

`ln -s /RAID_STORAGE2/mschedl/RNA-CASE/*.trim.fq.gz .`

I also want the reference genome

`ln -s /RAID_STORAGE2/mschedl/RNA-CASE/cvir_edited.fa .`

This should be what I need for running hisat2


Parameters from Manual:
building index file shouldn't be any different
-x is the index file just made
fast q is the default so I don't really need to specify it
-- rna-strandedness "Specify strand-specific information: the default is unstranded.
For single-end reads, use F or R. 'F' means a read corresponds to a transcript. 'R' means a read corresponds to the reverse complemented counterpart of a transcript. For paired-end reads, use either FR or RF.
With this option being used, every read alignment will have an XS attribute tag: '+' means a read belongs to a transcript on '+' strand of genome. '-' means a read belongs to a transcript on '-' strand of genome."

I think I want to use FR. I am assuming that the way the reads are labeled means the right stranded ness

After I want to sort the SAM file into a bam file and remove the SAM file

`nano 2020hisat.sh`

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
