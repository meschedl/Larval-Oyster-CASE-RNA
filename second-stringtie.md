
I think using the -A flag gets rid of a lot of information in the merged gtf file for some reason...

stringtie -e -G C_Vir_ST_merged.gtf -A gene_abund.tab -o CASE_J03_2.merge.gtf CASE_J03.F.trim.fq.gz.s.bam


make new directory for trying it without the -A flag
```
mkdir restring
cd restring

ln -s /RAID_STORAGE2/mschedl/RNA-CASE/ref_C_virginica-3.0_top_level.gff3 .
ln -s /home/mschedl/Working-CASE-RNA/histat/*.s.bam .
ln -s /home/mschedl/Working-CASE-RNA/histat/stringtie/*.bam.gtf .
```


Test the re-estimated after merge gtf files are/are not different using the -A?

try with the A

```
stringtie --merge -A -G ref_C_virginica-3.0_top_level.gff3 -o C_Vir_ST_merged_test.gtf CA_J06.F.trim.fq.gz.s.bam.gtf CA_J08.F.trim.fq.gz.s.bam.gtf  CA_J11.F.trim.fq.gz.s.bam.gtf CA_J18.F.trim.fq.gz.s.bam.gtf CASE_J03.F.trim.fq.gz.s.bam.gtf CASE_J09.F.trim.fq.gz.s.bam.gtf CASE_J12.F.trim.fq.gz.s.bam.gtf CASE_J13.F.trim.fq.gz.s.bam.gtf CON_J02.F.trim.fq.gz.s.bam.gtf CON_J05.F.trim.fq.gz.s.bam.gtf CON_J10.F.trim.fq.gz.s.bam.gtf SE_J01.F.trim.fq.gz.s.bam.gtf SE_J04.F.trim.fq.gz.s.bam.gtf SE_J07.F.trim.fq.gz.s.bam.gtf
```

Try without A
```
stringtie --merge -G ref_C_virginica-3.0_top_level.gff3 -o C_Vir_ST_merged_NO_A_test.gtf mergelist.txt
```
What is the difference between them/is there?

```
wc -l C_Vir_ST_merged_NO_A_test.gtf
```

1252378 C_Vir_ST_merged_NO_A_test.gtf
```
wc -l C_Vir_ST_merged_test.gtf
```
835503 C_Vir_ST_merged_test.gtf

why is the second one so much smaller? I made mergelist.txt by copy and pasting the file names from the "try with the A" code, so it's not like I missed a file....  

try making a new final gtf for a sample and see if those differ

```
stringtie -e -G C_Vir_ST_merged_NO_A_test.gtf -o SE_J07.test.merge.gtf SE_J07.F.trim.fq.gz.s.bam

stringtie -e -G C_Vir_ST_merged_test.gtf -o SE_J07.test.with.A.merge.gtf SE_J07.F.trim.fq.gz.s.bam
```
And are those different?
```
wc -l SE_J07.test.merge.gtf
```

1252378 SE_J07.test.merge.gtf
```
wc -l SE_J07.test.with.A.merge.gtf
```
835503 SE_J07.test.with.A.merge.gtf


ok those are the same numbers as above... is that right?
yes because some of the coverages will be 0 but it makes a line for everything in the gtf file?

well, the paper does not say to use the -A with the merge step, so I won't. I don't know what the difference is between them though... should I use gtf compare between both to the original annotation file and see what that says??

first make sure that everything in this directory is made with the No A file

```
nano NO_A_Stringtie.sh

#!/bin/bash
# In the same directory now with the BAM files and the annotation file link
F=/home/mschedl/Working-CASE-RNA/histat/stringtie/restring/

array1=($(ls $F/*.bam))

Re-estimate transcript abundance after merge step
	for i in ${array1[@]}; do
		stringtie -e -G $F/C_Vir_ST_merged_NO_A_test.gtf -o $(echo ${i}|sed "s/\..*//").NO.A.merge.gtf ${i}
		echo "${i}"
	done
	# input here is the original set of alignment bam files
	# here -G refers to the merged GTF files
	# -e creates more accurate abundance estimations with input transcripts, needed when converting to DESeq2 tables, says in the manual that this is recommended
echo "DONE" $(date)
```
And then compare?
```
gffcompare -r ref_C_virginica-3.0_top_level.gff3 -G -o c_vir_merged_compare_No_A C_Vir_ST_merged_NO_A_test.gtf
```

67891 reference transcripts loaded.
38 duplicate reference transcripts discarded.
108646 query transfrags loaded.

```
gffcompare -r ref_C_virginica-3.0_top_level.gff3 -G -o c_vir_merged_compare_With_A C_Vir_ST_merged_test.gtf
```

67891 reference transcripts loaded.
  38 duplicate reference transcripts discarded.
  79635 query transfrags loaded.
