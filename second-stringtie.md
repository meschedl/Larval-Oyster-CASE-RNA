
I think using the -A flag gets rid of a lot of information in the merged gtf file for some reason...

Because the -A is supposed to make a gene abundance table and it never did.


Make new directory for trying it without the -A flag
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

Why is the second one so much smaller? I made mergelist.txt by copy and pasting the file names from the "try with the A" code, so it's not like I missed a file.

Try making a new final gtf for a sample and see if those differ

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


Ok those are the same numbers as above... is that right?
I would say yes because some of the coverages will be 0 but it makes a line for everything in the gtf file?

Well, the paper does not say to use the -A with the merge step, so I won't. I don't know what the difference is between them though... should I use gtf compare between both to the original annotation file and see what that says??

First make sure that everything in this directory is made with the No A file

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

/but what if I look at the files that the compare made? The stats file seems pertinent.

```
less c_vir_merged_compare_No_A.stats
```

# gffcompare v0.11.2 | Command line was:
#gffcompare -r ref_C_virginica-3.0_top_level.gff3 -G -o c_vir_merged_compare_No_A C_Vir_ST_merged_NO_A_test.gtf
#

#= Summary for dataset: C_Vir_ST_merged_NO_A_test.gtf
#     Query mRNAs :  108646 in   45439 loci  (100665 multi-exon transcripts)
#            (16934 multi-transcript loci, ~2.4 transcripts per locus)
# Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
# Super-loci w/ reference transcripts:    37492
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    88.9    |
        Exon level:   100.0     |    85.8    |
      Intron level:    99.9     |    90.5    |
Intron chain level:   100.0     |    65.0    |
  Transcript level:   100.0     |    62.5    |
       Locus level:   100.0     |    82.5    |

     Matching intron chains:   65384
       Matching transcripts:   67853
              Matching loci:   39152

          Missed exons:       0/352731  (  0.0%)
           Novel exons:   27597/425154  (  6.5%)
        Missed introns:     258/310704  (  0.1%)
         Novel introns:   10444/342908  (  3.0%)
           Missed loci:       0/39152   (  0.0%)
            Novel loci:    7947/45439   ( 17.5%)

 Total union super-loci across all input datasets: 45439
108646 out of 108646 consensus transcripts written in c_vir_merged_compare_No_A.annotated.gtf (0 discarded as redundant)


Ok what about the one with the A though, there has to be a reason it's smaller?

less c_vir_merged_compare_With_A.stats

# gffcompare v0.11.2 | Command line was:
#gffcompare -r ref_C_virginica-3.0_top_level.gff3 -G -o c_vir_merged_compare_With_A C_Vir_ST_merged_test.gtf
#

#= Summary for dataset: C_Vir_ST_merged_test.gtf
#     Query mRNAs :   79635 in   42550 loci  (72680 multi-exon transcripts)
#            (11864 multi-transcript loci, ~1.9 transcripts per locus)
# Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
# Super-loci w/ reference transcripts:    34880
#-----------------| Sensitivity | Precision  |
        Base level:    89.8     |    86.7    |
        Exon level:    86.4     |    83.6    |
      Intron level:    87.5     |    89.6    |
Intron chain level:    59.8     |    53.8    |
  Transcript level:    60.4     |    51.5    |
       Locus level:    82.5     |    75.2    |

     Matching intron chains:   39117
       Matching transcripts:   41012
              Matching loci:   32315

          Missed exons:   32408/352731  (  9.2%)
           Novel exons:   26576/372535  (  7.1%)
        Missed introns:   18845/310704  (  6.1%)
         Novel introns:   10302/303560  (  3.4%)
           Missed loci:    2576/39152   (  6.6%)
            Novel loci:    7417/42550   ( 17.4%)

 Total union super-loci across all input datasets: 42297
79635 out of 79635 consensus transcripts written in c_vir_merged_compare_With_A.annotated.gtf (0 discarded as redundant)


Ah-ha! It missed way more exons, introns, and loci for some reason! Well, I certainly want that. So the merge step without the -A is the way to go.

While it doesn't say to use -A in the merge step in the paper, it also doesn't say not to. I wonder why it does things differently, or why there isn't an error.


make the tables for DESeq2 with these Files

`nano list.sh`
```
#!/bin/bash

F=/home/mschedl/Working-CASE-RNA/histat/stringtie/restring


array2=($(ls *NO.A.merge.gtf))

for i in ${array2[@]}; do
    echo "$(echo ${i}|sed "s/\..*//") $F/${i}" >> sample_list.txt
done
```

looks good then run the prepDE.py script

```
conda activate python27

python prepDE.py -i sample_list.txt 

```
