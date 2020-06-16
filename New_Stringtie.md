New stringtie

-f <0.0-1.0>	Sets the minimum isoform abundance of the predicted transcripts as a fraction of the most abundant transcript assembled at a given locus. Lower abundance transcripts are often artifacts of incompletely spliced precursors of processed transcripts. Default: 0.01

-j <float>	There should be at least this many spliced reads that align across a junction (i.e. junction coverage). This number can be fractional, since some reads align in more than one place. A read that aligns in n places will contribute 1/n to the junction coverage. Default: 1

-c <float>	Sets the minimum read coverage allowed for the predicted transcripts. A transcript with a lower coverage than this value is not shown in the output. Default: 1


basically I should increase all of these??

`mkdir 2020-stringtie`
`cd 2020-stringtie`
`ln -s /RAID_STORAGE2/mschedl/RNA-CASE/ref_C_virginica-3.0_top_level.gff3 .`
`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/*.s.bam .`

try to find literature on this??


from cufflinks, the pre-stringtie??http://cole-trapnell-lab.github.io/cufflinks/cufflinks/
"After calculating isoform abundance for a gene, Cufflinks filters out transcripts that it believes are very low abundance, because isoforms expressed at extremely low levels often cannot reliably be assembled, and may even be artifacts of incompletely spliced precursors of processed transcripts. This parameter is also used to filter out introns that have far fewer spliced alignments supporting them. The default is 0.1, or 10% of the most abundant isoform (the major isoform) of the gene"
this is a lot higher than the defualt of the stringtie one

maybe should set to 0.1

 -j minimum junction coverage (default: 1)

 maybe increase to 2

"--conservative : conservative transcriptome assembly, same as -t -c 1.5 -f 0.05"

-t disable trimming of predicted transcripts based on coverage
   (default: coverage trimming is enabled)

I don't understand how disabling trimming is conservative, I want to keep that I think

-s minimum reads per bp coverage to consider for single-exon transcript
    (default: 4.75)

how about increase this to 7

-c minimum reads per bp coverage to consider for multi-exon transcript
    (default: 1)

increase this to 3

should make -f .1 from the cufflinks example

and make sure to do strandedness

--rf "--rf	Assumes a stranded library fr-firststrand"



download new annotation file, might just be formatted differently https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/6565/100/GCF_002022765.2_C_virginica-3.0/

scp to KITT

scp -P 2292 /Users/maggieschedl/Desktop/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz mschedl@KITT.uri.edu:/home/mschedl/Working-CASE-RNA/2020-New-Analysis/


/home/mschedl/Working-CASE-RNA/2020-New-Analysis/2020-stringtie/

`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/GCF_002022765.2_C_virginica-3.0_genomic.gff .`

nano first-stringtie.sh
```
#!/bin/bash
# In the same directory now with the BAM files and the annotation file link
F=/home/mschedl/Working-CASE-RNA/2020-New-Analysis/2020-stringtie

# StringTie to assemble transcripts for each sample with the annotation file
array1=($(ls $F/*.bam))
for i in ${array1[@]}; do
	stringtie -G $F/GCF_002022765.2_C_virginica-3.0_genomic.gff --rf -f .1 -c 3 -s 7 -j 2 -o ${i}.gtf ${i}
	echo "${i}"
done
```


so far says:
Coverage saturation parameter is deprecated starting at version 1.0.5
there isn't a version 1.0.5????

still was able to run, so not sure how that error message messed things up

`mv nohup.out stringtie1.out`

Now to do the merge mode, this will make all the transcripts have the same name and be comparable to each other in each sample. need to make a list file to call the gtf files to use in it and need to use reference annotation file again

merge options to look at:

-i	keep merged transcripts with retained introns (default: these are not kept unless there is strong evidence for them) Don't think I want to do this

-f <min_iso>	minimum isoform fraction (default: 0.01) this should be 0.1 again

-c <min_cov>	minimum input transcript coverage to include in the merge (default: 0)
-F <min_fpkm>	minimum input transcript FPKM to include in the merge (default: 0)
-T <min_tpm>	minimum input transcript TPM to include in the merge (default: 0)

I actually don't think I want to change these, because if one sample has 0 I don't want it to be removed from the merged set? I do not know if it is 0 for 1 or 0 across all... across all doesn't seem like it would make sense so it must be 0 for one.

`nano merge-gtf-list.txt`
CA_J06.F.trim.fq.gz.s.bam.gtf
CA_J08.F.trim.fq.gz.s.bam.gtf
CA_J11.F.trim.fq.gz.s.bam.gtf
CA_J18.F.trim.fq.gz.s.bam.gtf
CASE_J03.F.trim.fq.gz.s.bam.gtf
CASE_J09.F.trim.fq.gz.s.bam.gtf
CASE_J12.F.trim.fq.gz.s.bam.gtf
CASE_J13.F.trim.fq.gz.s.bam.gtf
CON_J02.F.trim.fq.gz.s.bam.gtf
CON_J05.F.trim.fq.gz.s.bam.gtf
CON_J10.F.trim.fq.gz.s.bam.gtf
SE_J01.F.trim.fq.gz.s.bam.gtf
SE_J04.F.trim.fq.gz.s.bam.gtf
SE_J07.F.trim.fq.gz.s.bam.gtf


`stringtie --merge -G GCF_002022765.2_C_virginica-3.0_genomic.gff -o 2020-Cvir-merged.gtf merge-gtf-list.txt`

done very quickly, now want to use gff compare

https://ccb.jhu.edu/software/stringtie/gffcompare.shtml

-r	An optional “reference” annotation GFF file. Each sample is matched against this file, and sample isoforms are tagged as overlapping, matching, or novel where appropriate. See the .refmap and .tmap output file descriptions below.
unsure about any of the other options...


`gffcompare -r GCF_002022765.2_C_virginica-3.0_genomic.gff -o Cvir-compare 2020-Cvir-merged.gtf`

67891 reference transcripts loaded.
 38 duplicate reference transcripts discarded.
 132608 query transfrags loaded.

 .stats file has info

 interesting that I have HIGH sensitivity but low precision, this is "worse" than before

Sensitivity = TP / (TP+FN)

Precision = TP / (TP+FP)

Here TP means "true positives", which in this case are the query features (bases, exons, introns, transcripts, etc.) which agree with the corresponding reference annotation features
FN means "false negatives", i.e. features which are found in the reference annotation but missed (not present) in the "query" data. Accordingly, FP are features present in the input "query" data but not confirmed by any reference annnotation data. Notice that FP+TP amounts to the whole input set of "query" features (per sample/input file)

TP (true positives) number is the number of exon bases that are reported at the same coordinate on both the query (the assembled transfrag) and *any* reference transcript (that is, the overlap length of exons)
FN (false negatives) is the number of bases in reference data exons which are not covered at all by any of the predicted transcripts(transfrags) exons;
FP (false positives) then is the number of bases which are only covered by any predicted transcript's exons but not covered by any reference transcript exons


old stringtie (did not have any strand information and did not have any stringent parameters):
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

new stringtie:

#= Summary for dataset: 2020-Cvir-merged.gtf
#     Query mRNAs :  132608 in   67496 loci  (122412 multi-exon transcripts)
#            (21190 multi-transcript loci, ~2.0 transcripts per locus)
# Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
# Super-loci w/ reference transcripts:    38987
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    59.1    |
        Exon level:   100.0     |    56.5    |
      Intron level:   100.0     |    59.8    |
Intron chain level:   100.0     |    53.4    |
  Transcript level:   100.0     |    51.2    |
       Locus level:   100.0     |    57.8    |

     Matching intron chains:   65384
       Matching transcripts:   67853
              Matching loci:   39152

          Missed exons:       0/352731  (  0.0%)
           Novel exons:  266994/624279  ( 42.8%)
        Missed introns:       6/310704  (  0.0%)
         Novel introns:  206267/519190  ( 39.7%)
           Missed loci:       0/39152   (  0.0%)
            Novel loci:   28509/67496   ( 42.2%)

 Total union super-loci across all input datasets: 67496
132608 out of 132608 consensus transcripts written in Cvir-compare.annotated.gtf (0 discarded as redundant)


this is a huge decreseas in precision..

it's also interesting that there are MORE query mRNAs, MORE novel loci... I am very confused...

there just seems so be SO many more novel exons, which means thery're not in the reference, which may be bringing down the precision? because they are considered at false positive? are they false?


is there a difference with this annotation file? should i try gffcomapre with the other one??
ref_C_virginica-3.0_top_level.gff3

gffcompare -r ref_C_virginica-3.0_top_level.gff3 -o Cvir-top-level-compare 2020-Cvir-merged.gtf

gave same of this
67891 reference transcripts loaded.
  38 duplicate reference transcripts discarded.
  132608 query transfrags loaded.

Ok this gives absolutely the same statistics so there is not a problem with the reference anotation file
deleted those

http://seqanswers.com/forums/archive/index.php/t-66297.html
"The other ways to filter more of the transcripts are with the -f parameter just as mentioned before, or with the -F or -T parameters that filter out transcripts of very low abundance in the samples. We like filtering with -F and -T more than with the -f option, because -f filters transcripts that have a relative low abundance compared to the most abundant transcript in the bundle, even if sometimes the transcripts that are filtered out are highly expressed."

maybe I should be using the parameters in the merge step as well..

I also didn't use the -f parameter. moved those comparison files into a directory called Cvir-compare

First: re-do merge with -f parameter set and see what happens

`stringtie --merge -G GCF_002022765.2_C_virginica-3.0_genomic.gff -o 2020-Cvir-merged-f.gtf -f 0.1 merge-gtf-list.txt`

look at comparsion again

`gffcompare -r GCF_002022765.2_C_virginica-3.0_genomic.gff -o Cvir-compare-f 2020-Cvir-merged-f.gtf`

67891 reference transcripts loaded.
 38 duplicate reference transcripts discarded.
 112319 query transfrags loaded

 less query transfrags. now this is interesting, how are there less? Ohhhh there are less because the -f parameter is a proportion, so even though I had it as a proportion in the individual gtf files, this is a proportion accross all of them?

 ok so this increased the precision but not by a lot

 #= Summary for dataset: 2020-Cvir-merged-f.gtf
#     Query mRNAs :  112319 in   67362 loci  (102080 multi-exon transcripts)
#            (18285 multi-transcript loci, ~1.7 transcripts per locus)
# Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
# Super-loci w/ reference transcripts:    39018
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    60.2    |
        Exon level:   100.0     |    60.5    |
      Intron level:   100.0     |    62.4    |
Intron chain level:   100.0     |    64.1    |
  Transcript level:   100.0     |    60.4    |
       Locus level:   100.0     |    57.9    |

     Matching intron chains:   65384
       Matching transcripts:   67853
              Matching loci:   39152

          Missed exons:       0/352731  (  0.0%)
           Novel exons:  226528/583032  ( 38.9%)
        Missed introns:       6/310704  (  0.0%)
         Novel introns:  185925/498287  ( 37.3%)
           Missed loci:       0/39152   (  0.0%)
            Novel loci:   28344/67362   ( 42.1%)

 Total union super-loci across all input datasets: 67362
112319 out of 112319 consensus transcripts written in Cvir-compare-f.annotated.gtf (0 discarded as redundant)

Interesting that the novel loci is the same.. looks like novel exons and novel introns went down, but mostly exons.
Still the same number of novel loci, so maybe this didn't take care of the issue

Do I even know if this is the issue? Not really but it seems like too many new ones

stringtie man say -F and -T so I guess I should try those

-F <min_fpkm>	minimum input transcript FPKM to include in the merge (default: 0)
-T <min_tpm>	minimum input transcript TPM to include in the merge (default: 0)

I guess I should put 1 for each of these, these are very different metrics... hard to know what cut off to use but I think the should be the same

`stringtie --merge -G GCF_002022765.2_C_virginica-3.0_genomic.gff -o 2020-Cvir-merged-fFT.gtf -f 0.1 -F 1 -T 1 merge-gtf-list.txt`


`gffcompare -r GCF_002022765.2_C_virginica-3.0_genomic.gff -o Cvir-compare-fFT 2020-Cvir-merged-fFT.gtf`

67891 reference transcripts loaded.
  38 duplicate reference transcripts discarded.
  112319 query transfrags loaded.

ok so it is exactly the same. I will try with 5 for each of those parameters now?

`stringtie --merge -G GCF_002022765.2_C_virginica-3.0_genomic.gff -o 2020-Cvir-merged-fFT.gtf -f 0.1 -F 5 -T 5 merge-gtf-list.txt`


`gffcompare -r GCF_002022765.2_C_virginica-3.0_genomic.gff -o Cvir-compare-fFT 2020-Cvir-merged-fFT.gtf`

ok this wildly cut down
67891 reference transcripts loaded.
 38 duplicate reference transcripts discarded.
 85417 query transfrags loaded.


#= Summary for dataset: 2020-Cvir-merged-fFT.gtf
#     Query mRNAs :   85417 in   51155 loci  (81538 multi-exon transcripts)
#            (14145 multi-transcript loci, ~1.7 transcripts per locus)
# Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
# Super-loci w/ reference transcripts:    39108
#-----------------| Sensitivity | Precision  |
         Base level:   100.0     |    74.6    |
         Exon level:   100.0     |    74.4    |
       Intron level:   100.0     |    74.8    |
 Intron chain level:   100.0     |    80.2    |
   Transcript level:   100.0     |    79.4    |
        Locus level:   100.0     |    76.5    |

      Matching intron chains:   65384
        Matching transcripts:   67853
               Matching loci:   39152

           Missed exons:       0/352731  (  0.0%)
            Novel exons:  120058/474258  ( 25.3%)
         Missed introns:       0/310704  (  0.0%)
          Novel introns:  104083/415371  ( 25.1%)
            Missed loci:       0/39152   (  0.0%)
             Novel loci:   12047/51155   ( 23.5%)

  Total union super-loci across all input datasets: 51155
 85417 out of 85417 consensus transcripts written in Cvir-compare-fFT.annotated.gtf (0 discarded as redundant)

 what I'm thinking is precision is going up because novel loci is going down. And these are things with I guess low reads/counts so they might not be real. hard to say

 not sure if I should increase parameters again or not

### next steps
Re-run stringtie without strandness parameter
re-run stringtie with only strandedness parameter
re-run with only conservative flag
want to use -C with the one we go forward with?


1. w/out strandedness
`mkdir stringtie-no-rf`
`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/*.s.bam .`
`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/GCF_002022765.2_C_virginica-3.0_genomic.gff .`

new script for stringtie without the rf
going to use the stringent parameters I used before though
nano stringtie-no-rf.sh
```
#!/bin/bash
# In the same directory now with the BAM files and the annotation file link
F=/home/mschedl/Working-CASE-RNA/2020-New-Analysis/stringtie-no-rf

# StringTie to assemble transcripts for each sample with the annotation file
array1=($(ls $F/*.bam))
for i in ${array1[@]}; do
	stringtie -G $F/GCF_002022765.2_C_virginica-3.0_genomic.gff -f .1 -c 3 -s 7 -j 2 -o ${i}.gtf ${i}
	echo "${i}"
done
```
nohup bash stringtie-no-rf.sh


cd ..
mkdir stringtie-rf-only
`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/*.s.bam .`
`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/GCF_002022765.2_C_virginica-3.0_genomic.gff .`
nano stringtie-rf-only.sh
```
#!/bin/bash
# In the same directory now with the BAM files and the annotation file link
F=/home/mschedl/Working-CASE-RNA/2020-New-Analysis/stringtie-rf-only

# StringTie to assemble transcripts for each sample with the annotation file
array1=($(ls $F/*.bam))
for i in ${array1[@]}; do
	stringtie -G $F/GCF_002022765.2_C_virginica-3.0_genomic.gff --rf -o ${i}.gtf ${i}
	echo "${i}"
done
```
nohup bash stringtie-rf-only.sh

now to try different merging on these

cd stringtie-no-rf

`nano merge-gtf-list.txt`
CA_J06.F.trim.fq.gz.s.bam.gtf
CA_J08.F.trim.fq.gz.s.bam.gtf
CA_J11.F.trim.fq.gz.s.bam.gtf
CA_J18.F.trim.fq.gz.s.bam.gtf
CASE_J03.F.trim.fq.gz.s.bam.gtf
CASE_J09.F.trim.fq.gz.s.bam.gtf
CASE_J12.F.trim.fq.gz.s.bam.gtf
CASE_J13.F.trim.fq.gz.s.bam.gtf
CON_J02.F.trim.fq.gz.s.bam.gtf
CON_J05.F.trim.fq.gz.s.bam.gtf
CON_J10.F.trim.fq.gz.s.bam.gtf
SE_J01.F.trim.fq.gz.s.bam.gtf
SE_J04.F.trim.fq.gz.s.bam.gtf
SE_J07.F.trim.fq.gz.s.bam.gtf

Basic stringtie merge no parameters
`stringtie --merge -G GCF_002022765.2_C_virginica-3.0_genomic.gff -o Merged-no-rf.gtf merge-gtf-list.txt`
`gffcompare -r GCF_002022765.2_C_virginica-3.0_genomic.gff -o Compare-no-rf Merged-no-rf.gtf`
67891 reference transcripts loaded.
  38 duplicate reference transcripts discarded.
  132608 query transfrags loaded.

### Stringent stringtie parameters but no strandedness

#= Summary for dataset: Merged-no-rf.gtf
#     Query mRNAs :  132608 in   67496 loci  (122412 multi-exon transcripts)
#            (21190 multi-transcript loci, ~2.0 transcripts per locus)
# Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
# Super-loci w/ reference transcripts:    38987
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    59.1    |
        Exon level:   100.0     |    56.5    |
      Intron level:   100.0     |    59.8    |
Intron chain level:   100.0     |    53.4    |
  Transcript level:   100.0     |    51.2    |
       Locus level:   100.0     |    57.8    |

     Matching intron chains:   65384
       Matching transcripts:   67853
              Matching loci:   39152

          Missed exons:       0/352731  (  0.0%)
           Novel exons:  266994/624279  ( 42.8%)
        Missed introns:       6/310704  (  0.0%)
         Novel introns:  206267/519190  ( 39.7%)
           Missed loci:       0/39152   (  0.0%)
            Novel loci:   28509/67496   ( 42.2%)

 Total union super-loci across all input datasets: 67496
132608 out of 132608 consensus transcripts written in Compare-no-rf.annotated.gtf (0 discarded as redundant)

looks the same as the stringent stringtie with strandedness.... so the strandedness is not causing the increase in novel loci. This is strange

what about the only strandedness??

cd..
cd stringtie-rf-only/

`nano merge-gtf-list.txt`
CA_J06.F.trim.fq.gz.s.bam.gtf
CA_J08.F.trim.fq.gz.s.bam.gtf
CA_J11.F.trim.fq.gz.s.bam.gtf
CA_J18.F.trim.fq.gz.s.bam.gtf
CASE_J03.F.trim.fq.gz.s.bam.gtf
CASE_J09.F.trim.fq.gz.s.bam.gtf
CASE_J12.F.trim.fq.gz.s.bam.gtf
CASE_J13.F.trim.fq.gz.s.bam.gtf
CON_J02.F.trim.fq.gz.s.bam.gtf
CON_J05.F.trim.fq.gz.s.bam.gtf
CON_J10.F.trim.fq.gz.s.bam.gtf
SE_J01.F.trim.fq.gz.s.bam.gtf
SE_J04.F.trim.fq.gz.s.bam.gtf
SE_J07.F.trim.fq.gz.s.bam.gtf

Basic stringtie merge no parameters
`stringtie --merge -G GCF_002022765.2_C_virginica-3.0_genomic.gff -o Merged-rf-only.gtf merge-gtf-list.txt`
`gffcompare -r GCF_002022765.2_C_virginica-3.0_genomic.gff -o Compare-rf-only Merged-rf-only.gtf`
67891 reference transcripts loaded.
 38 duplicate reference transcripts discarded.
 135112 query transfrags loaded.

 already are more query transfrags, so maybe the stringent filters do something..

 #= Summary for dataset: Merged-rf-only.gtf
#     Query mRNAs :  135112 in   66482 loci  (125406 multi-exon transcripts)
#            (21748 multi-transcript loci, ~2.0 transcripts per locus)
# Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
# Super-loci w/ reference transcripts:    38969
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    59.1    |
        Exon level:   100.0     |    55.9    |
      Intron level:   100.0     |    59.3    |
Intron chain level:   100.0     |    52.1    |
  Transcript level:   100.0     |    50.2    |
       Locus level:   100.0     |    58.6    |

     Matching intron chains:   65384
       Matching transcripts:   67853
              Matching loci:   39152

          Missed exons:       0/352731  (  0.0%)
           Novel exons:  273999/631712  ( 43.4%)
        Missed introns:       6/310704  (  0.0%)
         Novel introns:  210653/523936  ( 40.2%)
           Missed loci:       0/39152   (  0.0%)
            Novel loci:   27513/66482   ( 41.4%)

 Total union super-loci across all input datasets: 66482
135112 out of 135112 consensus transcripts written in Compare-rf-only.annotated.gtf (0 discarded as redundant)

ok... this is not much better, percentage of novel loci went down a little but it may not have changed at all if there are more loci in general

what about if I use conservative flag? I guess at this point it doesn't seem to matter if I use rf or not

cd ..
mkdir conservative-stringtie

`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/*.s.bam .`
`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/GCF_002022765.2_C_virginica-3.0_genomic.gff .`

new script for stringtie without the rf
going to use the stringent parameters I used before though
nano stringtie-cons.sh
```
#!/bin/bash
# In the same directory now with the BAM files and the annotation file link
F=/home/mschedl/Working-CASE-RNA/2020-New-Analysis/conservative-stringtie

# StringTie to assemble transcripts for each sample with the annotation file
array1=($(ls $F/*.bam))
for i in ${array1[@]}; do
	stringtie -G $F/GCF_002022765.2_C_virginica-3.0_genomic.gff -t -c 1.5 -f 0.05 -o ${i}.gtf ${i}
	echo "${i}"
done
```
nohup bash stringtie-cons.sh

--conservative is not an option apparently
just use -t -c 1.5 -f 0.05 instead?


still seems like here the problem is not anything I am doing at this step?? Was doing the alignment stranded wrong?


`nano merge-gtf-list.txt`
CA_J06.F.trim.fq.gz.s.bam.gtf
CA_J08.F.trim.fq.gz.s.bam.gtf
CA_J11.F.trim.fq.gz.s.bam.gtf
CA_J18.F.trim.fq.gz.s.bam.gtf
CASE_J03.F.trim.fq.gz.s.bam.gtf
CASE_J09.F.trim.fq.gz.s.bam.gtf
CASE_J12.F.trim.fq.gz.s.bam.gtf
CASE_J13.F.trim.fq.gz.s.bam.gtf
CON_J02.F.trim.fq.gz.s.bam.gtf
CON_J05.F.trim.fq.gz.s.bam.gtf
CON_J10.F.trim.fq.gz.s.bam.gtf
SE_J01.F.trim.fq.gz.s.bam.gtf
SE_J04.F.trim.fq.gz.s.bam.gtf
SE_J07.F.trim.fq.gz.s.bam.gtf

Basic stringtie merge no parameters
`stringtie --merge -G GCF_002022765.2_C_virginica-3.0_genomic.gff -o Merged-cons.gtf merge-gtf-list.txt`
`gffcompare -r GCF_002022765.2_C_virginica-3.0_genomic.gff -o Compare-cons Merged-cons.gtf`

67891 reference transcripts loaded.
38 duplicate reference transcripts discarded.
137896 query transfrags loaded.

#= Summary for dataset: Merged-cons.gtf
#     Query mRNAs :  137896 in   65749 loci  (128562 multi-exon transcripts)
#            (22107 multi-transcript loci, ~2.1 transcripts per locus)
# Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
# Super-loci w/ reference transcripts:    38925
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    58.5    |
        Exon level:   100.0     |    55.2    |
      Intron level:   100.0     |    58.8    |
Intron chain level:   100.0     |    50.9    |
  Transcript level:   100.0     |    49.2    |
       Locus level:   100.0     |    59.2    |

     Matching intron chains:   65384
       Matching transcripts:   67853
              Matching loci:   39152

          Missed exons:       0/352731  (  0.0%)
           Novel exons:  280123/638736  ( 43.9%)
        Missed introns:       6/310704  (  0.0%)
         Novel introns:  214503/528064  ( 40.6%)
           Missed loci:       0/39152   (  0.0%)
            Novel loci:   26824/65749   ( 40.8%)

 Total union super-loci across all input datasets: 65749
137896 out of 137896 consensus transcripts written in Compare-cons.annotated.gtf (0 discarded as redundant)


run with only stringent merge parameters

mkdir merge-only

`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/*.s.bam .`
`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/GCF_002022765.2_C_virginica-3.0_genomic.gff .`

nano stringtie.sh
```
#!/bin/bash
# In the same directory now with the BAM files and the annotation file link
F=/home/mschedl/Working-CASE-RNA/2020-New-Analysis/merge-only

# StringTie to assemble transcripts for each sample with the annotation file
array1=($(ls $F/*.bam))
for i in ${array1[@]}; do
	stringtie -G $F/GCF_002022765.2_C_virginica-3.0_genomic.gff  -o ${i}.gtf ${i}
	echo "${i}"
done
```
`nano merge-gtf-list.txt`
CA_J06.F.trim.fq.gz.s.bam.gtf
CA_J08.F.trim.fq.gz.s.bam.gtf
CA_J11.F.trim.fq.gz.s.bam.gtf
CA_J18.F.trim.fq.gz.s.bam.gtf
CASE_J03.F.trim.fq.gz.s.bam.gtf
CASE_J09.F.trim.fq.gz.s.bam.gtf
CASE_J12.F.trim.fq.gz.s.bam.gtf
CASE_J13.F.trim.fq.gz.s.bam.gtf
CON_J02.F.trim.fq.gz.s.bam.gtf
CON_J05.F.trim.fq.gz.s.bam.gtf
CON_J10.F.trim.fq.gz.s.bam.gtf
SE_J01.F.trim.fq.gz.s.bam.gtf
SE_J04.F.trim.fq.gz.s.bam.gtf
SE_J07.F.trim.fq.gz.s.bam.gtf

stringtie merge stringent parameters
`stringtie --merge -G GCF_002022765.2_C_virginica-3.0_genomic.gff -f 0.1 -F 5 -T 5 -o Merged-stringent.gtf merge-gtf-list.txt`
`gffcompare -r GCF_002022765.2_C_virginica-3.0_genomic.gff -o Compare-merge-stringent Merged-stringent.gtf`

67891 reference transcripts loaded.
38 duplicate reference transcripts discarded.
85265 query transfrags loaded.

#= Summary for dataset: Merged-stringent.gtf
#     Query mRNAs :   85265 in   50818 loci  (81453 multi-exon transcripts)
#            (14198 multi-transcript loci, ~1.7 transcripts per locus)
# Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
# Super-loci w/ reference transcripts:    39105
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    74.8    |
        Exon level:   100.0     |    74.5    |
      Intron level:   100.0     |    75.0    |
Intron chain level:   100.0     |    80.3    |
  Transcript level:   100.0     |    79.6    |
       Locus level:   100.0     |    77.0    |

     Matching intron chains:   65384
       Matching transcripts:   67853
              Matching loci:   39152

          Missed exons:       0/352731  (  0.0%)
           Novel exons:  118989/473290  ( 25.1%)
        Missed introns:       0/310704  (  0.0%)
         Novel introns:  103155/414521  ( 24.9%)
           Missed loci:       0/39152   (  0.0%)
            Novel loci:   11713/50818   ( 23.0%)

 Total union super-loci across all input datasets: 50818
85265 out of 85265 consensus transcripts written in Compare-merge-stringent.annotated.gtf (0 discarded as redundant)
Compare-merge-stringent.stats (END)




also want to do no stringtie no merge parameters

mkdir no-parameters

`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/*.s.bam .`
`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/GCF_002022765.2_C_virginica-3.0_genomic.gff .`

nano stringtie.sh
```
#!/bin/bash
# In the same directory now with the BAM files and the annotation file link
F=/home/mschedl/Working-CASE-RNA/2020-New-Analysis/no-parameters

# StringTie to assemble transcripts for each sample with the annotation file
array1=($(ls $F/*.bam))
for i in ${array1[@]}; do
	stringtie -G $F/GCF_002022765.2_C_virginica-3.0_genomic.gff  -o ${i}.gtf ${i}
	echo "${i}"
done
```
`nano merge-gtf-list.txt`
CA_J06.F.trim.fq.gz.s.bam.gtf
CA_J08.F.trim.fq.gz.s.bam.gtf
CA_J11.F.trim.fq.gz.s.bam.gtf
CA_J18.F.trim.fq.gz.s.bam.gtf
CASE_J03.F.trim.fq.gz.s.bam.gtf
CASE_J09.F.trim.fq.gz.s.bam.gtf
CASE_J12.F.trim.fq.gz.s.bam.gtf
CASE_J13.F.trim.fq.gz.s.bam.gtf
CON_J02.F.trim.fq.gz.s.bam.gtf
CON_J05.F.trim.fq.gz.s.bam.gtf
CON_J10.F.trim.fq.gz.s.bam.gtf
SE_J01.F.trim.fq.gz.s.bam.gtf
SE_J04.F.trim.fq.gz.s.bam.gtf
SE_J07.F.trim.fq.gz.s.bam.gtf

stringtie merge no parameters
`stringtie --merge -G GCF_002022765.2_C_virginica-3.0_genomic.gff  -o Merged.gtf merge-gtf-list.txt`
`gffcompare -r GCF_002022765.2_C_virginica-3.0_genomic.gff -o Compare-merge Merged.gtf`
67891 reference transcripts loaded.
38 duplicate reference transcripts discarded.
135112 query transfrags loaded.

#= Summary for dataset: Merged.gtf
#     Query mRNAs :  135112 in   66482 loci  (125406 multi-exon transcripts)
#            (21748 multi-transcript loci, ~2.0 transcripts per locus)
# Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
# Super-loci w/ reference transcripts:    38969
#-----------------| Sensitivity | Precision  |
        Base level:   100.0     |    59.1    |
        Exon level:   100.0     |    55.9    |
      Intron level:   100.0     |    59.3    |
Intron chain level:   100.0     |    52.1    |
  Transcript level:   100.0     |    50.2    |
       Locus level:   100.0     |    58.6    |

     Matching intron chains:   65384
       Matching transcripts:   67853
              Matching loci:   39152

          Missed exons:       0/352731  (  0.0%)
           Novel exons:  273999/631712  ( 43.4%)
        Missed introns:       6/310704  (  0.0%)
         Novel introns:  210653/523936  ( 40.2%)
           Missed loci:       0/39152   (  0.0%)
            Novel loci:   27513/66482   ( 41.4%)

 Total union super-loci across all input datasets: 66482
135112 out of 135112 consensus transcripts written in Compare-merge.annotated.gtf (0 discarded as redundant)

also want to try with -e flag to limit

mkdir ref-transcripts-only

`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/*.s.bam .`
`ln -s /home/mschedl/Working-CASE-RNA/2020-New-Analysis/GCF_002022765.2_C_virginica-3.0_genomic.gff .`

nano stringtie-e.sh
```
#!/bin/bash
# In the same directory now with the BAM files and the annotation file link
F=/home/mschedl/Working-CASE-RNA/2020-New-Analysis/ref-transcripts-only

# StringTie to assemble transcripts for each sample with the annotation file
array1=($(ls $F/*.bam))
for i in ${array1[@]}; do
	stringtie -G $F/GCF_002022765.2_C_virginica-3.0_genomic.gff -e -o ${i}.gtf ${i}
	echo "${i}"
done
```
`nano merge-gtf-list.txt`
CA_J06.F.trim.fq.gz.s.bam.gtf
CA_J08.F.trim.fq.gz.s.bam.gtf
CA_J11.F.trim.fq.gz.s.bam.gtf
CA_J18.F.trim.fq.gz.s.bam.gtf
CASE_J03.F.trim.fq.gz.s.bam.gtf
CASE_J09.F.trim.fq.gz.s.bam.gtf
CASE_J12.F.trim.fq.gz.s.bam.gtf
CASE_J13.F.trim.fq.gz.s.bam.gtf
CON_J02.F.trim.fq.gz.s.bam.gtf
CON_J05.F.trim.fq.gz.s.bam.gtf
CON_J10.F.trim.fq.gz.s.bam.gtf
SE_J01.F.trim.fq.gz.s.bam.gtf
SE_J04.F.trim.fq.gz.s.bam.gtf
SE_J07.F.trim.fq.gz.s.bam.gtf

stringtie merge with -e

doesn't really need to be merged but need to see how many loci
`stringtie --merge -G GCF_002022765.2_C_virginica-3.0_genomic.gff -o Merged-e.gtf merge-gtf-list.txt`
`gffcompare -r GCF_002022765.2_C_virginica-3.0_genomic.gff -o Compare-merge-e Merged-e.gtf`

67891 reference transcripts loaded.
  38 duplicate reference transcripts discarded.
  67883 query transfrags loaded.

  #= Summary for dataset: Merged-e.gtf
  #     Query mRNAs :   67883 in   39152 loci  (65414 multi-exon transcripts)
  #            (11433 multi-transcript loci, ~1.7 transcripts per locus)
  # Reference mRNAs :   67853 in   39152 loci  (65384 multi-exon)
  # Super-loci w/ reference transcripts:    39152
  #-----------------| Sensitivity | Precision  |
          Base level:   100.0     |   100.0    |
          Exon level:   100.0     |   100.0    |
        Intron level:   100.0     |   100.0    |
  Intron chain level:   100.0     |   100.0    |
    Transcript level:   100.0     |   100.0    |
         Locus level:   100.0     |   100.0    |

       Matching intron chains:   65384
         Matching transcripts:   67853
                Matching loci:   39152

            Missed exons:       0/352731  (  0.0%)
             Novel exons:       0/352757  (  0.0%)
          Missed introns:       0/310704  (  0.0%)
           Novel introns:       0/310704  (  0.0%)
             Missed loci:       0/39152   (  0.0%)
              Novel loci:       0/39152   (  0.0%)

   Total union super-loci across all input datasets: 39152
  67883 out of 67883 consensus transcripts written in Compare-merge-e.annotated.gtf (0 discarded as redundant)
