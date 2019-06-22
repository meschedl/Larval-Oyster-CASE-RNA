
# CASE RNASeq Analysis : QC, Read Trimming, Read Alignment, and Transcript Assembly

**Author: Maggie Schedl**

_All analysis was done on our lab shared server, KITT, made by [J. Puritz](https://github.com/jpuritz)_  

Lines of code with this symbol `¯\_(ツ)_/¯` are used in my personal computer terminal window

Programs Installed/Needed for this Project:  
- HISAT2 `conda install hisat2`
- StringTie `conda install stringtie`
- gffcomare `conda install gffcompare`
- fastp, fastQC, multiqc
- samtools

File Naming and Information:
- **cvir_edited.fa**: Eastern Oyster genome, I received this from [Erin Roberts](https://github.com/erinroberts), from what I can tell it is the exact same as the [_C. virginica_ genome on NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=Crassostrea%20virginica) with just a space in the header line of each chromosome removed
- **ref_C_virginica-3.0_top_level.gff3**: annotation file for the Eastern Oyster, I also got this from Erin, it has the matching header line convention to work with the above genome
- **prepDE.py**: python script for converting read count information from StringTie into a matrix fo DESeq2, full code [here](http://ccb.jhu.edu/software/stringtie/dl/prepDE.py)
- **CA_J06**: this is an example of the naming convention of the samples, CA refers to the treatment (coastal acidification), and the J number refers to the replicate jar

I used Erin Robert's [RNASeq pipeline](https://github.com/erinroberts/apoptosis_data_pipeline/blob/master/Streamlined%20Pipeline%20Tutorial/Apoptosis_Pipeline_Tutorial_with_Reduced_Dataset.md) or Kevin Wong's [project](https://github.com/jpuritz/BIO_594_2018/blob/master/FinalAssignment/KevinWong_FinalAssignment/P.dam_DE_Analysis.md) in [J. Puritz's](https://github.com/jpuritz) Bio594 2018, for reference.  If other resources were used they should be linked in this markdown. If anything is not linked properly, not-sourced, or looks missing please contact me at meschedl@uri.edu.

-------

#### Quality Control and Read Trimming

Reads were already de-multiplexed and assigned to each individual sample by [J. Puritz](https://github.com/jpuritz), and linked to a directory called CASE_RNA.

 I made a second directory to work in, and then linked in the files again.
```
mkdir Working-CASE-RNA
cd Working-CASE-RNA

ln -s /home/mschedl/CASE_RNA/* .
```

The first thing I did was look at the read counts for each file. I used a code from [this website](http://www.sixthresearcher.com/list-of-helpful-linux-commands-to-process-fastq-files-from-ngs-experiments/) and made it into a for-loop that went through all the files.

```
for fq in *.fq.gz
> do
>	echo $fq
> zcat $fq | echo $((`wc -l`/4))
> done
```
**Output:**  
CA_J06.F.fq.gz  
20445176  
CA_J06.R.fq.gz  
20445176  
CA_J08.F.fq.gz  
21746189  
CA_J08.R.fq.gz  
21746189  
CA_J11.F.fq.gz  
25550864  
CA_J11.R.fq.gz  
25550864  
CA_J18.F.fq.gz  
37263541  
CA_J18.R.fq.gz  
37263541  
CASE_J03.F.fq.gz  
26925142  
CASE_J03.R.fq.gz  
26925142  
CASE_J09.F.fq.gz  
31720810  
CASE_J09.R.fq.gz  
31720810  
CASE_J12.F.fq.gz  
24582739  
CASE_J12.R.fq.gz  
24582739  
CASE_J13.F.fq.gz  
36132924  
CASE_J13.R.fq.gz  
36132924  
CON_J02.F.fq.gz  
28850301  
CON_J02.R.fq.gz  
28850301  
CON_J05.F.fq.gz  
27446573  
CON_J05.R.fq.gz  
27446573  
CON_J10.F.fq.gz  
35291136  
CON_J10.R.fq.gz  
35291136  
SE_J01.F.fq.gz  
28376966  
SE_J01.R.fq.gz  
28376966  
SE_J04.F.fq.gz  
24827716  
SE_J04.R.fq.gz  
24827716  
SE_J07.F.fq.gz  
30894132  
SE_J07.R.fq.gz  
30894132  

Then I looked at the quality of the reads in each file. I like the look of the reports that [MultiQC](https://multiqc.info/) makes, and to make a MulitQC report, you need to use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

```
mkdir fastqc-before

cd fastqc-before/

fastqc ../*fq.gz

mv *fastqc.* fastqc-before/

cd fastq-before/

multiqc .

```
Then the report.html can be copied from KITT to my computer, so I can then look at it in a web browser.

In a terminal window:  
```
¯\_(ツ)_/¯ scp -P zzzz mschedl@KITT.uri.edu:/home/mschedl/Working-CASE-RNA/fastqc-before/multiqc_report.html /Users/maggieschedl/Desktop/URI/Classes/Puritz/CASE-RNA/outputs
```

![image1](https://raw.githubusercontent.com/jpuritz/BIO_594_2019/master/Final_Assignment/Maggie_Final_Project/images/before-qual.png)  
These are the mean quality scores across all reads, they look really good (all above 30). There is a little dip in the beginning of some though.

![image2](https://raw.githubusercontent.com/jpuritz/BIO_594_2019/master/Final_Assignment/Maggie_Final_Project/images/before-per-seq-qual.png)  
These are the per-sequence quality scores. Again, most look pretty good, where a lot of them are above 30 or at 40. But there is a tail that goes low, so I would like to shrink that tail.

![image3](https://raw.githubusercontent.com/jpuritz/BIO_594_2019/master/Final_Assignment/Maggie_Final_Project/images/GC.png)  
This is an image of the GC content of all the reads. In theory it should look like a normal distribution. This is not very normally distributed, indicating some sort of bias or contamination.

**Results Specific to RNASeq Data**

![image4](https://raw.githubusercontent.com/jpuritz/BIO_594_2019/master/Final_Assignment/Maggie_Final_Project/images/seq-counts.png)  
This shows the number of unique and duplicate reads for each sample. Because this is RNASeq data, there should be a lot of duplicate reads.

![image5](https://raw.githubusercontent.com/jpuritz/BIO_594_2019/master/Final_Assignment/Maggie_Final_Project/images/over-rep-seq.png)  
Same with the over-represented sequences. This is to be expected with RNASeq data, but not good with other types of data.

----

I used [fastp](https://github.com/OpenGene/fastp) to trim the reads some, especially to get the per-sequance quality score up. The -F and -f flags trim the front of both reads. The -q flag is for the phred quality value that a base is qualified. The -u flag is how many percents of bases are allowed to be unqualified. There are many options for trimming and cleaning of reads with fastp, these may have not been the best way to do this.

```
fastp --in1 CA_J06.F.fq.gz --in2 CA_J06.R.fq.gz --out1 CA_J06.F.trim.fq.gz --out2 CA_J06.R.trim.fq.gz  -f 5 -F 5 -q 15 -u 50 -j CA_J06.json -h CA_J06.html
```

The same parameters were applied to each sample. Afterwards, another MultiQC report was generated.

![image6](https://raw.githubusercontent.com/jpuritz/BIO_594_2019/master/Final_Assignment/Maggie_Final_Project/images/after-qual.png)  
Importantly, the quality scores are still high.

![image7](https://raw.githubusercontent.com/jpuritz/BIO_594_2019/master/Final_Assignment/Maggie_Final_Project/images/after-per-seq-qual.png)  
And the tail of that distribution scooted up a little.

----
#### Alignment to the Eastern Oyster Genome


The new fasta files generated with fastp can now be mapped to the Eastern Oyster genome with [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml), that link is to the manual, which I found most informative.

I made a new directory to work in as to not get confused by my files. Then I linked the trimmed fq.gz files into that directory, then I linked in the genome which lives in a storage directory.

```
mkdir histat
cd histat   

ln -s /home/mschedl/Working-CASE-RNA/histat/*.trim.fq.gz .
ln -s /RAID_STORAGE2/mschedl/RNA-CASE/cvir_edited.fa .
```

Then I used a script to first use the genome to make an index for hisat, then align each .trim.fq.gz file to the genome. This takes both the forward and reverse reads and makes a sam file of the alignment. Then, because sam files are HUGE, and StringTie needs a sorted bam file anyways, the script sorts the sam file into a bam file, and then removes the sam file. The --dta flag is very important, the StringTie website has this phrase "NOTE: be sure to run HISAT2 with the --dta option for alignment, or your results will suffer." That flag makes the bam files configured in a way that works well with the next program, StringTie.


`nano CASE-HISAT2.sh`

```
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
```

The output of this step is a sorted bam file, that has been aligned to the _C. virginica_ genome. These files are then used with StringTie to generate annotation files with transcript abundances.  
**Alignment rates were between 80-85%.** The entire output from HISAT2 can be found in the outputs folder.

----
#### Transcript Assembly

I created a new directory stringtie to make things clearer, and then linked in the annotation file which lives in storage. And I linked the bam files from the directory above to this one.

```
mkdir stringtie
cd stringtie
ln -s /RAID_STORAGE2/mschedl/RNA-CASE/ref_C_virginica-3.0_top_level.gff3 .
ln -s /home/mschedl/Working-CASE-RNA/histat/*.s.bam .

```
Originally, I had tried to run StringTie with the annotation file downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=Crassostrea%20virginica), which did not work, but that may have been for a variety of reasons at the time.


Note that I did not do all of StringTie in one script because of some problems with the merge step, as well as I had to re-do some steps in this process (see below).

Via the [StringTie manual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual):
 - Run each sample with the -G flag if a reference annotation is available
 - Run StringTie with --merge in order to generate a non-redundant set of transcripts observed in all the samples assembled previously
 - For each sample, run StringTie using the -B/-b and -e options in order to estimate transcript abundances. The -e option is not required but recommended for this run in order to produce more accurate abundance estimations of the input transcripts. Each StringTie run in this step will take as input the sorted read alignments (BAM file) obtained in step 1 for the corresponding sample and the -G option with the merged transcripts (GTF file) generated by stringtie --merge in step 3. Please note that this is the only case where the -G option is not used with a reference annotation, but with the global, merged set of transcripts as observed across all samples. _Note that I did not use the -B or -b flags because those are for if you plan on using the program ballgown_



First StringTie uses the annotation file to make a .gtf file containing transcript abundances for each sample. It uses the sorted bam file from the previous step as the input, and then with the -G flag, uses the Eastern Oyster annotation file as a reference.

`nano CASE-Stringtie-1.sh`

```
#!/bin/bash
# In the same directory now with the BAM files and the annotation file link
F=/home/mschedl/Working-CASE-RNA/histat/stringtie/

# StringTie to assemble transcripts for each sample with the annotation file
array1=($(ls $F/*.bam))
for i in ${array1[@]}; do
	stringtie -G $F/ref_C_virginica-3.0_top_level.gff3 -o ${i}.gtf ${i}
	echo "${i}"
done
```
Some notes about the above script, I didn't use the -l flag here. You use it to specify the name of output transcripts in each individual .gtf file. It ends up not mattering if you do because the merged .gtf file uses the default labels (MSTRG) whether you label in a previous step or not.  

 **This is where the MSTRG label comes from, all transcripts get a MSTRG label, some have a gene # associated with it in the ref_gene_id column of the final .gtf files if there was one, but all new transcripts have only the MSTRG label.**  

**Not** including the -e flag allows for novel transcript discovery because it "Limits the processing of read alignments to only estimate and output the assembled transcripts matching the reference transcripts given with the -G option." This was very important, and many of the transcripts are new. This may be because the reference annotation file is lacking in transcripts specific to larvae.

That created .gtf files for each sample, these will be merged into a single .gtf file that is non-redundant, and contains the transcripts from ALL of the files.

Here I originally made the merged .gtf file incorrectly. As I found out, do NOT use the -A flag when running the --merge option for StringTie. For some reason this causes the merged .gtf to miss exons and introns, which are not missed if you run this step as below. See the two gffcompare folders to see the stats of how the merged files compared to the reference. I now ran things in a directory called restring.

This step requires a .txt file of a list of all the sample names with the full path to each file on each line.

```
stringtie --merge -G ref_C_virginica-3.0_top_level.gff3 -o C_Vir_ST_merged_NO_A_test.gtf mergelist.txt
```
This made one merged .gtf file, C_Vir_ST_merged_NO_A_test.gtf. Notice the naming because I was testing if it was right. It was.



Then, I made a second script that uses the function gffcomare to give me some stats on how the reference annotation file and the merged annotation file compare. And it re-runs StringTie with the -e option and used the merged annotation file as the reference for assembly to the original bam files.

`nano CASE-Stringtie-2.sh`

```

#!/bin/bash
# In the same directory now with the BAM files and the annotation file link
F=/home/mschedl/Working-CASE-RNA/histat/stringtie/restring

# Want to use those bam files again to RE-estimate the transcript abundances
array1=($(ls $F/*.bam))

# gffcompare to compare how transcripts compare to reference annotation for the No A merged file

gffcompare -r $F/ref_C_virginica-3.0_top_level.gff3 -G -o c_vir_merged_compare_No_A C_Vir_ST_merged_NO_A_test.gtf
	# -o specifies prefix to use for output files
	# -r followed by the annotation file to use as a reference
 	# merged.annotation.gtf tells you how well the predicted transcripts track to the reference annotation file
 	# merged.stats file shows the sensitivity and precision statistics and total number for different features (genes, exons, transcripts)

#compare the merged file with the A
gffcompare -r $F/ref_C_virginica-3.0_top_level.gff3 -G -o c_vir_merged_compare_With_A C_Vir_ST_merged_test.gtf

#Re-estimate transcript abundance after merge step using the No A file
	for i in ${array1[@]}; do
		stringtie -e -G $F/C_Vir_ST_merged_NO_A_test.gtf -o $(echo ${i}|sed "s/\..*//").merge.gtf ${i}
		echo "${i}"
	done
	# input here is the original set of alignment bam files
	# here -G refers to the merged GTF files
	# -e creates more accurate abundance estimations with input transcripts, needed when converting to DESeq2 tables, says in the manual that this is recommended
echo "DONE" $(date)
```


The re-estimated .gtf files are then the input for the script prepDE.py. It's a really long script so I'm not going to write it out here but I just nano-d a file called prepDE.py and pasted [this code](https://ccb.jhu.edu/software/stringtie/dl/prepDE.py) into it.  
Then I used `chmod u+x prepDE.py` to activate it.

**Important** this script requires version 2.7 of python to run (it will not be able to read your list of files if you try using a different one). To check your version of python type `python --version`  
KITT uses version 3.6.8 by default, but you can create an environment to run a different version

```
conda create -n python27 python=2.7 anaconda

conda activate python27

python --version

```
**Python 2.7.15**

The script needs a text file with the path to the merged.gtf files and the name of the samples as the input. A simple script can do that for you.

`nano prepDEtest.sh`
```
#!/bin/bash

F=/home/mschedl/Working-CASE-RNA/histat/stringtie/restring


array2=($(ls *merge.gtf))

for i in ${array2[@]}; do
    echo "$(echo ${i}|sed "s/\..*//") $F/${i}" >> sample_list.txt
done
```
The sample_list.txt now has all the names and file paths. Then I ran the prepDE script.

```
python prepDE.py -i sample_lst.txt
```

This outputs a transcript_count_matrix.csv and a gene_count_matrix.csv

In excel I made a .csv file with the metadata/treatment information for each sample. It looks like this:

| sample   | treatment | library | extraction |
|----------|-----------|---------|------------|
| CASE_J03 | CASE      | three   | two        |
| CASE_J09 | CASE      | four    | two        |
| CASE_J12 | CASE      | two     | three      |
| CASE_J13 | CASE      | two     | three      |
| CA_J06   | CA        | three   | two        |
| CA_J08   | CA        | one     | two        |
| CA_J11   | CA        | four    | three      |
| CA_J18   | CA        | two     | three      |
| CON_J02  | CON       | three   | one        |
| CON_J05  | CON       | one     | two        |
| CON_J10  | CON       | four    | two        |
| SE_J01   | SE        | one     | one        |
| SE_J04   | SE        | four    | three      |
| SE_J07   | SE        | three   | two        |

Then I uploaded it to KITT so I could use it in the same directory as the count matrixes.

```
¯\_(ツ)_/¯ scp -P zzzz /Users/maggieschedl/Desktop/treatment_data.csv mschedl@KITT.uri.edu:/home/mschedl/Working-CASE-RNA/histat/stringtie/restring
```

Just the gene_count_matrix.csv and the treatment_data.csv are needed for DESeq2. See the Rmarkdown file that is set for GitHub view for all Differential Expression analysis!
