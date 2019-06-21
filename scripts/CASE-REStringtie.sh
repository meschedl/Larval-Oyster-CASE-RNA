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
