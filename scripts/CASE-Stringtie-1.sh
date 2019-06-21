```
#!/bin/bash
# In the same directory now with the BAM files and the annotation file link
F=/home/mschedl/Working-CASE-RNA/histat/stringtie

# StringTie to assemble transcripts for each sample with the annotation file
array1=($(ls $F/*.bam))
for i in ${array1[@]}; do
	stringtie -G $F/ref_C_virginica-3.0_top_level.gff3 -o ${i}.gtf ${i}
	echo "${i}"
done
```
