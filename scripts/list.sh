#!/bin/bash

F=/home/mschedl/Working-CASE-RNA/histat/stringtie/restring


array2=($(ls *NO.A.merge.gtf))

for i in ${array2[@]}; do
    echo "$(echo ${i}|sed "s/\..*//") $F/${i}" >> sample_list.txt
done
