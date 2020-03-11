### Getting the LOC names that match the MSTRG and gene names that are in the greenyellow module

**Take the csv from R and get rid of all the things I don't want in it. First make it tab separated instead of commas, then only select the second column, then remove the ""**

```
sed -e 's/,/  /g'  GYModGeneNames.csv | cut -f2 | sed -e 's/"//g' > GYModGeneNames
```
Then I went in with nano to get rid of the words "v1" that was just in the first row. There are 57 genes. Now it looks like this:

MSTRG.13878  
MSTRG.25557  
MSTRG.13690  
MSTRG.16274  
MSTRG.27893  
MSTRG.2530  
MSTRG.22163  
MSTRG.20077  
MSTRG.15392  
MSTRG.8152  
MSTRG.2984  
MSTRG.2767  
MSTRG.11988  
MSTRG.766  
MSTRG.10053  
MSTRG.7790  
MSTRG.2705  
MSTRG.25599  
MSTRG.14393  
MSTRG.16437  
MSTRG.12661  
MSTRG.5955  
MSTRG.17269  
MSTRG.3135  
MSTRG.23860  
MSTRG.32832  
MSTRG.1768  
gene9995  
MSTRG.435  
MSTRG.31251  
gene29006  
MSTRG.6399  
MSTRG.28055  
MSTRG.2399  
MSTRG.16456  
MSTRG.2457  
MSTRG.16584  
gene15520  
MSTRG.31917  
gene15540  
MSTRG.16180  
MSTRG.19568  
MSTRG.2908  
MSTRG.25289  
gene17224  
MSTRG.16547  
gene13664  
gene30051  
MSTRG.15756  
MSTRG.14902  
MSTRG.25273  
MSTRG.8362  
gene14559  
MSTRG.14531  
MSTRG.26795  
MSTRG.29582  
MSTRG.28511  

**Now search the merged transcript file for each of these gene names. The -f flag (lowercase) means to obtain the pattern from the file. This gives all the lines that have those gene names.**
```
grep -f GYModGeneNames C_Vir_ST_merged_NO_A_test.gtf > GYMODS
```

Ok there is a lot of information in these lines, and I don't need all of it. The [columns of .gtf files](https://useast.ensembl.org/info/website/upload/gff.html) are these:

1 chromosome name  
2 source of file  
3 feature (transcript or exon)  
4 start position  
5 end position  
6 score  
7 strand (+ is forward, - is reverse)  
8 frame (not sure how this one applies)  
9 attributes separated by semicolons  


The 9th column contains what I want. But things are tab delimited and then semicolon delimited and there are "" everywhere and a lot of duplicated entries because each are from a different exon of the same gene. Below is step by step what I did to get only the LOCxxxx information only. The bottom is a long one line of code piped that does all the same things in one go.

**Take only the 9th column of that file and make it into a new one**
```
cut -f9 GYMODS > GYMODS9
```

Ok how to separate all this out... there are so many duplicates because there can be like 18 different exons for one gene and it gives you them all. The easiest thing to look for first is probably the gene name. Not all have them but they're like LOCxxxxx and that is something I can search in the annotation file _ref_C_virginica-3.0_top_level.gff3_


**Make these tab delimited instead of with a semicolon. To make a tab you have to type this into the terminal and press ctl v then tab**

```
sed -e 's/;/  /g'  GYMODS9 > GYMODStab
```

So now gene_name should be in the 4th column (tab delimitated).  
**Separate out the 4th column**
```
cut -f4 GYMODStab > GYMODSgene_name
```

Now the file has the words "gene_name" or "ref_gene_id" in it, but I don't need those  
**Separate the space between "gene_name" and "LOCxxx" into two tab delimited columns again and only take the column with the LOCxxx**
```
sed -e 's/ /  /g'  GYMODSgene_name > GYMODSgene

cut -f3 GYMODSgene > GYMODSLOC
```
Now the only think left to remove are the duplicates and the quotation marks.

**First only accept unique entries**
```
uniq GYMODSLOC > GYMODSLOCuniq
```
**Replace quotation marks with nothing**
```
sed -e 's/"//g' GYMODSLOCuniq > GYMODSLOCuniq2
```
I'm not sure if the lines with nothing in them will cause me a problem, but I bet there is a way to get rid of them.  
**[Remove empty lines](https://serverfault.com/questions/252921/how-to-remove-empty-blank-lines-from-a-file-in-unix-including-spaces) because ^$ means empty line in bash**
```
sed '/^$/d' GYMODSLOCuniq2 > GYMODSLOCuniq3
```
**Instead of removing the gene# I can find the LOC. For some reason the other way didn't work. And I think I have to use awk if I want to use * for the numbers**
```
awk '/LOC*/' GYMODSLOCuniq3 > GYMODSLOConly
```
**How many are there?**
```
 wc -l GYMODSLOConly
```
 257

 I don't know if that's weird or not.... there were only 57 genes from the module???


**This is every step in one. It produces the same file as GYMODSLOConly. You just have to make sure to type _ctrl v tab_ to make the tab when you are in the terminal.**

 ```
 cut -f9 GYMODS | sed -e 's/;/  /g' | cut -f4 | sed -e 's/ /  /g' | cut -f3 | uniq | sed -e 's/"//g' | sed '/^$/d' | awk '/LOC*/' > LOCSONLY
 ```
![meme](images/the_limit.png)
