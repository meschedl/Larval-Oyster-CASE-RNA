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

**Now search the merged transcript file for each of these gene names. The -f flag (lowercase) means to obtain the pattern from the file. The -w flag makes sure it's looking for the whole "word" and doesn't find any that are partially what I'm looking for. This gives all the lines that have those gene names.**
```
grep -w -f GYModGeneNames C_Vir_ST_merged_NO_A_test.gtf > GYMODS
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
sort -u GYMODSLOC | uniq > GYMODSLOCuniq
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
 47


**This is every step in one. It produces the same file as GYMODSLOConly. You just have to make sure to type _ctrl v tab_ to make the tab when you are in the terminal.**

 ```
 cut -f9 GYMODS | sed -e 's/;/  /g' | cut -f4 | sed -e 's/ /  /g' | cut -f3 | sort -u | uniq | sed -e 's/"//g' | sed '/^$/d' | awk '/LOC*/' > LOCSONLY
 ```
![meme](/images/the_limit.png)

LOC111099188  
LOC111099966  
LOC111100038  
LOC111100488  
LOC111102714  
LOC111102774  
LOC111104759  
LOC111105233  
LOC111106005  
LOC111106900  
LOC111107717  
LOC111108895  
LOC111110251  
LOC111110310  
LOC111110588  
LOC111112253  
LOC111112489  
LOC111112982  
LOC111116750  
LOC111121978  
LOC111122036  
LOC111123014  
LOC111124965  
LOC111124966  
LOC111125708  
LOC111126890  
LOC111126932  
LOC111128463  
LOC111128610  
LOC111128876  
LOC111128883  
LOC111128917  
LOC111128920  
LOC111129584  
LOC111130548  
LOC111130853  
LOC111132149  
LOC111132218  
LOC111133121  
LOC111134699  
LOC111135334  
LOC111135615  
LOC111136423  
LOC111136478  
LOC111137362  
LOC111137875  
LOC111138326  

**Search through the reference annotation file because this one contains GenbankIDs. There are other things but I'm not sure if those are helpful**

```
grep -w -f LOCSONLY ref_C_virginica-3.0_top_level.gff3 > LOCSREF
```

Manipulating the LOCSREF file to be only the list of GENBANK ID information. First cut to the 9th column, then separate by semi colons, then cut to the third column, then separate by commas, then cut to the second column, then separate from the work genbank because I don't really need it, then only take the second column. After that there are a lot of duplicates so it's sorted and then removed of duplicates. I think if there isn't a genbank ID associated with the transcript it then just says Name:LOCxxx. So those aren't helpful. For some reason using awk to search for X* didn't work to separate out the IDs, so I searched for everything but the Name: ones. This left one number sting in the file that I went in with nano to delete out.  
This is what I'm trying to separate out:
![genbank](/images/genbank.png)

```
cut -f9 LOCSREF | sed -e 's/;/  /g' | cut -f3 | sed -e 's/,/  /g' | cut -f2 | sed -e 's/:/ /g' | cut -f2 > GYGENBANK

sort -u GYGENBANK | uniq > GYGENBANKuniq

awk '!/Name*/' GYGENBANKuniq > GYGENBANKX

wc -l GYGENBANKX
```

This is 175 lines. I don't know if I did this right.

This is actually right, it's just that for genes with multiple transcript varients (alternate splicing) each one of those gets a different genbank # for like coding sequence (CDS), exon, mRNA etc.   
Also, it seems like the XMxxx numbers are the numbers for the exon, and the XPxxx numbers are the numbers for the CDS, but the annotations are the same. So I split up GYGENBANKX to be just the XP#s, which there are 84. And in those 84 it contains the multiple transcript variants of the actual 47 genes with annotations.

```
grep XP_* GYGENBANKX > GYGENBANKXP
```
I don't know why XP* is not enough information for it to grep properly but it's not.

This is 84. Making cxv file of these and putting them into Genbank I guess

XP_022286306.1
XP_022287213.1
XP_022287333.1
XP_022288148.1
XP_022291277.1
XP_022291285.1
XP_022291347.1
XP_022294576.1
XP_022294577.1
XP_022294578.1
XP_022295105.1
XP_022295106.1
XP_022295107.1
XP_022296222.1
XP_022297489.1
XP_022297490.1
XP_022297491.1
XP_022298750.1
XP_022300688.1
XP_022302378.1
XP_022302456.1
XP_022302860.1
XP_022305370.1
XP_022305379.1
XP_022305388.1
XP_022305711.1
XP_022306597.1
XP_022311450.1
XP_022319187.1
XP_022319188.1
XP_022319269.1
XP_022320799.1
XP_022324042.1
XP_022324043.1
XP_022324044.1
XP_022325476.1
XP_022327518.1
XP_022327519.1
XP_022327520.1
XP_022329783.1
XP_022330047.1
XP_022330505.1
XP_022330506.1
XP_022330507.1
XP_022330521.1
XP_022330522.1
XP_022330523.1
XP_022330524.1
XP_022330525.1
XP_022330526.1
XP_022330527.1
XP_022330528.1
XP_022330529.1
XP_022330530.1
XP_022330531.1
XP_022330533.1
XP_022330534.1
XP_022330535.1
XP_022330536.1
XP_022330537.1
XP_022330539.1
XP_022330540.1
XP_022330583.1
XP_022330584.1
XP_022330586.1
XP_022330588.1
XP_022331747.1
XP_022335596.1
XP_022335597.1
XP_022335715.1
XP_022336920.1
XP_022336921.1
XP_022336923.1
XP_022339695.1
XP_022339696.1
XP_022339697.1
XP_022341007.1
XP_022341552.1
XP_022342972.1
XP_022343063.1
XP_022345279.1
XP_022345944.1
XP_022345945.1
XP_022345946.1

 XP are "protein" in Genbank search. coding sequence would be the protein, makes sense.   
 From the page in Genbank for each protein, there are options to BLAST the sequence, look for conserved domains, highlight sequence features, and "find in this sequence" whatever that means


List of GO terms from the GY module
```
mkdir GOGO  
cd GOGO  
ln -s /home/mschedl/Working-CASE-RNA/histat/stringtie/restring/STR.GO.list .
ln -s /home/mschedl/Working-CASE-RNA/histat/stringtie/restring/GYModGeneNames .
grep -w -f GYModGeneNames STR.GO.list > GYmodGO

```

Do the same thing for the lists of sigDEGs

the files that are like "results_CA_F.txt" are all the results from DESeq2. The files that are ordered are subsetted to be only the significant ones padj<0.05. There's no harm in them being ordered. I kinda want to decrease that value...

Need to search only the first line of the txt DEG files. I don't know if I can do that in one but I can probably pipe it to be selecting the first column in the txt file, then using that to search in the second command




```
ln -s /home/mschedl/Working-CASE-RNA/histat/stringtie/restring/ordered_sig_CASE.txt .
ln -s /home/mschedl/Working-CASE-RNA/histat/stringtie/restring/ordered_sig_CA.txt .
ln -s /home/mschedl/Working-CASE-RNA/histat/stringtie/restring/ordered_sig_SE.txt .
cut -f1 ordered_sig_CASE.txt > CASElist
grep -w -f CASElist STR.GO.list > CASE_DEG_GO
cut -f1 ordered_sig_CA.txt > CAlist
grep -w -f CAlist STR.GO.list > CA_DEG_GO
cut -f1 ordered_sig_SE.txt > SElist
grep -w -f SElist STR.GO.list > SE_DEG_GO

```

Hopefully these will go into like TopGO?








# Ok I am going to do something different. Seeing what overlap with the differentially expressed genes from DESeq2

_R code_

This is the overlaping list
```
x
1       MSTRG.13690
2       MSTRG.20077
3       MSTRG.2530
4       MSTRG.27893
5       MSTRG.16274
6       MSTRG.11988
7       MSTRG.5955
8       MSTRG.25599
9       MSTRG.13878
10      MSTRG.25289
11      MSTRG.2984
12      MSTRG.10053
13      MSTRG.8152
14      MSTRG.14902
```


What are the genbank IDs for these? Hopefully they all have them



```
cut -f2 GYSECASEoverlap.txt > GYCASESEoverlap
grep -w -f GYCASESEoverlap C_Vir_ST_merged_NO_A_test.gtf > GYoverlapref

cut -f9 GYoverlapref| sed -e 's/;/  /g' | cut -f4 | sed -e 's/ /  /g' | cut -f3 | sort -u | uniq | sed -e 's/"//g' | sed '/^$/d' | awk '/LOC*/' > overlapLOCS

wc -l overlapLOCS
less overlapLOCS
```
LOC111104759  
LOC111106005  
LOC111112982  
LOC111122036  
LOC111124965  
LOC111124966  
LOC111125708  
LOC111126890  
LOC111128610  
LOC111128920  
LOC111130548  
LOC111130853  
LOC111134699  
LOC111135334  

```
grep -w -f overlapLOCS ref_C_virginica-3.0_top_level.gff3 > overlapLOCSREF
cut -f9 overlapLOCSREF | sed -e 's/;/  /g' | cut -f3 | sed -e 's/,/  /g' | cut -f2 | sed -e 's/:/ /g' | cut -f2 > overlapGB

sort -u overlapGB | uniq > overlapGBuniq

awk '!/Name*/' overlapGBuniq > overlapGBX

wc -l overlapGBX
less overlapGBX
```
43  
Dbxref=GeneID 111130548  
Genbank XM_022438868.1  
Genbank XM_022438869.1  
Genbank XM_022438870.1  
Genbank XM_022440514.1  
Genbank XM_022450889.1  
Genbank XM_022463561.1  
Genbank XM_022468334.1  
Genbank XM_022468335.1  
Genbank XM_022468336.1  
Genbank XM_022469768.1  
Genbank XM_022471810.1  
Genbank XM_022471811.1  
Genbank XM_022471812.1  
Genbank XM_022474339.1  
Genbank XM_022474880.1  
Genbank XM_022483987.1  
Genbank XM_022483988.1  
Genbank XM_022483989.1  
Genbank XM_022485299.1  
Genbank XP_022294576.1  
Genbank XP_022294577.1  
Genbank XP_022294578.1  
Genbank XP_022296222.1  
Genbank XP_022306597.1  
Genbank XP_022319269.1  
Genbank XP_022324042.1  
Genbank XP_022324043.1  
Genbank XP_022324044.1  
Genbank XP_022325476.1  
Genbank XP_022327518.1  
Genbank XP_022327519.1  
Genbank XP_022327520.1  
Genbank XP_022330047.1  
Genbank XP_022330588.1  
Genbank XP_022339695.1  
Genbank XP_022339696.1  
Genbank XP_022339697.1  
Genbank XP_022341007.1  
Genbank XR_002638988.1  
Genbank XR_002638989.1  
Genbank XR_002638990.1  
Genbank XR_002638991.1  

Again only the XP ones I am going to search in Genbank. These should all be ones from the other list.
