CASE-RNA-WGCNA
================
Maggie Schedl
1/30/2020

Setup
=====

``` r
getwd()
```

    ## [1] "/home/mschedl/Working-CASE-RNA/histat/stringtie/restring"

``` r
library(WGCNA)
```

    ## Loading required package: dynamicTreeCut

    ## Loading required package: fastcluster

    ##
    ## Attaching package: 'fastcluster'

    ## The following object is masked from 'package:stats':
    ##
    ##     hclust

    ##

    ##
    ## Attaching package: 'WGCNA'

    ## The following object is masked from 'package:stats':
    ##
    ##     cor

``` r
library(genefilter)
library(dplyr)
```

    ##
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ##
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ##
    ##     intersect, setdiff, setequal, union

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ##
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ##
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ##
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ##
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ##
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colMeans,
    ##     colnames, colSums, dirname, do.call, duplicated, eval, evalq,
    ##     Filter, Find, get, grep, grepl, intersect, is.unsorted, lapply,
    ##     lengths, Map, mapply, match, mget, order, paste, pmax, pmax.int,
    ##     pmin, pmin.int, Position, rank, rbind, Reduce, rowMeans, rownames,
    ##     rowSums, sapply, setdiff, sort, table, tapply, union, unique,
    ##     unsplit, which, which.max, which.min

    ##
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ##
    ##     first, rename

    ## The following object is masked from 'package:base':
    ##
    ##     expand.grid

    ## Loading required package: IRanges

    ##
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ##
    ##     collapse, desc, slice

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ##
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ##
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ##
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ##
    ##     count

    ## The following objects are masked from 'package:genefilter':
    ##
    ##     rowSds, rowVars

    ## Loading required package: BiocParallel

    ##
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ##
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following objects are masked from 'package:base':
    ##
    ##     aperm, apply

``` r
library(tidyr)
```

    ##
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:S4Vectors':
    ##
    ##     expand

``` r
library(edgeR)
```

    ## Loading required package: limma

    ##
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:DESeq2':
    ##
    ##     plotMA

    ## The following object is masked from 'package:BiocGenerics':
    ##
    ##     plotMA

Start with data input and cleaning <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf>

``` r
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#Read in the female liver data set
GeneCountData <- as.data.frame(read.csv("gene_count_matrix.csv", row.names="gene_id"))
# Take a quick look at what is in the data set:
dim(GeneCountData) # what are the dimentions?
```

    ## [1] 45481    14

``` r
names(GeneCountData) # what are the names
```

    ##  [1] "CASE_J03" "CASE_J09" "CASE_J12" "CASE_J13" "CA_J06"   "CA_J08"  
    ##  [7] "CA_J11"   "CA_J18"   "CON_J02"  "CON_J05"  "CON_J10"  "SE_J01"  
    ## [13] "SE_J04"   "SE_J07"

We first check for genes and samples with too many missing values:

``` r
gsg = goodSamplesGenes(GeneCountData, verbose = 3)
```

    ##  Flagging genes and samples with too many missing values...
    ##   ..step 1

``` r
gsg$allOK
```

    ## [1] TRUE

Pre-filtering
=============

Preliminary filtering of the counts with P over A, Meaning that I wanted to remove all rows that have less than 21.4% of the samples with less than 5 counts. 21.4% was chosen because that is 3/14, or the minimum ratio of samples per one treatment. These are the raw counts

``` r
###filtering values for PoverA
#set filter values for PoverA, P percent of the samples have counts over A
filt <- filterfun(pOverA(0.214,5))

#create filter for the counts data
tfil <- genefilter(GeneCountData, filt)

#identify transcripts to keep by count filter
keep <- GeneCountData[tfil,]

#identify transcript list
gn.keep <- rownames(keep)

#data filtered in PoverA, P percent of the samples have counts over A
GeneCountData_Filt <- as.data.frame(GeneCountData[which(rownames(GeneCountData) %in% gn.keep),])
head(GeneCountData_Filt,10)
```

    ##             CASE_J03 CASE_J09 CASE_J12 CASE_J13 CA_J06 CA_J08 CA_J11 CA_J18
    ## MSTRG.10383      198      306      198      235    190    133    180    176
    ## MSTRG.28362       18        7       18       18      3      4     11     30
    ## MSTRG.10380        5       14        9       19      6     23     12      7
    ## MSTRG.32256        4        4        4        8     15      0     10      7
    ## MSTRG.18313       92      169      165      206    120    129    134    236
    ## MSTRG.10381      143      146       90      124     65     42     60     91
    ## MSTRG.19848       42       74       21       55     51     32     29     59
    ## MSTRG.19849       50       38       37       75     21     33     43     52
    ## MSTRG.19846       26       94       79       78     44     59     31     59
    ## MSTRG.19847       12       27        4       11     26     11      3      2
    ##             CON_J02 CON_J05 CON_J10 SE_J01 SE_J04 SE_J07
    ## MSTRG.10383     131     139     276    179    104    259
    ## MSTRG.28362      16       7       8     13      5      2
    ## MSTRG.10380      26       5       3     11      2     47
    ## MSTRG.32256       5       5       0      7     14      5
    ## MSTRG.18313     226     340     374    159     86    313
    ## MSTRG.10381      69      76     109     45     63     46
    ## MSTRG.19848      49      24      24     33     18     45
    ## MSTRG.19849      49      30      60     37     50     39
    ## MSTRG.19846      80      22      12     46     30     79
    ## MSTRG.19847       0       3       0      4      5      7

The question of how to normalize this dataset came up. The creators of WGCNA insist that "properly normalized" data is required. It turns out transforming your data (ex vst or log) is a type of normalization method. Going forward in this script are two different normalization methods. One use the TMM method from EdgeR and then log2 transforms it from there. As far as I can tell, the reason to log2 (or any base log) transform any data is to make it more symetrical. It is likely true that my data has a skew in it.

The other method is to use the varience stablizing transformation on just the raw counts. "Many common statistical methods for exploratory analysis of multidimensional data, especially methods for clustering and ordination (e.g., principal-component analysis and the like), work best for (at least approximately) homoskedastic data; this means that the variance of an observable quantity (i.e., here, the expression strength of a gene) does not depend on the mean. In RNA-Seq data, however, variance grows with the mean." <https://www.biostars.org/p/109104/> <http://www.nxn.se/valent/2017/10/15/variance-stabilizing-scrna-seq-counts> The count data is in no way normally distributed "The process of counting implies that variation will propagate as the number of events increase. The effect of this is that there will be an inherent relation between mean (expected value) and variance of counts" From what I can gather the vst transformation removes the relationship between the mean and the varience.

TMM Normalization
=================

<https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html>

I think I should use EdgeR’s trimmed mean of M values (TMM) which is for "gene count comparisons between and within samples and for DE analysis" I thing WGCNA does both within and between sample comparisons

``` r
TMM <- calcNormFactors(GeneCountData_Filt, lib.size = NULL, method = "TMM")

print(TMM) #these are the normalization factors calculated by the TMM method
```

    ##  [1] 0.9784394 0.9638976 1.0826775 0.9558671 1.1030259 1.0265708 1.0260582
    ##  [8] 0.9599571 0.9395334 1.0321207 0.9280542 1.0602186 0.9412628 1.0228635

Now to appy these to the filtered raw counts dataset. I wasn't exactualy sure how to apply this

<https://support.bioconductor.org/p/92781/>

``` r
lib <- colSums(GeneCountData_Filt) #get the sum of all the counts per each sample. This is the library size
effectivelibsize <- TMM*lib # the effective library size is the normalization factor multiplied by the library size

TMMnormalized <- GeneCountData_Filt/effectivelibsize # then you divide the counts by the effective library size.



head(TMMnormalized)
```

    ##                 CASE_J03     CASE_J09     CASE_J12     CASE_J13       CA_J06
    ## MSTRG.10383 3.514085e-06 5.715366e-06 4.225733e-06 4.608969e-06 3.176090e-06
    ## MSTRG.28362 2.635926e-07 9.867294e-08 3.955495e-07 2.576519e-07 5.187845e-08
    ## MSTRG.10380 9.338833e-08 2.987892e-07 1.765137e-07 3.176090e-07 9.251642e-08
    ## MSTRG.32256 5.638453e-08 8.789988e-08 5.725599e-08 1.383425e-07 2.514815e-07
    ## MSTRG.18313 1.963472e-06 3.314535e-06 2.758183e-06 3.176397e-06 2.281876e-06
    ## MSTRG.10381 3.142421e-06 2.089844e-06 1.556354e-06 2.078914e-06 1.018430e-06
    ##                   CA_J08       CA_J11       CA_J18      CON_J02      CON_J05
    ## MSTRG.10383 2.050781e-06 3.422814e-06 3.123631e-06 2.446774e-06 2.966550e-06
    ## MSTRG.28362 6.706175e-08 1.723497e-07 4.393209e-07 2.255381e-07 1.538248e-07
    ## MSTRG.10380 4.373596e-07 2.129748e-07 1.307437e-07 5.548943e-07 9.806317e-08
    ## MSTRG.32256 0.000000e+00 1.464403e-07 9.867294e-08 1.098749e-07 7.156998e-08
    ## MSTRG.18313 2.289479e-06 2.502807e-06 5.036732e-06 4.432455e-06 5.683529e-06
    ## MSTRG.10381 6.150493e-07 8.457680e-07 1.999722e-06 9.876658e-07 1.314254e-06
    ##                  CON_J10       SE_J01       SE_J04       SE_J07
    ## MSTRG.10383 5.413087e-06 2.992211e-06 1.603618e-06 4.925050e-06
    ## MSTRG.28362 1.145120e-07 2.248066e-07 8.382718e-08 3.133631e-08
    ## MSTRG.10380 5.014879e-08 1.696134e-07 3.803127e-08 8.341514e-07
    ## MSTRG.32256 0.000000e+00 1.173581e-07 2.193541e-07 7.322015e-08
    ## MSTRG.18313 5.766857e-06 3.023486e-06 1.526320e-06 5.846109e-06
    ## MSTRG.10381 1.827433e-06 7.050669e-07 9.225739e-07 6.484221e-07

``` r
#Now to get the 5000 most variable genes from this

TMMnormalized <- as.matrix(TMMnormalized) # make into a matrix so the vars fucntion will take it as an input
vars <- rowVars(TMMnormalized) # gets the varience for each row (each row is a gene) and makes it into a vector list
vars <- as.data.frame(vars) # lists aren't helpful so make it into a data frame, one column basically

TMMnormalized <- as.data.frame(TMMnormalized) # turn the counts data back into a dataframe
GeneCountData_Filt$nNames <- rownames(GeneCountData_Filt) # create a new column in this dataframe that makes the row names an actual column
rownames(vars) <- GeneCountData_Filt$nNames # use that column as the template for the rownames of the vars dataframe

TMMnormalizedvar <- cbind(TMMnormalized, vars) #add the row varience to the counts dataframe , binds by having the same rownames
TMMnormalizedorderedvar <- TMMnormalizedvar[order(-vars),] # order this dataframe by the varience column and in decreasing order, so highest varience should be first
orderedTMMM <- head(TMMnormalizedorderedvar, 5000) # then only take the top 5000 amount for a scaled down dataset
head(orderedTMMM)
```

    ##                CASE_J03    CASE_J09    CASE_J12    CASE_J13      CA_J06
    ## MSTRG.1     0.180423851 0.126283647 0.051492168 0.162353437 0.058956911
    ## MSTRG.18894 0.024214215 0.037941031 0.022093576 0.052552413 0.001060638
    ## MSTRG.18893 0.027310166 0.029972794 0.013865019 0.033347113 0.001293406
    ## MSTRG.31291 0.017400327 0.022328608 0.020982857 0.024128248 0.015023490
    ## MSTRG.19512 0.021555557 0.016410511 0.016280185 0.019270444 0.012830776
    ## MSTRG.3     0.003741206 0.005782597 0.003992706 0.005417765 0.001233686
    ##                  CA_J08      CA_J11     CA_J18     CON_J02     CON_J05
    ## MSTRG.1     0.053673477 0.048803563 0.14191410 0.133536578 0.064135760
    ## MSTRG.18894 0.016010910 0.040852130 0.05686509 0.024460346 0.014284487
    ## MSTRG.18893 0.013296764 0.019155553 0.05714647 0.020023951 0.010783343
    ## MSTRG.31291 0.012163119 0.018132702 0.01914738 0.026905290 0.018432029
    ## MSTRG.19512 0.010863059 0.012707058 0.02274095 0.017663372 0.016197336
    ## MSTRG.3     0.002223299 0.006963495 0.00452162 0.005334834 0.002348519
    ##                 CON_J10      SE_J01       SE_J04      SE_J07         vars
    ## MSTRG.1     0.070504043 0.055972394 0.0881147290 0.116344229 2.082694e-03
    ## MSTRG.18894 0.037739693 0.032460149 0.0279067359 0.003797802 2.707578e-04
    ## MSTRG.18893 0.024217592 0.039839506 0.0253586639 0.002789846 2.200032e-04
    ## MSTRG.31291 0.022535839 0.014247955 0.0124403595 0.026465201 2.327045e-05
    ## MSTRG.19512 0.019379420 0.013580983 0.0120393702 0.017720954 1.317600e-05
    ## MSTRG.3     0.007256846 0.005447111 0.0006913449 0.009642485 6.141052e-06

``` r
#now need to remove that last column and transform
orderedTMMM$vars <- NULL #remove the last column because it's not needed anymore
torderedTMM <- as.data.frame(t(as.matrix(orderedTMMM))) # transpose it to be samples as rows for next steps


GeneCountData_Filt$nNames <- NULL

log2orderedTMMM <- log2(orderedTMMM + 1) #log 2 transformation on the TMM normalized data ordered and subsetted to 5000 most varience
head(log2orderedTMMM)
```

    ##                CASE_J03    CASE_J09    CASE_J12    CASE_J13      CA_J06
    ## MSTRG.1     0.239304977 0.171570206 0.072438104 0.217048816 0.082643888
    ## MSTRG.18894 0.034517488 0.053724482 0.031527286 0.073892075 0.001529366
    ## MSTRG.18893 0.038871827 0.042606230 0.019865592 0.047324953 0.001864784
    ## MSTRG.31291 0.024887463 0.031858998 0.029958642 0.034396390 0.021513115
    ## MSTRG.19512 0.030767666 0.023483199 0.023298204 0.027536894 0.018393149
    ## MSTRG.3     0.005387349 0.008318496 0.005748787 0.007795086 0.001778736
    ##                  CA_J08     CA_J11      CA_J18     CON_J02     CON_J05
    ## MSTRG.1     0.075427859 0.06874449 0.191454127 0.180830946 0.089682218
    ## MSTRG.18894 0.022915894 0.05776513 0.079791221 0.034864142 0.020462357
    ## MSTRG.18893 0.019056757 0.02737427 0.080175276 0.028603028 0.015473795
    ## MSTRG.31291 0.017441811 0.02592561 0.027362693 0.038303130 0.026349697
    ## MSTRG.19512 0.015587570 0.01821691 0.032440775 0.025260418 0.023180587
    ## MSTRG.3     0.003203981 0.01001138 0.006508615 0.007676082 0.003384225
    ##                CON_J10      SE_J01       SE_J04      SE_J07
    ## MSTRG.1     0.09829024 0.078572120 0.1218306799 0.158781956
    ## MSTRG.18894 0.05344460 0.046086097 0.0397093717 0.005468693
    ## MSTRG.18893 0.03452224 0.056360873 0.0361286435 0.004019292
    ## MSTRG.31291 0.03215141 0.020410394 0.0178369247 0.037684718
    ## MSTRG.19512 0.02769113 0.019461362 0.0172654146 0.025342048
    ## MSTRG.3     0.01043161 0.007837194 0.0009970552 0.013844525

``` r
tlog2orderedTMMM <- as.data.frame(t(as.matrix(log2orderedTMMM)))
```

Varience Stablization Transformation
====================================

Preparing to varience stablizing transform the filtered raw counts dataset. This function takes a DESeq2 dataset as the input so first have to make one of those.

``` r
# beginning part copy and pasted from DESeq2 rmarkdown
# Make sure things are in the correct format

CASE_treatment <- read.csv("treatment_data.csv", header=TRUE, sep=",") # need treatment info for formatting purposes
rownames(CASE_treatment) <- CASE_treatment$sample
colnames(GeneCountData_Filt) <- CASE_treatment$sample


#The row and column names for the two data frames need to be exactly the same for the rest of the analysis, so it is good to check
all(rownames(CASE_treatment) %in% colnames(GeneCountData_Filt))  #Should return TRUE
```

    ## [1] TRUE

``` r
all(rownames(CASE_treatment) == colnames(GeneCountData_Filt))    # should return TRUE
```

    ## [1] TRUE

``` r
# Create the matrix to be able to use as an input for normalizing, using just the filtered raw counts

CASE_deseq_Matrix <- DESeqDataSetFromMatrix(countData = GeneCountData_Filt,
                              colData = CASE_treatment,
                              design = ~ treatment )
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
CASE_deseq_Matrix
```

    ## class: DESeqDataSet
    ## dim: 34582 14
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(34582): MSTRG.10383 MSTRG.28362 ... gene6081 MSTRG.10385
    ## rowData names(0):
    ## colnames(14): CASE_J03 CASE_J09 ... SE_J04 SE_J07
    ## colData names(4): sample treatment library extraction

Now apply the varience stablizing transformation on the matrix, then subset it to the top 5000 most varying genes in this dataset

``` r
vCASE_deseq_Matrix <- vst(CASE_deseq_Matrix, blind = FALSE) # don't use normaized matrix, this transformation normalizes the matrix
head(assay(vCASE_deseq_Matrix), 3)
```

    ##             CASE_J03 CASE_J09 CASE_J12 CASE_J13   CA_J06   CA_J08   CA_J11
    ## MSTRG.10383 8.209232 8.468603 8.264432 8.160522 8.349202 8.026280 8.219437
    ## MSTRG.28362 6.551011 6.193528 6.573460 6.465736 6.090871 6.156474 6.406539
    ## MSTRG.10380 6.158933 6.389667 6.325759 6.485900 6.246055 6.758175 6.436777
    ##               CA_J18  CON_J02  CON_J05  CON_J10   SE_J01   SE_J04   SE_J07
    ## MSTRG.10383 7.897232 7.758752 7.854094 8.412348 8.068422 7.700919 8.349057
    ## MSTRG.28362 6.684721 6.480058 6.235637 6.239351 6.413471 6.181812 5.977497
    ## MSTRG.10380 6.188935 6.684046 6.155329 6.036274 6.358263 6.010130 6.959458

``` r
MostVaryGenes <- head(order(rowVars(assay(vCASE_deseq_Matrix)), decreasing = TRUE), 5000)

# make that information into a matrix
VMostVaryGenes_Mat <- assay(vCASE_deseq_Matrix)[ MostVaryGenes, ] #which in the DESeq2 style matrix are those genes?
# make that matrix a dataframe
VMostVaryGenes_df <- as.data.frame(VMostVaryGenes_Mat)

tVMostVaryGenes_df <- as.data.frame(t(as.matrix(VMostVaryGenes_df))) #transpose to right format
```

Now everything going forward works with both types of transformed data. If it had the vst transformation it's named with a "v" If it had the TMM and log2 transformation it's named with "log"

Exploritory Sample Clustering
=============================

Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outlier samples. This also lets you see where there are potential big splits in the data.

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
sampleTree = hclust(dist(tVMostVaryGenes_df), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
```
![clust](/images/vst-clust.png)

Not really any obvious outliers. There are two smaller clusters of two/three samples. Not sure how much that's a problem.

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
sampleTree = hclust(dist(tlog2orderedTMMM), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
```
![clust](/images/tmm-clust.png)

There is a large split between two groups here. This could be why there is a problem reaching scale free topology below for this method.

Soft Thresholding Power
=======================

The next thing to do is determine the power or soft threshold to apply for the adjacency matrix. This is because we want to make a weighted matrix, where each gene has a weighted connection, not a hard threshold of either connected or not. From my understanding, by raising connection/adjacencies to a power there in lies the weighted-ness of their connection.

I need to pick a power where scale free topolgy of the matrix is met. From some googling, this is similar to a scale free network. "A scale-free network can be constructed by progressively adding nodes to an existing network and introducing links to existing nodes with preferential attachment so that the probability of linking to a given node i is proportional to the number of existing links k\_i that node has" So there is a higher chance of being connected to a node gene if it has more connections. That makes intuitive sense, which is why it has to reach this criterion to be passable for further analysis. <https://mathworld.wolfram.com/Scale-FreeNetwork.html>

"If the scale-free topology fit index fails to reach values above 0.8 for reasonable powers (less than 15 for unsigned or signed hybrid networks, and less than 30 for signed networks) and the mean connectivity remains relatively high (in the hundreds or above), chances are that the data exhibit a strong driver that makes a subset of the samples globally different from the rest. The difference causes high correlation among large groups of genes which invalidates the assumption of the scale-free topology approximation."<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html>

"If the lack of scale-free topology fit turns out to be caused by an interesting biological variable that one does not want to remove (i.e., adjust the data for), the appropriate soft-thresholding power can be chosen based on the number of samples as in the table below. This table has been updated in December 2017 to make the resulting networks conservative." For a signed network with less than 15 samples, 18 is an appropreate power to use.

I am doing a signed network analysis. <https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/>

something to think about <https://peterlangfelder.com/2018/11/25/__trashed/>

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(tVMostVaryGenes_df, powerVector = powers, verbose = 5, networkType = "signed") #make sure to specify network type
```

    ## pickSoftThreshold: will use block size 5000.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 5000 of 5000

    ## Warning: executing %dopar% sequentially: no parallel backend registered

    ##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1   0.7450  2.510          0.915  2770.0    2930.0   3310
    ## 2      2   0.0319  0.226          0.777  1790.0    1890.0   2580
    ## 3      3   0.0494 -0.228          0.613  1270.0    1300.0   2140
    ## 4      4   0.2300 -0.403          0.611   959.0     936.0   1840
    ## 5      5   0.4020 -0.490          0.616   756.0     695.0   1620
    ## 6      6   0.5390 -0.540          0.662   615.0     525.0   1440
    ## 7      7   0.6470 -0.573          0.706   513.0     407.0   1300
    ## 8      8   0.7390 -0.596          0.780   436.0     317.0   1190
    ## 9      9   0.8170 -0.622          0.838   376.0     252.0   1090
    ## 10    10   0.8530 -0.647          0.856   329.0     204.0   1010
    ## 11    12   0.9050 -0.694          0.885   258.0     143.0    884
    ## 12    14   0.8790 -0.743          0.844   209.0     106.0    790
    ## 13    16   0.8630 -0.800          0.832   173.0      79.1    715
    ## 14    18   0.8580 -0.843          0.840   146.0      60.7    653
    ## 15    20   0.8410 -0.891          0.834   125.0      48.2    601
    ## 16    22   0.8260 -0.929          0.826   108.0      38.1    556
    ## 17    24   0.8180 -0.964          0.827    94.2      30.2    516
    ## 18    26   0.8200 -0.996          0.839    82.9      24.3    481
    ## 19    28   0.8180 -1.020          0.848    73.4      19.5    450
    ## 20    30   0.8190 -1.050          0.855    65.4      15.7    422

``` r
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```
![pwr](/images/vst-pwr.png)

12 has a .9 signed R2, so I think I want to pick this power.

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(tlog2orderedTMMM, powerVector = powers, verbose = 5, networkType = "signed") #make sure specify network type
```

    ## pickSoftThreshold: will use block size 5000.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 5000 of 5000
    ##    Power SFT.R.sq    slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1 0.954000  7.14000          0.971  3230.0    3320.0   3770
    ## 2      2 0.899000  3.43000          0.946  2240.0    2330.0   2940
    ## 3      3 0.769000  2.02000          0.889  1620.0    1690.0   2400
    ## 4      4 0.545000  1.22000          0.829  1220.0    1260.0   2010
    ## 5      5 0.264000  0.69700          0.775   948.0     966.0   1720
    ## 6      6 0.056000  0.28800          0.726   753.0     754.0   1490
    ## 7      7 0.000023 -0.00561          0.727   609.0     596.0   1300
    ## 8      8 0.043900 -0.24100          0.745   501.0     479.0   1140
    ## 9      9 0.152000 -0.49000          0.735   417.0     390.0   1020
    ## 10    10 0.235000 -0.60800          0.773   352.0     320.0    906
    ## 11    12 0.378000 -0.78600          0.833   258.0     222.0    733
    ## 12    14 0.483000 -0.94700          0.850   195.0     159.0    602
    ## 13    16 0.567000 -1.07000          0.871   152.0     118.0    501
    ## 14    18 0.629000 -1.14000          0.892   121.0      90.0    424
    ## 15    20 0.659000 -1.22000          0.897    97.5      70.9    365
    ## 16    22 0.689000 -1.28000          0.912    80.1      56.6    317
    ## 17    24 0.713000 -1.33000          0.920    66.6      45.7    278
    ## 18    26 0.732000 -1.37000          0.924    56.1      37.5    245
    ## 19    28 0.744000 -1.41000          0.925    47.7      31.2    218
    ## 20    30 0.749000 -1.44000          0.926    40.9      25.9    194

``` r
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```
![pwr](/images/tmm-pwr.png)

Does not reach scale free topology. I will have to chose 18 for this one.

Now create the network adjacency matrix, using the power picked and specifying the sign. The adjacency matrix is for all samples and is between each gene.

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
vsoftPower = 12;
vadjacency = adjacency(tVMostVaryGenes_df, power = vsoftPower, type= "signed");
```

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
logsoftPower = 18;
logadjacency = adjacency(tlog2orderedTMMM, power = logsoftPower, type= "signed");
```

The adjacencies are just correlations of one gene with another, until it does all of the comparissons. But to get a larger network view interconnectedness and similarities of one connection with another need to be taken into account. I need a topological overlap matrix. "More specifically, genes are said to have high topological overlap if they are connected to roughly the same group of genes in the network (i.e. they share the same neighborhood). To calculate the topological overlap for a pair of genes, their connections with all other genes in the network are compared"

To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity:

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
# Turn adjacency into topological overlap
vTOM = TOMsimilarity(vadjacency) # don't specify tomtype so it will use the adjacency ??
```

    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.

``` r
vdissTOM = 1-vTOM # this is the dissimilarity, which is 1 minus the similarity
```

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
# Turn adjacency into topological overlap
logTOM = TOMsimilarity(logadjacency);
```

    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.

``` r
logdissTOM = 1-logTOM
```

<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/WORKSHOP/2014/Langfelder-NetworkDay-clustering.pdf>

We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes. Note that we use the function hclust that provides a much faster hierarchical clustering routine than the standard hclust function.

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
# Call the hierarchical clustering function
vgeneTree = hclust(as.dist(vdissTOM), method = "average");

# Average linkage: average the dissimilarities between all objects
#single (minimum dissimilarity) and complete (maximum dissimilarity) are other options, but seem too harsh

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(vgeneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04);
```
![tree](/images/vst-tree.png)

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
# Call the hierarchical clustering function
loggeneTree = hclust(as.dist(logdissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(loggeneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04);
```
![tree](/images/tmm-tree.png)

Cutting the dendrogram into modules
===================================

"The clustering dendrogram plotted by the last command is shown in Figure 2. In the clustering tree (dendrogram), each leaf, that is a short vertical line, corresponds to a gene. Branches of the dendrogram group together densely interconnected, highly co-expressed genes."

"Module identification amounts to the identification of individual branches (”cutting the branches off the dendrogram”)"

Want to dynamically cut the tree height because no single cut height captures the branches "Branches are followed bottom to top. When two branches merge, they are evaluated using shape criteria such as minimum number of objects (genes), their core scatter and the gap between the branches" So the most similar are what are clustered first and then followed up the branches The core areas are the bottom ends of branches where the most similarly expressed genes are. The gap is the height difference between the cut height and the bottom of the branch. The core scatter is the amount of horizontal space between the branch where the tree is cut and the core section.

"DeepSplit controls how finely the branches should be split Higher values give more smaller modules, low values (0) give fewer larger modules" **I didn't change this, maybe I should have considered it. It's 2, which is low I guess**

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30; #must have 30 genes in it to be a module
# Module identification using dynamic tree cut:
logdynamicMods = cutreeDynamic(dendro = loggeneTree, distM = logdissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
```

    ##  ..cutHeight not given, setting it to 0.996  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.

``` r
table(logdynamicMods)
```

    ## logdynamicMods
    ##   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
    ##   1 511 191 173 156 146 145 139 128 126 125 118 117 114 103 100  96  95  93  91
    ##  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39
    ##  91  91  87  85  84  83  81  80  78  77  75  74  71  71  70  68  68  60  60  59
    ##  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55
    ##  58  58  52  52  50  48  48  48  47  41  39  37  37  36  35  33

55 modules.

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
# Convert numeric lables into colors
logdynamicColors = labels2colors(logdynamicMods)
table(logdynamicColors)
```

    ## logdynamicColors
    ##         bisque4           black            blue           brown          brown4
    ##              48             139             191             173              48
    ##            cyan       darkgreen        darkgrey     darkmagenta  darkolivegreen
    ##             103              87              84              70              71
    ##      darkorange     darkorange2         darkred   darkslateblue   darkturquoise
    ##              81              48              91              47              85
    ##     floralwhite           green     greenyellow            grey          grey60
    ##              50             146             118               1              95
    ##           ivory       lightcyan      lightcyan1      lightgreen lightsteelblue1
    ##              52              96              52              93              58
    ##     lightyellow         magenta          maroon   mediumpurple3    midnightblue
    ##              91             126              33              58             100
    ##    navajowhite2          orange      orangered4   paleturquoise  palevioletred3
    ##              35              83              59              74              36
    ##            pink           plum1           plum2          purple             red
    ##             128              60              41             125             145
    ##       royalblue     saddlebrown          salmon         salmon4         sienna3
    ##              91              77             114              37              68
    ##         skyblue        skyblue3       steelblue             tan        thistle1
    ##              78              60              75             117              37
    ##        thistle2       turquoise          violet           white          yellow
    ##              39             511              71              80             156
    ##     yellowgreen
    ##              68

``` r
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(loggeneTree, logdynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
```

![mods](/images/tmm-mods.png)

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
vdynamicMods = cutreeDynamic(dendro = vgeneTree, distM = vdissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
```

    ##  ..cutHeight not given, setting it to 0.99  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.

``` r
table(vdynamicMods)
```

    ## vdynamicMods
    ##    1    2    3    4    5    6    7    8    9   10   11   12   13
    ## 1907 1146  495  375  253  210  171  131   87   71   57   49   48

13 modules

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
# Convert numeric lables into colors
vdynamicColors = labels2colors(vdynamicMods)
table(vdynamicColors)
```

    ## vdynamicColors
    ##       black        blue       brown       green greenyellow     magenta
    ##         171        1146         495         253          57          87
    ##        pink      purple         red      salmon         tan   turquoise
    ##         131          71         210          48          49        1907
    ##      yellow
    ##         375

``` r
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(vgeneTree, vdynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
```

![mods](/images/vst-mods.png)

Correlating modules to each other and creating a threshold where they should be merged
======================================================================================

"The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation:"

Color "grey" is reserved for unassigned genes. grey is excluded from being a module with excludeGREY I don't think I have to include the soft power here "The power used in soft-thresholding the adjacency matrix. Only used when the hubgene approximation is necessary because the principal component calculation failed. It must be non-negative. The default value should only be changed if there is a clear indication that it leads to incorrect results" I don't think any calculations ever fail. trapErrors = "Controls handling of errors from that may arise when there are too many NA entries in expression data. If TRUE, errors from calling these functions will be trapped without abnormal exit. If FALSE, errors will cause the function to stop. Note, however, that subHubs takes precedence in the sense that if subHubs==TRUE and trapErrors==FALSE, an error will be issued only if both the principal component and the hubgene calculations have failed." don't think there are ever any errors with this but just in case

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
# Calculate eigengenes
logMEList = moduleEigengenes(tlog2orderedTMMM, colors = logdynamicColors, trapErrors = FALSE, excludeGrey = TRUE)
logMEs = logMEList$eigengenes
# Calculate dissimilarity of module eigengenes
logMEDiss = 1-cor(logMEs);
# Cluster module eigengenes
logMETree = hclust(as.dist(logMEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(logMETree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

logMEDissThres = 0.20 # mege modules if they are 80% correlated
# Plot the cut line into the dendrogram
abline(h=logMEDissThres, col = "red")
# Call an automatic merging function
logmerge = mergeCloseModules(tlog2orderedTMMM, logdynamicColors, cutHeight = logMEDissThres, verbose = 3)
```

    ##  mergeCloseModules: Merging modules whose distance is less than 0.2
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 56 module eigengenes in given set.
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 30 module eigengenes in given set.
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 25 module eigengenes in given set.
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 23 module eigengenes in given set.
    ##    Calculating new MEs...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 23 module eigengenes in given set.

``` r
# The merged module colors
logmergedColors = logmerge$colors;
# Eigengenes of the new merged modules:
logmergedMEs = logmerge$newMEs;
```

![modtree](/images/tmm-modtree.png)

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
# Calculate eigengenes
vMEList = moduleEigengenes(tVMostVaryGenes_df, colors = vdynamicColors, trapErrors = FALSE, excludeGrey = TRUE)
vMEs = vMEList$eigengenes
# Calculate dissimilarity of module eigengenes
vMEDiss = 1-cor(vMEs);
# Cluster module eigengenes
vMETree = hclust(as.dist(vMEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(vMETree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

vMEDissThres = 0.2 # this means to merge modlues if they are 99.7% correlated..... something is wrong here
# Plot the cut line into the dendrogram
abline(h=vMEDissThres, col = "red")
# Call an automatic merging function
vmerge = mergeCloseModules(tVMostVaryGenes_df, vdynamicColors, cutHeight = vMEDissThres, verbose = 3)
```

    ##  mergeCloseModules: Merging modules whose distance is less than 0.2
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 13 module eigengenes in given set.
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 12 module eigengenes in given set.
    ##    Calculating new MEs...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 12 module eigengenes in given set.

``` r
# The merged module colors
vmergedColors = vmerge$colors;
# Eigengenes of the new merged modules:
vmergedMEs = vmerge$newMEs;
```
![modtree](/images/vst-modtree.png)

Merging modules
===============

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(loggeneTree, cbind(logdynamicColors, logmergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
```
![mergemod](/images/tmm-mergemod.png)

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(vgeneTree, cbind(vdynamicColors, vmergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
```
![mergemod](/images/vst-mergemod.png)

Relating modules to external information and identifying important genes
========================================================================

<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf>

set up

``` r
#load in bianary treatment data
Traitdat <- as.data.frame(read.csv("Trait_Data.csv"))

Traitdat2 <- Traitdat[,-1] #make column 1 into row names
rownames(Traitdat2) <- Traitdat[,1]
print(Traitdat2)
```

    ##          pCO2 Effluent
    ## CASE_J03    1        1
    ## CASE_J09    1        1
    ## CASE_J12    1        1
    ## CASE_J13    1        1
    ## CA_J06      1        0
    ## CA_J08      1        0
    ## CA_J11      1        0
    ## CA_J18      1        0
    ## CON_J02     0        0
    ## CON_J05     0        0
    ## CON_J10     0        0
    ## SE_J01      0        1
    ## SE_J04      0        1
    ## SE_J07      0        1

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
#set up for using the module colors in further code
# Rename to moduleColors
logmoduleColors = logmergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
logmoduleLabels = match(logmoduleColors, colorOrder)-1;
logMEs = logmergedMEs;
```

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
#set up for using the module colors in further code
# Rename to moduleColors
vmoduleColors = vmergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
vmoduleLabels = match(vmoduleColors, colorOrder)-1;
vMEs = vmergedMEs;
```

"In this analysis we would like to identify modules that are significantly associated with the measured clinical traits. Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external traits and look for the most significant associations:"

"**Robust correlation.** The default correlation method in all functions in WGCNA is standard Pearson correlation. In general, unless there is good reason to believe that there are no outlier measurements, we recommend (and use ourselves) the biweight mid-correlation as a robust alternative. This is implemented in WGCNA function bicor. Many WGCNA functions take the argument corFnc that allows one to specify an alternative correlation function to the standard cor and bicor is one option. Additional arguments to the correlation function can be specified using the argument corOptions (depending on function, this argument may require one of two alternate forms, please see the help for each function for details). In certain functions, notably the of the blockwise family, correlation function cannot be specified directly as a function; rather, one must use the argument corType to specify either Pearson or biweight mid-correlation" "**Dealing with binary data.** When relating high-throughput data x to binary variable y such as sample traits, one can use argument robustY = FALSE to turn off the robust treatment for the y argment of bicor. This results in a hybrid robust-Pearson correlation as described in Langfelder and Horvath (2011). The hybrid correlation can also be used when one of the inputs is numeric but known to not have any outliers." <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html> <https://pubmed.ncbi.nlm.nih.gov/23050260-fast-r-functions-for-robust-correlations-and-hierarchical-clustering/>

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
# Define numbers of genes and samples
lognGenes = ncol(tlog2orderedTMMM);
lognSamples = nrow(tlog2orderedTMMM);
# Recalculate MEs with color labels
logMEs0 = moduleEigengenes(tlog2orderedTMMM, logmoduleColors)$eigengenes
logMEs = orderMEs(logMEs0)
logmoduleTraitCor = bicor(logMEs, Traitdat2, robustX = TRUE, robustY = FALSE, use = "all.obs"); #try using hybrid pearson robust correlation
logmoduleTraitPvalue = corPvalueStudent(logmoduleTraitCor, lognSamples)
```

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
# Define numbers of genes and samples
vnGenes = ncol(tVMostVaryGenes_df);
vnSamples = nrow(tVMostVaryGenes_df);
# Recalculate MEs with color labels
vMEs0 = moduleEigengenes(tVMostVaryGenes_df, vmoduleColors)$eigengenes
vMEs = orderMEs(vMEs0)
vmoduleTraitCor = bicor(vMEs, Traitdat2, robustX = TRUE, robustY = FALSE, use = "all.obs"); #try using hybrid pearson robust correlation
vmoduleTraitPvalue = corPvalueStudent(vmoduleTraitCor, vnSamples)
```

Now visualize the correlations of the treatments with the modules in heatmap style.

``` r
### TMM NORMALIZED THEN LOG2 TRANSFORMED
sizeGrWindow(10,6)
# Will display correlations and their p-values
logtextMatrix = paste(signif(logmoduleTraitCor, 2), "\n(",
signif(logmoduleTraitPvalue, 1), ")", sep = "");
dim(logtextMatrix) = dim(logmoduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = logmoduleTraitCor,
xLabels = names(Traitdat2),
yLabels = names(logMEs),
ySymbols = names(logMEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = logtextMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships TMMLog2"))
```

![cor](/images/tmm-cor.png)

There are no significant correlations with this method

``` r
# VST TRANSFORMED AND MOST VARIABLE FROM THOSE
sizeGrWindow(10,6)
# Will display correlations and their p-values
vtextMatrix = paste(signif(vmoduleTraitCor, 2), "\n(",
signif(vmoduleTraitPvalue, 1), ")", sep = "");
dim(vtextMatrix) = dim(vmoduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = vmoduleTraitCor,
xLabels = names(Traitdat2),
yLabels = names(vMEs),
ySymbols = names(vMEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = vtextMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
```
![cor](/images/vst-cor.png)

Ok, with this type of transformed data, the greenyellow module is highly negatively correlated with effluent treatment. Magenta is negatively correlated with pCO2 treatment (is that a 0.05 or 0.06?). And module purple is positively correlated with effluent treatment.

``` r
print(vmoduleTraitPvalue)
```

    ##                     pCO2     Effluent
    ## MEbrown       0.87094402 5.498015e-01
    ## MEmagenta     0.05080923 2.821362e-01
    ## MEblue        0.54234781 6.797402e-01
    ## MEpurple      0.11065918 3.246610e-02
    ## MEgreenyellow 0.72657064 2.431323e-05
    ## MEtan         0.94773974 7.431802e-01
    ## MEsalmon      0.70722052 7.512333e-01
    ## MEpink        0.62772826 2.330855e-01
    ## MEred         0.93423697 5.978912e-01
    ## MEyellow      0.20908817 8.874884e-01
    ## MEblack       0.63623631 5.567852e-01
    ## MEgreen       0.98723601 2.613448e-01

Not sure why it represents the effulent ones like that

I think I am going to continue on only with the vst transformed data for now. It gives me something to go on with!

Visualization of Networks
=========================

<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-05-Visualization.pdf>

Visualizing the gene network
============================

"One way to visualize a weighted network is to plot its heatmap, Fig. 1. Each row and column of the heatmap correspond to a single gene. The heatmap can depict adjacencies or topological overlaps, with light colors denoting low adjacency (overlap) and darker colors higher adjacency (overlap). In addition, the gene dendrograms and module colors are plotted along the top and left side of the heatmap. The package provides a convenient function to create such network plots; Fig. 1 was created using the following code. This code can be executed only if the network was calculated using a single-block approach (that is, using the 1-step automatic or the step-by-step tutorials). If the networks were calculated using the block-wise approach, the user will need to modify this code to perform the visualization in each block separately. The modification is simple and we leave it as an exercise for the interested reader."

This is like a heatmap representation of the dendrogram. I like that better maybe??

Lets do only the vst transformed data.

``` r
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
 # dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# I did do that!
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = vdissTOM^7; # how do you decide this??
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, vgeneTree, vmoduleColors, main = "Network heatmap plot, all genes")
```
![hmpa](/images/vst-network-hmap.png)

There is so much overlap... This makes me worry I guess. Is everything just correlated together?

Or is this not the merged modules? I think it is.. I don't fully understand this I guess.

Visualizing the network of eigengenes
=====================================

"It is often interesting to study the relationships among the found modules. One can use the eigengenes as representative profiles and quantify module similarity by eigengene correlation. The package contains a convenient function plotEigengeneNetworks that generates a summary plot of the eigengene network. It is usually informative to add a clinical trait (or multiple traits) to the eigengenes to see how the traits fit into the eigengene network"

``` r
Effluent = as.data.frame(Traitdat2$Effluent);
names(Effluent) = "effluent"
MET = orderMEs(cbind(vMEs, Effluent))
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
plotDendrograms = FALSE, xLabelsAngle = 90)
```
![dendro](/images/vst-eigen-dendro-ef.png)
![hmap](/images/vst-eigen-hmap-ef.png)

Ok so there are some groupings. They don't really go with the effluent treatment except for greenyellow

``` r
pCO2 = as.data.frame(Traitdat2$pCO2);
names(pCO2) = "pCO2"
MET = orderMEs(cbind(vMEs, pCO2))
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
plotDendrograms = FALSE, xLabelsAngle = 90)
```
![dendro](/images/vst-eigen-dendro-co.png)
![hmap](/images/vst-eigen-hmap-co.png)

Only the slight association with purple. This is another way of showing the groupings. Cool I guess.

<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf>

Gene relationship to trait and important modules: Gene Significance and Module Membership
=========================================================================================

"We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module."

First trying effluent because of the super high correlation it has with the greenyellow module

``` r
# Define variable effluent containing the weight column of datTrait
Effluent = as.data.frame(Traitdat2$Effluent);
names(Effluent) = "effluent"
# names (colors) of the modules
modNames = substring(names(vMEs), 3)
geneModuleMembership = as.data.frame(bicor(tVMostVaryGenes_df, vMEs, robustX = TRUE, robustY = FALSE, use = "pairwise.complete.obs"));
```

    ## Warning in bicor(tVMostVaryGenes_df, vMEs, robustX = TRUE, robustY = FALSE, :
    ## bicor: zero MAD in variable 'x'. Pearson correlation was used for individual
    ## columns with zero (or missing) MAD.

``` r
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), vnSamples)); # correlation of the module eigengene and the expression profile

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(bicor(tVMostVaryGenes_df, Effluent, robustX = TRUE, robustY = FALSE, use = "pairwise.complete.obs"));
```

    ## Warning in bicor(tVMostVaryGenes_df, Effluent, robustX = TRUE, robustY =
    ## FALSE, : bicor: zero MAD in variable 'x'. Pearson correlation was used for
    ## individual columns with zero (or missing) MAD.

``` r
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), vnSamples)); # absolute value of correlation between genes and the effluent treatment
names(geneTraitSignificance) = paste("GS.", names(Effluent), sep="");
names(GSPvalue) = paste("p.GS.", names(Effluent), sep="");
```

I think it used the normal pearson correlation for this. I don't know exactly why I can't do the other correlation...

Intramodular analysis: identifying genes with high GS and MM
============================================================

Using the GS and MM measures, we can identify genes that have a high significance for \[effluent\] as well as high module membership in interesting modules. As an example, we look at the \[greenyellow\] module that has the highest association with \[effluent\]. We plot a scatterplot of Gene Significance vs. Module Membership in the \[greenyellow\] module:

``` r
module = "greenyellow"
column = match(module, modNames);
moduleGenes = vmoduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for Effluent Treatment",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
```
![gy](/images/gy-cor.png)

High values on the y axis means high correlation between the gene and the treatment. High values of the y-axis are high correlation between the module and expression profile. I am starting to understand this maybe. If there is a correlation that suggests that there is a relationship. Compared to a cloud of points...

"The plot is shown in Fig. 2. Clearly, GS and MM are highly correlated, illustrating that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait. The reader is encouraged to try this code with other significance trait/module correlation (for example, the magenta, midnightblue, and red modules with weight)."

Do this for the other modules

Purple module

``` r
module = "purple"
column = match(module, modNames);
moduleGenes = vmoduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for Effluent Treatment",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
```
![purp](/images/purp-cor.png)

Magenta for pCO2

``` r
# Define variable effluent containing the weight column of datTrait
pCO2 = as.data.frame(Traitdat2$pCO2);
names(pCO2) = "pCO2"
# names (colors) of the modules
modNames = substring(names(vMEs), 3)
geneModuleMembership = as.data.frame(cor(tVMostVaryGenes_df, vMEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), vnSamples)); # correlation of the module eigengene and the expression profile

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(tVMostVaryGenes_df, pCO2, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), vnSamples)); # absolute value of correlation between genes and the pCO2 treatment
names(geneTraitSignificance) = paste("GS.", names(pCO2), sep="");
names(GSPvalue) = paste("p.GS.", names(pCO2), sep="");
```

``` r
module = "magenta"
column = match(module, modNames);
moduleGenes = vmoduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for pCO2 Treatment",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
```
![mag](/images/mag-cor.png)

Ok, so this isn't really correlated, even though this is supposed to be correlated significantly with pCO2 treatment. This may be because it's not the same robust correlation measure?

Ok now I want to get the genes that are in the greenyellow module. Who are they??

``` r
GYGeneNames <- names(tVMostVaryGenes_df)[vmoduleColors=="greenyellow"]
print(GYGeneNames)
```

    ##  [1] "MSTRG.13878" "MSTRG.25557" "MSTRG.13690" "MSTRG.16274" "MSTRG.27893"
    ##  [6] "MSTRG.2530"  "MSTRG.22163" "MSTRG.20077" "MSTRG.15392" "MSTRG.8152"
    ## [11] "MSTRG.2984"  "MSTRG.2767"  "MSTRG.11988" "MSTRG.766"   "MSTRG.10053"
    ## [16] "MSTRG.7790"  "MSTRG.2705"  "MSTRG.25599" "MSTRG.14393" "MSTRG.16437"
    ## [21] "MSTRG.12661" "MSTRG.5955"  "MSTRG.17269" "MSTRG.3135"  "MSTRG.23860"
    ## [26] "MSTRG.32832" "MSTRG.1768"  "gene9995"    "MSTRG.435"   "MSTRG.31251"
    ## [31] "gene29006"   "MSTRG.6399"  "MSTRG.28055" "MSTRG.2399"  "MSTRG.16456"
    ## [36] "MSTRG.2457"  "MSTRG.16584" "gene15520"   "MSTRG.31917" "gene15540"  
    ## [41] "MSTRG.16180" "MSTRG.19568" "MSTRG.2908"  "MSTRG.25289" "gene17224"  
    ## [46] "MSTRG.16547" "gene13664"   "gene30051"   "MSTRG.15756" "MSTRG.14902"
    ## [51] "MSTRG.25273" "MSTRG.8362"  "gene14559"   "MSTRG.14531" "MSTRG.26795"
    ## [56] "MSTRG.29582" "MSTRG.28511"

Ok, these are only 57 genes.

Save them for doing stuff in the terminal.

``` r
GYModGeneNames <- as.data.frame(GYGeneNames)
print(GYModGeneNames)
```

    ##    GYGeneNames
    ## 1  MSTRG.13878
    ## 2  MSTRG.25557
    ## 3  MSTRG.13690
    ## 4  MSTRG.16274
    ## 5  MSTRG.27893
    ## 6   MSTRG.2530
    ## 7  MSTRG.22163
    ## 8  MSTRG.20077
    ## 9  MSTRG.15392
    ## 10  MSTRG.8152
    ## 11  MSTRG.2984
    ## 12  MSTRG.2767
    ## 13 MSTRG.11988
    ## 14   MSTRG.766
    ## 15 MSTRG.10053
    ## 16  MSTRG.7790
    ## 17  MSTRG.2705
    ## 18 MSTRG.25599
    ## 19 MSTRG.14393
    ## 20 MSTRG.16437
    ## 21 MSTRG.12661
    ## 22  MSTRG.5955
    ## 23 MSTRG.17269
    ## 24  MSTRG.3135
    ## 25 MSTRG.23860
    ## 26 MSTRG.32832
    ## 27  MSTRG.1768
    ## 28    gene9995
    ## 29   MSTRG.435
    ## 30 MSTRG.31251
    ## 31   gene29006
    ## 32  MSTRG.6399
    ## 33 MSTRG.28055
    ## 34  MSTRG.2399
    ## 35 MSTRG.16456
    ## 36  MSTRG.2457
    ## 37 MSTRG.16584
    ## 38   gene15520
    ## 39 MSTRG.31917
    ## 40   gene15540
    ## 41 MSTRG.16180
    ## 42 MSTRG.19568
    ## 43  MSTRG.2908
    ## 44 MSTRG.25289
    ## 45   gene17224
    ## 46 MSTRG.16547
    ## 47   gene13664
    ## 48   gene30051
    ## 49 MSTRG.15756
    ## 50 MSTRG.14902
    ## 51 MSTRG.25273
    ## 52  MSTRG.8362
    ## 53   gene14559
    ## 54 MSTRG.14531
    ## 55 MSTRG.26795
    ## 56 MSTRG.29582
    ## 57 MSTRG.28511

``` r
write.csv(GYModGeneNames, "GYModGeneNames.csv")
```

Point-Biserial Correlation for Bianary Data
===========================================

Tried a different type of correlation specific to bianary data. It has really similar correlations to bicor (makes sense) but the direction is changed. I'm not sure what I think about that yet. And it did not find anything significant to the pCO2 treatment... Because I went forward only with the greenyellow module so far, I'm not too concerned because it is significantly correlated in both. The directionality though... I could have messed something up in this code though.

Try using Point-Biserial Correlation <https://www.rdocumentation.org/packages/ltm/versions/1.1-1/topics/biserial.cor>

TMM normalized then log2 transformed

``` r
library(ltm)
```

    ## Loading required package: MASS

    ##
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ##
    ##     select

    ## The following object is masked from 'package:genefilter':
    ##
    ##     area

    ## Loading required package: msm

    ## Loading required package: polycor

``` r
# Define numbers of genes and samples
lognGenes = ncol(tlog2orderedTMMM);
lognSamples = nrow(tlog2orderedTMMM);
# Recalculate MEs with color labels
logMEs0 = moduleEigengenes(tlog2orderedTMMM, logmoduleColors)$eigengenes
logMEs = orderMEs(logMEs0)

colnames(logMEs)
```

    ##  [1] "MEdarkgreen"      "MEdarkgrey"       "MEpalevioletred3" "MEblack"         
    ##  [5] "MEbrown"          "MEdarkolivegreen" "MElightgreen"     "MEdarkorange2"   
    ##  [9] "MEplum2"          "MEturquoise"      "MEdarkorange"     "MEskyblue"       
    ## [13] "MEyellowgreen"    "MEorange"         "MEgrey60"         "MEbisque4"       
    ## [17] "MEivory"          "MElightcyan1"     "MElightyellow"    "MEbrown4"        
    ## [21] "MEblue"           "MEorangered4"     "MEgrey"

``` r
# TMM NORMALIZED
# this type of correlation is clunky. I first have to change every module into a vector

Logdarkgreen <- as.vector(logMEs$MEdarkgreen)
Logdarkgrey <- as.vector(logMEs$MEdarkgrey)
Logpalevioletred3<- as.vector(logMEs$MEpalevioletred3)
Logblack <- as.vector(logMEs$MEblack)
Logbrown <- as.vector(logMEs$MEbrown)
Logdarkolivegreen <- as.vector(logMEs$MEdarkolivegreen)
Loglightgreen <- as.vector(logMEs$MElightgreen)
Logdarkorange2 <- as.vector(logMEs$MEdarkorange2)
Logplum2 <- as.vector(logMEs$MEplum2)
Logturquoise <- as.vector(logMEs$MEturquoise)
Logdarkorange <- as.vector(logMEs$MEdarkorange)
Logskyblue <- as.vector(logMEs$MEskyblue)
Logyellowgreen <- as.vector(logMEs$MEyellowgreen)
Logorange <- as.vector(logMEs$MEorange)
Loggrey60 <- as.vector(logMEs$MEgrey60)
Logbisque4 <- as.vector(logMEs$MEbisque4)
Logivory <- as.vector(logMEs$MEivory)
Loglightcyan1 <- as.vector(logMEs$MElightcyan1)
Loglightyellow <- as.vector(logMEs$MElightyellow)
Logbrown4 <- as.vector(logMEs$MEbrown4)
Logblue <- as.vector(logMEs$MEblue)
Logorangered4 <- as.vector(logMEs$MEorangered4)
Loggrey <- as.vector(logMEs$MEgrey) # I dont know why theres still a grey one??

# treatments need to be a vector too
effluent <- as.vector(Traitdat2$Effluent)
cO2 <- as.vector(Traitdat2$pCO2)

# then correlate each individually to each other

#darkgreen
LogdarkgreenCore = biserial.cor(Logdarkgreen, effluent, use = c("all.obs"), level = 1)
LogdarkgreenCorc = biserial.cor(Logdarkgreen, cO2, use = c("all.obs"), level = 1)
#darkgrey
LogdarkgreyCore = biserial.cor(Logdarkgrey, effluent, use = c("all.obs"), level = 1)
LogdarkgreyCorc = biserial.cor(Logdarkgrey, cO2, use = c("all.obs"), level = 1)
#palevioletred3
Logpalevioletred3Core = biserial.cor(Logpalevioletred3, effluent, use = c("all.obs"), level = 1)
Logpalevioletred3Corc = biserial.cor(Logpalevioletred3, cO2, use = c("all.obs"), level = 1)
#black
LogblackCore = biserial.cor(Logblack, effluent, use = c("all.obs"), level = 1)
LogblackCorc = biserial.cor(Logblack, cO2, use = c("all.obs"), level = 1)
#brown
LogbrownCore = biserial.cor(Logbrown, effluent, use = c("all.obs"), level = 1)
LogbrownCorc = biserial.cor(Logbrown, cO2, use = c("all.obs"), level = 1)
#darkolivegreen
LogdarkolivegreenCore = biserial.cor(Logdarkolivegreen, effluent, use = c("all.obs"), level = 1)
LogdarkolivegreenCorc = biserial.cor(Logdarkolivegreen, cO2, use = c("all.obs"), level = 1)
#lightgreen
LoglightgreenCore = biserial.cor(Loglightgreen, effluent, use = c("all.obs"), level = 1)
LoglightgreenCorc = biserial.cor(Loglightgreen, cO2, use = c("all.obs"), level = 1)
#darkorange2
Logdarkorange2Core = biserial.cor(Logdarkorange2, effluent, use = c("all.obs"), level = 1)
Logdarkorange2Corc = biserial.cor(Logdarkorange2, cO2, use = c("all.obs"), level = 1)
#plum2
Logplum2Core = biserial.cor(Logplum2, effluent, use = c("all.obs"), level = 1)
Logplum2Corc = biserial.cor(Logplum2, cO2, use = c("all.obs"), level = 1)
#turquoise
LogturquoiseCore = biserial.cor(Logturquoise, effluent, use = c("all.obs"), level = 1)
LogturquoiseCorc = biserial.cor(Logturquoise, cO2, use = c("all.obs"), level = 1)
#darkorange
LogdarkorangeCore = biserial.cor(Logdarkorange, effluent, use = c("all.obs"), level = 1)
LogdarkorangeCorc = biserial.cor(Logdarkorange, cO2, use = c("all.obs"), level = 1)
#skyblue
LogskyblueCore = biserial.cor(Logskyblue, effluent, use = c("all.obs"), level = 1)
LogskyblueCorc = biserial.cor(Logskyblue, cO2, use = c("all.obs"), level = 1)
#yellowgreen
LogyellowgreenCore = biserial.cor(Logyellowgreen, effluent, use = c("all.obs"), level = 1)
LogyellowgreenCorc = biserial.cor(Logyellowgreen, cO2, use = c("all.obs"), level = 1)
#orange
LogorangeCore = biserial.cor(Logorange, effluent, use = c("all.obs"), level = 1)
LogorangeCorc = biserial.cor(Logorange, cO2, use = c("all.obs"), level = 1)
#grey60
Loggrey60Core = biserial.cor(Loggrey60, effluent, use = c("all.obs"), level = 1)
Loggrey60Corc = biserial.cor(Loggrey60, cO2, use = c("all.obs"), level = 1)
#bisque4
Logbisque4Core = biserial.cor(Logbisque4, effluent, use = c("all.obs"), level = 1)
Logbisque4Corc = biserial.cor(Logbisque4, cO2, use = c("all.obs"), level = 1)
#ivory
LogivoryCore = biserial.cor(Logivory, effluent, use = c("all.obs"), level = 1)
LogivoryCorc = biserial.cor(Logivory, cO2, use = c("all.obs"), level = 1)
#lightcyan1
Loglightcyan1Core = biserial.cor(Loglightcyan1, effluent, use = c("all.obs"), level = 1)
Loglightcyan1Corc = biserial.cor(Loglightcyan1, cO2, use = c("all.obs"), level = 1)
#lightyellow
LoglightyellowCore = biserial.cor(Loglightyellow, effluent, use = c("all.obs"), level = 1)
LoglightyellowCorc = biserial.cor(Loglightyellow, cO2, use = c("all.obs"), level = 1)
#brown4
Logbrown4Core = biserial.cor(Logbrown4, effluent, use = c("all.obs"), level = 1)
Logbrown4Corc = biserial.cor(Logbrown4, cO2, use = c("all.obs"), level = 1)
#blue
LogblueCore = biserial.cor(Logblue, effluent, use = c("all.obs"), level = 1)
LogblueCorc = biserial.cor(Logblue, cO2, use = c("all.obs"), level = 1)
#orangered4
Logorangered4Core = biserial.cor(Logorangered4, effluent, use = c("all.obs"), level = 1)
Logorangered4Corc = biserial.cor(Logorangered4, cO2, use = c("all.obs"), level = 1)
#grey
LoggreyCore = biserial.cor(Loggrey, effluent, use = c("all.obs"), level = 1)
LoggreyCorc = biserial.cor(Loggrey, cO2, use = c("all.obs"), level = 1)

# then bring all those values together

Logeffluentcor <- list(LogdarkgreenCore, LogdarkgreyCore, Logpalevioletred3Core, LogblackCore, LogbrownCore, LogdarkolivegreenCore, LoglightgreenCore, Logdarkorange2Core, Logplum2Core, LogturquoiseCore, LogdarkorangeCore, LogskyblueCore, LogyellowgreenCore, LogorangeCore, Loggrey60Core, Logbisque4Core, LogivoryCore, Loglightcyan1Core, LoglightyellowCore, Logbrown4Core, LogblueCore, Logorangered4Core, LoggreyCore)

LogpCO2cor <- list(LogdarkgreenCorc, LogdarkgreyCorc, Logpalevioletred3Corc, LogblackCorc, LogbrownCorc, LogdarkolivegreenCorc, LoglightgreenCorc, Logdarkorange2Corc, Logplum2Corc, LogturquoiseCorc, LogdarkorangeCorc, LogskyblueCorc, LogyellowgreenCorc, LogorangeCorc, Loggrey60Corc, Logbisque4Corc, LogivoryCorc, Loglightcyan1Corc, LoglightyellowCorc, Logbrown4Corc, LogblueCorc, Logorangered4Corc, LoggreyCorc)

# make a dataframe from the lists
Logcor <- do.call(rbind, Map(data.frame, LogEffulentCor=Logeffluentcor, LogpCO2Cor=LogpCO2cor))
# give them rownames
rownames(Logcor) <- c("MEdarkgreen", "MEdarkgrey", "MEpalevioletred3", "MEblack", "MEbrown", "MEdarkolivegreen", "MElightgreen", "MEdarkorange2", "MEplum2", "MEturquoise", "MEdarkorange", "MEskyblue", "MEyellowgreen", "MEorange", "MEgrey60", "MEbisque4", "MivoryE", "MElightcyan1", "MElightyellow", "MEbrown4", "MEblue", "MEorangered4", "MEgrey")

# make it a matrix so I could do the same p value calculation as above
Logmatrixcor <- data.matrix(Logcor, rownames.force = NA)
LogcorPvalue = corPvalueStudent(Logmatrixcor, lognSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
LogtextMatrix = paste(signif(Logmatrixcor, 2), "\n(",
signif(LogcorPvalue, 1), ")", sep = "");
dim(LogtextMatrix) = dim(Logmatrixcor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = Logmatrixcor,
xLabels = names(Logcor),
yLabels = names(logMEs),
ySymbols = names(logMEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = LogtextMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships biserialcor"))
```
![biscor](/images/tmm-biscor.png)
Ok, I think this is the same as the other correlation method.

Vst transformed point baserial correlation
==========================================

``` r
# Define numbers of genes and samples
vnGenes = ncol(tVMostVaryGenes_df);
vnSamples = nrow(tVMostVaryGenes_df);
# Recalculate MEs with color labels
vMEs0 = moduleEigengenes(tVMostVaryGenes_df, vmoduleColors)$eigengenes
vMEs = orderMEs(vMEs0)
colnames(vMEs)
```

    ##  [1] "MEbrown"       "MEmagenta"     "MEblue"        "MEpurple"     
    ##  [5] "MEgreenyellow" "MEtan"         "MEsalmon"      "MEpink"       
    ##  [9] "MEred"         "MEyellow"      "MEblack"       "MEgreen"

``` r
# VST TRANSFORMED

vbrown <- as.vector(vMEs$MEbrown)
vmagenta <- as.vector(vMEs$MEmagenta)
vblue <- as.vector(vMEs$MEblue)
vpurple <- as.vector(vMEs$MEpurple)
vgreenyellow <- as.vector(vMEs$MEgreenyellow)
vtan <- as.vector(vMEs$MEtan)
vsalmon <- as.vector(vMEs$MEsalmon)
vpink <- as.vector(vMEs$MEpink)
vred <- as.vector(vMEs$MEred)
vyellow <- as.vector(vMEs$MEyellow)
vblack <- as.vector(vMEs$MEblack)
vgreen <- as.vector(vMEs$MEgreen)


effluent <- as.vector(Traitdat2$Effluent)
cO2 <- as.vector(Traitdat2$pCO2)

#brown
vbrownCore = biserial.cor(vbrown, effluent, use = c("all.obs"), level = 1)
vbrownCorc = biserial.cor(vbrown, cO2, use = c("all.obs"), level = 1)
#magenta
vmagentaCore = biserial.cor(vmagenta, effluent, use = c("all.obs"), level = 1)
vmagentaCorc = biserial.cor(vmagenta, cO2, use = c("all.obs"), level = 1)
#blue
vblueCore = biserial.cor(vblue, effluent, use = c("all.obs"), level = 1)
vblueCorc = biserial.cor(vblue, cO2, use = c("all.obs"), level = 1)
#purple
vpurpleCore = biserial.cor(vpurple, effluent, use = c("all.obs"), level = 1)
vpurpleCorc = biserial.cor(vpurple, cO2, use = c("all.obs"), level = 1)
#greenyellow
vgreenyellowCore = biserial.cor(vgreenyellow, effluent, use = c("all.obs"), level = 1)
vgreenyellowCorc = biserial.cor(vgreenyellow, cO2, use = c("all.obs"), level = 1)
#tan
vtanCore = biserial.cor(vtan, effluent, use = c("all.obs"), level = 1)
vtanCorc = biserial.cor(vtan, cO2, use = c("all.obs"), level = 1)
#salmon
vsalmonCore = biserial.cor(vsalmon, effluent, use = c("all.obs"), level = 1)
vsalmonCorc = biserial.cor(vsalmon, cO2, use = c("all.obs"), level = 1)
#pink
vpinkCore = biserial.cor(vpink, effluent, use = c("all.obs"), level = 1)
vpinkCorc = biserial.cor(vpink, cO2, use = c("all.obs"), level = 1)
#red
vredCore = biserial.cor(vred, effluent, use = c("all.obs"), level = 1)
vredCorc = biserial.cor(vred, cO2, use = c("all.obs"), level = 1)
#yellow
vyellowCore = biserial.cor(vyellow, effluent, use = c("all.obs"), level = 1)
vyellowCorc = biserial.cor(vyellow, cO2, use = c("all.obs"), level = 1)
#black
vblackCore = biserial.cor(vblack, effluent, use = c("all.obs"), level = 1)
vblackCorc = biserial.cor(vblack, cO2, use = c("all.obs"), level = 1)
#green
vgreenCore = biserial.cor(vgreen, effluent, use = c("all.obs"), level = 1)
vgreenCorc = biserial.cor(vgreen, cO2, use = c("all.obs"), level = 1)


Veffluentcor <- list(vbrownCore, vmagentaCore, vblueCore, vpurpleCore, vgreenyellowCore, vtanCore, vsalmonCore, vpinkCore, vredCore, vyellowCore, vblackCore, vgreenCore)

VpCO2cor <- list(vbrownCorc, vmagentaCorc, vblueCorc, vpurpleCorc, vgreenyellowCorc, vtanCorc, vsalmonCorc, vpinkCorc, vredCorc, vyellowCorc, vblackCorc, vgreenCorc)

Vcor <- do.call(rbind, Map(data.frame, vEffulentCor=Veffluentcor, vpCO2Cor=VpCO2cor))
rownames(Vcor) <- c("brown", "magenta", "blue", "purple", "greenyellow", "tan", "salmon", "pink", "red", "yellow", "black", "green")

Vmatrixcor <- data.matrix(Vcor, rownames.force = NA)
VcorPvalue = corPvalueStudent(Vmatrixcor, vnSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
vtextMatrix = paste(signif(Vmatrixcor, 2), "\n(",
signif(VcorPvalue, 1), ")", sep = "");
dim(vtextMatrix) = dim(Vmatrixcor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = Vmatrixcor,
xLabels = names(Vcor),
yLabels = names(vMEs),
ySymbols = names(vMEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = vtextMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships biserialcor"))
```
![biscor](/images/vst-biscor.png)

These are the same modules but slightly different. The greenyellow relationship is different...
