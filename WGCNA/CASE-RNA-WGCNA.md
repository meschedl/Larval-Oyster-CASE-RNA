CASE-RNA-WGCNA
================
Maggie Schedl
1/30/2020

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
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colMeans, colnames, colSums, dirname, do.call, duplicated,
    ##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
    ##     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
    ##     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
    ##     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which, which.max,
    ##     which.min

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

start with data input and cleaning <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf>

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

preliminary filtering of the counts with P over A, Meaning that I wanted to remove all rows that have less than 21.4% of the samples with less than 5 counts. 21.4% was chosen because that is 3/14, or the minimum ratio of samples per one treatment.

``` r
###filtering values for PoverA #set filter values for PoverA, P percent of the samples have counts over A 
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

    ##             CASE_J03 CASE_J09 CASE_J12 CASE_J13 CA_J06 CA_J08 CA_J11
    ## MSTRG.10383      198      306      198      235    190    133    180
    ## MSTRG.28362       18        7       18       18      3      4     11
    ## MSTRG.10380        5       14        9       19      6     23     12
    ## MSTRG.32256        4        4        4        8     15      0     10
    ## MSTRG.18313       92      169      165      206    120    129    134
    ## MSTRG.10381      143      146       90      124     65     42     60
    ## MSTRG.19848       42       74       21       55     51     32     29
    ## MSTRG.19849       50       38       37       75     21     33     43
    ## MSTRG.19846       26       94       79       78     44     59     31
    ## MSTRG.19847       12       27        4       11     26     11      3
    ##             CA_J18 CON_J02 CON_J05 CON_J10 SE_J01 SE_J04 SE_J07
    ## MSTRG.10383    176     131     139     276    179    104    259
    ## MSTRG.28362     30      16       7       8     13      5      2
    ## MSTRG.10380      7      26       5       3     11      2     47
    ## MSTRG.32256      7       5       5       0      7     14      5
    ## MSTRG.18313    236     226     340     374    159     86    313
    ## MSTRG.10381     91      69      76     109     45     63     46
    ## MSTRG.19848     59      49      24      24     33     18     45
    ## MSTRG.19849     52      49      30      60     37     50     39
    ## MSTRG.19846     59      80      22      12     46     30     79
    ## MSTRG.19847      2       0       3       0      4      5      7

I need to use normalized count data, so use the DESeq2 method to do this.

Ok to get normalized count data I have to make it into a DESeq dataset, use their function to normalize the counts, then extract out the normalized counts from the DESeq format, and write them into a table that will then be moved forward. So basically what I have done above is useless and I'll have to do over after using the normalized counts. Starting with GeneCountData\_Filt because that is the same actually. Using code from previous markdown file

<https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html>

``` r
# beginning part copy and pasted from DESeq2 rmarkdown
# Make sure things are in the correct format

CASE_treatment <- read.csv("treatment_data.csv", header=TRUE, sep=",") # need treatment info for formatting purposes
rownames(CASE_treatment) <- CASE_treatment$sample
colnames(GeneCountData_Filt) <- CASE_treatment$sample
head(CASE_treatment)
```

    ##            sample treatment library extraction
    ## CASE_J03 CASE_J03      CASE   three        two
    ## CASE_J09 CASE_J09      CASE    four        two
    ## CASE_J12 CASE_J12      CASE     two      three
    ## CASE_J13 CASE_J13      CASE     two      three
    ## CA_J06     CA_J06        CA   three        two
    ## CA_J08     CA_J08        CA     one        two

``` r
head(GeneCountData_Filt)
```

    ##             CASE_J03 CASE_J09 CASE_J12 CASE_J13 CA_J06 CA_J08 CA_J11
    ## MSTRG.10383      198      306      198      235    190    133    180
    ## MSTRG.28362       18        7       18       18      3      4     11
    ## MSTRG.10380        5       14        9       19      6     23     12
    ## MSTRG.32256        4        4        4        8     15      0     10
    ## MSTRG.18313       92      169      165      206    120    129    134
    ## MSTRG.10381      143      146       90      124     65     42     60
    ##             CA_J18 CON_J02 CON_J05 CON_J10 SE_J01 SE_J04 SE_J07
    ## MSTRG.10383    176     131     139     276    179    104    259
    ## MSTRG.28362     30      16       7       8     13      5      2
    ## MSTRG.10380      7      26       5       3     11      2     47
    ## MSTRG.32256      7       5       5       0      7     14      5
    ## MSTRG.18313    236     226     340     374    159     86    313
    ## MSTRG.10381     91      69      76     109     45     63     46

``` r
write.csv(GeneCountData_Filt,"filtered_genecount.csv")

write.table(GeneCountData_Filt,"filtered_genecount.txt",quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t") #these were written out just in case they were needed for later

#The row and column names for the two data frames need to be exactly the same for the rest of the analysis, so it is good to check
all(rownames(CASE_treatment) %in% colnames(GeneCountData_Filt))  #Should return TRUE
```

    ## [1] TRUE

``` r
all(rownames(CASE_treatment) == colnames(GeneCountData_Filt))    # should return TRUE
```

    ## [1] TRUE

``` r
# Create the matrix to be able to use as an input for normalizing, not really concerned about the design

CASE_deseq_Matrix <- DESeqDataSetFromMatrix(countData = GeneCountData_Filt,
                              colData = CASE_treatment,
                              design = ~ treatment ) # column name of the treatment information as teh design 
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

``` r
# To perform the median of ratios method of normalization, DESeq2 has a single estimateSizeFactors() function that will generate size factors for us
SizeFact_deseq_Matrix <- estimateSizeFactors(CASE_deseq_Matrix)

# We can take a look at the normalization factor applied to each sample using:
sizeFactors(SizeFact_deseq_Matrix)
```

    ##  CASE_J03  CASE_J09  CASE_J12  CASE_J13    CA_J06    CA_J08    CA_J11 
    ## 0.9881839 1.1900647 0.9358968 1.2312619 0.8273069 0.7993982 0.8893150 
    ##    CA_J18   CON_J02   CON_J05   CON_J10    SE_J01    SE_J04    SE_J07 
    ## 1.2133347 1.0530743 1.0044657 1.1316559 1.0299238 0.8934197 1.1279059

``` r
# Ok there doesn't seem to be any bias on the size factor normalization for any of the specific treatments. Not exactually sure what they do with these values though. 

# Now, to retrieve the normalized counts matrix from dds, we use the counts() function and add the argument normalized=TRUE
normalized_Genecounts <- counts(SizeFact_deseq_Matrix, normalized=TRUE)

# We can save this normalized data matrix to file
write.table(normalized_Genecounts, "normalized_Genecounts.txt", sep="\t", quote=F, col.names=NA)
```

I am pretty sure (but not completely) that SizeFact\_deseq\_Matrix is the same matrix, but with normalized count values instead of raw count values. DESeq2 documentation is not completely clear.

I also need to limit the number of genes I am working with. I ended up chosing the top 5000 most variable genes/transcripts. I also did a varience stablizing transformation on the normalized data, because that is what is suggested before looking at ranking most variable genes, and also can be used for input WGCNA data. The major caviate is that vst transformation is more for sample sizes over 30, the suggested one to use is rlog transformation. I can get the rlog transformation to give an ok R2 value if I chose 3000 most variable genes, but that's lower than I wanted and also just barely is over 0.8. Without any transformation it the power is only acceptable at 14, which I guess is high, and that the r suqared sometimes is negative which is weird. So I went with the vst but I'm not super happy about it.

``` r
vSizeFact_deseq_Matrix <- vst(SizeFact_deseq_Matrix, blind = FALSE) # use the normalized matrix.. but I'm not sure exactly what this matrix is
head(assay(vSizeFact_deseq_Matrix), 3)
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
MostVaryGenes3 <- head(order(rowVars(assay(vSizeFact_deseq_Matrix)), decreasing = TRUE), 5000)

# make that information into a matrix
MostVaryGenes3_Mat <- assay(vSizeFact_deseq_Matrix)[ MostVaryGenes3, ] # don't know how this works 
# make that matrix a dataframe 
MostVaryGenes3_df <- as.data.frame(MostVaryGenes3_Mat)

tMostVaryGenes3_df <- as.data.frame(t(as.matrix(MostVaryGenes3_df))) #transpose to right format

#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.

sampleTree = hclust(dist(tMostVaryGenes3_df), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
```

Great there are no obvious outlier samples. Now to pick the power.

``` r
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(tMostVaryGenes3_df, powerVector = powers, verbose = 5)
```

    ## pickSoftThreshold: will use block size 5000.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 5000 of 5000

    ## Warning: executing %dopar% sequentially: no parallel backend registered

    ##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1    0.555  0.593         0.4360  1960.0   2090.00   2930
    ## 2      2    0.152 -0.187        -0.0688  1090.0   1130.00   2100
    ## 3      3    0.527 -0.413         0.3920   708.0    689.00   1630
    ## 4      4    0.676 -0.520         0.5840   502.0    451.00   1320
    ## 5      5    0.789 -0.579         0.7300   377.0    307.00   1100
    ## 6      6    0.860 -0.624         0.8260   294.0    216.00    943
    ## 7      7    0.858 -0.686         0.8390   236.0    156.00    831
    ## 8      8    0.852 -0.738         0.8510   194.0    115.00    742
    ## 9      9    0.836 -0.797         0.8510   162.0     86.30    673
    ## 10    10    0.837 -0.839         0.8670   137.0     65.50    615
    ## 11    12    0.808 -0.938         0.8610   102.0     39.40    524
    ## 12    14    0.801 -1.000         0.8720    78.7     24.60    454
    ## 13    16    0.805 -1.060         0.8860    62.2     16.40    398
    ## 14    18    0.799 -1.120         0.8890    50.2     10.90    353
    ## 15    20    0.809 -1.150         0.9020    41.3      7.60    316
    ## 16    22    0.816 -1.190         0.9140    34.4      5.44    285
    ## 17    24    0.836 -1.210         0.9370    29.0      4.06    258
    ## 18    26    0.839 -1.230         0.9410    24.7      2.99    235
    ## 19    28    0.849 -1.240         0.9490    21.2      2.29    214
    ## 20    30    0.839 -1.270         0.9430    18.3      1.76    197

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

I want to chose 6 because the R2 is 0.86 which is pretty acceptable. The mean connectivity is also really high. I am actually worried about it being that high.

We now calculate the adjacencies, using the soft thresholding power 6:

``` r
softPower = 6;
adjacency = adjacency(tMostVaryGenes3_df, power = softPower);
```

The adjacencies are just correlations of one gene with another, until it does all of the comparissons. But to get a larger network view interconnectedness and similarities of one connection with another need to be taken into account. I need a topological overlap matrix. "More specifically, genes are said to have high topological overlap if they are connected to roughly the same group of genes in the network (i.e. they share the same neighborhood). To calculate the topological overlap for a pair of genes, their connections with all other genes in the network are compared"

To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity:

``` r
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
```

    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.

``` r
dissTOM = 1-TOM
```

<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/WORKSHOP/2014/Langfelder-NetworkDay-clustering.pdf>

We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes. Note that we use the function hclust that provides a much faster hierarchical clustering routine than the standard hclust function.

``` r
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Average linkage: average the dissimilarities between all objects
#single (minimum dissimilarity) and complete (maximum dissimilarity) are other options, but seem too harsh

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04);
```

"The clustering dendrogram plotted by the last command is shown in Figure 2. In the clustering tree (dendrogram), each leaf, that is a short vertical line, corresponds to a gene. Branches of the dendrogram group together densely interconnected, highly co-expressed genes."

"Module identification amounts to the identification of individual branches (”cutting the branches off the dendrogram”)"

In general, I can see really only 2 branches, although there is some separation in that big chunk. I wonder if the connectivity is too high. I just don't know. I'm going to start with the tutorial suggestion of module trimming. So I think this means the minimum size module will have at least 30 genes in it.

Want to dynamically cut the tree height because no single cut height captures the branches "Branches are followed bottom to top. When two branches merge, they are evaluated using shape criteria such as minimum number of objects (genes), their core scatter and the gap between the branches" So the most similar are what are clustered first and then followed up the branches The core areas are the bottom ends of branches where the most similarly expressed genes are. The gap is the height difference between the cut height and the bottom of the branch. The core scatter is the amount of horizontal space between the branch where the tree is cut and the core section.

"deepSplit controls how finely the branches should be split Higher values give more smaller modules, low values (0) give fewer larger modules"

``` r
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
```

    ##  ..cutHeight not given, setting it to 0.986  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.

``` r
table(dynamicMods)
```

    ## dynamicMods
    ##    0    1    2    3    4    5    6    7    8    9   10   11   12 
    ##   31 2204 1187  632  180  140  133  128  102   76   71   69   47

12 modules, 31 unassigned genes. Module 1 and 2 are huge, others are smaller

``` r
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
```

    ## dynamicColors
    ##       black        blue       brown       green greenyellow        grey 
    ##         128        1187         632         140          69          31 
    ##     magenta        pink      purple         red         tan   turquoise 
    ##          76         102          71         133          47        2204 
    ##      yellow 
    ##         180

``` r
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
```

Hmmm. Well I thought there would be 12 colors

"The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation:"

``` r
# Calculate eigengenes
MEList = moduleEigengenes(tMostVaryGenes3_df, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

MEDissThres = 0.3 # chose to cut here
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(tMostVaryGenes3_df, dynamicColors, cutHeight = MEDissThres, verbose = 3)
```

    ##  mergeCloseModules: Merging modules whose distance is less than 0.3
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 13 module eigengenes in given set.
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 11 module eigengenes in given set.
    ##    Calculating new MEs...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 11 module eigengenes in given set.

``` r
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
```

``` r
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
#dev.off()
```

Relating modules to external information and identifying important genes
========================================================================

<https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf>

set up

``` r
#load in trait data, same as for DESeq2 but also with mortality data from Melati
Traitdat <- as.data.frame(read.csv("Trait_Data.csv"))

Traitdat2 <- Traitdat[,-1] #make column 1 into row names
rownames(Traitdat2) <- Traitdat[,1]

#set up for using the module colors in further code

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
```

In this analysis we would like to identify modules that are significantly associated with the measured clinical traits. Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external traits and look for the most significant associations:

``` r
# Define numbers of genes and samples
nGenes = ncol(tMostVaryGenes3_df);
nSamples = nrow(tMostVaryGenes3_df);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(tMostVaryGenes3_df, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, Traitdat2, use = "pairwise.complete.obs"); #will only calculate correlation between two variables for those pairs where both are not-NA
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
```

``` r
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(Traitdat2),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
```

    ## Warning in greenWhiteRed(50): WGCNA::greenWhiteRed: this palette is not suitable for people
    ## with green-red color blindness (the most common kind of color blindness).
    ## Consider using the function blueWhiteRed instead.

Ok, so nothing really correlates with treatment. Because it's taking treatment in the wrong way??? OMG I should separate out the treatments into their actual conditions because assigning numbers to the treatments means nothing???/? Also not sure about having the library and extraction info in there.

``` r
#load in different trait data
Traitdat3 <- as.data.frame(read.csv("Trait_Data.csv"))

Traitdat4 <- Traitdat3[,-1] #make the 
rownames(Traitdat4) <- Traitdat3[,1]

head(Traitdat4)
```

    ##          library extraction T1_Larvae T24_Larvae Percent_Decrease pCO2
    ## CASE_J03       3          2     117.0      55.08         0.529231 2800
    ## CASE_J09       4          2        NA         NA               NA 2800
    ## CASE_J12       2          3     202.5      85.68         0.576889 2800
    ## CASE_J13       2          3        NA         NA               NA 2800
    ## CA_J06         3          2     127.0      44.64         0.648504 2800
    ## CA_J08         1          2     140.0      72.72         0.480571 2800
    ##          Effluent
    ## CASE_J03     0.05
    ## CASE_J09     0.05
    ## CASE_J12     0.05
    ## CASE_J13     0.05
    ## CA_J06       0.00
    ## CA_J08       0.00

"In this analysis we would like to identify modules that are significantly associated with the measured traits. Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external traits and look for the most significant associations:"

Correlate modules with the treatment and mortality data. I don't actually understand how this works at all, does it use the sample information or not?

``` r
# Define numbers of genes and samples
nGenes = ncol(tMostVaryGenes3_df);
nSamples = nrow(tMostVaryGenes3_df);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(tMostVaryGenes3_df, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, Traitdat4, use = "pairwise.complete.obs"); #will only calculate correlation between two variables for those pairs where both are not-NA
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
```

"Since we have a moderately large number of modules and traits, a suitable graphical representation will help in reading the table. We color code each association by the correlation value:"

``` r
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(Traitdat4),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
```

"Each row corresponds to a module eigengene, column to a trait. Each cell contains the corresponding correlation and p-value." ok so there's only one sig correlation for either of the treatments cool cool cool cool cool cool cool. and it's negative, what does that mean?

"We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module."

Ok I'm just going to do effluent anyways

``` r
# Define variable weight containing the weight column of datTrait
effluent = as.data.frame(Traitdat4$Effluent);
names(effluent) = "effluent"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(tMostVaryGenes3_df, MEs, use = "pairwise.complete.obs"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(tMostVaryGenes3_df, effluent, use = "pairwise.complete.obs"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(effluent), sep="");
names(GSPvalue) = paste("p.GS.", names(effluent), sep="")
```

"Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules. As an example, we look at the brown module that has the highest association with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:"

Ok I'll use the greenyellow module and effluent

``` r
module = "greenyellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for effluent treatment",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
```

correlation: low. but p value also low???

I think I should try others I guess. I will stop here for now.

notes
=====

ok with rlog 3000 its ok at 7

with vst 5000 its pretty good at 7

if I don't transform it its ok for 5000 its ok at 14 but it's weird and has negative values
