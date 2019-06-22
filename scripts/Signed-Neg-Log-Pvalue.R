#Preparing log pvalue with downregulated negative fo GO_MWU
#Maggie Schedl
#2019

#For preparing results from DESeq2 for use with GO_MWU https://github.com/z0on/GO_MWU
#The analysis requires a text file with one column of gene names and one column of "raw" p-values. 
#But they really mean p-values "in the form of "signed negative log p-values". These measures are negative decimal logarithms of the raw (uncorrected) p-value for each gene, multiplied by -1 if the gene was down-regulated." 
#We took that to mean the -log of the pvalue, and then to get the down-regulated genes to be negative, the solution was to multiply by the log2foldchange, then divide by the absolute value of the log2foldchange



setwd(/home/mschedl/Working-CASE-RNA/histat/stringtie/restring)

### CA
results_CA_F_text <- read.delim("results_CA_F.txt") #read in table of full results dataset of CA vs CON contrast
log_results_CA_F <- -log(results_CA_F$pvalue) #negative log of p value (NOT ADJUSTED)
logCA <- data.frame(matrix(unlist(log_results_CA_F), nrow=length(log_results_CA_F), byrow=T)) #make the list into a data frame
colnames(logCA) <- c("logpvalue") #give it a column name
results_CA_F_log <- cbind(results_CA_F_text, logCA) #add that column to the results dataframe
absvallogfoldchange <- abs(results_CA_F_log$log2FoldChange) #take the absolute value for the log2foldchange
abschange <- data.frame(matrix(unlist(absvallogfoldchange), nrow=length(absvallogfoldchange), byrow=T)) #make that list into a data frame

colnames(abschange) <- c("absolute_value_log_2change") #give it a column name

results_CA_F_log <- cbind(results_CA_F_log, abschange) # add it to the results dataframe
results_CA_F_log$multilogpvalue <- results_CA_F_log$logpvalue * results_CA_F_log$log2FoldChange #multiply the -log pvalue by the log2fold change
results_CA_F_log[, "neglogpvalue"] <- results_CA_F_log[, "multilogpvalue"] / results_CA_F_log[, "absolute_value_log_2change"] #divide the mulitpled -logpvalue by the absolute value of the log2fold change

write.table(results_CA_F_log, "results_CA_F_log.txt",quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t") #write into a text file for GO_MWU input 

###CASE
results_CASE_F_text <- read.delim("results_CASE_F.txt") #read in table of full results dataset of CASE vs CON contrast
log_results_CASE_F <- -log(results_CASE_F$pvalue) #negative log of p value (NOT ADJUSTED)
logCASE <- data.frame(matrix(unlist(log_results_CASE_F), nrow=length(log_results_CASE_F), byrow=T)) #make the list into a data frame
colnames(logCASE) <- c("logpvalue") #give it a column name
results_CASE_F_log <- cbind(results_CASE_F_text, logCASE) #add that column to the results dataframe
absvallogfoldchangeCASE <- abs(results_CASE_F_log$log2FoldChange) #take the absolute value for the log2foldchange
abschangeCASE <- data.frame(matrix(unlist(absvallogfoldchangeCASE), nrow=length(absvallogfoldchangeCASE), byrow=T)) #make that list into a data frame

colnames(abschangeCASE) <- c("absolute_value_log_2change") #give it a column name

results_CASE_F_log <- cbind(results_CASE_F_log, abschangeCASE) # add it to the results dataframe
results_CASE_F_log$multilogpvalue <- results_CASE_F_log$logpvalue * results_CASE_F_log$log2FoldChange #multiply the -log pvalue by the log2fold change
results_CASE_F_log[, "neglogpvalue"] <- results_CASE_F_log[, "multilogpvalue"] / results_CASE_F_log[, "absolute_value_log_2change"] #divide the mulitpled -logpvalue by the absolute value of the log2fold change

write.table(results_CASE_F_log, "results_CASE_F_log.txt",quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t") #write into a text file for GO_MWU input 


####SE
results_SE_F_text <- read.delim("results_SE_F.txt") #read in table of full results dataset of SE vs CON contrast
log_results_SE_F <- -log(results_SE_F$pvalue) #negative log of p value (NOT ADJUSTED)
logSE <- data.frame(matrix(unlist(log_results_SE_F), nrow=length(log_results_SE_F), byrow=T)) #make the list into a data frame
colnames(logSE) <- c("logpvalue") #give it a column name
results_SE_F_log <- cbind(results_SE_F_text, logSE) #add that column to the results dataframe
absvallogfoldchangeSE <- abs(results_SE_F_log$log2FoldChange) #take the absolute value for the log2foldchange
abschangeSE <- data.frame(matrix(unlist(absvallogfoldchangeSE), nrow=length(absvallogfoldchangeSE), byrow=T)) #make that list into a data frame

colnames(abschangeSE) <- c("absolute_value_log_2change") #give it a column name

results_SE_F_log <- cbind(results_SE_F_log, abschangeSE) # add it to the results dataframe
results_SE_F_log$multilogpvalue <- results_SE_F_log$logpvalue * results_SE_F_log$log2FoldChange #multiply the -log pvalue by the log2fold change
results_SE_F_log[, "neglogpvalue"] <- results_SE_F_log[, "multilogpvalue"] / results_SE_F_log[, "absolute_value_log_2change"] #divide the mulitpled -logpvalue by the absolute value of the log2fold change

write.table(results_SE_F_log, "results_SE_F_log.txt",quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t") #write into a text file for GO_MWU input 

