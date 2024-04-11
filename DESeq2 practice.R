#Load libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

install.packages("tidyverse")
library(tidyverse)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("airway")
library(airway)

#Preparing data--> Creating two CSV files (sample info and counts)
data("airway")
sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

#Reading in counts data
counts_data <- read.csv('counts_data.csv')

#Reading in sample info
colData <- read.csv('sample_info.csv')

#Making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

#Are they in the same order?
all(colnames(counts_data) == rownames(colData))

#Construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData  = counts_data,
                       colData = colData,
                       design = ~ dexamethasone)

#Pre-Filtering the data to omit rows with low gene counts (<10)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Setting the factor level (what is the baseline factor to compare the others to?)
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

#Collapse technical replicates if they exist in data
#N/A

#Running DESeq
dds <- DESeq(dds)
res <- results(dds)

#Viewing a summary of the results. Adjusting the alpha level to 0.01 (significance cut-off)
res0.01 <- results(dds, alpha=0.01)
summary(res0.01)

#MA plot (blue points are significant and differentially expressed genes)
plotMA(res)
