if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("biomaRt")
BiocManager::install("mixOmics")


library(edgeR) #https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
library(biomaRt)
library(mixOmics)

setwd("~/Desktop/RP 2/endometrium analysis")

Countdata<-read.delim("GSE132711_RNAse.txt")
row.names(Countdata)= Countdata$Gene_id
Countdata=Countdata[,-1]
Countdata<-Countdata[,-c(1:11)]

samples<-read.delim("Sample.txt")
