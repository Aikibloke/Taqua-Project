if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
aBiocManager::install("biomaRt")
BiocManager::install("mixOmics")


library(edgeR) #https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
library(biomaRt)
library(mixOmics)

setwd("~/Desktop/RP 2/Endometrium analysis/Chi data analysis ")

Countdata<-read.delim("GSE132711_RNAse.txt")
row.names(Countdata)= Countdata$Gene_id
Countdata=Countdata[,-1]
Countdata<-Countdata[,-c(1:11)]

samples<-read.delim("Sample.txt")

#Create DGE object
DGE<-DGEList(counts=Countdata,group=samples$Type)

#Filter low expression genes 
keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] #Remove low expression genes

#Calculate effective library sizes
DGE <- calcNormFactors(DGE)

#Identifiy DEGs
DGE <- estimateDisp(DGE) #Estimate dispersions
et <- exactTest(DGE) #Calculate DE
significant<-topTags(et,n = 6637,p.value = 0.05)

write.csv(significant$table)





