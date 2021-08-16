library(entropy)
library(edgeR)
library(ggplot2)
library(gplots)
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(DESeq2)
library(MAST)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#setwd("~/Desktop/RP 2/Single cell analysis  ")

#Importing data 
Var<-read.csv("var.csv")
Count <- read.csv("X.csv", col_names = FALSE)
Samples<-read.csv("Sample.csv")
#Samples<-read.csv("Samples2.csv")

Count<- t(Count)
rownames(Count)<- var$...1


#group 

Samples$des<-paste(Samples$Stage, Samples$Epithelial.celltype, sep=".")
Samples$group<-NA

Samples$group[which(Samples$Stage=="late.secretory")]<- "irrl"
Samples$group[which(Samples$Stage=="mid.secretory")]<- "irrl"
Samples$group[which(Samples$Stage=="proliferative")]<- "irrl"
Samples$group[which(Samples$Stage=="early.secretory")]<- "irrl"
Samples$group[which(Samples$Stage=="early.mid.secretory")]<- "irrl"

Samples$group[which(Samples$des=="proliferative.Ciliated")]<- "Y"
Samples$group[which(Samples$des=="proliferative.Lumenal1")]<- "Y"
Samples$group[which(Samples$des=="proliferative.Glandular")]<- "Y"
Samples$group[which(Samples$des=="mid.secretory.Glandular")]<- "Y"
Samples$group[which(Samples$des=="mid.secretory.Lumenal1")]<- "Y"
Samples$group[which(Samples$des=="mid.secretory.Ciliated")]<- "Y"
Samples$group[which(Samples$des=="early.mid.secretory.Glandular")]<- "Y"
Samples$group[which(Samples$des=="early.mid.secretory.Lumenal1")]<- "Y"
Samples$group[which(Samples$des=="early.mid.secretory.Ciliated")]<- "Y"



#extracting data 
Count<- as.data.frame(Count)
dat<-Count[,which(Samples$group=="Y")]
Samples<-Samples[which(Samples$group=="Y"),]

colnames(dat)<-Samples$...1



# Grouping factor

Samples$group[which(Samples$des=="proliferative.Ciliated")]<- "PC"
Samples$group[which(Samples$des=="proliferative.Lumenal1")]<- "PL"
Samples$group[which(Samples$des=="proliferative.Glandular")]<- "PG"
Samples$group[which(Samples$des=="mid.secretory.Glandular")]<- "MSG"
Samples$group[which(Samples$des=="mid.secretory.Lumenal1")]<- "MSL"
Samples$group[which(Samples$des=="mid.secretory.Ciliated")]<- "MSC"
Samples$group[which(Samples$des=="early.mid.secretory.Glandular")]<- "MSG"
Samples$group[which(Samples$des=="early.mid.secretory.Lumenal1")]<- "MSL"
Samples$group[which(Samples$des=="early.mid.secretory.Ciliated")]<- "MSC"


# Seurat

data<-CreateSeuratObject(counts = dat, project = "Single cell", min.cells = 3, min.features = 200)


Group<-Samples$group
names(Group)<-Samples$SampleID
data$group<-Group



orig_ident<-Samples$SampleID
names(orig_ident)<-Samples$X
data$orig.ident<-orig_ident

data$orig.ident<-as.factor(data$orig.ident)

data@active.ident<-data@meta.data$orig.ident

#data$group <- as.factor(data$group)
#data@active.ident <- data@meta.data$group


# Filter out low-quality cells

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),)

data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)



#Normalize data

data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 30000)


all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)

#differential gene expression 
DGE <- FindMarkers(data, group.by= "group", ident.1= "PC", ident.2 = "MSC")
Significant<-DGE[which(DGE$p_val_adj<0.05),]

write.csv(Significant,"Glandular 0.05.csv")

DimHeatmap(data, dims = 1, cells = 2493, balanced = TRUE)

data2<-GetAssayData(data)
data2<-as.data.frame(data2)

#differential gene expression 

DGE <- FindMarkers(data, group.by= "group", ident.1= "PL", ident.2 = "MSL", test.use = "DESeq2")

write.csv (Significant, ".csv")



