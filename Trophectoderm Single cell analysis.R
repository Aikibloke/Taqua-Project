

library(Seurat)

Count<- t(Count)
rownames(Count)<- var$...1


day6<-Count[,which(Samples$Age=="6")]
day6<-as.data.frame(day6)
day7<-Count[,which(Samples$Age=="7")]
day7<-as.data.frame(day7)


sday6<-CreateSeuratObject(counts = day6, project = "Single cell", min.cells = 3, min.features = 200)

sday7<-CreateSeuratObject(counts = day7, project = "Single cell", min.cells = 3, min.features = 200)


sday6[["percent.mt"]] <- PercentageFeatureSet(sday6, pattern = "^MT-")

VlnPlot(sday6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),)

sday6 <- subset(sday6, subset = nFeature_RNA > 7500 )

sday7[["percent.mt"]] <- PercentageFeatureSet(sday7, pattern = "^MT-")

VlnPlot(sday7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),)

sday7 <- subset(sday7, subset = nFeature_RNA < 14000 )

day6file<- sday6@assays$RNA@counts
day6file<- as.data.frame(day6file)

day7file<- sday7@assays$RNA@counts
day7file<- as.data.frame(day7file)

day6<-row.names(day6file)
day7<-row.names(day7file)


write.csv(day6,"day6 genelist")


write.csv(day7,"day7 genelist")






