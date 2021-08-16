# The object 'data' is a matrix genes x samples
# The object 'gene_list' is a list of genes (eg DEGs) which form the nodes of the hypernetwork
# The object 'rf_variable' is a variable separating the samples used for random forest


install.packages("gplots")
install.packages("Boruta")
install.packages("DMwR")
install.packages("data.table")

library(gplots) # heatmap.2 function required
library(Boruta) # Boruta function required
library(DMwR) # SMOTE function required
library(dplyr)
library(entropy)

setwd("~/Desktop/RP 2/Endometrium analysis/Chi")

cellsurface<-read.delim("cell surface .txt")


cellsurface<- merge(cellsurface,gene_IDs, by.x ="external_gene_name")
genelist<- cellsurface$ensembl_gene_id

samples<-read.delim("Sample.txt")
colnames(TMM)<-samples$cycle

# principal component analysis
BiocManager::install("mixOmics")

library(ggplot2)
library(mixOmics)
Y<- t(data)
Y<-as.data.frame(Y)
Y<-as.data.frame(cbind(samples$cycle, Y))
colnames(Y)[1]<- "Cycle"
Y$`Cycle`<-as.factor(Y$`Cycle`)
clusters<-as.factor(Y[,1])
data1<-Y[,-1]
data1<-as.data.frame(apply(data1,2,as.numeric))#make values numeric

#PCA method 1 - to assess general structure
pca_data<-prcomp(data1,scale. = F,center = T,retx = T)
pd<-as.data.frame(pca_data$x)
pd$clusters<-clusters
ggplot(pd,aes(x=PC1,y=PC2,group=clusters,col=clusters))+
  geom_point()

#PCA method 2 - to assess distribution of structure
pd2<-pca(data1,ncomp = 2)
plotIndiv(pd2,ind.names = F,group = clusters, ellipse = T, star=F,title="Secretory compared to proliferative")

#correlation matrix generation
cor_matrix<-cor(t(TMM[na.omit(match(genelist,rownames(TMM))),]), #correlate genes which match the list against all other genes omit NA matches for instances where genes not found in matrix)
                t(TMM[-na.omit(match(genelist,rownames(TMM))),]))  #NB the '-' in the square brackets on this line but absent from the previous line. This means exclude these entries (entries which match the gene list)

hist(cor_matrix) # examine correlation matrix values for normality

sd(cor_matrix)


binary<-abs(cor_matrix) # generate absolute values for all correlation r-values. This matrix will be binarized
thresh<-0.8 # set r-value threshold for binarization
binary[which(binary>thresh)]<-1 # set any values greater than the threshold to 1
binary[which(binary!=1)]<-0 # set any values which aren't 1s as 0s ('!=' means does not equal)

binary<-as.data.frame(binary)
binary<- tibble::rownames_to_column(binary, "ensembl_gene_id")
binary<- merge(binary,gene_IDs, by.x = "ensembl_gene_id", by.y ="ensembl_gene_id")
rownames(binary)<- binary[,9233]
binary<-binary[,-c(1,9233)]
binary<-as.matrix(binary)


#hypernetwork generation and visualisation
hyp<-binary%*%t(binary) # generate the hypernetwork by multiplying the binary matrix by the transpose of itself. ('%*%' is the operator for matrix multiplication )


hm<-heatmap.2(hyp,trace="none") # generate a heatmap from the hypernetwork matrix, save the heatmap object as it contains dendrograms.Exclude the trace which aims to separates rows and columns

dend<-as.hclust(hm$rowDendrogram)
ct<-cutree(dend,k=5) # Cut dendrogram to generate 2 groups (first dendrogram split)
table(ct)
cc1<-names(ct[which(ct==1)])
cc2<-names(ct[which(ct==2)])
cc3<-names(ct[which(ct==4|ct==5)]) # central_cluster_genes<-names(ct[which(ct==1)]) # extract names of the 1st cluster

list(cc1)

galoise1<-binary[match(cc1,rownames(binary)),]
galoise1<-galoise1[,which(colSums(galoise1)==nrow(galoise1))]
galoise1<-galoise1[,which(colSums(galoise1)>nrow(galoise1)*0.8)]

galoise2<-binary[match(cc2,rownames(binary)),]
galoise2<-galoise2[,which(colSums(galoise2)==nrow(galoise2))]
galoise2<-galoise2[,which(colSums(galoise2)>nrow(galoise2)*0.8)]

galoise3<-binary[match(cc3,rownames(binary)),]
galoise3<-galoise3[,which(colSums(galoise3)==nrow(galoise3))]
galoise3<-galoise3[,which(colSums(galoise3)>nrow(galoise3)*0.8)]

#list for central cluster
print(galoise1)
write.csv (galoise1, "galoise1.csv")
write.csv (galoise2, "galoise2.csv")
write.csv (galoise3, "galoise3.csv")



