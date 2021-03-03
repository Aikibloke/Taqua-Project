# The object 'data' is a matrix genes x samples
# The object 'gene_list' is a list of genes (eg DEGs) which form the nodes of the hypernetwork
# The object 'rf_variable' is a variable separating the samples used for random forest


install.packages("gplots")
install.packages("Boruta")
install.packages("DMwR")

library(gplots) # heatmap.2 function required
library(Boruta) # Boruta function required
library(DMwR) # SMOTE function required

#Data input - source is transcriptomic data
data<-readRDS("List.rds")

Expression<-data[[1]]

plot(Expression)

colnames(Expression)<-data[[3]]$`Pubertal Group`

gene_list1<-read.delim("List.txt")

full_list<-data[[2]]$`Gene Symbol`

gene_list<-gene_list1$Gene.Symbol

# principal component analysis
BiocManager::install("mixOmics")

library(ggplot2)
library(mixOmics)
Y<- t(Expression)
Y<-as.data.frame(Y)
Y<-as.data.frame(cbind(data[[3]]$`Pubertal Group`,Y))
colnames(Y)[1] <- "Pubertal Group"
Y$`Pubertal Group`<-as.factor(Y$`Pubertal Group`)
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
plotIndiv(pd2,ind.names = F,group = clusters, ellipse = T, star=F,title="pubertal compared to non-pubertal")

#correlation matrix generation
cor_matrix<-cor(t(Expression[na.omit(match(gene_list,rownames(Expression))),]), #correlate genes which match the list against all other genes (omit NA matches for instances where genes not found in matrix)
                t(Expression[-na.omit(match(gene_list,rownames(Expression))),])) #NB the '-' in the square brackets on this line but absent from the previous line. This means exclude these entries (entries which match the gene list)

hist(cor_matrix) # examine correlation matrix values for normality

binary<-abs(cor_matrix) # generate absolute values for all correlation r-values. This matrix will be binarized
thresh<-0.3 # set r-value threshold for binarization
binary[which(binary>thresh)]<-1 # set any values greater than the threshold to 1
binary[which(binary!=1)]<-0 # set any values which aren't 1s as 0s ('!=' means does not equal)

#hypernetwork generation and visualisation
hyp<-binary%*%t(binary) # generate the hypernetwork by multiplying the binary matrix by the transpose of itself. ('%*%' is the operator for matrix multiplication )

hm<-heatmap.2(hyp,trace="none") # generate a heatmap from the hypernetwork matrix, save the heatmap object as it contains dendrograms.Exclude the trace which aims to separates rows and columns

dendrogram<-as.hclust(hm$rowDendrogram) 

ct<-cutree(dendrogram,k = 3) # Cut dendrogram to generate 2 groups (first dendrogram split)
# central_cluster_genes<-names(ct[which(ct==1)]) # extract names of the 1st cluster (check cluster 1 is appropriate) 

#list for central cluster
central_cluster_genes<-names(ct[which(ct!=1)])
print(central_cluster_genes)
write.csv (central_cluster_genes, "CCG.csv")

#list for differentially expressed genes not in central cluster
Not_central_cluster_genes<-names(ct[which(ct==1)])

#random forest data foramt for genes in central cluster
rf_data<-t(Expression[na.omit(match(central_cluster_genes,row.names(Expression))),]) # refine original dataset to just genes from the hypernetwork central cluster
rf_variable<-as.factor(data[[3]]$`Pubertal Group`)
rf_data<-as.data.frame(cbind(rf_variable,rf_data))
rf_data$rf_variable<-as.factor(rf_data$rf_variable)


#sample genes to generate a comparison for analysis from the set not in the central cluster.
sampled_genes<-sample(Not_central_cluster_genes,size = 212,replace = FALSE)


#sample genes to generate a comparison for analysis from full list.
sampled_genes1<-sample(full_list,size = 212,replace = FALSE)

#random forest data foramt for genes NOT in central cluster
rf_data1<-t(Expression[na.omit(match(sampled_genes,row.names(Expression))),]) # refine original dataset to just genes from the hypernetwork central cluster
rf_variable<-as.factor(data[[3]]$`Pubertal Group`)
rf_data1<-as.data.frame(cbind(rf_variable,rf_data1))
rf_data1$rf_variable<-as.factor(rf_data1$rf_variable)

#random forest data foramt for genes NOT in central cluster from full list
rf_data2<-t(Expression[na.omit(match(sampled_genes1,row.names(Expression))),]) # refine original dataset to just genes from the hypernetwork central cluster
rf_variable<-as.factor(data[[3]]$`Pubertal Group`)
rf_data2<-as.data.frame(cbind(rf_variable,rf_data2))
rf_data2$rf_variable<-as.factor(rf_data2$rf_variable)

#SMOTE 
#for hypernetwork
rfSM<- SMOTE(rf_variable~., rf_data, perc.over = 200, k = 5, perc.under=200)

#for genes not in hypernetwork central cluster
rfSM1<- SMOTE(rf_variable~., rf_data1, perc.over = 200, k = 5, perc.under=200)

#for genes not in random sample
rfSM2<-SMOTE(rf_variable~., rf_data2, perc.over = 200, k = 5, perc.under=200)

#BORUTA implementation of random forest
#for hypernetwork smote data
Bor<-Boruta(formula=rf_variable~.,data=rfSM, doTrace=2)
stats1<-attStats(Bor)
print(stats1)
plot(Bor)
write.csv (stats1, "HNBOR.csv")

#for genes not in central cluster
Bor1<-Boruta(formula=rf_variable~.,data=rfSM1, doTrace=2)
stats2<-attStats(Bor1)
print(stats2)
plot(Bor1)
write.csv (stats2, "NotHNBOR.csv")

#for genes sampled from full set
Bor2<-Boruta(formula=rf_variable~.,data=rfSM2, doTrace=2)
stats3<-attStats(Bor2)
print(stats3)
plot(Bor2)
write.csv (stats3, "RanBOR.csv")

#load rattle GUI
library(rattle)

rattle()


