

install.packages("gplots")
install.packages("Boruta")
install.packages("DMwR")

library(gplots) # heatmap.2 function required
library(Boruta) # Boruta function required
library(DMwR) # SMOTE function required
library(edgeR)

#Data input - source is transcriptomic data



PVMSCiliated<-dat[,which(Samples$group=="PC"|Samples$group=="MSC")]

PVMSGlandular<-dat[,which(Samples$group=="PG"|Samples$group=="MSG")]

PVMSLuminal<-dat[,which(Samples$group=="PL"|Samples$group=="MSL")]

TMM<-cpm(PVMSLuminal, normalized.lib.sizes=TRUE)

TMM<- TMM[apply(TMM[,-1], 1, function(x) !all(x==0)),]


apply(TMM, 1, sd)

gene_list1<-read.delim("Ciliated day 6.txt")

gene_list1<-read.delim("Luminal day7.txt")

gene_list5<-gene_list1$Genes


#correlation matrix generation

cor_matrix<-cor(t(TMM[na.omit(match(gene_list5,rownames(TMM))),]), #correlate genes which match the list against all other genes (omit NA matches for instances where genes not found in matrix)
                t(TMM[-na.omit(match(gene_list5,rownames(TMM))),]),) #NB the '-' in the square brackets on this line but absent from the previous line. This means exclude these entries (entries which match the gene list)

hist(cor_matrix) # examine correlation matrix values for normality

# You need to generate a new object with the correlation matrix values converted to their absolute equivalents (i.e. remove any minus signs from negative values)
binary<-abs(cor_matrix)
sd_thresh<-sd(cor_matrix) # Next calculate sd of the original correlation matrix
binary[which(binary>sd_thresh)]<-1 #anything greater than the sd is converted to a 1
binary[which(binary!=1)]<-0 #anything that is not now a 1 becomes a zero


binary<-abs(cor_matrix) # generate absolute values for all correlation r-values. This matrix will be binarized
thresh<-0.3 # set r-value threshold for binarization
binary[which(binary>thresh)]<-1 # set any values greater than the threshold to 1
binary[which(binary!=1)]<-0 # set any values which aren't 1s as 0s ('!=' means does not equal)

#hypernetwork generation and visualisation
hyp<-binary%*%t(binary) # generate the hypernetwork by multiplying the binary matrix by the transpose of itself. ('%*%' is the operator for matrix multiplication )

hm<-heatmap.2(hyp,trace="none") 
 # generate a heatmap from the hypernetwork matrix, save the heatmap object as it contains dendrograms.Exclude the trace which aims to separates rows and columns

dendrogram<-as.hclust(hm$rowDendrogram) 

ct<-cutree(dendrogram,k = 4) # Cut dendrogram to generate 2 groups (first dendrogram split)
# central_cluster_genes<-names(ct[which(ct==1)]) # extract names of the 1st cluster (check cluster 1 is appropriate) 
table(ct)

#list for central cluster
cc1<-names(ct[which(ct==1)])
cc2<-names(ct[which(ct==4)])

galois1<-binary[match(cc1,rownames(binary)),]
galois1<-galois1[,which(colSums(galois1)==nrow(galois1))]
galois1<-galois1[,which(colSums(galois1)>nrow(galois1)*0.8)]

galois1<-binary[match(cc2,rownames(binary)),]
galois1<-galois1[,which(colSums(galois1)==nrow(galois1))]
galois1<-galois1[,which(colSums(galois1)>nrow(galois1)*0.8)]
galois1<-t(galois1)


write.csv (galois1, " .csv")



