
library(edgeR) #https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
library(biomaRt)
library(mixOmics)
library(dplyr)

setwd("~/Desktop/RP 2/Endometrium analysis/Chi")
samples<-read.delim("Sample.txt")
dat<-read.delim("Count.txt")
rownames(dat)<-dat[,1]
dat<-dat[,-c(1)]


genes<-rownames(dat) # select list of gene IDs
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),
                  values = genes, mart= ensembl)
dat<-dat[na.omit(match(gene_IDs$ensembl_gene_id,rownames(dat))),]


DGE<-DGEList(counts=dat,group=samples$cycle) 

keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] #Remove low expression genes

DGE <- calcNormFactors(DGE)

DGE <- estimateDisp(DGE) #Estimate dispersions
et <- exactTest(DGE) #Calculate DE
et$table$Gene_name<-gene_IDs$external_gene_name[na.omit(match(rownames(et$table),gene_IDs$ensembl_gene_id))]
significant<-topTags(et,n = 9346,p.value = 0.01)
Signfigenes<-significant$table$Gene_name

cpm_data<-cpm(DGE,log = T) #Extract normalised log2 CPM values
rownames(cpm_data)<-gene_IDs$external_gene_name[na.omit(match(rownames(et$table),gene_IDs$ensembl_gene_id))]

pd<-pca(t(cpm_data))

# png("LH_unsupervised_transcriptome_pca.png")
plotIndiv(pd,group=samples$cycle,ind.names = samples$ID,ellipse = T)
# dev.off()

cpm_data<- cpm_data[na.omit(match(Signfigenes,rownames(cpm_data))),]
pd<-pca(t(cpm_data[na.omit(match(rownames(cpm_data),topTags(et,n = 440)$table$Gene_name)),]))


plotIndiv(pd,group=samples$cycle,ind.names = samples$ID, ellipse = T)



print(genelist)


#genes <-dat$Geneid
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),
#values = genes, mart= mart)
#colnames(gene_IDs)<- c("ensembl_gene_id","Geneid")

#dat<- (left_join(dat, gene_IDs, by = c("Geneid"="ensembl_gene_id")))
#at<- dat[!(dat$Geneid.y == ""), ]
#rownames(dat)<-dat[,11]







