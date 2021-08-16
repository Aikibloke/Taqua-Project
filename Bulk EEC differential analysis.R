library(edgeR) #https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
library(biomaRt)
library(mixOmics)
library(dplyr)


samples<-read.delim("Sample.txt")
data<-read.delim("Count.txt")
rownames(data)<-data[,1]
data<-data[,-c(1)]


genes<-rownames(data) # select list of gene IDs
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),
                  values = genes, mart= ensembl)
data<-data[na.omit(match(gene_IDs$ensembl_gene_id,rownames(data))),]


DGE<-DGEList(counts=data,group=samples$cycle) 

keep <- filterByExpr(DGE) #ID low expression genes
DGE <- DGE[keep, , keep.lib.sizes=FALSE] #Remove low expression genes

DGE <- calcNormFactors(DGE)

DGE <- estimateDisp(DGE) #Estimate dispersions
et <- exactTest(DGE) #Calculate DE
et$table$Gene_name<-gene_IDs$external_gene_name[na.omit(match(rownames(et$table),gene_IDs$ensembl_gene_id))]
significant<-topTags(et,n = 9346,p.value = 0.01)
Signfigenes<-significant$table$Gene_name

#Convert into TMM file 
TMM<-cpm(DGE, normalized.lib.sizes=TRUE)
rownames(TMM)<-gene_IDs$external_gene_name[na.omit(match(rownames(et$table),gene_IDs$ensembl_gene_id))]



cpm_data<-cpm(DGE,log = T) #Extract normalised log2 CPM values
rownames(cpm_data)<-gene_IDs$external_gene_name[na.omit(match(rownames(et$table),gene_IDs$ensembl_gene_id))]



plotIndiv(pd,group=samples$cycle,ind.names = samples$ID, ellipse = T)

write.csv(significant$table,"~/Desktop/RP 2/Endometrium analysis/Chi\\genelist.csv")

write.csv (TMM,"TMM.csv")


P<-ggplot(data=et$table, aes(x=logFC, y=-log10(PValue))) + geom_point() + theme_minimal()

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")


et$table$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.01, set as "UP" 
et$table$diffexpressed[et$table$logFC > 0.6 & et$table$PValue < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.01, set as "DOWN"
et$table$diffexpressed[et$table$logFC < -0.6 & et$table$PValue < 0.05] <- "DOWN"

p <- ggplot(data=et$table, aes(x=et$table$logFC, y=-log10(PValue), col=diffexpressed)) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")

p3 <- p2 + scale_color_manual(values=c("green", "black", "blue"))


mycolors <- c("green", "black", "blue")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)


ggplot(data=et$table, aes(x=logFC, y=-log10(PValue), col=diffexpressed,)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("green", "black", "blue")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red")

library(plyr)

count(significant$table, "diffexpressed")









