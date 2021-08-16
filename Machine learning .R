Machine learning code

library(DMwR) # SMOTE function required
library(Boruta) # Boruta function required

Y<- read.table("RIF1B.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Y$group = as.factor(Y$group)


#SMOTE
DDS<- SMOTE(group ~ ., data = Y, perc.over = 200, k = 5, perc.under=200)
#BORUTA
library(mlbench)
Boruta(group~.,data=DDS,doTrace=2)->Bor1.hvo;

print(Bor1.hvo);
plot(Bor1.hvo);
stats1<-attStats(Bor1.hvo);
print(stats1);
write.csv (stats1, "GHD_data.csv")

rattle()
install.packages("RGtk2", depen=T)

