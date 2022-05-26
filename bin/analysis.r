setwd("./Ranalysis")

#install.packages("minfi")
#library(minfi)
library(cluster)
#BiocManager::install("WGCNA", lib=c("./lib/R"))
library(WGCNA)
install.packages("ggVennDiagram", lib=c("./lib/R"))
library(ggVennDiagram, lib=c("./lib/R"))

matrixn <- read.csv("inputs/matrixn.csv", row.names = 1)
matrixp <- read.csv("inputs/matrixp.csv", row.names = 1)
sum <- read.csv("inputs/Metadata.csv", row.names=1, sep= "\t", skip=1)
busco <- read.csv("inputs/busco_values.csv", row.names=1, sep= ",")
vfMap <- read.csv("inputs/VFmap.csv", row.names=1)
sum$Host <- tolower(sum$Host)
sum$Isolation_source <- tolower(sum$Isolation_source)
sum$Geo_loc_name <- tolower(sum$Geo_loc_name)


rownames(busco) <- substr(rownames(busco),17,31)
rownames(matrixp) <- substr(rownames(matrixp),9,23)
rownames(matrixn) <- substr(rownames(matrixn),9,23)


goodBuscoP <- busco[rownames(matrixp),]
goodBuscoN <- busco[rownames(matrixn),]

setEPS()
postscript("outputs/Completeness_hist.eps")
hist(busco$C, col="light blue", xlab="Pilnumo įvertis (%)", ylab="Dažnis", main="Busco genomų pilnumo įvertis")
dev.off()

cleanBuscoP <- goodBuscoP[goodBuscoP[,1]>90,]
cleanBuscoN <- goodBuscoN[goodBuscoN[,1]>90,]
cleanMatrixp <- matrixp[rownames(cleanBuscoP),]
cleanMatrixn <- matrixn[rownames(cleanBuscoN),]
cleanSumP <- sum[rownames(cleanBuscoP),]
cleanSumN <- sum[rownames(cleanBuscoN),]



binMatrixp <-cleanMatrixp;
binMatrixp[binMatrixp>0] <- 1;
binMatrixn <-cleanMatrixn;
binMatrixn[binMatrixn>0] <- 1;

# Skirtumų matrica, 'binary' metodas yra Tanimoto indeksas iš esmės
distancesp = dist( binMatrixp, method = 'binary' )
distancesn = dist( binMatrixn, method = 'binary' )



# Klasterizavimas bei pavaizdavimas
clustersp = hclust( distancesp )
clustersn = hclust( distancesn )
Pkmeans <- kmeans(distancesp,4)
PkmeansColors <- Pkmeans$cluster
PkmeansColors[PkmeansColors==1] <- "turquoise"
PkmeansColors[PkmeansColors==2] <- "blue"
PkmeansColors[PkmeansColors==3] <- "brown"
PkmeansColors[PkmeansColors==4] <- "yellow"

Nkmeans <- kmeans(distancesn,4)
NkmeansColors <- Nkmeans$cluster
NkmeansColors[NkmeansColors==1] <- "turquoise"
NkmeansColors[NkmeansColors==2] <- "blue"
NkmeansColors[NkmeansColors==3] <- "brown"
NkmeansColors[NkmeansColors==4] <- "yellow"


cleanSumP <- cbind(cleanSumP,Country=sub(":.*", "", cleanSumP$Geo_loc_name))
cleanSumN <- cbind(cleanSumN,Country=sub(":.*", "", cleanSumN$Geo_loc_name))
#countries <- sub(":.*", "", cleanSumP$Geo_loc_name)

colorp <- cbind(PkmeansColors, labels2colors(cleanSumP$Host),labels2colors(cleanSumP$Isolation_source), labels2colors(cleanSumP$Country), numbers2colors(cleanBuscoP$C))

colnames(colorp) <- c("Klasteris", "Nesiotojas", "Izoliatas", "Salis", "Genomo pilnumas");

setEPS()
postscript("outputs/pro_dendro.eps")
plotDendroAndColors(clustersp,colorp,cex.dendroLabels = 0.1, addGuide=TRUE, ylab="Atstumas", main="Klasterių dendrograma")
dev.off()

colorn <- cbind(labels2colors(Nkmeans$cluster), labels2colors(cleanSumN$Host),labels2colors(cleanSumN$Isolation_source), labels2colors(cleanSumN$Country), numbers2colors(cleanBuscoN$C))
colnames(colorn) <- c("Klasteris", "Nesiotojas", "Izoliatas", "Salis", "Genomo pilnumas");

setEPS()
postscript("outputs/nuc_dendro.eps")
plotDendroAndColors(clustersn,colorn,cex.dendroLabels = 0.1, addGuide=TRUE, ylab="Atstumas", main="Klasterių dendrograma")
dev.off()

#library(ggplot2)

MDSp = cmdscale( distancesp )
MDSn = cmdscale( distancesn )

setEPS()
postscript("outputs/pro_MDS_clusters.eps")
plot( MDSp , col=colorp[,1] , main="Multidimensinis masteliavimas. Klasteriai")
dev.off()
setEPS()
postscript("outputs/nuc_MDS_clusters.eps")
plot( MDSn , col=colorn[,1] , main="Multidimensinis masteliavimas. Klasteriai")
dev.off()

setEPS()
postscript("outputs/pro_MDS_hosts.eps")
plot( MDSp , col=colorp[,2] , main="Multidimensinis masteliavimas. Nešiotojai")
dev.off()
setEPS()
postscript("outputs/nuc_MDS_hosts.eps")
plot( MDSn , col=colorn[,2] , main="Multidimensinis masteliavimas. Nešiotojai")
dev.off()

setEPS()
postscript("outputs/pro_MDS_isolation_source.eps")
plot( MDSp , col=colorp[,3] , main="Multidimensinis masteliavimas. Izoliatai")
dev.off()
setEPS()
postscript("outputs/nuc_MDS_isolation_source.eps")
plot( MDSn , col=colorn[,3] , main="Multidimensinis masteliavimas. Izoliatai")
dev.off()

Nclust1vfg <- substr(labels(which ((colSums(binMatrixn[Nkmeans$cluster==1,])>length(which(Nkmeans$cluster==1))/3) == TRUE)),1,9)
Nclust2vfg <- substr(labels(which ((colSums(binMatrixn[Nkmeans$cluster==2,])>length(which(Nkmeans$cluster==2))/3) == TRUE)),1,9)
Nclust3vfg <- substr(labels(which ((colSums(binMatrixn[Nkmeans$cluster==3,])>length(which(Nkmeans$cluster==3))/3) == TRUE)),1,9)
Nclust4vfg <- substr(labels(which ((colSums(binMatrixn[Nkmeans$cluster==4,])>length(which(Nkmeans$cluster==4))/3) == TRUE)),1,9)

Pclust1vfg <- substr(labels(which ((colSums(binMatrixp[Pkmeans$cluster==1,])>length(which(Pkmeans$cluster==1))/3) == TRUE)),1,9)
Pclust2vfg <- substr(labels(which ((colSums(binMatrixp[Pkmeans$cluster==2,])>length(which (Pkmeans$cluster==2))/3) == TRUE)),1,9)
Pclust3vfg <- substr(labels(which ((colSums(binMatrixp[Pkmeans$cluster==3,])>length(which (Pkmeans$cluster==3))/3) == TRUE)),1,9)
Pclust4vfg <- substr(labels(which ((colSums(binMatrixp[Pkmeans$cluster==4,])>length(which (Pkmeans$cluster==4))/3) == TRUE)),1,9)


Pclust1vfg <- cbind(Pclust1vfg, vfMap[Pclust1vfg,1:2])
Pclust2vfg <- cbind(Pclust2vfg, vfMap[Pclust2vfg,1:2])
Pclust3vfg <- cbind(Pclust3vfg, vfMap[Pclust3vfg,1:2])
Pclust4vfg <- cbind(Pclust4vfg, vfMap[Pclust4vfg,1:2])

PclustersVFG <- list("Turkio"  = Pclust1vfg[,1] , "Melynas" = Pclust2vfg[,1] , "Rudas" = Pclust3vfg[,1] , "Geltonas" = Pclust4vfg[,1] )
PclustersVF <- list("Turkio"  = Pclust1vfg[,2] , "Melynas" = Pclust2vfg[,2] , "Rudas" = Pclust3vfg[,2] , "Geltonas" = Pclust4vfg[,2] )
PclustersVFC <- list("Turkio"  = Pclust1vfg[,3] , "Melynas" = Pclust2vfg[,3] , "Rudas" = Pclust3vfg[,3] , "Geltonas" = Pclust4vfg[,3] )

Nclust1vfg <- cbind(Nclust1vfg, vfMap[Nclust1vfg,1:2])
Nclust2vfg <- cbind(Nclust2vfg, vfMap[Nclust2vfg,1:2])
Nclust3vfg <- cbind(Nclust3vfg, vfMap[Nclust3vfg,1:2])
Nclust4vfg <- cbind(Nclust4vfg, vfMap[Nclust4vfg,1:2])
NclustersVFG <- list("Turkio"  = Nclust1vfg[,1] , "Melynas" = Nclust2vfg[,1] , "Rudas" = Nclust3vfg[,1] , "Geltonas" = Nclust4vfg[,1]  )
NclustersVF <- list("Turkio"  = Nclust1vfg[,2] , "Melynas" = Nclust2vfg[,2] , "Rudas" = Nclust3vfg[,2] , "Geltonas" = Nclust4vfg[,2]  )
NclustersVFC <- list("Turkio"  = Nclust1vfg[,3] , "Melynas" = Nclust2vfg[,3] , "Rudas" = Nclust3vfg[,3] , "Geltonas" = Nclust4vfg[,3]  )


setEPS()
postscript("outputs/pro_Venn_VFG.eps")
ggVennDiagram(PclustersVFG, label_color = "purple")
dev.off()
setEPS()
postscript("outputs/pro_Venn_VF.eps")
ggVennDiagram(PclustersVF, label_color = "purple")
dev.off()
setEPS()
postscript("outputs/pro_Venn_VFC.eps")
ggVennDiagram(PclustersVFC, label_color = "purple")
dev.off()

setEPS()
postscript("outputs/nuc_Venn_VFG.eps")
ggVennDiagram(NclustersVFG, label_color = "purple")
dev.off()
setEPS()
postscript("outputs/nuc_Venn_VF.eps")
ggVennDiagram(NclustersVF, label_color = "purple")
dev.off()
setEPS()
postscript("outputs/nuc_Venn_VFC.eps")
ggVennDiagram(NclustersVFC, label_color = "purple")
dev.off()


Turquoise<-sort(table(cleanSumP[Pkmeans$cluster==1,"Country"]))
Blue<-sort(table(cleanSumP[Pkmeans$cluster==2,"Country"]))
Brown<-sort(table(cleanSumP[Pkmeans$cluster==3,"Country"]))
Yellow<-sort(table(cleanSumP[Pkmeans$cluster==4,"Country"]))
write.table(Turquoise, "outputs/pro_turquoise_cluster_countries.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Blue, "outputs/pro_blue_cluster_countries.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Brown, "outputs/pro_brown_cluster_countries.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Yellow, "outputs/pro_yellow_cluster_countries.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))

Turquoise<-sort(table(cleanSumP[Pkmeans$cluster==1,"Host"]))
Blue<-sort(table(cleanSumP[Pkmeans$cluster==2,"Host"]))
Brown<-sort(table(cleanSumP[Pkmeans$cluster==3,"Host"]))
Yellow<-sort(table(cleanSumP[Pkmeans$cluster==4,"Host"]))
write.table(Turquoise, "outputs/pro_turquoise_cluster_hosts.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Blue, "outputs/pro_blue_cluster_hosts.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Brown, "outputs/pro_brown_cluster_hosts.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Yellow, "outputs/pro_yellow_cluster_hosts.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))

Turquoise<-sort(table(cleanSumP[Pkmeans$cluster==1,"Isolation_source"]))
Blue<-sort(table(cleanSumP[Pkmeans$cluster==2,"Isolation_source"]))
Brown<-sort(table(cleanSumP[Pkmeans$cluster==3,"Isolation_source"]))
Yellow<-sort(table(cleanSumP[Pkmeans$cluster==4,"Isolation_source"]))
write.table(Turquoise, "outputs/pro_turquoise_cluster_isolation_sources.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Blue, "outputs/pro_blue_cluster_isolation_sources.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Brown, "outputs/pro_brown_cluster_isolation_sources.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Yellow, "outputs/pro_yellow_cluster_isolation_sources.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))



Turquoise<-sort(table(cleanSumN[Nkmeans$cluster==1,"Country"]))
Blue<-sort(table(cleanSumN[Nkmeans$cluster==2,"Country"]))
Brown<-sort(table(cleanSumN[Nkmeans$cluster==3,"Country"]))
Yellow<-sort(table(cleanSumN[Nkmeans$cluster==4,"Country"]))
write.table(Turquoise, "outputs/nuc_turquoise_cluster_countries.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Blue, "outputs/nuc_blue_cluster_countries.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Brown, "outputs/nuc_brown_cluster_countries.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Yellow, "outputs/nuc_yellow_cluster_countries.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))

Turquoise<-sort(table(cleanSumN[Nkmeans$cluster==1,"Host"]))
Blue<-sort(table(cleanSumN[Nkmeans$cluster==2,"Host"]))
Brown<-sort(table(cleanSumN[Nkmeans$cluster==3,"Host"]))
Yellow<-sort(table(cleanSumN[Nkmeans$cluster==4,"Host"]))
write.table(Turquoise, "outputs/nuc_turquoise_cluster_hosts.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Blue, "outputs/nuc_blue_cluster_hosts.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Brown, "outputs/nuc_brown_cluster_hosts.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Yellow, "outputs/nuc_yellow_cluster_hosts.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))

Turquoise<-sort(table(cleanSumN[Nkmeans$cluster==1,"Isolation_source"]))
Blue<-sort(table(cleanSumN[Nkmeans$cluster==2,"Isolation_source"]))
Brown<-sort(table(cleanSumN[Nkmeans$cluster==3,"Isolation_source"]))
Yellow<-sort(table(cleanSumN[Nkmeans$cluster==4,"Isolation_source"]))
write.table(Turquoise, "outputs/nuc_turquoise_cluster_isolation_sources.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Blue, "outputs/nuc_blue_cluster_isolation_sources.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Brown, "outputs/nuc_brown_cluster_isolation_sources.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))
write.table(Yellow, "outputs/nuc_yellow_cluster_isolation_sources.csv", append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Country","Count"))

