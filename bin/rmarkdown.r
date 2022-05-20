
```{r}
install.packages("minfi")
library(minfi)
library(cluster)

matrixn <- read.csv("C:/Users/D/Desktop/studijos/bakis/opera-search/matrixn.csv", row.names = 1)
matrixp <- read.csv("C:/Users/D/Desktop/studijos/bakis/opera-search/matrixp.csv", row.names = 1)
sum <- read.csv("C:/Users/D/Desktop/studijos/bakis/opera-search/inputs/Sum.csv", row.names=1, sep= ";", skip=1)
busco <- read.csv(file="C:/Users/D/Desktop/studijos/bakis/opera-search/busco_values.csv", row.names=1, sep= ";")

```
```{r}

rownames(busco) <- substr(rownames(busco),17,31)
rownames(matrixp) <- substr(rownames(matrixp),9,23)
rownames(matrixn) <- substr(rownames(matrixn),9,23)


goodBuscoP <- busco[rownames(matrixp),]
goodBuscoN <- busco[rownames(matrixn),]

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
plot( clustersn )
plot( clustersp )

```
``` {r}
BiocManager::install("WGCNA", force=TRUE)
library("WGCNA")

cleanSumP <- cbind(cleanSumP,Country=sub(":.*", "", cleanSumP$Geo_loc_name))
countries <- sub(":.*", "", cleanSumP$Geo_loc_name)

colorp <- cbind(labels2colors(cleanSumP$Host),labels2colors(cleanSumP$Isolation_source), labels2colors(cleanSumP$Geo_loc_name), labels2colors(sub(":.*", "", cleanSumP$Geo_loc_name)), numbers2colors(cleanBuscoP$C))

colnames(colorp) <- c("Host", "Isolation_source", "Geo_loc_name", "Country", "Completeness");

plotDendroAndColors(clustersp,colorp,cex.dendroLabels = 0.1, addGuide=TRUE)

colorn <- cbind(labels2colors(cleanSumN$Host),labels2colors(cleanSumN$Isolation_source), labels2colors(cleanSumN$Geo_loc_name), labels2colors(sub(":.*", "", cleanSumN$Geo_loc_name)), numbers2colors(cleanBuscoN$C))
colnames(colorn) <- c("Host", "Isolation_source", "Geo_loc_name", "Country", "Completeness");

plotDendroAndColors(clustersn,colorn, cex.dendroLabels = 0.1, addGuide=TRUE, col.dendroLabels="red")
```
```{r}
library(ggplot2)


MDSp = cmdscale( distancesp )
MDS1p = cmdscale( distancesp, k=1 )
MDSn = cmdscale( distancesn )
MDS1n = cmdscale( distancesn, k=1 )
```
```{r}
plot( MDSp , col=colorp[,1] )
plot(MDS1p , col=colorp[,1])
plot( MDSn , col=colorn[,1] )
plot(MDS1n , col=colorn[,1])

plot( MDSp , col=colorp[,2] )
plot(MDS1p , col=colorp[,2])
plot( MDSn , col=colorn[,2] )
plot(MDS1n , col=colorn[,2])

plot( MDSp , col=colorp[,3] )
plot(MDS1p , col=colorp[,3])
plot( MDSn , col=colorn[,3] )
plot(MDS1n , col=colorn[,3])

plot( MDSp , col=colorp[,4] )
plot(MDS1p , col=colorp[,4])
plot( MDSn , col=colorn[,4] )
plot(MDS1n , col=colorn[,4])
```
```{r}
library(matrixStats)
Nclust1 <- MDSn[MDSn[,1]>0.1,]
length(binMatrixn[,1])/3
Nclust1vfg <- labels(which ((colSums(binMatrixn[rownames(Nclust1),])>50) == TRUE))
Nclust2 <- MDSn[MDSn[,2]>0.2,]
Nclust2vfg <- labels(which ((colSums(binMatrixn[rownames(Nclust2),])>50) == TRUE))
Nclust3 <- MDSn[MDSn[,2]<0.2 & MDSn[,2]>-0.05 & MDSn[,1]< -0.11,]
Nclust3vfg <- labels(which ((colSums(binMatrixn[rownames(Nclust3),])>50) == TRUE))
Nclust4 <- MDSn[MDSn[,2]<0.2 & MDSn[,2]>-0.07 & MDSn[,1]>-0.11 & MDSn[,1]<0.1,]
Nclust4vfg <- labels(which ((colSums(binMatrixn[rownames(Nclust4),])>50) == TRUE))
Nclust5 <- MDSn[MDSn[,2]< -0.05,]
Nclust5vfg <- labels(which ((colSums(binMatrixn[rownames(Nclust5),])>50) == TRUE))
plot(MDSn[MDSn[,2]< -0.05,])

Pclust1 <- MDSp[MDSp[,1]>0.01 & MDSp[,2]>0.01,]
Pclust1vfg <- labels(which ((colSums(binMatrixp[rownames(Pclust1),])>length(Pclust1[,1])/3) == TRUE))
colorp1 <- colorp[MDSp[,1]>0.01 & MDSp[,2]>0.01,]
plot(Pclust1, col=colorp1[,1])

Pclust2 <- MDSp[MDSp[,1]>-0.05 & MDSp[,2]>-0.01 & MDSp[,2]<0.01,]
Pclust2vfg <- labels(which ((colSums(binMatrixp[rownames(Pclust2),])>length(Pclust2[,1])/3) == TRUE))
colorp2 <- colorp[MDSp[,1]>-0.05 & MDSp[,2]>-0.01 & MDSp[,2]<0.01,]
plot(Pclust2, col=colorp2[,1])

Pclust3 <- MDSp[MDSp[,1]>-0.05 & MDSp[,2]< -0.01,]
Pclust3vfg <- labels(which ((colSums(binMatrixp[rownames(Pclust3),])>length(Pclust3[,1])/3) == TRUE))
colorp3 <- colorp[MDSp[,1]>-0.05 & MDSp[,2]< -0.01,]
plot(Pclust3, col=colorp3[,1])

Pclust4 <- MDSp[MDSp[,1]< -0.07 & MDSp[,2]< 0.03,]
Pclust4vfg <- labels(which ((colSums(binMatrixp[rownames(Pclust4),])>length(Pclust4[,1])/3) == TRUE))
colorp4 <- colorp[MDSp[,1]< -0.07 & MDSp[,2]< 0.03,]
plot(Pclust4, col=colorp4[,1])

Pclust5 <- MDSp[MDSp[,1]< -0.07 & MDSp[,2]> 0.03,]
Pclust5vfg <- labels(which ((colSums(binMatrixp[rownames(Pclust5),])>length(Pclust5[,1])/3) == TRUE))
colorp5 <- colorp[MDSp[,1]< -0.07 & MDSp[,2]> 0.03,]
plot(Pclust5, col=colorp5[,1])

```
```{r}
install.packages("ggVennDiagram")
library(ggVennDiagram)

# List of items
Pclusters <- list("Pclust1"  = Pclust1vfg , "Pclust2" = Pclust2vfg , "Pclust3" = Pclust3vfg , "Pclust4" = Pclust4vfg , "Pclust5" = Pclust5vfg )
Nclusters <- list("Nclust1"  = Nclust1vfg , "Nclust2" = Nclust2vfg , "Nclust3" = Nclust3vfg , "Nclust4" = Nclust4vfg , "Nclust5" = Nclust5vfg )

ggVennDiagram(Pclusters)
ggVennDiagram(Nclusters)

```
cleanSumP[cleanSumP[,"Country"]==",,,","Country"] <- "missing,,,"
cleanSumP[cleanSumP[,"Isolation_source"]=="","Isolation_source"] <- "missing"
sort(prop.table(table(cleanSumP[rownames(Pclust5),"Country"])))
sort(table(cleanSumP[rownames(Pclust1),"Country"]))
tab <- table(cleanSumP[rownames(Pclust1),"Isolation_source"])
sort(tab)
sort(prop.table(tab))

1. išsivalyti duom
2. Aprašymą clusters pagal virulent gen.
3. Aprašyti clusters pagal metadata