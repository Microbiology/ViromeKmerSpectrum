library(reshape2)
library(vegan)
library(ggplot2)
library(cluster)
library(ape)

INPUT <- read.delim("/Users/Hannigan/git/ViromeKmerSpectrum/data/kmerSpectrum.tsv", sep="\t", header=F)

CINPUT <- dcast(INPUT, V1~V3, value.var="V2")
rownames(CINPUT) <- CINPUT$V1
CINPUT <- (CINPUT[,-1])
CINPUT[is.na(CINPUT)] <- 0
TCINPUT <- t(CINPUT)
colnames(TCINPUT) <- gsub("\\;","_", colnames(TCINPUT))

DistMatrix <- vegdist(TCINPUT, method = "bray")
OrdMds <- metaMDS(DistMatrix)
OrdMdsFit = data.frame(MDS1 = OrdMds$points[,1], MDS2 = OrdMds$points[,2])
#NMDS_AND_MAP <- cbind(BRAY_ORD_FIT,INPUT_MAP)
ggplot(OrdMdsFit, aes(x=MDS1, y=MDS2)) + theme_classic() + geom_point(size=4) + scale_colour_brewer(palette="Set2")

wssplot <- function(data, nc, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}
KMEAN_DIST<-as.matrix(DistMatrix)
wssplot(as.matrix(KMEAN_DIST),nc=25)

set.seed(1234)
fit<-kmeans(KMEAN_DIST, 8)
clusplot(KMEAN_DIST, fit$cluster, shade=TRUE, color=TRUE)

hvar <- hclust(DistMatrix)
pdf(file="/Users/Hannigan/git/ViromeKmerSpectrum/data/kmerSpectrumTree.pdf", height=30, width=8)
plot(as.phylo(hvar), cex = 0.2)
dev.off()



