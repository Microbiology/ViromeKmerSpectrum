library(reshape2)
library(vegan)
library(ggplot2)

INPUT <- read.delim("/Users/Hannigan/git/ViromeKmerSpectrum/data/kmerSpectrum.tsv", sep="\t", header=F)

CINPUT <- dcast(INPUT, V1~V3, value.var="V2")
rownames(CINPUT) <- CINPUT$V1
CINPUT <- (CINPUT[,-1])
CINPUT[is.na(CINPUT)] <- 0
TCINPUT <- t(CINPUT)

DistMatrix <- vegdist(CINPUT, method = "bray")
OrdMds <- metaMDS(DistMatrix)
OrdMdsFit = data.frame(MDS1 = OrdMds$points[,1], MDS2 = OrdMds$points[,2])
#NMDS_AND_MAP <- cbind(BRAY_ORD_FIT,INPUT_MAP)
ggplot(OrdMdsFit, aes(x=MDS1, y=MDS2)) + theme_classic() + geom_point(size=4) + scale_colour_brewer(palette="Set2")
