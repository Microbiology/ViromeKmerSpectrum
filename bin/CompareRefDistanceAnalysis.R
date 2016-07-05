# CompareRedDistanceAnalysis.R
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#################
# Set Libraries #
#################

library("ggplot2")
library("optparse")
library("reshape2")
library("vegan")
library("RColorBrewer")
library("cluster")
library("flexclust")

option_list <- list(
  make_option(c("-d", "--distance"), type = "character", default = NULL,
      help = "Input distance file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
      help = "Output file for count summary", metavar = "character")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

inputtable <- read.delim(file = "../data/CompareRefs/RefCompareSkinFormat.tsv", sep = "\t", header = F)

dcasttable <- dcast(inputtable, formula = V1~V3, value.var = "V2")
rownames(dcasttable) <- dcasttable$V1
dcasttable <- dcasttable[,-1]
ordnmds <- metaMDS(dcasttable,k=2)
ordnmdsfit = data.frame(MDS1 = ordnmds$points[,1], MDS2 = ordnmds$points[,2])
ordnmdsfit$ID <- rownames(ordnmdsfit)
ordnmdsfit$ID <- sub("[Pp]hage_.*", "phage", ordnmdsfit$ID, perl=TRUE)
ordnmdsfit$ID <- sub("contig.*", "Contig_phage", ordnmdsfit$ID, perl=TRUE)

ordnmdsfit <- ordnmdsfit[ordnmdsfit$ID %in% names(table(ordnmdsfit$ID))[table(ordnmdsfit$ID) >= 0],]

plot1 <- ggplot(ordnmdsfit, aes(x=MDS1, y=MDS2, colour=ID)) +
		theme_classic() +
		theme(
			axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
			axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
		geom_point(size=2.5) +
		scale_colour_brewer(palette="Set1") +
		ggtitle("Kmer Bray-Curtis")

# Try out the clustering analysis
wssplot <- function(data, nc, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}

wssplot(dcasttable,nc=15)

set.seed(1234)
fit<-kmeans(dcasttable,3)
clusplot(as.matrix(dcasttable), fit$cluster,shade=TRUE, color=TRUE,main="K-means 3 Clustering- Virome - Time Points 2&3")

ct.km <- table(ordnmdsfit$ID, fit$cluster)
randIndex(ct.km)

# Merge in the classes
classes <- data.frame(fit$cluster)
classes$IDnew <- rownames(classes)
ordnmdsfit$IDnew <- rownames(ordnmdsfit)
mergedord <- merge(ordnmdsfit, classes, by="IDnew")
mergedord <- mergedord[mergedord$ID %in% names(table(mergedord$ID))[table(mergedord$ID) >= 3],]

plot2 <- ggplot(mergedord, aes(x=MDS1, y=MDS2, colour=ID, shape=factor(fit.cluster))) +
		theme_classic() +
		theme(
			axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
			axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
		geom_point(size=2.5) +
		scale_colour_brewer(palette="Set1") +
		ggtitle("Kmer Bray-Curtis")
plot2

pdf(file="../Figures/RefOrdination.pdf", width=6, height=4)
	plot2
dev.off()
