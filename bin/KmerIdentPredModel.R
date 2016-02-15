library(reshape2)
library(vegan)
library(ggplot2)
library(cluster)
library(ape)
library(caret)
library(C50)
library(randomForest)
library(ROCR)


INPUT <- read.delim("/Users/Hannigan/git/ViromeKmerSpectrum/data/kmerSpectrum.tsv", sep="\t", header=F)

CINPUT <- dcast(INPUT, V1~V3, value.var="V2")
rownames(CINPUT) <- CINPUT$V1
CINPUT <- (CINPUT[,-1])
CINPUT[is.na(CINPUT)] <- 0
TCINPUT <- as.data.frame(t(CINPUT))
colnames(CINPUT) <- gsub("\\;","_", colnames(CINPUT))
rownames(TCINPUT) <- gsub("\\;","_", rownames(TCINPUT))

PcaInput <- prcomp(CINPUT)
# Calculate the percent variance accounted for by each component
pcaPercentVar <- PcaInput$sd^2/sum(PcaInput$sd^2)*100
# Plot the variance explained by each component
screeplot(PcaInput, type = "lines")
PcaScree <- melt(pcaPercentVar)
PcaScree$name <- sequence(length(row.names(PcaScree)))
PcaScreePlot <- ggplot(PcaScree, aes(x=name, y=value)) + theme_classic() + geom_point() + geom_path() + xlab("PCA Component") + ylab("Percent Variance Explained")

pdf(file="/Users/Hannigan/git/ViromeKmerSpectrum/Figures/kmerSpectrumPcaScree.pdf", height=6, width=8)
PcaScreePlot
dev.off()

PCAsubset <- as.data.frame(PcaInput$rotation[,c(1:4)])
PCAsubset$ID <- sub("phage.*", "phage", rownames(PCAsubset))

TCINPUT$ID <- sub("[pP]hage.*", "phage", rownames(TCINPUT))
rownames(TCINPUT) <- NULL


crx <- TCINPUT[ sample( nrow( TCINPUT ) ), ]
# Create feature matrix and target vector
X <- crx[,1:256]
y <- as.factor(crx[,257])
trainX <- X[1:1000,]
trainy <- y[1:1000]
testX <- X[1001:2011,]
testy <- y[1001:2011]

#Boost
model <-  C50::C5.0( trainX, trainy, trials=10 )
summary( model )

p <- predict( model, testX, type="class" )
100 * sum( p == testy ) / length( p )


rf <-randomForest(factor(trainy)~.,data=trainX, ntree=1000, keep.forest=TRUE, importance=TRUE)

pred <- as.data.frame(predict(rf, testX, type="response"))
pred$ID <- as.numeric(rownames(pred))
colnames(pred) <- c("ID")
rownames(pred) <- NULL
pred <- as.factor(pred$ID)
pred <- factor(pred, levels=c(levels(testy)))

100 * sum( pred == testy ) / length( pred )


