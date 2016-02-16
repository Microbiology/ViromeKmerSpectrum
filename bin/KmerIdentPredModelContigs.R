library(reshape2)
library(vegan)
library(ggplot2)
library(cluster)
library(ape)
library(caret)
library(C50)
library(randomForest)
library(ROCR)

INPUT <- read.delim("/Users/Hannigan/git/ViromeKmerSpectrum/data/trainingKmer.tsv", sep="\t", header=F)
CONTIGS <- read.delim("/Users/Hannigan/git/ViromeKmerSpectrum/data/testContigKmer.tsv", sep="\t", header=F)


FormatForPredModel <- function(input) {
  CINPUT <- dcast(input, V1~V3, value.var="V2")
  rownames(CINPUT) <- CINPUT$V1
  CINPUT <- (CINPUT[,-1])
  CINPUT[is.na(CINPUT)] <- 0
  TCINPUT <- as.data.frame(t(CINPUT))
  rownames(TCINPUT) <- gsub("\\;","_", rownames(TCINPUT))

  TCINPUT$ID <- sub("[p]hage.*", "phage", rownames(TCINPUT))
  rownames(TCINPUT) <- NULL
  return(TCINPUT)
}

ReferenceDf <- as.data.frame(FormatForPredModel(INPUT))
ContigsDf <- as.data.frame(FormatForPredModel(CONTIGS))

# Create feature matrix and target vector
trainX <- ReferenceDf[,1:256]
trainX <- as.data.frame(sapply(trainX, function(x) x/sum(x)))
trainy <- as.factor(ReferenceDf[,257])

testX <- ContigsDf[,1:256]
testX <- as.data.frame(sapply(testX, function(x) x/sum(x)))
testy <- as.factor(ContigsDf[,257])

#Boost
#model <-  C50::C5.0(trainX, trainy, trials=100)
#summary(model)

tuned <- train(trainX, factor(trainy), method = "gbm", trControl = trainControl(method = "repeatedcv", number=5, repeats = 5))

pred <- predict(model, testX)
postResample(pred, testy)
