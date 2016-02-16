library(reshape2)
library(vegan)
library(ggplot2)
library(cluster)
library(ape)
library(caret)
library(C50)
library(randomForest)
library(ROCR)
#library(doMC)
#registerDoMC(4)


INPUT <- read.delim("/Users/Hannigan/git/ViromeKmerSpectrum/data/kmerSpectrum.tsv", sep="\t", header=F)

CINPUT <- dcast(INPUT, V1~V3, value.var="V2")
rownames(CINPUT) <- CINPUT$V1
CINPUT <- (CINPUT[,-1])
CINPUT[is.na(CINPUT)] <- 0
TCINPUT <- as.data.frame(t(CINPUT))
colnames(CINPUT) <- gsub("\\;","_", colnames(CINPUT))
rownames(TCINPUT) <- gsub("\\;","_", rownames(TCINPUT))

TCINPUT$ID <- sub("[p]hage.*", "phage", rownames(TCINPUT))
TCINPUT$ID <- sub("[p]hage.*", "phage", rownames(TCINPUT))
rownames(TCINPUT) <- NULL

########################
# C50 Prediction Model #
########################

crx <- TCINPUT[ sample( nrow( TCINPUT ) ), ]
# Create feature matrix and target vector
X <- crx[,1:256]
y <- as.factor(crx[,257])
trainX <- X[1:1000,]
trainX <- as.data.frame(sapply(trainX, function(x) x/sum(x)))
trainy <- y[1:1000]

testX <- X[1001:2011,]
testX <- as.data.frame(sapply(testX, function(x) x/sum(x)))
testy <- y[1001:2011]

# # Train
# tuned <- train(trainX, factor(trainy), method = "C5.0", tuneLength = 11, trControl = trainControl(method = "repeatedcv", repeats = 5), metric = "Kappa")

#Boost
model <-  C50::C5.0(trainX, trainy, trials=75, rules=TRUE)
summary(model)

p <- predict(model, testX, type="class")
C50Result <- 100 * sum(p == testy) / length(p)
C50Result

##################################
# Random Forest Prediction Model #
##################################

rf <-randomForest(factor(trainy)~.,data=trainX, ntree=1000, keep.forest=TRUE, importance=TRUE)
summary(rf)

pred <- as.data.frame(predict(rf, testX, type="response"))
pred$ID <- as.numeric(rownames(pred))
colnames(pred) <- c("ID")
rownames(pred) <- NULL
pred <- as.factor(pred$ID)
pred <- factor(pred, levels=c(levels(testy)))

RfResult <- 100 * sum( pred == testy ) / length( pred )
RfResult

