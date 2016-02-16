# ReturnModelResults.R
# Geoffrey Hannigan
# Pat Schoss Lab
# University of Michigan

library(reshape2)
library(vegan)
library(ape)
library(caret)
library(C50)
library(optparse)

# Parse the paramters from the command line.
option_list = list(
	make_option(c("-r", "--reference"), type="character", default=NULL, 
              help="reference kmer spectrum", metavar="character"),
	make_option(c("-t", "--test"), type="character", default=NULL, 
              help="test kmer spectrum", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

INPUT <- read.delim(opt$reference, sep="\t", header=F)
CONTIGS <- read.delim(opt$test, sep="\t", header=F)


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

# Get number of columns
ColCount <- ncol(ReferenceDf)
ColCountContigs <- ncol(ContigsDf)

# Create feature matrix and target vector
trainX <- ReferenceDf[,-ColCount]
trainX <- as.data.frame(sapply(trainX, function(x) x/sum(x)))
trainy <- as.factor(ReferenceDf[,ColCount])

testX <- ContigsDf[,-ColCount]
testX <- as.data.frame(sapply(testX, function(x) x/sum(x)))
testy <- as.factor(ContigsDf[,ColCount])

#Boost
model <-  C50::C5.0(trainX, trainy, trials=50)

pred <- predict(model, testX)
results <- postResample(pred, testy)

write(results, stdout())
