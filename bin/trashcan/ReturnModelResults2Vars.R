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
library(doMC)
registerDoMC(8)

# Parse the paramters from the command line.
option_list = list(
	make_option(c("-r", "--reference"), type="character", default=NULL, 
              help="short reference kmer spectrum", metavar="character"),
	make_option(c("-t", "--test"), type="character", default=NULL, 
              help="short test kmer spectrum", metavar="character"),
  make_option(c("-h", "--highreference"), type="character", default=NULL, 
              help="long test kmer spectrum", metavar="character"),
  make_option(c("-e", "--hightest"), type="character", default=NULL, 
              help="long test kmer spectrum", metavar="character"),
  make_option(c("-a", "--arow"), type="numeric", default=NULL, 
              help="row number for reference", metavar="numeric"),
  make_option(c("-b", "--brow"), type="numeric", default=NULL, 
              help="row number for test", metavar="numeric"),
  make_option(c("-c", "--crow"), type="numeric", default=NULL, 
              help="row number for long reference", metavar="numeric"),
  make_option(c("-d", "--drow"), type="numeric", default=NULL, 
              help="row number for long test", metavar="numeric")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

FormatForPredModel <- function(input) {
  CINPUT <- dcast(input, V1~V3, value.var="V2")
  rownames(CINPUT) <- CINPUT$V1
  CINPUT <- (CINPUT[,-1])
  CINPUT[is.na(CINPUT)] <- 0
  TCINPUT <- as.data.frame(t(as.data.frame(CINPUT)))
  remove(CINPUT)
  rownames(TCINPUT) <- gsub("\\;","_", rownames(TCINPUT))
  TCINPUT$ID <- sub("[p]hage.*", "phage", rownames(TCINPUT))
  rownames(TCINPUT) <- NULL
  return(TCINPUT)
}

INPUT <- read.delim(file=opt$reference, sep="\t", header=F, colClasses = c("factor", "integer", "factor"), nrows=opt$arow, comment.char="")
ReferenceDf <- as.data.frame(FormatForPredModel(INPUT))
remove(INPUT)

SECONDINPUT <- read.delim(file=opt$highreference, sep="\t", header=F, colClasses = c("factor", "integer", "factor"), nrows=opt$crow, comment.char="")
LongReferenceDf <- as.data.frame(FormatForPredModel(SECONDINPUT))
remove(SECONDINPUT)

ReferenceDf <- cbind(ReferenceDf, LongReferenceDf)
remove(LongReferenceDf)

CONTIGS <- read.delim(file=opt$test, sep="\t", header=F, colClasses = c("factor", "integer", "factor"), nrows=opt$brow, comment.char="")
ContigDf <- as.data.frame(FormatForPredModel(CONTIGS))
remove(CONTIGS)

SECONDCONTIGS <- read.delim(file=opt$hightest, sep="\t", header=F, colClasses = c("factor", "integer", "factor"), nrows=opt$drow, comment.char="")
LongContigDf <- as.data.frame(FormatForPredModel(SECONDCONTIGS))
remove(SECONDCONTIGS)

ReferenceDf <- cbind(ContigDf, LongContigDf)
remove(LongContigDf)


# Get number of columns
ColCount <- ncol(ReferenceDf)
ColCountContigs <- ncol(ContigsDf)

# Create feature matrix and target vector
trainX <- ReferenceDf[,-ColCount]
trainX <- as.data.frame(sapply(trainX, function(x) x/sum(x)))
trainy <- as.factor(ReferenceDf[,ColCount])
remove()

testX <- ContigsDf[,-ColCount]
testX <- as.data.frame(sapply(testX, function(x) x/sum(x)))
testy <- as.factor(ContigsDf[,ColCount])

#Boost
model <-  C50::C5.0(trainX, trainy, trials=25)

#print("Testing model", stderr())

pred <- predict(model, testX)
results <- postResample(pred, testy)

write(results, stdout())
