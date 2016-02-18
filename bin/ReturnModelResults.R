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
              help="reference kmer spectrum", metavar="character"),
	make_option(c("-t", "--test"), type="character", default=NULL, 
              help="test kmer spectrum", metavar="character"),
  make_option(c("-a", "--arow"), type="numeric", default=NULL, 
              help="row number for reference", metavar="numeric"),
  make_option(c("-b", "--brow"), type="numeric", default=NULL, 
              help="row number for test", metavar="numeric"),
  make_option(c("-s", "--scree"), type="character", default=NULL, 
              help="scree figure output", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("Reading in files", stderr())

INPUT <- read.delim(file=opt$reference, sep="\t", header=F, colClasses = c("factor", "integer", "factor"), nrows=opt$arow, comment.char="")
CONTIGS <- read.delim(file=opt$test, sep="\t", header=F, colClasses = c("factor", "integer", "factor"), nrows=opt$brow, comment.char="")


FormatForPredModel <- function(input) {
  CINPUT <- dcast(input, V1~V3, value.var="V2")
  rownames(CINPUT) <- CINPUT$V1
  CINPUT <- (CINPUT[,-1])
  CINPUT[is.na(CINPUT)] <- 0
  # PCA
  PcaInput <- prcomp(CINPUT)
  pcaPercentVar <- PcaInput$sd^2/sum(PcaInput$sd^2)*100
  PcaScree <- melt(pcaPercentVar)
  PcaScree$name <- sequence(length(row.names(PcaScree)))
  PcaScreePlot <- ggplot(PcaScree, aes(x=name, y=value)) + theme_classic() + geom_point() + geom_path() + xlab("PCA Component") + ylab("Percent Variance Explained")
  # End PCA 
  TCINPUT <- as.data.frame(t(CINPUT))
  rownames(TCINPUT) <- gsub("\\;","_", rownames(TCINPUT))
  TCINPUT$ID <- sub("[p]hage.*", "phage", rownames(TCINPUT))
  rownames(TCINPUT) <- NULL
  return(list(TCINPUT, PcaScreePlot))
}

ReferenceList <- FormatForPredModel(INPUT)
ReferenceDf <- as.data.frame(ReferenceList[1])
ReferencePlot <- ReferenceList[2]

ContigList <- FormatForPredModel(CONTIGS)
ContigsDf <- as.data.frame(ContigList[1])
ContigsPlot <- ContigList[2]

pdf(file=opt$scree, width=8, height=6)
  ReferencePlot
dev.off()

# # Get number of columns
# ColCount <- ncol(ReferenceDf)
# ColCountContigs <- ncol(ContigsDf)

# # Create feature matrix and target vector
# trainX <- ReferenceDf[,-ColCount]
# # trainX <- as.data.frame(sapply(trainX, function(x) x/sum(x)))
# trainy <- as.factor(ReferenceDf[,ColCount])

# testX <- ContigsDf[,-ColCount]
# # testX <- as.data.frame(sapply(testX, function(x) x/sum(x)))
# testy <- as.factor(ContigsDf[,ColCount])

# #print("Training model", stderr())

# #Boost
# model <-  C50::C5.0(trainX, trainy, trials=25)

# #print("Testing model", stderr())

# pred <- predict(model, testX)
# results <- postResample(pred, testy)

# write(results, stdout())
