# ProcessComparison.R
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#################
# Set Libraries #
#################

library("ggplot2")
library("optparse")

option_list <- list(
  make_option(c("-k", "--kmer"), type = "character", default = NULL,
      help = "Kmer data table", metavar = "character"),
  make_option(c("-b", "--blast"), type = "character", default = NULL,
      help = "Blast data table.", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
      help = "Output file for count summary", metavar = "character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

kmer <- read.delim(file=opt$kmer, sep="\t", header=F)
blast <- read.delim(file=opt$blast, sep="\t", header=F)

# I can write it this way beacuse this always has the max length
total <- length(kmer$V1)

kmercorrect <- 100 * sum(c(kmer$V1 == kmer$V2)) / total
blastcorrect <- 100 * sum(c(blast$V1 == blast$V2)) / total

values <- c(kmercorrect, blastcorrect)
labels <- c("Correct Kmer", "Correct Blast")

df <- data.frame(values, labels)

pdf(opt$out, height=4, width=3)
	ggplot(df, aes(x=labels, y=values)) +
		theme_classic() +
		theme(
			axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
			axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
		geom_bar(stat="identity", fill="tomato3") +
		ggtitle("Accuracy\nAssembling Genomes") +
		ylab("Correctly Identified Genomes\nAt Strain Level (%)") +
		xlab("Methods")
dev.off()