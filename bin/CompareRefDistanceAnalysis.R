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

option_list <- list(
  make_option(c("-d", "--distance"), type = "character", default = NULL,
      help = "Input distance file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
      help = "Output file for count summary", metavar = "character")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

inputtable <- read.delim(file = opt$distance, sep = "\t", header = F)

dcasttable <- dcast(inputtable, formula = V1~V3, value.var = "V2")
rownames(dcasttable) <- dcasttable$V1
dcasttable <- dcasttable[,-1]
ordnmds <- metaMDS(dcasttable,k=2)
ordnmdsfit = data.frame(MDS1 = ordnmds$points[,1], MDS2 = ordnmds$points[,2])
ordnmdsfit$ID <- rownames(ordnmdsfit)
ordnmdsfit$ID <- sub("[Pp]hage_.*", "phage", ordnmdsfit$ID, perl=TRUE)

ordnmdsfit <- ordnmdsfit[ordnmdsfit$ID %in% names(table(ordnmdsfit$ID))[table(ordnmdsfit$ID) >= 3],]

plot <- ggplot(ordnmdsfit, aes(x=MDS1, y=MDS2, colour=ID)) +
		theme_classic() +
		theme(
			axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
			axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
		geom_point(size=2.5) +
		scale_colour_brewer(palette="Set1") +
		ggtitle("Kmer Bray-Curtis")

pdf(file=opt$output, width=6, height=4)
	plot
dev.off()
