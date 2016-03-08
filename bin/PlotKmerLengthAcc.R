# PlotKmerLengthAcc.R
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# Read in libraries
library("optparse")
library("ggplot2")
library("RColorBrewer")

# Parse the paramters from the command line.
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Data file from kmer distance benchmark"),
	make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Destination for resulting figure (PDF)", metavar="character"),
	make_option(c("-p", "--png"), type="character", default=NULL, 
              help="Destination for resulting figure (PNG)", metavar="numeric"),
	make_option(c("-t", "--title"), type="character", default=NULL,
                        help="Title for the resulting plot", metavar="character"),
	make_option(c("-l", "--log"), action="store_true", default=FALSE,
			help="Should we use a log scale for y? [default %default]")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Read data into memory
RawKmerTime <- read.delim(file=opt$input, sep="\t", header=FALSE)

# Plot out the data
TimePlot <- ggplot(RawKmerTime, aes(x=factor(V2), y=V3, colour=V4)) +
	theme_classic() +
	geom_jitter() +
	ylab("Distance to Reference Genome (Inv Confidence)") +
	xlab("Kmer Length") +
	scale_colour_brewer(palette="Set2") +
	ggtitle(opt$title)

if (opt$log) {
   ComparePlot <- ComparePlot + scale_y_log10()
}

pdf(file=opt$out, width=8, height=6)
	TimePlot
dev.off()

png(file=opt$png, width=8, height=6, units="in", res=300)
	TimePlot
dev.off()
