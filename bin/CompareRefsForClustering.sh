#! /bin/bash
# CompareRefsForClustering.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

module load perl/5.22.1
module load perl-modules/5.22.1

# Set Variables
export Output=$1

export Genomes=$2
export OutputFile=$3 # Without file ext please

# Move and make output dir
mkdir ./data/${Output}

perl ./bin/remove_block_fasta_format.pl ${Genomes} ${Genomes}.tmp

perl ./bin/CalculateKmerDistancesPar.pl \
	-i ${Genomes}.tmp \
	-t ${Genomes}.tmp \
	-o ./data/${Output}/${OutputFile}.tsv \
	-f ./data/${Output}/${OutputFile}-format.tsv \
	-w 5 \
	-r \
	-p 4

# Clean up
rm ${Genomes}.tmp
