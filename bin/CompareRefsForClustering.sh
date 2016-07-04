#! /bin/bash
# CompareRefsForClustering.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# Set Variables
export WorkingDirectory=/Users/Hannigan/git/ViromeKmerSpectrum/data/
export Output='CompareRefs'

export LocalPath=/Users/Hannigan/git/ViromeKmerSpectrum/bin/
export FigurePath=/Users/Hannigan/git/ViromeKmerSpectrum/Figures/

export Genomes=/Users/Hannigan/git/ViromeKmerSpectrum/data/PhageSVAFormat.fa

# Load in the proper perl module
module load perl/5.22.1 

# Move and make output dir
cd ${WorkingDirectory} || exit

mkdir ./${Output}

perl ../bin/CalculateKmerDistances.pl \
	-i ${Genomes} \
	-t ${Genomes} \
	-o ./${Output}/RefCompare.tsv \
	-f ./${Output}/RefCompareFormat.tsv \
	-w 5 \
	-r
