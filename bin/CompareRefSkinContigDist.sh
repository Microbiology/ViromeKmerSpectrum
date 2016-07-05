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

export Genomes=/Users/Hannigan/git/ViromeKmerSpectrum/data/SkinMergeWithBacFungiNoBlock.fa

# Move and make output dir
cd ${WorkingDirectory} || exit

mkdir ./${Output}

perl ../bin/CalculateKmerDistances.pl \
	-i ${Genomes} \
	-t ${Genomes} \
	-o ./${Output}/RefCompareSkin.tsv \
	-f ./${Output}/RefCompareSkinFormat.tsv \
	-w 5 \
	-p 16 \
	-r

