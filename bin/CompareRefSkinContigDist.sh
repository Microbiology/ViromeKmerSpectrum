#! /bin/bash
# CompareRefsForClustering.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#PBS -N CompareRefsForClustering
#PBS -q first
#PBS -l nodes=1:ppn=1,mem=40gb
#PBS -l walltime=600:00:00
#PBS -j oe
#PBS -V
#PBS -A schloss_lab

# Set Variables
export WorkingDirectory=/Users/Hannigan/git/ViromeKmerSpectrum/data/
export Output='CompareRefs'

export LocalPath=/Users/Hannigan/git/ViromeKmerSpectrum/bin/
export FigurePath=/Users/Hannigan/git/ViromeKmerSpectrum/Figures/

export Genomes=/Users/Hannigan/git/ViromeKmerSpectrum/data/SkinMergeWithBacFungiNoBlockFormat.fa

# Move and make output dir
cd ${WorkingDirectory} || exit

mkdir ./${Output}

perl ../bin/CalculateKmerDistances.pl \
	-i ${Genomes} \
	-t ${Genomes} \
	-o ./${Output}/RefCompareSkin.tsv \
	-f ./${Output}/RefCompareSkinFormat.tsv \
	-w 5 \
	-r
