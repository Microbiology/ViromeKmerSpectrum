#! /bin/bash
# CompareRefsForClustering.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

module load perl/5.22.1
module load perl-modules/5.22.1

#PBS -N CompareRefsForClustering
#PBS -q first
#PBS -l nodes=1:ppn=1,mem=40gb
#PBS -l walltime=600:00:00
#PBS -j oe
#PBS -V
#PBS -A schloss_lab

# Set Variables
export WorkingDirectory=/mnt/EXT/Schloss-data/ghannig/ViromeKmerSpectrum/data/
export Output='CompareRefs'

export LocalPath=/mnt/EXT/Schloss-data/ghannig/ViromeKmerSpectrum/bin/
export FigurePath=/mnt/EXT/Schloss-data/ghannig/ViromeKmerSpectrum/Figures/

export Genomes=/mnt/EXT/Schloss-data/ghannig/ViromeKmerSpectrum/data/PhageSVAFormat.fa

# Move and make output dir
cd ${WorkingDirectory} || exit

mkdir ./${Output}

perl ../bin/CalculateKmerDistances.pl \
	-i ${Genomes} \
	-t ${Genomes} \
	-o ./${Output}/RefCompare.tsv \
	-f ./${Output}/RefCompareFormat.tsv \
	-w 5

