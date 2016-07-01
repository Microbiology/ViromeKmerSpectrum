#! /bin/bash
# BenchmarkGenomeSplits.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# Set Variables
export WorkingDirectory=/Users/Hannigan/git/ViromeKmerSpectrum/data/
export Output='LengthBenchmark'

export LocalPath=/Users/Hannigan/git/ViromeKmerSpectrum/bin/
export FigurePath=/Users/Hannigan/git/ViromeKmerSpectrum/Figures/
export ncbibin=/Users/Hannigan/Documents/ncbi-blast-2.3.0+/bin/

export Genome=/Users/Hannigan/git/ViromeKmerSpectrum/data/tmpbench/subsample.fa

# Move and make output dir
cd ${WorkingDirectory} || exit

mkdir ./${Output}

perl ../bin/GenerateHalfGenomes.pl \
	-i ${Genome} \
	-o ./Genome1.fa \
	-t ./Genome2.fa

perl ../bin/CalculateKmerDistancesPar.pl \
	-i ./Genome1.fa \
	-t ./Genome2.fa \
	-o ./Kmerout.tsv \
	-f ./KmerOutFormat.tsv \
	-w 5 \
	-p 64

${ncbibin}makeblastdb \
	-dbtype nucl \
	-in ./Genome1.fa \
	-out ./Genome1Ref

${ncbibin}blastn \
	-query ./Genome2.fa \
	-out ./BlastOut.tsv \
	-db ./Genome1Ref \
	-outfmt 6 \
	-num_threads 4 \
	-max_target_seqs 1 \
	-evalue 1e-5

# Final kmer format
sort -k 1,1 -k 2,2 ./KmerOutFormat.tsv \
	| awk -F"[. ]" '!a[$1]++' \
	| awk '{ print $1"\t"$3"\t"$2 }' \
	| awk {' print $0"\tKMER" '} \
	> ./KmerForAnalysis.tsv

cut -f 1,2,12 ./BlastOut.tsv \
	| sort -nrk 3,3 \
	| awk '!a[$1]++' \
	| awk {' print $0"\tBLAST" '} \
	> BlastForAnalysis.tsv

Rscript ../bin/ProcessComparison.R \
	-k KmerForAnalysis.tsv \
	-b BlastForAnalysis.tsv \
	-o ${FigurePath}ComparisonBarGraph
