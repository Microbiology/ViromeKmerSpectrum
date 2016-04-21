#! /usr/bin/
# Benchmark.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

#################
# Set Variables #
#################

# Paths
export WorkingDirectory='../data/'
export Output='Benchmark'
export seqtkprog='/Users/Hannigan/Documents/seqtk-1.1/seqtk'
export ncbibin='/Users/Hannigan/Documents/ncbi-blast-2.3.0+/bin/'

# Files
export PhageFa=./PhageSVAFormat.fa

# Move and make output dir
cd ${WorkingDirectory} || exit

mkdir ./${Output}

# Make tmp directory
mkdir ./tmpbench

# Randomly subsample reference fasta
${seqtkprog} sample ${PhageFa} 100 > ./tmpbench/subsample.fa

# ###################################
# # Confirm Accuracy Against Itself #
# ###################################
# # Kmer
# perl ../bin/CalculateKmerDistancesPar.pl \
# 	-i ./tmpbench/subsample.fa \
# 	-t ./tmpbench/subsample.fa \
# 	-o ./tmpbench/BenchPerlComplete.tsv \
# 	-f ./tmpbench/BenchPerlCompleteFormat.tsv \
# 	-w 8 \
# 	-p 64

# ${ncbibin}makeblastdb \
# 	-dbtype nucl \
# 	-in ./tmpbench/subsample.fa \
# 	-out ./tmpbench/subsampleblastref

# ${ncbibin}blastn \
# 	-query ./tmpbench/subsample.fa \
# 	-out ./tmpbench/BenchBlastnComplete.tsv \
# 	-db ./tmpbench/subsampleblastref \
# 	-outfmt 6 \
# 	-num_threads 64 \
# 	-max_target_seqs 1 \
# 	-evalue 1e-3

################################
# Accuracy With Contig Subsets #
################################
perl ../bin/RandomContigGenerator.pl \
	-i ./tmpbench/subsample.fa \
	-o ./tmpbench/subsampleContigs.fa

perl ../bin/CalculateKmerDistancesPar.pl \
	-i ./tmpbench/subsampleContigs.fa \
	-t ./tmpbench/subsampleContigs.fa \
	-o ./tmpbench/BenchPerlCompleteContigs.tsv \
	-f ./tmpbench/BenchPerlCompleteContigsFormat.tsv \
	-w 8 \
	-p 64

${ncbibin}blastn \
	-query ./tmpbench/subsampleContigs.fa \
	-out ./tmpbench/BenchBlastnCompleteContigs.tsv \
	-db ./tmpbench/subsampleblastref \
	-outfmt 6 \
	-num_threads 64 \
	-max_target_seqs 1 \
	-evalue 1e-3



