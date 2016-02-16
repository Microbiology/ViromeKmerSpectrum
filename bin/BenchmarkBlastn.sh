# BenchmarkBlastn.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# NOTE: I wrote this to be run locally, not on a server

# Set the variables to be used in this script
export WorkingDirectory=/Users/Hannigan/git/ViromeKmerSpectrum/data/
export Output='BlastnBenchmark'

export PhageReference=/Users/Hannigan/git/ViromeKmerSpectrum/data/PhageSVAFormat.fa

export LocalPath=/Users/Hannigan/git/ViromeKmerSpectrum/bin/
export BlastPath=/Users/Hannigan/Documents/ncbi-blast-2.3.0+/bin/
export OpenMetToolkit=/Users/Hannigan/git/OpenMetagenomeToolkit/

# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}

RandomBlast () {
	# 1 = Phage genome reference
	# 2 = Final Output Name

	# First split the fasta into a test and training dataset
	echo Splitting reference file...
	perl ${OpenMetToolkit}RandomFastaSplit.pl \
		-i ${1} \
		-c ./${Output}/trainingSet.fa \
		-r ./${Output}/testSetWhole.fa \
		-l 1000

	# To look at contigs, subset the test set
	perl ${LocalPath}RandomContigGenerator.pl \
		-i ./${Output}/testSetWhole.fa \
		-o ./${Output}/testSet.fa

	# Use the training set as the reference database
	echo Training database...
	${BlastPath}makeblastdb \
		-in ./${Output}/trainingSet.fa \
		-dbtype nucl \
		-out ./${Output}/TrainedDb

	# Run the test dataset against the trained database
	echo Testing blast against training database...
	${BlastPath}blastn \
			-query ./${Output}/testSet.fa \
			-db ./${Output}/TrainedDb \
			-evalue 1e-3 \
			-outfmt 6 \
			-num_threads 4 \
		| sort -k1,1 -k12,12nr -k11,11n \
		| sort -u -k1,1 --merge \
		> ./${Output}/BlastnTestResults.tsv

	cut -f 1,2 ./${Output}/BlastnTestResults.tsv \
		| gsed 's/phage\S*/phage/g' \
		| awk '$1 == $2 { print $0 }' \
		| awk 'END { print 100*NR/1011 }' \
		> ./${Output}/${2}
		# The 1011 comes from number of seqs in test fasta
}

export -f RandomBlast

RandomBlast \
	${PhageReference} \
	BlastnAccuracy.txt

