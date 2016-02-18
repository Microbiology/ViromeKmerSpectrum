# BenchmarkKmerLength.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# NOTE: I wrote this to be run locally, not on a server

# Set the variables to be used in this script
export WorkingDirectory=/Home/ghannig/git/ViromeKmerSpectrum/data/
export Output='LengthBenchmark'

export LocalPath=/Home/ghannig/git/ViromeKmerSpectrum/bin/

export TrainingSet=/Home/ghannig/git/ViromeKmerSpectrum/data/BlastnBenchmark/trainingSet.fa
export TestSet=/Home/ghannig/git/ViromeKmerSpectrum/data/BlastnBenchmark/testSet.fa

# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}

ModelWithLength () {
	# 1 = Window size
	# 2 = Input training set
	# 3 = Input testing set
	echo Window size is ${1}

	perl ${LocalPath}CalculateKmerSpectrum.pl \
		-i ${2} \
		-o ./${Output}/trainingKmer.tsv \
		-w ${1}

	perl ${LocalPath}CalculateKmerSpectrum.pl \
		-i ${3} \
		-o ./${Output}/testingKmer.tsv \
		-w ${1}

	arow=$(wc -l ./${Output}/trainingKmer.tsv | gsed 's/^ *//' | gsed 's/ .*//')
	brow=$(wc -l ./${Output}/testingKmer.tsv | gsed 's/^ *//' | gsed 's/ .*//')
	echo Train length is ${arow}
	echo Test length is ${brow}

	# Run kmer spectra through predictive model and collect output
	Rscript ${LocalPath}ReturnModelResults.R \
		-r ./${Output}/trainingKmer.tsv \
		-t ./${Output}/testingKmer.tsv \
		-a ${arow} \
		-b ${brow} \
		> ./${Output}/${1}-ModelAccuracyOutput.tsv
}

export -f ModelWithLength

for i in $(seq 4 4 40); do
	ModelWithLength \
		${i} \
		${TrainingSet} \
		${TestSet}
done
