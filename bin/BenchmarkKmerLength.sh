# BenchmarkKmerLength.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# NOTE: I wrote this to be run locally, not on a server

# Set the variables to be used in this script
export WorkingDirectory=/Users/Hannigan/git/ViromeKmerSpectrum/data/
export Output='LengthBenchmark'

export LocalPath=/Users/Hannigan/git/ViromeKmerSpectrum/bin/

export TrainingSet=/Users/Hannigan/git/ViromeKmerSpectrum/data/BlastnBenchmark/trainingSet.fa
export TestSet=/Users/Hannigan/git/ViromeKmerSpectrum/data/BlastnBenchmark/testSet.fa

# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}
# Appending output so ensure it is clean going in

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
		-w ${w}

	# Run kmer spectra through predictive model and collect output
	Rscript ${LocalPath}ReturnModelResults.R \
		-r ./${Output}/trainingKmer.tsv \
		-t ./${Output}/testingKmer.tsv \
		>> ./${Output}/ModelAccuracyOutput.tsv
}

export -f ModelWithLength

for iterator in $(seq 4 4 60); do
	ModelWithLength \
		${iterator} \
		${TrainingSet} \
		${TestSet}
done
