# BenchmarkKmerLength.sh
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# NOTE: I wrote this to be run locally, not on a server

# Set the variables to be used in this script
export WorkingDirectory=/Users/Hannigan/git/ViromeKmerSpectrum/data/
export Output='LengthBenchmark'

export LocalPath=/Users/Hannigan/git/ViromeKmerSpectrum/bin/
export FigurePath=/Users/Hannigan/git/ViromeKmerSpectrum/Figures/

export TrainingSet=/Users/Hannigan/git/ViromeKmerSpectrum/data/trainingSub.fa
export TestSet=/Users/Hannigan/git/ViromeKmerSpectrum/data/trainingSubContigs.fa

# Make the output directory and move to the working directory
echo Creating output directory...
cd ${WorkingDirectory}
mkdir ./${Output}

ModelWithLength () {
	# 1 = Window size
	# 2 = Input training set
	# 3 = Input testing set
	echo Window size is ${1}

	perl ${LocalPath}CalculateKmerDistances.pl \
		-i ${2} \
		-t ${3} \
		-o ./${Output}/${1}-Result.tsv \
		-f ./${Output}/${1}-ResultFormat.tsv \
		-w ${1} \
		-r \
		> ./${Output}/${1}-BenchmarkTime.tsv
}

export -f ModelWithLength

for i in $(seq 4 4 44); do
	ModelWithLength \
		${i} \
		${TrainingSet} \
		${TestSet}
done

# Combine the time files
cat \
	./${Output}/*-BenchmarkTime.tsv \
	> ./${Output}/TotalBenchTime.tsv

# Get only the RUNTIME lines
grep RUNTIME \
	./${Output}/TotalBenchTime.tsv \
	> ./${Output}/TotalBenchTimeFormat.tsv

# Plot out the time over kmers
Rscript ${LocalPath}PlotKmerTimeLength.R \
	-i ./${Output}/TotalBenchTimeFormat.tsv \
	-o ${FigurePath}TotalBenchTimeContigs.pdf \
	-p ${FigurePath}TotalBenchTimeContigs.png \
	-t "Kmer Length vs Time Benchmark"

rm ./${Output}/ResultingDistance.tsv

# Get confidence score in annotations by length
for i in $(seq 4 4 44); do
	echo Number is $i
	awk -v seq=$i \
		'{ if ($1 == $3) {print $1"\t"seq"\t"$2"\tTRUE_ID"} else {print $1"\t"seq"\t"$2"\tFALSE_ID"} }' \
		./${Output}/${i}-ResultFormat.tsv \
		>> ./${Output}/ResultingDistance.tsv
done

Rscript ${LocalPath}PlotKmerLengthAcc.R \
	-i ./${Output}/ResultingDistance.tsv \
	-o ${FigurePath}TotalBenchAccContigs.pdf \
	-p ${FigurePath}TotalBenchAccContigs.png \
	-t "Kmer Length vs Confidence"

