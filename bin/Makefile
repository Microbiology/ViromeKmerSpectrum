# Run a test of the distance script

./tmp/test.tsv ./tmp/testFormat.tsv : ../data/test.fa
	perl CalculateKmerDistancesPar.pl \
		-i $^ \
		-t $^ \
		-o ./tmp/test.tsv \
		-f ./tmp/testFormat.tsv \
		-w 24 \
		-p 64

