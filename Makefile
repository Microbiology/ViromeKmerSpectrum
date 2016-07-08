# Makefile
# Geoffrey Hannigan

# Allow for everything to be run with all

OBJECTS = \
	./data/phage.txt ./data/bacteria.txt ./data/virus.txt ./data/eukaryota-tmp.txt \
	./data/eukaryota.txt \
	./data/PhageRef.fa ./data/BacteriaRef.fa ./data/VirusRef.fa ./data/EukaryotaRef.fa \
	./data/allReferences.fa \
	./data/CompareRefs/comparerefs.tsv ./data/CompareRefs/comparerefs-format.tsv \
	./data/CompareRefs/comparephages.tsv ./data/CompareRefs/comparephages-format.tsv \
	./data/CompareRefs/compareviruses.tsv ./data/CompareRefs/compareviruses-format.tsv

all: $(OBJECTS)

#################
# Download Data #
#################
# Get the phage accession list from EBI
./data/phage.txt ./data/bacteria.txt ./data/virus.txt ./data/eukaryota-tmp.txt :
	wget http://www.ebi.ac.uk/genomes/phage.txt -O ./data/phage.txt
	wget http://www.ebi.ac.uk/genomes/bacteria.txt -O ./data/bacteria.txt
	wget http://www.ebi.ac.uk/genomes/virus.txt -O ./data/virus.txt
	wget http://www.ebi.ac.uk/genomes/eukaryota.txt -O ./data/eukaryota-tmp.txt

# Only get a subset of the eukaryotes since they are huge files
# I'm pulling a random subset here
./data/eukaryota.txt : ./data/eukaryota-tmp.txt
	shuf ./data/eukaryota-tmp.txt | head -n 25 > ./data/eukaryota.txt; rm ./data/eukaryota-tmp.txt

# Must be downloaded in chunks or else there will be errors
# due to it trying to download multiple file classes.
./data/PhageRef.fa ./data/BacteriaRef.fa ./data/VirusRef.fa ./data/EukaryotaRef.fa : ./data/phage.txt ./data/bacteria.txt ./data/virus.txt ./data/eukaryota.txt
	bash ./bin/DownloadReferences.sh

#####################################
# Cluster Phages and Other Microbes #
#####################################
# Get together the fasta files
./data/allReferences.fa : ./data/PhageRef.fa ./data/BacteriaRef.fa ./data/VirusRef.fa ./data/EukaryotaRef.fa
	cat ./data/PhageRef.fa ./data/BacteriaRef.fa ./data/VirusRef.fa ./data/EukaryotaRef.fa > ./data/allReferences.fa

./data/CompareRefs/comparerefs.tsv ./data/CompareRefs/comparerefs-format.tsv : ./data/allReferences.fa
	bash ./bin/CompareRefsForClustering.sh \
		"CompareRefs" \
		./data/allReferences.fa \
		"comparerefs"

##################
# Cluster Phages #
##################
./data/CompareRefs/comparephages.tsv ./data/CompareRefs/comparephages-format.tsv : ./data/PhageRef.fa
	bash ./bin/CompareRefsForClustering.sh \
		"CompareRefs" \
		./data/PhageRef.fa \
		"comparephages"

./data/CompareRefs/compareviruses.tsv ./data/CompareRefs/compareviruses-format.tsv : ./data/VirusRef.fa
	bash ./bin/CompareRefsForClustering.sh \
		"CompareRefs" \
		./data/VirusRef.fa \
		"compareviruses"

####################
# Write Manuscript #
####################

./doc/manuscript.pdf:./doc/manuscript.md
	./doc/rendermanuscript
