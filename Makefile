# Makefile
# Geoffrey Hannigan

# Allow for everything to be run with all

OBJECTS = \
	./data/phage.txt ./data/bacteria.txt ./data/virus.txt ./data/eukaryota-tmp.txt \
	./data/eukaryota.txt \
	./data/PhageRef.fa ./data/BacteriaRef.fa ./data/VirusRef.fa ./data/EukaryotaRef.fa \
	./data/PhageRefSub.fa ./data/BacteriaRefSub.fa ./data/VirusRefSub.fa ./data/EukaryotaRefSub.fa \
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

######################
# Subset Fasta Files #
######################
# Phage
./data/PhageRefSub.fa : ./data/PhageRef.fa
	./bin/seqtk sample -s 100 ./data/PhageRef.fa 100 > ./data/PhageRefSub.fa

# Virus
./data/VirusRefSub.fa : ./data/VirusRef.fa
	./bin/seqtk sample -s 100 ./data/VirusRef.fa 100 > ./data/VirusRefSub.fa

# Bacteria
./data/BacteriaRefSub.fa : ./data/BacteriaRef.fa
	./bin/seqtk sample -s 100 ./data/BacteriaRef.fa 25 > ./data/BacteriaRefSub.fa

# Eukaryotes
./data/EukaryotaRefSub.fa : ./data/EukaryotaRef.fa
	./bin/seqtk sample -s 100 ./data/EukaryotaRef.fa 5 > ./data/EukaryotaRefSub.fa

#####################################
# Cluster Phages and Other Microbes #
#####################################
# Get together the fasta files
./data/allReferences.fa : ./data/PhageRefSub.fa ./data/BacteriaRefSub.fa ./data/VirusRefSub.fa ./data/EukaryotaRefSub.fa
	cat ./data/PhageRefSub.fa ./data/BacteriaRefSub.fa ./data/VirusRefSub.fa ./data/EukaryotaRefSub.fa > ./data/allReferences.fa

./data/CompareRefs/comparerefs.tsv ./data/CompareRefs/comparerefs-format.tsv : ./data/allReferences.fa
	bash ./bin/CompareRefsForClustering.sh \
		"CompareRefs" \
		./data/allReferences.fa \
		"comparerefs"

##################
# Cluster Phages #
##################
./data/CompareRefs/comparephages.tsv ./data/CompareRefs/comparephages-format.tsv : ./data/PhageRefSub.fa
	bash ./bin/CompareRefsForClustering.sh \
		"CompareRefs" \
		./data/PhageRefSub.fa \
		"comparephages"

./data/CompareRefs/compareviruses.tsv ./data/CompareRefs/compareviruses-format.tsv : ./data/VirusRefSub.fa
	bash ./bin/CompareRefsForClustering.sh \
		"CompareRefs" \
		./data/VirusRefSub.fa \
		"compareviruses"

####################
# Write Manuscript #
####################

./doc/manuscript.pdf:./doc/manuscript.md
	./doc/rendermanuscript
