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
	./data/CompareRefs/comparephages.tsv ./data/CompareRefs/comparephages-format.tsv ./data/CompareRefs/comparephages-formatfinal.tsv

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

# Bacteria
./data/BacteriaRefSub.fa : ./data/BacteriaRef.fa
	./bin/seqtk sample -s 100 ./data/BacteriaRef.fa 10 > ./data/BacteriaRefSub.fa

# Eukaryotes
./data/EukaryotaRefSub.fa : ./data/EukaryotaRef.fa
	./bin/seqtk sample -s 100 ./data/EukaryotaRef.fa 3 > ./data/EukaryotaRefSub.fa

#####################################
# Cluster Phages and Other Microbes #
#####################################
# Get together the fasta files
./data/allReferences.fa : ./data/PhageRefSub.fa ./data/BacteriaRefSub.fa ./data/EukaryotaRefSub.fa
	cat ./data/PhageRefSub.fa ./data/BacteriaRefSub.fa ./data/EukaryotaRefSub.fa > ./data/allReferences.fa

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

# Format the phages
./data/CompareRefs/comparephages-formatfinal.tsv : ./data/CompareRefs/comparephages-format.tsv
	sed 's/[Pp]hage[^\t]*\t/Phage\t/' ./data/CompareRefs/comparephages-format.tsv \
		| sed 's/[Pp]hage[^\t]*$/Phage/' \
		| sed 's/ENA|\(\S*\)|\S* /\1 /g' \
		| sed 's/\(\S* \S* \S* \)[^\t]*Phage/\1Phage/g' \
		| sed 's/\(\S*\) \([^\t]*\)/\2 \1/g' \
		| sed 's/\(\S*\) \(\S*\) Phage/\1__\2 Phage/g' \
		| sed 's/ /_/g' \
		| sed '$d' \
		> ./data/CompareRefs/comparephages-formatfinal.tsv

####################
# Write Manuscript #
####################

./doc/manuscript.pdf:./doc/manuscript.md
	./doc/rendermanuscript
