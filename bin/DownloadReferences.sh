#! /bin/bash
# Geoffrey Hannigan

DownloadGenomes () {
	mkdir ./tmp
	for ID in $(cut -c1-2 ${1} | sort | uniq); do
		export AccString=$(cut -f 1 ${1} | egrep "^${ID}" | tr '\n' ',' | gsed 's/,$//')
		wget "http://www.ebi.ac.uk/ena/data/view/${AccString}&display=fasta" -O ./tmp/Ref_${ID}.fa
	done
	cat ./tmp/Ref_*.fa > ${2}
	rm -r ./tmp
}

export -f DownloadGenomes

# Download the phages
DownloadGenomes \
	./data/phage.txt \
	./data/PhageRef.fa
# Download the bacteria
DownloadGenomes \
	./data/bacteria.txt \
	./data/BacteriaRef.fa
# Download the viruses
DownloadGenomes \
	./data/virus.txt \
	./data/VirusRef.fa
# Download the eukaryotes
DownloadGenomes \
	./data/eukaryota.txt \
	./data/EukaryotaRef.fa
