#!/usr/bin/perl
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# Set use
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
# Timer
my $start_run = time();

# Set variables
my %FastaHash;
my %fastaHash;
my %windowHash;
my %kmerTableHash;
my %KmerReadHash;
my $window = 4;
my $opt_help;
my $input;
my $output;
my $FastaErrorCounter = 0;
my $FastaID;
my $length = 0;
my $correct = 0;
my $windowValue;
my $key;
my $skip = 0;
my $maxCount = 1000;
my $Mean;
my $MeanCounter;
my $MeanSum;

# Set the options
GetOptions(
	'h|help' => \$opt_help,
	'i|input=s' => \$input,
	'o|output=s' => \$output,
	'w|window=n' => \$window,
	'm|max=n' => \$maxCount
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

open(IN, "<$input") || die "Unable to read $input: $!";
open(OUT, ">$output") || die "Unable to write to $output: $!";

sub ReadInFasta {
	print "Reading fasta\n";
	# Set the variable for the fasta input file
	my $fastaInput = shift;
	# Setup fasta hash to return at the end
	while (my $line = <$fastaInput>) {
		if ($line =~ /\>/ && $FastaErrorCounter == 0) {
			# print "Made it to ID!\n";
			chomp $line;
			$FastaID = $line;
			# Get rid of the arrow from the ID
			$FastaID =~ s/\>//;
			$FastaID =~ s/_//g;
			$FastaErrorCounter = 1;
		} elsif ($line =~ /\>/ && $FastaErrorCounter == 1) {
			die "KILLED BY BAD FASTA! There was no sequence before ID $line: $!";
		} elsif ($line !~ /\>/ && $FastaErrorCounter == 0) {
			print STDERR "Yikes, is this in block format? That is totally not allowed!\n";
			die "KILLED BY BAD FASTA! There was a missing ID: $!";
		} elsif ($line !~ /\>/ && $FastaErrorCounter == 1) {
			chomp $line;
			# Change out the lower case letters so they match the codon hash
			$line =~ s/g/G/g;
			$line =~ s/a/A/g;
			$line =~ s/c/C/g;
			$line =~ s/t/T/g;
			$FastaHash{$FastaID} = $line;
			$FastaErrorCounter = 0;
		}
	}
	return %FastaHash;
}

sub SlideForKmerSpectrum {
	print "Running sliding window\n";
	my $fastaHash = shift;
	while (my ($fastaKey, $fastaSeq) = each(%{$fastaHash})) {
		my $TestVal = 0;
		# Reset counter and hash
		undef %windowHash;
		$windowValue = 0;
		$length = length($fastaSeq) - 1;
		$correct = $length - $window;
		foreach my $interation (0 .. $correct) {
			# print "Interating $interation\n";
			$windowValue = substr $fastaSeq, $interation, $window;
			# Count the occurnace of that window sequence
			$windowHash{$windowValue}++;
		}
		# For now print it out
		foreach my $key (sort keys %windowHash){
			$TestVal = %windowHash -> {$key};
			$skip = 1 if ($TestVal >= $maxCount);
			# mean
			$MeanSum = $MeanSum + $TestVal;
			++$MeanCounter;
			# \mean
		}
		next if ($skip == 1);
		foreach my $key (sort keys %windowHash){
			$TestVal = %windowHash -> {$key};
			print OUT "$fastaKey\t$TestVal\t$key\n";
		}
	}
	$Mean = $MeanSum / $MeanCounter;
	print STDERR "Mean unfiltered $window mer frequency is $Mean\n";
}

sub FilterKmerSpectrum {
	print "Filtering the spectrum\n";
	my $KmerReadHash = shift;
	foreach my $KmerKey (sort keys %$KmerReadHash) {
		my $TestVal = 0;
		$skip = 0;
		# Iterate through each hash
		foreach my $IdKey (sort keys %{ $KmerReadHash -> {$KmerKey} }) {
			$TestVal = $KmerReadHash -> {$KmerKey}{$IdKey};
			$skip = 1 if ($TestVal >= $maxCount);
			# mean
			$MeanSum = $MeanSum + $TestVal;
			++$MeanCounter;
			# \mean
		}
		next if ($skip == 1);
		foreach my $IdKey (sort keys %{ $KmerReadHash -> {$KmerKey} }) {
			$TestVal = $KmerReadHash -> {$KmerKey}{$IdKey};
			print OUT "$KmerKey\t$TestVal\t$IdKey\n";
		}
	}
	$Mean = $MeanSum / $MeanCounter;
	print STDERR "Mean unfiltered $window mer frequency is $Mean\n";
}

my %Fasta = ReadInFasta(\*IN);
# print Dumper \%Fasta;
SlideForKmerSpectrum(\%Fasta);

close(IN);
close(OUT);

# Get the toal time to run the script
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDERR "\nCalculated kmer spectrum in $run_time seconds.\n";
