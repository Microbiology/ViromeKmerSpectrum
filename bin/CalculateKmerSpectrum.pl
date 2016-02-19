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
my %TestHash;
my %testHash;
my %FinalHash;
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
my $maxCount = 0;
my $minCount = 0;
my $Mean;
my $MeanCounter = 0;
my $MeanSum = 0;
my $progress;
my $test;
my $testOuput;

# Set the options
GetOptions(
	'h|help' => \$opt_help,
	'i|input=s' => \$input,
	't|test=s' => \$test,
	'o|output=s' => \$output,
	'r|testOuput=s' => \$testOuput,
	'w|window=n' => \$window,
	'm|max=n' => \$maxCount,
	'q|min=n' => \$minCount
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

open(IN, "<$input") || die "Unable to read $input: $!";
open(TEST, "<$test") || die "Unable to read $test: $!";
open(OUT, ">$output") || die "Unable to write to $output: $!";
open(TOUT, ">$testOuput") || die "Unable to write to $testOuput: $!";

sub ReadInFasta {
	print "Reading fasta\n";
	# Ensure the hash is empty before use
	undef %FastaHash;
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
	my ($fastaHash, $testHash) = @_;
	my $TestKeyCount = keys %{$testHash};
	my $TrainKeyCount = keys %{$fastaHash};
	my $ProgressCounter = 1;
	# Get hash of kmer patterns within test
	while (my ($fastaKey, $fastaSeq) = each(%{$testHash})) {
		# Use 50 because this is half of the percent progress
		$progress = 100 * $ProgressCounter / $TestKeyCount;
		print STDERR "\rFilter progress (1/3): $progress \%";
		++$ProgressCounter;
		$skip = 0;
		undef %TestHash;
		my $TestVal = 0;
		$windowValue = 0;
		$length = length($fastaSeq) - 1;
		$correct = $length - $window;
		foreach my $interation (0 .. $correct) {
			# print "Interating $interation\n";
			$windowValue = substr $fastaSeq, $interation, $window;
			# Count the occurnace of that window sequence
			# print "Window value is $windowValue\n";
			$TestHash{$windowValue}++;
		}
		foreach my $key (sort keys %TestHash){
			$TestVal = %TestHash -> {$key};
			$skip = 1 if ($TestVal >= $maxCount && $maxCount != 0 && $minCount == 0);
			$skip = 1 if ($TestVal < $minCount && $minCount != 0 && $maxCount == 0);
		}
		next if ($skip == 1);
		%FinalHash = (%FinalHash, %TestHash);
	}
	undef %TestHash;
	$ProgressCounter = 1;
	# Print the test
	while (my ($fastaKey, $fastaSeq) = each(%{$testHash})) {
		$progress = 100 * $ProgressCounter / $TestKeyCount ;
		print STDERR "\rWrite test progress (2/3): $progress \%";
		++$ProgressCounter;
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
		# Only get the kmers present in the test set
		foreach my $key (sort keys %FinalHash){
			unless (exists %windowHash -> {$key}){
				print TOUT "$key\t0\t$fastaKey\n";
				next;
			}
			$TestVal = %windowHash -> {$key};
			# mean
			$MeanSum = $MeanSum + $TestVal;
			++$MeanCounter;
			# \mean
			print TOUT "$key\t$TestVal\t$fastaKey\n";
		}
	}
	$Mean = $MeanSum / $MeanCounter unless ($MeanCounter == 0);
	print STDERR "\nMean unfiltered $window mer frequency in test set is $Mean\n";
	undef %testHash;
	$ProgressCounter = 1;
	# Print the training set
	while (my ($fastaKey, $fastaSeq) = each(%{$fastaHash})) {
		$progress = 100 * $ProgressCounter / $TrainKeyCount ;
		print STDERR "\rWrite train progress (3/3): $progress \%";
		++$ProgressCounter;
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
		# Only get the kmers present in the test set
		foreach my $key (sort keys %FinalHash){
			unless (exists %windowHash -> {$key}){
				print OUT "$key\t0\t$fastaKey\n";
				next;
			}
			$TestVal = %windowHash -> {$key};
			# mean
			$MeanSum = $MeanSum + $TestVal;
			++$MeanCounter;
			# \mean
			print OUT "$key\t$TestVal\t$fastaKey\n";
		}
	}
	$Mean = $MeanSum / $MeanCounter unless ($MeanCounter == 0);
	print STDERR "\nMean unfiltered $window mer frequency is $Mean\n";
}

my %Fasta = ReadInFasta(\*IN);
my %FilterFasta = ReadInFasta(\*TEST);
# print Dumper \%FilterFasta;
SlideForKmerSpectrum(\%Fasta, \%FilterFasta);

close(IN);
close(OUT);
close(TEST);
close(TOUT);

# Get the toal time to run the script
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDERR "\nCalculated kmer spectrum in $run_time seconds.\n";


