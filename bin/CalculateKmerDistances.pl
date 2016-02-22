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
my %windowHash;
my %ReferenceHash;
my %referenceHash;
my %BrayCurtisHash;
my %frequency;
my %WindowResult;
my $windowValue;
my $opt_help;
my $input;
my $test;
my $window = 4;
my $output;
my $FastaErrorCounter = 0;
my $length = 0;
my $correct = 0;
my $skip = 0;
my $FastaID;
my $ReferenceCount;
my $referenceHash;
my $HitID;
my $outformat;

# Set the options
GetOptions(
	'h|help' => \$opt_help,
	'i|input=s' => \$input,
	't|test=s' => \$test,
	'o|output=s' => \$output,
	'f|outformat=s' => \$outformat,
	'w|window=n' => \$window
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

open(IN, "<$input") || die "Unable to read $input: $!";
open(TEST, "<$test") || die "Unable to read $test: $!";
open(OUT, ">$output") || die "Unable to write to $output: $!";
open(OUTFMT, ">$outformat") || die "Unable to write to $outformat: $!";

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

sub ReturnSlidingWindow {
	# Get fasta sequence string as variable
	my $fastaSeq = shift;
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
	# Return frequencies as hash
	return %windowHash;
}

sub GetFrequencyCount {
	my $CounterHash = shift;
	my $TotalCounter = 0;
	while (my ($key, $value) = each(%{$CounterHash})) {
		$TotalCounter = $TotalCounter + $value;
	}
	return $TotalCounter;
}

sub CreateKmerHash {
	my $fastaHash = shift;
	while (my ($fastaKey, $fastaSeq) = each(%{$fastaHash})) {
		my $TestVal = 0;
		my %WindowResult = &ReturnSlidingWindow($fastaSeq);
		# my $TotalCount = &GetFrequencyCount(\%WindowResult);
		# Add individual kmer counts to master hash
		while (my ($key, $value) = each(%WindowResult)) {
			$ReferenceHash{$fastaKey}{$key} = $value;
		}
	}
	return %ReferenceHash;
}

sub BcDistanceFromReference {
	my ($queryHash, $referenceHash) = @_;
	my $ProgressCounter = 1;
	my $TestKeyCount = keys %{$queryHash};
	# Get distance for each query sequence
	while (my ($fastaKey, $fastaSeq) = each(%{$queryHash})) {
		my $progress = 100 * $ProgressCounter / $TestKeyCount;
		print STDERR "\rProgress: $progress \%";
		++$ProgressCounter;
		my $TotalReferenceCount = 0;
		my $TotalCount = 0;
		# Get kmer frequency for this sequence
		undef %WindowResult;
		%WindowResult = &ReturnSlidingWindow($fastaSeq);
		$TotalCount = &GetFrequencyCount(\%WindowResult);
		# Calculate distance from each reference
		foreach my $referenceID (sort keys %{$referenceHash}) {
			my $BrayCurtis = 0;
			undef %frequency;
			while (my ($key, $value) = each(%{$referenceHash} -> {$referenceID})) {
				$frequency{$key} = $value;
			}
			$TotalReferenceCount = &GetFrequencyCount(\%frequency);
			# print "Ref Count is $TotalReferenceCount\n";
			# print "Total count is $TotalCount\n";
			# Sum less of shared kmer sequences
			my $LesserValueSum = 0;
			while (my ($KmerSeq, $KmerCount) = each(\%WindowResult)) {
				$ReferenceCount = 0;
				# Here we are iterating through query kmers
				next unless (exists %{$referenceHash} -> {$referenceID}{$KmerSeq});
				my $ReferenceCount = %{$referenceHash} -> {$referenceID}{$KmerSeq};
				# print "Reference is $ReferenceCount\n";
				# print "Kmer is $KmerCount\n";
				# Add on lesser of the two shared counts
				if ($ReferenceCount >= $KmerCount) {
					$LesserValueSum = $LesserValueSum + $KmerCount;
				} elsif ($ReferenceCount < $KmerCount) {
					$LesserValueSum = $LesserValueSum + $ReferenceCount;
				}
			}
			# Calculate Bray Curtis Distance
			next if ($TotalCount == 0 || $TotalReferenceCount == 0);
			$BrayCurtis = 1 - ( (2 * $LesserValueSum) / ($TotalCount + $TotalReferenceCount) );
			# print "$fastaKey\t$referenceID\t$BrayCurtis\n";
			# print "$LesserValueSum\t$TotalCount\t$TotalReferenceCount\n";
			# Save the result into a hash
			$BrayCurtisHash{$fastaKey}{$referenceID} = $BrayCurtis;
		}
	}
	return %BrayCurtisHash;
}

sub IdentityScores {
	my $ResultHash = shift @_;
	# Return the results for each query sequence
	foreach my $queryID (sort keys %{$ResultHash}) {
		my $counter = 0;
		print OUT "Query: $queryID\n";
		foreach my $BCvalue (sort keys %{$ResultHash -> {$queryID}}) {
			$HitID = %{$ResultHash} -> {$queryID}{$BCvalue};
			print OUTFMT "$queryID\t$HitID\t$BCvalue\n";
		}
		# Get the top 5 distances
		foreach my $BCvalue (sort { %{$ResultHash -> {$queryID}} -> {$a} <=> %{$ResultHash -> {$queryID}} -> {$b} } keys %{$ResultHash -> {$queryID}}) {
			last if ($counter == 6);
			$HitID = %{$ResultHash} -> {$queryID}{$BCvalue};
			print OUT "\t$BCvalue: $HitID\n";
			$counter++;
		}
	}
}

my %Fasta = &ReadInFasta(\*IN);
my %TestFasta = &ReadInFasta(\*TEST);
my %referencekmer = &CreateKmerHash(\%Fasta);
my %Results = &BcDistanceFromReference(\%TestFasta, \%referencekmer);
# print Dumper \%Results;
&IdentityScores(\%Results);

close(IN);
close(OUT);
close(TEST);

# Get the toal time to run the script
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDERR "\nCalculated kmer distances in $run_time seconds.\n";


