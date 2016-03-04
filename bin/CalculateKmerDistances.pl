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
my %InHash;
my %FastaHash;
my %CounterHash;
my %windowHash;
my %ReferenceHash;
my %referenceHash;
my %BrayCurtisHash;
my %frequency;
my %frequencyint;
my %WindowResult;
my %NewWindowResult;
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
my $reverse = '';
my $key;

# Set the options
GetOptions(
	'h|help' => \$opt_help,
	'i|input=s' => \$input,
	't|test=s' => \$test,
	'o|output=s' => \$output,
	'f|outformat=s' => \$outformat,
	'w|window=n' => \$window,
	'r|reverse' => \$reverse
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

sub ReverseCompliment {
	# Take in a string of nucleotides
	my $fastaSeq = '';
	$fastaSeq = shift;
	# Get the compliment sequences of the gene
	$fastaSeq =~ tr/ACGT/TGCA/;
	# Reverse to finish getting reverse compliment
	$fastaSeq = reverse $fastaSeq;
	return $fastaSeq;
}

sub ReturnSlidingWindow {
	# Get fasta sequence string as variable
	my $fastaSeq = '';
	$fastaSeq = shift;
	my $reverseSeq = '';
	$reverseSeq = &ReverseCompliment($fastaSeq);
	undef %windowHash;
	$windowValue = 0;
	$length = length($fastaSeq) - 1;
	$correct = $length - $window;
	foreach my $interation (0 .. $correct) {
		# print "Interating $interation\n";
		$windowValue = substr $fastaSeq, $interation, $window;
		# Count the occurnace of that window sequence
		$windowHash{$windowValue}++;
		# Also run reverse compliment if specified
		next unless ($reverse);
		$windowValue = 0;
		$windowValue = substr $reverseSeq, $interation, $window;
		$windowHash{$windowValue}++;
	}
	# Return frequencies as hash
	return %windowHash;
}

sub GetFrequencyCount {
	undef %CounterHash;
	my $CounterHash = shift;
	my $TotalCounter = 0;
	while (my ($key, $value) = each(%{$CounterHash})) {
		$TotalCounter = $TotalCounter + $value;
	}
	return $TotalCounter;
}

sub CreateKmerHash {
	undef %ReferenceHash;
	undef %FastaHash;
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

sub HashRandomSubsample {
	undef %InHash;
	my ($InHash, $subcount) = @_;
	# print Dumper \%{$InHash};
	# print STDERR "Subsampling to depth $subcount\n";
	my $totalHashCount = &GetFrequencyCount(\%{$InHash});
	# print STDERR "input hash count is $totalHashCount\n";
	# Input is the resulting subsample, so calculate how much needs
	# to be randomly subtracted from the hash.
	my $correctCount = $totalHashCount - $subcount;
	# print STDERR "corrected count is $correctCount\n";
	# print STDERR "corrected count is $correctCount\n";
	foreach my $subiter (1..$correctCount) {
		# print STDERR "iteration is $subiter\n";
		$key = (keys \%{$InHash})[rand keys \%{$InHash}];
		# print STDERR "KEY IS $key\n";
		my $valueint = $InHash -> {$key};
		# print STDERR "VALUE IS $valueint\n";
		$InHash -> {$key}--;
		my $valueafter = $InHash -> {$key};
		$valueafter = $valueafter + 1;
		# print STDERR "VALUE IS $valueint\n";
		# Some quality control
		die "The subsampling is not subtracting properly: $!" unless ($valueint == $valueafter);
		die "Subsampling is giving negative counts and that makes no sense: $!" if ($valueint < 0);
		delete($InHash -> {$key}) if ($InHash -> {$key} == 0);
	}
	# print Dumper \%{$InHash};
	# print "$InHash\n";
	return %{$InHash};
}

sub BcDistanceFromReference {
	my ($queryHash, $referenceHash) = @_;
	my $ProgressCounter = 1;
	my $TestKeyCount = keys %{$queryHash};
	# Get distance for each query sequence
	while (my ($fastaKey, $fastaSeq) = each(%{$queryHash})) {
		# print STDERR "Reading sample $fastaKey\n";
		my $progress = 100 * $ProgressCounter / $TestKeyCount;
		print STDERR "\rProgress: $progress \%";
		++$ProgressCounter;
		# Get kmer frequency for this sequence
		undef %WindowResult;
		%WindowResult = &ReturnSlidingWindow($fastaSeq);
		# print Dumper \%WindowResult;
		# Calculate distance from each reference
		foreach my $referenceID (sort keys %{$referenceHash}) {
			# print STDERR "Referencing $referenceID\n";
			my $BrayCurtis = 0;
			my $TotalReferenceCount = 0;
			my $TotalCount = 0;
			undef %frequency;
			undef %frequencyint;
			undef %NewWindowResult;

			# Get the count for the query
			$TotalCount = &GetFrequencyCount(\%WindowResult);
			# print "Total query frequency is $TotalCount\n";

			# Make hash frequency which will contain subsampled kmer
			# frequency of the current loop reference.
			while (my ($key, $value) = each( $referenceHash -> {$referenceID}) ) {
				$frequencyint{$key} = $value;
			}
			# print STDERR "FREQINT\n";
			# print Dumper \%frequencyint;
			# Subample the reference freq hash
			$TotalReferenceCount = &GetFrequencyCount(\%frequencyint);
			# print "Total reference frequency is $TotalReferenceCount\n";
			# print Dumper \%frequencyint;
			if ($TotalCount < $TotalReferenceCount) {
				# print "Reference is greater.\n";
				%frequency = &HashRandomSubsample(\%frequencyint, $TotalCount);
				%NewWindowResult = %WindowResult;
			} elsif ($TotalCount > $TotalReferenceCount) {
				# print "Query is greater.\n";
				%NewWindowResult = &HashRandomSubsample(\%WindowResult, $TotalReferenceCount);
				%frequency = %frequencyint;
			} elsif ($TotalCount == $TotalReferenceCount) {
				# print "Reference same length as query.\n";
				%frequency = %frequencyint;
				%NewWindowResult = %WindowResult;
			}
			# print Dumper \%frequency;
			$TotalReferenceCount = &GetFrequencyCount(\%frequency);
			my $NewQueryCount = &GetFrequencyCount(\%NewWindowResult);
			# print STDERR "RESULT: Query = $NewQueryCount\nReference = $TotalReferenceCount\n";
			# Make sure the subsampling has both hashes equal
			die "Subsampling did not provide equal frequency counts: $!" unless ($NewQueryCount == $TotalReferenceCount);

			# Sum less of shared kmer sequences
			my $LesserValueSum = 0;
			while (my ($KmerSeq, $KmerCount) = each(\%NewWindowResult)) {
				$ReferenceCount = 0;
				# Here we are iterating through query kmers
				next unless (exists my $frequency -> {$KmerSeq});
				my $ReferenceCount = $frequency -> {$KmerSeq};
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
			$HitID = $ResultHash -> {$queryID}{$BCvalue};
			print OUTFMT "$queryID\t$HitID\t$BCvalue\n";
		}
		# Get the top 5 distances
		foreach my $BCvalue (sort { $ResultHash -> {$queryID} -> {$a} <=> $ResultHash -> {$queryID} -> {$b} } keys $ResultHash -> {$queryID}) {
			last if ($counter == 6);
			$HitID = $ResultHash -> {$queryID}{$BCvalue};
			print OUT "\t$BCvalue: $HitID\n";
			$counter++;
		}
	}
}

print STDERR "As requested, running additional reverse compliment.\n" if ($reverse);
my %Fasta = &ReadInFasta(\*IN);
my %TestFasta = &ReadInFasta(\*TEST);
print "Creating reference hash.\n";
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


