#! /usr/bin/perl
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# Set use
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use List::Util 'shuffle';
use Parallel::Loops;
# Timer
my $start_run = time();

# Set variables
my %MidHash;
my %InHash;
my %OutHash;
my %FastaHash;
my %CounterHash;
my %windowHash;
my %ReferenceHash;
my %referenceHash;
my %BrayCurtisHash;
my %frequency;
my $frequency;
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
my $processors = 1;

# Set the options
GetOptions(
	'h|help' => \$opt_help,
	'i|input=s' => \$input,
	't|test=s' => \$test,
	'o|output=s' => \$output,
	'f|outformat=s' => \$outformat,
	'w|window=n' => \$window,
	'p|processors=n' => \$processors,
	'r|reverse' => \$reverse
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

open(my $IN, "<", "$input") || die "Unable to read $input: $!";
open(my $TEST, "<", "$test") || die "Unable to read $test: $!";
open(my $OUT, ">", "$output") || die "Unable to write to $output: $!";
open(my $OUTFMT, ">", "$outformat") || die "Unable to write to $outformat: $!";

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
	my ($CounterHash) = @_;
	my $TotalCounter = 0;
	while (my ($key, $value) = each(%{$CounterHash})) {
		$TotalCounter = $TotalCounter + $value;
	}
	return $TotalCounter;
}

sub CreateKmerHash {
	undef %ReferenceHash;
	undef %FastaHash;
	my ($fastaHash) = @_;
	while (my ($fastaKey, $fastaSeq) = each(%{$fastaHash})) {
		my $TestVal = 0;
		my %WindowResult = &ReturnSlidingWindow($fastaSeq);
		# Add individual kmer counts to master hash
		while (my ($key, $value) = each(%WindowResult)) {
			$ReferenceHash{$fastaKey}{$key} = $value;
		}
	}
	return %ReferenceHash;
}

my $pl = Parallel::Loops->new($processors);

sub HashRandomSubsample {
	undef %InHash;
	undef %OutHash;
	my ($InHash, $subcount) = @_;
	my $totalHashCount = &GetFrequencyCount(\%{$InHash});
	my @InArray=(0)x$totalHashCount;
	# Convert hash to array
	@InArray = map { ($_) x $InHash -> {$_}} keys \%{$InHash};
	@InArray = shuffle(@InArray);
	@InArray = splice(@InArray, 0, $subcount);
	foreach my $interation (@InArray) {
		# Count the occurnace of that window sequence
		$OutHash{$interation}++;
	}
	# print Dumper \%OutHash;
	return %OutHash;
}

sub BcDistanceFromReference {
	my ($queryHash, $referenceHash) = @_;
	my $ProgressCounter = 1;
	my $TestKeyCount = keys %{$queryHash};
	my $RefKeysCount = keys %{$referenceHash};
	# Get distance for each query sequence
	while (my ($fastaKey, $fastaSeq) = each(%{$queryHash})) {
		# print "Fasta key is $fastaKey\n";
		my $progress = 100 * $ProgressCounter / $TestKeyCount;
		my $refprogcounter = 1;
		print STDERR "\rProgress: $progress \%";
		++$ProgressCounter;
		# Get kmer frequency for this sequence
		undef %WindowResult;
		%WindowResult = &ReturnSlidingWindow($fastaSeq);
		# print Dumper \%WindowResult;
		# Calculate distance from each reference
		my @arraysub;
		foreach my $linekey (sort keys %{$referenceHash}) {
			# print "$linekey\n";
			push @arraysub, $linekey;
		}
		undef %MidHash;
		# print Dumper \%MidHash;
		$pl->share( \%MidHash );
		$pl->foreach( \@arraysub, sub {
			my $referenceID = $_;
			# print "ID is $referenceID\n";
			my $refProgress = 100 * $refprogcounter / $RefKeysCount;
			# print STDERR "\rProgress: $progress % -- $refProgress %";
			++$refprogcounter;
			my $BrayCurtis = 0;
			my $TotalReferenceCount = 0;
			my $TotalCount = 0;
			undef %frequency;
			undef %frequencyint;
			undef %NewWindowResult;

			# Get the count for the query
			# print Dumper \%WindowResult;
			$TotalCount = &GetFrequencyCount(\%WindowResult);
			# print "Total count is $TotalCount\n";

			# Make hash frequency which will contain subsampled kmer
			# frequency of the current loop reference.
			while (my ($key, $value) = each( $referenceHash -> {$referenceID}) ) {
				$frequencyint{$key} = $value;
			}
			# Subample the reference freq hash
			$TotalReferenceCount = &GetFrequencyCount(\%frequencyint);
			# print "Total reference count $TotalReferenceCount\n";
			if ($TotalCount < $TotalReferenceCount) {
				# print "Problem 1\n";
				%frequency = &HashRandomSubsample(\%frequencyint, $TotalCount);
				%NewWindowResult = %WindowResult;
			} elsif ($TotalCount > $TotalReferenceCount) {
				# print "Problem 2\n";
				%NewWindowResult = &HashRandomSubsample(\%WindowResult, $TotalReferenceCount);
				%frequency = %frequencyint;
			} elsif ($TotalCount == $TotalReferenceCount) {
				# print "Problem 3\n";
				%frequency = %frequencyint;
				%NewWindowResult = %WindowResult;
			}
			# print Dumper \%frequency;
			$TotalReferenceCount = &GetFrequencyCount(\%frequency);
			# print "New total reference count $TotalReferenceCount\n";
			my $NewQueryCount = &GetFrequencyCount(\%NewWindowResult);
			# Make sure the subsampling has both hashes equal
			die "Subsampling did not provide equal frequency counts: $!" unless ($NewQueryCount == $TotalReferenceCount);

			# Sum less of shared kmer sequences
			my $LesserValueSum = 0;
			while (my ($KmerSeq, $KmerCount) = each(\%NewWindowResult)) {
				# print "Kmer is $KmerSeq\n";
				$ReferenceCount = 0;
				my $frequency = \%frequency;
				# Here we are iterating through query kmers
				next unless (exists $frequency -> {$KmerSeq});
				my $ReferenceCount = $frequency -> {$KmerSeq};
				# Add on lesser of the two shared counts
				if ($ReferenceCount >= $KmerCount) {
					$LesserValueSum = $LesserValueSum + $KmerCount;
				} elsif ($ReferenceCount < $KmerCount) {
					$LesserValueSum = $LesserValueSum + $ReferenceCount;
				}
				# print "Lesser value sum is $LesserValueSum\n";
			}
			# Calculate Bray Curtis Distance
			# print "Final total count is $TotalCount\n";
			# print "Final total reference count $TotalReferenceCount\n";
			unless ($TotalCount == 0 || $TotalReferenceCount == 0) {
				$BrayCurtis = 1 - ( (2 * $LesserValueSum) / ($TotalCount + $TotalReferenceCount) );
				# print "Resulting BC: $BrayCurtis\n";
				# Save the result into a hash
				$MidHash{$referenceID} = $BrayCurtis;
			}
			# print "Outside BC: $BrayCurtis\n";
			# print Dumper \%BrayCurtisHash;
		});
		while (my ($key, $value) = each( %MidHash ) ) {
				$BrayCurtisHash{$fastaKey}{$key} = $value;
		}
		# print "Fasta Key is $fastaKey\n";
		# print Dumper \%MidHash;
		# print Dumper \%BrayCurtisHash;
	}
	# print Dumper \%BrayCurtisHash;
	return %BrayCurtisHash;
}

sub IdentityScores {
	my $ResultHash = shift @_;
	# Return the results for each query sequence
	foreach my $queryID (sort keys %{$ResultHash}) {
		my $counter = 0;
		print $OUT "Query: $queryID\n";
		foreach my $BCvalue (sort keys %{$ResultHash -> {$queryID}}) {
			$HitID = $ResultHash -> {$queryID}{$BCvalue};
			print $OUTFMT "$queryID\t$HitID\t$BCvalue\n";
		}
		# Get the top 5 distances
		foreach my $BCvalue (sort { $ResultHash -> {$queryID} -> {$a} <=> $ResultHash -> {$queryID} -> {$b} } keys $ResultHash -> {$queryID}) {
			last if ($counter == 6);
			$HitID = $ResultHash -> {$queryID}{$BCvalue};
			print $OUT "\t$BCvalue: $HitID\n";
			$counter++;
		}
	}
}

# Now run the subroutines defined above

# Remind the users that the reverse compliment was requested
print STDERR "As requested, running additional reverse compliment.\n" if ($reverse);
# Read in the querry fasta
my %Fasta = &ReadInFasta(\*$IN);
# Read in the reference fasta
my %TestFasta = &ReadInFasta(\*$TEST);
print "Creating reference hash.\n";
# Create a kmer hash from the reference file
my %referencekmer = &CreateKmerHash(\%Fasta);
undef %Fasta;
# Run Bray Curtis Distance
print "Getting Distances.\n";
my %Results = &BcDistanceFromReference(\%TestFasta, \%referencekmer);
# Format the scores and print to output
&IdentityScores(\%Results);

# Close out
close($IN);
close($OUT);
close($TEST);
close($OUTFMT);

# Get the toal time to run the script
my $end_run = time();
my $run_time = $end_run - $start_run;
print STDOUT "\nRUNTIME\t$window\t$run_time\n";

# Add a help menu for the user

=head1 NAME

FilterFasta.pl

=head1 SYNOPSIS

Maintained by Geoffrey Hannigan (ghannig@umich.edu)

Given two fasta files, this will calculate and compare the kmer frequency profiles of each sequence to every other sequence.

=head1 OPTIONS

CalculateKmerDistancesPar.pl 
	-i <input fasta 1>
	-t <input fasta 2>
	-o <output file>
	-f <formatted output file>
	-w <kmer window size>
	-p <processors>
	-r <reverse compliment option>
	-h <print helpful help menu> 

=head1 ARGUMENTS

-h | --help	Print this helpful help menu.

-i | --input	Input fasta file for kmer comparison. This should be the smaller file.

-t | --test Second fasta file for kmer comparison. This should be the larger file.

-o | --output	Output file name. This will be a table showing the top five closest kmer profiles for each test sequence.

-f | --outformat	Output file name. This will be a table including the distance for every comparison.

-w | --window	The kmer length [number] to use for the analysis. Default is four (tetramer).

-p | --processors	The [number] of processors to use for the analysis.

-r | --reverse	Add this flag if you would like to include the reverse compliments of each sequence in the analysis.

=cut
