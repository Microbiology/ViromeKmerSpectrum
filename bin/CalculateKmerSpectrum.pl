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

# Set the options
GetOptions(
	'h|help' => \$opt_help,
	'i|input=s' => \$input,
	'o|output=s' => \$output
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

open(IN, "<$input") || die "Unable to read $input: $!";
open(OUT, ">$output") || die "Unable to write to $output: $!";

sub ReadInFasta {
	print "Made it to fasta\n";
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
			print "Added $FastaID\n";
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

sub SlideForContigSpectrum {
	my $fastaHash = shift;
	print "Reading $fastaHash\n";
	while (my ($fastaKey, $fastaSeq) = each(%{$fastaHash})) {
		# Reset counter and hash
		undef %windowHash;
		$windowValue = 0;
		print "The key is $fastaKey\n";
		$length = length($fastaSeq) - 1;
		$correct = $length - $window;
		foreach my $interation (0 .. $correct) {
			# print "Interating $interation\n";
			$windowValue = substr $fastaSeq, $interation, $window;
			# Count the occurnace of that window sequence
			$windowHash{$windowValue}++;
		}
		# For now print it out
		foreach $key (sort keys %windowHash){
			print OUT "$key\t$windowHash{$key}\t$fastaKey\n";
		}
	}
}

my %Fasta = ReadInFasta(\*IN);
# print Dumper \%Fasta;
SlideForContigSpectrum(\%Fasta);


