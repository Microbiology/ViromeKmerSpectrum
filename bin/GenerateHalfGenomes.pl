#!/usr/bin/perl
# Geoffrey Hannigan
# Pat Schloss Lab
# University of Michigan

# Set use
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Set Variables
my %FastaHash;
my $opt_help;
my $input;
my $output1;
my $output2;
my $FastaErrorCounter;
my $FastaID;


# Set the options
GetOptions(
	'h|help' => \$opt_help,
	'i|input=s' => \$input,
	'o|output1=s' => \$output1,
	't|output2=s' => \$output2
);

pod2usage(-verbose => 1) && exit if defined $opt_help;

open(my $IN, "<", "$input") || die "Unable to read $input: $!";
open(my $ONE, ">", "$output1") || die "Unable to write to $output1: $!";
open(my $TWO, ">", "$output2") || die "Unable to write to $output2: $!";

while (my $line = <$IN>) {
	#If the line is a fasta title line (starts with arrow)
	if ($line =~ /\>/) {
		print $ONE $line;
		print $TWO $line;
	}
	else {
		chomp $line;
		my $length = length $line;
		my $halfway = $length / 2;
		my $genome1 = substr $line, 0, $halfway;
		my $genome2 = substr $line, $halfway, $length;
		print $ONE "$genome1\n";
		print $TWO "$genome2\n";
	}
}

# Print the final newline

#Close out files and print completion note to STDOUT
close($IN);
close($ONE);
close($TWO);
print "Completed\n"

