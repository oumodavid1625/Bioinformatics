#!/usr/bin/perl
use strict;
use warnings;

# Open the FASTA file
my $filename = 'tubby.fasta';
open(my $fh, '<', $filename) or die "Could not open file '$filename' $!";

# Variables for storing sequence and header
my $header = "";
my $sequence = "";

# Read FASTA file
while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ /^>/) {  
        $header = $line;
    } else {
        $sequence .= $line;
    }
}
close $fh;

# Calculate sequence length
my $length = length($sequence);

# Count nucleotide frequencies
my $A = ($sequence =~ tr/Aa//);
my $T = ($sequence =~ tr/Tt//);
my $G = ($sequence =~ tr/Gg//);
my $C = ($sequence =~ tr/Cc//);

# Calculate GC content
my $GC_content = (($G + $C) / $length) * 100;

# Print results
print "Header: $header\n";
print "Sequence Length: $length\n";
print "A: $A\nT: $T\nG: $G\nC: $C\n";
print "GC Content: $GC_content%\n";
