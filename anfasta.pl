#!/usr/bin/perl
use strict;
use warnings;

# ===== Subroutines =====

# Reverse complement
sub revcomp {
    my $seq = shift;
    $seq = reverse($seq);
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}

# GC content
sub gc_content {
    my $seq = shift;
    my $gc = ($seq =~ tr/GCgc//);
    return ($gc / length($seq)) * 100;
}

# Transcribe DNA to RNA
sub transcribe {
    my $seq = shift;
    $seq =~ tr/T/U/;
    return $seq;
}

# Print codons
sub print_codons {
    my $seq = shift;
    print "Codons:\n";
    for (my $i = 0; $i <= length($seq) - 3; $i += 3) {
        my $codon = substr($seq, $i, 3);
        print "$codon ";
    }
    print "\n";
}

# ===== Main Script =====

my $fasta_file = "tubby_mrna.fasta";
open(my $fh, '<', $fasta_file) or die "Cannot open file '$fasta_file': $!";

my ($header, $sequence) = ("", "");

while (<$fh>) {
    chomp;
    if (/^>(.*)/) {
        # Process the previous sequence if any
        if ($sequence ne "") {
            process_sequence($header, $sequence);
            $sequence = "";
        }
        $header = $1;
    } else {
        $sequence .= $_;
    }
}
# Process the last sequence
process_sequence($header, $sequence) if $sequence ne "";

close($fh);

# ===== Sequence Processor =====

sub process_sequence {
    my ($hdr, $seq) = @_;
    print "\n>$hdr\n";
    print "Original Sequence:\n$seq\n";

    my $rev = revcomp($seq);
    print "Reverse Complement:\n$rev\n";

    my $gc = gc_content($seq);
    printf "GC Content: %.2f%%\n", $gc;

    my $rna = transcribe($seq);
    print "Transcribed RNA:\n$rna\n";

    print_codons($seq);
}
