#!/usr/bin/perl
use strict;
use warnings;

sub calculate_gc_content {
    my ($sequence) = @_;
    my $gc_count = ($sequence =~ tr/GC//);
    my $total_length = length($sequence);
    return $total_length > 0 ? ($gc_count / $total_length) * 100 : 0.0;
}

sub parse_fasta {
    my ($file) = @_;
    my %sequences;
    my $current_id;
    my $current_sequence = '';
    open my $fh, '<', $file or die "Could not open file '$file': $!";
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/^\s+|\s+$//g; # Trim whitespace
        if ($line =~ /^>/) {
            if (defined $current_id) {
                $sequences{$current_id} = $current_sequence;
                $current_sequence = '';
            }
            $current_id = substr($line, 1);
            $current_id =~ s/^\s+|\s+$//g; # Trim whitespace
        } else {
            $current_sequence .= $line if $line;
        }
    }
    close $fh;
    if (defined $current_id) {
        $sequences{$current_id} = $current_sequence;
    }
    return %sequences;
}

sub main {
    my $file = 'rosalind.txt';
    my %sequences = parse_fasta($file);
    unless (%sequences) {
        die "Error: No valid sequences found in the input file '$file'.\n";
    }
    my ($max_gc_id, $max_gc) = (undef, -1.0);
    while (my ($seq_id, $sequence) = each %sequences) {
        my $gc = calculate_gc_content($sequence);
        if ($gc > $max_gc) {
            $max_gc = $gc;
            $max_gc_id = $seq_id;
        }
    }
    if (defined $max_gc_id) {
        printf "%s\n%.6f\n", $max_gc_id, $max_gc;
    } else {
        die "Error: Could not determine the sequence with the highest GC-content.\n";
    }
}

main();
