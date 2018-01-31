#!/usr/bin/perl

# this script takes in a target file and a background file
# and compares the level of low complexity and composition
# to random datasets of the same size from the background file

use strict;
use warnings;
use JFR::Fasta;
use File::Temp;
use List::Util qw/shuffle/;
use Data::Dumper;

our $VERSION = 0.01;
our $SRAND   = 420;
our $REPS    = 1000;
our $SEG     = 'segmasker';

MAIN: {
    srand($SRAND);
    my $targ_fa = $ARGV[0] or die "usage: $0 TARGET_FA BACKGROUND_FA\n";
    my $bg_fa   = $ARGV[1] or die "usage: $0 TARGET_FA BACKGROUND_FA\n";
    my $ra_targ = get_seqs_from_fasta($targ_fa);
    my $ra_bg   = get_seqs_from_fasta($bg_fa);
    my $targ_complexity = get_complexity($ra_targ);
    my $rh_targ_composition = get_composition($ra_targ);
    my $num_per_set = scalar(@{$ra_targ});
    my $gte_complexity = 0;
    my %gte_composition = ();
    for (1 .. $REPS) {
        my $ra_rand = get_rand($ra_bg,$num_per_set);
        my $rcpx    = get_complexity($ra_rand);
        $gte_complexity++ if ($rcpx >= $targ_complexity);
        my $rh_comp = get_composition($ra_rand);
        compare_composition($rh_targ_composition,$rh_comp,\%gte_composition);
    }
    my $rh_pvals = convert_to_pvals(\%gte_composition,$REPS);
    print "compcompmcmc.pl version $VERSION\n";
    print "TARGET FASTA: $targ_fa\n";
    print "BACKGROUND FASTA: $bg_fa\n\n";
    print "Amino acid composition:\n  low p-values are over represented; ";
    print "high p-values are under represented\n";
    foreach my $aa (sort {$rh_pvals->{$a} <=> $rh_pvals->{$b}} keys %{$rh_pvals}) {
        print "  $aa: $rh_pvals->{$aa}\n";
    }
    my $complexity_pval = $gte_complexity / $REPS;
    print "\nComplexity p-value: $complexity_pval\n";
    print "$gte_complexity of our $REPS random sets had complexity higher than sequences in $targ_fa\n";
}

sub convert_to_pvals {
    my $rh_comp = shift;
    my $reps = shift;
    my %pvals = ();
    foreach my $aa (keys %{$rh_comp}) {
        $pvals{$aa} = $rh_comp->{$aa} / $reps;
    }
    return \%pvals;
}
sub compare_composition {
    my $rh_t = shift;
    my $rh_r = shift;
    my $rh_c = shift;
    foreach my $aa (keys %{$rh_t}) {
        next if ($aa eq 'X');
        $rh_c->{$aa}++ if ($rh_r->{$aa} >= $rh_t->{$aa});
    }
}

sub get_rand {
    my $ra_bg = shift;
    my $num   = shift;
    my @rand  = ();
    my @shuf  = List::Util::shuffle(@{$ra_bg});
    for (my $i = 0; $i < $num; $i++) {
        push @rand, $shuf[$i]; 
    }
    return \@rand;
}

sub get_complexity {
    my $ra_seqs = shift;
    my ($fh, $file) = tmpnam();    
    my $count = 0;
    foreach my $seq (@{$ra_seqs}) {
        $count++;
        print $fh ">$count\n$seq\n";
    }
    my @out = `$SEG -outfmt interval -in $file`;
    my $num_aas = get_num_aas($ra_seqs); 
    my $num_seqs = 0;
    my $num_lc   = 0;
    foreach my $line (@out) {
        next if ($line =~ m/^\s*$/);
        if ($line =~ m/^>/) {
            $num_seqs++;
            next;
        } else {
            $line =~ m/^(\d+)\s+-\s+(\d+)/ or die "unexpected:$line";
            $num_lc += ($2 - $1);
        }
    }
    my $ratio = $num_lc / $num_aas;
    return $ratio;
}

sub get_num_aas {
    my $ra_seqs = shift;
    my $num = 0;
    foreach my $seq (@{$ra_seqs}) {
        $num += length($seq);
    }
    return $num;
}

sub get_composition {
    my $ra_seqs = shift;
    my %counts  = ();
    my %comp    = ();
    my $sum     = 0;
    foreach my $seq (@{$ra_seqs}) {
        my @aas = split /|/, $seq;
        $sum += scalar(@aas);
        foreach my $aa (@aas) {
            $counts{$aa}++;
        }
    }
    foreach my $key (sort keys %counts) {
        $comp{$key} = $counts{$key} / $sum;
    }
    return \%comp;
}

sub get_seqs_from_fasta {
    my $file = shift;
    my @seqs = ();
    my $fp   = JFR::Fasta->new($file);
    while (my $rec = $fp->get_record()) {
        push @seqs, $rec->{'seq'};
    }
    return \@seqs;
}

sub composition {

}

