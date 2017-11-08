#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;
use List::Util;
use Data::Dumper;

our $SEED  = 42;
our $BLAST = 'ML2.2_v_meta_minusML.blastp';
our $FASTA = '/bwdata1/jfryan/00-DATA/ML2.2.aa';
our $REPS  = 10000;

MAIN: {
    srand($SEED);
    my $file        = $ARGV[0] or die "usage: $0 FILE\n";
    my $ra_seq_ids  = get_seq_ids($file);
    my $ra_all_ids  = get_seq_ids($FASTA);
    my $rh_bl       = get_bl($BLAST);

    my $num         = get_num_of_non_orphs($rh_bl,$ra_seq_ids,scalar(@{$ra_seq_ids}));
  
    my $total = scalar(@{$ra_seq_ids});
    my $num_orphs = scalar(@{$ra_seq_ids}) - $num;
    print "out of $total, $num_orphs have no hits to our 11-taxa animal database with E-Vals at or below 0.01\n";

    my $less_orphs = 0;
    for (my $i = 0; $i < $REPS; $i++) {
        my @shuf = List::Util::shuffle(@{$ra_all_ids});
        my $rnum = get_num_of_non_orphs($rh_bl,\@shuf,scalar(@{$ra_seq_ids}));
        $less_orphs++ if ($rnum <= $num);
    }
    my $pval = $less_orphs / $REPS;
    print "p-value: $pval ($less_orphs / $REPS)\n";
}

sub get_num_of_non_orphs {
    my $rh_bl = shift;
    my $ra_s  = shift;
    my $num   = shift;
    my $hits  = 0;
    for (my $i = 0; $i < $num; $i++) {
        $hits++ if ($rh_bl->{$ra_s->[$i]});    
print "$ra_s->[$i]\n" unless ($rh_bl->{$ra_s->[$i]}); # if want to print ids
    }
exit; # if want to print ids
    return $hits;
}

sub get_seq_ids {
    my $fa = shift;
    my $ns = 0;
    my @ids = ();
    my $fp = JFR::Fasta->new($fa);
    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        push @ids, $id;
    }
    return \@ids;
}

sub get_bl {
    my $file = shift;
    my %bl = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        chomp $line;
        my @f = split /\t/, $line;
        $bl{$f[0]} = 1;
    }
    return \%bl;
}
