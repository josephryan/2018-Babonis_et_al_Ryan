#!/usr/bin/perl

our $SEED = 420;

use strict;
use warnings;
use List::Util qw(shuffle);
use Data::Dumper;

our $REPS = 10000;
our $SEBE_FILE = 'cell_clusters_sebe-pedros.csv';
our @TEST_IDS = qw(189.id 165.id);

MAIN: {
    srand($SEED);
    foreach my $test_file (@TEST_IDS) {
        my $ra_ids = get_ids($test_file);
        my $rh_sebe = get_data($SEBE_FILE);
        my $num_present_in_sebe = get_num_present_in_sebe($rh_sebe,$ra_ids);
        my $test_stat = get_top_cluster_count($rh_sebe,$ra_ids);
        my $pval = mc($rh_sebe,$num_present_in_sebe,$test_stat,$REPS);
        print_pvalue($test_file,$pval);
    }
}

sub print_pvalue {
    my $file = shift;
    my $pval = shift;
    if ($pval == 0) {
        my $pv = 1/$REPS;
        print "PVAL: <$pv ($file)\n";
    } else {
        print "PVAL: = $pval ($file)\n";
    }
}

sub get_ids {
    my $file = shift;
    my @ids  = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        chomp $line;
        push @ids, $line unless ($line =~ m/^\s*$/);
    }
    return \@ids;
}

sub mc {
    my $rh_sebe = shift;
    my $num     = shift;
    my $test    = shift;
    my $reps    = shift;
    my @ids = keys %{$rh_sebe};
    my $gte = 0;
    for (my $i = 0; $i < $reps; $i++) {
        my @shuf = shuffle(@ids);
        my @randset = @shuf[0..$num];
        my $tcc = get_top_cluster_count($rh_sebe,\@randset);
        $gte++ if ($tcc >= $test);
    } 
    my $pval = $gte / $reps;
    return $pval;
}

sub get_num_present_in_sebe {
    my $rh_d = shift;
    my $ra_test = shift;
    my $count = 0;
    foreach my $id (@{$ra_test}) {
        $count++ if ($rh_d->{$id});
    }
    return $count;
}

sub get_top_cluster_count {
    my $rh_dat = shift;
    my $ra_test = shift;
    my %cl_count = ();
    foreach my $id (@{$ra_test}) {
        foreach my $c_id (@{$rh_dat->{$id}}) {
            $cl_count{$c_id}++; 
        }
    }
    my @sorted_vals = sort {$b <=> $a} values %cl_count;
    my $tcc = $sorted_vals[0];
    return $tcc;
}

sub get_data {
    my $file = shift;
    my %data = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        chomp $line;
        my @f = split /,/, $line;
        push @{$data{$f[1]}},  $f[0] if ($f[1] =~ m/^ML/);
    }
    return \%data;
}
