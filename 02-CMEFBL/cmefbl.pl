#!/usr/bin/perl

use strict;
use warnings;
use List::Util 'shuffle';
use Data::Dumper;

our $SRAND = 420;
srand($SRAND);

our $FILE = 'ML_all_datapoints.tab';
our $LAST_EARLY_TIMEPT = 8.5;
our $FIRST_LATE_TIMEPT = 13;
our $REPS = 10000;

our $FA_DIR = './fasta';

# the following lines were used to test the effect of our 70% cutoff of
# presence in other ctenophores. We found little effect so chose 70%
#our @NULL_FA = qw(ml.0.50.fa ml.0.55.fa ml.0.60.fa ml.0.65.fa ml.0.70.fa
#            ml.0.75.fa ml.0.80.fa ml.0.85.fa ml.0.90.fa ml.0.95.fa ml.1.00.fa);
#
#our @H_TARG_FA = qw(ml_haeck_others.0.50.fa ml_haeck_others.0.55.fa
#                    ml_haeck_others.0.60.fa ml_haeck_others.0.65.fa
#                    ml_haeck_others.0.70.fa ml_haeck_others.0.75.fa
#                    ml_haeck_others.0.80.fa ml_haeck_others.0.85.fa
#                    ml_haeck_others.0.90.fa ml_haeck_others.0.95.fa
#                    ml_haeck_others.1.00.fa);
#
#our @O_TARG_FA = qw(ml_others.0.50.fa ml_others.0.55.fa ml_others.0.60.fa
#                    ml_others.0.65.fa ml_others.0.70.fa ml_others.0.75.fa
#                    ml_others.0.80.fa ml_others.0.85.fa ml_others.0.90.fa
#                    ml_others.0.95.fa ml_others.1.00.fa);

# comment below 3 lines if you uncomment above
our @NULL_FA = qw(ml.0.70.fa);
our @H_TARG_FA = qw(ml_haeck_others.0.70.fa);
our @O_TARG_FA = qw(ml_others.0.70.fa);

MAIN: {
    print_header();
    for (my $i = 0; $i < @NULL_FA; $i++) {
        my $rh_null = get_ids_from_fasta("$FA_DIR/$NULL_FA[$i]");
        my $rh_h    = get_ids_from_fasta("$FA_DIR/$H_TARG_FA[$i]");
        my $rh_o    = get_ids_from_fasta("$FA_DIR/$O_TARG_FA[$i]");

# (ratio > 1 means higher expression after tentacle development)
# uncomment next 9 lines to print expr ratios 
#my $rh_b = get_data($FILE,$rh_h);
#foreach my $key (sort {$rh_b->{$b} <=> $rh_b->{$a}} keys %{$rh_b}) {
#    print "$key,$rh_b->{$key}\n";
#}
#my $rh_c = get_data($FILE,$rh_o);
#foreach my $key (sort {$rh_b->{$b} <=> $rh_b->{$a}} keys %{$rh_b}) {
#    print "$key,$rh_b->{$key}\n";
#}
#exit;
        my $rh_dat  = get_data($FILE,$rh_null);

        my $h_out = process_data($H_TARG_FA[$i],$rh_null,$rh_h,$rh_dat);
        my $h_rand = 0;
        for (my $i = 0; $i < $REPS; $i++) {
            my $h_ml = mcmc(scalar(keys(%{$rh_h})),$rh_dat);
            $h_rand++ if ($h_ml >= $h_out);
        }
        my $h_pval = $h_rand / $REPS;
        print "$h_pval\n";

        my $o_out = process_data($O_TARG_FA[$i],$rh_null,$rh_o,$rh_dat);
        my $o_rand = 0;
        for (my $i = 0; $i < $REPS; $i++) {
            my $o_ml = mcmc(scalar(keys(%{$rh_o})),$rh_dat);
            $o_rand++ if ($o_ml >= $o_out);
        }
        my $o_pval = $o_rand / $REPS;
        print "$o_pval\n";
    } 
}

sub mcmc {
    my $num = shift;
    my $rh_dat = shift;
    my @rand = shuffle(keys %{$rh_dat});
    my %rand_subset = ();
    for (my $i = 0; $i < $num; $i++) {
        $rand_subset{$rand[$i]}++;
    }
    my ($more_later,$ra_vals) = get_more_later($rh_dat,\%rand_subset);
    return $more_later;
}    

sub print_header {
    print "#file,null_total,null_more_later,targ_total,targ_more_later,pval\n";
}

sub get_more_later {
    my $rh_dat = shift;
    my $rh_targ = shift;
    my $targ = 0;

    foreach my $id (keys %{$rh_targ}) {
        die "unexpected:$id" unless (defined($rh_dat->{$id}));
        $targ++ if ($rh_dat->{$id} > 1);
    }
    return $targ;
}

sub process_data {
    my $file    = shift;
    my $rh_null = shift;
    my $rh_targ = shift;
    my $rh_dat  = shift;

    my @vals = values %{$rh_dat};
    my $more_later = 0;
    foreach my $v (@vals) {
        $more_later++ if ($v > 1);
    }
 
    my $total = scalar(@vals);
    print "$file,$total,$more_later";

    my $targ_total = scalar(keys %{$rh_targ});
    my $targ = get_more_later($rh_dat,$rh_targ);
    print ",$targ_total,$targ,";
    return $targ;
}

sub get_percent {
    my $ra_expr = shift;
    my $percent = 0;


    return $percent;
}

sub get_ids_from_fasta {
    my $file = shift;
    my %ids = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        next unless ($line =~ m/^>(\S+)/);
        $ids{$1}++;
    }
    return \%ids;
}


# input = timecoursefile
# input = null ids  #ids from which to build distribution
# output = %data
# %data = ('ML0001a' => ['avg_time0', 'avg_time1'];
sub get_data {
    my $file = shift;
    my $rh_ids = shift;
    my %data = ();
    open IN, $file or die "cannot open $file:$!";
    my $devnull = <IN>;
    my $i_line = <IN>;
    chomp $i_line;
    my @indices = split /\t/, $i_line;
    shift @indices;

    my $let = 0;
    my $flt = 0;
    for (my $i = 0; $i < @indices; $i++) {
        $let = $i if $indices[$i] == $LAST_EARLY_TIMEPT;
        if ($indices[$i] == $FIRST_LATE_TIMEPT) {
            $flt = $i;
            last;
        }
    }
    
    while (my $line = <IN>) {
        chomp $line;
        my @f = split /\t/, $line;
        my $id = shift @f;
        $id =~ s/'//g;
        next unless ($rh_ids->{$id});
        my $pre = 0;
        my $pre_ct = 0;
        for (my $i = 0; $i < $let; $i++) {
            $pre += $f[$i];
            $pre_ct++;
        }
        my $post = 0;
        my $post_ct = 0;
        for (my $i = $flt; $i < @indices; $i++) {
            $post += $f[$i];
            $post_ct++;
        }
        $data{$id} = ($post / $post_ct) / ($pre / $pre_ct + 1);
    }
    return \%data;
}

