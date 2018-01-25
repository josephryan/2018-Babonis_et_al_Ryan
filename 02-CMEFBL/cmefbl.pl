#!/usr/bin/perl

use strict;
use warnings;
use List::Util 'shuffle';
use Data::Dumper;

our $VERSION = 0.03;  # this version ignores transcripts with zero (does better job at this than v0.02;

our $SRAND = 420;
srand($SRAND);

our $FILE = 'ML_all_datapoints.tab';
our $LAST_EARLY_TIMEPT = 8.5;
our $FIRST_LATE_TIMEPT = 13;
our $REPS = 10000;
our $MINIMUM_EXPR_TO_BE_COUNTED = 1;

our $FA_DIR = '.';
our @NULL_FA = qw(ml.0.50.fa ml.0.55.fa ml.0.60.fa ml.0.65.fa ml.0.70.fa
            ml.0.75.fa ml.0.80.fa ml.0.85.fa ml.0.90.fa ml.0.95.fa ml.1.00.fa);

our @H_TARG_FA = qw(ml_haeck_others.0.50.fa ml_haeck_others.0.55.fa
                    ml_haeck_others.0.60.fa ml_haeck_others.0.65.fa
                    ml_haeck_others.0.70.fa ml_haeck_others.0.75.fa
                    ml_haeck_others.0.80.fa ml_haeck_others.0.85.fa
                    ml_haeck_others.0.90.fa ml_haeck_others.0.95.fa
                    ml_haeck_others.1.00.fa);

our @O_TARG_FA = qw(ml_others.0.50.fa ml_others.0.55.fa ml_others.0.60.fa
                    ml_others.0.65.fa ml_others.0.70.fa ml_others.0.75.fa
                    ml_others.0.80.fa ml_others.0.85.fa ml_others.0.90.fa
                    ml_others.0.95.fa ml_others.1.00.fa);

# to just run .70 uncomment the following 3 lines
@NULL_FA = qw(ml.0.70.fa);
@H_TARG_FA = qw(ml_haeck_others.0.70.fa);
@O_TARG_FA = qw(ml_others.0.70.fa);

MAIN: {
    print_header();
    for (my $i = 0; $i < @NULL_FA; $i++) {
        my $rh_expr = get_data($FILE);  # get all expr ratios
        my $rh_null = get_ids_from_fasta("$FA_DIR/$NULL_FA[$i]",$rh_expr);
        my $rh_h    = get_ids_from_fasta("$FA_DIR/$H_TARG_FA[$i]",$rh_expr);
        my $rh_o    = get_ids_from_fasta("$FA_DIR/$O_TARG_FA[$i]",$rh_expr);
        my $rh_dat  = get_data($FILE,$rh_null); # get only exprratios in our set

        # uncomment the following for a quick look at expr ratios
        # print_expr_ratios($FILE,$rh_h); exit;
        # print_expr_ratios($FILE,$rh_o); exit;

        my $h_out = process_data($H_TARG_FA[$i],$rh_h,$rh_dat);
        my $h_rand = 0;
        for (my $i = 0; $i < $REPS; $i++) {
            my $h_ml = mcmc(scalar(keys(%{$rh_h})),$rh_dat);
            $h_rand++ if ($h_ml >= $h_out);
        }
        my $h_pval = $h_rand / $REPS;
        print "$h_pval\n";

        my $o_out = process_data($O_TARG_FA[$i],$rh_o,$rh_dat);
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
    print "file,null_total,null_more_later,targ_total,targ_more_later,pval\n";
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

# input: $file = FASTA file - this is for reporting only
# input: $rh_targ = TARGET set of expr values (ids same as $file)
# input: $rh_dat  = 
# output = %data
# %data = ('ML0001a' => ['avg_time0', 'avg_time1'];
sub process_data {
    my $file    = shift;
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
    my $rh_dat = shift;
    my %ids = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        next unless ($line =~ m/^>(\S+)/);
        my $id = $1;
        $ids{$id}++ if ($rh_dat->{$id}); # no expr info if not in rh_dat
    }
    return \%ids;
}


# input = timecoursefile
# input = null ids  #ids from which to build distribution
# return value is a ratio of expression:
#     reads after time point / # samples after the time point
#     over reads before time point / # samples before the time point
#     (values > 1 are higher after time point; values < 1 are higher before)
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
        next if expr_less_than_min(\@f,$MINIMUM_EXPR_TO_BE_COUNTED);
        next if ($rh_ids && !$rh_ids->{$id});
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
        next unless ($post);
        $data{$id} = ($post / $post_ct) / ($pre / $pre_ct + 1);
    }
    return \%data;
}

sub expr_less_than_min {
    my $ra_expr = shift;
    my $min = shift;
    my $expr = 0;
    foreach my $val (@{$ra_expr}) {
        $expr += $val;
    }
    return 1 if ($expr < $min);
    return 0;
}

# this routine is just for debugging
sub print_expr_ratios {
    my $file = shift;
    my $rh_d = shift;
    my $rh_b = get_data($FILE,$rh_d);
    foreach my $key (sort {$rh_b->{$b} <=> $rh_b->{$a}} keys %{$rh_b}) {
        print "$key,$rh_b->{$key}\n";
    }
}

