#!/usr/bin/perl

use strict;
use warnings;
use IO::Uncompress::Gunzip qw($GunzipError);
use Data::Dumper;

# Author: Joseph Ryan <joseph.ryan@whitney.ufl.edu>

# note: in calculating 70% for genes absent in Beroe but present in at
#       least one of the two Haeckalia species, this script treats the
#       two Haeckalias as a single taxon.

our $VERSION = 0.01;
our $FILE = 'ogs.txt.gz';
our $T_OUT = 'ml_tentacle_candidates_ids.txt';
our $C_OUT = 'ml_colloblast_candidates_ids.txt';
our $PERC  = 0.70;

MAIN: {
    my $rh_ogs     = get_ogs_w_no_beroe($FILE);
    my $total_orgs = get_total_orgs($rh_ogs);
    my $cutoff     = int($PERC * $total_orgs);
    my $hcutoff     = int($PERC * ($total_orgs + 1));
    open T_OUT, ">$T_OUT" or die "cannot open $T_OUT:$!";
    open C_OUT, ">$C_OUT" or die "cannot open $C_OUT:$!";
    foreach my $og (keys %{$rh_ogs}) {
        next unless ($rh_ogs->{$og}->{'ml'});
        if ($rh_ogs->{$og}->{'haeck'}) {
            next unless (scalar(keys(%{$rh_ogs->{$og}})) > $hcutoff);
            foreach my $ml (@{$rh_ogs->{$og}->{'ml'}}) {
                print T_OUT "$ml\n";
            }
        } else {
            next unless (scalar(keys(%{$rh_ogs->{$og}})) > $cutoff);
            foreach my $ml (@{$rh_ogs->{$og}->{'ml'}}) {
                print C_OUT "$ml\n";
            }
        }
    }
}

sub get_total_orgs {
    my $rh_og = shift;
    my %all = ();
    foreach my $og (keys %{$rh_og}) {
        foreach my $org (keys %{$rh_og->{$og}}) {
            next if ($org eq 'ml');
            next if ($org eq 'haeck');
            $all{$org} = 1;
        }
    }
    return scalar(keys(%all));
}

sub get_ogs_w_no_beroe {
    my $file = shift;
    my %data = ();
    my $fh = IO::Uncompress::Gunzip->new($file) 
        or die "IO::Uncompress::Gunzip of $file failed: $GunzipError\n";
    LINE: while (my $line = $fh->getline() ) {
        chomp $line;
        my @fields = split /,/, $line;
        my $og = shift @fields;
        foreach my $f (@fields) {
            next LINE if ($f =~ m/^Bero/ || $f =~ m/^beroe/ || $f =~ m/^Bova/);
            next if ($f =~ m/^Ml/);
            if ($f =~ m/^(haeckelia[^_]+_)/) {
                $data{$og}->{haeck} = 1;
            } elsif ($f =~ m/ML2.2.aa_(.*)/) {
                push @{$data{$og}->{'ml'}}, $1;
            } else {
                $f =~ m/^([^_]+)_/ || die "unexpected: $f";
                $data{$og}->{$1} = 1;
            }
        }
    }
    return \%data;
}
