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
our $H_OUT = 'present_in_ml_haeck_70perc.out';
our $O_OUT = 'present_in_ml_70perc.out';
our $PERC  = 0.70;

MAIN: {
    my $rh_ogs     = get_ogs_w_no_beroe($FILE);
    my $total_orgs = get_total_orgs($rh_ogs);
    my $cutoff     = int($PERC * $total_orgs);
    my $hcutoff     = int($PERC * ($total_orgs + 1));
    open H_OUT, ">$H_OUT" or die "cannot open $H_OUT:$!";
    open O_OUT, ">$O_OUT" or die "cannot open $O_OUT:$!";
    foreach my $og (keys %{$rh_ogs}) {
        next unless ($rh_ogs->{$og}->{'ml'});
        if ($rh_ogs->{$og}->{'haeck'}) {
            next unless (scalar(keys(%{$rh_ogs->{$og}})) > $hcutoff);
            foreach my $ml (@{$rh_ogs->{$og}->{'ml'}}) {
                print H_OUT "$ml\n";
            }
        } else {
            next unless (scalar(keys(%{$rh_ogs->{$og}})) > $cutoff);
            foreach my $ml (@{$rh_ogs->{$og}->{'ml'}}) {
                print O_OUT "$ml\n";
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
