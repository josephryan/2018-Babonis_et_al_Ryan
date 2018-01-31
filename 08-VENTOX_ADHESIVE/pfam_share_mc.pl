#!/usr/bin/perl

use strict;
use warnings;
use autodie;
use List::Util;
use Data::Dumper;

our $SEED = 420;
our $REPS = 10000;

our $VENTOX_INTERPRO     = 'ventox.fa.tsv';
our $ADHESIVE_INTERPRO   = 'adhesives.fa.tsv';
our $COLLOBLAST_INTERPRO = 'colloblast_candidates.fa.tsv';
our $TENTACLE_INTERPRO   = 'tentacle_candidates.fa.tsv';
our $ML_INTERPRO         = 'ML2.2.aa.tsv';
our $ML_IDS              = 'ML2.2.ids.txt';
our $NUM_COLLOBLAST      = 189; # number of candidates
our $NUM_TENTACLE        = 165; # number of candidates

MAIN: {
    srand($SEED);
    my $rh_ventox_pf    = get_pfams_from_interpro($VENTOX_INTERPRO);
    my $rh_adhesive_pf  = get_pfams_from_interpro($ADHESIVE_INTERPRO);
    my $rh_collo_pf     = get_pfams_from_interpro($COLLOBLAST_INTERPRO);
    my $rh_tent_pf      = get_pfams_from_interpro($TENTACLE_INTERPRO);
    my $rh_ml_pf        = get_ml_pfs($ML_INTERPRO);
    my $ra_ml_ids       = get_ml_ids($ML_IDS); 
    my $coll_ventox_num = pfams_shared($rh_ventox_pf,$rh_collo_pf);
    my $coll_adhes_num  = pfams_shared($rh_adhesive_pf,$rh_collo_pf);
    my $tent_ventox_num = pfams_shared($rh_ventox_pf,$rh_tent_pf);
    my $tent_adhes_num  = pfams_shared($rh_adhesive_pf,$rh_tent_pf);
    my %gte = ();
    for (1..$REPS) {
        my $ra_c_rand_ids = get_rand_ids($NUM_COLLOBLAST,$ra_ml_ids);
        my $rh_c_rand_pf  = get_rand_pf($ra_c_rand_ids,$rh_ml_pf);
        my $ra_t_rand_ids = get_rand_ids($NUM_TENTACLE,$ra_ml_ids);
        my $rh_t_rand_pf  = get_rand_pf($ra_t_rand_ids,$rh_ml_pf);
        my $vc_count     = pfams_shared($rh_ventox_pf,$rh_c_rand_pf);
        my $vt_count     = pfams_shared($rh_ventox_pf,$rh_t_rand_pf);
        my $ac_count     = pfams_shared($rh_adhesive_pf,$rh_c_rand_pf);
        my $at_count     = pfams_shared($rh_adhesive_pf,$rh_t_rand_pf);
        $gte{v_coll}++ if ($vc_count >= $coll_ventox_num);
        $gte{v_tent}++ if ($vt_count >= $tent_ventox_num);
        $gte{a_coll}++ if ($ac_count >= $coll_adhes_num);
        $gte{a_tent}++ if ($at_count >= $tent_adhes_num);
    }
    foreach my $key (sort keys %gte) {
        my $pval = $gte{$key} / $REPS;
        print "$key: p-value = $pval\n";
    }
    print_pfam_ids($rh_ventox_pf,$rh_adhesive_pf,$rh_collo_pf,$rh_tent_pf);
}

sub print_pfam_ids {
    my $rh_ventox = shift;
    my $rh_adhesive = shift;
    my $rh_collo = shift;
    my $rh_tent = shift;

    print "PFAM ids present in Ventox database:\n";
    foreach my $pf (sort keys %{$rh_ventox}) { print "  $pf\n"; }
    print "\nPFAM ids present in adhesive dataset:\n";
    foreach my $pf (sort keys %{$rh_adhesive}) { print "  $pf\n"; }
    print "\nPFAM ids present in colloblast candidates:\n";
    foreach my $pf (sort keys %{$rh_collo}) { print "  $pf\n"; }
    print "\nPFAM ids present in tentacle candidates:\n";
    foreach my $pf (sort keys %{$rh_tent}) { print "  $pf\n"; }
}

sub get_ml_pfs {
    my $file = shift;
    my %dat  = ();
    open(my $fh, "<", $file);
    while (my $line = <$fh>) {
        my @fields = split /\t/, $line;
        next unless ($fields[3] eq 'Pfam');
        $dat{$fields[0]}->{$fields[4]}++;
    }
    return \%dat;
}

sub get_rand_pf {
    my $ra_rand_ids = shift;
    my $rh_ml_pf    = shift;
    my %pf = ();
    foreach my $id (@{$ra_rand_ids}) {
        foreach my $pfid (keys %{$rh_ml_pf->{$id}}) {
            $pf{$pfid}++;
        }
    }
    return \%pf;
}

sub get_rand_ids {
    my $n = shift;
    my $ra_mlids = shift;
    my @rands = ();
    my @shuf = List::Util::shuffle(@{$ra_mlids});
    for (my $i = 0; $i < $n; $i++) {
        push @rands, $shuf[$i];
    }
    return \@rands;
}

sub pfams_shared {
    my $rh_pf      = shift;
    my $rh_test_pf = shift;
    my $shared     = 0;
    foreach my $pf  (keys %{$rh_pf}) {
        $shared++ if ($rh_test_pf->{$pf});
    }
    return $shared;
}

sub get_ml_ids {
    my $file = shift;
    my @ids  = ();
    open(my $fh, "<", $file);
    while (my $line = <$fh>) {
        chomp $line;
        next if ($line =~ m/^\s*$/); # skip blanks
        push @ids, $line; 
    }
    return \@ids;
}

sub get_pfams_from_interpro {
    my $file = shift;
    my %dat  = ();
    open(my $fh, "<", $file);
    while (my $line = <$fh>) {
        my @fields = split /\t/, $line;
        next unless ($fields[3] eq 'Pfam');
        $dat{$fields[4]}++;
    }
    return \%dat;
}
