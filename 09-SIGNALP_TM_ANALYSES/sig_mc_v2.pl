#!/usr/bin/perl

use strict;
use warnings; 
use JFR::Fasta;
use List::Util 'shuffle';
use Data::Dumper;

our $SIGNALP_OUT = './ML2.2.signalp.out';
our $TMDATA = './ML2.2_TM.out';
our $FASTA = 'ML2.2.aa'; # download from https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz
our $REPS  = 10000;
our $SEED  = 42;

MAIN: {
    srand($SEED);
    my $fa = $ARGV[0] or die "usage: $0 FASTA\n";
    my $ra_seq_ids  = get_seq_ids($fa);
    my $ra_all_ids  = get_seq_ids($FASTA);
    my $rh_sig      = get_signalp_data($SIGNALP_OUT);
    get_tm_data($rh_sig,$TMDATA);
    my $num_sig     = get_num($rh_sig,'sig',$ra_seq_ids);
    my $num_tm      = get_num($rh_sig,'tm' ,$ra_seq_ids);
    my $num_both    = get_num_overlap($rh_sig,$ra_seq_ids);
    my $all_sig     = get_num($rh_sig,'sig',$ra_all_ids);
    my $all_tm      = get_num($rh_sig,'tm' ,$ra_all_ids);
    my $all_both    = get_num_overlap($rh_sig,$ra_all_ids);
print "\$num_sig = $num_sig | \$num_tm = $num_tm | \$num_both = $num_both\n";
print "\$all_sig = $all_sig | \$all_tm = $all_tm | \$all_both = $all_both\n";
    my $sig_hit = 0;
    my $tm_hit = 0;
    for (1..$REPS) {
        my @shuf = shuffle(@{$ra_all_ids});
        my $num_sig_mc = 0;
        my $num_tm_mc  = 0;
        for (my $i = 0; $i < scalar(@{$ra_seq_ids}); $i++) { 
            $num_sig_mc++ if ($rh_sig->{$ra_all_ids->[$i]}->{'sig'} eq 'Y');
            $num_tm_mc++  if ($rh_sig->{$ra_all_ids->[$i]}->{'tm'});
        } 
        $sig_hit++ if ($num_sig_mc >= $num_sig);
        $tm_hit++ if ($num_tm_mc   >= $num_tm);
    }
    my $sig_pval = $sig_hit / $REPS;
    my $tm_pval  = $tm_hit / $REPS;
    print "sig p-val = $sig_pval ($sig_hit / $REPS)\n";
    print "tm p-val  = $tm_pval  ($tm_hit / $REPS)\n";
}

sub get_tm_data {
    my $rh_sigdat = shift;
    my $tmfile      = shift;

    open IN, $tmfile or die "cannot open $tmfile:$!";
    while (my $line = <IN>) {
        chomp $line;
        my @f = split /\s+/, $line;
        $rh_sigdat->{$f[0]}->{'tm'} = 'Y';
    }
}

sub get_num_overlap {
    my $rh_sig = shift;
    my $ra_ids = shift;
    my $num    = 0;
    foreach my $id (@{$ra_ids}) {
        $num++ if ($rh_sig->{$id}->{'sig'} eq 'Y' && $rh_sig->{$id}->{'tm'});
    }
    return $num;
}
sub get_num {
    my $rh_sig = shift;
    my $feat   = shift;
    my $ra_ids = shift;
    my $num    = 0;
    foreach my $id (@{$ra_ids}) {
        $num++ if ($rh_sig->{$id}->{$feat} && $rh_sig->{$id}->{$feat} eq 'Y');
    }
    return $num;
}

sub get_signalp_data {
    my $sigfile = shift;
    my %sigdat  = ();
    open IN, $sigfile or die "cannot open $sigfile:$!";
    my $devnull_header = <IN>;
    my $devnull_header2 = <IN>;
    while (my $line = <IN>) {
        chomp $line;
        my @f = split /\s+/, $line;
        $sigdat{$f[0]}->{'sig'} = $f[9];
    }
    return \%sigdat;
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

