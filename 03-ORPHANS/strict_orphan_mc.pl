#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;
use IO::Uncompress::Gunzip qw($GunzipError);
use List::Util;
use Data::Dumper;

our $SEED  = 42;
our $REPS  = 10000;
our $BLAST = 'ML2.2_v_meta_minusML.blastp.gz';

our $FASTA = 'ML2.2.aa'; 
# ML2.2.aa can be downloade from:
#   https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz

#our @INFILES = qw(../01-LOST/ml_haeck_others.0.70.fa ../01-LOST/ml_others.0.70.fa ../02-CMEFBL/MLHO_87.fa ../02-CMEFBL/MLO_120.fa);
#our @OUTFILES = qw( MLHO_165_orphs.txt MLO_189_orphs.txt MLHO_120_orphs.txt MLO_87_orphs.txt );

our @INFILES = qw(../01-LOST/ml_haeck_others.0.70.fa ../01-LOST/ml_others.0.70.fa);
our @OUTFILES = qw( MLHO_165_orphs.txt MLO_189_orphs.txt);

MAIN: {
    srand($SEED);
    foreach (my $i = 0; $i < @INFILES; $i++) {
        my $file = $INFILES[$i];
        my $outfile = $OUTFILES[$i];
        my $ra_seq_ids  = get_seq_ids($file);
        my $ra_all_ids  = get_seq_ids($FASTA);
        my $rh_bl       = get_bl($BLAST);

        my $num_orphs   = get_num_orphs($rh_bl,$ra_seq_ids,scalar(@{$ra_seq_ids}),$outfile);
        my $total = scalar(@{$ra_seq_ids});
        my $div = '#' x 80;
        print "$div\nINFILE: $file\n  out of $total, $num_orphs have no hits to our 11-taxa animal database (E-Vals <= 0.01)\n";

        my $more_orphs = 0;
        for (my $j = 0; $j < $REPS; $j++) {
            my @shuf = List::Util::shuffle(@{$ra_all_ids});
            my $rnum = get_num_orphs($rh_bl,\@shuf,scalar(@{$ra_seq_ids}));
            $more_orphs++ if ($rnum >= $num_orphs);
        }
        my $pval = $more_orphs / $REPS;
        print "  p-value: $pval ($more_orphs / $REPS)\n";
    }
}

sub get_num_orphs {
    my $rh_bl = shift;
    my $ra_s  = shift;
    my $num   = shift;
    my $out   = shift;
    my $hits  = 0;
    open OUT, ">$out" or die "cannot open $out:$!" if ($out);
    for (my $i = 0; $i < $num; $i++) {
        $hits++ unless ($rh_bl->{$ra_s->[$i]});
        next unless ($out);
        print OUT "$ra_s->[$i]\n" unless ($rh_bl->{$ra_s->[$i]});
    }
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
    my $fh = IO::Uncompress::Gunzip->new($file)
            or die "IO::Uncompress::Gunzip of $file failed: $GunzipError\n";
    while (my $line = $fh->getline()) {
        chomp $line;
        my @f = split /\t/, $line;
        $bl{$f[0]} = 1;
    }
    return \%bl;
}
