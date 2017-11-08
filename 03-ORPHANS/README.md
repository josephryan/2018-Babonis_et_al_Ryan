### How to build the meta_minusML.fa BLAST database and run the BLAST
lwp-download 'http://ryanlab.whitney.ufl.edu/downloads/alien_index/meta.fa.gz'

gzip -d meta.fa.gz 

head -n 1668882 meta.fa > meta.fa.head.1668882.fa

tail -n 921539 meta.fa > meta.fa.tail.921539

cat meta.fa.head.1668882.fa meta.fa.tail.921539 > meta_minusML.fa

rm meta.fa*

makeblastdb -dbtype prot -in meta_minusML.fa

lwp-download https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz

gzip -d ML2.2.aa

blastp -num_threads 60 -query ML2.2.aa -db meta_minusML.fa -outfmt 6 -max_target_seqs 10 -seg yes -evalue 0.01 -out ML2.2_v_meta_minusML.blastp 2> blastp.err &

#### strict_orphan_mc.pl 

1. uses the BLASTP above to determine if a gene is an orphan (those w/o hits)

2. determines how many of our set of interest are orphans

3. randomly generates 10,000 sets of Mnemiopsis genes the same size as 
   our set of interest and determines how many of these randoms sets 
   have the same number of orphans. From that generates a P-Value

#### identify orphans in the whole set of 189 absent from Beroe and Haeckalia
perl strict_orphan_mc.pl ../ml_others.0.70.fa
    out of 189, 79 have no hits to our 11-taxa animal database with E-Vals at or below 0.01
    p-value: 0 (0 / 10000)

#### identify orphans in the whole set of 165 absent from Beroe
perl strict_orphan_mc.pl ../ml_haeck_others.0.70.fa
    out of 165, 46 have no hits to our 11-taxa animal database with E-Vals at or below 0.01
    p-value: 0.6431 (6431 / 10000)
#################################################
# this is on those expressed higher after onset of tentacle dev

#### identify orphans in those absent from Beroe and Haeckalia and expressed higher after onset of tentacle development
perl strict_orphan_mc.pl ../11-CMEFBL_RESULTS/MLH_87.fa
    out of 87, 21 have no hits with E-Vals at or below 0.01
    p-value: 0.8599 (8599 / 10000)

#### identify orphans in those absent from Beroe and expressed higher after onset of tentacle development
perl strict_orphan_mc.pl ../11-CMEFBL_RESULTS/MLO_120.fa
    out of 120, 52 have no hits to our 11-taxa animal database with E-Vals at or below 0.01
    p-value: 0.0007 (7 / 10000)

#################################################
# there are 2 lines in the strict_orphan_mc.pl script that can be uncommented
# to generate ids of orphans. Those files have been generated:
189_orphs.txt (110 orphs)
165_orphs.txt  (119 orphs)
120_orphs.txt (68 orphs)
87_orphs.txt (66 orphs)

