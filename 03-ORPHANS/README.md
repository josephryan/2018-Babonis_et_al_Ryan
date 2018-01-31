### How to build the meta_minusML.fa BLAST database and run the BLAST

    ```download 'http://ryanlab.whitney.ufl.edu/downloads/alien_index/meta.fa.gz'

    gzip -d meta.fa.gz 

    head -n 1668882 meta.fa > meta.fa.head.1668882.fa

    tail -n 921539 meta.fa > meta.fa.tail.921539

    cat meta.fa.head.1668882.fa meta.fa.tail.921539 > meta_minusML.fa

    rm meta.fa*

    makeblastdb -dbtype prot -in meta_minusML.fa

    download https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz

    gzip -d ML2.2.aa

    blastp -num_threads 60 -query ML2.2.aa -db meta_minusML.fa -outfmt 6 -max_target_seqs 10 -seg yes -evalue 0.01 -out ML2.2_v_meta_minusML.blastp 2> blastp.err &```

#### strict_orphan_mc.pl 

1. uses the BLASTP above to determine if a gene is an orphan (those w/o hits)

2. determines how many of our set of interest are orphans

3. randomly generates 10,000 sets of Mnemiopsis genes the same size as 
   our set of interest and determines how many of these randoms sets 
   have the same number of orphans. From that generates a P-Value

#### identify orphans in colloblast candidates

    perl strict_orphan_mc.pl ../01-LOST/ml_colloblast_candidates.fa

out of 189, 79 have no hits to our 11-taxa animal database with E-Vals at or below 0.01

    p-value: 0 (0 / 10000)

#### identify orphans in tentacle candidates

    perl strict_orphan_mc.pl ../01-LOST/ml_tentacle_candidates.fa

out of 165, 46 have no hits to our 11-taxa animal database with E-Vals at or below 0.01

    p-value: 0.6431 (6431 / 10000)

#################################################
# there are 2 lines in the strict_orphan_mc.pl script that can be uncommented
# to generate ids of orphans.

