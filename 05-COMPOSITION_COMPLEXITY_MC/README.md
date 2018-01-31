### Generate complexity and composition stats
#### script requires segmaster binary to be in your path
#### script requires ML2.2.aa, which can be downloaded here:
    https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz

#### run on tentacle candidate genes
perl compcomp_mc.pl ../01-LOST/ml_tentacle_candidates.fa ML2.2.aa > compcomp_mc.tentacles.out

#### run on colloblast candidate genes
perl compcomp_mc.pl ../01-LOST/ml_colloblast_candidates.fa ML2.2.aa > compcomp_mc.colloblasts.out


## Files

#### compcomp_mc.pl
the script

#### compcomp_mc.tentacles.out
output for tentacle candidate genes

#### compcomp_mc.colloblasts.out
output for colloblast candidate genes 


