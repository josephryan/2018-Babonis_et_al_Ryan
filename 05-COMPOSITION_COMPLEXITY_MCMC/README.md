### Generate complexity and composition stats
#### script requires segmaster binary to be in your path
#### script requires ML2.2.aa, which can be downloaded here:
    https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz

#### run on tentacle candidate genes
perl compcompmcmc.pl ../01-LOST/ml_others.0.70.fa ML2.2.aa > compcompmcmc.others.out

#### run on colloblast candidate genes
perl compcompmcmc.pl ../01-LOST/ml_haeck_others.0.70.fa ML2.2.aa > compcompmcmc.haeck_others.out


## Files

#### compcompmcmc.pl
the script

#### compcompmcmc.others.out
output for tentacle candidate genes

#### compcompmcmc.haeck_others.out
output for colloblast candidate genes 


