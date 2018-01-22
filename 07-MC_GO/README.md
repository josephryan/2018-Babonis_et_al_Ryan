### mcgo.pl
script that runs the GO Monte Carlo analysis
    `./mcgo.pl ML_GO.txt colloblast_candidates.txt MLIDS.txt > mcgo.colloblast_candidates.out`
    `./mcgo.pl ML_GO.txt tentacle_candidates.txt MLIDS.txt > mcgo.tentacle_candidates.out`

### MLIDS.txt
the ids of all the Mnemiopsis leidyi gene models

### ML_GO.txt
GO analysis of all M. leidyi gene models performed in Levin et al. 2016

### colloblast_candidates.txt
list of candidate collablast genes

### mcgo.colloblast_candidates.out
output of `./mcgo.pl ML_GO.txt colloblast_candidates.txt MLIDS.txt`

### tentacle_candidates.txt
list of candidate tentacle genes

### mcgo.tentacle_candidates.out
output of `./mcgo.pl ML_GO.txt tentacle_candidates.txt MLIDS.txt`

## References
Levin, Michal, et al. "The mid-developmental transition and the evolution of animal body plans." Nature 531.7596 (2016): 637-641.
