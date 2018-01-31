### get_tentacle_colloblast_candidates.pl
Script identifies orthogroups that contain no Beroe gene, but genes from 70%
of other species (including Mnemiopsis leidyi). If a group has a Haeckalia gene
the Mnemiopsis leidyi gene ids in the group are printed to ml_tentacle_candidates_ids.txt. Otherwise, they are printed to ml_colloblast_candidates_ids.txt.

### ml_tentacle_candidates_ids.txt
output of get_tentacle_colloblast_candidates.pl. 

### ml_tentacle_candidates.fa
These are the ML tentacle candidate sequences in FASTA format (ids from ml_tentacle_candidates_ids.txt)

### ml_colloblast_candidates_ids.txt
output of get_tentacle_colloblast_candidates.pl.

### ml_colloblast_candidates.fa
These are the ML colloblast candidate sequences in FASTA format (ids from ml_colloblast_candidates_ids.txt)

### ogs.txt.gz
compressed orthogroups file from orthofinder run

