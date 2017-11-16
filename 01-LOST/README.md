### ogs.txt.gz
compressed orthogroups file from orthofinder run

### get_absent_from_beroe_ids.pl
Program looks for orthogroups where no Beroe gene is present but genes
from 70% of other species are present. It prints Mnemiopsis leidyi IDs from 
orthogroups to the file present_in_ml_70perc.out. 

The script also identifies orthogroups where no Beroe gene and no Haeckalia
genes are present, but genes from 70% of other species are present. It
prints Mnemiopsis leidyi IDS from these orthogroups to the file present_in_ml_haeck_70perc.out

### ml_others.0.70.fa
These are the ML sequences in FASTA format corresponding with the ids in
the present_in_ml_70perc.out

### ml_haeck_others.0.70.fa
These are the ML sequences in FASTA format corresponding with the ids in 
the present_in_ml_haeck_70perc.out file

