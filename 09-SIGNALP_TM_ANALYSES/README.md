### download ML2.2.aa from NHGRI
wget https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz

### uncompress ML2.2.aa.gz
gzip -d ML2.2.aa.gz

### run signalp on ML2.2.aa (ML2.2.signalp.out is included so don't need to run)
signalp -n ML.signalp.gff -l ML2.2.signalp.err ML2.2.aa > ML2.2.signalp.out

### run tmhmm on ML2.2.aa (ML2.2_TMHMM.out is included so don't need to run)
tmhmm -short ML2.2.aa > ML2.2_TMHMM.out

### grep only those rows that have transmembranes (i.e., PredHel > 0)
grep -v PredHel=0 ML2.2_TMHMM.out > ML2.2_TM.out

### script to count signal-peps and tms, then run Monte Carlo from full ML2.2
perl sig_mc_v2.pl ../02-CMEFBL/fasta/ml_colloblast_candidates.0.70.fa

### same but on tentacle candidates rather than colloblast candidates
perl sig_mc_v2.pl ../02-CMEFBL/fasta/ml_tentacle_candidates.0.70.fa

