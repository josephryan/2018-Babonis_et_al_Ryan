# CMEFBL (Check Mnemiopsis expression for Beroe losses)

### cmefbl.pl

to run:  perl cmefbl.pl

output:
file,null_total,null_more_later,targ_total,targ_more_later,pval
ml_tentacle_candidates.0.70.fa,12646,4691,155,87,0
ml_colloblast_candidates.0.70.fa,12646,4691,182,120,0

Perl script generates CSV data which includes the following columns

1. file — FASTA file with genes present in 70% of ctenophore transcriptomes and present in Mnemiopsis leidyi, but absent in Beroe (./fasta/ml_tentacle_candidates.0.70.fa) or absent in Beroe and Haeckalia (./fasta/ml_colloblast_candidates.0.70.fa)

2. null_total — number of genes in all Mnemiopsis leidyi gene models (ML2.2) with total expression in time course > 1

3. null_more_later - number of genes from null_total  with higher expression at 13-20-hour timepoints than at 0-8.5-hour timepoints.

4. targ_total - number of genes in target file (from first column)

5. targ_more_later - number of genes in target file (from first column) with higher expression at 13-20-hour timepoints than at 0-8.5-hour timepoints.

6. pval - p-value based on 10,000 monte carlo simulations

### ML_all_datapoints.tab
time course expression data for Mnemiopsis leidyi gene models (ML2.2)

