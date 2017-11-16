# CMEFBL (Check Mnemiopsis expression for Beroe losses)

### cmefbl.pl

Perl script generates CSV data which includes the following columns

1. file — input FASTA file which are genes present in N% of ctenophore transcriptomes and present in Mnemiopsis leidyi, but absent in Beroe (ml_haeck_others) or absent in Beroe and Haeckalia (ml_others)

2. null_total — number of genes in "file" (file from column 1)

3. null_more_later - number of genes in "file" with higher expression at 13-20-hour timepoints than at 0-8.5-hour timepoints.

4. targ_total - number of genes in target file with 

5. targ_more_later - number of genes in "file" with higher expression at 13-20-hour timepoints than at 0-8.5-hour timepoints.

6. pval - p-value based on 10,000 monte carlo simulations

### ML_all_datapoints.tab

