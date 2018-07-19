# colloblast and tentacle candidates in single-cell data

### sebe-pedros_etal_mc.pl
Script counts how many candidates (input id file) are in the Sebé-Pedrós et al. data (the sample count = X), then determines the greatest number of genes from these data that occur in a single cluster (this is the test statistic). It then randomly draws X identifiers from the Sebé-Pedrós data and compares the greatest number of genes from this sampled data set occur in a single cluster. It does this 10,000 times and determines a p value by dividing the number of the sampled data that are greater than or equal to the test statistic by 10,000. 

### 165.id
identifiers from tentacle candidates

### 189.id
identifiers from colloblast candidates

### cell_clusters_sebe-pedros.csv
data file generated from supplemental table 4 from Sebe-Pedros et al. 2018  

Sebé-Pedrós, Arnau, et al. "Early metazoan cell type diversity and the evolution of multicellular gene regulation." Nature Ecology & Evolution (2018): 1.  

#### LINK to orig xls file:   
https://static-content.springer.com/esm/art%3A10.1038%2Fs41559-018-0575-6/MediaObjects/41559_2018_575_MOESM7_ESM.xlsx

