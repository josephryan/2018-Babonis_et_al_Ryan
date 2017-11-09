# 18S tree

### run trees

#### run iqtree
iqtree-omp -nt 15 -m GTR+G4 -s ../00-DATA/cteno_names_w_sp.fa -pre ct_iq

#### run raxml with 10 random starting trees
raxmlHPC-PTHREADS-SSE3 -T 15 -p 1234 -# 10 -m GTRGAMMA --no-bfgs -s ../00-DATA/cteno_names_w_sp.fa -n ct_mp

#### run raxml with 10 parsimony starting trees
raxmlHPC-PTHREADS-SSE3 -T 15 -d -p 1234 -# 10 -m GTRGAMMA --no-bfgs -s ../00-DATA/cteno_names_w_sp.fa -n ct_rt

### evaluate trees

#### use raxml to generate a likelihood score of iqtree
raxmlHPC -f e -m GTRGAMMA --no-bfgs -t ct_iq.treefile -s ../00-DATA/cteno_names_w_sp.fa -n iq_eval

#### use raxml to generate a likelihood score of raxml tree (w/parsimony starts)
raxmlHPC -f e -m GTRGAMMA --no-bfgs -t RAxML_bestTree.ct_mp -s ../00-DATA/cteno_names_w_sp.fa -n mp_eval

#### use raxml to generate a likelihood score of raxml tree (w/random starts)
raxmlHPC -f e -m GTRGAMMA --no-bfgs -t RAxML_bestTree.ct_rt -s ../00-DATA/cteno_names_w_sp.fa -n rt_eval

#### iqtree produced best tree
cat RAxML_log.iq_eval  # 0.950715 -4633.053553
cat RAxML_log.mp_eval  # 0.964952 -4633.298470
cat RAxML_log.rt_eval  # 1.041750 -4633.298502

### bootstrap best tree

#### generate 100 bootstraps
raxmlHPC-PTHREADS-SSE3 -T 10 -p 1234 -m GTRGAMMA -b 12345 -# 100 -s ../00-DATA/cteno_names_w_sp.fa -n bs

#### apply bootstraps to iqtree
raxmlHPC -m GTRCAT -p 12345 -f b -z /bwdata1/jfryan/07-BEROE/07-18S/03-BS/RAxML_bootstrap.bs -n applybs -t ct_iq.treefile


#### clone AfterPhylo for clade collapsing
git clone https://github.com/qiyunzhu/AfterPhylo

#### collapse nodes with less than 50% bs support
perl AfterPhylo/AfterPhylo.pl -format=newick -collapse=50 ../03-BS/RAxML_bipartitions.applybs



