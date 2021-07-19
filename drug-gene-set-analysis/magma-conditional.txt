### ----- MAGMA conditional analysis ----- ###

# Run using MAGMA: https://ctg.cncr.nl/software/magma

## ----- FORMAT INPUT ----- ##

# for drug gene sets for step 3 
sed 's/"/ /g' cond_gene_sets.txt > tmp
sed 's/^ //g' tmp > tmp1           
sed 's/NA, //g' tmp1 > tmp2 
sed 's/, /\t/g' tmp2 > tmp3
awk 'NF>=3' tmp3 > conditional.sets
rm tmp*

# for moa gene sets
sed 's/"/ /g' moa_cond_gsets.txt > tmp
sed 's/\tNA\t/ /g' tmp > tmp1
sed 's/NA / /g' tmp1 > moa.sets
rm tmp*

# for moa gene sets
sed 's/"/ /g' ind_cond_gsets.txt > tmp
sed 's/\tNA\t/ /g' tmp > tmp1
sed 's/NA / /g' tmp1 > ind.sets
rm tmp*

# step 3b: drug gene-sets 
./magma --gene-results t2d-step2.genes.raw --set-annot conditional.sets --model analyse=file,drug.signif.txt condition=cond_drug --out t2d.drug.step3b

# step 3b: moa gene-sets 
./magma --gene-results t2d-step2.genes.raw --set-annot moa.sets --model analyse=file,moa.signif.txt condition=cond_drug --out t2d.moa.step3b

# step 3b: indication gene-sets 
./magma --gene-results t2d-step2.genes.raw --set-annot ind.sets --model analyse=file,ind.signif.txt condition=cond_drug --out t2d.ind.step3b

### --- END --- ###
