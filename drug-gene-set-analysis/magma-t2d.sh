### ---- Drug gene-set analysis (DGSA): Type 2 diabetes ---- ###

# DGSA is run in MAGMA: https://ctg.cncr.nl/software/magma 

## ----- FORMAT INPUT DATA ----- ##

# format .bim file for step 1
awk '{print $1, $10, $5=0, $2, $3, $4}' T2D_Mahajan2018.txt > tmpfile
sed '1d' tmpfile > tmpfile2
sed 's/ /\t/g' tmpfile2 > t2d_magma_input37.bim
rm tmpfile*

# format p value file for step 2
awk '{print $10, $1, $2, $8, $9}' T2D_Mahajan2018.txt > tmpfile 
sed '1s/pos/BP/g; 1s/p/P/g; 1s/rsid/SNP/g; 1s/n_samPle/NOBS/g; 1s/chromosome/CHR/g' tmpfile > tmpfile2
sed 's/ /\t/g' tmpfile2 > t2d_pval37.txt
rm tmpfile*
# this file will cause errors because the p-values are too small
# need to go in and manually change p-values less than 1.0e-305 to 1.0e-303. Roughly 10 SNPs that need to be changed. 

# for drug gene sets for step 3 
sed 's/"/ /g' entrez_genesets.txt > tmp
sed 's/^ //g' tmp > tmp1           
sed 's/NA, //g' tmp1 > tmp2 
sed 's/, /\t/g' tmp2 > tmp3
awk 'NF>=3' tmp3 > entrez_drugs.sets

# format lincs gene sets with expression data for step 4
sed 's/"//g' lincs_entrez_gsets.txt > lincs_sets.cov

## ----- RUN ANALYSIS ----- ##

# step 1: annotate gwas summary stats
./magma --annotate window=1,0.5 --snp-loc t2d_magma_input37.bim --gene-loc NCBI37.3.gene.loc --out t2d-step1

# step 2: gene analysis 
./magma --bfile g1000_eur synonyms=0 --pval t2d_pval37.txt ncol=NOBS --gene-annot t2d-step1.genes.annot \
--out t2d-step2       # error here if the p-values are too small (i.e., > 1.0e-303) 

# step 3: drug-gene-set analysis 
./magma --gene-results t2d-step2.genes.raw --set-annot entrez_drugs.sets --settings gene-info --out t2d-step3a-targ

# step 3: drug-gene-set analysis for moa drug gene sets
./magma --gene-results t2d-step2.genes.raw --set-annot moa_sets_magma.sets --settings gene-info --out t2d-moa-step3a

# step 3: drug-gene-set analysis for indication drug gene sets
./magma --gene-results t2d-step2.genes.raw --set-annot ind_sets_magma.sets --settings gene-info --out t2d-ind-step3a

# step 4: lincs drugs pathway analysis
./magma --gene-results t2d-step2.genes.raw --gene-covar lincs_sets.cov --model direction-covar=positive --out t2d-step4a-all


### --- END --- ###
