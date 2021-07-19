### --- Drug gene-set analysis (DGSA): Parkinson's --- ### 

# DGSA is run with MAGMA: https://ctg.cncr.nl/software/magma

## ----- FORMAT INPUT DATA ----- ##

# get HARMONIZED + LIFTED-OVER (build 38) gwas input ready (cd to .../gwas_folder/) (don't use)
awk '{print $1, $3, $4, $8}' PD2019_harmonized.txt > tmp1
sed 's/chr//g' tmp1 > tmp2
sed 's/omosome/CHR/g' tmp2 > tmp3 
sed 's/ /\t/g' tmp3 > pd_magma_input38.txt
rm tmp*

# format .bim file for step 1
awk '{print $1, $12, $5=0, $2, $3, $4}' NALLS2019-PD-ALL.txt > tmpfile
sed '1d' tmpfile > tmpfile2
sed 's/ /\t/g' tmpfile2 > pd_magma_input37.bim
rm tmpfile*

# format p value file for step 2
awk '{print $12, $1, $2, $8, $7 = $9+$10}' NALLS2019-PD-ALL.txt > tmpfile 
sed '1s/bp_pos/BP/g; 1s/p/P/g; 1s/rsid/SNP/g; 1s/0/NOBS/g; 1s/chromosome/CHR/g' tmpfile > tmpfile2
sed 's/ /\t/g' tmpfile2 > pd_pval37.txt
rm tmpfile*

# for drug gene sets for step 3 moa / indication (change file names)
sed 's/"/ /g' moa_uniq_gsets.txt > tmp
sed 's/^ //g' tmp > tmp1          
sed 's/NA, //g' tmp1 > tmp2
sed 's/NA\t//g' tmp2 > tmp3
sed 's/NA//g' tmp3 > tmp4
awk 'NF>=3' tmp4 > moa_sets_magma.sets
rm tmp*    

## ----- RUN ANALYSIS --- ##

# step 1: annotate gwas summary stats
./magma --annotate window=1,0.5 --snp-loc pd_magma_input37.bim --gene-loc NCBI37.3.gene.loc --out pd-step1

# step 2: gene analysis 
./magma --bfile g1000_eur synonyms=0 --pval pd_pval37.txt ncol=NOBS --gene-annot pd-step1.genes.annot \
--out pd-step2

# step 3: drug-gene-set analysis for moa drug gene sets
./magma --gene-results pd-step2.genes.raw --set-annot moa_sets_magma.sets --settings gene-info --out pd-moa-step3a

# step 3: drug-gene-set analysis for indication drug gene sets
./magma --gene-results pd-step2.genes.raw --set-annot ind_sets_magma.sets --settings gene-info --out pd-ind-step3a

# step 3: drug-gene-set analysis for top n gene sets from LINCS
./magma --gene-results pd-step2.genes.raw --set-annot lincs_top_200.sets --settings gene-info --out pd-topn-step3a

# step 4: lincs drugs pathway analysis
./magma --gene-results pd-step2.genes.raw --gene-covar lincs_sets.cov --model direction-covar=positive --out pd-step4a

### --- END --- ###
