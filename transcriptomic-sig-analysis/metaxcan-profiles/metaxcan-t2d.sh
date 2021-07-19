### --- Impute gene expression signature of T2D using MetaXcan --- ###

# NOTE: To download MetaXcan, see here: https://github.com/hakyimlab/MetaXcan
# MetaXcan tutorial, see here: https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS

# ------------------------------------------------------------------------- #
### ----- STEP 1: run full harmonization for T2D summary statistics ----- ###
# ------------------------------------------------------------------------- #

python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_parsing.py \
-gwas_file ~/genetics-tools/metaxcan-for-local/data/gwas_folder/T2D_Mahajan2018.txt.gz \
-liftover ~/genetics-tools/metaxcan-for-local/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map A2 non_effect_allele \
-output_column_map A1 effect_allele \
-output_column_map b effect_size \
-output_column_map chromosome chromosome \
-output_column_map se standard_error \
--chromosome_format \
-output_column_map pos position \
-output_column_map freq frequency \
--insert_value sample_size 898130 --insert_value n_cases 74124 \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size n_cases \
-output ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz

# ---------------------------------------------------------------------------------- #
### ------ STEP 2: Run imputation for each chromosome (should for loop this) ----- ###
# ---------------------------------------------------------------------------------- #

# imputation chromosome 1
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr1.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 1 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr1_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 2
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr2.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 2 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr2_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 3
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr3.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 3 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr3_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 4
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr4.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 4 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr4_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 5
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr5.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 5 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr5_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 6
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr6.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 6 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr6_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 7
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr7.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 7 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr7_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 8
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr8.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 8 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr8_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 9
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr9.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 9 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr9_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 10
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr10.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 10 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr10_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 11
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr11.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 11 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr11_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 12
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr12.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 12 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr12_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 13
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr13.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 13 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr13_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 14
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr14.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 14 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr14_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 15
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr15.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 15 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr15_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 16
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr16.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 16 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr16_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 17
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr17.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 17 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr17_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 18
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr18.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 18 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr18_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 19
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr19.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 19 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr19_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 20
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr20.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 20 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr20_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 21
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr21.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 21 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr21_sb0_reg0.1_ff0.01_by_region.txt.gz

# imputation chromosome 22
python ~/genetics-tools/metaxcan-for-local/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file ~/genetics-tools/metaxcan-for-local/data/eur_ld.bed.gz \
-gwas_file ~/genetics-tools/metaxcan-for-local/output/t2d-maha_harmonized.txt.gz \
-parquet_genotype ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/chr22.variants.parquet \
-parquet_genotype_metadata ~/genetics-tools/metaxcan-for-local/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 22 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output ~/genetics-tools/metaxcan-for-local/output/summary_imputation-t2d/t2d_chr22_sb0_reg0.1_ff0.01_by_region.txt.gz

# ----------------------------------------- #
### ----- STEP 3: Merge chromosomes ----- ###
# ----------------------------------------- #

python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation_postprocess.py \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-folder /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation \
-pattern PD_chr* \
-parsimony 7 \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/imputed_PD2019.txt.gz

# ------------------------------------------------------------------------------ #
### ----- STEP 4: Run S-PrediXcan for every disease relevant tissue (11) ----- ###
# ------------------------------------------------------------------------------ #

# Adipose_Subcutaneous
/Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Adipose_Subcutaneous.db \
--covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Adipose_Subcutaneous.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D/mashr_t2d_Adipose_Subcutaneous.csv


# mashr_Adipose_Visceral_Omentum 
/Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Adipose_Visceral_Omentum.db \
--covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Adipose_Visceral_Omentum.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D/mashr_t2d_Adipose_Visceral_Omentum.csv

# Liver
/Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Liver.db \
--covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Liver.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D/mashr_t2d_Liver.csv

# Muscle_Skeletal
/Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Muscle_Skeletal.db \
--covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Muscle_Skeletal.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D/mashr_t2d_Muscle_Skeletal.csv


# Pancreas 
/Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Pancreas.db \
--covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Pancreas.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D/mashr_t2d_Pancreas.csv

# hypothalamus
/Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Brain_Hypothalamus.db \
--covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Brain_Hypothalamus.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D/mashr_t2d_Brain_Hypothalamus.csv

# Colon_Sigmoid
/Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Colon_Sigmoid.db \
--covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Colon_Sigmoid.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D/mashr_t2d_Colon_Sigmoid.csv

# Colon_Transverse
/Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Colon_Transverse.db \
--covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Colon_Transverse.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D/mashr_t2d_Colon_Transverse.csv

# Kidney_Cortex
/Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Kidney_Cortex.db \
--covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Kidney_Cortex.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D/mashr_t2d_Kidney_Cortex.csv

# Small_Intestine_Terminal_Ileum
/Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Small_Intestine_Terminal_Ileum.db \
--covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Small_Intestine_Terminal_Ileum.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D/mashr_t2d_Small_Intestine_Terminal_Ileum.csv

# Whole_Blood
/Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Whole_Blood.db \
--covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use/mashr_Whole_Blood.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D/mashr_t2d_Whole_Blood.csv

# --------------------------------------- #
### ----- STEP 5: Run S-MulTiXcan ----- ###
# --------------------------------------- #

python /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/MetaXcan-master/software/SMulTiXcan.py \
--models_folder /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/eqtl/mashr/models-to-use \
--models_name_pattern "mashr_(.*).db" \
--snp_covariance /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/data/models/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/s-predi-output/T2D \
--metaxcan_filter "mashr_t2d_(.*).csv" \
--metaxcan_file_name_parse_pattern "(.*)_t2d_(.*).csv" \
--gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/t2d-maha_imputed.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--keep_non_rsid \
--model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output /Volumes/nyb_mac_passport/gwas2drug/metaxcan-for-local/output/T2D2018-smulti-output.txt


### --- END --- ###





