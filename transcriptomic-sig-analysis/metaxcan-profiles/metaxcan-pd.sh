### --- Impute gene expression signature of PD using MetaXcan --- ###

# NOTE: To download MetaXcan, see here: https://github.com/hakyimlab/MetaXcan
# MetaXcan tutorial, see here: https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS

# ----------------------------------------------------------------------------- #
### ----- STEP 1: Run full harmonization for 2019 PD summary statistics ----- ###
# ----------------------------------------------------------------------------- #

python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_parsing.py \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/gwas_folder/NALLS2019-PD-ALL.txt.gz \
-liftover /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/liftover/hg19ToHg38.over.chain.gz \
-snp_reference_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.txt.gz METADATA \
-output_column_map rsid variant_id \
-output_column_map A2 non_effect_allele \
-output_column_map A1 effect_allele \
-output_column_map b effect_size \
-output_column_map p pvalue \
-output_column_map chromosome chromosome \
--chromosome_format \
-output_column_map bp_pos position \
-output_column_map freq frequency \
-output_column_map n_cases n_cases \
-output_column_map n_sample sample_size \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error n_cases sample_size \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz

# ----------------------------------- #
### ------ STEP 2: Imputation ----- ###
# ----------------------------------- #

# chromosome 1
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr1.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 1 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr1_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 2
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr2.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 2 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr2_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 3
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr3.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 3 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr3_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 4
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr4.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 4 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr4_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 5
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr5.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 5 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr5_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 6
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr6.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 6 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr6_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 7
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr7.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 7 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr7_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 8
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr8.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 8 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr8_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 9
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr9.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 9 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr9_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 10
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr10.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 10 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr10_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 11
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr11.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 11 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr11_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 12
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr12.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 12 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr12_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 13
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr13.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 13 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr13_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 14
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr14.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 14 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr14_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 15
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr15.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 15 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr15_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 16
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr16.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 16 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr16_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 17
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr17.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 17 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr17_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 18
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr18.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 18 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr18_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 19
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr19.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 19 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr19_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 20
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr20.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 20 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr20_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 21
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr21.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 21 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr21_sb0_reg0.1_ff0.01_by_region.txt.gz

# chromosome 22
python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation.py \
-by_region_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/eur_ld.bed.gz \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-parquet_genotype /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/chr22.variants.parquet \
-parquet_genotype_metadata /Volumes/nyb_mac_passport/gwas2drug/metaxcan/data/reference_panel_1000G/variant_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome 22 \
-regularization 0.1 \
-frequency_filter 0.01 \
-sub_batches 10 \
-sub_batch 0 \
--standardise_dosages \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation/PD_chr22_sb0_reg0.1_ff0.01_by_region.txt.gz

# ------------------------------------------------- #
### ----- STEP 3: Merge imputed chromosomes ----- ###
# ------------------------------------------------- #

python /Volumes/nyb_mac_passport/gwas2drug/metaxcan/summary-gwas-imputation-master/src/gwas_summary_imputation_postprocess.py \
-gwas_file /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/PD2019_harmonized.txt.gz \
-folder /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/summary_imputation \
-pattern PD_chr* \
-parsimony 7 \
-output /Volumes/nyb_mac_passport/gwas2drug/metaxcan/output/imputed_PD2019.txt.gz

# ----------------------------------------------------------------- #
### ----- STEP 4: Run S-PrediXcan for each brain model (13) ----- ###
# ----------------------------------------------------------------- #

# NOTE: Need to change the output file names to match the GTEx file names, other wise will get error for step 5
# (did this manually first time round)

# # amygdala
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Amygdala.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Amygdala.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-amygdala.csv


# anterior cingulate cortex 
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Anterior_cingulate_cortex_BA24.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Anterior_cingulate_cortex_BA24.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-ant-cingulate-cortex.csv

# caudate basal ganglia
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Caudate_basal_ganglia.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Caudate_basal_ganglia.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-caud-basal-ganglia.csv

# cerebellar hemisphere
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Cerebellar_Hemisphere.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Cerebellar_Hemisphere.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-cere-hemi.csv


# cerebellum
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Cerebellum.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Cerebellum.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-cerebellum.csv

# cortex
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Cortex.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Cortex.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-cortex.csv

# frontal cortex
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Frontal_Cortex_BA9.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Frontal_Cortex_BA9.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-frontal-cortex.csv

# hippocampus
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Hippocampus.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Hippocampus.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-hippocampus.csv

# hypothalamus
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Hypothalamus.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Hypothalamus.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-hypothalamus.csv

# nucleus accumbens basal ganglia
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Nucleus_accumbens_basal_ganglia.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Nucleus_accumbens_basal_ganglia.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-nuc-accumbens-basalg.csv

# putamen basal ganglia
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Putamen_basal_ganglia.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Putamen_basal_ganglia.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-putamen-basalg.csv

# spinal cord cervical
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Spinal_cord_cervical_c-1.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Spinal_cord_cervical_c-1.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-spinal-cord.csv

# substantia nigra
/Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SPrediXcan.py \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--model_db_path /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Substantia_nigra.db \
--covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr/mashr_Brain_Substantia_nigra.txt.gz \
--keep_non_rsid \
--additional_output \
--model_db_snp_key varID \
--throw \
--output_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output/mashr-pd-sub-nigra.csv

# ---------------------------------------- #
### ----- STEP 5: Run S-MulTiXcan ------ ###
# ---------------------------------------- #

python /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/MetaXcan-master/software/SMulTiXcan.py \
--models_folder /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/eqtl/mashr \
--models_name_pattern "mashr_Brain_(.*).db" \
--snp_covariance /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/data/models/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz \
--metaxcan_folder /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/s-predi-output \
--metaxcan_filter "mashr_pd_(.*).csv" \
--metaxcan_file_name_parse_pattern "(.*)_pd_(.*).csv" \
--gwas_file /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/imputed_PD2019.txt.gz \
--snp_column panel_variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--keep_non_rsid \
--model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output /Users/nyb-macbook/genetics-tools/so_methods/metaxcan/output/multixcan-out/PD2019-smulti-output.txt

### --- END --- ###





