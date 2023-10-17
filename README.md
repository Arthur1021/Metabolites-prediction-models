# Metabolites-prediction-models
 Metabolites prediction models using TwinsUK dataset

## Example code for association test:
Rscript scripts/FUSION.assoc_test.R --sumstats Alzheimer_GWAS_summary.txt --weights pos_file/M45415.pos --weights_dir models/ --ref_ld_chr LD_ref/M45415.ld.1000g. --chr Z --out M45415_on_AD_association.txt