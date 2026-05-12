#!/bin/bash
#SBATCH --job-name=vcftools_filt   # Job name
#SBATCH --ntasks=16                   # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=30G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=vcftools_filt_chrx.log   # Standard output and error log
pwd; hostname; date


#Input files
VCF_CHR1=OrmiaHawaiiAllSitesBcftools.chr1.vcf.gz
VCF_CHR2=OrmiaHawaiiAllSitesBcftools.chr2.vcf.gz
VCF_CHRX=OrmiaHawaiiAllSitesBcftools.CHRX.vcf.gz

#Output files
OUT_CHR1=FilteredOrmiaAllSitesBcftools.CHR1.vcf.gz
OUT_CHR2=FilteredOrmiaAllSitesBcftools.CHR2.vcf.gz
OUT_CHRX=FilteredOrmiaAllSitesBcftools.CHRX.vcf.gz

#Now specify the filters.

MISS=0.8
QUAL=30
MIN_DEPTH=10

activate vcftools_env

bcftools view --threads 16 -M 2 -i 'MIN(INFO/MQ)>30 && MIN(INFO/DP)>10 && F_MISSING<0.1' $VCF_CHRX | sed 's/ENA|OY727110|OY727110.1/chrX/g' | bgzip -c > $OUT_CHRX
