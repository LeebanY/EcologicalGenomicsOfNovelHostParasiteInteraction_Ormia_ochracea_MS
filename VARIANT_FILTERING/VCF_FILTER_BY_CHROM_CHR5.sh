#!/bin/bash
#SBATCH --job-name=vcftools_filt   # Job name
#SBATCH --ntasks=16                   # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=30G                     # Job memory request
#SBATCH --partition=short
#SBATCH --output=vcftools_filt_chr4.log   # Standard output and error log
pwd; hostname; date


#Input files
VCF_CHR1=OrmiaHawaiiAllSitesBcftools.chr1.vcf.gz
VCF_CHR2=OrmiaHawaiiAllSitesBcftools.chr2.vcf.gz
VCF_CHRX=OrmiaHawaiiAllSitesBcftools.CHRX.vcf.gz
VCF_CHR4=OrmiaHawaiiAllSitesBcftools.CHR4.vcf.gz
VCF_CHR5=OrmiaHawaiiAllSitesBcftools.CHR5.vcf.gz

#Output files
OUT_CHR1=FilteredOrmiaAllSitesBcftools.CHR1.vcf.gz
OUT_CHR2=FilteredOrmiaAllSitesBcftools.CHR2.vcf.gz
OUT_CHRX=FilteredOrmiaAllSitesBcftools.CHRX.vcf.gz
OUT_CHR4=FilteredOrmiaAllSitesBcftools.CHR4.vcf.gz
OUT_CHR5=FilteredOrmiaAllSitesBcftools.CHR5.vcf.gz

#Now specify the filters.

MISS=0.8
QUAL=30
MIN_DEPTH=10

activate vcftools_env

bcftools view --threads 16 -M 2 -i 'MIN(INFO/MQ)>30 && MIN(INFO/DP)>10' $VCF_CHR5 | sed 's/ENA|OY727109|OY727109\.1/chr5/g' | bgzip -c > $OUT_CHR5

