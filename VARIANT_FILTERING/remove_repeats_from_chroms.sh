#!/bin/bash
#SBATCH --job-name=repeat_filter   # Job name
#SBATCH --ntasks=16                   # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=30G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=vcftools_remove_reps.log   # Standard output and error log
pwd; hostname; date

# REMOVE REPEATS FROM CHROMS

activate bcftools_env

VCF_IN_1=FilteredOrmiaAllSitesBcftools.CHR1.vcf.gz
VCF_IN_2=FilteredOrmiaAllSitesBcftools.CHR2.vcf.gz
VCF_IN_3=FilteredOrmiaAllSitesBcftools.CHR3.vcf.gz
VCF_IN_4=FilteredOrmiaAllSitesBcftools.CHR4.vcf.gz
VCF_IN_5=FilteredOrmiaAllSitesBcftools.CHR5.vcf.gz
VCF_IN_X=FilteredOrmiaAllSitesBcftools.CHRX.vcf.gz


VCF_OUT_1=FilteredOrmiaAllSitesBcftools.CHR1.norepeats.vcf.gz
VCF_OUT_2=FilteredOrmiaAllSitesBcftools.CHR2.norepeats.vcf.gz
VCF_OUT_3=FilteredOrmiaAllSitesBcftools.CHR3.norepeats.vcf.gz
VCF_OUT_4=FilteredOrmiaAllSitesBcftools.CHR4.norepeats.vcf.gz
VCF_OUT_5=FilteredOrmiaAllSitesBcftools.CHR5.norepeats.vcf.gz
VCF_OUT_X=FilteredOrmiaAllSitesBcftools.CHRX.norepeats.vcf.gz

repeat_bed=../OrmiaRepeatOutput.bed

#vcftools --gzvcf $VCF_IN_1 --exclude-bed $repeat_bed --recode --stdout | bgzip -c > \
#$VCF_OUT_1

#vcftools --gzvcf $VCF_IN_2 --exclude-bed $repeat_bed --recode --stdout | bgzip -c > \
#$VCF_OUT_2

#vcftools --gzvcf $VCF_IN_3 --exclude-bed $repeat_bed --recode --stdout | bgzip -c > \
#$VCF_OUT_3

vcftools --gzvcf $VCF_IN_4 --exclude-bed $repeat_bed --recode --stdout | bgzip -c > \
$VCF_OUT_4

#vcftools --gzvcf $VCF_IN_5 --exclude-bed $repeat_bed --recode --stdout | bgzip -c > \
#$VCF_OUT_5

#vcftools --gzvcf $VCF_IN_X --exclude-bed $repeat_bed --recode --stdout | bgzip -c > \
#$VCF_OUT_X

for all in *.norepeats.vcf.gz;do
tabix -p vcf $all
done
