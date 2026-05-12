#!/bin/bash
#SBATCH --job-name=tadj   # Job name
#SBATCH --ntasks=8                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=30G                     # Job memory request
#SBATCH --partition=short
#SBATCH --output=tajd.log   # Standard output and error log
pwd; hostname; date


ARIZONA=PBStest/populations/Arizona.txt
VCF1=VCF_by_Chrom/FilteredOrmiaAllSitesBcftools.CHR1.norepeats.vcf.gz
VCF2=VCF_by_Chrom/FilteredOrmiaAllSitesBcftools.CHR2.norepeats.vcf.gz
VCF3=VCF_by_Chrom/FilteredOrmiaAllSitesBcftools.CHR3.norepeats.vcf.gz
VCF4=VCF_by_Chrom/FilteredOrmiaAllSitesBcftools.CHR4.norepeats.vcf.gz
VCF5=VCF_by_Chrom/FilteredOrmiaAllSitesBcftools.CHR5.norepeats.vcf.gz

outdir=TajimasDOutput

activate vcftools_env

vcftools --gzvcf $VCF1 --keep $ARIZONA --TajimaD 10000 --out $outdir/Arizona_CHR1
vcftools --gzvcf $VCF2 --keep $ARIZONA --TajimaD 10000 --out $outdir/Arizona_CHR2
vcftools --gzvcf $VCF3 --keep $ARIZONA --TajimaD 10000 --out $outdir/Arizona_CHR3
vcftools --gzvcf $VCF4 --keep $ARIZONA --TajimaD 10000 --out $outdir/Arizona_CHR4
vcftools --gzvcf $VCF5 --keep $ARIZONA --TajimaD 10000 --out $outdir/Arizona_CHR5
