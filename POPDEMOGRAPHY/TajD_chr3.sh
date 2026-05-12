#!/bin/bash
#SBATCH --job-name=tadj   # Job name
#SBATCH --ntasks=8                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=30G                     # Job memory request
#SBATCH --partition=short
#SBATCH --output=tajd.log   # Standard output and error log
pwd; hostname; date


KAUAI=kauai.txt
OAHU=oahu.txt
HILO=hilo.txt
VCF=VCF_by_Chrom/FilteredOrmiaAllSitesBcftools.CHR3.vcf.gz
outdir=TajimasDOutput

activate vcftools_env

vcftools --gzvcf $VCF --keep $KAUAI --TajimaD 10000 --out $outdir/KAUAI_CHR3

vcftools --gzvcf $VCF --keep $OAHU --TajimaD 10000 --out $outdir/OAHU_CHR3

vcftools --gzvcf $VCF --keep $HILO --TajimaD 10000 --out $outdir/HILO_CHR3
