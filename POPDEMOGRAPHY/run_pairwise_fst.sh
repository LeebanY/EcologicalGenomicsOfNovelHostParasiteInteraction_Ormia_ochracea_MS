#!/bin/bash
#SBATCH --job-name=pairwise_fst   # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=30G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=pairwisefst.log   # Standard output and error log
pwd; hostname; date

activate vcftools_env

VCF_chr2=../VCF_by_Chrom/FilteredOrmiaAllSitesBcftools.CHR2.norepeats.vcf.gz
VCF_chr3=../VCF_by_Chrom/FilteredOrmiaAllSitesBcftools.CHR3.norepeats.vcf.gz
VCF_chr4=../VCF_by_Chrom/FilteredOrmiaAllSitesBcftools.CHR4.norepeats.vcf.gz
VCF_chr5=../VCF_by_Chrom/FilteredOrmiaAllSitesBcftools.CHR5.norepeats.vcf.gz
VCF_chrX=../VCF_by_Chrom/FilteredOrmiaAllSitesBcftools.CHRX.norepeats.vcf.gz
popfile=PBS_test.IDs

# calculate pairwise fst takes a population file with two ta-separated columns (samples and pops) and an uncompressed vcf file.

bash calculate_pairwise_fst.v3.sh $popfile $VCF_chr2
mv pairwise_fst_output/*.fst pairwise_fst_output/chr2/

bash calculate_pairwise_fst.v3.sh $popfile $VCF_chr3
mv pairwise_fst_output/*.fst pairwise_fst_output/chr3/

bash calculate_pairwise_fst.v3.sh $popfile $VCF_chr4
mv pairwise_fst_output/*.fst pairwise_fst_output/chr4/

bash calculate_pairwise_fst.v3.sh $popfile $VCF_chr5
mv pairwise_fst_output/*.fst pairwise_fst_output/chr5/

bash calculate_pairwise_fst.v3.sh $popfile $VCF_chrX
mv pairwise_fst_output/*.fst pairwise_fst_output/chrX/
