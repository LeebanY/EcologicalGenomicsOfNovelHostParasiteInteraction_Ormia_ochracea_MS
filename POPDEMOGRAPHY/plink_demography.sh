#!/bin/bash
#SBATCH --job-name=plink    # Job name
#SBATCH --ntasks=4                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=20G                     # Job memory request
#SBATCH --partition=short
#SBATCH --output=plink.log   # Standard output and error log
pwd; hostname; date

#VCF=~/scratch/OrmiaPopGenAnalysis/VCF_by_Chrom/Autosomal_Variants.filtered.norepeats.pruned.vcf.gz
#VCF=../VCF_by_Chrom/NoDepletaVCFforPca/Autosomal_Variants.filtered.norepeats.pruned.nodepleta.vcf.gz
VCF=../VCF_by_Chrom/Autosomal_Variants.filtered.norepeats.pruned.variants.Hawaiian.vcf.gz
mainland=../Tree/MainlandAutosomal_Variantsonly.filtered.norepeats.pruned.vcf.gz

activate hapflk_env

#plink2 --vcf $VCF --double-id --allow-extra-chr --bp-space 100 --set-missing-var-ids @:# \
#--make-bed --pca --out HawaiiOrmia_pruned_input

plink2 --vcf $mainland --double-id --allow-extra-chr --bp-space 100 --set-missing-var-ids @:# \
--make-bed --pca --out Mainland_Ormia_pruned_input
