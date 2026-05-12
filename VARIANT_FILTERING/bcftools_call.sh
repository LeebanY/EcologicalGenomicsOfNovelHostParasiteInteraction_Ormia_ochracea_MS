#!/bin/bash
#SBATCH --job-name=variantcalling    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=60G                     # Job memory request
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --output=variantcalling.log   # Standard output and error log
pwd; hostname; date


conda activate bcftools_env

bamlist=bamlist.txt
ref=Ormia_ochracea-GCA_963402855.1-softmasked.headeredited.fa
output=calls/FinalBatchOrmiaSequences.vcf.gz


bcftools mpileup --threads 16 -f $ref -b $bamlist | bcftools call --threads 16 -m -Oz -f GQ -o $output
