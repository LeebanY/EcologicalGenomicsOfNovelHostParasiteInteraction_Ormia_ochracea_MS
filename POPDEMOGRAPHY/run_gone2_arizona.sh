#!/bin/bash
#SBATCH --job-name=gone_k   # Job name
#SBATCH --ntasks=8                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=40G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=gone_kauai.log   # Standard output and error log
pwd; hostname; date

VCF=../SelectiveSweepHaplotypeScans/Arizona/Arizona.merged.Allchrs.vcf

./GONE2/gone2 -g 2 -r 1.5 $VCF

#./GONE2/gone2 -g 2 -r 1.1 $VCF_1
#./GONE2/gone2 -g 2 -r 1.1 $VCF_2
#./GONE2/gone2 -g 2 -r 1.1 $VCF_3
#./GONE2/gone2 -g 2 -r 1.1 $VCF_4
#./GONE2/gone2 -g 2 -r 1.1 $VCF_5
