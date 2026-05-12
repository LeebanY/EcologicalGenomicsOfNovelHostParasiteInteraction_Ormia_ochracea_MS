#!/bin/bash
#SBATCH --job-name=pixy   # Job name
#SBATCH --ntasks=8                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=30G                     # Job memory request
#SBATCH --partition=short
#SBATCH --output=pixynucdiv.log   # Standard output and error log
pwd; hostname; date

activate pixy_env

popfile=populations_file.txt

#chr1
pixy --stats pi \
--vcf VCF_by_Chrom/CHR1.varsonly.chrrenamed.vcf.gz \
--populations $popfile \
--window_size 100000 \
--n_cores 8 \
--output_folder pixy_nucleotide_diversity \
--output_prefix pixy_nucdiv_per_pop
