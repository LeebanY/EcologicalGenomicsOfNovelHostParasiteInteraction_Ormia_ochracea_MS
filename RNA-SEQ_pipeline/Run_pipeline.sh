#!/bin/bash
#SBATCH --job-name=speccontin    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=20G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=specconti.log   # Standard output and error log
pwd; hostname; date

activate snakemake_dros
#REFERENCE should the only one with .fas as ending.

#snakemake -s Kallisto.snakefile --unlock --cluster "sbatch --job-name=snakemake_test --nodes=1 --ntasks=16 --partition=medium --mem=20G" -j 1 --conda-frontend conda
snakemake -s fastp.snakefile --cluster "sbatch --job-name=snakemake_test --nodes=1 --ntasks=16 --partition=medium --mem=20G" -j 1 --conda-frontend conda
#mkdir calls

