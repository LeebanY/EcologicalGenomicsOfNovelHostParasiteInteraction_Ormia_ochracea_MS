#!/bin/bash
#SBATCH --job-name=gadma   # Job name
#SBATCH --ntasks=3                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=250G                     # Job memory request
#SBATCH --partition=himem
#SBATCH --output=gadma.log   # Standard output and error log
pwd; hostname; date

activate gadma_env2

gadma -p params.file -o gadma_PRUNED-SYNONYMOUS_OrmiaInference_momentsdrawing
