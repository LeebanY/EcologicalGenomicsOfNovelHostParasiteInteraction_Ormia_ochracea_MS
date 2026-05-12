#!/bin/bash
#SBATCH --job-name=FUBAR    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=60G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=FUBAR.log   # Standard output and error log
pwd; hostname; date

conda activate HYPHY

for tree in FUBAR_TREES/*;do
nuc=$(basename "$tree" .output.trees | cut -d'_' -f1)
ali=FUBAR_INPUT/$nuc.fa
HYPHYMPI CPU=16 fubar --alignment $ali --tree $tree --output FUBAR_RESULTS/$nuc.fubar.txt
done
