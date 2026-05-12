#!/bin/bash
#SBATCH --job-name=absrel    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=60G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=absrel.log   # Standard output and error log
pwd; hostname; date

conda activate HYPHY

for tree in ABSREL_TREES_LABELLED/*;do
nuc=$(basename "$tree" _tree.txt.txt | cut -d'_' -f1)
ali=ABSREL_PBS_INPUT/$nuc.fa
HYPHYMPI CPU=16 absrel --alignment $ali --tree $tree --output ABSREL_PBS_RESULTS/$nuc.absrel.txt
done
