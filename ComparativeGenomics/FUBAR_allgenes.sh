#!/bin/bash
#SBATCH --job-name=FUBAR    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=60G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=FUBAR.log   # Standard output and error log
pwd; hostname; date

conda activate HYPHY

for tree in ALL_TREES_LABELLED/*;do
nuc=$(basename "$tree" _tree.txt.txt | cut -d'_' -f1)
ali=ALL_ALIGNMENTS_LABELLED/$nuc.fa
HYPHYMPI CPU=16 fubar --alignment $ali --tree $tree --output FUBAR_ALL_GENES_RESULTS/$nuc.fubar.txt
done
