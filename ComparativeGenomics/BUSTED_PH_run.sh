#!/bin/bash
#SBATCH --job-name=bustedph    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=60G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=bustedph.log   # Standard output and error log
pwd; hostname; date

conda activate HYPHY

for tree in ALL_TREES_LABELLED/*;do
nuc=$(basename "$tree" _tree.txt.txt | cut -d'_' -f1)
ali=ALL_ALIGNMENTS_LABELLED/$nuc.fa
HYPHYMPI CPU=16 scripts/hyphy-analyses/BUSTED-PH/BUSTED-PH.bf --alignment $ali --srv No --tree $tree --branches foreground --output BUSTEDPH_ALL_GENES_RESULTS/$nuc.bustedph.txt
done
