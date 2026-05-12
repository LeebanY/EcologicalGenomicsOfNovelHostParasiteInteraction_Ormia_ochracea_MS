#!/bin/bash
#SBATCH --job-name=label_trees    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=20G                     # Job memory request
#SBATCH --partition=short
#SBATCH --output=label_trees.log   # Standard output and error log
pwd; hostname; date

conda activate HYPHY

# Input and output directories
INPUT_DIR="KeyPartsToSave/GeneTreesToUse"
OUTPUT_DIR="ALL_TREES_LABELLED"

# Make sure output dir exists
mkdir -p "$OUTPUT_DIR"

# Loop over all Newick files in input dir
for treefile in "$INPUT_DIR"/*.txt; do
    # Get the base filename (no directory)
    base=$(basename "$treefile")
    # Construct output path
    outfile="$OUTPUT_DIR/${base}.txt"
    # Run HyPhy label-tree
    hyphy scripts/hyphy-analyses/LabelTrees/label-tree.bf \
        --tree "$treefile" \
        --regexp 'OrmOch' \
        --label foreground \
        --output "$outfile"

    echo "Processed $treefile -> $outfile"
done
