#!/bin/bash
#SBATCH --job-name=matchpbs    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=20G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=matchpbs.log   # Standard output and error log
pwd; hostname; date

# Your input files
idfile="../Ormia_annotation_data/PBS.transcripts.id"
indir="PAL2NAL"
outfile="PBS_gene_files.txt"

# Clear the output file first
> "$outfile"

# Loop through each identifier line
while IFS= read -r id; do
    # Skip empty lines
    [[ -z "$id" ]] && continue
    grep -rlx -- "$id" "$indir" >> "$outfile"
done < "$idfile"

# Remove duplicate filenames
sort -u "$outfile" -o "$outfile"
