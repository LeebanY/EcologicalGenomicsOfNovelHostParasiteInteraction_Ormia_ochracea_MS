#!/bin/bash
#SBATCH --job-name=orthofinder    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=250G                     # Job memory request
#SBATCH --partition=himem
#SBATCH --output=orthofinder.log   # Standard output and error log
pwd; hostname; date

activate orthofinder

#orthofinder -t 16 -o Orthofinder_Output_NUCLEOTIDES -d -z -f protein_coding_genes_nucleotides

orthofinder -t 16 -o Orthofinder_Output_SEXCHROM_divergence -d -z -f ormia_outgroup_orthogroup
