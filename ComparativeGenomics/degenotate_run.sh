#!/bin/bash
#SBATCH --job-name=degenot    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=40G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=degenotate_attempt.log   # Standard output and error log
pwd; hostname; date

mamba activate degenotate_2025_env

annotation=gene_annotation.final.gff3
fasta=Ormia_genome.FINAL.fas
vcf=AllBatches.Filtered.Variants.NOWREADY.WITHTACHINAFERA.vcf.gz
outgroup=outgroup.txt
#exclude=Mainland_Depleta.ids
#output_dir=MK_TEST_withTACHINA
exclude=Hawaii_Depleta.id
output_dir=MK_TEST_withTACHINA_MAINLAND

# Hawaii
degenotate.py -a $annotation -g $fasta -v $vcf -u $outgroup -e $exclude -o $output_dir -sfs

# Mainland
#degenotate.py -a $annotation -g $fasta -v $vcf -u $outgroup -e $exclude -o $output_dir -sfs
