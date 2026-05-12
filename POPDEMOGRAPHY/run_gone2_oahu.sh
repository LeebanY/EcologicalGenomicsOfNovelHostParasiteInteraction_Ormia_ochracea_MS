#!/bin/bash
#SBATCH --job-name=gone_o   # Job name
#SBATCH --ntasks=8                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=40G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=gone_oahu.log   # Standard output and error log
pwd; hostname; date

VCF_1=../SelectiveSweepHaplotypeScans/phased_vcfs_hawaiian_samples/Oahu/Autosomal_Variants.filtered.norepeats.pruned.phased.chr1.Hawaii.vcf
VCF_2=../SelectiveSweepHaplotypeScans/phased_vcfs_hawaiian_samples/Oahu/Autosomal_Variants.filtered.norepeats.pruned.phased.chr2.Hawaii.vcf
VCF_3=../SelectiveSweepHaplotypeScans/phased_vcfs_hawaiian_samples/Oahu/Autosomal_Variants.filtered.norepeats.pruned.phased.chr3.Hawaii.vcf
VCF_4=../SelectiveSweepHaplotypeScans/phased_vcfs_hawaiian_samples/Oahu/Autosomal_Variants.filtered.norepeats.pruned.phased.chr4.Hawaii.vcf
VCF_5=../SelectiveSweepHaplotypeScans/phased_vcfs_hawaiian_samples/Oahu/Autosomal_Variants.filtered.norepeats.pruned.phased.chr5.Hawaii.vcf
VCF=../SelectiveSweepHaplotypeScans/phased_vcfs_hawaiian_samples/Oahu/merge.vcf

./GONE2/gone2 -g 2 -r 1.5 $VCF

#./GONE2/gone2 -g 2 -r 1.1 $VCF_1
#./GONE2/gone2 -g 2 -r 1.1 $VCF_2
#./GONE2/gone2 -g 2 -r 1.1 $VCF_3
#./GONE2/gone2 -g 2 -r 1.1 $VCF_4
#./GONE2/gone2 -g 2 -r 1.1 $VCF_5
