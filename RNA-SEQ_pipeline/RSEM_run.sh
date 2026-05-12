#!/bin/bash
#SBATCH --job-name=rsem    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=20G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=RSEM.log   # Standard output and error log
pwd; hostname; date

activate rsemenv

rsemindex=star_output/OrmOch.rsem.index

AS2=star_output/AS2/AS2_test_Aligned.toTranscriptome.out.bam
CS2=star_output/CS2/CS2_test_Aligned.toTranscriptome.out.bam
DC3=star_output/DC3/DC3_test_Aligned.toTranscriptome.out.bam
AS3=star_output/AS3/AS3_test_Aligned.toTranscriptome.out.bam
CS3=star_output/CS3/CS3_test_Aligned.toTranscriptome.out.bam
DC6=star_output/DC6/DC6_test_Aligned.toTranscriptome.out.bam
B2R1=star_output/B2R1/B2R1_test_Aligned.toTranscriptome.out.bam
CS5=star_output/CS5/CS5_test_Aligned.toTranscriptome.out.bam
DC7=star_output/DC7/DC7_test_Aligned.toTranscriptome.out.bam
B2R2=star_output/B2R2/B2R2_test_Aligned.toTranscriptome.out.bam
CS6=star_output/CS6/CS6_test_Aligned.toTranscriptome.out.bam
DC8=star_output/DC8/DC8_test_Aligned.toTranscriptome.out.bam
B2R3=star_output/B2R3/B2R3_test_Aligned.toTranscriptome.out.bam
DC1=star_output/DC1/DC1_test_Aligned.toTranscriptome.out.bam
CS1=star_output/CS1/CS1_test_Aligned.toTranscriptome.out.bam
DC2=star_output/DC2/DC2_test_Aligned.toTranscriptome.out.bam



#rsem-calculate-expression -p 16 --paired-end --alignments \
#--estimate-rspd \
#--append-names \
#--seed 123456 \
#--no-bam-output \
#$AS2 $rsemindex rsem_output/AS2

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--append-names \
--seed 123456 \
--no-bam-output \
$AS3 $rsemindex rsem_output/AS3

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--append-names \
--seed 123456 \
--no-bam-output \
$B2R1 $rsemindex rsem_output/BSR1

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--append-names \
--seed 123456 \
--no-bam-output \
$B2R2 $rsemindex rsem_output/B2R2

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$B2R3 $rsemindex rsem_output/B2R3

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$CS1 $rsemindex rsem_output/CS1

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$CS2 $rsemindex rsem_output/CS2

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$CS3 $rsemindex rsem_output/CS3

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$CS5 $rsemindex rsem_output/CS5

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$CS6 $rsemindex rsem_output/CS6

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$DC1 $rsemindex rsem_output/DC1

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$DC2 $rsemindex rsem_output/DC2

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$DC3 $rsemindex rsem_output/DC3

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$DC6 $rsemindex rsem_output/DC6

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$DC7 $rsemindex rsem_output/DC7

rsem-calculate-expression -p 16 --paired-end --alignments \
--estimate-rspd \
--seed 123456 \
--append-names \
--no-bam-output \
$DC8 $rsemindex rsem_output/DC8
