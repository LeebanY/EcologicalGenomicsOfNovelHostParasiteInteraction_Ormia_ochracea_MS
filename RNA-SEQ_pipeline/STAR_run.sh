#!/bin/bash
#SBATCH --job-name=star    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=20G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=STARtest.log   # Standard output and error log
pwd; hostname; date

activate STARenv

fasta_genome=star_output/OrmiaOchracea.genome.newheads.fasta
gff=star_output/OrmOch.rsem.index.gtf

#STAR --runThreadN 16 \
#--runMode genomeGenerate \
#--genomeDir ormia_star_index \
#--genomeFastaFiles $fasta_genome \
#--genomeSAindexNbases 13 \
#--sjdbGTFfile $gff \
#--sjdbOverhang 149


AS2_1=cleaned/AS2_EKRN240015870-1A_223JWCLT4_L8_1.cleaned.fastq
AS2_2=cleaned/AS2_EKRN240015870-1A_223JWCLT4_L8_2.cleaned.fastq

AS3_1=cleaned/AS3_EKRN240015871-1A_223JWCLT4_L8_1.cleaned.fastq
AS3_2=cleaned/AS3_EKRN240015871-1A_223JWCLT4_L8_2.cleaned.fastq

B2R1_1=cleaned/B2R1_EKRN240015865-1A_223MFNLT4_L2_1.cleaned.fastq
B2R1_2=cleaned/B2R1_EKRN240015865-1A_223MFNLT4_L2_2.cleaned.fastq

B2R2_1=cleaned/B2R2_EKRN240015866-1A_223JWCLT4_L8_1.cleaned.fastq
B2R2_2=cleaned/B2R2_EKRN240015866-1A_223JWCLT4_L8_2.cleaned.fastq

B2R3_1=cleaned/B2R3_EKRN240015867-1A_22JV35LT3_L5_1.cleaned.fastq
B2R3_2=cleaned/B2R3_EKRN240015867-1A_22JV35LT3_L5_2.cleaned.fastq

DC1_1=cleaned/DC1_EKRN240015877-1A_223JWCLT4_L8_1.cleaned.fastq
DC1_2=cleaned/DC1_EKRN240015877-1A_223JWCLT4_L8_2.cleaned.fastq

DC2_1=cleaned/DC2_EKRN240015878-1A_22JMTHLT3_L1_1.cleaned.fastq
DC2_2=cleaned/DC2_EKRN240015878-1A_22JMTHLT3_L1_2.cleaned.fastq

DC3_1=cleaned/DC3_EKRN240015879-1A_223JWCLT4_L8_1.cleaned.fastq
DC3_2=cleaned/DC3_EKRN240015879-1A_223JWCLT4_L8_2.cleaned.fastq

DC6_1=cleaned/DC6_EKRN240015885-1A_222JKJLT4_L4_1.cleaned.fastq
DC6_2=cleaned/DC6_EKRN240015885-1A_222JKJLT4_L4_2.cleaned.fastq
DC6_3=cleaned/DC6_EKRN240015885-1A_223JWCLT4_L8_1.cleaned.fastq
DC6_4=cleaned/DC6_EKRN240015885-1A_223JWCLT4_L8_2.cleaned.fastq

DC7_1=cleaned/DC7_EKRN240015886-1A_222JKJLT4_L4_1.cleaned.fastq
DC7_2=cleaned/DC7_EKRN240015886-1A_222JKJLT4_L4_2.cleaned.fastq
DC7_3=cleaned/DC7_EKRN240015886-1A_223JWCLT4_L8_1.cleaned.fastq
DC7_4=cleaned/DC7_EKRN240015886-1A_223JWCLT4_L8_2.cleaned.fastq

DC8_1=cleaned/DC8_EKRN240015887-1A_223JWCLT4_L8_1.cleaned.fastq
DC8_2=cleaned/DC8_EKRN240015887-1A_223JWCLT4_L8_2.cleaned.fastq

CS1_1=cleaned/CS1_EKRN240015873-1A_223JWCLT4_L8_1.cleaned.fastq
CS1_2=cleaned/CS1_EKRN240015873-1A_223JWCLT4_L8_2.cleaned.fastq

CS2_1=cleaned/CS2_EKRN240015874-1A_223JWCLT4_L8_1.cleaned.fastq
CS2_2=cleaned/CS2_EKRN240015874-1A_223JWCLT4_L8_2.cleaned.fastq

CS3_1=cleaned/CS3_EKRN240015875-1A_223JWCLT4_L8_1.cleaned.fastq
CS3_2=cleaned/CS3_EKRN240015875-1A_223JWCLT4_L8_2.cleaned.fastq

CS5_1=cleaned/CS5_EKRN240015881-1A_223MFNLT4_L2_1.cleaned.fastq
CS5_2=cleaned/CS5_EKRN240015881-1A_223MFNLT4_L2_2.cleaned.fastq

CS6_1=cleaned/CS6_EKRN240015882-1A_222JKJLT4_L4_1.cleaned.fastq
CS6_2=cleaned/CS6_EKRN240015882-1A_222JKJLT4_L4_2.cleaned.fastq
CS6_3=cleaned/CS6_EKRN240015882-1A_223JWCLT4_L8_1.cleaned.fastq
CS6_4=cleaned/CS6_EKRN240015882-1A_223JWCLT4_L8_2.cleaned.fastq
CS6_5=cleaned/CS6_EKRN240015882-1A_22HM5FLT3_L2_1.cleaned.fastq
CS6_6=cleaned/CS6_EKRN240015882-1A_22HM5FLT3_L2_2.cleaned.fastq

#AS2
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $AS2_1 $AS2_2 \
--outFileNamePrefix star_output/AS2/AS2_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#AS3
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $AS3_1 $AS3_2 \
--outFileNamePrefix star_output/AS3/AS3_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#B2R1
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $B2R1_1 $B2R1_2 \
--outFileNamePrefix star_output/B2R1/B2R1_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#B2R2
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $B2R2_1 $B2R2_2 \
--outFileNamePrefix star_output/B2R2/B2R2_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#B2R3
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $B2R3_1 $B2R3_2 \
--outFileNamePrefix star_output/B2R3/B2R3_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#CS1
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $CS1_1 $CS1_2 \
--outFileNamePrefix star_output/CS1/CS1_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#CS2
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $CS2_1 $CS2_2 \
--outFileNamePrefix star_output/CS2/CS2_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#CS3
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $CS3_1 $CS3_2 \
--outFileNamePrefix star_output/CS3/CS3_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#CS5
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $CS5_1 $CS5_2 \
--outFileNamePrefix star_output/CS5/CS5_test_ \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMunmapped Within \
--outSAMattributes Standard

#CS6
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn cleaned/CS6_1.cleaned.fastq cleaned/CS6_2.cleaned.fastq \
--outFileNamePrefix star_output/CS6/CS6_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#DC1
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $DC1_1 $DC1_2 \
--outFileNamePrefix star_output/DC1/DC1_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#DC2
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $DC2_1 $DC2_2 \
--outFileNamePrefix star_output/DC2/DC2_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#DC3
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $DC3_1 $DC3_2 \
--outFileNamePrefix star_output/DC3/DC3_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#DC6
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn cleaned/DC6_1.cleaned.fastq cleaned/DC6_2.cleaned.fastq \
--outFileNamePrefix star_output/DC6/DC6_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#DC7
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn cleaned/DC7_1.cleaned.fastq cleaned/DC7_2.cleaned.fastq \
--outFileNamePrefix star_output/DC7/DC7_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM

#DC8
STAR --genomeDir star_output/ormia_star_index \
--runThreadN 16 \
--readFilesIn $DC8_1 $DC8_2 \
--outFileNamePrefix star_output/DC8/DC8_test_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--quantMode TranscriptomeSAM
