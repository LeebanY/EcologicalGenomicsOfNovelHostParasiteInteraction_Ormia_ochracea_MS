#!/bin/bash
#SBATCH --job-name=clues_array
#SBATCH --array=1-54               # number of lines in positions.txt
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --partition=long
#SBATCH --output=logs/clues_%A_%a.out
#SBATCH --error=logs/clues_%A_%a.err

pwd; hostname; date

# Activate conda environment
conda activate CLUES2_env

# Input file
POSITION_FILE=positions.txt

# Extract the line for this task (tab-delimited: position<TAB>freq)
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $POSITION_FILE)
POS=$(echo $LINE | awk '{print $1}')
FREQ=$(echo $LINE | awk '{print $2}')

# Set chromosome manually (adjust if needed)
CHR=2

# Prefix for output files
OUTPREFIX="AllFlies_Chr${CHR}_POS_${POS}"

echo "Processing position $POS on chromosome $CHR with popFreq=$FREQ"

# Step 1: Extract tree data (SingerToCLUES.py)
python ../CLUES2/SingerToCLUES.py \
--position $POS \
--tree_path Trees/ \
--output $OUTPREFIX

# Step 2: Run inference (inference.py) with allele frequency
python ../CLUES2/inference.py \
--times ${OUTPREFIX}_times.txt \
--popFreq $FREQ \
--N 100000 \
--tCutoff 10000 \
--CI 0.95 \
--out RESULTS_CHR${CHR}_${POS}.outfile_1epochmodel

python ../CLUES2/inference.py \
--times ${OUTPREFIX}_times.txt \
--popFreq $FREQ \
--N 100000 \
--timeBins 1000 \
--tCutoff 10000 \
--CI 0.95 \
--out RESULTS_CHR${CHR}_${POS}.outfile_2epochmodel

python ../CLUES2/inference.py \
--times ${OUTPREFIX}_times.txt \
--popFreq $FREQ \
--N 100000 \
--timeBins 200 1000 \
--tCutoff 10000 \
--CI 0.95 \
--out RESULTS_CHR${CHR}_${POS}.outfile_3epochmodel
