#!/bin/bash
#SBATCH --job-name=compute_trace  # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=20G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=computetrace.log   # Standard output and error log
pwd; hostname; date

conda activate SINGER_env

python ../../SINGER/SINGER/compute_trace.py -prefix  -m 2.8e-9 -start_index 0 -end_index 99 -output_filename TRACE_OUTPUTS
