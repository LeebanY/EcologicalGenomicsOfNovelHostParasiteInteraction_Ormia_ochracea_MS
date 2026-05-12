#!/bin/bash
#SBATCH --job-name=sweed   # Job name
#SBATCH --ntasks=1                   # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=20G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=sweed_sim.log   # Standard output and error log
pwd; hostname; date

SimulationVCF=TRUE_BOTTLENECK_ORMIA.vcf

./../../../SweeD/SweeD_v3.2.1_Linux/SweeD -name SimulationBottleneck_Ormia -input $SimulationVCF -grid 10000 -strictPolymorphic -folded
