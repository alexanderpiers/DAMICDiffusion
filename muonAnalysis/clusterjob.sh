#!/bin/bash
## Job Name
#SBATCH --job-name=muon
## Allocation Definition
#SBATCH --account=damic
#SBATCH --partition=damic
## Resources
## Nodes
#SBATCH --nodes=1#
###SBATCH --ntasks-per-node=12
#SBATCH --time=2:00:00
#SBATCH --mem=20G
#SBATCH --workdir=/gscratch/home/apiers/DAMICDiffusion/muonAnalysis

## Execute script
module load icc_18-impi_2018
mpirun -np 1 /gscratch/home/apiers/DAMICDiffusion/muonAnalysis/build/muonanalysis
