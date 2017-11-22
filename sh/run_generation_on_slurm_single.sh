#!/bin/bash -l
#SBATCH -p standard
#SBATCH --time=2:00:00
#SBATCH --job-name=five_generation
#SBATCH -o /scratch/tyapici/slurm/logs/mock_generation.log
#SBATCH -e /scratch/tyapici/slurm/errors/mock_generation.err
#SBATCH --nodes=1

source activate python3

cd /home/tyapici/scratch/gcatgen

srun -N 1 python generate.py

