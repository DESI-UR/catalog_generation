#!/bin/bash -l
#SBATCH -p standard
#SBATCH --mem=10gb
#SBATCH --time=8:00:00
#SBATCH --job-name=five_generation
#SBATCH -o /scratch/tyapici/slurm/logs/mock_generation.log
#SBATCH -e /scratch/tyapici/slurm/errors/mock_generation.err
#SBATCH --mail-user=tyapici@ur.rochester.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1

source activate python3
export LD_LIBRARY_PATH="/home/tyapici/scratch/software/anaconda3/envs/python3/lib"

cd /home/tyapici/scratch/catalog_generation

srun -N 1 python py/generate.py
