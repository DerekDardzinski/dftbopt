#!/bin/bash -l
#SBATCH -J atomic_dft_test # Job name
#SBATCH -N 4  # Number of nodes
#SBATCH -o stdout # File to which STDOUT will be written %j is the job #
#SBATCH -t 30
#SBATCH -q debug
#SBATCH -A m3578
#SBATCH --constraint=knl
export OMP_NUM_THREADS=1

srun -n 272 -c 4 --cpu_bind=cores python atomic_dft.py > log.out
