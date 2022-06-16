#!/bin/bash -l
#SBATCH -J test # Job name
#SBATCH -N 1  # Number of nodes
#SBATCH -o stdout # File to which STDOUT will be written %j is the job #
#SBATCH -t 3
#SBATCH -q debug
#SBATCH -A m3578
#SBATCH --constraint=knl
export OMP_NUM_THREADS=1

srun -n 272 -c 1 --cpu_bind=cores python ./test.py > log2.out
