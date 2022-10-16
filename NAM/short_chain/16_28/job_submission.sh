#! /bin/bash -l
#SBATCH -J NAM_sc16
#SBATCH -o out.NAM_PL
#SBATCH -n 96
#SBATCH -t 200:00:00

./NAM_parallel.py
