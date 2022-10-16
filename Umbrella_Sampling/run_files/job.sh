#!/bin/bash
###SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=100:00:00
#SBATCH --job-name=3_2_fib_usample
##SBATCH --nodelist=compute-0-[0-10]
./a.out

