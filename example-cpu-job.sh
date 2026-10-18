#!/bin/bash -l
#SBATCH --job-name=examplejob
#SBATCH --output=examplejob.o%j      # Name of stdout output file
#SBATCH --error=examplejob.e%j       # Name of stderr error file
#SBATCH --partition=small
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks=1               # Total number of mpi tasks
#SBATCH --mem=0                      # Allocate all the memory on each node
#SBATCH --time=00:10:00              # Run time (d-hh:mm:ss)
#SBATCH --account=<your_project>   # Project for billing

# Launch MPI code 
srun build/part1_ns --N 64 --Re 10 --T 10.0 --eps 1e-6 --dt_factor 0.4
