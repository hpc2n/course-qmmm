#!/bin/bash
#SBATCH -A Project_ID
#Asking for 10 min.
#SBATCH -t 00:50:00
#Number of nodes
#SBATCH -N 1
#Ask for 28 processes
#SBATCH -n 28
#SBATCH --output=job_o.out
#SBATCH --error=job_o.err
#SBATCH --reservation=*FIXME*

ml purge  > /dev/null 2>&1 
ml GCC/9.3.0  OpenMPI/4.0.3 
ml NAMD/2.14-mpi

mpirun -np 28 namd2 4ake_eq.conf > logfile.txt
