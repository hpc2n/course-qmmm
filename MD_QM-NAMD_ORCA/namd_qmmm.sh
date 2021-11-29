#!/bin/bash
#SBATCH -A Project_ID
#SBATCH -J qmmm
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -n 14

#Load modules necessary for running NAMD
ml purge  > /dev/null 2>&1 
ml GCC/10.3.0  OpenMPI/4.1.1
ml NAMD/2.14-nompi
ml ORCA/5.0.1

namd2 +p4 QMMM-Min.conf   > output_minimization.dat
#namd2 +p4 QMMM-Ann.con    > output_annealing.dat
#namd2 +p4 QMMM-Equi.conf  > output_equilibration.dat 
#namd2 +p4 QMMM.conf       > output_production.dat 

