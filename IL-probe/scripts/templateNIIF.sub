#!/bin/bash
#SBATCH --partition=prod
#SBATCH --job-name=XXXXXX
#SBATCH --time=12:00:00
#SBATCH --ntasks-per-node=24
#SBATCH -o XXXXXX.out
module load openmpi
module load gromacs/5.1.4-single
mpirun gmx_mpi mdrun -deffnm XXXXXX

