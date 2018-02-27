#!/bin/bash
#SBATCH -J g_methyl_lam0.65.pf                     # Job name
#SBATCH -o g_methyl_lam0.65.%j.out                 # Name of stdout output file (%j expands to jobId)
#SBATCH -p gpu                           # Queue name
#SBATCH -N 1                             # Total number of nodes requested (16 cores/node)
#SBATCH -n 16                            # Total number of mpi tasks requested
#SBATCH --mail-user=tud16919@temple.edu
#SBATCH --mail-type=end                  # email me when the job finishes
#SBATCH -t 24:00:00                      # Run time (hh:mm:ss) 
#SBATCH -A TG-MCB140270                  # acount name
#SBATCH -e errors.%j.out                 # error file

module load boost
module load gromacs/5.0.4
module load cuda/6.0
export OMP_NUM_THREADS=16

mdrun_gpu -s topol.tpr -maxh 23.5
