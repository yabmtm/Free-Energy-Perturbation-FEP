import os, sys, string
import numpy as np


lambdas = np.arange(0.0, 1.05, 0.050)

# write SLURM job batch script 

for i in range(len(lambdas)):

    fout = open('lambda_%d/gmx5_FEP_lam%2.2f.sh'%(i,lambdas[i]), 'w')

    # write job script header

    fout.write("""#!/bin/bash
#SBATCH -J g_methyl_lam%2.2f.pf                     # Job name
#SBATCH -o g_methyl_lam%2.2f.%%j.out                 # Name of stdout output file (%%j expands to jobId)
#SBATCH -p gpu                           # Queue name
#SBATCH -N 1                             # Total number of nodes requested (16 cores/node)
#SBATCH -n 16                            # Total number of mpi tasks requested
#SBATCH --mail-user=tud16919@temple.edu
#SBATCH --mail-type=end                  # email me when the job finishes
#SBATCH -t 24:00:00                      # Run time (hh:mm:ss) 
#SBATCH -A TG-MCB140270                  # acount name
#SBATCH -e errors.%%j.out                 # error file

"""%(lambdas[i],lambdas[i])  )

    fout.write('module load boost\n')
    fout.write('module load gromacs/5.0.4\n')
    fout.write('module load cuda/6.0\n')
    fout.write('export OMP_NUM_THREADS=16\n\n')
    fout.write('mdrun_gpu -s topol.tpr -maxh 23.5\n')

    fout.close()


