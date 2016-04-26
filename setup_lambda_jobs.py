import os, sys, string


biglist = [ x / 100. for x in range(0, 105, 5) ]
lambdas = []

for i in biglist:
    lambdas.append(float(format(i, '.2f')))

for i in range(len(lambdas)):

    fout = open('lambda_%d/gmx5_FEP_lam%2.2f.sh'%(i,lambdas[i]), 'w')

######################################

    if sys.argv[1] == "stampede-normal": # Write this batch script if running on stampede-cpu

        fout.write("""#!/bin/bash
#SBATCH -J ic_everything_lam%2.2f.pf                     
#SBATCH -o ic_everything_lam%2.2f.%%j.out                 
#SBATCH -p normal
#SBATCH -N 1                             
#SBATCH -n 16                            
#SBATCH -t 12:00:00                      
#SBATCH -A TG-MCB140270                  
#SBATCH -e errors.%%j.out                 

module load boost cxx11 gromacs

ibrun mdrun_mpi -notunepme -dlb yes -npme -1 -s topol.tpr -o out_prod.gro -maxh 11.9

"""%(lambdas[i],lambdas[i])  )

######################################

    if sys.argv[1] == "cb2rr-normal": # Write this batch script if running on cb2rr-cpu

        fout.write("""#!/bin/bash
#PBS -N ic_ev_nm_f_lam%2.2f.pf
#PBS -o ic_ev_nm_f_lam%2.2f.out
#PBS -q normal
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS

cd $PBS_O_WORKDIR

module use -a /home/tuf10875/pkg/modulefiles/
module load gromacs/5.1.2 openmpi

mpirun -np 12 mdrun_mpi -s topol.tpr -c solvent_ions_afterprod.gro -maxh 11.9

"""%(lambdas[i],lambdas[i])  )

######################################

    if sys.argv[1] == "cb2rr-gpu": # Write this batch script if running on cb2rr-gpu

	fout.write("""#!/bin/bash

#PBS -N ic_ev_nm_f_lam%2.2f.pf
#PBS -o ic_ev_nm_f_lam%2.2f.out
#PBS -q gpu
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS

cd $PBS_O_WORKDIR

module use -a /home/tuf10875/pkg/modulefiles/
module load gromacs/5.1.2 openmpi

OMP_NUM_THREADS=4

mpirun -np 4 mdrun_mpi -ntomp 4 -gpu_id 0123 -s topol.tpr -c solvent_ions_afterprod.gro -maxh 11.9

"""%(lambdas[i],lambdas[i])  )

###########################################

    if sys.argv[1] == "owlsnest-normal": # Write this batch script if running on owlsnest-cpu

	fout.write("""#!/bin/bash

#PBS -N ic_ev_nm_f_lam%2.2f.pf
#PBS -o ic_ev_nm_f_lam%2.2f.out
#PBS -q normal
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00
#PBS

cd $PBS_O_WORKDIR

. /home/tue91994/Programs/gromacs-5.0.4-installed/bin/GMXRC

module load openmpi

mpirun -np 12 mdrun_mpi -s topol.tpr -c solvent_ions_afterprod.gro -maxh 11.9

"""%(lambdas[i],lambdas[i])  )

    fout.close()



