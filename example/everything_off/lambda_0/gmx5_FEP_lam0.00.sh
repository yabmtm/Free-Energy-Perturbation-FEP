#!/bin/bash
#PBS -N ic_ev_nm_f_lam0.00.pf
#PBS -o ic_ev_nm_f_lam0.00.out
#PBS -q normal
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS

cd $PBS_O_WORKDIR

module use -a /home/tuf10875/pkg/modulefiles/
module load gromacs/5.1.2 openmpi

mpirun -np 12 mdrun_mpi -s topol2.tpr -cpi state.cpt -c prod2.gro -maxh 11.9

