#!/bin/bash

# script for extending fep simulations that are organized into 20 lambda values
# should be run in one of {everything_off, nonmethyl_off, methyl_off}

PREPARE=TRUE
SUBMIT=TRUE

steps=1750000

if [ "$PREPARE" == "TRUE" ]; then

    for i in {0..20}; do # range of lambdas to prepare extended run
	cd lambda_$i
	j=2 # denotes the extension number

	for extend in {2..15}; do # doing up to 15 extensions for 50ns
            if [ -e topol$extend.tpr ];  then # check how many extensions have been done
		j=$((j+1))
            fi
	done

	if (( j == 2 )); then # start first extension
            gmx convert-tpr -s topol.tpr -nsteps $(($steps*$j)) -o topol2.tpr
	fi

	if (( j > 2 )); then # start subsequent extensions
            last=$(($j-1))
            gmx convert-tpr -s topol$last.tpr -nsteps $(($steps*$j)) -o topol$j.tpr
	fi	# edit the submission script to point to the correct tpr file
	sed "s/mpirun.*/mpirun -np 12 mdrun_mpi -s topol$j.tpr -cpi state.cpt -c prod$j.gro/" gmx5_FEP_lam* > gmx5_FEP_extended.sh 
	cd ..
    done
fi

if [ "$SUBMIT" == "TRUE" ]; then # submit jobs with mpi (PBS)
    for i in {0..20}; do
        cd lambda_$i
        qsub gmx5_FEP_extended.sh
        cd ..
    done
fi

# SLIMY JOBS
# for i in {1..20}; do
#     cd lambda_$i
#     sbatch gmx5_FEP_extended.sh
#     cd ..
