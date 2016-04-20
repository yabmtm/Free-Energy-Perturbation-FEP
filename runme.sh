#!/bin/bash

#######################################
# This script will take input parameters to setup a series of alchemical free-energy based on parameters defined in run.mdp
# and topologies defined in solvated2.top solvent_ions_equilibrated.gro, index.ndx, and any *.itp files.
# The alchemical shift will occur over 20 lambda values.
# Valid $HOST names are owlsnest, cb2rr, stampede.
# Valid $QUEUE names are currently normal and gpu.
# Remember to adjust nsteps in run.mdp and ion charge in topol_Ions*.itp
#######################################


PREPARE=FALSE
SUBMIT=TRUE
HOST=cb2rr
QUEUE=normal


### source the gromacs installation directory
if [ "$HOST" == "cb2rr" ]; then
    module use -a /home/tuf10875/pkg/modulefiles
    module load gromacs/5.1.2

elif [ "$HOST" == "owlsnest" ]; then
    . /home/tue91994/Programs/gromacs-5.0.4-installed/bin/GMXRC

elif [ "$HOST" == "stampede" ]; then
    module load boost cxx11 gromacs

fi

### prepare simulations for 20 lambda values
if [ "$PREPARE" == "TRUE" ]; then

    # create and populate lambda directories
    python mklambdas.py run.mdp solvated2.top solvent_ions_equilibrated.gro index.ndx *.itp

    # generate simulation binaries for all lambda directories
    for i in {0..20}; do
        echo "Generating binary for lambda_"$i
        cd lambda_$i
        gmx grompp -f grompp.mdp -c solvent_ions_equilibrated.gro -p solvated2.top -n index.ndx -maxwarn 1
	echo SOL | gmx genion -s topol.tpr -n index.ndx -p solvated2.top -o solvent_ions_equilibrated.gro -pname NA -nname CL -neutral
        cd ..
    done

    # setup job files for submission
    python setup_lambda_jobs.py $HOST-$QUEUE

fi

### submit simulations with the relevant batch scheduler
if [ "$SUBMIT" == "TRUE" ]; then

    for i in {0..20}; do
        number=$(echo "scale=2; $i/20" | bc | awk '{printf "%.2f\n", $0}')
        cd lambda_$i

	if [ "$HOST" == "owlsnest" -o "$HOST" == "cb2rr" ]; then
            qsub gmx5_FEP_lam$number.sh
	elif [ "$HOST" == "stampede" ]; then
            sbatch gmx5_FEP_lam$number.sh
	fi
            cd ..
    done
fi
