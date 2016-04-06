#!/bin/bash

# Remember to adjust parameters in setup_lambda_pbsjobs.py and run.mdp

PREPARE=FALSE
SUBMIT=TRUE

# source the gromacs installation directory
source /home/tue91994/Programs/gromacs-5.0.4-installed/bin/GMXRC

if [ "$PREPARE" == "TRUE" ]; then # prepare simulations for 20 lambda values

    # create and populate lambda directories
    python mklambdas.py run.mdp solvated2.top solvent_ions_equilibrated.gro index.ndx *.itp

    # generate simulation binaries for all lambda directories
    for i in {0..20}; do
        echo "Generating binary for lambda_"$i
        cd lambda_$i
        grompp -f grompp.mdp -c solvent_ions_equilibrated.gro -p solvated2.top -n index.ndx -maxwarn 1
	echo SOL | genion -s topol.tpr -n index.ndx -p solvated2.top -o solvent_ions_equilibrated.gro -pname NA -nname CL -neutral
        cd ..
    done

    # setup job files for pbs
    python setup_lambda_pbsjobs.py
fi


if [ "$SUBMIT" == "TRUE" ]; then # submit simulations in accordance with parameters in setup_lambda_pbsjobs.py

    ## submit jobs with mpi (PBS)
    for i in {11..20}; do
        number=$(echo "scale=2; $i/20" | bc | awk '{printf "%.2f\n", $0}')
        cd lambda_$i
        qsub gmx5_FEP_lam$number.sh
        cd ..
    done
fi


# SLIME
# for i in {1..20}; do
#     number=$(echo "scale=2; $i/20" | bc | awk '{printf "%.2f\n", $0}')
#     cd lambda_$i
#     sbatch gmx5_FEP_lam$number.sh
#     cd ..

