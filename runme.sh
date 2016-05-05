#!/bin/bash

#######################################
# This script will take input parameters to setup a series of alchemical free-energy based on parameters defined in run.mdp
# and topologies defined in solvated2.top solvent_ions_equilibrated.gro, index.ndx, and any *.itp files.
# The alchemical shift will occur over 20 lambda values.
# Valid $HOST names are owlsnest, cb2rr, stampede.
# Valid $QUEUE names are currently normal and gpu.
# The number of steps is set for a normal queue run for each host, but may be altered by uncommenting it below.
#######################################

# parameters may be set here, or as arguments when running the script.

PREPARE=FALSE
EXTEND=FALSE # this feature currently only works for normal queue
SUBMIT=FALSE
HOST=FALSE
QUEUE=normal # assumes normal queue unless --gpu flag is set.
NSTEPS=FALSE

# parse parameters for input
while test $# -gt 0
do
    case "$1" in
        --prepare) PREPARE=TRUE
            ;;
        --extend) EXTEND=TRUE
            ;;
        --submit) SUBMIT=TRUE
            ;;
	--gpu) QUEUE=gpu
            ;;
        --performance) for i in lambda*; do tail $i/md.log | grep Performance; done
            ;;
        *) echo "Possible inputs are: --prepare, --extend, --submit, --performance,"
                echo 'and --gpu.'; exit
            ;;
    esac
    shift
done

### host discovery
if [[ $(hostname) =~ "stampede" ]]; then
    HOST=stampede
elif [[ $(hostname) =~ "cb2rr" ]]; then
    HOST=cb2rr
else
    HOST=owlsnest
fi

### load gromacs and define steps (host-dependent)
if [ "$HOST" == "cb2rr" ]; then
    module use -a /home/tuf10875/pkg/modulefiles
    module load gromacs/5.1.2
    STEPS=1750000
    if [ "$QUEUE" == "gpu" ]; then
	STEPS=4750000
    fi

elif [ "$HOST" == "owlsnest" ]; then
    . /home/tue91994/Programs/gromacs-5.0.4-installed/bin/GMXRC
    STEPS=1750000
    if [ "$QUEUE" == "gpu" ]; then
        STEPS=4750000
    fi

elif [ "$HOST" == "stampede" ]; then
    module load boost cxx11 gromacs
    STEPS=2000000
    if [ "$QUEUE" == "gpu" ]; then
        STEPS=4750000
    fi
fi


### prepare simulations for 20 lambda values
if [ "$PREPARE" == "TRUE" ]; then

    # set correct number of steps in run.mdp, overriding default is NSTEPS is set above.
    case $NSTEPS in
        ''|*[!0-9]*) echo Using $STEPS steps.;;
        *) STEPS=$NSTEPS; echo Using $STEPS steps. ;;
    esac
    sed -i "/^nsteps/s/[^ ]*/$STEPS/17" run.mdp

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


### extend simulations by $STEPS steps
if [ "$EXTEND" == "TRUE" ]; then
        for i in {0..20}; do # range of lambdas to prepare extended run
        cd lambda_$i
        j=2 # denotes the extension number

        for extend in {2..15}; do # doing up to 15 extensions.
            if [ -e topol$extend.tpr ];  then
                j=$((j+1))
            fi
        done

        if (( j == 2 )); then
            gmx convert-tpr -s topol.tpr -nsteps $(($STEPS*$j)) -o topol2.tpr
        fi

        if (( j > 2 )); then
            last=$(($j-1))
            echo "creating topol$j.tpr"
            gmx convert-tpr -s topol$last.tpr -nsteps $(($STEPS*$j)) -o topol$j.tpr
        fi

	if [ "$HOST" == "owlsnest" -o "$HOST" == "cb2rr" ]; then
            sed -i "s/mpirun.*/mpirun -np 12 mdrun_mpi -s topol$j.tpr -cpi state.cpt -c prod$j.gro -maxh 11.9/" gmx5_FEP_lam*

	elif [ "$HOST" == "stampede" ]; then
            sed -i "s/ibrun.*/ibrun mdrun_mpi -notunepme -dlb yes -npme -1 -s topol$j.tpr -cpi state.cpi -c out_prod.gro -maxh 11.9/" gmx5_FEP_lam*
	fi
	
        cd ..
    done
fi


### submit simulations with the relevant batch scheduler
if [ "$SUBMIT" == "TRUE" ]; then
    for i in {0..20}; do
        number=$(echo "scale=2; $i/20" | bc | awk '{printf "%.2f\n", $0}')
        cd lambda_$i

	if [ "$HOST" == "owlsnest" -o "$HOST" == "cb2rr" ]; then
            qsub gmx5_FEP_lam* > ../jobIDs
	elif [ "$HOST" == "stampede" ]; then
            sbatch gmx5_FEP_lam* > ../jobIDs
	fi
            cd ..
    done

    sed -i "s/\..*//" jobIDs

    if [ "$HOST" == "stampede" ]; then
	cat jobIDs | grep Submitted | sed 's/Submitted batch job //' > jobIDs
    fi
fi
