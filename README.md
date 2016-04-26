# Free-Energy-Perturbation-FEP
Scripts for running free energy perturbation simulations using Gromacs.

# runme.sh
Allows specification of host/queue and contains commands to prepare and submit simulations

# mbar_analysis.sh
Copies all dhdl.xvg files into the respective mbar_analysis folder for each thermodynamic leg and runs mbar calculations on it from the alchemical-analysis repository. The resulting free energies are organized into a plottable output file with error bars, and the overall energies are also reported.
