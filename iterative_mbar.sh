#!/bin/bash

for i in {0..15}; do
    echo $i >> count
    mkdir results_$i
    python ~/research/repos/alchemical-analysis/alchemical_analysis/alchemical_analysis.py\
	-s $(($i*1000)) -o results_$i
    awk '/TOTAL/{print $17}' results_$i/results.txt >> energies.txt
    awk '/TOTAL/{print $19}' results_$i/results.txt >> uncertainties.txt
    echo "Finished processing last $((20-$i)) ns."
done
paste energies.txt uncertainties.txt > energy_errors
paste count energy_errors > post_eq_energies.txt
rm -f energies.txt uncertainties.txt count energy_errors
