#!/bin/bash

# this script will analyze dhdl.xvg data from a set of 20 lambdas in each of three different thermodynamic cycles
# to be run in parent file that has {everything_off, methyl_off, and non_methyl_off}
# output will be saved to final_results.txt in this directory and may be graphed as a series of overall lambdas vs. overall energy

cp /home/tug27224/research/scripts/delg_to_g.py .
echo "0" > del_g.txt

# analyze first two legs of cycle

for file in non_methyl_off everything_off; do
    cd $file
    python copy_dhdl_files.py
    cd mbar_analysis
    /home/tug27224/research/repos/alchemical-analysis/alchemical_analysis/alchemical_analysis.py .
    awk '{print $19}' results.txt | sed '$ d' >> ../../del_g.txt
    cd ../..
done

# treat methyl_off separately because we have to reverse the lambdas

cd methyl_off
python copy_dhdl_files.py
cd mbar_analysis
python /home/tug27224/research/repos/alchemical-analysis/alchemical_analysis/alchemical_analysis.py . # runs mbar analysis
awk '{print $19}' results.txt | sed '$ d' > ../../del_g_neg.txt
cd ../..

sed -i '/^$/d' del_g.txt del_g_neg.txt # removes blank lines
sed -i 's/^/-/' del_g_neg.txt # switches delG to negative
tac del_g_neg.txt >> del_g.txt # flips lambdas and adds on methyl_off del_Gs

python delg_to_g.py > energies.txt # converts del_G's to energies
paste /home/tug27224/testing/0-60.txt energies.txt > final_results.txt # adds x-values (overall lambda numbers)

# cleanup
rm -f energies.txt del_g.txt delg_to_g.py del_g_neg.txt
