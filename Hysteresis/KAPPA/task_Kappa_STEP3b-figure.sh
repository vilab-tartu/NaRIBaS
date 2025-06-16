#!/bin/bash

# Place to add code
set -euo pipefail

# Transfer list entries to bash variables
electrode_name=$(echo ${current_electrode[0]} | awk '{print $1}')

ion_name=$(echo ${current_ion[0]} | awk '{print $1}')

#numberofions_name=$(echo ${current_numberofions[0]} | awk '{print $1}') 
temperature_name=$(echo ${current_temperature[0]} | awk '{print $1}')
version_name=$(echo ${current_version[0]} | awk '{print $1}')


####################################
mkdir -p $dir_analysis/$ion_name

# Define path for storing simulation data
analysispath=$electrode_name/$ion_name/$temperature_name/$version_name/
echo $dir_experiments
cp $dir_systempreparation/scripts/plot_theta_sigma_potential.py $dir_experiments/$analysispath/
cp $dir_systempreparation/scripts/distance_extract.py $dir_experiments/$analysispath/
cp $dir_systempreparation/scripts/plot_distance.py $dir_experiments/$analysispath/


####################################
cd $dir_experiments/$analysispath/
pwd

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mdanalysis

python plot_theta_sigma_potential.py 
python distance_extract.py
python plot_distance.py

cp theta_sigma_plot.png $dir_analysis/$ion_name
cp sigma_potential_plot.png $dir_analysis/$ion_name
cp distance_graph.png $dir_analysis/$ion_name
cp uM.txt $dir_analysis/$ion_name


