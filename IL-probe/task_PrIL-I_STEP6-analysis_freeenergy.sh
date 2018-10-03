#!/bin/bash

# Place to add code

# Some basic functions to access the data

#echo ${concentration[*]} # All elements of the array
#echo ${concentration[9]} # Element number 10 (counter starts at 0)
#numberofitems=${#concentration[*]} # calculate the number of elements in array
#echo $numberofitems
#echo ${concentration[1]} | awk '{print $1}' # Access the data that is stored within an element array
#echo '-----'
#print_current_setup

# Transfer list entries to bash variables
electrode_name=$(echo ${current_electrode[0]} | awk '{print $1}')
xbox_nm=$(echo ${current_electrode[0]} | awk '{print $2}') #in nm
ybox_nm=$(echo ${current_electrode[0]} | awk '{print $3}') #in nm
r_wall_nm=$(echo ${current_electrode[0]} | awk '{print $4}') #in nm
electrodeatoms=$(echo ${current_electrode[0]} | awk '{print $5}')

anion_name=$(echo ${current_anion[0]} | awk '{print $1}')
r_anion_nm=$(echo ${current_anion[0]} | awk '{print $2}') #in nm
anion_r_name=$(echo ${current_anion[0]} | awk '{print $3}')
cation_name=$(echo ${current_cation[0]} | awk '{print $1}')
r_cation_nm=$(echo ${current_cation[0]} | awk '{print $2}') #in nm
cation_r_name=$(echo ${current_cation[0]} | awk '{print $3}')
additive_name=$(echo ${current_additive[0]} | awk '{print $1}')
r_additive_nm=$(echo ${current_additive[0]} | awk '{print $2}') #in nm
additive_r_name=$(echo ${current_additive[0]} | awk '{print $3}')
probe_name=$(echo ${current_probe[0]} | awk '{print $1}')
probe_charge_name=$(echo ${current_probe[0]} | awk '{print $2}')
xprobepos=$(echo ${current_orientation[0]} | awk '{print $1}') #in A
yprobepos=$(echo ${current_orientation[0]} | awk '{print $2}') #in A
zprobepos=$(echo ${current_orientation[0]} | awk '{print $3}') #in A
xproberot=$(echo ${current_orientation[0]} | awk '{print $4}')
yproberot=$(echo ${current_orientation[0]} | awk '{print $5}')
zproberot=$(echo ${current_orientation[0]} | awk '{print $6}')

mol_fraction_name=$(echo ${current_mol_fraction[0]} | awk '{print $1}')
temperature_name=$(echo ${current_temperature[0]} | awk '{print $1}')
version_name=$(echo ${current_version[0]} | awk '{print $1}')
surfacecharge_name=$(echo ${current_surfacecharge[0]} | awk '{print $1}')

echo

# Define path for storing configurations / simulation data
fullpath=$electrode_name/$cation_name$anion_name/$additive_name$mol_fraction/$probe_name/$temperature_name/$version_name/1.20/$surfacecharge_name/

# Copy scripts
cp $dir_systempreparation/scripts/average-xvg.py $dir_experiments/$fullpath
cp $dir_systempreparation/scripts/pmf_plot.par $dir_experiments/$fullpath
cp $dir_systempreparation/scripts/pmf-s-fullstd-sc.py $dir_experiments/$fullpath
cp $dir_systempreparation/scripts/pmf-s-sc.py $dir_experiments/$fullpath
cp -r $dir_systempreparation/scripts/Mara $dir_experiments/$fullpath


cd $dir_experiments/$fullpath
pwd

python pmf-s-sc.py
python pmf-s-fullstd-sc.py
xmgrace -p pmf_plot.par -settype xydy output-pmf-sc.dat -settype xydy output-pmf-fullstd-sc.dat -world 0 -20 2 100 -saveall $cation_name$anion_name$additive_name$mol_fraction$probe_name$temperature_name$version_name\_$surfacecharge_name.dat -hardcopy
cp $cation_name$anion_name$additive_name$mol_fraction$probe_name$temperature_name$version_name\_$surfacecharge_name.dat ~/Desktop/
for a in *.png; do convert -trim "$a" "$a"; done
