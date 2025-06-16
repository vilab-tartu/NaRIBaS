#!/bin/bash

# Place to add code
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mdanalysis

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
zbox_nm=$(echo ${current_electrode[0]} | awk '{print $4}') #in nm
r_wall_nm=$(echo ${current_electrode[0]} | awk '{print $5}') #in nm
electrodeatoms=$(echo ${current_electrode[0]} | awk '{print $6}')

ion_name=$(echo ${current_ion[0]} | awk '{print $1}')
r_ion_nm=$(echo ${current_ion[0]} | awk '{print $2}') #in nm
ion_r_name=$(echo ${current_ion[0]} | awk '{print $5}')
ion_charge_name=$(echo ${current_ion[0]} | awk '{print $4}')

#numberofions_name=$(echo ${current_numberofions[0]} | awk '{print $1}') 
temperature_name=$(echo ${current_temperature[0]} | awk '{print $1}')
version_name=$(echo ${current_version[0]} | awk '{print $1}')
replica_name=$(echo ${current_replica[0]} | awk '{print $1}')

####################################
# Calculations
xbox=$(awk "BEGIN {print "$xbox_nm*10.0"}" /dev/null)  #in A
ybox=$(awk "BEGIN {print "$ybox_nm*10.0"}" /dev/null)  #in A
A_box=$(awk "BEGIN {print "$xbox*$ybox"}" /dev/null)	#box area in sqA
A_ion=$(awk "BEGIN {print "4*$r_ion_nm*$r_ion_nm*100"}" /dev/null)	#ion area in sqA
maxnumberofions_name=$(awk "BEGIN {print int("$A_box/$A_ion")}" /dev/null)

####################################

# Define path for storing simulation data
fullpath_start=$electrode_name/$ion_name/$temperature_name/$version_name/$replica_name/
echo $dir_experiments
cp $dir_systempreparation/scripts/script.py $dir_experiments/$fullpath_start/
#cp $dir_systempreparation/scripts/script_theta.py $dir_experiments/$fullpath_start/
cp $dir_systempreparation/scripts/script_potential.py $dir_experiments/$fullpath_start/
cp $dir_systempreparation/scripts/script_distance.py $dir_experiments/$fullpath_start/

####################################
cd $dir_experiments/$fullpath_start/
pwd
#python script.py
echo $maxnumberofions_name
#python script_theta.py 
python script_potential.py
python script_distance.py "$ion_r_name" "$ion_charge_name" "$A_box"
