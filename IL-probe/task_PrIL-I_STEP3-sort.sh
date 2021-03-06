#!/bin/bash

# Place to add code
export LC_NUMERIC="en_US.UTF-8"
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
pull_distance_start_name=$(echo ${current_pull_distance[0]} | awk '{print $1}')
pull_distance_end_name=$(echo ${current_pull_distance[1]} | awk '{print $1}')
replica_name=$(echo ${current_replica[0]} | awk '{print $1}')
temperature_annealing_name=$(echo ${current_replica[0]} | awk '{print $2}')

numberofionpairs_name=$(echo ${current_numberofionpairs[0]} | awk '{print $1}')
anion_num_name=$(awk "BEGIN {print "$numberofionpairs_name+$probe_charge_name"}" /dev/null)
cation_num_name=$(awk "BEGIN {print "$(($numberofionpairs_name*$((100-$mol_fraction_name))/100))"}" /dev/null)
additive_num_name=$(awk "BEGIN {print "$numberofionpairs_name*$mol_fraction_name/100"}" /dev/null)

echo
echo SYSMEM COMPOSITION:
echo Number of anions $anion_num_name
echo Number of cations $cation_num_name
echo Number of additives $additive_num_name
echo Distance $pull_distance_start_name

# Define path for storing configurations / simulation data
fullpath_pull=$electrode_name/$cation_name$anion_name/$additive_name$mol_fraction/$probe_name/$temperature_name/$version_name/configuration/$pull_distance_start_name/$surfacecharge_name/pull

cd $dir_experiments/$fullpath_pull
pwd

#create gro-s
gmx trjconv -f NVT_pull.trr -o NVT_pull.gro -s NVT_pull.tpr -sep << EOF
0
q
EOF

count=0
while [ $count -lt 21 ]; do
  name=NVT_pull$count

  exp=1.2-$count/20
  dist=$(awk "BEGIN {print $exp}" /dev/null)

  dist=$(printf "%.2f" $dist)

  mkdir -p ../$dist
  mv $name.gro ../$dist/NVT_pull.gro
  rm *#
  count=$((count+1))
done
