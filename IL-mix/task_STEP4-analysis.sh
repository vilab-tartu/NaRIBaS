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
anion1_name=$(echo ${current_anion1[0]} | awk '{print $1}')
anion1_size_factor=$(echo ${current_anion1[0]} | awk '{print $2}')
anion1_r_name=$(echo ${current_anion1[0]} | awk '{print $3}')
anion2_name=$(echo ${current_anion2[0]} | awk '{print $1}')
anion2_size_factor=$(echo ${current_anion2[0]} | awk '{print $2}')
anion2_r_name=$(echo ${current_anion2[0]} | awk '{print $3}')
cation_name=$(echo ${current_cation[0]} | awk '{print $1}')
cation_size_factor=$(echo ${current_cation[0]} | awk '{print $2}')
cation_r_name=$(echo ${current_cation[0]} | awk '{print $3}')
temperature_name=$(echo ${current_temperature[0]} | awk '{print $1}')
version_name=$(echo ${current_version[0]} | awk '{print $1}')
numberofionpairs_name=$(echo ${current_numberofionpairs[0]} | awk '{print $1}')
initial_box_size=$(echo ${current_numberofionpairs[0]} | awk '{print $2}') #in A
molefraction_name=$(echo ${current_numberofionpairs[0]} | awk '{print $3}')

exp="$numberofionpairs_name*(100-$molefraction_name)/100"
anion1_n_name=$(awk "BEGIN {print $exp}" /dev/null)

exp="$numberofionpairs_name-$anion1_n_name"
anion2_n_name=$(awk "BEGIN {print $exp}" /dev/null)

echo
echo SYSMEM COMPOSITION:
echo Number of anions 1 $anion1_n_name
echo Number of anions 2 $anion2_n_name
echo Number of cations $numberofionpairs_name

# Define path for storing configurations / simulation data
fullpath_start=$cation_name$anion1_name/$anion2_name$molefraction_name/$temperature_name/$version_name/
fullpath_temp=$temperature_name$version_name

mkdir -p $dir_experiments/$fullpath_start
mkdir -p $dir_analysis/

cd $dir_experiments/$fullpath_start
pwd

####################################
gmx energy -f NPT.edr -s NPT.tpr -o energy  << EOF
13 21
EOF

gmx msd -f NPT.xtc -s NPT.tpr -n index.ndx  -trestart 10 -b 10000 -o msd_cat << EOF
1
EOF
gmx msd -f NPT.xtc -s NPT.tpr -n index.ndx  -trestart 10 -b 10000 -o msd_an1 << EOF
2
EOF
gmx msd -f NPT.xtc -s NPT.tpr -n index.ndx  -trestart 10 -b 10000 -o msd_an2 << EOF
3
EOF

cp energy.xvg $dir_analysis/energy_$cation_name$anion1_name\_$anion2_name$molefraction_name\_$temperature_name$version_name.xvg

cp msd_cat.xvg $dir_analysis/msd_cat_$molefraction_name.xvg
cp msd_an1.xvg $dir_analysis/msd_an1_$molefraction_name.xvg
cp msd_an2.xvg $dir_analysis/msd_an2_$molefraction_name.xvg
