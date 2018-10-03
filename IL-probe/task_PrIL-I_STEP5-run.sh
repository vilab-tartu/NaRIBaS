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
distance_name=$(echo ${current_distance[0]} | awk '{print $1}')
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

# Define path for storing configurations / simulation data
fullpath=$electrode_name/$cation_name$anion_name/$additive_name$mol_fraction/$probe_name/$temperature_name/$version_name/1.20/$surfacecharge_name/$distance_name/$replica_name
fullpath_charge=$electrode_name/$cation_name$anion_name/$additive_name$mol_fraction/$probe_name/$temperature_name/$version_name/1.20/$surfacecharge_name/$distance_name/$replica_name/ch
fullpath_temp=$cation_name$anion_name$additive_name$mol_fraction$probe_name$temperature_name$version_name\_$surfacecharge_name

mkdir -p $dir_experiments/$fullpath
mkdir -p $dir_temp/$fullpath_temp

# Calculate the electrode positions
calc_electrode_positions_in_slab

cd $dir_experiments/$fullpath
pwd

electrode_left=$electrode_name'l'
electrode_right=$electrode_name'r'

exp="$surfacecharge_name * 10^-6 / (1.602 * 10^-19 * $electrodeatoms/($xbox_nm*$ybox_nm) * 10^14)"
charge=$(awk "BEGIN {print $exp}" /dev/null)

exp="-1 * $surfacecharge_name * 10^-6 / (1.602 * 10^-19 * $electrodeatoms/($xbox_nm*$ybox_nm) * 10^14)"
minus_charge=$(awk "BEGIN {print $exp}" /dev/null)

echo 'Start converting the topology files ...'
cp $dir_systempreparation/top/$electrode_name.itp electrode_local.itp
sed -i 's/SED_electrode_charge_SED/'$charge'/g' electrode_local.itp
sed -i 's/SED_minus_electrode_charge_SED/'$minus_charge'/g' electrode_local.itp

sed 's+SED_dir_systempreparation_SED+'$dir_systempreparation'+g' $dir_systempreparation/topol_local_PrIL-S.top > topol_local.top
sed -i 's/SED_anion_name_SED/'$anion_name'/g' topol_local.top
sed -i 's/SED_anion_r_name_SED/'$anion_r_name'/g' topol_local.top
sed -i 's/SED_anion_num_SED/'$anion_num_name'/g' topol_local.top
sed -i 's/SED_cation_name_SED/'$cation_name'/g' topol_local.top
sed -i 's/SED_cation_r_name_SED/'$cation_r_name'/g' topol_local.top
sed -i 's/SED_cation_num_SED/'$cation_num_name'/g' topol_local.top
sed -i 's/SED_additive_name_SED/'$additive_name'/g' topol_local.top
sed -i 's/SED_additive_r_name_SED/'$additive_r_name'/g' topol_local.top
sed -i 's/SED_additive_num_SED/'$additive_num_name'/g' topol_local.top
sed -i 's/SED_probe_name_SED/'$probe_name'/g' topol_local.top
sed -i 's/SED_electrode_SED/'$electrode_name'/g' topol_local.top
sed -i 's/SED_electrode_left_SED/'$electrode_left'/g' topol_local.top
sed -i 's/SED_electrode_right_SED/'$electrode_right'/g' topol_local.top
sed -i 's/SED_electrodeatoms_SED/'$electrodeatoms'/g' topol_local.top

sed 's/SED_temperature_name_SED/'$temperature_name'/g' $dir_systempreparation/5_NVT_probe_production_varTemp.mdp > 5_NVT_probe_run.mdp
sed -i 's/SED_temperature_replica_name_SED/'$temperature_annealing_name'/g' 5_NVT_probe_run.mdp

#Move everything to the rundirectory
mv 5_NVT_probe_run.mdp electrode_local.itp topol_local.top $dir_experiments/$fullpath/

####################################

#Prepare the index file
gmx make_ndx -f NVT_pull.gro -o index.ndx << EOF
keep 0
r $electrode_left
name 1 Cathode
r $electrode_right
name 2 Anode
r $cation_r_name
name 3 Cation
r $anion_r_name
name 4 Anion
r $probe_name
name 5 Probe
q
EOF

##Edit the box size to insert the vacuum slab
#exp="$pos_right_electrode_nm+2*$r_wall_nm"
#zbox_nm=$(awk "BEGIN {print $exp}" /dev/null)  #in nm

##Edit the box size to insert the vacuum slab
#cp NVT.gro NVT_anneal.gro
#sed -i '$d' NVT_anneal.gro
#echo $xbox_nm $ybox_nm $zbox_nm >> NVT_anneal.gro

####################################

#grompp and initiate production run
gmx grompp -f 5_NVT_probe_run.mdp -c NVT_anneal.gro -p topol_local.top -n index.ndx -o NVT -maxwarn 1
#gmx mdrun -deffnm NVT
rm *#
