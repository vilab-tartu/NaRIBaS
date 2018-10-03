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
replica_name=$(echo ${current_replica[0]} | awk '{print $1}')
pull_distance_start_name=$(echo ${current_pull_distance[0]} | awk '{print $1}')

numberofionpairs_name=$(echo ${current_numberofionpairs[0]} | awk '{print $1}')

electrode_left=$electrode_name'l'
electrode_right=$electrode_name'r'
anion_num_name=$(awk "BEGIN {print "$numberofionpairs_name+$probe_charge_name"}" /dev/null)
cation_num_name=$(awk "BEGIN {print "$(($numberofionpairs_name*$((100-$mol_fraction_name))/100))"}" /dev/null)
additive_num_name=$(awk "BEGIN {print "$numberofionpairs_name*$mol_fraction_name/100"}" /dev/null)

echo
echo SYSMEM COMPOSITION:
echo Number of anions $anion_num_name
echo Number of cations $cation_num_name

# Define path for storing configurations / simulation data
fullpath_start=$electrode_name/$cation_name$anion_name/$additive_name$mol_fraction/$probe_name/$temperature_name/$version_name/configuration/$pull_distance_start_name/
fullpath_temp=$cation_name$anion_name$additive_name$mol_fraction$probe_name$temperature_name$version_name

mkdir -p $dir_experiments/$fullpath_start
mkdir -p $dir_temp/$fullpath_temp
mkdir -p $dir_systempreparation

# Calculate the electrode positions
calc_electrode_positions_in_slab

# How many slices to take into account? deltaz is the width
#deltaz=0.015 #nm
#calc_slicing_of_slab
# Two properties should now be available: $numberofslices and $currentdens

####################################

# Start with packmol: Add all particle names and calculated values concerning the box size to the packmpl script.

cp -rp sysprep/* $dir_systempreparation

cd $dir_systempreparation
echo 'Start converting the packmol input file packmol.inp ...'
sed 's+SED_dir_systempreparation_SED+'$dir_systempreparation'+g' packmol_PrIL-I.inp > packmol.inp

seeds=( 1983757 3542638 6328548 6679881 7552618 9305733 9361537 39571623 )
seed=${seeds[$replica_name]}
echo $replica_name
sed -i 's/SED_seed_SED/'$seed'/g' packmol.inp

sed -i 's/SED_anion_name_SED/'$anion_name'/g' packmol.inp
sed -i 's/SED_anion_num_SED/'$anion_num_name'/g' packmol.inp
sed -i 's/SED_cation_name_SED/'$cation_name'/g' packmol.inp
sed -i 's/SED_cation_num_SED/'$cation_num_name'/g' packmol.inp
sed -i 's/SED_additive_name_SED/'$additive_name'/g' packmol.inp
sed -i 's/SED_additive_num_SED/'$additive_num_name'/g' packmol.inp
sed -i 's/SED_probe_name_SED/'$probe_name'/g' packmol.inp

sed -i 's/SED_xbox_SED/'$xbox'/g' packmol.inp
sed -i 's/SED_ybox_SED/'$ybox'/g' packmol.inp
sed -i 's/SED_zbox_left_SED/'$z1_ion'/g' packmol.inp
sed -i 's/SED_zbox_right_SED/'$z2_ion'/g' packmol.inp

zboxal=$(awk "BEGIN {print "$z1_ion+5"}" /dev/null)
zboxar=$(awk "BEGIN {print "$z2_ion-5"}" /dev/null)

sed -i 's/SED_zbox_left_ad_SED/'$zboxal'/g' packmol.inp
sed -i 's/SED_zbox_right_ad_SED/'$zboxar'/g' packmol.inp

xbox12=$(awk "BEGIN {print "$xbox*0.5"}" /dev/null)
ybox12=$(awk "BEGIN {print "$ybox*0.5"}" /dev/null)

sed -i 's/SED_xbox12_SED/'$xbox12'/g' packmol.inp
sed -i 's/SED_ybox12_SED/'$ybox12'/g' packmol.inp

zprobepos_r=$(awk "BEGIN {print "$zprobepos+$pos_left_electrode+10*$pull_distance_start_name"}" /dev/null)

sed -i 's/SED_xpos_SED/'$xprobepos'/g' packmol.inp
sed -i 's/SED_ypos_SED/'$yprobepos'/g' packmol.inp
sed -i 's/SED_zpos_SED/'$zprobepos_r'/g' packmol.inp
sed -i 's/SED_xrot_SED/'$xproberot'/g' packmol.inp
sed -i 's/SED_yrot_SED/'$yproberot'/g' packmol.inp
sed -i 's/SED_zrot_SED/'$zproberot'/g' packmol.inp

sed -i 's/SED_electrode_name_left_SED/'$electrode_left'/g' packmol.inp
sed -i 's/SED_electrode_name_right_SED/'$electrode_right'/g' packmol.inp

sed -i 's/SED_pos_left_electrode_SED/'$pos_left_electrode'/g' packmol.inp
sed -i 's/SED_pos_right_electrode_SED/'$pos_right_electrode'/g' packmol.inp

echo 'Run packmol and convert packmol output to gromacs input ...'
$dir_packmol/packmol < packmol.inp
gmx editconf -f packmol.pdb -o packmol.gro
rm packmol.pdb packmol.inp

#read -p "Press enter to continue..."

#+++++++++++++++++++++++

#Edit the box size to insert the vacuum slab
sed -i '$d' packmol.gro
echo $xbox_nm $ybox_nm $zbox_vacuum >> packmol.gro

#Prepare the index file
gmx make_ndx -f packmol.gro -o index.ndx << EOF
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

#Add all necessary .mdp files
sed 's/SED_temperature_name_SED/'$temperature_name'/g' 0_STEEP_probe_varTemp.mdp > 0_STEEP_probe.mdp
sed 's/SED_temperature_name_SED/'$temperature_name'/g' 1_NVT_probe_lowtimestep_varTemp.mdp > 1_NVT_probe_lowtimestep.mdp

echo 'Start converting the topology files ...'
cp $dir_systempreparation/top/$electrode_name.itp electrode_local.itp
sed -i 's/SED_electrode_charge_SED/'0.0000'/g' electrode_local.itp

sed 's+SED_dir_systempreparation_SED+'$dir_systempreparation'+g' $dir_systempreparation/topol_local_PrIL-I.top > topol_local.top
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
#Move everything to the rundirectory
mv 0_STEEP_probe.mdp 1_NVT_probe_lowtimestep.mdp electrode_local.itp topol_local.top packmol.gro index.ndx $dir_experiments/$fullpath_start/

####################################
cd $dir_experiments/$fullpath_start/

#grompp and run the energy minimization
gmx grompp -f 0_STEEP_probe.mdp -c packmol.gro -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
rm mdout.mdp
gmx mdrun -deffnm STEEP

#grompp the first equilibration step
gmx grompp -f 1_NVT_probe_lowtimestep.mdp -c STEEP.gro -p topol_local.top -n index.ndx -o NVT_lowtimestep -maxwarn 1
rm mdout.mdp
rm *#
gmx mdrun -deffnm NVT_lowtimestep
