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
fullpath_temp=$temperature_name$molefraction_name$version_name

mkdir -p $dir_experiments/$fullpath_start
mkdir -p $dir_temp/$fullpath_temp

####################################

# Start with packmol: Add all particle names and calculated values concerning the box size to the packmpl script.

cd $dir_systempreparation
echo 'Start converting the packmol input file packmol.inp ...'
sed 's+SED_dir_systempreparation_SED+'$dir_systempreparation'+g' packmol_IL-bulk.inp > packmol.inp

seeds=( 6328547 1983757 3542638 6679881 7552618 9305733 9361537 39571623 )
seed=${seeds[0]}
sed -i 's/SED_seed_SED/'$seed'/g' packmol.inp

sed -i 's/SED_anion1_name_SED/'$anion1_name'/g' packmol.inp
sed -i 's/SED_anion1_num_SED/'$anion1_n_name'/g' packmol.inp
sed -i 's/SED_anion2_name_SED/'$anion2_name'/g' packmol.inp
sed -i 's/SED_anion2_num_SED/'$anion2_n_name'/g' packmol.inp
sed -i 's/SED_cation_name_SED/'$cation_name'/g' packmol.inp
sed -i 's/SED_cation_num_SED/'$numberofionpairs_name'/g' packmol.inp

exp="$initial_box_size*$cation_size_factor*($anion1_size_factor*(100-$molefraction_name)+$anion2_size_factor*$molefraction_name)/100"
box=$(awk "BEGIN {print $exp}" /dev/null)
sed -i 's/SED_box_SED/'$box'/g' packmol.inp

echo 'Run packmol and convert packmol output to gromacs input ...'
$dir_packmol/packmol < packmol.inp

#+++++++++++++++++++++++

#Add all necessary .mdp files
sed 's/SED_temperature_name_SED/'$temperature_name'/g' 0_STEEP_varTemp.mdp > 0_STEEP.mdp
sed 's/SED_temperature_name_SED/'$temperature_name'/g' 1_NVT_lowtimestep_varTemp.mdp > 1_NVT_lowtimestep.mdp

echo 'Start converting the topology files ...'

sed 's+SED_dir_systempreparation_SED+'$dir_systempreparation'+g' $dir_systempreparation/topol_local_IL-bulk.top > topol_local.top
sed -i 's/SED_anion1_name_SED/'$anion1_name'/g' topol_local.top
sed -i 's/SED_anion1_r_name_SED/'$anion1_r_name'/g' topol_local.top
sed -i 's/SED_anion1_num_SED/'$anion1_n_name'/g' topol_local.top
sed -i 's/SED_anion2_name_SED/'$anion2_name'/g' topol_local.top
sed -i 's/SED_anion2_r_name_SED/'$anion2_r_name'/g' topol_local.top
sed -i 's/SED_anion2_num_SED/'$anion2_n_name'/g' topol_local.top
sed -i 's/SED_cation_name_SED/'$cation_name'/g' topol_local.top
sed -i 's/SED_cation_r_name_SED/'$cation_r_name'/g' topol_local.top
sed -i 's/SED_cation_num_SED/'$numberofionpairs_name'/g' topol_local.top
#Move everything to the rundirectory
mv 0_STEEP.mdp 1_NVT_lowtimestep.mdp topol_local.top packmol.pdb packmol.inp index.ndx $dir_experiments/$fullpath_start/

####################################
cd $dir_experiments/$fullpath_start/

gmx editconf -f packmol.pdb -o packmol.gro
rm packmol.pdb

#Prepare the index file
gmx make_ndx -f packmol.gro -o index.ndx << EOF
keep 0
r $cation_r_name
name 1 Cation
r $anion1_r_name
name 2 Anion1
r $anion2_r_name
name 3 Anion2
q
EOF

#grompp and run the energy minimization
gmx grompp -f 0_STEEP.mdp -c packmol.gro -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
rm mdout.mdp
gmx mdrun -deffnm STEEP
mv STEEP.gro packmol.gro

gmx grompp -f 0_STEEP.mdp -c packmol.gro -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
rm mdout.mdp
gmx mdrun -deffnm STEEP
mv STEEP.gro packmol.gro

gmx grompp -f 0_STEEP.mdp -c packmol.gro -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
rm mdout.mdp
gmx mdrun -deffnm STEEP
mv STEEP.gro packmol.gro

gmx grompp -f 0_STEEP.mdp -c packmol.gro -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
rm mdout.mdp
gmx mdrun -deffnm STEEP
mv STEEP.gro packmol.gro

gmx grompp -f 0_STEEP.mdp -c packmol.gro -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
rm mdout.mdp
gmx mdrun -deffnm STEEP

#grompp the first equilibration step
gmx grompp -f 1_NVT_lowtimestep.mdp -c STEEP.gro -p topol_local.top -n index.ndx -o NVT_lowtimestep -maxwarn 1
gmx mdrun -deffnm NVT_lowtimestep
rm *xtc *trr *cpt mdout.mdp *#
<<comment-multijob
#copy tpr files to run them somewhere as multijobs
runningnumber=$((totalnumberofsetups - 1))
cp NVT_lowtimestep.tpr $dir_temp/$fullpath_temp/NVT_lowtimestep$runningnumber.tpr
if [ $totalnumberofsetups -eq 1 ]; then
rm $dir_temp/$fullpath_temp/copyback.txt
fi
pathtocopyback=$(pwd)
echo "mv NVT_lowtimestep$runningnumber.tpr $pathtocopyback/NVT_lowtimestep.tpr" >> $dir_temp/$fullpath_temp/copyback.txt
comment-multijob
