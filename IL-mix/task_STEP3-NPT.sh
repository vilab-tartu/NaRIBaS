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
mkdir -p $dir_temp/

#read -p "Press enter to continue..."

#+++++++++++++++++++++++

#Add all necessary .mdp files
cd $dir_systempreparation

sed 's/SED_temperature_name_SED/'$temperature_name'/g' 3_NPT_varTemp.mdp > 3_NPT.mdp

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
mv 3_NPT.mdp topol_local.top $dir_experiments/$fullpath_start/

####################################
cd $dir_experiments/$fullpath_start/

#gmx grompp the first equilibration step

gmx grompp -f 3_NPT.mdp -c NPT_eq.gro -p topol_local.top -n index.ndx -o NPT -maxwarn 1
rm mdout.mdp
#gmx mdrun -deffnm NPT
#rm *xtc *trr *#
rm *#

runningnumber=$((totalnumberofsetups - 1))
# copy tpr files to run them somewhere as multijobs
cp NPT.tpr $dir_temp/NPT$runningnumber.tpr
if [ $totalnumberofsetups -eq 1 ]; then
rm $dir_temp/copyback.txt
fi
pathtocopyback=$(pwd)
echo "mv NPT$runningnumber.tpr $pathtocopyback/NPT.tpr" >> $dir_temp/copyback.txt
echo "mv NPT$runningnumber.log $pathtocopyback/NPT.log" >> $dir_temp/copyback.txt
echo "mv NPT$runningnumber.edr $pathtocopyback/NPT.edr" >> $dir_temp/copyback.txt
echo "mv NPT$runningnumber.xtc $pathtocopyback/NPT.xtc" >> $dir_temp/copyback.txt
echo "mv NPT$runningnumber.out $pathtocopyback/NPT.out" >> $dir_temp/copyback.txt
echo "mv NPT$runningnumber.gro $pathtocopyback/NPT.gro" >> $dir_temp/copyback.txt
