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
anion_name=$(echo ${current_anion[0]} | awk '{print $1}')
cation_name=$(echo ${current_cation[0]} | awk '{print $1}')
functional_name=$(echo ${current_functional_opt[0]} | awk '{print $1}' | sed 's,\/,\\/,g')
basisset_name=$(echo ${current_basisset[0]} | awk '{print $1}' | sed 's,\/,\\/,g')
version_name=$(echo ${current_version[0]} | awk '{print $1}')

#These variables are now accessed as in the following example:
echo SYSMEM PARAMETERS
echo Anion $anion_name
echo Cation $cation_name
echo Functional $functional_name
echo Basis-set $basisset_name

# Define path for storing computational data
fullpath=$(echo $cation_name$anion_name/opt/ | sed 's,\\\/,,g')
mkdir -p $dir_experiments/$fullpath

# Make new input file by substituting in variables to the template input file
# Then substitute in the current variables: cation, anion, functional, basis set
cd $dir_systempreparation
sed 's+SED_FUNCTIONAL_SED+'$functional_name'+g' ionic-pair-opt.inp > $cation_name$anion_name.inp
sed -i 's/SED_ANION_SED/'$anion_name'/g' $cation_name$anion_name.inp
sed -i 's/SED_CATION_SED/'$cation_name'/g' $cation_name$anion_name.inp
sed -i 's/SED_BASISSET_SED/'$basisset_name'/g' $cation_name$anion_name.inp
sed -i 's/___/ /g' $cation_name$anion_name.inp

# Move the resulting input file and pre-existing geometry to the calculation directory
cp $dir_systempreparation/xyz/$cation_name$anion_name.xyz $dir_experiments/$fullpath/
mv $cation_name$anion_name.inp $dir_experiments/$fullpath/
cd $dir_experiments/$fullpath/

# Execute the calculation
#$dir_orca/orca $cation_name$anion_name.inp > $cation_name$anion_name.out

