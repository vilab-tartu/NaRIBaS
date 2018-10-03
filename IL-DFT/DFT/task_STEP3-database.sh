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
cation_atom_num=$(echo ${current_cation[0]} | awk '{print $2}')
functional_name=$(echo ${current_functional_sp[0]} | awk '{print $1}' | sed 's,\/,\\/,g')
basisset_name=$(echo ${current_basisset[0]} | awk '{print $1}' | sed 's,\/,\\/,g')
version_name=$(echo ${current_version[0]} | awk '{print $1}')

echo SYSMEM PARAMETERS
echo Anion $anion_name
echo Cation $cation_name
echo Functional $functional_name
echo Basis-set $basisset_name

# Define path for storing configurations / computational data
fullpath=$(echo $cation_name$anion_name/$functional_name/ | sed 's,\\\/,,g')
cd $dir_experiments/$fullpath/

# Copy the database file and its creation script to current folder
mv $dir_experiments/db.json $dir_experiments/$fullpath/ # Unified results text moved to data.
cp $dir_experiments/db.py $dir_experiments/$fullpath/ # db default.

# Make modifications on the database generation script
sed -i 's/SED_ANION_SED/'$anion_name'/g' $dir_experiments/$fullpath/db.py
sed -i 's/SED_CATION_SED/'$cation_name'/g' $dir_experiments/$fullpath/db.py
sed -i 's/SED_FUNCTIONAL_SED/'$functional_name'/g' $dir_experiments/$fullpath/db.py
sed -i 's/SED_BASISSET_SED/'$basisset_name'/g' $dir_experiments/$fullpath/db.py
sed -i 's/SED_VERSION_SED/'$version_name'/g' $dir_experiments/$fullpath/db.py
sed -i 's/___/ /g' $dir_experiments/$fullpath/db.py

# Execute the script, which indexes the properties of interest. That script template is found in $dir_experiments
python $dir_experiments/$fullpath/db.py 

# Clean up: remove the modified version of the script and overwrite the old database with the new one
rm $dir_experiments/$fullpath/db.py
mv $dir_experiments/$fullpath/db.json $dir_experiments/ 
