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

# Define suffix for task
task_name="clean"

# Transfer list entries to bash variables
ion_name=$(echo ${current_ion[0]} | awk '{print $1}')
ion_size=$(echo ${current_ion[0]} | awk '{print $2}')
ion_dist=$(echo ${current_ion[0]} | awk '{print $3}')
slab_name=$(echo ${current_slab[0]} | awk '{print $1}')
slab_face=$(echo ${current_slab[0]} | awk '{print $2}')
slab_size=$(echo ${current_slab[0]} | awk '{print $3}')
slab_vacm=$(echo ${current_slab[0]} | awk '{print $4}')
method_soft=$(echo ${current_method[0]} | awk '{print $1}')
method_mode=$(echo ${current_method[0]} | awk '{print $2}')
method_func=$(echo ${current_method[0]} | awk '{print $3}')
method_disp=$(echo ${current_method[0]} | awk '{print $4}')
method_ctff=$(echo ${current_method[0]} | awk '{print $5}')
method_grid=$(echo ${current_method[0]} | awk '{print $6}')
method_kden=$(echo ${current_method[0]} | awk '{print $7}')
version_name=$(echo ${current_version[0]} | awk '{print $1}')

#These variables are now accessed as in the following example:
echo
echo SYSMEM PARAMETERS
echo ion name $ion_name
echo ion size $ion_size
echo ion dist $ion_dist
echo slab name $slab_name
echo slab face $slab_face
echo slab size $slab_size
echo slab vacuum $slab_vacm
echo soft $method_soft
echo mode $method_mode
echo functional $method_func
echo dispersion $method_disp
echo PW cut-off $method_ctff
echo grid-space $method_grid
echo kp-density $method_kden
echo version $version_name

file_name=$(echo ${slab_name}_${slab_face}_${ion_size},${slab_size}_${ion_name})

# Define path for storing computational data
fullpath=$(echo $slab_name/$slab_face/${ion_size},${slab_size}/$ion_name/$method_soft-$method_mode/$method_func-$method_disp/$method_ctff-$method_grid-$method_kden/$version_name/)
mkdir -p $dir_experiments/$fullpath
echo $dir_experiments/$fullpath
cd $dir_experiments/$fullpath/

#shopt -s extglob
#rm !(*.traj)

#rm bader.p
#rm *.den
#rm *.pot
#rm *_ddec.xyz
#rm *.nrg
#rm core.*
#rm *_sys.gpw

echo CLEANING DONE!
