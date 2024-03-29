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
task_name="prop"

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

# Make new input file by substituting in variables to the template input file
# Then substitute in the current variables: molecule, vacuum, functional, dfttype, mode, param, version
cd $dir_systempreparation
sed 's+SED_ION_SED+'$ion_name'+g' scripts/${task_name}.py > ${file_name}_${task_name}.py
sed -i 's/SED_DISTANCE_SED/'$ion_dist'/g' ${file_name}_${task_name}.py
sed -i 's/SED_XY_SED/'$ion_size'/g' ${file_name}_${task_name}.py
sed -i 's/SED_SLAB_SED/'$slab_name'/g' ${file_name}_${task_name}.py
sed -i 's/SED_Z_SED/'$slab_size'/g' ${file_name}_${task_name}.py
sed -i 's/SED_FACE_SED/'$slab_face'/g' ${file_name}_${task_name}.py
sed -i 's/SED_VACUUM_SED/'$slab_vacm'/g' ${file_name}_${task_name}.py
sed -i 's/SED_FUNCTIONAL_SED/'$method_func'/g' ${file_name}_${task_name}.py
sed -i 's/SED_VDW_SED/'$method_disp'/g' ${file_name}_${task_name}.py
sed -i 's/SED_CUTOFF_SED/'$method_ctff'/g' ${file_name}_${task_name}.py
sed -i 's/SED_GRIDSPACE_SED/'$method_grid'/g' ${file_name}_${task_name}.py
sed -i 's/SED_KPDENSITY_SED/'$method_kden'/g' ${file_name}_${task_name}.py

#sed 's+SED_ION_SED+'$ion_name'+g' scripts/${task_name}2.py > ${file_name}_${task_name}2.py
#sed -i 's/SED_DISTANCE_SED/'$ion_dist'/g' ${file_name}_${task_name}2.py
#sed -i 's/SED_XY_SED/'$ion_size'/g' ${file_name}_${task_name}2.py
#sed -i 's/SED_SLAB_SED/'$slab_name'/g' ${file_name}_${task_name}2.py
#sed -i 's/SED_Z_SED/'$slab_size'/g' ${file_name}_${task_name}2.py
#sed -i 's/SED_FACE_SED/'$slab_face'/g' ${file_name}_${task_name}2.py
#sed -i 's/SED_VACUUM_SED/'$slab_vacm'/g' ${file_name}_${task_name}2.py
#sed -i 's/SED_FUNCTIONAL_SED/'$method_func'/g' ${file_name}_${task_name}2.py
#sed -i 's/SED_VDW_SED/'$method_disp'/g' ${file_name}_${task_name}2.py
#sed -i 's/SED_CUTOFF_SED/'$method_ctff'/g' ${file_name}_${task_name}2.py
#sed -i 's/SED_GRIDSPACE_SED/'$method_grid'/g' ${file_name}_${task_name}2.py
#sed -i 's/SED_KPDENSITY_SED/'$method_kden'/g' ${file_name}_${task_name}2.py

sed 's+SED_PATH_SED+'$dir_systempreparation/atomic_densities/'+g' scripts/job_control.txt > job_control.txt

# Move the resulting input script to the calculation directory
mv ${file_name}_${task_name}.py $dir_experiments/$fullpath/
#mv ${file_name}_${task_name}2.py $dir_experiments/$fullpath/
mv job_control.txt $dir_experiments/$fullpath/
cd $dir_experiments/$fullpath/

# Execute the calculation
mpirun gpaw python ${file_name}_${task_name}.py

# Run DDEC analysis
if [ -f "${file_name}_sys.cube" ]; then
    cp ${file_name}_sys.cube total_density.cube
    chargemol
    mv DDEC6_even_tempered_net_atomic_charges.xyz ${file_name}_ddec.xyz
    rm total_density.cube
    rm total_cube_DDEC_analysis.output
    rm overlap_populations.xyz
    rm DDEC*
fi

if [ -f "${file_name}_cdft_sys.cube" ]; then
    cp ${file_name}_cdft_sys.cube total_density.cube
    chargemol
    mv DDEC6_even_tempered_net_atomic_charges.xyz ${file_name}_cdft_ddec.xyz
    rm total_density.cube
    rm total_cube_DDEC_analysis.output
    rm overlap_populations.xyz
    rm DDEC*
fi

if [ -f "${file_name}_lcao_sys.cube" ]; then
    cp ${file_name}_lcao_sys.cube total_density.cube
    chargemol
    mv DDEC6_even_tempered_net_atomic_charges.xyz ${file_name}_lcao_ddec.xyz
    rm total_density.cube
    rm total_cube_DDEC_analysis.output
    rm overlap_populations.xyz
    rm DDEC*
fi
