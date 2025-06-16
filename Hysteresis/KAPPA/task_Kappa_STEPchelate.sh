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
zbox_nm=$(echo ${current_electrode[0]} | awk '{print $4}') #in nm
r_wall_nm=$(echo ${current_electrode[0]} | awk '{print $5}') #in nm
electrodeatoms=$(echo ${current_electrode[0]} | awk '{print $6}')

ion_name=$(echo ${current_ion[0]} | awk '{print $1}')
r_ion_nm=$(echo ${current_ion[0]} | awk '{print $2}') #in nm
ion_r_name=$(echo ${current_ion[0]} | awk '{print $3}')
ion_charge_name=$(echo ${current_ion[0]} | awk '{print $4}')
ion_centre_name=$(echo ${current_ion[0]} | awk '{print $5}')

numberofions_name=$(echo ${current_numberofions[0]} | awk '{print $1}') 
temperature_name=$(echo ${current_temperature[0]} | awk '{print $1}')
version_name=$(echo ${current_version[0]} | awk '{print $1}')
replica_name=$(echo ${current_replica[0]} | awk '{print $1}')

####################################
# Calculations
xbox=$(awk "BEGIN {print "$xbox_nm*10.0"}" /dev/null)  #in A
ybox=$(awk "BEGIN {print "$ybox_nm*10.0"}" /dev/null)  #in A
zbox=$(awk "BEGIN {print "$zbox_nm*10.0"}" /dev/null)  #in A
xbox12=$(awk "BEGIN {print "$xbox*0.5"}" /dev/null)
ybox12=$(awk "BEGIN {print "$ybox*0.5"}" /dev/null)
zbox12=$(awk "BEGIN {print "$zbox*0.5"}" /dev/null)
r_wall=$(awk "BEGIN {print "$r_wall_nm*10.0"}" /dev/null)   #in A
z1_ion=$(awk "BEGIN {print "$zbox12"}" /dev/null)   #in A
z2_ion=$(awk "BEGIN {print "$zbox12+4.0*$r_ion_nm*10.0"}" /dev/null) #in A

A_box=$(awk "BEGIN {print "$xbox*$ybox"}" /dev/null)	#box area in sqA
A_ion=$(awk "BEGIN {print "4*$r_ion_nm*$r_ion_nm*100"}" /dev/null)	#ion area in sqA
maxnumberofions_name=$(awk "BEGIN {print int("$A_box/$A_ion")}" /dev/null)
varnumberofions_name=$(awk "BEGIN {print "$maxnumberofions_name+1.0*$numberofions_name"}" /dev/null)


electrode_atom_charge=$(awk "BEGIN {printf \"%.20f\n\","-1.0*$ion_charge_name*$varnumberofions_name/$electrodeatoms"}" /dev/null)

####################################

# Define path for storing simulation data
fullpath_start=$electrode_name/$ion_name/$temperature_name/$version_name/$replica_name/$varnumberofions_name

if [ -s $dir_experiments/$fullpath_start/NVT_lowtimestep.gro ]
then
echo $fullpath_start
echo $fullpath_start
else
echo $fullpath_start
if [ -s $dir_experiments/$fullpath_start/STEEP.gro ]
then
grompp -f 1_NVT_lowtimestep.mdp -c STEEP.gro -p topol_local.top -n index.ndx -o NVT_lowtimestep -maxwarn 1
rm mdout.mdp

mpirun gmx mdrun -ntmpi 3 -deffnm NVT_lowtimestep
rm *trr
else
mkdir -p $dir_experiments/$fullpath_start

cd $dir_systempreparation

# Start with packmol: Add all particle names and calculated values concerning the box size to the packmpl script.

echo 'Start converting the packmol input file packmol.inp ...'
sed 's+SED_dir_systempreparation_SED+'$dir_systempreparation'+g' packmol_Kappa.inp > packmol.inp

seeds=( 1983757 3542638 6328548 6679881 7552618 9305733 9361537 39571623 )
seed=${seeds[$replica_name]}
sed -i 's/SED_seed_SED/'$seed'/g' packmol.inp

sed -i 's/SED_ion_name_SED/'$ion_name'/g' packmol.inp
sed -i 's/SED_ion_num_SED/'$varnumberofions_name'/g' packmol.inp

sed -i 's/SED_xbox_SED/'$xbox'/g' packmol.inp
sed -i 's/SED_ybox_SED/'$ybox'/g' packmol.inp
sed -i 's/SED_zbox_left_SED/'$z1_ion'/g' packmol.inp
sed -i 's/SED_zbox_right_SED/'$z2_ion'/g' packmol.inp

sed -i 's/SED_xbox12_SED/'$xbox12'/g' packmol.inp
sed -i 's/SED_ybox12_SED/'$ybox12'/g' packmol.inp
sed -i 's/SED_zbox12_SED/'$zbox12'/g' packmol.inp

sed -i 's/SED_electrode_name_SED/'$electrode_name'/g' packmol.inp


echo 'Run packmol and convert packmol output to gromacs input ...'
$dir_packmol/packmol < packmol.inp
editconf -f packmol.pdb -o packmol.gro
rm packmol.pdb packmol.inp

#read -p "Press enter to continue..."

#+++++++++++++++++++++++
#Edit the box size to insert the vacuum slab
sed -i '$d' packmol.gro
echo $xbox_nm $ybox_nm $zbox_nm >> packmol.gro

#Prepare the index file
make_ndx -f packmol.gro -o index.ndx << EOF
keep 0
r $electrode_name
name 1 Electrode
r $ion_r_name
name 2 Ion
a $ion_centre_name
name 3 Centre
q
EOF

#Add all necessary .mdp files
sed 's/SED_temperature_name_SED/'$temperature_name'/g' 0_STEEP_varTemp.mdp > 0_STEEP.mdp
sed 's/SED_temperature_name_SED/'$temperature_name'/g' 1_NVT_lowtimestep_varTemp.mdp > 1_NVT_lowtimestep.mdp

echo 'Start converting the topology files ...'
cp $dir_systempreparation/top/$electrode_name.itp electrode_local.itp
sed -i 's/SED_electrode_charge_SED/'$electrode_atom_charge'/g' electrode_local.itp

sed 's+SED_dir_systempreparation_SED+'$dir_systempreparation'+g' $dir_systempreparation/topol_local_chelate.top > topol_local.top
sed -i 's/SED_ion_name_SED/'$ion_name'/g' topol_local.top
sed -i 's/SED_ion_r_name_SED/'$ion_r_name'/g' topol_local.top
sed -i 's/SED_ion_num_SED/'$varnumberofions_name'/g' topol_local.top
sed -i 's/SED_electrode_SED/'$electrode_name'/g' topol_local.top
sed -i 's/SED_electrodeatoms_SED/'$electrodeatoms'/g' topol_local.top
#Move everything to the rundirectory
mv 0_STEEP.mdp 1_NVT_lowtimestep.mdp electrode_local.itp topol_local.top packmol.gro index.ndx $dir_experiments/$fullpath_start/

####################################
cd $dir_experiments/$fullpath_start/

#grompp and run the energy minimization
grompp -f 0_STEEP.mdp -c packmol.gro -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
rm mdout.mdp
mpirun gmx mdrun -ntmpi 3 -deffnm STEEP

#grompp -f 0_STEEP.mdp -c STEEP.gro -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
#rm mdout.mdp
#mdrun -deffnm STEEP

#grompp -f 0_STEEP.mdp -c STEEP.gro -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
#rm mdout.mdp
#mdrun -deffnm STEEP

#grompp -f 0_STEEP.mdp -c STEEP.gro -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
#rm mdout.mdp
#mdrun -deffnm STEEP

#grompp -f 0_STEEP.mdp -c STEEP.gro -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
#rm mdout.mdp
#mdrun -deffnm STEEP
rm *#

#grompp the first equilibration step

grompp -f 1_NVT_lowtimestep.mdp -c STEEP.gro -p topol_local.top -n index.ndx -o NVT_lowtimestep -maxwarn 1
rm mdout.mdp

mpirun gmx mdrun -ntmpi 3 -deffnm NVT_lowtimestep
rm *trr
fi
#rm *xtc
echo $electrode_atom_charge

echo $r_wall
fi
